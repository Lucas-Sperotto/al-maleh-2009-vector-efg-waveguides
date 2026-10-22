#!/usr/bin/env python3
"""Al Maleh (2009), Chapter 5 Cartesian EFG: source-anchored benchmark.

Independent current-method replication, NOT a claim of numerical identity:
MLS weighting, support and exact author's geometry partition parameters beyond
the quoted Figure 9 are not fully specified by the supplied pages.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import subprocess
import sys
from m2_ch5_quiver_svg import render

DEFAULT_REFERENCE = Path("docs/results/m2_author_ch5_tables_2_3.csv")
MODES = {"TM11": (1,1), "TM21": (2,1)}
EXACT = {"TM11": 1.25*math.pi**2, "TM21": 2*math.pi**2}

def reference_rows(file):
    with file.open(encoding="utf-8",newline="") as stream:
        rows=list(csv.DictReader(stream))
    expected={"mode","nodes_y","nodes_x","gauss_points_per_axis",
              "author_signed_error_percent","author_table","source_equation",
              "material_source"}
    if not rows or not expected.issubset(rows[0]):
        raise ValueError("Incomplete original Chapter 5 table reference")
    seen=set()
    for r in rows:
        key=(r["mode"],int(r["nodes_y"]),int(r["nodes_x"]),
             int(r["gauss_points_per_axis"]))
        if key in seen or r["mode"] not in MODES or key[-1] not in range(1,6):
            raise ValueError("Invalid/duplicate Al Maleh reference row")
        seen.add(key)
        if not math.isfinite(float(r["author_signed_error_percent"])):
            raise ValueError("Nonfinite author error value")
        if int(r["author_table"])!=(2 if r["mode"]=="TM11" else 3):
            raise ValueError("Al Maleh table index does not match mode")
    return rows

def read_kv(path):
    return dict(line.split("=",1) for line in
                path.read_text(encoding="utf-8").splitlines() if "=" in line)

def read_modes(path):
    with path.open(encoding="utf-8",newline="") as stream:
        rows=list(csv.DictReader(stream))
    result={}
    for row in rows:
        key=(int(row["m"]),int(row["n"]))
        if key in MODES.values():
            result[key]=row
    if set(result)!=set(MODES.values()):
        raise ValueError("Missing TM11 or TM21 diagnostics")
    return result

def signed_error_pct(lam,exact):
    return 100*(lam/exact-1)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--binary",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    p.add_argument("--reference",default=DEFAULT_REFERENCE,type=Path)
    p.add_argument("--support-factor",default=3.3,type=float)
    p.add_argument("--subdivision",default=2,type=int)
    # Default keeps the historical profile; quartic_spline is Eq. (4-2).
    p.add_argument("--weight",default="wendland_c2",
                   choices=("wendland_c2","quartic_spline"))
    p.add_argument("--timeout",default=240,type=int)
    args=p.parse_args()
    if not args.binary.is_file() or not args.reference.is_file():
        p.error("Missing executable or original table CSV")
    if not math.isfinite(args.support_factor) or args.support_factor<=0:
        p.error("Invalid MLS support")
    if args.subdivision!=2:
        p.error("The Chapter 5 Figure 9 experimental layout uses subdivision=2")
    if args.timeout<1:
        p.error("Invalid timeout")
    ref=reference_rows(args.reference)
    combinations=sorted({(int(r["nodes_x"]),int(r["nodes_y"]),
                          int(r["gauss_points_per_axis"])) for r in ref})
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
    results={}
    for nx,ny,gauss in combinations:
        name=f"{ny}x{nx}_gauss{gauss}"
        folder=out/name;folder.mkdir(parents=True,exist_ok=True)
        cmd=[str(args.binary.resolve()),str(nx),str(ny),"1",
             str(args.support_factor),str(gauss),str(folder),
             str(args.subdivision),args.weight]
        print("M2 CH5",name,"command:",*cmd,flush=True)
        proc=subprocess.run(cmd,text=True,capture_output=True,
                            check=False,timeout=args.timeout)
        (folder/"stdout.txt").write_text(proc.stdout,encoding="utf-8")
        (folder/"stderr.txt").write_text(proc.stderr,encoding="utf-8")
        if proc.returncode:
            raise RuntimeError(f"Failed benchmark {name}: {proc.stderr}")
        meta=read_kv(folder/"run_metadata.txt")
        expected_cells=(nx-1)*(ny-1)*args.subdivision**2
        expected_points=expected_cells*gauss**2
        if (int(meta["background_cells"])!=expected_cells or
            int(meta["quadrature_points"])!=expected_points or
            int(meta["subcells_per_node_interval"])!=2 or
            int(meta["dofs"])!=2*nx*ny):
            raise RuntimeError(f"Incorrect Ch5 integration topology: {name}")
        modes=read_modes(folder/"tm_matches.csv")
        for fig_mode in ("tm11","tm21"):
            render(folder/(fig_mode+"_vector_field.csv"),
                   folder/(fig_mode+"_vector_field.svg"),
                   f"EFG numérico - {fig_mode.upper()} - malha {ny}x{nx} - Gauss {gauss}x{gauss}")
        results[(nx,ny,gauss)]=(modes,meta)
    summary=[]
    for r in ref:
        mode=r["mode"]; nx=int(r["nodes_x"]);ny=int(r["nodes_y"])
        gauss=int(r["gauss_points_per_axis"])
        eig=results[(nx,ny,gauss)][0][MODES[mode]]
        exact=EXACT[mode]; author=float(r["author_signed_error_percent"])
        closest=float(eig["nearest_lambda"])
        selected=float(eig["strongest_lambda"])
        nearest_err=signed_error_pct(closest,exact)
        shape_err=signed_error_pct(selected,exact)
        summary.append(dict(
            mode=mode,nodes_y=ny,nodes_x=nx,gauss_points_per_axis=gauss,
            author_lambda_rounded=12.3370 if mode=="TM11" else 19.7392,
            analytic_lambda_current_precision=exact,
            author_signed_error_percent=author,
            current_nearest_lambda=closest,
            current_nearest_signed_error_percent=nearest_err,
            current_strongest_curl_lambda=selected,
            current_strongest_curl_signed_error_percent=shape_err,
            current_strongest_curl_overlap=float(eig["strongest_overlap"]),
            current_nearest_curl_overlap=float(eig["nearest_curl_overlap"]),
            current_divergence_ratio=float(eig["divergence_ratio"]),
            current_boundary_normal_ratio=float(eig["normal_boundary_ratio"]),
            absolute_difference_nearest_error_percentage_points=abs(nearest_err-author),
            absolute_difference_strongest_error_percentage_points=abs(shape_err-author),
            reproduction_status="different_MLS_and_quadrature_details_unverified"))
    fields=list(summary[0])
    with (out/"comparison.csv").open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=fields)
        writer.writeheader();writer.writerows(summary)
    setup=dict(source="Hussam Al Maleh (2009), Chapter 5, equations (5-1)-(5-11), "
                      "Figures 9-15, Tables 2 and 3, supplied screenshots",
               source_reference=str(args.reference),source_numeric_rows=len(ref),
               experiments=len(combinations),
               domain={"a_m":2,"b_m":1},node_naming="author_nodes_y_x",
               current_grid_nx_is_x_not_author_first_dimension=True,
               subcells_per_nodal_interval=args.subdivision,
               integration_cells_per_node_rectangle=args.subdivision**2,
               polynomial_degree=1,wendland_c2=args.weight=="wendland_c2",
               mls_weight=args.weight,
               support_factor=args.support_factor,
               mode_selection="closest_lambda_and_max_curl_overlap_recorded_separately",
               comparison_kind="exploratory_not_exact_replication")
    (out/"metadata.json").write_text(json.dumps(setup,indent=2,ensure_ascii=False)+"\n",
                                     encoding="utf-8")
    report=[
        "# M2 — Capítulo 5 de Al Maleh: comparação documentada (cartesiano)",
        "",
        "**Os percentuais do autor foram transcritos das imagens fornecidas.**",
        f"Os cálculos locais usam MLS linear com peso {args.weight}, fator de suporte",
        f"{args.support_factor}, geometria 2×1 e quatro células por retângulo nodal.",
        "Sem conhecimento dos parâmetros de suporte, pesos MLS e estratégias",
        "originais completos, diferenças NÃO significam automaticamente bugs.",
        "A coincidência numérica deve ser demonstrada, não assumida.",
        "",
        "O formato 3×5 do autor significa três posições no eixo y e cinco no x;",
        "no programa, os argumentos dessa nuvem são (nx=5, ny=3).",
        "Os percentuais são erros RELATIVOS ASSINADOS em λ=k²; os modos",
        "mais próximos espectralmente podem ser espúrios. Por isso também",
        "exportamos o candidato com maior correlação espacial curl(H).",
        "",
        "| Modo | Grade autor | Gauss | Autor (%) | Nosso mais próximo (%) | Nosso curl mais forte (%) | Δ pp mais próximo |",
        "|---|---|---:|---:|---:|---:|---:|",
        *[f"| {r['mode']} | {r['nodes_y']}×{r['nodes_x']} | "
          f"{r['gauss_points_per_axis']}×{r['gauss_points_per_axis']} | "
          f"{r['author_signed_error_percent']:+.3f} | "
          f"{r['current_nearest_signed_error_percent']:+.3f} | "
          f"{r['current_strongest_curl_signed_error_percent']:+.3f} | "
          f"{r['absolute_difference_nearest_error_percentage_points']:.3f} |"
          for r in summary],
        "",
        "Referências analíticas: TM11 = 1.25*pi² ≈ 12.33700550;",
        "TM21 = 2*pi² ≈ 19.73920880.",
        "Na dissertação, Tabela 2 é TM11 e Tabela 3 é TM21.",
        "Onde há espaços vazios nas tabelas, não se inventaram resultados.",
        "",
        "A Figura 15 mostra um candidato espúrio com k² ≈ 12.1456",
        "(grade 6×11, quadratura 2×2). Essa proximidade ao TM11",
        "não deve ser confundida com correspondência física de campo.",
        "",
        "Arquivos: comparison.csv, metadata.json e uma pasta por combinação",
        "com espectro completo, matrizes, diagnósticos e campos TM11/TM21",
        "em CSV e SVG originais, além de stdout/stderr.",
        "Nenhum espectro histórico M6a–M6f foi alterado.",
        ""
    ]
    (out/"REPORT.md").write_text("\n".join(report),encoding="utf-8")
    print("\n".join(report),flush=True)
    print(f"PASS: {len(combinations)} independent configurations, "
          f"{len(summary)} source-documented comparison rows",flush=True)
    return 0

if __name__=="__main__":
    raise SystemExit(main())
