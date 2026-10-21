#!/usr/bin/env python3
"""Reproduce Al Maleh 2009 Ch6 Tables 5/6 as an independent EFG experiment.

The original author numbers and the solver results remain separate. Signed
errors are measured against k²(TM11)=5pi²/4 and k²(TM21)=2pi². A shape-based
candidate and an eigenvalue-nearest candidate are both reported; neither is
arbitrarily declared physically validated.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
REFERENCE=ROOT/"docs/results/m5_author_ch6_tables_5_6.csv"
CASES={
    (5,3):[1,2,3,4],
    (7,4):[1,2,3,4],
    (9,5):[1,2],
    (11,6):[1],
    (13,7):[1],
}
EXACT={"TM11":1.25*math.pi**2,"TM21":2*math.pi**2}

def kv(path):
    return dict(s.split("=",1) for s in path.read_text(encoding="utf-8").splitlines()
                if "=" in s)

def read_reference(path=REFERENCE):
    with path.open(newline="",encoding="utf-8") as f:
        data=list(csv.DictReader(f))
    if len(data)!=24:
        raise ValueError("Ch6 Tables 5/6 require 24 reported values")
    seen=set()
    for r in data:
        key=(r["mode"],int(r["nodes_x"]),int(r["nodes_y"]),
             int(r["gauss_points_per_axis"]))
        if key in seen or r["mode"] not in EXACT or key[3] not in CASES.get((key[1],key[2]),()):
            raise ValueError("Unexpected/duplicate primary reference row")
        seen.add(key)
        if int(r["author_table"]) != (5 if r["mode"]=="TM11" else 6):
            raise ValueError("Mismatched source table number")
        if not math.isfinite(float(r["author_signed_error_percent"])):
            raise ValueError("Invalid author reference number")
    return data

def render_field(svg_tool,base,output,name):
    # SVG renderer imports no science plotting dependency; output is our
    # original numerical field and never reuses a thesis figure.
    proc=subprocess.run([sys.executable,str(svg_tool),
                         "--input",str(base),
                         "--output",str(output),"--title",name],
                        capture_output=True,text=True,check=False)
    if proc.returncode:
        raise RuntimeError("SVG creation failed: "+proc.stderr)

def run(binary,out,tol,support,timeout):
    ref=read_reference()
    out.mkdir(parents=True,exist_ok=True)
    measured={}
    svg_tool=ROOT/"scripts/m2_ch5_quiver_svg.py"
    for (nx,ny),orders in CASES.items():
        for order in orders:
            name=f"{ny}x{nx}_gauss{order}"
            path=out/name
            path.mkdir(parents=True,exist_ok=True)
            cmd=[str(binary.resolve()),str(nx),str(ny),str(order),
                 str(support),"2",str(path),str(tol)]
            print("M5 CH6",name," ".join(cmd),flush=True)
            proc=subprocess.run(cmd,check=False,capture_output=True,text=True,
                                timeout=timeout)
            (path/"stdout.txt").write_text(proc.stdout,encoding="utf-8")
            (path/"stderr.txt").write_text(proc.stderr,encoding="utf-8")
            if proc.returncode:
                raise RuntimeError(f"Failed {name}: {proc.stderr}")
            metadata=kv(path/"run_metadata.txt")
            assert int(metadata["dofs"])==2*nx*ny
            assert int(metadata["retained_rank"])+int(metadata["removed_mass_nullity"])==2*nx*ny
            assert int(metadata["integration_points"])==4*(nx-1)*(ny-1)*order*order
            with (path/"tm_matches.csv").open(newline="",encoding="utf-8") as f:
                modal={(int(r["m"]),int(r["n"])):r for r in csv.DictReader(f)}
            if set(modal)!={(1,1),(2,1)}:
                raise RuntimeError("Missing TM11/TM21 results")
            measured[(nx,ny,order)]=(metadata,modal)
            for key in ("tm11_vector_field","tm21_vector_field",
                        "tm11_irrot_field","tm11_rot_field"):
                render_field(svg_tool,path/(key+".csv"),path/(key+".svg"),
                             f"M5 {key} {name}")
            if nx==9 and ny==5:
                for key in ("irrotational_shape_central","rotational_shape_central"):
                    render_field(svg_tool,path/(key+".csv"),path/(key+".svg"),
                                 f"M5 Ch6 central MLS shape {key} {name}")
    rows=[]
    for r in ref:
        nx,ny,g=[int(r[n]) for n in ("nodes_x","nodes_y","gauss_points_per_axis")]
        mode=r["mode"]
        data=measured[(nx,ny,g)][1][(1,1) if mode=="TM11" else (2,1)]
        author=float(r["author_signed_error_percent"])
        near=float(data["nearest_signed_error_pct"])
        shape=float(data["best_curl_signed_error_pct"])
        rows.append(dict(
            mode=mode,nodes_y=ny,nodes_x=nx,gauss_points_per_axis=g,
            author_signed_error_percent=author,
            current_nearest_signed_error_percent=near,
            current_best_curl_signed_error_percent=shape,
            current_nearest_lambda=data["nearest_lambda"],
            current_best_curl_lambda=data["best_curl_lambda"],
            current_best_curl_rho=data["best_curl_rho"],
            difference_shape_author_percentage_points=shape-author,
            absolute_difference_shape_author_percentage_points=abs(shape-author),
            mass_nullity=int(measured[(nx,ny,g)][0]["removed_mass_nullity"]),
            retained_rank=int(measured[(nx,ny,g)][0]["retained_rank"]),
            numerical_curl_kernel=int(measured[(nx,ny,g)][0]["near_zero_eigenvalues"]),
            negative_eigenvalues=int(measured[(nx,ny,g)][0]["negative_eigenvalues"]),
            status="experimental_different_MLS_parameters_not_source_identical"))
    with (out/"comparison.csv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    summary=[]
    for (nx,ny,g),(meta,modal) in measured.items():
        summary.append(dict(nodes_y=ny,nodes_x=nx,gauss=g,
            retained_mass_rank=int(meta["retained_rank"]),
            removed_mass_nullity=int(meta["removed_mass_nullity"]),
            numerical_curl_kernel=int(meta["near_zero_eigenvalues"]),
            negative_eigenvalues=int(meta["negative_eigenvalues"]),
            mass_cutoff=float(meta["mass_cutoff"]),
            mass_orthogonality_error=float(meta["M_orthogonality_error"]),
            tm11_best_curl_rho=float(modal[(1,1)]["best_curl_rho"]),
            tm21_best_curl_rho=float(modal[(2,1)]["best_curl_rho"])))
    with (out/"spectrum_summary.csv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(summary[0]))
        writer.writeheader();writer.writerows(summary)
    provenance=dict(
        reference="Al Maleh 2009 Ch6 Eqs 6-1..6-6, Tables 5 and 6",
        reference_rows=len(ref),executed_configs=len(measured),
        alpha="gradient MLS",beta="z_cross_gradient MLS",
        prior_ch5_subdivision=2,original_author_grid_indexing="ny_x_nx",
        support_factor=support,mls_weight="Wendland C2",
        mass_relative_rank_cutoff=tol,
        negative_eigenvalues_retained_in_raw_exports=True,
        observed_not_an_exact_numerical_replication=True)
    (out/"metadata.json").write_text(json.dumps(provenance,indent=2)+"\n",
                                      encoding="utf-8")
    report=[
      "# M5 — Al Maleh (2009) Capítulo 6: campanha vetorial rotacional/irrotacional",
      "",
      "**Valores publicados e valores reproduzidos são dados distintos.**",
      "Eqs. (6-1)..(6-6) implementadas: grad Ni e z×grad Ni;",
      "o autoproblema usa as formas fracas curl-curl/massa do Capítulo 5.",
      "A massa singular é reduzida por autovetores positivos acima de",
      "um limiar declarado; o espectro retido inclui todos os autovalores",
      "negativos e quase zero, sem reclassificação seletiva.",
      "",
      "| Modo | Grade | Gauss | Autor (%) | Candidato λ próximo (%) | Melhor forma curl (%) | Diferença forma/autor (pp) |",
      "|---|---|---|---:|---:|---:|---:|",
      *[
        f"| {r['mode']} | {r['nodes_y']}×{r['nodes_x']} | {r['gauss_points_per_axis']}×{r['gauss_points_per_axis']} | "
        f"{r['author_signed_error_percent']:+.3f} | "
        f"{r['current_nearest_signed_error_percent']:+.3f} | "
        f"{r['current_best_curl_signed_error_percent']:+.3f} | "
        f"{r['difference_shape_author_percentage_points']:+.3f} |"
        for r in rows
      ],
      "",
      "## Rank e número de modos",
      "",
      "| Malha | Gauss | Rank M | Nullity M | Zeros curl | Autovalores negativos |",
      "|---|---:|---:|---:|---:|---:|",
      *[
        f"| {r['nodes_y']}×{r['nodes_x']} | {r['gauss']} | {r['retained_mass_rank']} | "
        f"{r['removed_mass_nullity']} | {r['numerical_curl_kernel']} | {r['negative_eigenvalues']} |"
        for r in summary
      ],
      "",
      "Os zeros de massa (dependências entre funções) NÃO são os mesmos",
      "que os autovalores zero de rigidez no problema reduzido.",
      "Erro pequeno em λ sozinho não identifica um modo físico.",
      "",
      "## Referência da dissertação",
      "",
      "Tabela 4 descreve muitos autovalores zero no método rotacional/irrotacional;",
      "e resultados negativos na quadratura 2×2, reduzidos por ordens mais altas.",
      "Para TM11, as três raízes destacadas são aproximadamente 16,8, 15,8 e",
      "15,5 (Gauss 2×2, 3×3, 4×4) no caso 3×5, conforme Tabela 4.",
      "As Tabelas 5 e 6 apresentam erros sem convergência sistemática;",
      "as Figuras 18–26 exibem campos irregulares e modos espúrios.",
      "",
      "Nenhum parâmetro foi ajustado para reproduzir o resultado negativo",
      "da dissertação. Diferenças quantitativas exigem cotejo com o capítulo 4",
      "(MLS original e integrais), seleção modal e sensibilidade da redução",
      "de posto. O driver preserva no CSV bruto o espectro retido completo",
      "e todo o conjunto de coeficientes dos dois potenciais.",
      "",
      "As saídas tm11_{irrot,rot}_field.csv/.svg mostram separadamente",
      "os dois termos da Eq. (6-6), nunca copiam as Figuras 22–23.",
      "",
    ]
    (out/"REPORT.md").write_text("\n".join(report),encoding="utf-8")
    print("\n".join(report),flush=True)
    print(f"PASS: {len(measured)} configurations and {len(rows)} source points",flush=True)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--binary",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--support-factor",type=float,default=3.3)
    p.add_argument("--mass-rank-relative-cutoff",type=float,default=1e-10)
    p.add_argument("--timeout",type=int,default=600)
    a=p.parse_args()
    if not a.binary.is_file():p.error("Missing rot_irrot_rect compiled binary")
    run(a.binary,a.out.resolve(),a.mass_rank_relative_cutoff,a.support_factor,a.timeout)

if __name__=="__main__":
    main()
