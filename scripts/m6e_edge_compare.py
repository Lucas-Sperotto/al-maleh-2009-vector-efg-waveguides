#!/usr/bin/env python3
"""M6e: reproducible EFG x independent Whitney edge FEM TM comparison.

No numerical eigenpairs are modified; the external edge_rect solver is invoked
at a pinned source revision and its raw eigenspectrum is kept intact.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess

FEM_REPO = "Lucas-Sperotto/fem-eigen-em-2d3d-whitney"
PINNED_FEM_SHA = "dc1d9545c1238554c2292aa0f64a3c1b883d4118"
EFG_CASES = (
    ("9x7", 8, 4, (54, 60, 63, 66, 70, 71), 49),
    ("13x11", 12, 6, (100, 108, 117, 123, 127, 128), 94),
    ("17x15", 16, 8, (153, 163, 178, 185, 192, 193), 144),
)
MODE_NAMES = ("TM11", "TM21", "TM31", "TM12", "TM22_A", "TM41_B")
PI = math.pi
ANALYTIC = (
    1.25 * PI**2, 2 * PI**2, 3.25 * PI**2, 4.25 * PI**2,
    5 * PI**2, 5 * PI**2,
)
EFG_GEOMETRY = (2.0, 1.0)
EDGE_GEOMETRY = (1.0, 0.5)
GEOMETRIC_LAMBDA_MULTIPLIER = (EDGE_GEOMETRY[0]/EFG_GEOMETRY[0])**2
# Lambda on a geometrically enlarged guide scales with inverse length squared:
# lambda_EFG = lambda_edge * (a_edge/a_EFG)^2 = lambda_edge / 4.

def read_csv(path):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))

def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b""):
            h.update(chunk)
    return h.hexdigest()

def scaled_lambda(lambda_edge):
    if not math.isfinite(lambda_edge):
        raise ValueError("Nonfinite FEM edge eigenvalue")
    return lambda_edge * GEOMETRIC_LAMBDA_MULTIPLIER

def validate_edge_exports(directory):
    directory=Path(directory)
    tm=directory/"linop"/"edge_rect_tm_eigenvalues.csv"
    te=directory/"linop"/"edge_rect_te_eigenvalues.csv"
    matched=directory/"csv"/"edge_rect_modes.csv"
    timing=directory/"run_timing.csv"
    if not all(p.is_file() for p in (tm,te,matched,timing)):
        raise RuntimeError("Missing one or more native FEM edge reference exports: "+str(directory))
    rows=read_csv(tm)
    tm_lambdas=[float(r["lambda"]) for r in rows]
    if not tm_lambdas or any(int(r["ordered_rank"]) != i+1 for i,r in enumerate(rows)):
        raise RuntimeError("Incomplete or reordered native edge TM spectrum")
    if any(not math.isfinite(x) for x in tm_lambdas) or any(
            a>b+1e-12 for a,b in zip(tm_lambdas,tm_lambdas[1:])):
        raise RuntimeError("Unsorted or nonfinite edge TM eigenvalues")
    return rows, tm_lambdas, read_csv(matched), read_csv(timing)

def match_near_analytic(spectrum, window=0.25):
    """Stable unique matching of six physical roots by proximity to analytic λ.

    Degenerate TM22/TM41 is matched as a 2-D spectral PAIR, not by eigenvector
    orientation. Nearby eigenvalues are not assumed to have correct fields.
    """
    usable=[(i,lam) for i,lam in enumerate(spectrum) if lam>1e-8]
    selected={}
    occupied=set()
    for mode in range(4):
        exact=ANALYTIC[mode]
        candidates=sorted(((abs(v-exact)/exact,i,v) for i,v in usable
                            if i not in occupied))
        if candidates and candidates[0][0]<=window:
            err,i,lam=candidates[0]
            selected[mode]=(i,lam,err)
            occupied.add(i)
        else:
            selected[mode]=None
    exact=ANALYTIC[4]
    doublet=sorted(((abs(v-exact)/exact,i,v) for i,v in usable
                    if i not in occupied))[:2]
    if len(doublet)!=2 or any(e>window for e,_,_ in doublet):
        selected[4]=None
        selected[5]=None
    else:
        doublet.sort(key=lambda x:x[2])
        for k,(_,i,lam) in enumerate(doublet):
            selected[4+k]=(i,lam,abs(lam-exact)/exact)
    return selected

def run_edge(binary, nx, ny, directory, timeout):
    env=os.environ.copy()
    # Explicit isolated experiment out root; the FEM executable honors this.
    env["TP3485_OUT_DIR"]=str(directory.resolve())
    env["OPENBLAS_NUM_THREADS"]="1"
    env["OMP_NUM_THREADS"]="1"
    cmd=[str(binary.resolve()),"--nx",str(nx),"--ny",str(ny),
         "--nmodos","12","--backend","closed-form"]
    proc=subprocess.run(cmd,env=env,capture_output=True,text=True,
                        check=False,timeout=timeout)
    directory.mkdir(parents=True,exist_ok=True)
    (directory/"edge_stdout.txt").write_text(proc.stdout,encoding="utf-8")
    (directory/"edge_stderr.txt").write_text(proc.stderr,encoding="utf-8")
    if proc.returncode:
        raise RuntimeError("Whitney FEM edge_rect failed: "+proc.stderr[-2500:])
    return directory/"helmvec"/"rect",cmd

def compare_case(edge_dir,efg_dir,mode_indices,expected_below):
    edge_rows,edge_native,matched,timing=validate_edge_exports(edge_dir)
    efg_rows=read_csv(efg_dir/"eigenvalues.csv")
    efg_values=[float(r["lambda"]) for r in efg_rows]
    if max(mode_indices)>=len(efg_values):
        raise RuntimeError("Archived EFG benchmark indices outside spectrum")
    if any(int(r["index"])!=i for i,r in enumerate(efg_rows)):
        raise RuntimeError("Incomplete EFG original eigenvalue indexes")
    scaled=[scaled_lambda(v) for v in edge_native]
    positive_cutoff_edge=max(1e-9,1e-10*max(edge_native))
    # State cutoff explicitly. Values slightly below the analytic first root
    # may approximate TM11; call them "below reference", not "spurious".
    below_analytic=sum(positive_cutoff_edge<v<ANALYTIC[0] for v in scaled)
    deeply_below=sum(positive_cutoff_edge<v<0.5*ANALYTIC[0] for v in scaled)
    original_count=sum(
        r["status"]=="positive_candidate" and 0<float(r["lambda"])<ANALYTIC[0]
        for r in efg_rows
    )
    # "near zero" EFG depends on the original published classification.
    if original_count!=expected_below:
        raise RuntimeError("Historical EFG low-positive counts changed")
    assignments=match_near_analytic(scaled)
    comparisons=[]
    for m,name in enumerate(MODE_NAMES):
        i=mode_indices[m]
        efg=efg_values[i]
        exact=ANALYTIC[m]
        selected=assignments[m]
        comparisons.append({
            "mode":name,"analytic_lambda_on_2x1":exact,
            "efg_index":i,"efg_lambda":efg,
            "efg_relative_error":abs(efg-exact)/exact,
            "edge_index":"" if selected is None else selected[0],
            "edge_native_lambda_on_1x0_5":"" if selected is None else edge_native[selected[0]],
            "edge_scaled_lambda_on_2x1":"" if selected is None else selected[1],
            "edge_relative_error":"" if selected is None else selected[2],
            "efg_edge_relative_difference":"" if selected is None else abs(efg-selected[1])/exact,
            "matching_method":"no_match" if selected is None else (
                "pair_by_eigenvalue_only_not_field" if m>=4 else
                "nearest_distinct_eigenvalue_not_field")
        })
    if assignments[0] is None:
        raise RuntimeError("Independent Whitney FEM failed to recover TM11 vicinity")
    records=[r for r in matched if r.get("family")=="TM"]
    field_modes=len(records)
    return {
        "comparisons":comparisons,
        "summary":{
            "efg_total_modes":len(efg_values),
            "edge_total_tm_modes":len(edge_native),
            "edge_tm_near_zero_cutoff_native":positive_cutoff_edge,
            "edge_tm_near_zero_cutoff_scaled":scaled_lambda(positive_cutoff_edge),
            "edge_positive_below_analytic_TM11":below_analytic,
            "edge_positive_below_half_TM11":deeply_below,
            "efg_positive_below_TM11_original":original_count,
            "edge_matched_csv_tm_modes":field_modes,
            "edge_matched_csv_min_rho":min((float(r["rho_abs"]) for r in records),default=None),
        },
        "native_tm":edge_rows, "native_matched":records, "timing":timing
    }

def export_csv(path,rows):
    with Path(path).open("w",newline="",encoding="utf-8") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)

def execute(edge_binary,out,baseline,timeout,edge_sha):
    if len(edge_sha)!=40 or not all(ch in "0123456789abcdef" for ch in edge_sha.lower()):
        raise ValueError("Full 40-character FEM edge source SHA is required")
    if edge_sha!=PINNED_FEM_SHA:
        raise ValueError("M6e protocol pins FEM edge at "+PINNED_FEM_SHA)
    out.mkdir(parents=True,exist_ok=True)
    summary,comparisons=[],[]
    for label,nx,ny,indices,expected in EFG_CASES:
        case=out/label
        edge_dir,command=run_edge(edge_binary,nx,ny,case/"edge_native",timeout)
        provenance=compare_case(edge_dir,baseline/label,indices,expected)
        entry=dict(case=label,edge_mesh_nx=nx,edge_mesh_ny=ny,
                   **provenance["summary"])
        summary.append(entry)
        comparisons.extend(dict(case=label,**r) for r in provenance["comparisons"])
        export_csv(case/"comparison.csv",provenance["comparisons"])
        info=dict(reference_repository=FEM_REPO,reference_commit=edge_sha,
                  edge_command=command,edge_native_rectangle=list(EDGE_GEOMETRY),
                  target_EFG_rectangle=list(EFG_GEOMETRY),
                  lambda_native_to_target_multiplier=GEOMETRIC_LAMBDA_MULTIPLIER,
                  edge_case_output=str(edge_dir.relative_to(out)),
                  comparison_category="TM_transverse_H_only",
                  note="Eigenvalue-only assignment; vector fields NOT independently cross-correlated.",
                  edge_source_binary=str(edge_binary),
                  efg_original_case=str((baseline/label)),
                  edge_input_sha256={str(p.relative_to(case)):sha256(p) for p in
                      sorted(edge_dir.rglob("*.csv"))},
                  efg_eigenvalues_sha256=sha256(baseline/label/"eigenvalues.csv"))
        (case/"provenance.json").write_text(
            json.dumps(info,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
        print(label,"EFG low-positive",expected,
              "edge deeply below TM11",entry["edge_positive_below_half_TM11"],
              "edge below analytical TM11",entry["edge_positive_below_analytic_TM11"],
              flush=True)
    export_csv(out/"summary.csv",summary)
    export_csv(out/"comparison.csv",comparisons)
    lines=["# M6e — EFG cartesiano versus FEM de aresta Whitney",
           "",
           "Referência FEM independente: "+FEM_REPO+" @ "+edge_sha,
           "Ramo usado: **TM com campo transversal H**, jamais TE.",
           "FEM edge no retângulo 1×0,5; EFG no retângulo 2×1.",
           "**Fator de conversão λ_edge para λ_EFG: 1/4.**",
           "Os meshes FEM utilizam células 8×4, 12×6 e 16×8;",
           "a nuvem EFG utiliza nós 9×7, 13×11 e 17×15.",
           "",
           "| Nuvem EFG | Células edge | EFG positivos abaixo TM11 | Edge positivos abaixo TM11 analítico | Edge abaixo de 0,5×TM11 |",
           "|---|---|---:|---:|---:|",
           *[f"| {s['case']} | {s['edge_mesh_nx']}×{s['edge_mesh_ny']} | "
             f"{s['efg_positive_below_TM11_original']} | "
             f"{s['edge_positive_below_analytic_TM11']} | "
             f"{s['edge_positive_below_half_TM11']} |" for s in summary],
           "",
           "Abaixo da referência analítica pode incluir a aproximação",
           "numérica do próprio TM11, portanto NÃO equivale automaticamente",
           "a espúrios. O cutoff numérico do FEM é publicado no summary.csv.",
           "",
           "Os pares TM22/TM41 são tratados como duas raízes distintas",
           "próximas de 5π²; **não** há validação de ângulos principais",
           "ou comparação direta dos campos H entre discretizações.",
           "O CSV de modos nativo do FEM pode trazer matching próprio,",
           "porém o protocolo comparativo usa os autovalores não filtrados.",
           "",
           "| Caso | Modo | Analítico | EFG λ | Edge λ reescalado | Erro EFG % | Erro edge % |",
           "|---|---|---:|---:|---:|---:|---:|",
           *[("| {case} | {mode} | {analytic_lambda_on_2x1:.7f} | "
             "{efg_lambda:.7f} | {edge} | {eerr:.5f} | {werr} |").format(
               **r,edge="n/a" if r["edge_scaled_lambda_on_2x1"]=="" else
               f"{r['edge_scaled_lambda_on_2x1']:.7f}",
               eerr=100*r["efg_relative_error"],
               werr="n/a" if r["edge_relative_error"]=="" else
               f"{100*r['edge_relative_error']:.5f}")
             for r in comparisons],
           "",
           "Dados brutos completos do FEM edge, matrizes S/T,",
           "autovalores/autovetores TE/TM, CSV de campos e logs,",
           "são preservados nas subpastas edge_native.",
           "Os baselines M6a-M6d não são alterados nem ocultados.",
           "",
           "## Limitações",
           "",
           "Este estágio é **comparação de espectros** e geometrias normalizadas.",
           "É indispensável uma campanha posterior com reconstrução",
           "vetorial de campos no mesmo grid físico, análise do dubleto",
           "por ângulos principais, divergência, fluxo normal e",
           "erro de orientação, antes de concluir conformidade física.",
           ""]
    (out/"REPORT.md").write_text("\n".join(lines),encoding="utf-8")
    print("\n".join(lines),flush=True)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--edge-binary",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--efg-baseline",type=Path,default=Path(
        "docs/results/independent_validation/2026-10-21_dc1d98b/efg"))
    parser.add_argument("--edge-source-sha",default=PINNED_FEM_SHA)
    parser.add_argument("--timeout",type=int,default=600)
    args=parser.parse_args()
    if args.timeout<1:parser.error("--timeout must be positive")
    binary=args.edge_binary.resolve()
    if not binary.is_file():parser.error("Edge FEM binary not found: "+str(binary))
    baseline=args.efg_baseline.resolve()
    if not baseline.is_dir():parser.error("EFG M6a baseline not found: "+str(baseline))
    execute(binary,args.out.resolve(),baseline,args.timeout,args.edge_source_sha)

if __name__=="__main__":
    main()
