#!/usr/bin/env python3
"""Reproduce M6d Cartesian gradient-space diagnostic without filtering eigenpairs."""
import argparse
import csv
import json
from pathlib import Path
import statistics
import subprocess

CASES=(("9x7",9,7,49),("13x11",13,11,94),("17x15",17,15,144))
SENSITIVITY=(("sensitive_1e-6",1e-6),("sensitive_1e-10",1e-10))

def load_csv(path):
    with path.open(newline="",encoding="utf-8") as stream:
        return list(csv.DictReader(stream))

def metadata(path):
    return dict(line.split("=",1) for line in
                path.read_text(encoding="utf-8").splitlines() if "=" in line)

def run(binary,out,nx,ny,rank_cutoff,timeout):
    out.mkdir(parents=True,exist_ok=True)
    cmd=[str(binary),str(nx),str(ny),str(out),str(rank_cutoff)]
    execution=subprocess.run(cmd,text=True,capture_output=True,
                             timeout=timeout,check=False)
    (out/"stdout.txt").write_text(execution.stdout,encoding="utf-8")
    (out/"stderr.txt").write_text(execution.stderr,encoding="utf-8")
    if execution.returncode:
        raise RuntimeError("M6d execution failed in "+str(out)+": "+
                           execution.stderr)
    return metadata(out/"metadata.txt")

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--timeout",type=int,default=600)
    args=parser.parse_args()
    binary=args.binary.resolve()
    if not binary.is_file() or args.timeout<1:
        parser.error("Binary does not exist or timeout invalid")
    out=args.out.resolve()
    out.mkdir(parents=True,exist_ok=True)
    summaries=[]
    for name,nx,ny,expected in CASES:
        target=out/name
        print("M6d",name,"default projection rank tolerance=1e-8",flush=True)
        m=run(binary,target,nx,ny,1e-8,args.timeout)
        modes=load_csv(target/"eigenmode_gradient.csv")
        selected=load_csv(target/"selected_modes.csv")
        low=[row for row in modes if row["original_class"]=="positive_below_TM11"]
        if len(modes)!=2*nx*ny or len(low)!=expected or int(m["positive_below_TM11"])!=expected:
            raise RuntimeError(name+": historical eigenmode count changed")
        if len(selected)!=6:
            raise RuntimeError("Missing selected modes "+name)
        fraction=[float(row["gradient_energy_fraction"]) for row in low]
        if not all(0<=f<=1 for f in fraction):
            raise RuntimeError("Gradient fraction out of bounds")
        summaries.append(dict(
            case=name,vector_dofs=2*nx*ny,
            low_positive=len(low),gradient_rank=int(m["projected_gradient_rank"]),
            gram_cutoff=float(m["gradient_rank_absolute_cutoff"]),
            low_grad_fraction_mean=statistics.mean(fraction),
            low_grad_fraction_median=statistics.median(fraction),
            low_grad_fraction_min=min(fraction),
            low_grad_fraction_max=max(fraction),
            low_grad_fraction_ge_0_9=sum(f>=0.9 for f in fraction),
            TM11_grad_fraction=float(selected[0]["gradient_energy_fraction"]),
            TM21_grad_fraction=float(selected[1]["gradient_energy_fraction"]),
            TM31_grad_fraction=float(selected[2]["gradient_energy_fraction"]),
            TM12_grad_fraction=float(selected[3]["gradient_energy_fraction"]),
            doublet_A_grad_fraction=float(selected[4]["gradient_energy_fraction"]),
            doublet_B_grad_fraction=float(selected[5]["gradient_energy_fraction"]),
            relative_gradient_fit_matrix_discrepancy=
                float(m["gradient_energy_matrix_relative_discrepancy"])
        ))
        print(name,"positive below TM11",expected,
              "projected gradient rank",summaries[-1]["gradient_rank"],
              "mean projected gradient",summaries[-1]["low_grad_fraction_mean"],
              flush=True)

    sensitivity=[]
    for name,threshold in SENSITIVITY:
        target=out/"rank_sensitivity"/name
        m=run(binary,target,17,15,threshold,args.timeout)
        modes=load_csv(target/"eigenmode_gradient.csv")
        low=[r for r in modes if r["original_class"]=="positive_below_TM11"]
        if len(low)!=144:
            raise RuntimeError("Rank sensitivity unexpectedly changed K/M spectrum")
        sensitivity.append(dict(case=name,relative_cutoff=threshold,
                                rank=int(m["projected_gradient_rank"]),
                                low_grad_mean=statistics.mean(
                                    float(r["gradient_energy_fraction"]) for r in low),
                                TM11_grad_fraction=float(
                                    modes[153]["gradient_energy_fraction"]),
                                doublet_A_grad_fraction=float(
                                    modes[192]["gradient_energy_fraction"]),
                                doublet_B_grad_fraction=float(
                                    modes[193]["gradient_energy_fraction"])))

    with (out/"summary.csv").open("w",newline="",encoding="utf-8") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(summaries[0]))
        writer.writeheader();writer.writerows(summaries)
    with (out/"rank_sensitivity.csv").open("w",newline="",encoding="utf-8") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(sensitivity[0]))
        writer.writeheader();writer.writerows(sensitivity)
    description=[
        "# M6d — Projeção L2 de gradientes escalares no EFG vetorial",
        "",
        "**Diagnóstico, não intervenção na solução.** Os autovalores e",
        "autovetores originais de K h = λ M h permanecem inalterados.",
        "A projeção não prova uma sequência discreta gradiente-rotacional",
        "compatível; a fração gradiente por si só não determina quais",
        "autovetores são físicos ou espúrios.",
        "",
        "| Nuvem | Abaixo de TM11 | Rank gradiente | Fração grad. média | Mediana | >=90% | TM11 | TM22/TM41 A | TM22/TM41 B |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
        *["| {case} | {low_positive} | {gradient_rank} | {low_grad_fraction_mean:.6f} | "
          "{low_grad_fraction_median:.6f} | {low_grad_fraction_ge_0_9} | "
          "{TM11_grad_fraction:.6f} | {doublet_A_grad_fraction:.6f} | "
          "{doublet_B_grad_fraction:.6f} |".format(**r) for r in summaries],
        "",
        "Rank definido por λ(Gram projetada)>cutoff*λmax, por padrão cutoff=1e-8.",
        "O estudo 17×15 com 1e-6 e 1e-10 está em rank_sensitivity.csv.",
        "A norma do erro de representação de gradientes escalares também",
        "foi salva; não se assumiu que gradientes MLS pertencem ao",
        "espaço vetorial MLS cartesiano.",
        "",
        "Arquivos por caso: eigenmode_gradient.csv (TODO o espectro),",
        "selected_modes.csv, gradient_fit_probes.csv, metadados e",
        "matrizes Matrix Market em precisão 17 algarismos:",
        "coupling H, escalar S, Gram projetada, base Q, massa vetorial",
        "M e autovetores V completos.",
        "",
        "**Não foi aplicado filtro, restrição solenoidal, penalização ou",
        "recalibração de limiar de autovalores.** A campanha histórica",
        "M6a e sua classificação permanecem separadas.",
        ""
    ]
    (out/"REPORT.md").write_text("\n".join(description),encoding="utf-8")
    print("\n".join(description),flush=True)

if __name__=="__main__":
    main()
