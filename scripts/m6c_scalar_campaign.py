#!/usr/bin/env python3
"""Run a distinct scalar MLS/EFG TM control, compare only grounded mode values."""
import argparse
import csv
import math
from pathlib import Path
import subprocess

CASES = (("9x7", 9, 7), ("13x11", 13, 11), ("17x15", 17, 15))
SIMPLE = ((1,1),(2,1),(3,1),(1,2))


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()
    binary=args.binary.resolve()
    if not binary.is_file() or args.timeout<=0:
        parser.error("Missing binary or invalid timeout")
    out=args.out.resolve()
    out.mkdir(parents=True,exist_ok=True)
    summary=[]
    for label,nx,ny in CASES:
        dest=out/label
        command=[str(binary),str(nx),str(ny),str(dest)]
        print("M6c:", " ".join(command),flush=True)
        job=subprocess.run(command,capture_output=True,text=True,
                           timeout=args.timeout,check=False)
        dest.mkdir(parents=True,exist_ok=True)
        (dest/"stdout.txt").write_text(job.stdout,encoding="utf-8")
        (dest/"stderr.txt").write_text(job.stderr,encoding="utf-8")
        if job.returncode:
            raise RuntimeError(label+": "+job.stderr)
        values=read_csv(dest/"eigenvalues.csv")
        matches=read_csv(dest/"tm_reference.csv")
        if len(values)<4 or not matches:
            raise RuntimeError("Incomplete scalar eigenpairs or matching table")
        if any(int(row["index"])!=i for i,row in enumerate(values)):
            raise RuntimeError("Missing/reordered scalar eigenpair indices")
        if not all(math.isfinite(float(row["lambda"])) and float(row["lambda"])>0
                   for row in values):
            raise RuntimeError("Nonpositive/nonfinite scalar spectrum")
        if any(float(row["lambda"])>float(values[i+1]["lambda"])
               for i,row in enumerate(values[:-1])):
            raise RuntimeError("Scalar eigenvalues are not ordered")
        if any(float(row["projected_relative_residual"])>1e-7 or
               float(row["pec_collocation_relative_residual"])>1e-7
               for row in values):
            raise RuntimeError("Scalar constrained problem failed algebraic test")
        for row in matches:
            if (int(row["m"]),int(row["n"])) not in SIMPLE:
                continue
            summary.append(dict(
                case=label,nodes=nx*ny,reduced_dofs=len(values),
                mode="TM"+row["m"]+row["n"],
                lambda_analytic=row["lambda_analytic"],
                lambda_scalar_efg=row["lambda_nearest"],
                relative_lambda_error=row["relative_lambda_error"],
                overlap_scalar_shape=row["scalar_overlap"],
                boundary_trace_ratio=row["boundary_trace_ratio"],
                normal_flux_ratio=row["normal_flux_ratio"],
                matching="exploratory_nearest_not_unique",
            ))
        print(label,"positive scalar roots",len(values),
              "lambda_first",values[0]["lambda"],flush=True)
    if len(summary)!=12:
        raise RuntimeError("Missing simple TM comparisons")
    with (out/"summary.csv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)
    md=["# M6c — Controle escalar MLS/EFG para TM (experimental)", "",
        "**Formulação nova e separada:** potencial longitudinal Ez aproximado "
        "por MLS; contorno PEC imposto por colocalização nos pontos de "
        "fronteira. Não é a formulação vetorial de Al Maleh e não corrige "
        "automaticamente o espectro do EFG cartesiano.",
        "",
        "A condição Ez=0 é exata nos pontos de colocalização, **não ao longo de "
        "toda a parede**. A identidade div(z×grad Ez)=0 vale para campos "
        "diferenciáveis; o fluxo normal remanescente revela o erro entre "
        "pontos de contorno.",
        "",
        "| Nuvem | DOFs escalares reduzidos | Modo | Erro em lambda (%) | Overlap escalar | Traço PEC | Fluxo normal |",
        "|---|---:|---|---:|---:|---:|---:|",
        *["| {case} | {reduced_dofs} | {mode} | {err:.6f} | {overlap:.8f} | {trace:.6f} | {normal:.6f} |".format(
            **row,err=100*float(row["relative_lambda_error"]),
            overlap=float(row["overlap_scalar_shape"]),
            trace=float(row["boundary_trace_ratio"]),
            normal=float(row["normal_flux_ratio"])) for row in summary],
        "",
        "Os índices dos autovalores são atribuídos aqui apenas por proximidade "
        "espectral exploratória. O par degenerado TM22/TM41 consta dos CSVs "
        "individuais, mas **não recebe certificação** sem análise conjunta "
        "do subespaço.",
        "",
        "Cada pasta tem os autovalores completos, coeficientes autovetoriais "
        "MLS escalares em Matrix Market, matrizes escalares S/M, matriz de "
        "colocalização C, base Q do núcleo de C e metadados.",
        "",
        "O resíduo algébrico só é nulo **no espaço admissível**, i.e., "
        "Qᵀ(Sv−lambda Mv). O resíduo sem projeção em geral não é zero, "
        "pois aparecem multiplicadores para as restrições.",
        ""]
    (out/"REPORT.md").write_text("\n".join(md),encoding="utf-8")
    print("\n".join(md),flush=True)


if __name__=="__main__":
    main()
