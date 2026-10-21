#!/usr/bin/env python3
"""Run independent M6b constraint audits; never rewrites historical M6a data."""
import argparse
import csv
from pathlib import Path
import subprocess

CASES = (("9x7", 9, 7, 70, 71, 49), ("13x11", 13, 11, 127, 128, 94),
         ("17x15", 17, 15, 192, 193, 144))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()
    if args.timeout <= 0:
        parser.error("Timeout must be positive")
    binary = args.binary.resolve()
    if not binary.exists():
        parser.error("Missing compiled M6b executable: " + str(binary))
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for name, nx, ny, i, j, expected_count in CASES:
        destination = out / name
        command = [str(binary), str(nx), str(ny), str(i), str(j),
                   str(destination)]
        print("M6b:", " ".join(command), flush=True)
        result = subprocess.run(command, text=True, capture_output=True,
                                timeout=args.timeout, check=False)
        destination.mkdir(parents=True, exist_ok=True)
        (destination / "stdout.txt").write_text(result.stdout, encoding="utf-8")
        (destination / "stderr.txt").write_text(result.stderr, encoding="utf-8")
        if result.returncode != 0:
            raise RuntimeError(name + " solver failed: " + result.stderr)
        with (destination / "subspace_constraints.csv").open(newline="",encoding="utf-8") as source:
            data=list(csv.DictReader(source))
        if len(data)!=1:
            raise RuntimeError("Expected exactly one subspace record")
        row=data[0]
        with (destination / "positive_below_tm11.csv").open(newline="",encoding="utf-8") as source:
            low=list(csv.DictReader(source))
        if len(low)!=expected_count:
            raise RuntimeError("Published M6a baseline counts changed for "+name)
        for key in ("divergence_ratio_min","divergence_ratio_max",
                    "normal_ratio_min","normal_ratio_max"):
            val=float(row[key])
            if val<0 or not (val<1e6):
                raise RuntimeError("Invalid M6b metric "+key)
        if float(row["divergence_ratio_min"])>float(row["divergence_ratio_max"]):
            raise RuntimeError("Divergence range reversed")
        if float(row["normal_ratio_min"])>float(row["normal_ratio_max"]):
            raise RuntimeError("Boundary range reversed")
        if name=="17x15" and float(row["divergence_ratio_max"]) < 0.415:
            raise RuntimeError("M6b invariant maximum fails historical lower bound")
        rows.append(dict(case=name,low_positive_count=len(low),**row))
        print(name,"invariant divergence",row["divergence_ratio_min"],
              row["divergence_ratio_max"],"low-positive",len(low),flush=True)
    with (out / "summary.csv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    md=["# M6b — Subespaço TM22/TM41: restrições invariantes",
        "",
        "**Diagnóstico independente sem alterar K/M, autovetores ou dados históricos de M6a.**",
        "As referências são lambda(TM22)=lambda(TM41)=5π²; valores mínimo/máximo "
        "resultam de problemas generalizados das formas quadráticas D/B projetadas "
        "no mesmo espaço de duas dimensões.",
        "",
        "| Nuvem | Positivos abaixo de TM11 | Divergência mín | Divergência máx | Normal mín | Normal máx |",
        "|---|---:|---:|---:|---:|---:|",
        *["| {} | {} | {:.9f} | {:.9f} | {:.9f} | {:.9f} |".format(
            row["case"],row["low_positive_count"],
            *[float(row[key]) for key in (
                "divergence_ratio_min","divergence_ratio_max",
                "normal_ratio_min","normal_ratio_max")])
          for row in rows],
        "",
        "O valor máximo é invariante a qualquer troca de base não singular do "
        "subespaço. As direções de máximo são coeficientes nas colunas "
        "originais do Eigen; os autovetores de entrada não foram alterados.",
        "Uma reprovação com máximo > 0,25 não pode ser eliminada por rotação "
        "dos autovetores. A classificação continua exploratória e nenhum "
        "candidato positivo adicional foi automaticamente identificado como físico.",
        "",
        "Arquivos de cada caso: subspace_constraints.csv, positive_below_tm11.csv, "
        "run_metadata.txt, stdout.txt e stderr.txt.",
        ""]
    (out/"REPORT.md").write_text("\n".join(md),encoding="utf-8")
    print("\n".join(md),flush=True)


if __name__ == "__main__":
    main()
