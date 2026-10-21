#!/usr/bin/env python3
"""M4 EFG benchmark campaign. Python stdlib only; no physical-mode claims.

Executes C++ tm_validation independently for each configuration. Distinguishes
numerical/conditioning failures from completed runs and retains run logs.
"""
import argparse
import csv
import os
from pathlib import Path
import platform
import subprocess
import sys

CASES = (
    ("refine_5x5", 5, 5, 1, 3.3, 3),
    ("refine_7x7", 7, 7, 1, 3.3, 3),
    ("refine_9x7", 9, 7, 1, 3.3, 3),
    ("radius_7x7_narrow", 7, 7, 1, 2.8, 3),
    ("radius_7x7_wide", 7, 7, 1, 3.8, 3),
    ("quadratic_7x7", 7, 7, 2, 3.8, 3),
    ("gauss2_7x7", 7, 7, 1, 3.3, 2),
    ("gauss4_7x7", 7, 7, 1, 3.3, 4),
)
MATCH_FIELDS = (
    "m", "n", "lambda_exact", "kc_exact",
    "nearest_index", "nearest_lambda", "nearest_relative_error",
    "nearest_curl_overlap", "nearest_cluster_overlap",
    "nearest_cluster_size", "strongest_index", "strongest_lambda",
    "strongest_overlap", "strongest_relative_error",
    "divergence_ratio", "normal_boundary_ratio",
    "algebraic_residual", "status",
)

def read_metadata(path):
    lines = path.read_text(encoding="utf-8").splitlines()
    return dict(line.split("=", 1) for line in lines if "=" in line)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path("out/m4"))
    parser.add_argument("--timeout", type=int, default=240)
    parser.add_argument("--profile", choices=("ci", "smoke"), default="ci")
    args = parser.parse_args()
    if args.timeout <= 0:
        parser.error("timeout must be positive")
    binary = args.binary.resolve()
    if not binary.is_file():
        parser.error("compiled executable not found: " + str(binary))
    output = args.out.resolve()
    output.mkdir(parents=True, exist_ok=True)
    cases = CASES if args.profile == "ci" else CASES[:1]
    run_rows = []
    match_rows = []
    for label, nx, ny, degree, radius, gauss in cases:
        destination = output / label
        destination.mkdir(exist_ok=True)
        cmd = [str(binary), str(nx), str(ny), str(degree),
               str(radius), str(gauss), str(destination)]
        print("CASE", label, "command:", " ".join(cmd), flush=True)
        status = "success"
        message = ""
        try:
            task = subprocess.run(cmd, timeout=args.timeout,
                                  text=True, capture_output=True, check=False)
            (destination / "stdout.txt").write_text(task.stdout, encoding="utf-8")
            (destination / "stderr.txt").write_text(task.stderr, encoding="utf-8")
            if task.returncode != 0:
                status = "solver_failed"
                message = task.stderr.strip().replace("\n", "; ")[:500]
            elif not (destination / "tm_matches.csv").exists():
                status = "missing_artifacts"
                message = "tm_matches.csv was not created"
            else:
                with (destination / "tm_matches.csv").open(
                    encoding="utf-8", newline=""
                ) as f:
                    for row in csv.DictReader(f):
                        if set(MATCH_FIELDS) - set(row):
                            raise ValueError("Incomplete matching table")
                        match_rows.append({"case": label, "nx": nx, "ny": ny,
                                           "degree": degree, "support": radius,
                                           "gauss": gauss, **row})
        except subprocess.TimeoutExpired:
            status, message = "timeout", "Execution time limit exceeded"
        except (OSError, ValueError, csv.Error) as exc:
            status, message = "report_failed", str(exc)[:500]
        metadata = {}
        if (destination / "run_metadata.txt").exists():
            metadata = read_metadata(destination / "run_metadata.txt")
        run_rows.append({"case": label, "nx": nx, "ny": ny,
                         "degree": degree, "support": radius, "gauss": gauss,
                         "status": status, "message": message,
                         "numerical_kernel": metadata.get("numerical_kernel", ""),
                         "mass_condition": metadata.get("mass_condition", ""),
                         "max_mls_condition": metadata.get("mls_max_condition", "")})
        print("RESULT", label, status, message, flush=True)
    with (output / "runs.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(run_rows[0]))
        writer.writeheader()
        writer.writerows(run_rows)
    with (output / "tm_matches_all.csv").open(
        "w", newline="", encoding="utf-8"
    ) as f:
        writer = csv.DictWriter(f, fieldnames=[
            "case", "nx", "ny", "degree", "support", "gauss", *MATCH_FIELDS
        ])
        writer.writeheader()
        writer.writerows(match_rows)

    def summary_line(case):
        matches = [r for r in match_rows if r["case"] == case and
                   r["m"] == "1" and r["n"] == "1"]
        if not matches:
            return "| " + case + " | falha | — | — | — |"
        x = matches[0]
        return ("| {} | {} | {:.5f} | {:.5f} | {} |".format(
            case, x["nearest_index"],
            float(x["nearest_relative_error"]),
            float(x["nearest_cluster_overlap"]), x["status"]))

    report = [
        "# M4 — Campanha espectral EFG cartesiano (diagnóstico exploratório)",
        "",
        "**Não constitui certificação de reprodução da dissertação nem "
        "identificação definitiva de modos TM físicos.**",
        "",
        "A geometria é um retângulo 2×1. O autovalor analítico TM11 é "
        "1.25π² ≈ 12.3370055. Cada caso monta K/M novamente e "
        "resolve o autoproblema.",
        "",
        "| Caso | Índice próximo TM11 | Erro relativo λ11 | Sobreposição curl TM11 | Diagnóstico |",
        "|---|---:|---:|---:|---|",
        *[summary_line(row["case"]) for row in run_rows],
        "",
        "## Correspondência dos modos TM21 e TM12 (refinamentos)",
        "",
        "| Caso | Referência | Erro relativo lambda | Projeção de curl | Divergência | Fluxo normal |",
        "|---|---|---:|---:|---:|---:|",
        *[
            "| {} | TM{}{} | {:.8f} | {:.8f} | {:.6f} | {:.6f} |".format(
                case, m, n, float(row["nearest_relative_error"]),
                float(row["nearest_cluster_overlap"]),
                float(row["divergence_ratio"]),
                float(row["normal_boundary_ratio"]))
            for case in (r["case"] for r in run_rows[:3])
            for m, n in (("2", "1"), ("1", "2"))
            for row in match_rows
            if row["case"] == case and row["m"] == m and row["n"] == n
        ],
        "",
        "## Diagnóstico do núcleo e massa",
        "",
        "| Caso | Estado | Núcleo numérico | Condição da massa |",
        "|---|---|---:|---:|",
        *[
            "| {} | {} | {} | {} |".format(
                r["case"], r["status"],
                r["numerical_kernel"], r["mass_condition"])
            for r in run_rows
        ],
        "",
        "Interpretação: a sobreposição mede a correlação normalizada entre "
        "curl(H) e sin(mπx/a) sin(nπy/b), com tratamento de clusters "
        "espectrais degenerados. Erro pequeno no autovalor **sozinho** não "
        "é evidência de um modo físico. Falhas de condicionamento ou "
        "de posto não são ocultadas; consulte runs.csv e os stderr.txt.",
        "",
        "A convergência deve ser estudada em sequência com diagnóstico "
        "de divergência, fluxo normal, estabilidade de suporte e "
        "convergência de autofunções.",
        "",
        "## Ambiente",
        "",
        "- Python: " + sys.version.split()[0],
        "- Plataforma: " + platform.platform(),
        "- Commit CI: " + os.environ.get("GITHUB_SHA", "não registrado"),
        "- Execução: " + args.profile,
        "",
        "Arquivos-fonte: tm_matches_all.csv, runs.csv e pasta por caso.",
        "",
    ]
    (output / "REPORT.md").write_text("\n".join(report), encoding="utf-8")
    succeeded = sum(r["status"] == "success" for r in run_rows)
    required = 3 if args.profile == "ci" else 1
    success_refinements = sum(r["status"] == "success" for r in run_rows[:3])
    print("\n".join(report), flush=True)
    print("M4:", succeeded, "of", len(run_rows),
          "completed;", success_refinements, "refinement cases succeeded")
    # The sweep records sensitivity failures without hiding them. The
    # three baseline refinements must complete for CI to pass.
    if succeeded < required or (
        args.profile == "ci" and success_refinements < 3
    ):
        return 1
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
