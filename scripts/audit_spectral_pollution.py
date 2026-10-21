#!/usr/bin/env python3
"""Evidence-only audit of TM-positive spectral pollution in archived M6a data.

Stdlib-only. Does not relabel eigenpairs, change numerical cutoffs, rewrite
scientific artifacts or conflate low curl-curl values with physical TM modes.
"""
import argparse
import csv
import json
import math
from pathlib import Path

TM11 = 1.25 * math.pi * math.pi
TM_DEGENERATE = 5.0 * math.pi * math.pi
BOUND = 0.25


def read_rows(path):
    with open(path, newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def audit_efg_case(directory):
    directory = Path(directory)
    metadata = json.loads((directory / "run_metadata.json").read_text(encoding="utf-8"))
    eigenvalues = read_rows(directory / "eigenvalues.csv")
    diagnostics = read_rows(directory / "mode_diagnostics.csv")
    if not eigenvalues or len(eigenvalues) != len(diagnostics):
        raise ValueError("Incomplete eigenpairs/diagnostics: " + str(directory))
    n = len(eigenvalues)
    if n != metadata["spectrum"]["dofs"]:
        raise ValueError("Exported eigenpair count disagrees with metadata")
    vals = []
    for i, (v, d) in enumerate(zip(eigenvalues, diagnostics)):
        if int(v["index"]) != i or int(d["index"]) != i:
            raise ValueError("Eigenvector index/diagnostic row mismatch")
        lam = float(v["lambda"])
        if not math.isfinite(lam):
            raise ValueError("Nonfinite eigenvalue")
        if abs(lam - float(d["lambda"])) > max(1e-11, 1e-12 * abs(lam)):
            raise ValueError("Eigenvalues and diagnostics disagree")
        vals.append(lam)
    if any(b < a for a, b in zip(vals, vals[1:])):
        raise ValueError("Spectrum is not ascending")
    cutoff = float(metadata["spectrum"]["zero_cutoff"])
    classified = sum(lam <= cutoff for lam in vals)
    if classified != metadata["spectrum"]["numerical_kernel_dimension"]:
        raise ValueError("Near-zero dimension disagrees with metadata")
    low = [i for i, lam in enumerate(vals) if cutoff < lam < TM11]
    ranges = [(0, 1e-4), (1e-4, 1e-2), (1e-2, 1), (1, TM11)]
    bands = {
        f"({lower:.6g},{upper:.6g})": sum(lower < vals[i] < upper for i in low)
        for lower, upper in ranges
    }
    bad_divergence = sum(
        float(diagnostics[i]["divergence_ratio"]) > BOUND
        for i in low
    )
    return {
        "case": directory.name, "dofs": n,
        "numerical_kernel": classified, "zero_cutoff": cutoff,
        "low_positive_below_tm11": len(low),
        "low_positive_with_divergence_above_025": bad_divergence,
        "low_positive_lambda_min": vals[low[0]] if low else None,
        "low_positive_lambda_max": vals[low[-1]] if low else None,
        "lambda_intervals": bands,
        "tm11_lambda": TM11,
    }


def audit_degenerate_pair(directory, indices):
    if len(indices) != 2 or len(set(indices)) != 2 or min(indices) < 0:
        raise ValueError('The selected degenerate subspace needs two distinct valid indices')
    directory = Path(directory)
    diagnostics = read_rows(directory / "mode_diagnostics.csv")
    if max(indices) >= len(diagnostics):
        raise ValueError("Requested degenerate indices not found")
    directions = []
    for i in indices:
        row = diagnostics[i]
        if int(row["index"]) != i:
            raise ValueError("Degenerate column/diagnostic index mismatch")
        norm = float(row["field_norm_l2"])
        div_norm = float(row["div_norm_l2"])
        if norm <= 0 or not all(map(math.isfinite, (norm, div_norm))):
            raise ValueError("Invalid degenerate field norms")
        directions.append({
            "index": i, "lambda": float(row["lambda"]),
            "divergence_ratio_at_analytic_lambda": div_norm / (math.sqrt(TM_DEGENERATE) * norm),
            "native_divergence_ratio": float(row["divergence_ratio"]),
            "native_normal_ratio": float(row["normal_ratio"]),
        })
    # For ANY chosen basis representing this exact two-dimensional span,
    # the supremum of the divergence Rayleigh quotient on the span is
    # at least the maximum of its values on these two vectors.
    # This is a mathematically justified lower bound, NOT the exact
    # rotation-invariant operator norm, which requires cross terms.
    bound = max(row["divergence_ratio_at_analytic_lambda"] for row in directions)
    return {
        "pair": "TM22/TM41",
        "reference_lambda": TM_DEGENERATE,
        "directions": directions,
        "lower_bound_on_basis_invariant_worst_case_divergence": bound,
        "threshold": BOUND,
        "subspace_cannot_pass_divergence_threshold": bound > BOUND,
        "exact_subspace_operator_norm_not_computed": True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", required=True, type=Path)
    parser.add_argument("--output", type=Path, help="Optional output path; never writes to source unless explicitly requested")
    args = parser.parse_args()
    root = args.campaign
    cases = [audit_efg_case(root / "efg" / key) for key in ("9x7", "13x11", "17x15")]
    cross = read_rows(root / "comparison.csv")
    pair_rows = [row for row in cross if row["mode"] in ("TM22", "TM41")]
    if len(pair_rows) != 2 or pair_rows[0]["efg_indices"] != pair_rows[1]["efg_indices"]:
        raise ValueError("Published degenerate pair assignments are inconsistent")
    indices = tuple(int(x) for x in pair_rows[0]["efg_indices"].split(";"))
    pair = audit_degenerate_pair(root / "efg" / "17x15", indices)
    results = {
        "purpose": "physical TM lower-spectrum and degenerate-subspace evidence, not a physical-mode certification",
        "case_audits": cases, "degenerate_pair": pair,
        "note": "Positive entries below TM11 are unconfirmed for the rectangular TM spectrum; direct numerical operator details and near-zero mode identification require additional analysis.",
    }
    encoded = json.dumps(results, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    print(encoded, end="")
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8")


if __name__ == "__main__":
    main()
