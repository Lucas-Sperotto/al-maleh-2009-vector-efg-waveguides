#!/usr/bin/env python3
"""Independently check an affine/quartic M5 export with midpoint quadrature.

Run with PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1. NumPy/SciPy are required.
The input directory is read-only. JSON goes to stdout, or --out (new file).
No repository C++ solver, MLS derivatives or campaign helper is imported.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh


def read_rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def check_export(folder):
    metadata = dict(line.split("=", 1) for line in
                    (folder / "run_metadata.txt").read_text().splitlines() if "=" in line)
    if (metadata["mls_degree"], metadata["mls_weight"], metadata["quad_order"]) != (
            "1", "quartic_spline", "1"):
        raise ValueError("This independent check requires degree 1, quartic weight, Gauss 1")
    nx, ny = int(metadata["nx"]), int(metadata["ny"])
    n, sub = nx * ny, int(metadata["subdivision"])
    k = np.zeros((2 * n, 2 * n))
    m = np.zeros_like(k)
    for row in read_rows(folder / "matrix_elements.csv"):
        i, j = int(row["row"]), int(row["col"])
        k[i, j], m[i, j] = float(row["K"]), float(row["M"])
    saved_rows = read_rows(folder / "eigenvalues.csv")
    saved = np.array([float(row["lambda"]) for row in saved_rows])
    v = np.zeros((2 * n, len(saved)))
    for row in read_rows(folder / "eigenvectors.csv"):
        v[int(row["dof"]), int(row["mode_index"])] = float(row["coefficient"])
    mu, u = eigh(m)
    cutoff = max(1e-12, float(metadata["mass_rank_relative_cutoff"]) * mu[-1])
    positive = mu > cutoff
    z = u[:, ~positive]
    q = u[:, positive] / np.sqrt(mu[positive])
    spectrum, _ = eigh(q.T @ k @ q)
    residue = k @ v - (m @ v) * saved[np.newaxis, :]
    residuals = np.linalg.norm(residue, axis=0) / (
        (np.linalg.norm(k) + abs(saved) * np.linalg.norm(m)) * np.linalg.norm(v, axis=0))
    report = dict(
        rank=int(positive.sum()), nullity=int((~positive).sum()),
        mass_eigen_max=float(mu[-1]), largest_removed=float(mu[~positive][-1]),
        smallest_retained=float(mu[positive][0]),
        KZ_relative=float(np.linalg.norm(k @ z) / np.linalg.norm(k)),
        max_full_residual=float(max(residuals)),
        M_orthogonality_error=float(np.max(np.abs(v.T @ m @ v - np.eye(v.shape[1])))),
        scipy_max_spectrum_difference=float(np.max(np.abs(saved - spectrum))), fd=[])

    # Independent affine weighted least squares and quartic radial polynomial.
    nodes = np.array([(x, y) for y in np.linspace(0, 1, ny) for x in np.linspace(0, 2, nx)])
    radius = float(metadata["support_factor"]) * max(2 / (nx - 1), 1 / (ny - 1))

    def phi(point):
        d = (nodes - point) / radius
        r = np.linalg.norm(d, axis=1)
        w = np.maximum(0, 1 - 6 * r * r + 8 * r**3 - 3 * r**4)
        w[r >= 1] = 0
        p = np.column_stack([np.ones(n), d])
        a = p.T @ (w[:, None] * p)
        return np.linalg.solve(a, np.array([1., 0, 0])) @ (p.T * w)

    # The exported driver uses a 2x1 rectangle and positive cell-midpoint weights.
    cx, cy = (nx - 1) * sub, (ny - 1) * sub
    points = np.array([(2 * (i + .5) / cx, (j + .5) / cy)
                       for j in range(cy) for i in range(cx)])
    weights = np.full(len(points), 2 / (cx * cy))
    expected = read_rows(folder / "tm_matches.csv")
    zero_cutoff = max(1e-9, 1e-9 * max(saved))
    for step in [.001, .0005, .00025]:
        laplacian = []
        for x in points:
            center, lap = phi(x), np.zeros(n)
            for direction in np.eye(2):
                lap += (-phi(x + 2 * step * direction) + 16 * phi(x + step * direction)
                        - 30 * center + 16 * phi(x - step * direction)
                        - phi(x - 2 * step * direction)) / (12 * step * step)
            laplacian.append(lap)
        laplacian = np.array(laplacian)
        curl = laplacian @ v[n:]
        result = dict(step=step, K_relative_error=float(np.linalg.norm(
            laplacian.T @ (weights[:, None] * laplacian) - k[n:, n:]) /
            np.linalg.norm(k[n:, n:])), selections=[])
        for mm in [1, 2]:
            analytic = np.sin(mm * np.pi * points[:, 0] / 2) * np.sin(np.pi * points[:, 1])
            cross = curl.T @ (weights * analytic)
            energy = (weights[:, None] * curl**2).sum(axis=0)
            rho = abs(cross) / np.sqrt(np.maximum(energy, 1e-300) * (weights * analytic**2).sum())
            rho[saved <= zero_cutoff] = 0
            best = int(np.argmax(rho))
            distances = abs(saved - ((mm / 2)**2 + 1) * np.pi**2)
            distances[saved <= zero_cutoff] = np.inf
            nearest = int(np.argmin(distances))
            result["selections"].append(dict(
                mode=f"TM{mm}1", best_curl_index=best, rho=float(rho[best]),
                best_lambda=float(saved[best]), nearest_index=nearest,
                nearest_lambda=float(saved[nearest])))
            source = expected[mm - 1]
            if best != int(source["best_curl_index"]) or nearest != int(source["nearest_index"]):
                raise AssertionError("Independent curl/nearest mode selection differs")
            if abs(rho[best] - float(source["best_curl_rho"])) > 1e-7:
                raise AssertionError("Independent curl correlation differs by more than 1e-7")
        if result["K_relative_error"] > 1e-7:
            raise AssertionError("Independent finite-difference stiffness differs by more than 1e-7")
        report["fd"].append(result)
    if report["nullity"] != 4 or report["rank"] != 2 * n - 4:
        raise AssertionError("Affine common-kernel count differs for this experiment")
    for key, tolerance in [("KZ_relative", 1e-12), ("max_full_residual", 1e-12),
                           ("M_orthogonality_error", 1e-10),
                           ("scipy_max_spectrum_difference", 1e-9)]:
        if not np.isfinite(report[key]) or report[key] > tolerance:
            raise AssertionError(f"{key} exceeds {tolerance}")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("export_dir", type=Path)
    parser.add_argument("--out", type=Path, help="write a new JSON file without overwriting")
    args = parser.parse_args()
    result = json.dumps(check_export(args.export_dir), indent=2)
    if args.out:
        with args.out.open("x") as stream:
            stream.write(result + "\n")
    else:
        print(result)


if __name__ == "__main__":
    main()
