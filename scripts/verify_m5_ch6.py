#!/usr/bin/env python3
"""Independent stdlib archive verification of original Al Maleh Ch6 M5 runs."""
import argparse
import csv
import hashlib
import math
from pathlib import Path
from m5_ch6_campaign import CASES,read_reference

def rows(path):
    with path.open(encoding="utf-8",newline="") as stream:
        return list(csv.DictReader(stream))

def kv(path):
    return dict(line.split("=",1) for line in path.read_text(encoding="utf-8").splitlines()
                if "=" in line)

def verify(root):
    manifest=rows(root/"manifest_sha256.csv")
    filenames=set()
    for entry in manifest:
        rel=entry["file"]
        if rel in filenames or rel.startswith("/") or ".." in Path(rel).parts:
            raise ValueError("Repeated/unsafe hash manifest path")
        filenames.add(rel)
        file=root/rel
        if not file.is_file() or hashlib.sha256(file.read_bytes()).hexdigest()!=entry["sha256"]:
            raise ValueError("SHA256 mismatch: "+rel)
    compared=rows(root/"comparison.csv")
    reference=read_reference()
    if len(compared)!=24:
        raise ValueError("M5 campaign comparison incomplete")
    for a,b in zip(reference,compared):
        if any(a[k]!=b[k] for k in ("mode","nodes_y","nodes_x","gauss_points_per_axis")):
            raise ValueError("M5 source/observed rows misaligned")
        if not math.isclose(float(a["author_signed_error_percent"]),
                            float(b["author_signed_error_percent"]),abs_tol=1e-12):
            raise ValueError("Changed thesis reference number")
        for k in ("current_nearest_signed_error_percent",
                  "current_best_curl_signed_error_percent",
                  "current_best_curl_rho"):
            if not math.isfinite(float(b[k])):
                raise ValueError("Nonfinite M5 modal comparison")
        if not 0<=float(b["current_best_curl_rho"])<=1+1e-8:
            raise ValueError("Unphysical curl correlation")
    for (nx,ny),orders in CASES.items():
        for order in orders:
            path=root/f"{ny}x{nx}_gauss{order}"
            meta=kv(path/"run_metadata.txt")
            dofs=2*nx*ny
            rank=int(meta["retained_rank"])
            mass_null=int(meta["removed_mass_nullity"])
            zero=int(meta["near_zero_eigenvalues"])
            nquad=int(meta["integration_points"])
            if rank+mass_null!=dofs or not 0<rank<dofs or zero<0 or zero>rank:
                raise ValueError("Invalid M5 mass/nullspace counts")
            if nquad!=4*(nx-1)*(ny-1)*order*order:
                raise ValueError("Incorrect four-background-cell Ch6 integration")
            if mass_null!=4:
                raise ValueError("Changed verified affine MLS mass nullity")
            if (nx,ny)==(5,3) and zero!=14:
                raise ValueError("Changed 3x5 Ch6 near-zero curl count")
            if float(meta["M_orthogonality_error"])>1e-5:
                raise ValueError("Mass-orthogonality failure")
            eig=rows(path/"eigenvalues.csv")
            if len(eig)!=rank:
                raise ValueError("Incomplete reduced eigenvalue spectrum")
            values=[float(r["lambda"]) for r in eig]
            if any(not math.isfinite(x) for x in values):
                raise ValueError("Nonfinite EVP spectrum")
            if any(values[i]>values[i+1]+1e-9 for i in range(rank-1)):
                raise ValueError("M5 eigenspectrum not sorted")
            if any(float(r["reduced_residual"])>1e-7 for r in eig):
                raise ValueError("Projected eigenpair fails algebraic residual test")
            vectors=rows(path/"eigenvectors.csv")
            if len(vectors)!=rank*dofs:
                raise ValueError("Missing M5 coefficient columns")
            if {int(r["mode_index"]) for r in vectors}!=set(range(rank)):
                raise ValueError("Missing eigenvector")
            if len(rows(path/"tm_matches.csv"))!=2:
                raise ValueError("Missing TM11 or TM21")
            for name in ("tm11_vector_field","tm21_vector_field",
                         "tm11_irrot_field","tm11_rot_field"):
                if len(rows(path/(name+".csv")))!=25*13:
                    raise ValueError("Incomplete M5 vector field: "+name)
                if not (path/(name+".svg")).is_file():
                    raise ValueError("Missing M5 SVG plot")
            if (nx,ny)==(9,5):
                for name in ("irrotational_shape_central",
                             "rotational_shape_central"):
                    if len(rows(path/(name+".csv")))!=25*13:
                        raise ValueError("Missing original Ch6 central shape function")
                    if not (path/(name+".svg")).is_file():
                        raise ValueError("Missing original Ch6 shape diagram")
            print(f"PASS M5 {ny}x{nx} gauss={order} rank={rank} mass_null={mass_null} "
                  f"curl_zero={zero} nquad={nquad}")
    print(f"PASS M5 full manifest ({len(filenames)} SHA256 hashes), "
          "12 configurations, 24 source comparisons")

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--campaign",type=Path,required=True)
    a=p.parse_args()
    verify(a.campaign.resolve())

if __name__=="__main__":
    main()
