#!/usr/bin/env python3
"""Independently reopen every M6c Matrix Market export and check constrained EVP."""
import argparse
import csv
from pathlib import Path
import numpy as np
from scipy.io import mmread


def csv_rows(path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def verify(path):
    path=Path(path)
    S=np.asarray(mmread(path/"stiffness_matrix.mtx"))
    M=np.asarray(mmread(path/"mass_matrix.mtx"))
    C=np.asarray(mmread(path/"boundary_collocation.mtx"))
    Q=np.asarray(mmread(path/"boundary_nullspace.mtx"))
    V=np.asarray(mmread(path/"eigenvectors.mtx"))
    values=csv_rows(path/"eigenvalues.csv")
    dof=csv_rows(path/"dof_map.csv")
    n=len(dof)
    k=len(values)
    if any(A.ndim!=2 for A in (S,M,C,Q,V)):
        raise ValueError("Matrix Market source must be 2D")
    if S.shape!=(n,n) or M.shape!=(n,n) or Q.shape!=(n,k) or V.shape!=(n,k):
        raise ValueError("M6c complete spectral matrix dimensions inconsistent")
    if C.shape[1]!=n or C.shape[0]!=sum(r["is_boundary_collocation_point"]=="true" for r in dof):
        raise ValueError("M6c boundary row count or DOF map inconsistent")
    if any(int(row["index"])!=j for j,row in enumerate(values)):
        raise ValueError("M6c eigenvalue index mismatch")
    lam=np.array([float(row["lambda"]) for row in values])
    if (not all(np.all(np.isfinite(A)) for A in (S,M,C,Q,V))
            or not np.all(np.isfinite(lam))
            or not np.all(lam>0)
            or not np.all(np.diff(lam)>=0)):
        raise ValueError("Nonfinite/nonpositive/unordered scalar EFG spectrum")
    if not np.allclose(S,S.T,atol=1e-10,rtol=0) or not np.allclose(M,M.T,atol=1e-10,rtol=0):
        raise ValueError("Scalar MLS stiffness/mass not symmetric")
    if not np.allclose(Q.T@Q,np.eye(k),atol=1e-7,rtol=0):
        raise ValueError("SVD boundary nullspace not Euclidean orthonormal")
    if not np.allclose(C@Q,0,atol=1e-8,rtol=0):
        raise ValueError("SVD nullspace violates PEC collocation")
    if not np.allclose(C@V,0,atol=1e-7,rtol=0):
        raise ValueError("Complete eigencoefficients violate PEC collocation")
    gram=V.T@M@V
    gram_err=float(np.max(np.abs(gram-np.eye(k))))
    if gram_err>1e-6:
        raise ValueError(f"Scalar M orthogonality error: {gram_err}")
    residual=Q.T@(S@V-(M@V)*lam)
    denominator=max(1e-300,(np.linalg.norm(Q.T@S@Q)+np.max(np.abs(lam))*np.linalg.norm(Q.T@M@Q))*np.max(np.linalg.norm(Q.T@V,axis=0)))
    relative=float(np.linalg.norm(residual)/denominator)
    if relative>1e-7:
        raise ValueError(f"Projected scalar EVP residual: {relative}")
    published=np.array([float(r["projected_relative_residual"]) for r in values])
    if not np.all(np.isfinite(published)) or published.max()>1e-7:
        raise ValueError("Published projected eigenpair residuals exceed bounds")
    return dict(case=path.name,scalar_dofs=n,reduced_dofs=k,
                orthogonality_max=gram_err,projected_residual=relative,
                collocation_max=float(np.max(np.abs(C@V))))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign",type=Path,required=True)
    args=parser.parse_args()
    for name in ("9x7","13x11","17x15"):
        row=verify(args.campaign/name)
        print("PASS",row,flush=True)


if __name__=="__main__":
    main()
