#!/usr/bin/env python3
"""Reopen all M6d numerical exports and independently audit L2 gradient projections."""
import argparse
import csv
from pathlib import Path
import numpy as np
from scipy.io import mmread

SHAPES=(("9x7",9,7,49),("13x11",13,11,94),("17x15",17,15,144))

def read(path):
    with Path(path).open(newline="",encoding="utf-8") as stream:
        return list(csv.DictReader(stream))

def matrix(path):
    return np.asarray(mmread(path),dtype=float)

def check_case(directory,nx,ny,expected,baseline):
    M=matrix(directory/"vector_mass.mtx")
    H=matrix(directory/"gradient_coupling.mtx")
    S=matrix(directory/"scalar_gradient_energy.mtx")
    G=matrix(directory/"projected_scalar_gram.mtx")
    Q=matrix(directory/"gradient_projector_basis.mtx")
    V=matrix(directory/"eigenvectors.mtx")
    rows=read(directory/"eigenmode_gradient.csv")
    singular=read(directory/"gradient_gram_eigenvalues.csv")
    n=nx*ny
    d=2*n
    if (M.shape!=(d,d) or H.shape!=(d,n) or S.shape!=(n,n)
        or G.shape!=(n,n) or V.shape!=(d,d) or Q.shape[0]!=d or
        len(rows)!=d or len(singular)!=n):
        raise ValueError("Matrix dimensions do not match full M6d EFG spectrum: "+str(directory))
    if not all(np.all(np.isfinite(a)) for a in (M,H,S,G,Q,V)):
        raise ValueError("Nonfinite matrices")
    if (not np.allclose(M,M.T,rtol=0,atol=1e-10)
        or not np.allclose(G,G.T,rtol=0,atol=1e-10)):
        raise ValueError("Asymmetric mass/projected scalar Gram")
    np.linalg.cholesky(M)
    reconstructed=H.T@np.linalg.solve(M,H)
    err=float(np.linalg.norm(G-reconstructed)/max(1,np.linalg.norm(G)))
    if err>1e-8:
        raise ValueError("Projected gradient Gram does not reconstruct: "+str(err))
    if (np.linalg.eigvalsh(S-G).min() <
          -1e-7*max(1,np.linalg.norm(S,ord=2))):
        raise ValueError("Best-fit vector gradient energy exceeds exact scalar gradient energy")
    gram=Q.T@M@Q
    max_m_orth=float(np.max(np.abs(gram-np.eye(Q.shape[1]))))
    if max_m_orth>5e-6:
        raise ValueError("Gradient basis not M-orthonormal")
    fraction=np.sum((Q.T@(M@V))**2,axis=0)/np.sum(V*(M@V),axis=0)
    expected_fractions=np.array([float(row["gradient_energy_fraction"]) for row in rows])
    if not np.all(np.isfinite(fraction)) or np.any((fraction<0)|(fraction>1+1e-6)):
        raise ValueError("Projected EFG modal energy is invalid")
    if np.max(np.abs(fraction-expected_fractions))>1e-6:
        raise ValueError("Archived modal gradient fractions fail reconstruction")
    if not all(int(row["index"])==i for i,row in enumerate(rows)):
        raise ValueError("Spectral indices do not align with eigenvectors")
    eigenvalues=np.array([float(row["lambda"]) for row in rows])
    if np.any(np.diff(eigenvalues)<0):
        raise ValueError("Eigenvalues not sorted")
    low=[r for r in rows if r["original_class"]=="positive_below_TM11"]
    if len(low)!=expected:
        raise ValueError("Historical count of low positive EFG modes changed")
    eps=1e-9
    baseline_rows=read(baseline/name_from(nx,ny)/"eigenvalues.csv")
    if len(baseline_rows)!=len(rows):
        raise ValueError("M6a original spectrum length changed")
    original=np.array([float(r["lambda"]) for r in baseline_rows])
    if np.max(np.abs(eigenvalues-original)/(1+np.abs(original)))>eps:
        raise ValueError("M6d K/M eigenvalues differ materially from published M6a baseline")
    g_values=np.array([float(r["sigma"]) for r in singular])
    g_kept=sum(r["kept"]=="true" for r in singular)
    if g_kept!=Q.shape[1] or np.max(np.abs(np.linalg.eigvalsh(G)-g_values))>1e-7:
        raise ValueError("Projected Gram rank or eigenvalues inconsistent")
    print("PASS",directory.name,"dofs",d,"gradient_rank",Q.shape[1],
          "Gram_error",f"{err:.3e}","M_orth_error",f"{max_m_orth:.3e}",
          "fraction_error",f"{np.max(np.abs(fraction-expected_fractions)):.3e}",
          "low_positive",len(low),flush=True)

def name_from(nx,ny):
    return str(nx)+"x"+str(ny)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign",type=Path,required=True)
    parser.add_argument("--baseline",type=Path,default=Path(
      "docs/results/independent_validation/2026-10-21_dc1d98b/efg"))
    args=parser.parse_args()
    for name,nx,ny,low in SHAPES:
        check_case(args.campaign/name,nx,ny,low,args.baseline)
    for sub in ("sensitive_1e-6","sensitive_1e-10"):
        check_case(args.campaign/"rank_sensitivity"/sub,17,15,144,args.baseline)

if __name__=="__main__":
    main()
