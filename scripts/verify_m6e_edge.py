#!/usr/bin/env python3
"""Independently reopen and audit every Whitney TM edge eigenpair/S/T CSV.

Reads native FEM edge_rect exports without modifying either solver.
Use PYTHONNOUSERSITE=1 on systems mixing NumPy 2 and distro SciPy.
"""
import argparse
import csv
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix

CASES=("9x7","13x11","17x15")


def csv_rows(path):
    with Path(path).open(newline="",encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def read_crs(path):
    sections={}
    for row in csv_rows(path):
        sections.setdefault(row["section"],[]).append(
            (int(row["index"]),row["value"]))
    def scalar(key):
        if key not in sections or len(sections[key])!=1:
            raise ValueError("Malformed CRS scalar field "+key)
        return int(sections[key][0][1])
    nrows=scalar("nrows")
    ncols=scalar("ncols")
    nnz=scalar("nnz")
    def indexed(key,count,cast):
        if len(sections.get(key,[]))!=count:
            raise ValueError("Unexpected CRS "+key+" length")
        rows=sorted(sections[key])
        if [r[0] for r in rows]!=list(range(count)):
            raise ValueError("Missing/duplicate CRS entries in "+key)
        return np.array([cast(r[1]) for r in rows])
    row_ptr=indexed("row_ptr",nrows+1,int)
    col_idx=indexed("col_idx",nnz,int)
    values=indexed("values",nnz,float)
    if row_ptr[0]!=0 or row_ptr[-1]!=nnz or np.any(np.diff(row_ptr)<0):
        raise ValueError("Invalid CRS row offsets")
    if np.any(col_idx<0) or np.any(col_idx>=ncols):
        raise ValueError("Invalid CRS column indexes")
    matrix=csr_matrix((values,col_idx,row_ptr),shape=(nrows,ncols)).toarray()
    if not np.all(np.isfinite(matrix)):
        raise ValueError("Nonfinite CRS data")
    return matrix


def read_modes(path,n):
    Z=np.zeros((n,n),dtype=float)
    counts=np.zeros(n,dtype=int)
    with Path(path).open(newline="",encoding="utf-8") as stream:
        reader=csv.DictReader(stream)
        for row in reader:
            rank=int(row["ordered_rank"])-1
            index=int(row["dof_index"])
            if not (0<=rank<n and 0<=index<n):
                raise ValueError("Whitney eigenvector DOF/index out of range")
            im=float(row["component_imag"])
            if abs(im)>1e-12:
                raise ValueError("Expected real symmetric FEM edge vectors")
            Z[index,rank]=float(row["component_real"])
            counts[rank]+=1
    if not np.all(counts==n):
        raise ValueError("Missing/duplicate Whitney eigenvector coefficient")
    if not np.all(np.isfinite(Z)):
        raise ValueError("Nonfinite Whitney eigenvector coefficients")
    return Z


def verify(case_path):
    edge=Path(case_path)/"edge_native"/"helmvec"/"rect"/"linop"
    prefix="edge_rect_tm"
    values=csv_rows(edge/(prefix+"_eigenvalues.csv"))
    eigen=np.array([float(row["lambda"]) for row in values])
    n=len(eigen)
    if n<10 or any(int(r["ordered_rank"])!=i+1 for i,r in enumerate(values)):
        raise ValueError("Missing/misaligned FEM edge TM eigenvalues")
    if np.any(np.diff(eigen)<-1e-10):
        raise ValueError("Whitney TM eigenvalues not sorted")
    S=read_crs(edge/(prefix+"_S_crs.csv"))
    T=read_crs(edge/(prefix+"_T_crs.csv"))
    V=read_modes(edge/(prefix+"_eigenvectors.csv"),n)
    if S.shape!=(n,n) or T.shape!=(n,n):
        raise ValueError("Inconsistent edge TM S/T/eigenvector sizes")
    if not np.allclose(S,S.T,rtol=0,atol=1e-9) or not np.allclose(T,T.T,rtol=0,atol=1e-9):
        raise ValueError("Whitney S/T symmetry violated")
    np.linalg.cholesky(T)
    gram=V.T@(T@V)
    ortho=float(np.max(np.abs(gram-np.eye(n))))
    if ortho>1e-6:
        raise ValueError("Edge FEM eigenvectors not mass-orthonormal")
    residual=S@V-(T@V)*eigen
    den=(np.linalg.norm(S)+np.max(np.abs(eigen))*np.linalg.norm(T))*max(
        1,np.max(np.linalg.norm(V,axis=0)))
    rel=float(np.max(np.linalg.norm(residual,axis=0))/max(1e-300,den))
    if rel>1e-7:
        raise ValueError("Independent Whitney EVP residual too large: "+str(rel))
    if np.any(eigen < -1e-7*np.max(np.abs(eigen))):
        raise ValueError("Materially negative FEM edge curl-curl eigenvalue")
    return dict(case=case_path.name,modes=n,M_orthogonality_max=ortho,
                relative_EVP_residual_max=rel)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--campaign",type=Path,required=True)
    a=p.parse_args()
    for label in CASES:
        r=verify(a.campaign/label)
        print("PASS",r,flush=True)


if __name__=="__main__":
    main()
