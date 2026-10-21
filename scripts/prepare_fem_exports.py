#!/usr/bin/env python3
"""Package exact-double TM FEM captures; never solve or modify their eigenpairs.

The three run directories must contain capture/call_02_*.mtx/csv and native
helm10_rect output in raw/helm10/rect/.  Call 01 is TE and is excluded here.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import io
import json
import platform
from pathlib import Path
import shutil
import subprocess
import sys
from datetime import datetime, timezone
import numpy as np

CONFIGS = {"coarse": (16, 8), "medium": (30, 15), "fine": (40, 20)}


def read_array(path):
    with open(path) as stream:
        for line in stream:
            if line.strip() and not line.startswith("%"):
                shape = tuple(map(int, line.split()))
                break
        values = np.loadtxt(stream, dtype=float)
    return values.reshape(shape, order="F")


def gzip_copy(src, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(src, "rb") as source, open(dest, "wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0) as target:
            shutil.copyfileobj(source, target)


def write_csv(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        raw = open(path, "wb")
        stream = io.TextIOWrapper(gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0), newline="")
    else:
        raw = None
        stream = open(path, "w", newline="")
    try:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(header)
        for row in rows:
            writer.writerow([format(float(x), ".17g") if isinstance(x, (float, np.floating)) else x for x in row])
    finally:
        stream.close()
        if raw:
            raw.close()


def version(command):
    result = subprocess.run(command, capture_output=True, text=True)
    return (result.stdout or result.stderr).strip().splitlines()[0] if result.returncode == 0 else None


def read_native_vtk(path):
    """Read actual executable mesh provenance; VTK display floats are rounded."""
    tokens=path.read_text().split()
    start=tokens.index("POINTS")
    count=int(tokens[start+1])
    points=np.array(tokens[start+3:start+3+count*3],dtype=float).reshape(count,3)
    start=tokens.index("CELLS")
    count=int(tokens[start+1])
    cells=np.array(tokens[start+3:start+3+count*4],dtype=int).reshape(count,4)
    if not np.all(cells[:,0]==3): raise ValueError("Nontriangular VTK cell")
    return points[:,:2],cells[:,1:]


def read_native_crs(path):
    sections={}
    with open(path) as stream:
        for row in csv.DictReader(stream):
            sections.setdefault(row["section"],[]).append(float(row["value"]))
    n=int(sections["nrows"][0])
    matrix=np.zeros((n,n))
    ptr=np.array(sections["row_ptr"],dtype=int)
    col=np.array(sections["col_idx"],dtype=int)
    values=np.array(sections["values"])
    for i in range(n): matrix[i,col[ptr[i]:ptr[i+1]]]=values[ptr[i]:ptr[i+1]]
    return matrix


def prepare(staging, out, name, nx, ny, repo):
    source = staging / name
    capture = source / "capture"
    native = source / "raw/helm10/rect"
    dest = out / name
    dest.mkdir(parents=True, exist_ok=True)
    k = read_array(capture / "call_02_stiffness_matrix.mtx")
    m = read_array(capture / "call_02_mass_matrix.mtx")
    v = read_array(capture / "call_02_eigenvectors.mtx")
    eigen = np.genfromtxt(capture / "call_02_eigenvalues.csv", delimiter=",", names=True)
    lam = eigen["lambda"]
    n = (nx - 1) * (ny - 1)
    if k.shape != (n,n) or m.shape != (n,n) or v.shape != (n,n) or len(lam) != n:
        raise ValueError(f"{name}: inconsistent dimensions")
    kv, mv = k @ v, m @ v
    absolute = np.linalg.norm(kv - mv * lam, axis=0)
    denom = np.linalg.norm(kv, axis=0) + np.abs(lam) * np.linalg.norm(mv, axis=0)
    residual = absolute / np.maximum(denom, np.finfo(float).tiny)
    gram = v.T @ mv
    orthogonality = float(np.max(np.abs(gram - np.eye(n))))
    chol = np.linalg.cholesky(m)
    mass_min = float(np.linalg.eigvalsh(m)[0])
    scale = max(1.0, float(np.max(np.abs(lam))))
    threshold = 1e-10 * scale
    classes = np.where(lam < -threshold, "negative", np.where(np.abs(lam) <= threshold,"near_zero","positive"))
    nodes, dofs, triangles = [], [], []
    node_dof = np.full((nx+1)*(ny+1), -1, dtype=int)
    for j in range(ny+1):
        for i in range(nx+1):
            node = j*(nx+1)+i
            x, y = 2.0*i/nx, 1.0*j/ny
            boundary = int(i == 0 or i == nx or j == 0 or j == ny)
            nodes.append((node,x,y,boundary))
            if not boundary:
                node_dof[node] = len(dofs)
                dofs.append((len(dofs),node,x,y))
    for j in range(ny):
        for i in range(nx):
            n00=j*(nx+1)+i
            n10,n01,n11=n00+1,n00+nx+1,n00+nx+2
            triangles.extend([(len(triangles),n00,n10,n11),(len(triangles)+1,n00,n11,n01)])
    points=np.array([[r[1],r[2]] for r in nodes])
    vtk_path=native/"vtk/tm11_rect_sv.vtk"
    vtk_points,vtk_triangles=read_native_vtk(vtk_path)
    mesh_vtk_coordinates_error=float(np.max(np.abs(vtk_points-points)))
    mesh_vtk_connectivity_identical=np.array_equal(vtk_triangles,np.array(triangles)[:,1:])
    native_k=read_native_crs(native/"linop/helm10_rect_tm_S_crs.csv")
    native_m=read_native_crs(native/"linop/helm10_rect_tm_T_crs.csv")
    native_k_error=float(np.max(np.abs(native_k-k)))
    native_m_error=float(np.max(np.abs(native_m-m)))
    # Independently reconstruct P1 element matrices from the exported mesh/map.
    reconstructed_k=np.zeros_like(k)
    reconstructed_m=np.zeros_like(m)
    element_gradients=[]
    for _, *indices in triangles:
        xy=points[indices]
        twice_area=np.linalg.det(np.array([[xy[1,0]-xy[0,0],xy[2,0]-xy[0,0]],
                                           [xy[1,1]-xy[0,1],xy[2,1]-xy[0,1]]]))
        area=abs(twice_area)/2
        grad=np.array([[xy[1,1]-xy[2,1],xy[2,0]-xy[1,0]],
                       [xy[2,1]-xy[0,1],xy[0,0]-xy[2,0]],
                       [xy[0,1]-xy[1,1],xy[1,0]-xy[0,0]]])/twice_area
        element_gradients.append(grad)
        ke=area*(grad@grad.T)
        me=(area/12)*(np.ones((3,3))+np.eye(3))
        local_dofs=node_dof[indices]
        for ii,di in enumerate(local_dofs):
            if di<0: continue
            for jj,dj in enumerate(local_dofs):
                if dj>=0:
                    reconstructed_k[di,dj]+=ke[ii,jj]
                    reconstructed_m[di,dj]+=me[ii,jj]
    mesh_k_error=float(np.max(np.abs(reconstructed_k-k)))
    mesh_m_error=float(np.max(np.abs(reconstructed_m-m)))
    checks = {
        "all_finite": all(np.isfinite(a).all() for a in (k,m,v,lam)),
        "complete_spectrum_dimensions": len(lam)==n and v.shape==(n,n),
        "stiffness_symmetric": bool(np.max(np.abs(k-k.T))<1e-13),
        "mass_symmetric": bool(np.max(np.abs(m-m.T))<1e-13),
        "mass_positive_definite": bool(np.all(np.diag(chol)>0) and mass_min>0),
        "all_eigenpair_residuals": bool(np.max(residual)<1e-10),
        "mass_orthonormality": orthogonality<1e-10,
        "mesh_stiffness_reconstruction": mesh_k_error<1e-12,
        "mesh_mass_reconstruction": mesh_m_error<1e-14,
        "positive_dirichlet_spectrum": bool(np.all(lam>0)),
        "native_vtk_connectivity": bool(mesh_vtk_connectivity_identical),
        "native_vtk_coordinates": mesh_vtk_coordinates_error<1e-5,
        "native_export_stiffness": native_k_error<1e-13,
        "native_export_mass": native_m_error<1e-15,
    }
    if not all(checks.values()):
        raise AssertionError(f"{name}: {checks}")
    for component in ("stiffness_matrix", "mass_matrix", "eigenvectors"):
        gzip_copy(capture / f"call_02_{component}.mtx",dest / f"{component}.mtx.gz")
    gzip_copy(vtk_path,dest/"native_exports/tm11_rect_sv.vtk.gz")
    write_csv(dest/"eigenvalues.csv",["index","lambda","kc","classification","relative_residual"],
              [(i,lam[i],float(np.sqrt(max(0.0,lam[i]))),classes[i],residual[i]) for i in range(n)])
    write_csv(dest/"mesh_nodes.csv",["node_index","x","y","is_boundary"],nodes)
    write_csv(dest/"mesh_elements.csv",["element_index","node0","node1","node2"],triangles)
    write_csv(dest/"dof_map.csv",["dof","node_index","x","y"],dofs)
    write_csv(dest/"mode_diagnostics.csv",["index","relative_residual","absolute_residual","mass_norm_squared","gradient_energy","dirichlet_trace_max","classification"],
              [(i,residual[i],absolute[i],gram[i,i],v[:,i]@kv[:,i],0.0,classes[i]) for i in range(n)])
    # All coefficients remain available above; physical Ez and exact element
    # gradients are additionally saved for the first 12 computed TM modes.
    count=min(12,n)
    nodal=np.zeros((len(nodes),count))
    free=np.flatnonzero(node_dof>=0)
    nodal[free]=v[node_dof[free],:count]
    write_csv(dest/"fields_nodes.csv.gz",["node_index","x","y"]+[f"Ez_mode_{i:03d}" for i in range(count)],
              [(i,*points[i],*nodal[i]) for i in range(len(nodes))])
    field_rows=[]
    for (element,*indices),grad in zip(triangles,element_gradients):
        values=grad.T@nodal[indices]
        field_rows.append((element,*points[indices].mean(axis=0),*values.T.ravel()))
    write_csv(dest/"fields_elements.csv.gz",["element_index","x","y"]+[f"dEz_d{axis}_mode_{i:03d}" for i in range(count) for axis in ("x","y")],field_rows)
    # Preserve native TM CSVs with their original 16-digit writer precision;
    # the authoritative 17-digit Matrix Market files above come from LAPACK.
    for original in sorted((native/"linop").glob("helm10_rect_tm_*.csv")):
        gzip_copy(original,dest/"native_exports"/(original.name+".gz"))
    for original in sorted((native/"csv").glob("helm10_rect_fields_TM_*.csv")):
        gzip_copy(original,dest/"native_exports"/(original.name+".gz"))
    for original in [native/"csv/helm10_rect_modes.csv",native/"run_timing.csv"]:
        gzip_copy(original,dest/"native_exports"/(original.name+".gz"))
    for log in ("stdout.log","stderr.log"):
        if (source/log).exists(): shutil.copy2(source/log,dest/log)
    git_commit=subprocess.check_output(["git","-C",str(repo),"rev-parse","HEAD"],text=True).strip()
    source_status=subprocess.check_output(["git","-C",str(repo),"status","--porcelain"],text=True).splitlines()
    tracked_modified=any(not item.startswith("??") for item in source_status)
    transparency=None
    if name=="coarse" and (staging/"uninstrumented/raw/helm10/rect/linop").exists():
        transparency={"independent_uninstrumented_command": ["helm10_rect","--ar-m","2.0","--nx","16","--ny","8","--nmodos","12","--backend","closed-form"],
                      "native_outputs_byte_identical":True,"sha256":{}}
        for filename in ("helm10_rect_tm_eigenvalues.csv","helm10_rect_tm_eigenvectors.csv",
                         "helm10_rect_tm_S_crs.csv","helm10_rect_tm_T_crs.csv"):
            original=(native/"linop"/filename).read_bytes()
            control=(staging/"uninstrumented/raw/helm10/rect/linop"/filename).read_bytes()
            if original!=control: raise AssertionError(f"Interceptor altered {filename}")
            transparency["sha256"][filename]=hashlib.sha256(original).hexdigest()
    metadata={
        "method":"independent scalar P1 triangular FEM", "family":"TM", "boundary_condition":"homogeneous Dirichlet Ez=0",
        "geometry":{"a":2.0,"b":1.0,"unit":"m"},"nx":nx,"ny":ny,
        "nodes":len(nodes),"triangles":len(triangles),"dofs":n,"matrix_dimensions":[n,n],"eigenpairs":n,
        "source_repository":str(repo.resolve()),"source_remote":"https://github.com/Lucas-Sperotto/fem-eigen-em-2d3d-whitney.git",
        "source_commit":git_commit,"source_modified":tracked_modified,
        "source_worktree_status_preserved":source_status,
        "command":["helm10_rect","--ar-m","2.0","--nx",str(nx),"--ny",str(ny),"--nmodos","12","--backend","closed-form"],
        "material":{"eps_r":1.0,"mu_r":1.0},"backend":"closed-form",
        "solver":{"library":"LAPACKE_dsygv","itype":1,"layout":"column-major","jobz":"V","uplo":"U","info":0,
                  "precision":"IEEE 754 binary64", "algorithm":"Cholesky reduction to standard symmetric EVP, DSYEV implicit QL/QR eigensolve, back transformation",
                  "normalization":"V.T @ T @ V = I; LAPACK normalization is preserved",
                  "input":"S u = lambda T u; lambda = kc squared", "capture_call":2},
        "export":{"matrix_market":"array real general; full matrix; column-major text order",
                  "text_significant_digits":17,"native_csv_significant_digits":16,
                  "index_base":0,"eigenvector_column":"eigenvalues.csv index",
                  "capture":"transparent LD_PRELOAD forwards real LAPACKE_dsygv; records S/T before and eigenpairs after without numerical substitution",
                  "coefficient_order":"interior nodes in row-major mesh node order; one Ez DOF per node",
                  "field_convention":"Ez is mass-normalized nodal scalar; dEz gradients are exact P1 element gradients; first 12 modes saved, all modes reconstructible from coefficients",
                  "degenerate_modes":"individual vectors may rotate; compare physical two-dimensional subspaces",
                  "raw_native_outputs":"native_exports/*.csv.gz; unmodified executable exports"},
        "diagnostics":{"maximum_relative_residual":float(np.max(residual)),"mass_orthogonality_max_abs_error":orthogonality,
                       "mass_smallest_eigenvalue":mass_min,"stiffness_symmetry_max_abs_error":float(np.max(np.abs(k-k.T))),
                       "mass_symmetry_max_abs_error":float(np.max(np.abs(m-m.T))),
                       "mesh_stiffness_reconstruction_max_abs_error":mesh_k_error,"mesh_mass_reconstruction_max_abs_error":mesh_m_error,
                       "native_vtk_coordinates_max_abs_error":mesh_vtk_coordinates_error,"native_vtk_connectivity_identical":bool(mesh_vtk_connectivity_identical),
                       "native_export_stiffness_max_abs_error":native_k_error,"native_export_mass_max_abs_error":native_m_error,
                       "near_zero_threshold":threshold,"near_zero_count":int(np.sum(classes=="near_zero")),"negative_count":int(np.sum(classes=="negative"))},
        "checks":{"passed":sum(checks.values()),"failed":sum(not item for item in checks.values()),"details":checks},
        "tools":{"python":sys.version.split()[0],"numpy":np.__version__,"cmake":version(["cmake","--version"]),"compiler":version(["c++","--version"]),"platform":platform.platform()},
        "executed_at_utc":datetime.fromtimestamp((capture/"call_02_eigenvalues.csv").stat().st_mtime,timezone.utc).isoformat(),
        "threads":{"OPENBLAS_NUM_THREADS":"1"},
        "interceptor_transparency_control":transparency,
    }
    (dest/"run_metadata.json").write_text(json.dumps(metadata,indent=2,allow_nan=False)+"\n")
    print(f"{name}: n={n}, residual={np.max(residual):.3e}, orthogonality={orthogonality:.3e}, checks={sum(checks.values())}/{len(checks)}")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--staging",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--fem-repo",type=Path,default=Path(__file__).resolve().parents[2]/"fem-eigen-em-2d3d-whitney")
    args=parser.parse_args()
    for name,(nx,ny) in CONFIGS.items(): prepare(args.staging,args.out,name,nx,ny,args.fem_repo)

if __name__=="__main__": main()
