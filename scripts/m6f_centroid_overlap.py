#!/usr/bin/env python3
"""M6f exploratory physical-vector cross-check on the SAME Whitney triangle centroids.

The centroid rule is NOT a converged common quadrature. It is a first,
reproducible, low-cost comparison of Hx,Hy. Degenerate pairs are compared
as invariant 2D L2 subspaces, not arbitrary one-to-one eigenvectors.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import numpy as np

CASES=(("9x7",9,7,8,4),("13x11",13,11,12,6),("17x15",17,15,16,8))
SIMPLE=("TM11","TM21","TM31","TM12")
ALL=SIMPLE+("TM22_A","TM41_B")
MNS={"TM11":(1,1),"TM21":(2,1),"TM31":(3,1),"TM12":(1,2),
     "TM22_A":(2,2),"TM41_B":(4,1)}

def rows(path):
    with Path(path).open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))

def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for piece in iter(lambda:f.read(1<<20),b""):
            h.update(piece)
    return h.hexdigest()

def vector_norm(v,w):
    energy=float(np.dot(w,v*v))
    if not math.isfinite(energy) or energy<=1e-25:
        raise ValueError("Zero/nonfinite vector field energy")
    return math.sqrt(energy)

def unit(v,w):
    return np.asarray(v,dtype=float)/vector_norm(np.asarray(v,dtype=float),w)

def correlation(a,b,w):
    a=unit(a,w);b=unit(b,w)
    r=float(np.dot(w,a*b))
    if not math.isfinite(r) or abs(r)>1+1e-8:
        raise ValueError("Impossible vector L2 correlation")
    return min(1.0,abs(r))

def principal_correlations(A,B,w):
    """Coordinates are columns; report basis-invariant L2 principal cosines."""
    w=np.asarray(w,dtype=float)
    if A.ndim!=2 or B.ndim!=2 or A.shape!=B.shape or A.shape[1]!=2 or len(w)!=A.shape[0]:
        raise ValueError("Expected two vector modes in common sample points")
    def orthonormal(mat):
        scaled=np.sqrt(w)[:,None]*mat
        q,r=np.linalg.qr(scaled,mode="reduced")
        if np.min(np.linalg.svd(r,compute_uv=False))<1e-9*np.max(np.linalg.svd(r,compute_uv=False)):
            raise ValueError("Degenerate pair has numerically rank-deficient fields")
        return q
    U=orthonormal(A)
    V=orthonormal(B)
    sv=np.linalg.svd(U.T@V,compute_uv=False)
    if not np.all(np.isfinite(sv)) or np.any(sv>1+1e-8):
        raise ValueError("Invalid principal cosines")
    return np.clip(np.sort(sv)[::-1],0,1)

def analytic_H(x,y,m,n):
    kx=m*math.pi/2
    ky=n*math.pi
    return np.column_stack((
        -ky*np.sin(kx*x)*np.cos(ky*y),
         kx*np.cos(kx*x)*np.sin(ky*y))).ravel()

def native_H(path):
    rec=rows(path)
    if not rec:raise RuntimeError("Missing FEM native centroid field "+str(path))
    cells=np.array([int(r["cell_id"]) for r in rec],dtype=int)
    if not np.array_equal(cells,np.arange(len(rec))):
        raise RuntimeError("Discontinuous/reordered native FEM cells")
    x=np.array([float(r["xc_m"]) for r in rec])
    y=np.array([float(r["yc_m"]) for r in rec])
    field=np.array([(float(r["Hx"]),float(r["Hy"])) for r in rec])
    if not np.all(np.isfinite(field)) or not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise RuntimeError("Nonfinite native FEM-H vector field")
    return cells,x,y,field.ravel()

def efg_H(path,mode_indices):
    rec=rows(path)
    cells=np.array([int(r["cell_id"]) for r in rec],dtype=int)
    x=np.array([float(r["x_efg"]) for r in rec])
    y=np.array([float(r["y_efg"]) for r in rec])
    F={}
    for mode,index in mode_indices.items():
        F[mode]=np.array([
            (float(r["mode"+str(index)+"_hx"]),
             float(r["mode"+str(index)+"_hy"])) for r in rec]).ravel()
        if not np.all(np.isfinite(F[mode])):
            raise RuntimeError("Nonfinite EFG samples")
    return cells,x,y,F

def compare_case(campaign,binary,out,label,nx,ny,cx,cy):
    native_dir=campaign/label/"edge_native"/"helmvec"/"rect"/"csv"
    matched=rows(native_dir/"edge_rect_modes.csv")
    indexed={int(r["eig_index"]):r for r in matched if r["family"]=="TM"}
    comp=rows(campaign/label/"comparison.csv")
    if len(comp)!=6 or {r["mode"] for r in comp}!=set(ALL):
        raise RuntimeError("Missing M6e six-mode comparison")
    edge_modes={}
    efg_idx={}
    for r in comp:
        mode=r["mode"]
        if r["edge_index"]=="" or int(r["edge_index"]) not in indexed:
            raise RuntimeError("M6e selected Whitney field not in native 12-mode export: "+mode)
        path=native_dir/indexed[int(r["edge_index"])]["fields_csv_file"]
        if not path.is_file():
            raise RuntimeError("Missing native FEM TM field "+str(path))
        edge_modes[mode]=path
        efg_idx[mode]=int(r["efg_index"])
    # Native FEM centroids serve as the common exploratory set.
    _,xf,yf,_=native_H(edge_modes["TM11"])
    sample_file=edge_modes["TM11"]
    out.mkdir(parents=True,exist_ok=True)
    cmd=[str(binary.resolve()),str(nx),str(ny),str(sample_file),
         str(out/"efg_samples.csv")]+[str(efg_idx[key]) for key in ALL]
    proc=subprocess.run(cmd,text=True,capture_output=True,check=False,timeout=600)
    (out/"efg_stdout.txt").write_text(proc.stdout,encoding="utf-8")
    (out/"efg_stderr.txt").write_text(proc.stderr,encoding="utf-8")
    if proc.returncode:
        raise RuntimeError("EFG centroid sampler error: "+proc.stderr)
    cell,xe,ye,fields=efg_H(out/"efg_samples.csv",efg_idx)
    if len(cell)!=2*cx*cy or not np.array_equal(cell,np.arange(len(cell))) or \
       not np.allclose(xe,2*xf,atol=1e-12,rtol=0) or \
       not np.allclose(ye,2*yf,atol=1e-12,rtol=0):
        raise RuntimeError("Nonmatching FEM and EFG physical sample points")
    # Area of each original triangular cell = 0.5*(1/cx)*(0.5/cy)
    # Physical scaling x2 in both axes multiplies area by 4.
    weight_per_cell=4*0.5*(1/cx)*(0.5/cy)
    if abs(weight_per_cell*len(cell)-2)>1e-12:
        raise RuntimeError("Common-grid area mismatch")
    weights=np.repeat(weight_per_cell,2*len(cell))
    fem_fields={}
    hashes={}
    for mode,path in edge_modes.items():
        ids,x,y,F=native_H(path)
        if not np.array_equal(ids,cell) or not np.allclose(x,xf) or not np.allclose(y,yf):
            raise RuntimeError("Whitney field grids differ between modes")
        fem_fields[mode]=F
        hashes[mode]=dict(native_file=str(path.relative_to(campaign)),
                         sha256=digest(path),efg_index=efg_idx[mode],
                         edge_index=int(next(r["edge_index"] for r in comp if r["mode"]==mode)))
    simple=[]
    for mode in SIMPLE:
        m,n=MNS[mode]
        ana=analytic_H(xe,ye,m,n)
        rho=correlation(fields[mode],fem_fields[mode],weights)
        simple.append(dict(case=label,mode=mode,
                           rho_efg_whitney=rho,
                           normalized_l2_difference=math.sqrt(max(0,2-2*rho)),
                           rho_efg_analytic=correlation(fields[mode],ana,weights),
                           rho_whitney_analytic=correlation(fem_fields[mode],ana,weights),
                           centroid_cells=len(cell)))
    A=np.column_stack([fields[m] for m in ALL[4:]])
    B=np.column_stack([fem_fields[m] for m in ALL[4:]])
    ref=np.column_stack([analytic_H(xe,ye,*MNS[m]) for m in ALL[4:]])
    s=principal_correlations(A,B,weights)
    sa=principal_correlations(A,ref,weights)
    sb=principal_correlations(B,ref,weights)
    pair=dict(case=label,centroid_cells=len(cell),
              min_cross_subspace_correlation=float(min(s)),
              max_cross_principal_angle_deg=float(np.degrees(np.arccos(min(s)))),
              min_efg_analytic_subspace_correlation=float(min(sa)),
              min_whitney_analytic_subspace_correlation=float(min(sb)),
              first_cross_principal_correlation=float(max(s)),
              second_cross_principal_correlation=float(min(s)))
    provenance=dict(
        m6e_source=str(campaign.resolve()),edge_mesh=[cx,cy],efg_cloud=[nx,ny],
        native_edge_domain=[1.0,0.5],efg_domain=[2.0,1.0],
        edge_to_efg_coordinates_multiplier=2,
        cell_weight_efg=weight_per_cell,
        sample_rule="native_Whitney_triangle_centroids_one_point_not_converged_quadrature",
        vec_field="TM_H_t",
        field_status="Whitney_native_per_mode_unit_peak_EFG_original_mass_normalized",
        amplitude_status="Each vector field renormalized in discrete centroid L2",
        matching="eigenvalue_selected_from_M6e_and_native_edge_index_not_shape-certified",
        reference_fields=hashes,
        efg_sample_sha256=digest(out/"efg_samples.csv")
    )
    (out/"provenance.json").write_text(json.dumps(provenance,ensure_ascii=False,indent=2)+"\n",
                                       encoding="utf-8")
    return simple,pair

def write_csv(path,data):
    with Path(path).open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(data[0]))
        writer.writeheader();writer.writerows(data)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--efg-binary",type=Path,required=True)
    p.add_argument("--m6e-campaign",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    binary=a.efg_binary.resolve()
    campaign=a.m6e_campaign.resolve()
    dest=a.out.resolve()
    if not binary.is_file() or not campaign.is_dir():
        p.error("Missing compiled EFG sampler or prior M6e raw FEM exports")
    dest.mkdir(parents=True,exist_ok=True)
    simple,pairs=[],[]
    for case,nx,ny,cx,cy in CASES:
        print("M6f common centroid vector comparison:",case,flush=True)
        s,d=compare_case(campaign,binary,dest/case,case,nx,ny,cx,cy)
        simple+=s;pairs.append(d)
    write_csv(dest/"simple_modes.csv",simple)
    write_csv(dest/"doublet_subspaces.csv",pairs)
    lines=[
        "# M6f — Diagnóstico exploratório de campos TM-H nos centroides Whitney",
        "",
        "**Atenção:** uma amostra por triângulo não é quadratura espacial",
        "convergida e NÃO constitui certificação física M6f integral.",
        "Campos EFG avaliados nos MESMOS pontos, com x,y mapeados por fator",
        "2 (guia FEM 1×0,5 para guia EFG 2×1). Valores Whitney de origem",
        "são normalizados por pico; ambos são renormalizados em L² de",
        "centroides com pesos das áreas físicas dos triângulos.",
        "Sinais e escalas não alteram as correlações.",
        "",
        "| Caso | Modo | |rho| EFG×Whitney | Erro L2 normalizado | |rho| EFG×analítico | |rho| Whitney×analítico |",
        "|---|---|---:|---:|---:|---:|",
        *[f"| {r['case']} | {r['mode']} | {r['rho_efg_whitney']:.8f} | {r['normalized_l2_difference']:.6f} | {r['rho_efg_analytic']:.8f} | {r['rho_whitney_analytic']:.8f} |" for r in simple],
        "",
        "## Dubleto TM22/TM41 (comparação conjunta, sem rotular colunas)",
        "",
        "| Caso | menor correlação principal EFG×Whitney | maior ângulo (°) | menor EFG×analítico | menor Whitney×analítico |",
        "|---|---:|---:|---:|---:|",
        *[f"| {r['case']} | {r['min_cross_subspace_correlation']:.8f} | {r['max_cross_principal_angle_deg']:.6f} | {r['min_efg_analytic_subspace_correlation']:.8f} | {r['min_whitney_analytic_subspace_correlation']:.8f} |" for r in pairs],
        "",
        "Os ângulos principais foram calculados pela SVD das bases",
        "ortonormalizadas no L² dos centroides. Rotação interna de",
        "autovetores do dubleto não altera os valores singulares.",
        "A falha de divergência do M6b deve continuar registrada, mesmo",
        "que o ângulo de subespaço seja pequeno.",
        "",
        "## Próximos testes indispensáveis",
        "",
        "1. Avaliar Whitney e EFG em quadratura Gauss na MESMA malha",
        "  e fazer refinamento de resolução do espaço de amostragem.",
        "2. Conferir H.normal e div(H) entre métodos, especialmente",
        "  a direção física pior do dubleto EFG.",
        "3. Rastrear o matching dos modos sem depender apenas de autovalores",
        "  próximos; a primeira campanha é explicitamente exploratória.",
        ""
    ]
    (dest/"REPORT.md").write_text("\n".join(lines),encoding="utf-8")
    print("\n".join(lines),flush=True)

if __name__=="__main__":
    main()
