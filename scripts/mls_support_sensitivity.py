#!/usr/bin/env python3
"""Small, controlled diagnostic sweep of MLS support factor; not thesis reproduction.

Current basis: degree-1 polynomial MLS, Wendland C2, support=hmax*factor.
Only parameter changed between runs is factor. No mode filtering; both
spectral-nearest and strongest curl overlap are retained side-by-side.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import subprocess
from audit_m2_m5_reproduction import read_original, M2_ORIGINAL, M5_ORIGINAL

GRIDS=((5,3),(7,4))
SUPPORTS=(3.0,3.3,3.6)
GAUSS=2
SUBDIVISION=2

def read_csv(path):
    with path.open(encoding="utf-8",newline="") as f:
        return list(csv.DictReader(f))

def parse_kv(path):
    return dict(line.split("=",1) for line in path.read_text(encoding="utf-8").splitlines() if "=" in line)

def run_case(binary,args,folder):
    folder.mkdir(parents=True,exist_ok=True)
    cmd=[str(binary.resolve()),*map(str,args),str(folder)]
    p=subprocess.run(cmd,text=True,capture_output=True,check=False,timeout=300)
    (folder/"stdout.txt").write_text(p.stdout,encoding="utf-8")
    (folder/"stderr.txt").write_text(p.stderr,encoding="utf-8")
    if p.returncode:
        raise RuntimeError(f"Failed: {cmd}\n{p.stderr}")
    return p

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cartesian-binary",required=True,type=Path)
    parser.add_argument("--rot-irrot-binary",required=True,type=Path)
    parser.add_argument("--out",required=True,type=Path)
    a=parser.parse_args()
    if not a.cartesian_binary.is_file() or not a.rot_irrot_binary.is_file():
        parser.error("Build the Cartesian and M5 executables first.")
    source={ "M2":read_original(M2_ORIGINAL), "M5":read_original(M5_ORIGINAL)}
    out=a.out.resolve()
    data=[]
    for nx,ny in GRIDS:
        for factor in SUPPORTS:
            # Cartesian tm_validation: nx ny degree support gauss out sub
            tag=f"{ny}x{nx}_h{factor:.1f}".replace(".","p")
            directory=out/"M2"/tag
            cmd=[str(a.cartesian_binary.resolve()),str(nx),str(ny),"1",
                 str(factor),str(GAUSS),str(directory),str(SUBDIVISION)]
            directory.mkdir(parents=True,exist_ok=True)
            p=subprocess.run(cmd,capture_output=True,text=True,timeout=300)
            (directory/"stdout.txt").write_text(p.stdout,encoding="utf-8")
            (directory/"stderr.txt").write_text(p.stderr,encoding="utf-8")
            if p.returncode:raise RuntimeError(f"M2 {tag}: {p.stderr}")
            modes={(int(r["m"]),int(r["n"])):r for r in read_csv(directory/"tm_matches.csv")}
            meta=parse_kv(directory/"run_metadata.txt")
            if int(meta["subcells_per_node_interval"])!=SUBDIVISION:
                raise ValueError("Wrong M2 background integration setting")
            for mode,mn in (("TM11",(1,1)),("TM21",(2,1))):
                r=modes[mn]
                author=source["M2"][mode,ny,nx,GAUSS]
                shape=100*(float(r["strongest_lambda"])/(1.25*math.pi**2 if mode=="TM11" else 2*math.pi**2)-1)
                nearest=100*(float(r["nearest_lambda"])/(1.25*math.pi**2 if mode=="TM11" else 2*math.pi**2)-1)
                data.append(dict(formulation="M2",mode=mode,nodes_y=ny,nodes_x=nx,
                    gauss=GAUSS,support_factor=factor,author_error_pct=author,
                    nearest_error_pct=nearest,field_error_pct=shape,
                    field_author_gap_pp=shape-author,
                    field_curl_rho=float(r["strongest_overlap"]),
                    mass_rank="",mass_nullity="",near_zero_curl=""))
            # Rot/irrot: nx ny gauss support sub out tolerance=default
            directory=out/"M5"/tag
            cmd=[str(a.rot_irrot_binary.resolve()),str(nx),str(ny),
                 str(GAUSS),str(factor),str(SUBDIVISION),str(directory)]
            directory.mkdir(parents=True,exist_ok=True)
            p=subprocess.run(cmd,capture_output=True,text=True,timeout=300)
            (directory/"stdout.txt").write_text(p.stdout,encoding="utf-8")
            (directory/"stderr.txt").write_text(p.stderr,encoding="utf-8")
            if p.returncode:raise RuntimeError(f"M5 {tag}: {p.stderr}")
            modes={(int(r["m"]),int(r["n"])):r for r in read_csv(directory/"tm_matches.csv")}
            meta=parse_kv(directory/"run_metadata.txt")
            for mode,mn in (("TM11",(1,1)),("TM21",(2,1))):
                r=modes[mn]
                author=source["M5"][mode,ny,nx,GAUSS]
                data.append(dict(formulation="M5",mode=mode,nodes_y=ny,nodes_x=nx,
                    gauss=GAUSS,support_factor=factor,author_error_pct=author,
                    nearest_error_pct=float(r["nearest_signed_error_pct"]),
                    field_error_pct=float(r["best_curl_signed_error_pct"]),
                    field_author_gap_pp=float(r["best_curl_signed_error_pct"])-author,
                    field_curl_rho=float(r["best_curl_rho"]),
                    mass_rank=int(meta["retained_rank"]),
                    mass_nullity=int(meta["removed_mass_nullity"]),
                    near_zero_curl=int(meta["near_zero_eigenvalues"])))
            print(f"PASS M2+M5 grid={ny}x{nx} factor={factor}",flush=True)
    file=out/"support_sensitivity.csv"
    with file.open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(data[0]))
        writer.writeheader();writer.writerows(data)
    lines=[
        "# Auditoria de sensibilidade ao suporte MLS (M2/M5)",
        "",
        "O parâmetro alterado é SOMENTE o multiplicador do raio de suporte MLS.",
        "Nuvens 3×5 e 4×7, Gauss 2×2, quatro subcélulas por retângulo,",
        "MLS linear e Wendland C2. Fatores 3,0; 3,3; 3,6 são",
        "**escolhas experimentais nossas, não valores atribuídos a Al Maleh**.",
        "Erros percentuais são ASSINADOS em lambda e diferenças em pp.",
        "A melhor correlação de curl não garante um campo fisicamente correto.",
        "",
        "| Formulação | TM | Nós | Suporte | Autor % | Nosso melhor curl % | Δ autor pp | rho curl | Rank da massa | Zeros curl |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in data:
        rank=str(r["mass_rank"]) if r["mass_rank"]!="" else "—"
        zero=str(r["near_zero_curl"]) if r["near_zero_curl"]!="" else "—"
        lines.append(f"| {r['formulation']} | {r['mode']} | {r['nodes_y']}×{r['nodes_x']} | {r['support_factor']:.1f} | {r['author_error_pct']:+.3f} | {r['field_error_pct']:+.3f} | {r['field_author_gap_pp']:+.3f} | {r['field_curl_rho']:.6f} | {rank} | {zero} |")
    lines+=["","Para cada caso foram guardados espectros completos, matrizes e saídas nativas;",
            "nenhum parâmetro foi otimizado para obter a referência do autor.",
            "Esta varredura NÃO estuda a troca de peso MLS ou base polinomial.",
            ""]
    (out/"REPORT.md").write_text("\n".join(lines),encoding="utf-8")
    (out/"metadata.json").write_text(json.dumps(dict(
        factors=SUPPORTS,grids=GRIDS,gauss=GAUSS,subdivision=SUBDIVISION,
        MLS_polynomial_degree=1,weight="Wendland_C2",
        interpretation="controlled_support_sensitivity_not_original_source_replication"),
        indent=2)+"\n",encoding="utf-8")
    print("\n".join(lines),flush=True)
    print(f"PASS {len(data)} measured mode rows from 12 independent solver runs",flush=True)

if __name__=="__main__":
    main()
