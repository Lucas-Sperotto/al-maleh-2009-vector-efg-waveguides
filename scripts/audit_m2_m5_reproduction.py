#!/usr/bin/env python3
"""Recompute transparent M2/M5 source-vs-simulation discrepancies (stdlib only).

Inputs are IMMUTABLE source fixtures and the archived 2026-10-21 measured
tables. No solver results are fabricated; all reported percentages are signed.
This tool audits reported statistics, not the physical identity of eigenmodes.
"""
import argparse
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
M2 = ROOT / "docs/results/m2_ch5_observed_2026-10-21.csv"
M5 = ROOT / "docs/results/M5_CH6_MEASURED_2026-10-21.md"
M5_ORIGINAL = ROOT / "docs/results/m5_author_ch6_tables_5_6.csv"
M2_ORIGINAL = ROOT / "docs/results/m2_author_ch5_tables_2_3.csv"

def csv_rows(path):
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))

def read_original(path):
    found = {}
    for row in csv_rows(path):
        key = (row["mode"], int(row["nodes_y"]), int(row["nodes_x"]),
               int(row["gauss_points_per_axis"]))
        if key in found:
            raise ValueError("Repeated author reference "+str(key))
        found[key] = float(row["author_signed_error_percent"])
    return found

def m2():
    result = []
    for r in csv_rows(M2):
        ny, nx = map(int, r["nodes_author_y_x"].split("×"))
        order = int(r["gauss_per_axis"].split("×")[0])
        result.append(dict(mode=r["mode"], ny=ny, nx=nx, order=order,
          author=float(r["author_signed_error_percent"]),
          nearest=float(r["observed_nearest_signed_error_percent"]),
          shape=float(r["observed_best_curl_signed_error_percent"])))
    return result

def m5():
    rows = []
    for line in M5.read_text(encoding="utf-8").splitlines():
        if not line.startswith("| TM11 |") and not line.startswith("| TM21 |"):
            continue
        fields = [v.strip() for v in line.strip().split("|")[1:-1]]
        if len(fields) != 7:
            raise ValueError("M5 measured row format changed")
        ny, nx = map(int, fields[1].split("×"))
        order = int(fields[2].split("×")[0])
        rows.append(dict(mode=fields[0], ny=ny, nx=nx, order=order,
                         author=float(fields[3]), nearest=float(fields[4]),
                         shape=float(fields[5])))
    return rows

def summarize(data):
    def metrics(rows):
        return dict(n=len(rows),
          MAE_nearest_pp=sum(abs(r["nearest"]-r["author"]) for r in rows)/len(rows),
          MAE_field_pp=sum(abs(r["shape"]-r["author"]) for r in rows)/len(rows),
          matching_ambiguous=sum(abs(r["nearest"]-r["shape"])>1 for r in rows),
          within_one_pp=sum(abs(r["shape"]-r["author"])<1 for r in rows))
    return {name:metrics(subset) for name,subset in
       [("all",data),("TM11",[r for r in data if r["mode"]=="TM11"]),
        ("TM21",[r for r in data if r["mode"]=="TM21"]),
        ("3x5",[r for r in data if (r["ny"],r["nx"])==(3,5)]),
        ("larger",[r for r in data if (r["ny"],r["nx"])!=(3,5)])]}

def validate(measured, reference, total):
    if len(measured) != total or len(reference) != total:
        raise ValueError("Expected exactly %d measured and source rows"%total)
    seen=set()
    for r in measured:
        key=(r["mode"],r["ny"],r["nx"],r["order"])
        if key in seen or key not in reference:
            raise ValueError("Duplicated or ungrounded measured comparison "+str(key))
        seen.add(key)
        if not all(math.isfinite(r[k]) for k in ("author","nearest","shape")):
            raise ValueError("Nonfinite measured outcome")
        if abs(reference[key]-r["author"]) > 0.000001:
            raise ValueError("Original-author fixture disagrees with measured table")
    if seen!=set(reference):
        raise ValueError("Missing comparisons")

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path)
    args=parser.parse_args()
    two=m2();five=m5()
    validate(two,read_original(M2_ORIGINAL),28)
    validate(five,read_original(M5_ORIGINAL),24)
    report=dict(
      evidence="Archived source tables and separately recorded measured results",
      precision="Percent errors rounded to 0.001 percentage point in archived reports",
      interpretation="Descriptive reproduction gaps; NOT proof of physical eigenmode identity",
      M2=summarize(two),M5=summarize(five),
      rank_observed=[dict(nodes=n,full_dofs=2*n,mass_rank=2*n-4,
                          mass_nullity=4,near_zero_curl=n-1)
                     for n in (15,28,45,66,91)],
      rank_scope="Observed for these five linear-MLS clouds; no general-rank theorem established")
    text=json.dumps(report,indent=2,ensure_ascii=False)+"\n"
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        args.out.write_text(text,encoding="utf-8")
    print(text,end="")

if __name__=="__main__":
    main()
