#!/usr/bin/env python3
"""Independent theory/code audit campaign (2026-10-22), fully reproducible.

Runs NEW experiments only; historical M2/M4/M5/M6 data are read, never
rewritten. Stdlib only. Every number in the generated REPORT.md comes from
CSV files written by the compiled executables of this repository.

Studies:
  1. Inverse identification of the Chapter-4 MLS radius against Table 4
     (quartic spline Eq. 4-2 is from the thesis; the radius is NOT stated
     there and is inferred here, explicitly labelled as an inference).
  2. Tables 2/3/5/6 with the identified profile (M2 and M5 drivers).
  3. Laplacian (second-order) consistency defect of the MLS shapes.
  4. M5 refinement with production and finer reference quadrature.
  5. M5 singular-mass rank-threshold sweep with full/projected residuals.
  6. Cartesian TM22/TM41 doublet: invariant subspace analysis.
  7. M7 mode tracking of the M5 TM21 3x5 support-factor sweep.
  8. Regression of archived M6a spectra after the MLS basis change.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TABLE4 = ROOT / "docs/results/author_ch6_table4_spectra.csv"
M6A = ROOT / "docs/results/independent_validation/2026-10-21_dc1d98b"
THESIS_SHA256 = "ba356899e516894301a907f5ead06392a2132dd15e5c7b3811733dbae06cc506"
PI = math.pi


def run(cmd, log=None, timeout=3600):
    print("+", " ".join(str(c) for c in cmd), flush=True)
    p = subprocess.run([str(c) for c in cmd], capture_output=True, text=True,
                       timeout=timeout, check=False)
    if log is not None:  # keep archived logs free of machine-specific paths
        log.write_text((p.stdout + p.stderr).replace(str(ROOT), "<repo>"), encoding="utf-8")
    if p.returncode:
        raise RuntimeError(f"Command failed ({p.returncode}): {cmd}\n{p.stderr[-2000:]}")
    return p.stdout


def read_csv(path):
    with Path(path).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, fields=None):
    if not rows:
        raise ValueError(f"No rows for {path}")
    with Path(path).open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields or list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def gzip_file(src, dst):
    with Path(src).open("rb") as fi, gzip.GzipFile(dst, "wb", mtime=0) as fo:
        shutil.copyfileobj(fi, fo)
    Path(src).unlink()


def git(*args):
    try:
        return subprocess.check_output(["git", "-C", str(ROOT), *args], text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


# ------------------------------------------------------------------ study 1
def author_table4():
    rows = read_csv(TABLE4)
    out = {}
    for r in rows:
        if r["value"] == "":
            continue
        out.setdefault((r["method"], int(r["gauss_points_per_axis"])), []).append(
            (int(r["row"]), float(r["value"])))
    for key in out:
        out[key].sort()
    if len(out[("cartesian", 2)]) != 20 or len(out[("rot_irrot", 2)]) != 5:
        raise ValueError("Unexpected Table 4 transcription size")
    return out


def score_scan(scan_rows, table):
    """Mean/max |log10(ours/author)| for Table-4 nonzero entries."""
    grouped = {}
    for r in scan_rows:
        if r["status"] != "ok":
            continue
        key = (r["weight"], int(r["degree"]), float(r["factor"]), r["formulation"],
               int(r["gauss"]))
        grouped.setdefault(key, {})[int(r["position"])] = float(r["lambda"])
    summary = []
    for (weight, degree, factor, formulation, gauss), values in sorted(grouped.items()):
        if formulation == "cartesian":
            ref = [(row, v) for row, v in table[("cartesian", gauss)] if row > 5]
            ours = [values.get(row) for row, _ in ref]          # rows 6..20
        else:
            # Table 4 row 1 is "14 zeros"; align the first nonzero value of
            # OUR reduced spectrum with row 2 instead of assuming 14 zeros.
            zeros = sum(abs(v) < 1e-8 for v in values.values())
            ref = table[("rot_irrot", gauss)]                  # rows 2..6
            ours = [values.get(zeros + row - 1) for row, _ in ref]
        # A nonpositive/missing value at a Table-4 position means the
        # positional comparison is meaningless for that profile: flag it
        # explicitly instead of writing a non-finite score.
        comparable = all(o is not None and o > 0 for o in ours)
        errs = [abs(math.log10(o / v)) for o, (_, v) in zip(ours, ref)] if comparable else []
        summary.append(dict(weight=weight, degree=degree, factor=f"{factor:.4f}",
                            formulation=formulation, gauss=gauss, compared=len(ref),
                            comparable=comparable,
                            mean_abs_log10_ratio=sum(errs) / len(errs) if comparable else "",
                            max_abs_log10_ratio=max(errs) if comparable else ""))
    return summary


def within_printed_rounding(value, printed):
    """Half a unit of the last printed decimal, including scientific notation."""
    reference = Decimal(printed)
    half_unit = Decimal(10) ** reference.as_tuple().exponent / 2
    return abs(Decimal(str(value)) - reference) <= half_unit


def table4_precision(out):
    """Check each nonzero Table-4 entry at the inferred profile, without fitting."""
    with gzip.open(out / "scan_quartic_spline_deg1.csv.gz", "rt", encoding="utf-8") as f:
        scan = [r for r in csv.DictReader(f)
                if r["status"] == "ok" and float(r["factor"]) == 1.5]
    rows = []
    for formulation in ("cartesian", "rot_irrot"):
        for gauss in (2, 3, 4):
            values = {int(r["position"]): float(r["lambda"]) for r in scan
                      if r["formulation"] == formulation and int(r["gauss"]) == gauss}
            zeros = sum(abs(v) < 1e-8 for v in values.values())
            for r in read_csv(TABLE4):
                if (r["method"] != formulation or int(r["gauss_points_per_axis"]) != gauss
                        or not r["value"] or (formulation == "cartesian" and int(r["row"]) <= 5)):
                    continue
                position = int(r["row"]) if formulation == "cartesian" else zeros + int(r["row"]) - 1
                value = values[position]
                printed = r["value_as_printed"]
                rows.append(dict(formulation=formulation, gauss=gauss, author_row=r["row"],
                                 computed_position=position, value_as_printed=printed,
                                 computed_lambda=value,
                                 within_printed_rounding=within_printed_rounding(value, printed)))
    write_csv(out / "table4_printed_precision.csv", rows)
    return rows


def study_ch4(exe, out):
    table = author_table4()
    scans = [("quartic_spline", 1, 1.20, 3.50), ("wendland_c2", 1, 1.20, 3.50),
             ("quartic_spline", 2, 2.20, 3.50)]
    summary = []
    for weight, degree, f0, f1 in scans:
        raw = out / f"scan_{weight}_deg{degree}.csv"
        run([exe, "ch4-scan", weight, degree, f"{f0:.2f}", f"{f1:.2f}", "0.01", raw])
        summary += score_scan(read_csv(raw), table)
        gzip_file(raw, raw.with_suffix(".csv.gz"))
    write_csv(out / "ch4_scan_vs_table4.csv", summary)
    best = []
    keys = sorted({(r["weight"], r["degree"], r["formulation"]) for r in summary})
    for weight, degree, formulation in keys:
        per_factor = {}
        for r in summary:
            if (r["weight"], r["degree"], r["formulation"]) == (weight, degree, formulation):
                per_factor.setdefault(r["factor"], []).append(r)
        complete = {f: rs for f, rs in per_factor.items()
                    if len(rs) == 3 and all(x["comparable"] for x in rs)}
        if not complete:
            best.append(dict(weight=weight, degree=degree, formulation=formulation,
                             best_factor="none", mean_abs_log10_ratio_gauss2="",
                             mean_abs_log10_ratio_gauss3="", mean_abs_log10_ratio_gauss4="",
                             worst_abs_log10_ratio="", worst_relative_error_pct=""))
            continue
        f_best = min(complete, key=lambda f: sum(x["mean_abs_log10_ratio"] for x in complete[f]))
        rs = complete[f_best]
        best.append(dict(weight=weight, degree=degree, formulation=formulation,
                         best_factor=f_best,
                         mean_abs_log10_ratio_gauss2=rs[0]["mean_abs_log10_ratio"],
                         mean_abs_log10_ratio_gauss3=rs[1]["mean_abs_log10_ratio"],
                         mean_abs_log10_ratio_gauss4=rs[2]["mean_abs_log10_ratio"],
                         worst_abs_log10_ratio=max(x["max_abs_log10_ratio"] for x in rs),
                         worst_relative_error_pct=100 * (10 ** max(x["max_abs_log10_ratio"] for x in rs) - 1)))
    write_csv(out / "ch4_identification_best.csv", best)
    table4_precision(out)
    return best


# ------------------------------------------------------------------ study 2
def classify(author, curl, nearest):
    if abs(curl - author) <= 0.01:
        return "match_best_curl"
    if abs(nearest - author) <= 0.01:
        return "match_nearest_lambda_only"
    if abs(abs(curl) - abs(author)) <= 0.01 and curl * author < 0:
        return "magnitude_match_opposite_sign"
    return "unexplained_difference"


def prune_case_outputs(folder, keep):
    for item in sorted(folder.rglob("*"), reverse=True):
        if item.is_file() and item.name not in keep:
            item.unlink()


def study_tables(build, out):
    m2 = out / "m2_identified_profile"
    m5 = out / "m5_identified_profile"
    env = ["env", "PYTHONNOUSERSITE=1", "OPENBLAS_NUM_THREADS=1"]
    run(env + [sys.executable, ROOT / "scripts/m2_ch5_benchmark.py", "--binary",
               build / "tm_validation", "--out", m2, "--support-factor", "1.5",
               "--weight", "quartic_spline"], log=out / "m2_identified_profile.log")
    run(env + [sys.executable, ROOT / "scripts/m5_ch6_campaign.py", "--binary",
               build / "rot_irrot_rect", "--out", m5, "--support-factor", "1.5",
               "--weight", "quartic_spline"], log=out / "m5_identified_profile.log")
    rows = []
    for name, folder, curl_key in (("M2", m2, "current_strongest_curl_signed_error_percent"),
                                   ("M5", m5, "current_best_curl_signed_error_percent")):
        reference = ROOT / "docs/results" / ("m2_author_ch5_tables_2_3.csv" if name == "M2"
                                             else "m5_author_ch6_tables_5_6.csv")
        printed_errors = {(r["mode"], r["nodes_y"], r["nodes_x"], r["gauss_points_per_axis"]):
                          r["author_signed_error_percent"] for r in read_csv(reference)}
        for r in read_csv(folder / "comparison.csv"):
            printed_error = printed_errors[(r["mode"], r["nodes_y"], r["nodes_x"], r["gauss_points_per_axis"])]
            author = float(r["author_signed_error_percent"])
            if author != float(printed_error):
                raise ValueError("Campaign changed a published reference value")
            curl = float(r[curl_key])
            nearest = float(r["current_nearest_signed_error_percent"])
            rows.append(dict(formulation=name, mode=r["mode"], nodes_y=r["nodes_y"],
                             nodes_x=r["nodes_x"], gauss=r["gauss_points_per_axis"],
                             author_signed_error_pct=author,
                             best_curl_signed_error_pct=curl,
                             nearest_lambda_signed_error_pct=nearest,
                             best_curl_minus_author_pp=curl - author,
                             nearest_minus_author_pp=nearest - author,
                             best_curl_within_printed_rounding=within_printed_rounding(
                                 curl, printed_error),
                             classification=classify(author, curl, nearest)))
    write_csv(out / "tables_2_3_5_6_identified_profile.csv", rows)
    # Keep compact, reviewable outputs; matrices and fields stay reproducible.
    keep = {"comparison.csv", "metadata.json", "REPORT.md", "spectrum_summary.csv",
            "eigenvalues.csv", "tm_matches.csv", "run_metadata.txt",
            "candidate_spectrum.csv"}
    for folder in (m2, m5):
        prune_case_outputs(folder, keep)
    return rows


# --------------------------------------------------------------- studies 3-6
def study_numerics(exe, out):
    run([exe, "laplacian", out / "laplacian_consistency.csv"])
    run([exe, "convergence", out / "m5_convergence.csv"], timeout=7200)
    run([exe, "rank-sweep", out / "m5_rank_sweep.csv"])
    for nx, ny in ((9, 7), (13, 11), (17, 15)):
        run([exe, "doublet", nx, ny, out / f"doublet_{nx}x{ny}.csv"])


# ------------------------------------------------------------------ study 7
def study_tracking(build, out):
    factors = [f"{3.0 + 0.02 * k:.2f}" for k in range(31)]
    folder = out / "m7_tm21_3x5_support_sweep"
    run([build / "mode_tracking", "m5", 5, 3, 2, 2, 1, "wendland_c2", 40, folder, *factors])
    sel = [r for r in read_csv(folder / "selection.csv") if r["target"] == "TM21"]
    trk = read_csv(folder / "tracking.csv")
    low = [r for r in trk if float(r["lambda_from"]) < 30]
    switches = [f"{float(r['factor']):.2f}" for r in sel
                if r["best_curl_from_tracked_previous_best"] == "0"]
    return dict(factors=len(factors), transitions_below_30=len(low),
                ambiguous_below_30=sum(r["status"] != "tracked" for r in low),
                min_tracked_overlap_below_30=min(float(r["h_overlap"]) for r in low),
                best_curl_selection_switches=switches,
                tm21_best_curl_error_range_pct=(min(float(r["best_curl_signed_error_pct"]) for r in sel),
                                                max(float(r["best_curl_signed_error_pct"]) for r in sel)),
                tm21_best_curl_rho_range=(min(float(r["best_curl_rho"]) for r in sel),
                                          max(float(r["best_curl_rho"]) for r in sel)))


# ------------------------------------------------------------------ study 8
def study_archive_regression(build, out):
    rows = []
    for nx, ny in ((9, 7), (13, 11), (17, 15)):
        folder = out / f"archive_regression_{nx}x{ny}"
        run([build / "tm_validation", nx, ny, 1, 3.3, 3, folder, 1])
        old = [float(r["lambda"]) for r in read_csv(M6A / f"efg/{nx}x{ny}/eigenvalues.csv")]
        new = [float(r["lambda"]) for r in read_csv(folder / "eigenvalues.csv")]
        if len(old) != len(new):
            raise RuntimeError("Archived and recomputed spectra differ in size")
        phys = [(a, b) for a, b in zip(old, new) if 1.0 < a < 60.0]
        meta_old = json.loads((M6A / f"efg/{nx}x{ny}/run_metadata.json").read_text())
        meta_new = json.loads((folder / "run_metadata.json").read_text())
        rows.append(dict(cloud=f"{nx}x{ny}", modes=len(old),
                         max_abs_lambda_difference=max(abs(a - b) for a, b in zip(old, new)),
                         max_relative_difference_1_lt_lambda_lt_60=max(abs(a - b) / a for a, b in phys),
                         archived_max_moment_condition_global_basis=meta_old["quadrature"]["max_moment_condition"],
                         current_max_moment_condition_local_basis=meta_new["quadrature"]["max_moment_condition"],
                         archived_kernel=meta_old["spectrum"]["numerical_kernel_dimension"],
                         current_kernel=meta_new["spectrum"]["numerical_kernel_dimension"]))
        shutil.rmtree(folder)
    write_csv(out / "archive_regression_m6a.csv", rows)
    return rows


# ------------------------------------------------------------------- report
def fmt(x, nd=3):
    return f"{float(x):+.{nd}f}"


def report(out, best, tables, track, regression):
    conv = read_csv(out / "m5_convergence.csv")
    lap = read_csv(out / "laplacian_consistency.csv")
    sweep = read_csv(out / "m5_rank_sweep.csv")
    doubl = {g: read_csv(out / f"doublet_{g}.csv") for g in ("9x7", "13x11", "17x15")}
    counts = {}
    for r in tables:
        counts.setdefault(r["formulation"], {}).setdefault(r["classification"], 0)
        counts[r["formulation"]][r["classification"]] += 1
    L = ["# Auditoria independente teoria × código — campanha 2026-10-22", "",
         "Campanha NOVA. Nenhum dado histórico foi regravado. Números de Al Maleh",
         "vêm das Tabelas 2–6 da dissertação (PDF McGill, SHA-256",
         f"`{THESIS_SHA256}`); números do programa vêm dos CSV desta pasta.", "",
         "## 1. Identificação inversa do perfil MLS do Capítulo 4 (inferência)", "",
         "A Eq. (4-2) fixa a spline quártica; o valor de r NÃO consta na tese.",
         "Varredura do fator r/h contra os 15 autovalores não nulos cartesianos e",
         "os 5 rotacionais/irrotacionais da Tabela 4 (malha 3×5, Gauss 2–4):", "",
         "| Peso | Grau | Formulação | Melhor r/h | média |log10| G2/G3/G4 | pior erro relativo |",
         "|---|---:|---|---:|---|---:|"]
    for b in best:
        if b["best_factor"] == "none":
            L.append(f"| {b['weight']} | {b['degree']} | {b['formulation']} | — | "
                     "sem fator comparável (autovalores não positivos nas posições da Tabela 4) | — |")
            continue
        L.append(f"| {b['weight']} | {b['degree']} | {b['formulation']} | {b['best_factor']} | "
                 f"{b['mean_abs_log10_ratio_gauss2']:.4f}/{b['mean_abs_log10_ratio_gauss3']:.4f}/"
                 f"{b['mean_abs_log10_ratio_gauss4']:.4f} | "
                 + (f"{b['worst_relative_error_pct']:.2f}% |" if b['worst_relative_error_pct'] < 1000
                    else ">1000% (sem correspondência) |"))
    precision = read_csv(out / "table4_printed_precision.csv")
    L += ["", "A conferência usa meia unidade da última casa de cada valor impresso, não um corte relativo uniforme.",
          f"Tabela 4: {sum(r['within_printed_rounding'] == 'True' for r in precision)}/{len(precision)} valores não nulos dentro do intervalo de arredondamento; ver table4_printed_precision.csv.",
          "As entradas fora do intervalo são listadas no CSV: a reprodução pode ser apenas aproximada.",
          "Melhor r/h no limite da faixa varrida indica ausência de mínimo interno, não identificação.",
          "", "## 2. Tabelas 2, 3, 5 e 6 com o perfil identificado", "",
          "| Formulação | Classificação | Linhas |", "|---|---|---:|"]
    for f in sorted(counts):
        for c in sorted(counts[f]):
            L.append(f"| {f} | {c} | {counts[f][c]} |")
    L += ["", "Classificação match: tolerância operacional de 0,01 ponto percentual; não significa igualdade na precisão impressa (0,001 pp).",
          "O perfil r/h=1,5 e grau 1 foi fixado pela Tabela 4 antes das comparações das Tabelas 2/3/5/6; não há novo ajuste nessas tabelas.",
          "F6 permanece hipótese de sinal; F8 permanece discrepância não explicada.",
          "Arquivo compacto: espectros, correlações, CSV de diagnósticos e logs; matrizes, coeficientes e campos completos dos drivers M2/M5 são regeneráveis e foram removidos desta campanha compacta.",
          "Os relatórios internos dos drivers descrevem a saída completa anterior à compactação; o manifesto de nível superior inventaria os arquivos realmente presentes."]
    L += ["", "| Form. | Modo | Malha | Gauss | Autor % | Melhor curl % | λ próximo % | Δ curl (pp) | Classe |",
          "|---|---|---|---:|---:|---:|---:|---:|---|"]
    for r in tables:
        L.append(f"| {r['formulation']} | {r['mode']} | {r['nodes_y']}×{r['nodes_x']} | {r['gauss']} | "
                 f"{fmt(r['author_signed_error_pct'])} | {fmt(r['best_curl_signed_error_pct'],4)} | "
                 f"{fmt(r['nearest_lambda_signed_error_pct'],4)} | {fmt(r['best_curl_minus_author_pp'],4)} | "
                 f"{r['classification']} |")
    L += ["", "## 3. Consistência de 2ª ordem (defeito de Δ para x²+y², valor exato 4)", "",
          "| Peso | Grau | r/h | Malha | cond. máx. (base local) | defeito máx. | RMS interior |",
          "|---|---:|---:|---|---:|---:|---:|"]
    for r in lap:
        if r["status"] == "ok":
            interior = f"{float(r['interior_laplacian_defect_rms']):.3e}" if r["interior_laplacian_defect_rms"] else "—"
            L.append(f"| {r['weight']} | {r['degree']} | {float(r['factor']):.2f} | {r['nodes_y']}×{r['nodes_x']} | "
                     f"{float(r['max_moment_condition']):.1f} | {float(r['laplacian_defect_max']):.3e} | {interior} |")
    L += ["", "## 4. Convergência M5 (TM11, melhor curl; quadratura de produção G2 e referência 16 subcélulas×G5)", "",
          "| Peso | Grau | r/h | Malha | erro G2 (%) | erro referência (%) | ρ curl | ‖KZ‖/‖K‖ | resíduo completo |",
          "|---|---:|---:|---|---:|---:|---:|---:|---:|"]
    index = {(r["weight"], r["degree"], r["factor"], r["nodes_y"], r["nodes_x"], r["subdivision"], r["gauss"], r["mode"]): r
             for r in conv}
    for r in conv:
        if r["mode"] != "TM11" or (r["subdivision"], r["gauss"]) != ("2", "2"):
            continue
        c = index[(r["weight"], r["degree"], r["factor"], r["nodes_y"], r["nodes_x"], "4", "5", "TM11")]
        L.append(f"| {r['weight']} | {r['degree']} | {float(r['factor']):.1f} | {r['nodes_y']}×{r['nodes_x']} | "
                 f"{fmt(r['best_curl_signed_error_pct'])} | {fmt(c['best_curl_signed_error_pct'])} | "
                 f"{float(c['best_curl_rho']):.4f} | {float(c['mass_kernel_stiffness_coupling']):.1e} | "
                 f"{float(c['max_full_residual']):.1e} |")
    L += ["", "## 5. Limiar de posto da massa M5", "",
          "| Caso | limiar | posto | nulidade | resíduo projetado | resíduo completo | ‖KZ‖/‖K‖ |",
          "|---|---:|---:|---:|---:|---:|---:|"]
    for r in sweep:
        if r["status"] != "ok":
            L.append(f"| {r['weight']} {r['degree']} {float(r['factor']):.1f} {r['nodes_y']}×{r['nodes_x']} G{r['gauss']} | "
                     f"{float(r['relative_cutoff']):.0e} | falhou | — | — | — | — |")
            continue
        L.append(f"| {r['weight']} {r['degree']} {float(r['factor']):.1f} {r['nodes_y']}×{r['nodes_x']} G{r['gauss']} | "
                 f"{float(r['relative_cutoff']):.0e} | {r['retained_rank']} | {r['mass_nullity']} | "
                 f"{float(r['max_reduced_residual']):.1e} | {float(r['max_full_residual']):.1e} | "
                 f"{float(r['mass_kernel_stiffness_coupling']):.1e} |")
    L += ["", "## 6. Dubleto TM22/TM41 cartesiano (grandezas invariantes de subespaço)", "",
          "Projeções L² e combinações da janela com autovalores distintos são aproximações geométricas, não novos autopares exatos. A hibridização é uma hipótese apoiada pelas correlações, sem prova causal.", "",
          "| Nuvem | Objeto | Índices | cos. principais c/ analítico | div mín–máx | fração gradiente máx |",
          "|---|---|---|---|---|---:|"]
    for g, rows in doubl.items():
        for r in rows:
            L.append(f"| {g} | {r['kind']} | {r['indices']} | {float(r['analytic_cos_1']):.5f}, "
                     f"{float(r['analytic_cos_2']):.5f} | {float(r['div_ratio_min']):.4f}–"
                     f"{float(r['div_ratio_max']):.4f} | {float(r['gradient_fraction_max']):.4f} |")
    L += ["", "## 7. Rastreamento M7 (M5, 3×5, Gauss 2×2, Wendland, r/h 3,00→3,60)", "",
          f"- Fatores: {track['factors']}; transições com λ<30: {track['transitions_below_30']};",
          f"  ambíguas: {track['ambiguous_below_30']}; menor sobreposição L² rastreada: "
          f"{track['min_tracked_overlap_below_30']:.4f}.",
          f"- Trocas do candidato 'melhor curl' TM21 (não sucessor rastreado): {', '.join(track['best_curl_selection_switches'])}.",
          f"- Erro TM21 do candidato melhor curl varia de {track['tm21_best_curl_error_range_pct'][0]:+.3f}% a "
          f"{track['tm21_best_curl_error_range_pct'][1]:+.3f}%; ρ entre {track['tm21_best_curl_rho_range'][0]:.4f} e "
          f"{track['tm21_best_curl_rho_range'][1]:.4f}.", "",
          "## 8. Regressão dos espectros arquivados M6a após a mudança de base MLS", "",
          "| Nuvem | modos | máx |Δλ| | máx Δλ/λ (1<λ<60) | cond. MLS arquivado (base global) | atual (base local) | núcleo arq./atual |",
          "|---|---:|---:|---:|---:|---:|---|"]
    for r in regression:
        L.append(f"| {r['cloud']} | {r['modes']} | {float(r['max_abs_lambda_difference']):.2e} | "
                 f"{float(r['max_relative_difference_1_lt_lambda_lt_60']):.2e} | "
                 f"{float(r['archived_max_moment_condition_global_basis']):.1f} | "
                 f"{float(r['current_max_moment_condition_local_basis']):.1f} | "
                 f"{r['archived_kernel']}/{r['current_kernel']} |")
    L += ["", "Interpretação científica: ver docs/results/INDEPENDENT_THEORY_CODE_AUDIT.md.", ""]
    (out / "REPORT.md").write_text("\n".join(L), encoding="utf-8")


def inventory(root):
    return [{"path": p.relative_to(root).as_posix(), "bytes": p.stat().st_size,
             "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
            for p in sorted(root.rglob("*")) if p.is_file() and p != root / "manifest.json"]


def prepare_output(out):
    """Create an empty campaign directory, preserving every existing file."""
    if out.exists():
        if not out.is_dir():
            raise ValueError(f"Output path {out} is not a directory; never overwrite a file")
        if any(out.iterdir()):
            raise ValueError(f"Output directory {out} is not empty; never overwrite a campaign")
    out.mkdir(parents=True, exist_ok=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build-dir", type=Path, default=ROOT / "build")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    build = a.build_dir.resolve()
    exe = build / "independent_audit"
    for tool in ("independent_audit", "tm_validation", "rot_irrot_rect", "mode_tracking"):
        if not (build / tool).is_file():
            ap.error(f"Missing executable {build / tool}; build the project first")
    out = a.out.resolve()
    try:
        prepare_output(out)
    except (ValueError, OSError) as error:
        ap.error(str(error))
    best = study_ch4(exe, out)
    tables = study_tables(build, out)
    study_numerics(exe, out)
    track = study_tracking(build, out)
    regression = study_archive_regression(build, out)
    provenance = dict(
        campaign="independent_theory_code_audit_2026-10-22",
        created_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        base_commit=git("rev-parse", "HEAD"), branch=git("branch", "--show-current"),
        working_tree_dirty=bool(git("status", "--porcelain")),
        tracked_diff_sha256=hashlib.sha256(git("diff", "HEAD", "--").encode()).hexdigest(),
        executable_sha256={name: hashlib.sha256((build / name).read_bytes()).hexdigest()
                           for name in ("independent_audit", "tm_validation", "rot_irrot_rect", "mode_tracking")},
        reference_csv_sha256={path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                              for path in (TABLE4, ROOT / "docs/results/m2_author_ch5_tables_2_3.csv",
                                           ROOT / "docs/results/m5_author_ch6_tables_5_6.csv")},
        python=sys.version.split()[0], platform=platform.platform(),
        thesis_pdf="https://mcgill.scholaris.ca/bitstreams/4d79b125-540f-4b38-8c06-722d14eb60bb/download",
        thesis_pdf_sha256=THESIS_SHA256,
        identified_profile=dict(weight="quartic_spline (Eq. 4-2)", degree=1,
                                radius_over_spacing="1.5 (inferred, not stated in thesis)"),
        m7=track)
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2, ensure_ascii=False) + "\n",
                                          encoding="utf-8")
    report(out, best, tables, track, regression)
    manifest = dict(schema_version=1, hash_policy="SHA-256 of every campaign file except manifest.json itself",
                    files=inventory(out))
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
                                        encoding="utf-8")
    print(f"PASS: audit campaign written to {out} ({len(manifest['files'])} files)")


if __name__ == "__main__":
    main()
