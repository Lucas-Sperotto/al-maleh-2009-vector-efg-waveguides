"""Archived independent audit campaign (2026-10-22): integrity and findings.

Stdlib only. Reads the versioned CSV outputs; never re-runs or rewrites them.
Each assertion pins a conclusion stated in
docs/results/INDEPENDENT_THEORY_CODE_AUDIT.md, so a future change that
alters the evidence fails loudly instead of silently changing the report.
"""
import csv
import hashlib
import json
import math
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN = ROOT / "docs/results/independent_audit_2026-10-22"
TABLE4 = ROOT / "docs/results/author_ch6_table4_spectra.csv"


def rows(path):
    with Path(path).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


class IndependentAuditCampaign(unittest.TestCase):
    def test_manifest_hashes(self):
        manifest = json.loads((CAMPAIGN / "manifest.json").read_text(encoding="utf-8"))
        listed = {r["path"]: r for r in manifest["files"]}
        actual = {p.relative_to(CAMPAIGN).as_posix() for p in CAMPAIGN.rglob("*")
                  if p.is_file() and p.name != "manifest.json"}
        self.assertEqual(set(listed), actual)
        for path, r in listed.items():
            data = (CAMPAIGN / path).read_bytes()
            self.assertEqual(hashlib.sha256(data).hexdigest(), r["sha256"], path)
            self.assertEqual(len(data), r["bytes"], path)

    def test_table4_transcription(self):
        t = rows(TABLE4)
        self.assertEqual(len(t), 78)
        self.assertEqual(sum(r["method"] == "cartesian" for r in t), 60)
        bold = {(r["method"], r["gauss_points_per_axis"]): r["value_as_printed"]
                for r in t if r["role"] == "TM11_bold"}
        self.assertEqual(bold[("cartesian", "2")], "1.45E+01")
        self.assertEqual([bold[("rot_irrot", g)] for g in "234"],
                         ["1.68E+01", "1.58E+01", "1.55E+01"])

    def test_chapter4_radius_identification(self):
        best = {(r["weight"], r["degree"], r["formulation"]): r
                for r in rows(CAMPAIGN / "ch4_identification_best.csv")}
        for formulation in ("cartesian", "rot_irrot"):
            r = best[("quartic_spline", "1", formulation)]
            self.assertEqual(r["best_factor"], "1.5000")
            # within the 3-significant-digit rounding of the printed table
            self.assertLess(float(r["worst_relative_error_pct"]), 0.5)
        # The project's historical Wendland profile cannot reproduce Table 4.
        self.assertGreater(float(best[("wendland_c2", "1", "cartesian")]["worst_relative_error_pct"]), 10)

    def test_tables_2_3_5_6_reproduction(self):
        t = rows(CAMPAIGN / "tables_2_3_5_6_identified_profile.csv")
        m2 = [r for r in t if r["formulation"] == "M2"]
        m5 = [r for r in t if r["formulation"] == "M5"]
        self.assertEqual(len(m2), 28)
        self.assertTrue(all(r["classification"] == "match_best_curl" for r in m2))
        self.assertLess(max(abs(float(r["best_curl_minus_author_pp"])) for r in m2), 0.001)
        classes = {}
        for r in m5:
            classes[r["classification"]] = classes.get(r["classification"], 0) + 1
        self.assertEqual(classes, {"match_best_curl": 21, "match_nearest_lambda_only": 1,
                                   "magnitude_match_opposite_sign": 1,
                                   "unexplained_difference": 1})
        odd = {(r["mode"], r["nodes_x"], r["gauss"]): r for r in m5
               if r["classification"] != "match_best_curl"}
        self.assertAlmostEqual(float(odd[("TM11", "5", "1")]["nearest_lambda_signed_error_pct"]), -3.572, delta=0.001)
        self.assertAlmostEqual(float(odd[("TM21", "5", "1")]["best_curl_signed_error_pct"]), 1.509, delta=0.001)
        self.assertAlmostEqual(float(odd[("TM11", "7", "1")]["best_curl_minus_author_pp"]), -0.227, delta=0.001)

    def test_m5_reduction_certificates(self):
        conv = rows(CAMPAIGN / "m5_convergence.csv")
        self.assertEqual(len(conv), 240)
        self.assertLess(max(float(r["mass_kernel_stiffness_coupling"]) for r in conv), 1e-12)
        self.assertLess(max(float(r["max_full_residual"]) for r in conv), 1e-12)

    def converged(self, weight, degree, factor):
        out = {}
        for r in rows(CAMPAIGN / "m5_convergence.csv"):
            if (r["weight"], r["degree"], float(r["factor"]), r["subdivision"], r["gauss"],
                    r["mode"]) == (weight, degree, factor, "4", "5", "TM11"):
                out[(int(r["nodes_y"]), int(r["nodes_x"]))] = r
        return out

    def test_author_profile_does_not_converge(self):
        conv = self.converged("quartic_spline", "1", 1.5)
        for (ny, nx), r in conv.items():
            if ny >= 5:
                self.assertGreater(float(r["best_curl_signed_error_pct"]), 20.0)
                self.assertLess(float(r["best_curl_rho"]), 0.91)

    def test_quadratic_basis_converges(self):
        conv = self.converged("wendland_c2", "2", 3.3)
        errors = [float(conv[g]["best_curl_signed_error_pct"])
                  for g in ((5, 9), (7, 13), (9, 17), (13, 25))]
        self.assertTrue(all(a > b for a, b in zip(errors, errors[1:])), errors)
        self.assertLess(errors[-1], 0.1)

    def test_laplacian_defect_is_refinement_independent_for_linear_basis(self):
        lap = {(r["weight"], r["degree"], float(r["factor"]), r["nodes_y"]): r
               for r in rows(CAMPAIGN / "laplacian_consistency.csv")}
        a = float(lap[("quartic_spline", "1", 1.5, "9")]["interior_laplacian_defect_rms"])
        b = float(lap[("quartic_spline", "1", 1.5, "17")]["interior_laplacian_defect_rms"])
        self.assertGreater(a, 1.0)
        self.assertAlmostEqual(a, b, delta=1e-6)
        q = float(lap[("wendland_c2", "2", 3.3, "17")]["laplacian_defect_max"])
        self.assertLess(q, 1e-9)

    def test_rank_threshold_plateau(self):
        for r in rows(CAMPAIGN / "m5_rank_sweep.csv"):
            if r["status"] != "ok":
                continue
            tol = float(r["relative_cutoff"])
            if 1e-14 <= tol <= 1e-6:
                self.assertLess(float(r["max_full_residual"]), 1e-12)
                self.assertLess(float(r["mass_kernel_stiffness_coupling"]), 1e-12)
            if tol == 1e-2:  # over-reduction: projected residual stays blind
                self.assertLess(float(r["max_reduced_residual"]), 1e-12)
                self.assertGreater(float(r["max_full_residual"]), 1e-4)

    def test_doublet_hybridization(self):
        d = {r["kind"]: r for r in rows(CAMPAIGN / "doublet_17x15.csv") if r["kind"] != "eigenvector"}
        eig = {r["indices"]: r for r in rows(CAMPAIGN / "doublet_17x15.csv") if r["kind"] == "eigenvector"}
        self.assertGreater(float(d["computed_doublet_span"]["div_ratio_max"]), 0.4)
        self.assertLess(float(d["l2_projection_of_analytic_doublet"]["div_ratio_max"]), 0.01)
        best = d["best_2d_subspace_of_window"]
        self.assertLess(float(best["div_ratio_max"]), 0.02)
        self.assertGreater(float(best["analytic_cos_2"]), 0.9999)
        self.assertGreater(float(eig["191"]["analytic_cos_1"]), 0.1)        # spurious partner
        self.assertGreater(float(eig["191"]["gradient_fraction_max"]), 0.9)

    def test_m7_tracking(self):
        prov = json.loads((CAMPAIGN / "provenance.json").read_text(encoding="utf-8"))["m7"]
        self.assertEqual(prov["ambiguous_below_30"], 0)
        self.assertGreater(prov["min_tracked_overlap_below_30"], 0.95)
        self.assertEqual(len(prov["best_curl_selection_switches"]), 3)

    def test_archived_m6a_regression(self):
        for r in rows(CAMPAIGN / "archive_regression_m6a.csv"):
            self.assertLess(float(r["max_relative_difference_1_lt_lambda_lt_60"]), 1e-8)
            self.assertEqual(r["archived_kernel"], r["current_kernel"])
            self.assertLess(float(r["current_max_moment_condition_local_basis"]),
                            float(r["archived_max_moment_condition_global_basis"]))
            self.assertTrue(math.isfinite(float(r["max_abs_lambda_difference"])))


if __name__ == "__main__":
    unittest.main()
