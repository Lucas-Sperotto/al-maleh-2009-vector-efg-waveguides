"""Regression for reproducible, non-destructive spectral-pollution evidence."""
import importlib.util
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "spectral_pollution", ROOT / "scripts" / "audit_spectral_pollution.py")
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)
LATEST = ROOT / "docs" / "results" / "independent_validation" / "LATEST"


class SpectralPollutionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.campaign = LATEST.parent / LATEST.read_text(encoding="utf-8").strip()

    def test_archived_counts_and_ordering(self):
        for shape, dofs, expected in (
            ("9x7", 126, 49), ("13x11", 286, 94), ("17x15", 510, 144)
        ):
            with self.subTest(shape=shape):
                row = audit.audit_efg_case(self.campaign / "efg" / shape)
                self.assertEqual(row["dofs"], dofs)
                self.assertEqual(row["low_positive_below_tm11"], expected)
                self.assertGreater(row["low_positive_with_divergence_above_025"], 0)
                self.assertLess(row["low_positive_lambda_max"], audit.TM11)

    def test_degenerate_pair_worst_direction_blocks_physical_claim(self):
        # The published pair has a near-solenoidal direction and a polluted
        # direction. A rotation can change *individual column metrics*, but
        # cannot make the supremum over their 2D span disappear.
        physical = audit.audit_degenerate_pair(
            self.campaign / "efg" / "17x15", (192, 193))
        self.assertTrue(physical["subspace_cannot_pass_divergence_threshold"])
        self.assertGreater(
            physical["lower_bound_on_basis_invariant_worst_case_divergence"], 0.41)
        self.assertTrue(physical["exact_subspace_operator_norm_not_computed"])
        self.assertLess(physical["directions"][0]["native_divergence_ratio"], 0.01)
        self.assertGreater(physical["directions"][1]["native_divergence_ratio"], 0.41)

    def test_invalid_or_single_direction_is_rejected(self):
        with self.assertRaises(ValueError):
            audit.audit_degenerate_pair(self.campaign / "efg" / "17x15", (192, 192))
        with self.assertRaises(ValueError):
            audit.audit_degenerate_pair(self.campaign / "efg" / "17x15", (192,))


if __name__ == "__main__":
    unittest.main()
