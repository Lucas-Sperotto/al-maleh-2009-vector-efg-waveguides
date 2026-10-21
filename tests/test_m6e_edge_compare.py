"""Offline scientific regression tests for independent Whitney TM edge comparison."""
import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest

MODULE=Path(__file__).resolve().parents[1]/"scripts"/"m6e_edge_compare.py"
spec=importlib.util.spec_from_file_location("m6e_edge_compare",MODULE)
m6e=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m6e)

class M6eEdgeComparisonTests(unittest.TestCase):
    def test_geometry_eigenvalue_scaling(self):
        self.assertAlmostEqual(m6e.GEOMETRIC_LAMBDA_MULTIPLIER,0.25)
        self.assertAlmostEqual(m6e.scaled_lambda(4*m6e.ANALYTIC[0]),
                               m6e.ANALYTIC[0],places=12)
        self.assertAlmostEqual(m6e.ANALYTIC[4],m6e.ANALYTIC[5])

    def test_unique_nearest_and_two_dimensional_degeneracy(self):
        modes=[m6e.ANALYTIC[0]*0.02, m6e.ANALYTIC[0]*1.02,
               m6e.ANALYTIC[1]*0.99,
               m6e.ANALYTIC[2]*1.01,m6e.ANALYTIC[3]*1.01,
               m6e.ANALYTIC[4]*0.999,m6e.ANALYTIC[4]*1.001]
        ans=m6e.match_near_analytic(modes)
        self.assertEqual([ans[i][0] for i in range(6)],list(range(1,7)))
        self.assertEqual(len(set(ans[i][0] for i in range(6))),6)

    def test_missing_pair_must_remain_missing(self):
        modes=[m6e.ANALYTIC[i] for i in range(4)]
        ans=m6e.match_near_analytic(modes)
        self.assertIsNone(ans[4])
        self.assertIsNone(ans[5])

    def test_nonfinite_rejected(self):
        for bad in (float("nan"), float("inf"), -float("inf")):
            with self.assertRaises(ValueError):
                m6e.scaled_lambda(bad)

    def test_tm_branch_only_on_native_exports(self):
        with tempfile.TemporaryDirectory() as name:
            d=Path(name); (d/"linop").mkdir(); (d/"csv").mkdir()
            with (d/"linop"/"edge_rect_tm_eigenvalues.csv").open("w") as f:
                f.write("ordered_rank,solver_index,lambda\n")
                f.write("1,0,0\n2,1,49.34802200544679\n")
            with (d/"linop"/"edge_rect_te_eigenvalues.csv").open("w") as f:
                f.write("ordered_rank,solver_index,lambda\n1,0,9\n")
            with (d/"csv"/"edge_rect_modes.csv").open("w") as f:
                f.write("family,rho_abs\nTE,1\nTM,0.999\n")
            (d/"run_timing.csv").write_text("case_label\nedge_rect\n")
            rows,lam,match,timing=m6e.validate_edge_exports(d)
            self.assertEqual(len(rows),2)
            self.assertEqual(lam[1],49.34802200544679)
            self.assertEqual(len([r for r in match if r["family"]=="TM"]),1)
            self.assertEqual(len(timing),1)
            (d/"linop"/"edge_rect_tm_eigenvalues.csv").unlink()
            with self.assertRaises(RuntimeError):
                m6e.validate_edge_exports(d)

if __name__=="__main__":
    unittest.main()
