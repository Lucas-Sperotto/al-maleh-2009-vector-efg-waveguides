"""Exact fixtures from user-provided Al Maleh (2009), Ch6 Tables 5/6."""
import csv
import math
from pathlib import Path
import sys
import unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from m5_ch6_campaign import read_reference,CASES,EXACT

class OriginalChapterSixTests(unittest.TestCase):
    def test_source_rows_and_only_reported_cells(self):
        rows=read_reference()
        self.assertEqual(len(rows),24)
        actual={(r["mode"],int(r["nodes_x"]),int(r["nodes_y"]),
                 int(r["gauss_points_per_axis"])):float(r["author_signed_error_percent"])
                 for r in rows}
        self.assertEqual(len(actual),24)
        for mode in ("TM11","TM21"):
            for (nx,ny),gausses in CASES.items():
                for order in range(1,5):
                    self.assertEqual((mode,nx,ny,order) in actual,
                                     order in gausses)
        self.assertAlmostEqual(actual[("TM11",5,3,1)],-3.572)
        self.assertAlmostEqual(actual[("TM11",5,3,2)],+36.183)
        self.assertAlmostEqual(actual[("TM11",5,3,3)],+28.258)
        self.assertAlmostEqual(actual[("TM11",5,3,4)],+26.008)
        self.assertAlmostEqual(actual[("TM11",7,4,2)],+34.120)
        self.assertAlmostEqual(actual[("TM11",13,7,1)],-13.257)
        self.assertAlmostEqual(actual[("TM21",5,3,1)],-1.509)
        self.assertAlmostEqual(actual[("TM21",5,3,2)],+38.224)
        self.assertAlmostEqual(actual[("TM21",7,4,4)],+23.227)
        self.assertAlmostEqual(actual[("TM21",13,7,1)],-12.831)
    def test_no_false_monotone_convergence_claim(self):
        rows=read_reference()
        selected={(r["mode"],int(r["nodes_x"]),int(r["nodes_y"]),
                   int(r["gauss_points_per_axis"])):float(r["author_signed_error_percent"])
                  for r in rows}
        # Increasing 1x1 to 2x2 actually increases the reported signed
        # error dramatically for this rotational/irrotational method.
        for mode in ("TM11","TM21"):
            self.assertLess(selected[mode,5,3,1],0)
            self.assertGreater(selected[mode,5,3,2],20)
        self.assertGreater(selected["TM11",9,5,2],
                           selected["TM11",7,4,2])
    def test_analytic_tm_guide(self):
        self.assertAlmostEqual(EXACT["TM11"],12.337005501,places=8)
        self.assertAlmostEqual(EXACT["TM21"],19.739208802,places=8)

if __name__=="__main__":
    unittest.main()
