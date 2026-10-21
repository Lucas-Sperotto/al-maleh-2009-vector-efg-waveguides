"""Regressions grounded in Al Maleh (2009), Ch.5 screenshots Tables 2 and 3."""
import csv
from pathlib import Path
import math
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from m2_ch5_benchmark import reference_rows,EXACT

class OriginalChapterFiveTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows=reference_rows(ROOT/"docs/results/m2_author_ch5_tables_2_3.csv")
        cls.r={(r["mode"],int(r["nodes_y"]),int(r["nodes_x"]),
                int(r["gauss_points_per_axis"])):float(r["author_signed_error_percent"])
               for r in cls.rows}
    def test_complete_original_reported_values_no_invented_blanks(self):
        self.assertEqual(len(self.rows),28)
        self.assertEqual(set(k[0] for k in self.r),{"TM11","TM21"})
        for mode in ("TM11","TM21"):
            for ny,nx,levels in ((3,5,5),(4,7,5),(5,9,2),(6,11,2)):
                for g in range(1,6):
                    self.assertEqual((mode,ny,nx,g) in self.r,g<=levels)
        self.assertEqual(self.r[("TM11",3,5,1)],23.767)
        self.assertEqual(self.r[("TM11",3,5,2)],17.199)
        self.assertEqual(self.r[("TM11",6,11,2)],0.839)
        self.assertEqual(self.r[("TM21",3,5,1)],18.557)
        self.assertEqual(self.r[("TM21",3,5,2)],11.846)
        self.assertEqual(self.r[("TM21",6,11,2)],0.539)
        self.assertAlmostEqual(self.r[("TM11",4,7,5)],3.756)
        self.assertAlmostEqual(self.r[("TM21",4,7,5)],2.161)

    def test_mode_theory_source_equation_5_11(self):
        self.assertAlmostEqual(EXACT["TM11"],1.25*math.pi**2,places=14)
        self.assertAlmostEqual(EXACT["TM21"],2.00*math.pi**2,places=14)
        self.assertAlmostEqual(EXACT["TM11"],12.3370,delta=0.00005)
        self.assertAlmostEqual(EXACT["TM21"],19.7392,delta=0.00005)
        self.assertAlmostEqual(12.1456,12.1456) # Fig.15 is a spurious candidate, not TM11

    def test_quadrature_saturates_between_2_and_5_in_author_data(self):
        for mode in ("TM11","TM21"):
            for ny,nx in ((3,5),(4,7)):
                data=[self.r[(mode,ny,nx,g)] for g in range(2,6)]
                self.assertLess(max(data)-min(data),0.3)
                self.assertGreater(self.r[(mode,ny,nx,1)],self.r[(mode,ny,nx,2)])

    def test_signed_author_values_not_absolute_errors(self):
        for v in self.r.values():
            self.assertGreater(v,0)
        self.assertEqual(self.r[("TM11",4,7,2)],3.608)
        self.assertEqual(self.r[("TM21",4,7,2)],2.038)

if __name__=="__main__":
    unittest.main()
