"""Audits signed thesis reference errors separately from observed errors."""
import importlib.util
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
source=ROOT/"scripts/audit_m2_m5_reproduction.py"
spec=importlib.util.spec_from_file_location("audit_m2m5",source)
audit=importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)

class ComparativeAuditTests(unittest.TestCase):
    def test_complete_source_citations_and_no_missing_rows(self):
        a=audit.m2(); b=audit.m5()
        audit.validate(a,audit.read_original(audit.M2_ORIGINAL),28)
        audit.validate(b,audit.read_original(audit.M5_ORIGINAL),24)
        self.assertEqual({r["mode"] for r in a},{"TM11","TM21"})
        self.assertEqual({r["mode"] for r in b},{"TM11","TM21"})
    def test_objective_signed_metrics(self):
        a=audit.summarize(audit.m2())
        b=audit.summarize(audit.m5())
        self.assertAlmostEqual(a["all"]["MAE_nearest_pp"],4.997964285714284,places=11)
        self.assertAlmostEqual(a["all"]["MAE_field_pp"],1.8594642857142858,places=11)
        self.assertAlmostEqual(b["all"]["MAE_nearest_pp"],20.76245833333334,places=11)
        self.assertAlmostEqual(b["all"]["MAE_field_pp"],18.138791666666666,places=11)
        self.assertEqual(a["all"]["matching_ambiguous"],5)
        self.assertEqual(b["all"]["matching_ambiguous"],5)
        self.assertEqual(a["all"]["within_one_pp"],11)
        self.assertEqual(b["all"]["within_one_pp"],0)
        # Coarse-grid agreement does NOT imply convergence agreement.
        self.assertGreater(a["larger"]["MAE_field_pp"],a["3x5"]["MAE_field_pp"])
        self.assertGreater(b["larger"]["MAE_field_pp"],b["3x5"]["MAE_field_pp"])
    def test_empirical_rank_table_separate_mass_and_curl_nullspaces(self):
        path=ROOT/"docs/results/M5_CH6_MEASURED_2026-10-21.md"
        rank_section=path.read_text(encoding="utf-8").split("## Posto e quase-núcleo")[1].split("## Comparações")[0]
        obs=[]
        for line in rank_section.splitlines():
            if not line.startswith("| ") or "×" not in line:
                continue
            p=[x.strip() for x in line.split("|")[1:-1]]
            if len(p)!=6 or not p[1].isdigit():
                continue
            ny,nx=map(int,p[0].split("×"))
            n=ny*nx
            rank=int(p[2]); mass_null=int(p[3]);curl_null=int(p[4])
            self.assertEqual((rank,mass_null,curl_null),(2*n-4,4,n-1))
            obs.append((n,p[1]))
        self.assertEqual(len(obs),12)
        self.assertEqual({n for n,_ in obs},{15,28,45,66,91})

if __name__=="__main__": unittest.main()
