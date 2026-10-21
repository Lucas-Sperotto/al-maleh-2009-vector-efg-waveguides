"""Integrity checks must detect changed and omitted scientific artifacts."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location("manifest", Path(__file__).resolve().parents[1] /
                                                "scripts/validation_manifest.py")
manifest = importlib.util.module_from_spec(spec)
spec.loader.exec_module(manifest)


class ManifestTests(unittest.TestCase):
    def test_tampering_and_missing_file(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            artifact = root / "spectrum.csv"
            artifact.write_text("index,lambda\n0,0\n1,12.337005501361698\n")
            original = artifact.read_text()
            (root / "manifest.json").write_text(json.dumps({"files": manifest.inventory(root)}))
            self.assertEqual(manifest.verify(root), 1)
            artifact.write_text(original.replace("12.337", "13.337"))
            with self.assertRaisesRegex(ValueError, "spectrum.csv"):
                manifest.verify(root)
            artifact.write_text(original)
            self.assertEqual(manifest.verify(root), 1)
            artifact.unlink()
            with self.assertRaisesRegex(ValueError, "spectrum.csv"):
                manifest.verify(root)


if __name__ == "__main__":
    unittest.main()
