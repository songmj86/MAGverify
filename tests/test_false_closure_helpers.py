import importlib
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from false_closure_scg_utils import load_ribosomal_families  # noqa: E402


SCG_TSV = """contig\thmm_name\tgene_callers_id\tsource
c1\tBacteria_71|Ribosomal_L2\t1\tBacteria_71
c1\tBacteria_71|Ribosomal_L3\t2\tBacteria_71
c1\tBacteria_71|Ribosomal_L3\t3\tBacteria_71
c1\tBacteria_71|ribosomal_L24\t4\tBacteria_71
c1\tBacteria_71|ADK\t5\tBacteria_71
c1\tArchaea_76|Ribosomal_L2\t1\tArchaea_76
c2\tArchaea_76|Ribosomal_S8\t7\tArchaea_76
c3\tBacteria_71|gyrB\t8\tBacteria_71
"""


class FalseClosureHelperTests(unittest.TestCase):
    def test_ribosomal_families_per_contig(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "scg.tsv"
            path.write_text(SCG_TSV)
            families = load_ribosomal_families(str(path))
        self.assertEqual(families["c1"], ("Ribosomal_L2", "Ribosomal_L3", "ribosomal_L24"))
        self.assertEqual(families["c2"], ("Ribosomal_S8",))
        self.assertNotIn("c3", families)

    def test_false_closure_scripts_import(self):
        for name in (
            "false_closure_detect_false_closure",
            "false_closure_resolve_false_closure_contig_circularity",
            "false_closure_summarize_false_closure_by_contig",
        ):
            importlib.import_module(name)


if __name__ == "__main__":
    unittest.main()
