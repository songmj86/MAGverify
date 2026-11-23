import importlib
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from false_closure import load_ribosomal_families  # noqa: E402


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
        for name in ("false_closure", "circularity_finalize"):
            importlib.import_module(name)


if __name__ == "__main__":
    unittest.main()


class BreakpointDetailsTests(unittest.TestCase):
    def test_every_breakpoint_in_one_cell(self):
        import false_closure
        rows = [
            {'Region_type': 'internal', 'Breakpoint_1based': '500', 'Distance_to_nearest_end_bp': '400',
             'Left_cross_contig_connected_QNAMEs': '3', 'Left_clipped_QNAMEs': '5',
             'Left_target_contigs': 'ctgA:3;ctgB:1', 'Right_cross_contig_connected_QNAMEs': '0',
             'Right_clipped_QNAMEs': '2', 'Right_target_contigs': 'NA', 'Coverage': '20',
             'Clipping_count': '7', 'Clipping_fraction': '0.35', 'Validated_clean_spanning_QNAMEs': '1',
             'Validated_all_clipped_QNAMEs': '7', 'False_closure_evidence': 'N', 'False_closure_reason': 'R1'},
            {'Region_type': 'TERMINAL', 'Breakpoint_1based': '10', 'Distance_to_nearest_end_bp': '9',
             'Left_cross_contig_connected_QNAMEs': '2', 'Left_clipped_QNAMEs': '2',
             'Left_target_contigs': 'ctgC:2', 'Right_cross_contig_connected_QNAMEs': '2',
             'Right_clipped_QNAMEs': '3', 'Right_target_contigs': 'ctgC:2', 'Coverage': '9',
             'Clipping_count': '3', 'Clipping_fraction': '0.3', 'Validated_clean_spanning_QNAMEs': '0',
             'Validated_all_clipped_QNAMEs': '3', 'False_closure_evidence': 'Y', 'False_closure_reason': 'R2'},
        ]
        self.assertEqual(
            false_closure.breakpoint_details(rows),
            'TERMINAL@10(9)|L=2/2[ctgC:2]|R=2/3[ctgC:2]|cov=9,clip=3,frac=0.3|spanning=0/3|evidence=Y|R2;'
            'INTERNAL@500(400)|L=3/5[ctgA:3,ctgB:1]|R=0/2|cov=20,clip=7,frac=0.35|spanning=1/7|evidence=N|R1')
        self.assertEqual(false_closure.breakpoint_details([]), 'NA')
