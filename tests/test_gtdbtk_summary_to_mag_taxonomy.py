"""GTDB-Tk summary conversion keeps exactly the MAGs given in -m."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / 'scripts' / 'prep_gtdbtk_summary_to_mag_taxonomy.py'
HEADER = 'user_genome\tclassification\tnote\n'


class GtdbtkConversionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        root = Path(self.temp.name)
        self.gtdbtk = root / 'gtdbtk'
        self.gtdbtk.mkdir()
        (self.gtdbtk / 'gtdbtk.bac120.summary.tsv').write_text(
            HEADER
            + 'Bin1L\td__Bacteria;p__Nitrospirota;c__Nitrospiria;o__Nitrospirales;f__X;g__;s__\tx\n'
            + 'Bin2L\td__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__;f__;g__;s__\tx\n'
            + 'Bin3L\tUnclassified Bacteria\tx\n')
        (self.gtdbtk / 'gtdbtk.ar53.summary.tsv').write_text(
            HEADER + 'Bin9L\td__Archaea;p__Halobacteriota;c__Methanomicrobia;o__Methanomicrobiales;f__;g__;s__\tx\n')
        self.mags = root / 'MAGs'
        self.mags.mkdir()
        self.out = root / 'MAG_taxonomy.tsv'

    def run_tool(self, *mags):
        for mag in mags:
            (self.mags / f'{mag}.fa').write_text('>c\nACGT\n')
        return subprocess.run([sys.executable, '-I', str(SCRIPT), '--gtdbtk-dir', str(self.gtdbtk),
                               '-m', str(self.mags), '--output', str(self.out)],
                              text=True, capture_output=True)

    def test_subset_of_gtdbtk_genomes(self):
        result = self.run_tool('Bin1L', 'Bin9L')
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = [line.split('\t') for line in self.out.read_text().splitlines()]
        self.assertEqual(rows[0], ['ID', 'GTDB_Domain', 'GTDB_Phylum', 'GTDB_Class', 'GTDB_Order'])
        self.assertEqual([r[0] for r in rows[1:]], ['Bin1L', 'Bin9L'])
        self.assertEqual(rows[2][1], 'd__Archaea')

    def test_unassigned_rank_keeps_prefix(self):
        result = self.run_tool('Bin2L')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.out.read_text().splitlines()[1].split('\t')[4], 'o__')

    def test_missing_mag_is_an_error(self):
        result = self.run_tool('Bin1L', 'Bin5L')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Bin5L', result.stderr)
        self.assertFalse(self.out.exists())

    def test_unclassified_mag_in_input_is_an_error(self):
        result = self.run_tool('Bin3L')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Bin3L', result.stderr)


if __name__ == '__main__':
    unittest.main()
