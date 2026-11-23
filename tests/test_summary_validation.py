"""Tests for summary-table validation."""
import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
VALIDATOR = ROOT / 'scripts/summary_tables.py'


def write_tsv(path, header, rows):
    with path.open('w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(header)
        writer.writerows(rows)


class SummaryValidationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        d = Path(self.tmp.name)
        self.files = [d / name for name in ('contig.tsv', 'mag.tsv')]
        write_tsv(self.files[0], ['MAG_ID','contig','length'],
                  [['Bin1L','c1','1200'],['Bin1L','c2','4000'],['NA','c3','200']])
        self.mag_header = ['MAG_ID','n_contigs','MAG_integrated_QC_call','Any_detected_QC_signal_Y_N']
        write_tsv(self.files[1], self.mag_header, [['Bin1L','2','No_major_signal_in_completed_modules','N']])

    def check(self):
        return subprocess.run([sys.executable, str(VALIDATOR), 'validate', *map(str,self.files)],
                              text=True, capture_output=True)

    def test_consistent_tables(self):
        p = self.check()
        self.assertEqual(p.returncode, 0, p.stdout+p.stderr)
        self.assertIn('1 MAGs, 3 contigs',p.stdout)

    def test_disagreement_is_error(self):
        write_tsv(self.files[1], self.mag_header, [['Bin1L','3','QC_signal_detected','N']])
        p = self.check()
        self.assertNotEqual(p.returncode,0)
        self.assertIn('Contig count mismatch',p.stderr)

    def test_mag_missing_from_contigs_is_error(self):
        write_tsv(self.files[1], self.mag_header, [['Bin2L','2','QC_signal_detected','N']])
        p = self.check()
        self.assertNotEqual(p.returncode,0)
        self.assertIn('MAG membership differs',p.stderr)


if __name__=='__main__':
    unittest.main()
