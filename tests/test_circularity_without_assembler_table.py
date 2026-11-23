"""Circularity results without an assembler-reported circularity table."""
import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / 'scripts' / 'circularity_finalize.py'


class CircularityWithoutAssemblerTableTests(unittest.TestCase):
    def test_assembler_calls_are_na(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            canonical = root / 'canonical.tsv'
            canonical.write_text(
                'MAG_ID\tNumber_of_contigs\tCIRCULARITY_ANALYZABLE\tHIFI_END_START_SUPPORTED\n'
                'Bin1L\t1\tY\tY\n'
                'Bin2L\t3\tN\tNA\n')
            final = root / 'final.tsv'
            result = subprocess.run([sys.executable, '-I', str(SCRIPT), 'refresh', '--canonical', str(canonical),
                                     '--assembler-circularity', '', '--final', str(final)],
                                    text=True, capture_output=True, cwd=root)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with final.open() as handle:
                rows = {row['ID']: row for row in csv.DictReader(handle, delimiter='\t')}
            self.assertEqual(rows['Bin1L']['Assembler_reported_circularity'], 'NA')
            self.assertEqual(rows['Bin1L']['Read_supported_circular_topology'], 'Y')
            self.assertEqual(rows['Bin1L']['High_confidence_circular'], 'NA')
            self.assertEqual(rows['Bin1L']['Additional_circular_candidate'], 'NA')
            self.assertEqual(rows['Bin1L']['Circular'], 'Y')
            self.assertEqual(rows['Bin2L']['Read_supported_circular_topology'], 'NA')


if __name__ == '__main__':
    unittest.main()
