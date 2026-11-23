"""Circularity requires an assembler-reported circularity table."""
import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / 'scripts' / 'circularity_finalize.py'


class CircularityWithoutAssemblerTableTests(unittest.TestCase):
    def test_assembler_calls_are_required(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            canonical = root / 'canonical.tsv'
            canonical.write_text(
                'MAG_ID\tNumber_of_contigs\tCIRCULARITY_ANALYZABLE\tHIFI_END_START_SUPPORTED\n'
                'Bin1L\t1\tY\tY\n'
                'Bin2L\t3\tN\tNA\n')
            final = root / 'final.tsv'
            result = subprocess.run([sys.executable, '-I', str(SCRIPT), 'refresh', '--canonical', str(canonical),
                                     '--final', str(final)],
                                    text=True, capture_output=True, cwd=root)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('the following arguments are required: --assembler-circularity', result.stderr)
            self.assertFalse(final.exists())


if __name__ == '__main__':
    unittest.main()
