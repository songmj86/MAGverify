"""A partial QC run must keep unavailable circularity calls unassessed."""
import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'scripts/summary_finalize.py'

class PartialSummaryTests(unittest.TestCase):
    def test_unassessed_circularity_is_na(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'summary.tsv'
            output = Path(directory) / 'validity.tsv'
            fields = ['MAG_ID', 'Circular', 'Taxonomic_chimera_Y_N',
                      'Phantom_sequence_Y_N', 'Excessive_repeat_Y_N',
                      'False_duplication_Y_N', 'False_closure_Y_N',
                      'Premature_circularization_Y_N']
            with source.open('w', newline='') as handle:
                w=csv.writer(handle,delimiter='\t')
                w.writerow(fields)
                w.writerow(['Bin1L','Disabled','N','N','N','N','N','N'])
                w.writerow(['Bin2L','Y','N','N','N','N','N','N'])
            p=subprocess.run([sys.executable,str(SCRIPT),'complete-genome',str(source),str(output)],
                             text=True,capture_output=True)
            self.assertEqual(p.returncode,0,p.stderr)
            with output.open() as handle:
                rows=list(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual([r['Circular'] for r in rows],['NA','Y'])

if __name__=='__main__':
    unittest.main()
