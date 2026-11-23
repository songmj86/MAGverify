"""phantom_mag.py --contig-output writes one row per assembly contig."""
import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parent.parent / 'scripts' / 'phantom_mag.py'
HEADER = ('region_id\tcontig\tstart_0based\tend_0based\tlength_bp\traw_read_blast_hsps\tbest_hit_pident\t'
          'best_hit_aligned_bp\traw_read_blast_query_covered_bp\traw_read_blast_query_coverage\t'
          'phantom_sequence_call\n')
CALL = '0\tNA\t0\t0\t0.000000\tNo_raw_read_BLAST_hit_phantom_candidate\n'


class PhantomContigTests(unittest.TestCase):
    def test_contig_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / 'phantom.tsv').write_text(HEADER + f'r1\tc1\t10\t60\t50\t{CALL}' + f'r2\tc1\t100\t120\t20\t{CALL}'
                                             + f'r3\tu1\t0\t30\t30\t{CALL}')
            (tmp / 'map.tsv').write_text('MAG_ID\tContig\nBinA\tc1\nBinA\tc2\n')
            (tmp / 'all.fai').write_text('c1\t1000\t0\t60\t61\nc2\t500\t0\t60\t61\nu1\t300\t0\t60\t61\n')
            out = tmp / 'phantom.contig.tsv'
            subprocess.run([sys.executable, str(SCRIPT), '--phantom', str(tmp / 'phantom.tsv'),
                            '--mag-contig-map', str(tmp / 'map.tsv'), '--fai', str(tmp / 'all.fai'),
                            '--site-output', str(tmp / 'sites.tsv'), '--summary-output', str(tmp / 'mag.tsv'),
                            '--contig-output', str(out)], check=True, capture_output=True)
            with out.open() as handle:
                rows = {row['Contig']: row for row in csv.DictReader(handle, delimiter='\t')}
            self.assertEqual(list(rows), ['c1', 'c2', 'u1'])
            self.assertEqual((rows['c1']['MAG_ID'], rows['c1']['n_phantom_regions'], rows['c1']['phantom_bp'],
                              rows['c1']['phantom_regions_0based'], rows['c1']['LR_PHANTOM']),
                             ('BinA', '2', '70', '10-60,100-120', 'Y'))
            self.assertEqual((rows['c2']['LR_PHANTOM'], rows['c2']['phantom_regions_0based']), ('N', 'NA'))
            self.assertEqual((rows['u1']['MAG_ID'], rows['u1']['phantom_fraction']), ('NA', '0.100000'))


if __name__ == '__main__':
    unittest.main()
