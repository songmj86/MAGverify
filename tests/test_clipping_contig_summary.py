"""Contig- and MAG-level clipping summaries from one all-region site table."""
import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / 'scripts' / 'clipping_summarize.py'

SITE_FIELDS = [
    'contig', 'length', 'pos_1based', 'coverage', 'clipping_count', 'clipping_ratio',
    'default_candidate', 'stringent_candidate', 'Clipping_ge1', 'Clipping_analyzable',
    'Clipping', 'Stringent_clipping', 'region_type', 'n_unique_all_clipped_qnames',
    'n_unique_exact_spanning_qnames',
    'exact_spanning_min_flank_bp_min', 'exact_spanning_min_flank_bp_median',
    'exact_spanning_min_flank_bp_max',
    'n_unique_spanning_with_I_qnames', 'spanning_I_length_bp_min',
    'spanning_I_length_bp_median', 'spanning_I_length_bp_max',
    'n_unique_spanning_with_D_qnames', 'spanning_D_length_bp_min',
    'spanning_D_length_bp_median', 'spanning_D_length_bp_max',
    'n_unique_spanning_with_N_qnames', 'spanning_N_length_bp_min',
    'spanning_N_length_bp_median', 'spanning_N_length_bp_max',
    'n_unique_spanning_with_any_IDN_qnames',
]


def site(contig, length, pos, coverage, clipped, clipping, region='INTERNAL'):
    row = {field: 'NA' for field in SITE_FIELDS}
    row.update(contig=contig, length=length, pos_1based=pos, coverage=coverage,
               clipping_count=clipped, clipping_ratio=f'{clipped / coverage:.6f}',
               default_candidate='Y', stringent_candidate='N',
               Clipping_ge1='Y' if clipping == 'Y' else 'N',
               Clipping_analyzable='Y', Clipping=clipping, Stringent_clipping='N',
               region_type=region, n_unique_all_clipped_qnames=clipped)
    for field in SITE_FIELDS:
        if field.startswith('n_unique_') and field != 'n_unique_all_clipped_qnames':
            row[field] = 0
    return row


def read_tsv(path):
    with open(path, newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


class ClippingContigSummaryTests(unittest.TestCase):
    def test_contig_and_mag_tables(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / 'assembly.fa.fai').write_text(
                'c1\t1000\t0\t60\t61\nc2\t2000\t0\t60\t61\nc3\t500\t0\t60\t61\n')
            (root / 'map.tsv').write_text('MAG_ID\tContig\nBinA\tc1\nBinA\tc2\n')
            with open(root / 'sites.tsv', 'w', newline='') as handle:
                writer = csv.DictWriter(handle, SITE_FIELDS, delimiter='\t', lineterminator='\n')
                writer.writeheader()
                writer.writerow(site('c1', 1000, 500, 10, 9, 'Y'))
                writer.writerow(site('c2', 2000, 100, 10, 9, 'N', 'TERMINAL'))
                writer.writerow(site('c3', 500, 250, 10, 9, 'Y'))
            out = {name: root / f'{name}.tsv' for name in ('mag', 'contig')}
            result = subprocess.run([
                sys.executable, '-I', str(SCRIPT), 'mag',
                '--clipping', str(root / 'sites.tsv'),
                '--mag-contig-map', str(root / 'map.tsv'),
                '--fai', str(root / 'assembly.fa.fai'),
                '--site-output', str(root / 'sites.tsv'),
                '--summary-output', str(out['mag']),
                '--contig-summary-output', str(out['contig']),
            ], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

            contigs = {row['contig']: row for row in read_tsv(out['contig'])}
            self.assertEqual(sorted(contigs), ['c1', 'c2', 'c3'])
            self.assertEqual(list(read_tsv(out['contig'])[0])[:3], ['contig', 'MAG_ID', 'length_bp'])
            self.assertEqual(contigs['c1']['MAG_ID'], 'BinA')
            self.assertEqual(contigs['c3']['MAG_ID'], 'NA')
            self.assertEqual(contigs['c1']['length_bp'], '1000')
            self.assertEqual((contigs['c1']['Clipping'], contigs['c1']['Internal_clipping']), ('Y', 'Y'))
            self.assertEqual((contigs['c2']['Clipping'], contigs['c2']['n_candidate_sites']), ('N', '1'))
            self.assertEqual(contigs['c3']['Clipping'], 'Y')
            self.assertNotIn('clipping_contigs', contigs['c1'])

            sites = read_tsv(root / 'sites.tsv')
            self.assertEqual(len(sites), 3)
            self.assertEqual([row['MAG_ID'] for row in sites], ['BinA', 'BinA', 'NA'])
            self.assertEqual(list(sites[0]).count('MAG_ID'), 1)

            mags = read_tsv(out['mag'])
            self.assertEqual([row['MAG_ID'] for row in mags], ['BinA'])
            self.assertEqual((mags[0]['Clipping'], mags[0]['n_contigs'], mags[0]['clipping_contigs']),
                             ('Y', '2', 'c1'))


if __name__ == '__main__':
    unittest.main()
