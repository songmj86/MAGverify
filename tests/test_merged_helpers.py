"""Regression checks for consolidated internal chimera and summary helpers."""
import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SELECTOR = ROOT / 'scripts/chimera_select_sequences.py'
FINALIZER = ROOT / 'scripts/summary_finalize.py'


class MergedHelperTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)
        self.assembly = self.root / 'assembly.fa'
        self.assembly.write_text('>c1 mag\nATGC\n>c2 unbinned\nACGTAC\n')
        self.map = self.root / 'map.tsv'
        self.map.write_text('MAG_ID\tcontig\nBin1L\tc1\n')
        self.genes = self.root / 'genes.faa'
        self.genes.write_text('>gene1\nMK\n>gene2\nMV\n')
        self.coords = self.root / 'genes.tsv'
        self.coords.write_text('gene1\tc1\t1\t4\ngene2\tc2\t1\t6\n')

    def run_script(self, script, *args):
        p = subprocess.run([sys.executable, str(script), *map(str, args)],
                           capture_output=True, text=True)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        return p

    def test_fasta_include_and_exclude(self):
        for mode, output in (('include', '>c1 mag\nATGC\n'),
                             ('exclude', '>c2 unbinned\nACGTAC\n')):
            with self.subTest(mode=mode):
                fasta = self.root / f'{mode}.fa'
                stats = self.root / f'{mode}.json'
                self.run_script(SELECTOR, 'fasta', mode,
                                self.assembly, self.map, fasta, stats)
                self.assertEqual(fasta.read_text(), output)
                self.assertIn('"mapped_contigs": 1', stats.read_text())

    def test_gene_include_and_exclude(self):
        for mode, faa_expected, coord_expected in (
            ('include', '>gene1\nMK\n', 'gene1\tc1\t1\t4\n'),
            ('exclude', '>gene2\nMV\n', 'gene2\tc2\t1\t6\n'),
        ):
            with self.subTest(mode=mode):
                faa = self.root / f'{mode}.faa'
                coord = self.root / f'{mode}.tsv'
                stats = self.root / f'genes_{mode}.json'
                self.run_script(SELECTOR, 'genes', mode, self.assembly, self.map,
                                self.genes, self.coords, faa, coord, stats)
                self.assertEqual(faa.read_text(), faa_expected)
                self.assertEqual(coord.read_text(), coord_expected)

    def test_missing_contig_raises(self):
        self.map.write_text('MAG_ID\tcontig\nBin1L\tmissing\n')
        p = subprocess.run([sys.executable, str(SELECTOR), 'fasta', 'include',
                            str(self.assembly), str(self.map),
                            str(self.root / 'out.fa'), str(self.root / 'out.json')],
                           capture_output=True, text=True)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn('mapped contigs are absent', p.stderr)

    def test_summary_complete_genome(self):
        source = self.root / 'integrated.tsv'
        output = self.root / 'complete.tsv'
        fields = ['MAG_ID', 'Circular', 'Taxonomic_chimera_Y_N',
                  'Phantom_sequence_Y_N', 'Excessive_repeat_Y_N',
                  'False_duplication_Y_N', 'False_closure_Y_N',
                  'Premature_circularization_Y_N']
        with source.open('w') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(fields)
            writer.writerow(['Bin1L', 'Disabled', 'N', 'N', 'N', 'N', 'Y', 'N'])
        self.run_script(FINALIZER, 'complete-genome', source, output)
        with output.open() as handle:
            row = next(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(row['Circular'], 'NA')
        self.assertEqual(row['False_closure_or_premature_circularization'], 'Y')


if __name__ == '__main__':
    unittest.main()
