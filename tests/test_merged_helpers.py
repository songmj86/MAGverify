"""Regression checks for consolidated internal chimera and summary helpers."""
import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SELECTOR = ROOT / 'scripts/chimera_gene_taxonomy.py'
FINALIZER = ROOT / 'scripts/summary_tables.py'


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
                self.run_script(SELECTOR, 'select', 'fasta', mode,
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
                self.run_script(SELECTOR, 'select', 'genes', mode, self.assembly, self.map,
                                self.genes, self.coords, faa, coord, stats)
                self.assertEqual(faa.read_text(), faa_expected)
                self.assertEqual(coord.read_text(), coord_expected)

    def test_missing_contig_raises(self):
        self.map.write_text('MAG_ID\tcontig\nBin1L\tmissing\n')
        p = subprocess.run([sys.executable, str(SELECTOR), 'select', 'fasta', 'include',
                            str(self.assembly), str(self.map),
                            str(self.root / 'out.fa'), str(self.root / 'out.json')],
                           capture_output=True, text=True)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn('mapped contigs are absent', p.stderr)

    def test_summary_arrange_and_validate(self):
        mag = self.root / 'summary.MAG.tsv'
        contig = self.root / 'summary.contig.tsv'
        fields = ['MAG_ID', 'n_contigs', 'Clipping_Y_N', 'Any_detected_QC_signal_Y_N',
                  'Premature_circularization_candidate', 'Premature_circularization',
                  'MAG_integrated_QC_call']
        with mag.open('w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(fields)
            writer.writerow(['Bin1L', '2', 'Y', 'Y', 'N', 'N', 'QC_signal_detected'])
        contig.write_text('MAG_ID\tcontig\tlength\nBin1L\tc1\t10\nBin1L\tc2\t20\nNA\tc3\t30\n')
        self.run_script(FINALIZER, 'arrange', mag)
        self.run_script(FINALIZER, 'validate', contig, mag)
        with mag.open() as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            row = next(reader)
        self.assertEqual(reader.fieldnames[:3], ['MAG_ID', 'n_contigs', 'MAG_integrated_QC_call'])
        self.assertNotIn('Premature_circularization_candidate', reader.fieldnames)
        self.assertEqual(row['MAG_integrated_QC_call'], 'QC_signal_detected')

if __name__ == '__main__':
    unittest.main()


class MergedDataclassTests(unittest.TestCase):
    """The merged scripts must keep the @dataclass decorators of the original scripts."""

    def test_dataclasses_survive_merge(self):
        import dataclasses
        import importlib.util
        expected = {
            'chimera_call': ('Gene', 'SideBlock', 'GeneIndex', 'PositionEvaluation',
                             'ClipRecord', 'RawShift', 'ShiftCluster'),
            'false_closure': ('Alignment', 'MagInfo', 'SourceMapping', 'Breakpoint',
                              'LargeInsertionSignal'),
            'false_duplication': ('ClippingSite', 'MappingSource', 'DuplicationPair',
                                  'CigarEvent', 'BoundaryEvidence'),
        }
        scripts = Path(__file__).resolve().parent.parent / 'scripts'
        sys.path.insert(0, str(scripts))
        for module_name, classes in expected.items():
            spec = importlib.util.spec_from_file_location('merged_' + module_name, scripts / f'{module_name}.py')
            module = importlib.util.module_from_spec(spec)
            sys.modules[spec.name] = module
            spec.loader.exec_module(module)
            for name in classes:
                with self.subTest(script=module_name, cls=name):
                    self.assertTrue(dataclasses.is_dataclass(getattr(module, name)))
        clip = sys.modules['merged_chimera_call'].ClipRecord(10, 5, 2, 0.4, True)
        self.assertEqual(clip.position, 10)
