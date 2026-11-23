"""prep-taxonomy writes OUTPUT/00_prep/gtdbtk_taxonomy.MAG.tsv, also when -o points at 00_prep."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

CLI = Path(__file__).resolve().parent.parent / 'bin' / 'MAGverify'


class PrepTaxonomyTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.home = Path(self.temp.name)
        (self.home / 'gtdbtk').mkdir()
        (self.home / 'gtdbtk/gtdbtk.bac120.summary.tsv').write_text(
            'user_genome\tclassification\nBinA\td__Bacteria;p__X;c__Y;o__Z;f__F;g__G;s__S\n')
        (self.home / 'mags').mkdir()
        (self.home / 'mags/BinA.fa').write_text('>c1\nACGT\n')

    def run_cli(self, output):
        return subprocess.run([sys.executable, str(CLI), 'prep-taxonomy', '-d', str(self.home / 'gtdbtk'),
                               '-g', str(self.home / 'mags'), '-o', str(output)], text=True, capture_output=True)

    def test_output_folder_and_00_prep_give_the_same_table(self):
        for output in (self.home / 'out', self.home / 'out2' / '00_prep'):
            with self.subTest(output=output):
                p = self.run_cli(output)
                self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
                root = output.parent if output.name == '00_prep' else output
                self.assertTrue((root / '00_prep/gtdbtk_taxonomy.MAG.tsv').is_file())
                self.assertFalse((root / '00_prep/00_prep').exists())

    def test_file_path_is_refused(self):
        p = self.run_cli(self.home / 'out' / 'x.tsv')
        self.assertNotEqual(p.returncode, 0)
        self.assertIn('output folder', p.stderr)
        self.assertFalse((self.home / 'out' / 'x.tsv').exists())

    def test_old_assembler_links_are_removed(self):
        links = self.home / 'out/00_prep/assembler_links'
        links.mkdir(parents=True)
        (links / 'a.fa').symlink_to(self.home / 'mags/BinA.fa')
        (links / 'a.fa.fai').write_text('c1\t4\n')
        p = self.run_cli(self.home / 'out')
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertFalse(links.exists())


if __name__ == '__main__':
    unittest.main()
