"""Final tables stay in the module directories; no OUTPUT/results/ copies are made."""
import contextlib
import io
import tempfile
import unittest
from pathlib import Path

from scripts.core_magverify_cli import FINAL_TABLES, remove_generated_results_dir, report_final_tables


class FinalTableTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_final_tables_follow_module_level_names(self):
        for stage, tables in FINAL_TABLES.items():
            for relative in tables:
                name = Path(relative).name
                if stage != 'map':
                    self.assertRegex(name, r'^[a-z_]+\.(MAG|contig|site|region|dataset|HSP)\.tsv$', name)

    def test_report_lists_existing_tables_only(self):
        (self.root / '06_repeats').mkdir()
        (self.root / '06_repeats/repeats.MAG.tsv').write_text('MAG_ID\nBin1L\n')
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            report_final_tables('repeats', self.root)
        self.assertIn('06_repeats/repeats.MAG.tsv', out.getvalue())
        self.assertNotIn('repeats.contig.tsv', out.getvalue())
        self.assertFalse((self.root / 'results').exists())

    def test_generated_results_dir_is_removed(self):
        results = self.root / 'results'
        results.mkdir()
        (results / 'repeats.MAG.tsv').write_text('x')
        (results / 'circularity_mag.tsv').symlink_to('../nowhere.tsv')
        remove_generated_results_dir(self.root)
        self.assertFalse(results.exists())

    def test_results_dir_with_user_files_is_kept(self):
        results = self.root / 'results'
        results.mkdir()
        (results / 'repeats.MAG.tsv').write_text('x')
        (results / 'my_notes.txt').write_text('keep')
        remove_generated_results_dir(self.root)
        self.assertTrue((results / 'my_notes.txt').exists())


if __name__ == '__main__':
    unittest.main()
