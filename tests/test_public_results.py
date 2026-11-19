import tempfile
import unittest
from pathlib import Path

from scripts.core_magverify_cli import PUBLIC_RESULTS, publish_results


class PublicResultNamesTest(unittest.TestCase):
    def test_public_results_are_short_links(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / '03_clipping' / 'MAG_clipping.final.tsv'
            source.parent.mkdir(parents=True)
            source.write_text('MAG_ID\tClipping\nBin1L\tN\n', encoding='utf-8')
            publish_results('clipping', root)
            result = root / 'results' / 'clipping_mag.tsv'
            self.assertTrue(result.is_symlink())
            self.assertEqual(result.read_text(), source.read_text())
            self.assertEqual(result.resolve(), source.resolve())
            publish_results('clipping', root)  # Idempotent
            self.assertFalse((root / 'results' / 'clipping_sites.tsv').exists())

    def test_summary_alias_is_mag_level(self):
        self.assertEqual(PUBLIC_RESULTS['summary']['summary_mag.tsv'],
                         '10_summary/integrated_QC_summary.by_MAG.tsv')

    def test_no_fake_result_when_stage_has_no_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            publish_results('map', root)
            self.assertFalse((root / 'results').exists())


if __name__ == '__main__':
    unittest.main()
