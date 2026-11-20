import json
import tempfile
import unittest
from pathlib import Path

from scripts.core_magverify_cli import PUBLIC_RESULTS, publish_results


class PublicResultTablesTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        prep = self.root / '00_prep'
        prep.mkdir()
        (prep / 'MAG_contig_map.tsv').write_text('MAG_ID\tContig\nBin1L\tctg1\n')
        asm = self.root / 'asm'
        asm.mkdir()
        (asm / 'metaMDBG.fa').write_text('>ctg1\nACGT\n>ctg2\nACGT\n')
        (prep / 'all_assemblies.fasta.sources.json').write_text(
            json.dumps([{'path': str(asm / 'metaMDBG.fa'), 'size': 1, 'mtime_ns': 1}]))

    def write(self, relative, text):
        path = self.root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)

    def test_mag_and_contig_tables_get_standard_leading_columns(self):
        self.write('06_repeats/MAG_repeat_breadth.final.tsv', 'repeat_fraction\tMAG_ID\n0.1\tBin1L\n')
        self.write('06_repeats/contig_repeat_breadth.final.tsv',
                   'contig\tsource_assembler\tfraction\nctg1\tmetaMDBG\t0.2\nctg2\tmetaMDBG\t0.3\n')
        publish_results('repeats', self.root)
        mag = (self.root / 'results/repeats.MAG.tsv').read_text().splitlines()
        self.assertEqual(mag, ['MAG_ID\tAssembler\trepeat_fraction', 'Bin1L\tmetaMDBG\t0.1'])
        contig = (self.root / 'results/repeats.contig.tsv').read_text().splitlines()
        self.assertEqual(contig, ['Contig\tMAG_ID\tAssembler\tfraction',
                                  'ctg1\tBin1L\tmetaMDBG\t0.2', 'ctg2\tNA\tmetaMDBG\t0.3'])

    def test_circularity_id_column_becomes_mag_id_and_old_links_go(self):
        results = self.root / 'results'
        results.mkdir()
        (results / 'circularity_mag.tsv').symlink_to('../nowhere.tsv')
        self.write('02_circularity/MAG_circularity.final.tsv', 'ID\tCircular\nBin1L\tY\n')
        publish_results('circularity', self.root)
        self.assertEqual((results / 'circularity.MAG.tsv').read_text().splitlines()[0],
                         'MAG_ID\tAssembler\tCircular')
        self.assertFalse((results / 'circularity_mag.tsv').is_symlink())

    def test_names_follow_module_level_pattern(self):
        for stage, tables in PUBLIC_RESULTS.items():
            for name, (_, level) in tables.items():
                if level:
                    self.assertTrue(name.endswith(f'.{level}.tsv'), name)

    def test_no_fake_result_when_stage_has_no_outputs(self):
        publish_results('map', self.root)
        self.assertFalse((self.root / 'results').exists())


if __name__ == '__main__':
    unittest.main()
