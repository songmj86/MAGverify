"""Result tables are written from the internal tables in .cache/."""
from pathlib import Path
import csv
import re
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'scripts'))
import core_public_tables as pt  # noqa: E402


def write(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text('\t'.join(header) + '\n' + ''.join('\t'.join(r) + '\n' for r in rows))


def read(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def summary_header(start):
    src = (ROOT / 'scripts/summary_tables.py').read_text()
    i = src.index(start)
    j = src.index("\\n'", i) + 3
    return ''.join(re.findall(r"'([^']*)'", src[i:j])).replace('\\n', '').split('\\t')


class PublicTableTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_names_and_values(self):
        self.assertEqual(pt.plain_name('HIFI_PRIMARY_SUPPLEMENTARY_END_START_SUPPORTED_Y_N'),
                         'HiFi_END_START_supported')
        self.assertEqual(pt.plain_name('LR_COPY1_LEFT_SUPPORTED'), 'Copy1_left_supported')
        self.assertEqual(pt.plain_name('SOURCE_CONTIG_ASSEMBLER_CIRCULARITY'), 'Assembler_reported_circularity')
        self.assertEqual(pt.plain_name('Paper_small_genome_proxy'), 'Small_genome_proxy')
        self.assertEqual(pt.plain_name('n_contigs'), 'n_contigs')
        self.assertNotIn('ource', ' '.join(pt.plain_name(n) for n in pt.RENAME))

    def test_mapping_table_has_no_paths(self):
        write(self.root / '01_mapping/.cache/mapping_manifest.tsv',
              ['Assembler', 'Assembly_FASTA', 'HiFi_all_BAM', 'HiFi_primary_BAM', 'HiFi_primary_supplementary_BAM',
               'Output_directory', 'Reads', 'Primary', 'Mapped_percent'],
              [['metaflye', '/a.fa', '/b.bam', '/c.bam', '/d.bam', '/o', '/r.fq', '10', '100.00']])
        pt.publish(self.root)
        rows = read(self.root / '01_mapping/mapping.assembler.tsv')
        self.assertEqual(list(rows[0]), ['Assembler', 'Mapped_reads'])

    def test_circularity_contig(self):
        write(self.root / '02_circularity/.cache/circularity.contig.tsv',
              ['Contig', 'MAG_ID', 'Assembler', 'Length_bp', 'MAG_n_contigs', 'Circularity_eligibility',
               'Circularity_evaluable', 'Circular'],
              [['c1', 'B1', 'flye', '10', '1', 'Y_single_contig', 'Y', 'Y'],
               ['c2', 'B2', 'flye', '10', '2', 'N_not_single_contig', 'N', 'NA']])
        pt.publish(self.root)
        rows = read(self.root / '02_circularity/circularity.contig.tsv')
        self.assertNotIn('Circularity_eligibility', rows[0])
        self.assertEqual(rows[1]['n_contigs'], '2')

    def test_summary_tables_go_to_the_top(self):
        contig_header = summary_header("'MAG_ID\\tcontig\\tlength")
        mag_header = summary_header("'MAG_ID\\tn_contigs\\t'")
        contig = dict.fromkeys(contig_header, 'NA')
        contig.update(MAG_ID='B1', contig='c1', length='100', validated_clipping_breakpoints='2',
                      Taxonomic_chimera='Not_run', phantom_sequence_call='Phantom_sequence_candidate',
                      assembly_allele_support_call='No_assembly_allele_disagreement',
                      excessive_repeat_call='No_excessive_repeat_signal',
                      false_duplication_call='No_false_duplication_signal', n_flagged_modules='2',
                      integrated_interpretation='QC_signal_detected')
        mag = dict.fromkeys(mag_header, 'NA')
        mag.update(MAG_ID='B1', n_contigs='1', Clipping_Y_N='Y', Taxonomic_chimera_Y_N='N',
                   Phantom_sequence_Y_N='Y', Circular='Not_run', False_closure_evidence='Not_run',
                   MAG_integrated_QC_call='QC_signal_detected')
        write(self.root / '10_summary/summary.contig.tsv', contig_header, [[contig[k] for k in contig_header]])
        write(self.root / '10_summary/summary.MAG.tsv', mag_header, [[mag[k] for k in mag_header]])
        pt.migrate(self.root)
        self.assertFalse((self.root / '10_summary').exists())
        pt.publish(self.root)
        c = read(self.root / 'summary.contig.tsv')[0]
        self.assertEqual((c['High_confidence_clipping'], c['Chimera'], c['Phantom_sequence']),
                         ('NA', 'NA', 'Y'))
        detail_c = read(self.root / 'detail.contig.tsv')[0]
        self.assertEqual((detail_c['Clipping_candidate'], detail_c['Checks_not_run']),
                         ('Y', 'taxonomic_chimera'))
        for removed in ('QC_call', 'MAG_integrated_QC_call', 'integrated_interpretation', 'excessive_repeat_call',
                        'false_duplication_clipping_summary', 'Circularity_eligibility'):
            self.assertNotIn(removed, c)
        m = read(self.root / 'summary.MAG.tsv')[0]
        self.assertEqual(list(m), pt.SUMMARY_MAG_FIELDS)
        self.assertEqual(list(c), pt.SUMMARY_CONTIG_FIELDS)
        for removed in ('n_contigs', 'n_flagged_checks', 'n_checks_not_run', 'Clipping', 'Phantom'):
            self.assertNotIn(removed, m)
        detail_m = read(self.root / 'detail.MAG.tsv')[0]
        for removed in ('n_flagged_checks', 'n_checks_not_run'):
            self.assertNotIn(removed, detail_m)
        self.assertEqual(detail_m['Clipping_candidate'], 'Y')
        for removed in ('QC_call', 'MAG_integrated_QC_call', 'Circularity_eligibility', 'excessive_repeat_call',
                        'false_duplication_clipping_summary'):
            self.assertNotIn(removed, m)
        self.assertEqual((m['Chimera'], m['Read_supported_circularity'], m['Premature_circularization']), ('NA', 'NA', 'NA'))
        self.assertEqual(list(detail_m)[:len(pt.DETAIL_MAG_LEAD)], pt.DETAIL_MAG_LEAD)
        self.assertIn('circularity', read(self.root / 'detail.MAG.tsv')[0]['Checks_not_run'])

    def test_summary_has_yes_no_and_full_versions(self):
        self.test_summary_tables_go_to_the_top()
        write(self.root / '03_clipping/.cache/clipping.MAG.tsv', ['MAG_ID', 'n_contigs', 'Clipping', 'n_clipping_sites'],
              [['B1', '1', 'Y', '2']])
        write(self.root / '03_clipping/.cache/clipping.contig.tsv', ['contig', 'MAG_ID', 'length_bp', 'Clipping'],
              [['c1', 'B1', '100', 'Y'], ['u1', 'NA', '50', 'N']])
        write(self.root / '05_phantom/.cache/phantom.contig.tsv', ['Contig', 'MAG_ID', 'length_bp', 'LR_PHANTOM'],
              [['c1', 'B1', '100', 'Y']])
        pt.publish(self.root)
        m = read(self.root / 'summary.MAG.tsv')[0]
        for detail in ('n_clipping_sites', 'Repeat_fraction', 'False_closure_partner_contigs'):
            self.assertNotIn(detail, m)
        self.assertTrue(all(m[k] in ('Y', 'N', 'NA')
                            for k in ('High_confidence_clipping', 'Phantom_sequence',
                                      'Read_supported_circularity')))
        full = read(self.root / 'detail.MAG.tsv')[0]
        self.assertEqual(list(full)[:len(pt.DETAIL_MAG_LEAD) + len(pt.SUMMARY_MAG_DETAILS)],
                         pt.DETAIL_MAG_LEAD + pt.SUMMARY_MAG_DETAILS)
        self.assertEqual((full['n_clipping_sites'], full['clipping.Clipping'], full['clipping.n_clipping_sites']),
                         ('NA', 'Y', '2'))
        self.assertNotIn('clipping.MAG_ID', full)
        contigs = read(self.root / 'detail.contig.tsv')
        self.assertEqual([row['Contig'] for row in contigs], ['c1'])
        self.assertEqual((contigs[0]['clipping.Clipping'], contigs[0]['phantom.Phantom'], contigs[0]['n_clipping_sites']),
                         ('Y', 'Y', '2'))
        self.assertNotIn('Chimera_ranks', read(self.root / 'summary.contig.tsv')[0])

    def test_module_order_and_summary_counts(self):
        self.assertEqual(pt.SUMMARY_MAG_FIELDS, [
            'MAG_ID', 'Assembler_reported_circularity',
            'Read_supported_circularity', 'High_confidence_circularity',
            'High_confidence_clipping', 'Premature_circularization',
            'Phantom_sequence', 'Excessive_repeat', 'False_duplication',
            'Chimera', 'Assembly_variant'])
        self.assertEqual(pt.SUMMARY_CONTIG_FIELDS, [
            'Contig', 'MAG_ID', 'Length_bp', 'High_confidence_clipping',
            'Phantom_sequence', 'Excessive_repeat', 'False_duplication',
            'Chimera', 'Assembly_variant'])
        self.assertEqual([name for name, _ in pt.FULL_MODULES], [
            'circularity', 'clipping', 'false_closure', 'phantom', 'repeats',
            'false_duplication', 'chimera', 'variants'])
        self.assertEqual(pt.DETAIL_MODULE_FLAGS,
                         {'chimera': ('Chimera',), 'variants': ('Assembly_variant',)})
        # Changing the column order does not change which flags enter the count.
        self.assertEqual(len(pt.MAG_CHECKS), 7)

    def test_existing_summary_republished_when_header_order_changes(self):
        self.test_summary_tables_go_to_the_top()
        old_path = self.root / 'summary.MAG.tsv'
        with open(old_path, encoding='utf-8', newline='') as h:
            rows = list(csv.DictReader(h, delimiter='\t'))
        old_fields = ['MAG_ID', 'n_contigs', 'n_flagged_checks', 'Taxonomic_chimera',
                      'Clipping', 'False_closure', 'Premature_circularization', 'Circular',
                      'Checks_not_run']
        pt.write(old_path, old_fields, rows)
        self.assertFalse(pt.summary_headers_current(self.root))
        # All cached tables are unchanged. Publishing must still update the headers.
        pt.publish(self.root)
        self.assertTrue(pt.summary_headers_current(self.root))
        self.assertEqual(list(read(old_path)[0]), pt.SUMMARY_MAG_FIELDS)

    def test_three_circularity_evidence_columns_preserve_na_and_independent_support(self):
        self.test_summary_tables_go_to_the_top()
        path = self.root / '02_circularity/.cache/circularity.MAG.tsv'
        write(path, ['ID', 'Assembler_reported_circularity',
                     'Read_supported_circular_topology', 'High_confidence_circularity'],
              [['B1', 'N', 'Y', 'N']])
        pt.publish(self.root)
        row = read(self.root / 'summary.MAG.tsv')[0]
        self.assertEqual((row['Assembler_reported_circularity'],
                          row['Read_supported_circularity'],
                          row['High_confidence_circularity']), ('N', 'Y', 'N'))
        # The old final Circular could be Y when assembler evidence was unavailable.
        # The new both-supported column must remain NA, not copy that old alias.
        write(path, ['ID', 'Assembler_reported_circularity',
                     'Read_supported_circular_topology', 'High_confidence_circularity'],
              [['B1', 'NA', 'Y', 'NA']])
        pt.publish(self.root, force=True)
        row = read(self.root / 'summary.MAG.tsv')[0]
        self.assertEqual((row['Assembler_reported_circularity'],
                          row['Read_supported_circularity'],
                          row['High_confidence_circularity']), ('NA', 'Y', 'NA'))
        write(path, ['ID', 'Assembler_reported_circularity',
                     'Read_supported_circular_topology', 'High_confidence_circularity'],
              [['B1', 'Y', 'Y', 'Y']])
        pt.publish(self.root, force=True)
        row = read(self.root / 'summary.MAG.tsv')[0]
        self.assertEqual((row['Assembler_reported_circularity'],
                          row['Read_supported_circularity'],
                          row['High_confidence_circularity']), ('Y', 'Y', 'Y'))

    def test_basic_summary_does_not_conflate_distinct_diagnostic_calls(self):
        self.test_summary_tables_go_to_the_top()
        public = read(self.root / 'summary.MAG.tsv')[0]
        detailed = read(self.root / 'detail.MAG.tsv')[0]
        self.assertNotIn('False_closure', public)
        self.assertNotIn('Premature_circularization_screen', public)
        self.assertNotIn('False_duplication_with_clipping', public)
        self.assertIn('False_closure', detailed)
        self.assertIn('Premature_circularization_screen', detailed)
        self.assertIn('False_duplication_with_clipping', detailed)
        self.assertIn('Checks_not_run', detailed)

    def test_previous_full_summary_names_migrate_to_detail(self):
        for level in ('contig', 'MAG'):
            (self.root / f'summary_full.{level}.tsv').write_text('old\n')
        pt.migrate(self.root)
        for level in ('contig', 'MAG'):
            self.assertFalse((self.root / f'summary_full.{level}.tsv').exists())
            self.assertEqual((self.root / f'detail.{level}.tsv').read_text(), 'old\n')

    def test_summary_detail_names_migrate_to_detail(self):
        for level in ('contig', 'MAG'):
            (self.root / f'summary_detail.{level}.tsv').write_text('older\n')
        pt.migrate(self.root)
        for level in ('contig', 'MAG'):
            self.assertFalse((self.root / f'summary_detail.{level}.tsv').exists())
            self.assertEqual((self.root / f'detail.{level}.tsv').read_text(), 'older\n')

    def test_earlier_tables_move_to_cache(self):
        write(self.root / '06_repeats/repeats.dataset.tsv', ['scope_id'], [['ALL']])
        write(self.root / '06_repeats/repeats.MAG.tsv', ['MAG_ID', 'INTRA_REPEAT'], [['B1', 'Y']])
        pt.migrate(self.root)
        pt.publish(self.root)
        self.assertTrue((self.root / '06_repeats/.cache/repeats.dataset.tsv').is_file())
        self.assertFalse((self.root / '06_repeats/repeats.dataset.tsv').exists())
        self.assertEqual(read(self.root / '06_repeats/repeats.MAG.tsv')[0], {'MAG_ID': 'B1', 'Intra_repeat': 'Y'})


if __name__ == '__main__':
    unittest.main()
