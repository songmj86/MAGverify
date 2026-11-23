#!/usr/bin/env python3
"""Write the result tables of each module folder from the tables in .cache/.

Subcommands:
  migrate ROOT   Move result tables of the module folders into .cache/
  publish ROOT   Write the result tables (column names, Y/N/NA values, no paths)
"""
from __future__ import annotations

import csv
import os
import re
import shutil
import sys
from pathlib import Path

TABLES = (
    ('01_mapping/mapping.assembler.tsv', '01_mapping/.cache/mapping_manifest.tsv', 'mapping'),
    ('02_circularity/circularity.MAG.tsv', '02_circularity/.cache/circularity.MAG.tsv', 'circularity_mag'),
    ('02_circularity/circularity.contig.tsv', '02_circularity/.cache/circularity.contig.tsv', 'circularity_contig'),
    ('03_clipping/clipping.site.tsv', '03_clipping/.cache/clipping.site.tsv', 'clipping_site'),
    ('03_clipping/clipping.contig.tsv', '03_clipping/.cache/clipping.contig.tsv', 'contig'),
    ('03_clipping/clipping.MAG.tsv', '03_clipping/.cache/clipping.MAG.tsv', 'mag'),
    ('03_clipping/clipping.assembler.tsv', '03_clipping/.cache/clipping.assembler.tsv', 'clipping_assembler'),
    ('04_false_closure/false_closure.contig.tsv', '04_false_closure/.cache/false_closure.contig.tsv', 'false_closure'),
    ('04_false_closure/false_closure.MAG.tsv', '04_false_closure/.cache/false_closure.MAG.tsv', 'false_closure'),
    ('05_phantom/phantom.region.tsv', '05_phantom/.cache/phantom.region.tsv', 'phantom_region'),
    ('05_phantom/phantom.contig.tsv', '05_phantom/.cache/phantom.contig.tsv', 'contig'),
    ('05_phantom/phantom.MAG.tsv', '05_phantom/.cache/phantom.MAG.tsv', 'mag'),
    ('05_phantom/phantom.dataset.tsv', '05_phantom/.cache/phantom.dataset.tsv', None),
    ('06_repeats/repeats.contig.tsv', '06_repeats/.cache/repeats.contig.tsv', 'contig'),
    ('06_repeats/repeats.MAG.tsv', '06_repeats/.cache/repeats.MAG.tsv', 'mag'),
    ('06_repeats/repeats.region.tsv', '06_repeats/.cache/repeats.region.tsv', 'contig'),
    ('06_repeats/repeats.HSP.tsv', '06_repeats/.cache/repeats.HSP.tsv', 'hsp'),
    ('06_repeats/repeats.dataset.tsv', '06_repeats/.cache/repeats.dataset.tsv', None),
    ('07_chimera/chimera.region.tsv', '07_chimera/.cache/chimera.region.tsv', 'chimera_region'),
    ('07_chimera/chimera.contig.tsv', '07_chimera/.cache/chimera.contig.tsv', 'chimera_contig'),
    ('07_chimera/chimera.MAG.tsv', '07_chimera/.cache/chimera.MAG.tsv', 'chimera_mag'),
    ('08_false_duplication/false_duplication.pair.tsv', '08_false_duplication/.cache/false_duplication.pair.tsv', 'false_dup_pair'),
    ('08_false_duplication/false_duplication.contig.tsv', '08_false_duplication/.cache/false_duplication.contig.tsv', 'false_dup'),
    ('08_false_duplication/false_duplication.MAG.tsv', '08_false_duplication/.cache/false_duplication.MAG.tsv', 'false_dup'),
    ('09_assembly_variants/variants.site.tsv', '09_assembly_variants/.cache/variants.site.tsv', 'contig'),
    ('09_assembly_variants/variants.contig.tsv', '09_assembly_variants/.cache/variants.contig.tsv', 'contig'),
    ('09_assembly_variants/variants.MAG.tsv', '09_assembly_variants/.cache/variants.MAG.tsv', 'mag'),
    ('summary.contig.tsv', '.cache/summary/summary.contig.tsv', 'summary_contig'),
    ('summary.MAG.tsv', '.cache/summary/summary.MAG.tsv', 'summary_mag'),
)
MOVED = (
    ('01_mapping/mapping_manifest.tsv', '01_mapping/.cache/mapping_manifest.tsv'),
    ('10_summary/summary.contig.tsv', '.cache/summary/summary.contig.tsv'),
    ('10_summary/summary.MAG.tsv', '.cache/summary/summary.MAG.tsv'),
)

DROP = {
    'ROTATION_VALIDATED', 'REPEAT_ANALYZABLE', 'k', 'meryl_interpretation',
    'Assembly_FASTA', 'HiFi_all_BAM', 'HiFi_primary_BAM', 'HiFi_primary_supplementary_BAM',
    'Output_directory', 'Reads', 'MAG_FASTA', 'File', 'taxonomy_source', 'taxonomy_release',
    'final_taxonomic_chimera_policy', 'high_confidence_policy', 'depth_bam', 'depth_mapq_filter',
    'maximum_search_genes_per_side', 'target_classified_genes_per_side', 'minimum_classified_genes_per_side',
    'minimum_block_dominance', 'unclassified_gene_policy', 'final_call_ranks', 'hifi_anchor_bp',
    'hifi_min_support_qnames', 'event_id',
    'Large_insertion_sensitivity_min_bp', 'Large_insertion_sensitivity_min_QNAMEs',
    'Min_connected_QNAMEs_per_side',
    'MAG_integrated_QC_call', 'integrated_interpretation', 'excessive_repeat_call',
    'false_duplication_clipping_summary', 'FALSE_DUPLICATION_CLIPPING_SUMMARY',
    'Circularity_eligibility', 'SINGLE_CONTIG', 'CIRCULARITY_ANALYZABLE', 'HIFI_END_START_SUPPORTED',
    'CIRCULARITY_SUPPORTED', 'IN_FINAL_MAG_Y_N', 'MAG_level_input_assembler_circularity',
    'MAG_LEVEL_INPUT_ASSEMBLER_CIRCULARITY', 'HIFI_PRIMARY_SUPPLEMENTARY_FALSE_CLOSURE_Y_N',
    'HIFI_PRIMARY_SUPPLEMENTARY_FALSE_CLOSURE_BREAKPOINT_COUNT',
}
RENAME_BY_KIND = {
    'hsp': {'qseqid': 'Query_contig', 'sseqid': 'Subject_contig', 'length': 'Alignment_length_bp',
            'qstart': 'Query_start', 'qend': 'Query_end', 'sstart': 'Subject_start', 'send': 'Subject_end',
            'qlen': 'Query_length_bp', 'slen': 'Subject_length_bp'},
    'phantom_region': {'length_bp': 'Region_length_bp'},
    'false_dup': {'n_fully_evaluable_pairs': 'n_analyzable_pairs',
                  'n_false_duplication_candidates': 'n_false_duplication_pairs',
                  'n_no_false_duplication_signal_pairs': 'n_no_false_duplication_pairs',
                  'n_insufficient_evidence_pairs': 'n_not_analyzable_pairs'},
    'chimera_region': {'assessment_class': 'Structural_assessment',
                       'final_taxonomic_chimera_status': 'Taxonomic_assessment'},
    'chimera_contig': {'GTDB_assessment_status': 'Taxonomic_assessment'},
    'chimera_mag': {'GTDB_assessment_status': 'Taxonomic_assessment'},
    'mapping': {'Primary': 'Mapped_reads', 'Secondary': 'Secondary_alignments',
                'Supplementary': 'Supplementary_alignments'},
}
ASSEMBLER_METRIC_NAMES = {'Clipping_ge1_Y_sites': 'clipping_ge1_sites', 'Clipping_Y_sites': 'clipping_sites'}
RENAME_BY_KIND['clipping_assembler'] = {
    f'{prefix}{old}{suffix}': f'{prefix.capitalize() if prefix else "Terminal_to_internal_"}{new}{suffix}'
    for old, new in ASSEMBLER_METRIC_NAMES.items()
    for prefix, suffix in (('INTERNAL_', ''), ('INTERNAL_', '_per_Mbp'), ('TERMINAL_', ''),
                           ('TERMINAL_', '_per_Mbp'))
}
RENAME_BY_KIND['clipping_assembler'].update({
    f'terminal_to_internal_{old}_density_ratio': f'Terminal_to_internal_{new}_density_ratio'
    for old, new in ASSEMBLER_METRIC_NAMES.items()})

DROP_BY_KIND = {
    'clipping_assembler': {f'{prefix}validated_clipping_sites{suffix}'
                           for prefix in ('INTERNAL_', 'TERMINAL_')
                           for suffix in ('', '_per_Mbp')} | {'terminal_to_internal_validated_clipping_sites_density_ratio'},
    'mapping': {'Duplicates', 'Mapped', 'Primary_mapped', 'Mapped_percent', 'Primary_mapped_percent'},
    'clipping_site': {'default_candidate', 'valid_spanning_reads', 'support_class', 'evidence_class', 'validated_clipping',
                      'candidate_clipping', 'Clipping_analyzable', 'Clipping_analysis_reason'},
    'false_dup': {'FALSE_DUPLICATION', 'TANDEM_ORIENTATION', 'false_duplication_call',
                  'false_duplication_candidate_detected', 'false_duplication_assessment_call',
                  'structural_assessment_complete', 'structural_assessment_status'},
    'circularity_mag': {'Contig_length_bp'},
    'false_dup_pair': {'FALSE_DUPLICATION', 'TANDEM_ORIENTATION', 'LR_ANALYZABLE_Y_N'},
    'chimera_region': {'CHIMERA_EVIDENCE', 'final_taxonomic_chimera_evidence', 'HIFI_PRIMARY_SUPP_STRUCTURAL_CALL',
                       'default_clipping_support'},
    'chimera_contig': {'CHIMERA_EVIDENCE', 'n_taxonomic_shift_events', 'n_taxonomic_chimera_with_clipping_support_events',
                       'n_clipping_independent_taxonomic_shift_events', 'taxonomic_chimera_with_clipping_support_ranks',
                       'taxonomic_chimera_with_clipping_support_positions_1based',
                       'clipping_independent_taxonomic_shift_positions_1based',
                       'Taxonomic_chimera_with_clipping_support', 'Clipping_independent_taxonomic_shift',
                       'Clipping_independent_taxonomic_chimera'},
    'chimera_mag': {'CHIMERA_EVIDENCE', 'n_taxonomic_shift_events'},
    'phantom_region': {'raw_read_blast_hsps', 'best_hit_pident', 'best_hit_aligned_bp',
                       'raw_read_blast_query_covered_bp', 'raw_read_blast_query_coverage',
                       'phantom_sequence_call'},
}

RENAME = {
    'ID': 'MAG_ID', 'contig': 'Contig', 'Contig_ID': 'Contig', 'Source_assembler': 'Assembler',
    'source_assembler': 'Assembler', 'assembler': 'Assembler',
    'length': 'Length_bp', 'Contig_length_bp': 'Length_bp', 'length_bp': 'Length_bp',
    'total_length_bp': 'MAG_length_bp', 'Genome_length_bp': 'MAG_length_bp',
    'Number_of_contigs': 'n_contigs', 'MAG_contig_count': 'n_contigs', 'Contig_count': 'n_contigs',
    'MAG_CONTIG_COUNT': 'n_contigs', 'MAG_n_contigs': 'n_contigs',
    'Longest_exact_terminal_overlap_bp_up_to_100000': 'Terminal_overlap_bp',
    'SOURCE_CIRCULAR_CONTIG_COUNT': 'n_circular_contigs', 'SOURCE_CIRCULAR_CONTIGS': 'Circular_contigs',
    'SOURCE_CONTIG_ASSEMBLER_CIRCULARITY': 'Assembler_reported_circularity',
    'ASSEMBLER_REPORTED_CIRCULARITY': 'Assembler_reported_circularity',
    'SOURCE_CONTIG_CIRCULARITY_RESOLVED_COUNT': 'n_contig_circularity_resolved',
    'SOURCE_CONTIG_CIRCULARITY_UNRESOLVED_COUNT': 'n_contig_circularity_unresolved',
    'SOURCE_CONTIG_CIRCULARITY_COMPLETE_Y_N': 'Contig_circularity_complete',
    'Candidate_additional_target_contigs': 'Partner_contigs',
    'Candidate_additional_target_MAGs': 'Partner_MAGs',
    'GTDB_assessment_status': 'GTDB_assessment',
    'pos_1based': 'Position', 'anvio_pos_0based': 'Anvio_position_0based', 'pident': 'Identity_percent',
    'source_table': 'Anvio_table', 'SCG_HMM_source': 'SCG_HMM_collection', 'k': 'Kmer_size',
    'qname_clipping_ratio': 'Read_clipping_ratio', 'relative_pos': 'Relative_position', 'assembly_allele_support_call': 'Assembly_variant_class',
    'MAG_assembly_allele_support_call': 'Assembly_variant_class',
}
REPLACE = (
    ('HIFI_PRIMARY_SUPPLEMENTARY_', 'HiFi_'), ('HIFI_PRIMARY_SUPP_', 'HiFi_'), ('hifi_primary_supp_', 'HiFi_'),
    ('HIFI_ALL_MAPPED_', 'HiFi_all_'),
    ('HIFI_PRIMARY_', 'HiFi_primary_'), ('HIFI_SUPPLEMENTARY_', 'HiFi_supplementary_'),
    ('HIFI_SECONDARY_', 'HiFi_secondary_'), ('HIFI_', 'HiFi_all_'), ('hifi_', 'HiFi_'),
    ('_UNIQUE_READS', '_reads'), ('_QNAMEs', '_reads'), ('_qnames', '_reads'), ('_pos_', '_position_'),
    ('Paper_', ''), ('PAPER_', ''),
)
KEEP_UPPER = {'MAG', 'MAGs', 'ID', 'GTDB', 'SCG', 'HMM', 'END', 'START', 'HSP', 'HSPs',
              'SNV', 'GC', 'IDN', 'MAPQ', 'QC', 'I', 'D', 'N', 'HiFi', 'Mbp', 'bp'}


def plain_name(name: str) -> str:
    name = RENAME.get(name, name)
    for old, new in REPLACE:
        if old.startswith('_'):
            name = name.replace(old, new)
        elif name.startswith(old):
            name = new + name[len(old):]
    if name.startswith('LR_') and len(name) > 3:
        name = name[3:]
    name = re.sub(r'_Y_N(_NA)?$', '', name)
    tokens = []
    for token in name.split('_'):
        if token in KEEP_UPPER or not token.isupper() or len(token) == 1:
            tokens.append(token)
        else:
            tokens.append(token.lower())
    name = '_'.join(tokens)
    if not name.startswith('n_') and name[:1].islower():
        name = name[:1].upper() + name[1:]
    return RENAME.get(name, name)


VALUES = {
    'NOT_ANALYZABLE': 'NA', 'NOT_ASSESSABLE': 'NA', 'Not_assessable': 'NA',
    'NOT_EVALUABLE_CONTIG_END': 'NA', 'Not_applicable': 'NA', 'Not_applicable_non_MAG_contig': 'NA',
    'Not_applicable_no_candidate_pairs': 'NA', 'SUPPORTED': 'Y', 'UNSUPPORTED': 'N',
    'Y_single_contig': 'Y', 'N_not_single_contig': 'N',
    'Minority_assembly_variant': 'Minority', 'Unsupported_assembly_variant': 'Unsupported',
    'No_assembly_allele_disagreement': 'None',
}

LEAD_CONTIG = ('Contig', 'MAG_ID', 'Assembler')
LEAD_MAG = ('MAG_ID',)


def read(path):
    with open(path, encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        return list(reader.fieldnames or []), list(reader)


def write(path, fields, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(f'.{path.name}.tmp.{os.getpid()}')
    with open(temp, 'w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n',
                                extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temp, path)


def column_plan(fields, kind):
    """Pairs of (internal column, output column) in output order."""
    drop = DROP | DROP_BY_KIND.get(kind, set())
    plan, used = [], set()
    for field in fields:
        if field in drop:
            continue
        name = RENAME_BY_KIND.get(kind, {}).get(field) or plain_name(field)
        if name in used:
            name = field
        used.add(name)
        plan.append((field, name))
    lead = LEAD_MAG if kind in ('mag', 'circularity_mag', 'chimera_mag') else LEAD_CONTIG
    first = [item for key in lead for item in plan if item[1] == key]
    return first + [item for item in plan if item not in first]


def generic(internal, public, kind):
    public = Path(public)
    public.parent.mkdir(parents=True, exist_ok=True)
    temp = public.with_name(f'.{public.name}.tmp.{os.getpid()}')
    with open(internal, encoding='utf-8-sig', newline='') as source, \
            open(temp, 'w', encoding='utf-8', newline='') as target:
        reader = csv.reader(source, delimiter='\t')
        fields = next(reader, [])
        plan = column_plan(fields, kind)
        index = {field: i for i, field in enumerate(fields)}
        columns = [index[field] for field, _name in plan]
        writer = csv.writer(target, delimiter='\t', lineterminator='\n')
        writer.writerow([name for _field, name in plan])
        for row in reader:
            writer.writerow([VALUES.get(row[i], row[i]) if i < len(row) else '' for i in columns])
    os.replace(temp, public)


NOT_RUN = {'Disabled', 'Not_run'}


def to_yn(value, yes=(), no=()):
    if value in yes or value == 'Y':
        return 'Y'
    if value in no or value == 'N' or value.startswith('No_'):
        return 'N'
    return 'NA'


CONTIG_CALLS = (
    ('Chimera', 'Taxonomic_chimera', (), ()),
    ('Assembly_variant', 'assembly_allele_support_call',
     ('Unsupported_assembly_variant', 'Minority_assembly_variant'), ('Assembly_allele_supported',)),
    ('Phantom', 'phantom_sequence_call', ('Phantom_sequence_candidate',), ()),
    ('Excessive_repeat', 'excessive_repeat_call', ('Excessive_repeat_candidate',), ()),
    ('False_duplication', 'false_duplication_call', ('False_duplication_candidate',), ()),
)


def not_run_modules(row, columns):
    return [label for label, column in columns if row.get(column) in NOT_RUN]


def summary_contig(internal, public):
    _fields, rows = read(internal)
    out = []
    for row in rows:
        clip = row['validated_clipping_breakpoints']
        new = {
            'Contig': row['contig'], 'MAG_ID': row['MAG_ID'], 'Length_bp': row['length'],
            'n_flagged_checks': row['n_flagged_modules'],
            'Clipping': 'NA' if clip == 'NA' else ('Y' if int(clip) > 0 else 'N'),
            'n_clipping_sites': clip,
        }
        for name, column, yes, no in CONTIG_CALLS:
            new[name] = to_yn(row[column], yes, no)
        new['Chimera_ranks'] = row['GTDB_taxonomic_chimera_ranks']
        new['Repeat_fraction'] = row['repeat_fraction']
        new['n_duplicated_pairs'] = row['n_false_duplication_candidate_pairs']
        new['n_false_duplication_pairs'] = row['n_false_duplication_candidates']
        new['n_unresolved_duplicated_pairs'] = row['n_false_duplication_insufficient_pairs']
        new['False_duplication_with_clipping'] = row['false_duplication_with_validated_clipping_Y_N_NA']
        missing = not_run_modules(row, (
            ('taxonomic_chimera', 'Taxonomic_chimera'), ('assembly_variant', 'assembly_allele_support_call'),
            ('phantom', 'phantom_sequence_call'), ('excessive_repeat', 'excessive_repeat_call'),
            ('false_duplication', 'false_duplication_call')))
        if clip == 'NA':
            missing.append('clipping')
        new['Checks_not_run'] = ','.join(sorted(set(missing))) or 'NA'
        new['n_checks_not_run'] = str(len(set(missing)))
        out.append(new)
    write(public, SUMMARY_CONTIG_FIELDS, out)
    return {row['contig']: row for row in rows}, out


# The public check count includes only the diagnostic flags in summary.MAG.tsv.
# False_closure remains available in summary_detail.MAG.tsv and in the
# false-closure module, because it is not synonymous with premature circularization.
MAG_CHECKS = (
    ('Clipping', 'Clipping_Y_N', None),
    ('Premature_circularization', 'Premature_circularization_Y_N', None),
    ('Phantom', 'Phantom_sequence_Y_N', 'phantom_sequence_call'),
    ('Excessive_repeat', 'Excessive_repeat_Y_N', 'excessive_repeat_call'),
    ('False_duplication', 'False_duplication_Y_N', 'false_duplication_call'),
    ('Assembly_variant', 'Assembly_allele_disagreement_Y_N', 'assembly_allele_support_call'),
    ('Chimera', 'Taxonomic_chimera_Y_N', 'Taxonomic_chimera'),
)


def circularity_summary_calls(root):
    """Use existing circularity evidence; never infer read support from a final flag."""
    root = Path(root)
    for path in (root / '02_circularity/circularity.MAG.tsv',
                 root / '02_circularity/.cache/circularity.MAG.tsv'):
        if not path.is_file():
            continue
        fields, rows = read(path)
        key = 'MAG_ID' if 'MAG_ID' in fields else 'ID' if 'ID' in fields else None
        if not key or not {'Assembler_reported_circularity',
                           'Read_supported_circular_topology',
                           'High_confidence_circular'} <= set(fields):
            continue
        calls = {}
        for row in rows:
            mag = row[key]
            if mag in calls:
                raise ValueError(f'Duplicate circularity MAG ID: {mag}')
            assembler = row['Assembler_reported_circularity']
            read_support = row['Read_supported_circular_topology']
            both = row['High_confidence_circular']
            if any(v not in ('Y', 'N', 'NA') for v in (assembler, read_support, both)):
                raise ValueError(f'Invalid circularity flag for MAG {mag}')
            calls[mag] = (assembler, read_support, both)
        return calls
    return {}


def summary_mag(internal, public, contigs):
    _fields, rows = read(internal)
    by_mag = {}
    for contig in contigs.values():
        by_mag.setdefault(contig['MAG_ID'], []).append(contig)
    circularity = circularity_summary_calls(Path(public).parent)
    out = []
    for row in rows:
        members = by_mag.get(row['MAG_ID'], [])
        new = {'MAG_ID': row['MAG_ID'], 'n_contigs': row['n_contigs']}
        not_run = []
        for name, column, contig_column in MAG_CHECKS:
            value = row[column]
            ran = True
            if contig_column and members and all(m.get(contig_column) in NOT_RUN for m in members):
                ran = False
            if name == 'Clipping' and members and all(m.get('validated_clipping_breakpoints') == 'NA'
                                                      for m in members):
                ran = False
            if name == 'Premature_circularization' and row['False_closure_evidence'] in NOT_RUN:
                ran = False
            if not ran:
                not_run.append(name.lower())
            new[name] = value if ran else 'NA'
        new['n_flagged_checks'] = str(sum(new[name] == 'Y' for name, _c, _r in MAG_CHECKS))
        assembly, read_support, both = circularity.get(row['MAG_ID'], ('NA', 'NA', 'NA'))
        new.update({
            'Assembler_reported_circularity': assembly,
            'Read_supported_circularity': read_support,
            'High_confidence_circular': both,
        })
        if row['Circular'] in NOT_RUN:
            not_run.append('circularity')
        new.update({
            'False_closure': row['False_closure_Y_N'] if row['False_closure_Y_N'] in ('Y', 'N') else 'NA',
            'Small_genome_proxy': row['Paper_small_genome_proxy'] if row['Paper_small_genome_proxy'] in ('Y', 'N') else 'NA',
            'Premature_circularization_screen': row['Premature_circularization_screen']
            if row['Premature_circularization_screen'] in ('Y', 'N') else 'NA',
            'False_closure_partner_contigs': row['False_closure_candidate_additional_contigs']
            if row['False_closure_candidate_additional_contigs'] not in NOT_RUN | {'Not_target'} else 'NA',
            'n_clipping_sites': row['n_validated_clipping_breakpoints'],
            'n_chimera_contigs': row['n_taxonomic_chimera_flags'],
            'n_variant_contigs': row['n_assembly_allele_disagreement_flags'],
            'n_phantom_contigs': row['n_phantom_sequence_flags'],
            'Repeat_fraction': row['repeat_fraction'],
            'n_excessive_repeat_contigs': row['n_excessive_repeat_contig_flags'],
            'n_false_duplication_contigs': row['n_false_duplication_flags'],
            'n_duplicated_pairs': row['n_false_duplication_candidate_pairs'],
            'n_unresolved_duplicated_pairs': row['n_false_duplication_insufficient_pairs'],
            'False_duplication_with_clipping': row['false_duplication_with_validated_clipping_Y_N_NA'],
            'Checks_not_run': ','.join(sorted(set(not_run))) or 'NA',
            'n_checks_not_run': str(len(set(not_run))),
        })
        out.append(new)
    write(public, SUMMARY_MAG_FIELDS, out)
    return out


# Public flags follow end-to-end module order. Counters are the last two columns.
SUMMARY_CONTIG_FIELDS = [
    'Contig', 'MAG_ID', 'Length_bp', 'Clipping', 'Phantom', 'Excessive_repeat',
    'False_duplication', 'Assembly_variant', 'Chimera',
    'n_flagged_checks', 'n_checks_not_run',
]
SUMMARY_CONTIG_DETAILS = [
    'Assembler', 'n_clipping_sites', 'Repeat_fraction', 'n_duplicated_pairs',
    'n_false_duplication_pairs', 'n_unresolved_duplicated_pairs', 'Chimera_ranks',
    'False_duplication_with_clipping', 'Checks_not_run',
]
SUMMARY_MAG_FIELDS = [
    'MAG_ID', 'n_contigs',
    'Assembler_reported_circularity', 'Read_supported_circularity',
    'High_confidence_circular', 'Clipping', 'Premature_circularization',
    'Phantom', 'Excessive_repeat', 'False_duplication',
    'Assembly_variant', 'Chimera', 'n_flagged_checks', 'n_checks_not_run',
]
SUMMARY_MAG_DETAILS = [
    'False_closure', 'Small_genome_proxy', 'Premature_circularization_screen',
    'False_duplication_with_clipping', 'Checks_not_run',
    'n_clipping_sites', 'False_closure_partner_contigs', 'n_phantom_contigs',
    'Repeat_fraction', 'n_excessive_repeat_contigs',
    'n_false_duplication_contigs', 'n_duplicated_pairs',
    'n_unresolved_duplicated_pairs', 'n_variant_contigs', 'n_chimera_contigs',
]
FULL_MODULES = (
    ('circularity', '02_circularity/circularity'),
    ('clipping', '03_clipping/clipping'),
    ('false_closure', '04_false_closure/false_closure'),
    ('phantom', '05_phantom/phantom'),
    ('repeats', '06_repeats/repeats'),
    ('false_duplication', '08_false_duplication/false_duplication'),
    ('variants', '09_assembly_variants/variants'),
    ('chimera', '07_chimera/chimera'),
)
FULL_SKIP = {'MAG': {'MAG_ID', 'n_contigs', 'MAG_length_bp', 'Assembler'},
             'contig': {'Contig', 'MAG_ID', 'Assembler', 'Length_bp', 'n_contigs'}}


def full_table(root, level, rows, key, lead):
    """Summary rows followed by every column of each module table of this level, named MODULE.COLUMN."""
    root = Path(root)
    wanted = {row[key] for row in rows}
    fields = list(lead)
    merged = {name: {} for name in wanted}
    for module, stem in FULL_MODULES:
        path = root / f'{stem}.{level}.tsv'
        if not path.is_file():
            continue
        with open(path, encoding='utf-8', newline='') as handle:
            reader = csv.reader(handle, delimiter='\t')
            header = next(reader, [])
            if key not in header:
                continue
            index = header.index(key)
            assembler = header.index('Assembler') if 'Assembler' in header else None
            keep = [(i, f'{module}.{name}') for i, name in enumerate(header) if name not in FULL_SKIP[level]]
            fields += [name for _i, name in keep]
            for values in reader:
                target = merged.get(values[index] if index < len(values) else None)
                if target is not None:
                    if assembler is not None and assembler < len(values) and values[assembler] not in ('', 'NA'):
                        target.setdefault('Assembler', values[assembler])
                    for i, name in keep:
                        target[name] = values[i] if i < len(values) else ''
    out = []
    for row in rows:
        extra = merged.get(row[key], {})
        full = {name: row.get(name, extra.get(name, 'NA')) for name in lead}
        for name in fields[len(lead):]:
            full[name] = extra.get(name, 'NA')
        out.append(full)
    return fields, out


def write_full_tables(root, contig_rows, mag_rows):
    root = Path(root)
    mag_contigs = [row for row in contig_rows if row.get('MAG_ID') not in ('', 'NA', None)]
    fields, rows = full_table(root, 'contig', mag_contigs, 'Contig',
                              SUMMARY_CONTIG_FIELDS[:3] + ['Assembler'] + SUMMARY_CONTIG_FIELDS[3:]
                              + SUMMARY_CONTIG_DETAILS[1:])
    write(root / 'summary_detail.contig.tsv', fields, rows)
    fields, rows = full_table(root, 'MAG', mag_rows, 'MAG_ID', SUMMARY_MAG_FIELDS + SUMMARY_MAG_DETAILS)
    write(root / 'summary_detail.MAG.tsv', fields, rows)
    return [str(root / 'summary_detail.contig.tsv'), str(root / 'summary_detail.MAG.tsv')]


def header_matches(path, expected):
    """Check the public table's column order without reading its data rows."""
    path = Path(path)
    if not path.is_file():
        return False
    with open(path, encoding='utf-8-sig', newline='') as handle:
        header = next(csv.reader(handle, delimiter='\t'), [])
    return header[:len(expected)] == expected


def summary_headers_current(root):
    root = Path(root)
    return (header_matches(root / 'summary.contig.tsv', SUMMARY_CONTIG_FIELDS)
            and header_matches(root / 'summary.MAG.tsv', SUMMARY_MAG_FIELDS))


def full_tables_current(root):
    root = Path(root)
    targets = [root / 'summary_detail.contig.tsv', root / 'summary_detail.MAG.tsv']
    if not all(path.is_file() for path in targets):
        return False
    detailed_contig_prefix = (SUMMARY_CONTIG_FIELDS[:3] + ['Assembler']
                              + SUMMARY_CONTIG_FIELDS[3:] + SUMMARY_CONTIG_DETAILS[1:])
    if not (header_matches(targets[0], detailed_contig_prefix)
            and header_matches(targets[1], SUMMARY_MAG_FIELDS + SUMMARY_MAG_DETAILS)):
        return False
    oldest = min(path.stat().st_mtime_ns for path in targets)
    sources = [root / '.cache/summary/summary.contig.tsv', root / '.cache/summary/summary.MAG.tsv']
    sources += [root / f'{stem}.{level}.tsv' for _m, stem in FULL_MODULES for level in ('contig', 'MAG')]
    return all(path.stat().st_mtime_ns <= oldest for path in sources if path.is_file())


def migrate(root):
    root = Path(root)
    # Migrate older public names without leaving duplicate summary tables.
    for level in ('contig', 'MAG'):
        old = root / f'summary_full.{level}.tsv'
        new = root / f'summary_detail.{level}.tsv'
        if old.is_file():
            if not new.exists():
                os.replace(old, new)
            else:
                old.unlink()
    for old, internal in MOVED:
        if (root / old).is_file() and not (root / internal).exists():
            (root / internal).parent.mkdir(parents=True, exist_ok=True)
            os.replace(root / old, root / internal)
    for public, internal, _kind in TABLES:
        if (root / public).is_file() and not (root / internal).exists():
            (root / internal).parent.mkdir(parents=True, exist_ok=True)
            os.replace(root / public, root / internal)
    old_summary = root / '10_summary'
    if old_summary.is_dir() and not any(old_summary.iterdir()):
        old_summary.rmdir()


def publish(root, force=False):
    root = Path(root)
    written = []
    contigs = None
    summary_internal = (root / '.cache/summary/summary.contig.tsv', root / '.cache/summary/summary.MAG.tsv')
    force_summary = None
    contig_rows = mag_rows = None
    for public, internal, kind in TABLES:
        source, target = root / internal, root / public
        if kind is None:
            if target.is_file() and source.is_file():
                target.unlink()
            continue
        if not source.is_file():
            continue
        if kind in ('summary_contig', 'summary_mag') and force_summary is None:
            force_summary = (all(path.is_file() for path in summary_internal)
                             and (not full_tables_current(root) or not summary_headers_current(root)))
        if kind == 'summary_mag' and contigs is None:
            contig_internal = root / '.cache/summary/summary.contig.tsv'
            contigs = {r['contig']: r for r in read(contig_internal)[1]} if contig_internal.is_file() else {}
        current = target.is_file() and target.stat().st_mtime_ns >= source.stat().st_mtime_ns
        if current and not force and not (force_summary and kind in ('summary_contig', 'summary_mag')):
            continue
        if kind == 'summary_contig':
            contigs, contig_rows = summary_contig(source, target)
        elif kind == 'summary_mag':
            mag_rows = summary_mag(source, target, contigs or {})
        else:
            generic(source, target, kind)
        written.append(str(target))
    if contig_rows is not None and mag_rows is not None:
        written += write_full_tables(root, contig_rows, mag_rows)
    return written


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) < 2 or argv[0] not in ('migrate', 'publish'):
        sys.exit(__doc__)
    if argv[0] == 'migrate':
        migrate(argv[1])
    else:
        for path in publish(argv[1], force='--force' in argv):
            print(path)
    return 0


if __name__ == '__main__':
    sys.exit(main())
