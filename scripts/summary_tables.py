#!/usr/bin/env python3
"""Integrated QC summary: build the MAG and contig tables, put the overall call first and validate them.

Subcommands:
  build              Build summary.MAG.tsv and summary.contig.tsv from the module tables
  arrange            Put the overall call next to the MAG ID and drop duplicated columns
  validate           Check that the MAG and contig summaries agree
"""

import csv
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path


def available(enabled, path):
    return enabled == 'Y' and os.path.isfile(path) and os.path.getsize(path) > 0


def cmd_build(argv):
    """Build summary.MAG.tsv and summary.contig.tsv from the module tables"""
    argv = argv

    if len(argv) != 25:
        raise SystemExit(f'ERROR: expected 25 split-summary arguments, found {len(argv)}')

    (
        fai, clipping_file, gtdb_chimera_file, allele_file,
        phantom_file, repeat_file, repeat_mag_file, false_dup_file,
        false_dup_clipping_contig_file, false_dup_clipping_mag_file, mag_map_file,
        circularity_file, circularity_audit_file,
        false_closure_file, false_closure_audit_file, out_contig, out_mag,
        run_clipping, run_taxonomic_chimera, run_allele, run_phantom, run_repeat,
        run_false_dup, run_circularity, run_false_closure,
    ) = argv

    expected_contig_modules = 5 + (1 if run_false_dup == 'Y' else 0)

    module_available = {
        'clipping': available(run_clipping, clipping_file),
        'gtdb_chimera': available(run_taxonomic_chimera, gtdb_chimera_file),
        'assembly_variants': available(run_allele, allele_file),
        'phantom': available(run_phantom, phantom_file),
        'repeat': available(run_repeat, repeat_file) and available(run_repeat, repeat_mag_file),
        'false_duplication': available(run_false_dup, false_dup_file),
        'false_duplication_clipping': (
            available(run_false_dup, false_dup_clipping_contig_file)
            and available(run_false_dup, false_dup_clipping_mag_file)
        ),
        'circularity': (
            available(run_circularity, circularity_file)
            and available(run_circularity, circularity_audit_file)
        ),
        'false_closure': (
            available(run_false_closure, false_closure_file)
            and available(run_false_closure, false_closure_audit_file)
        ),
    }

    def unavailable_label(name, enabled):
        if enabled != 'Y':
            return 'Disabled'
        return 'Not_run' if not module_available[name] else None

    lengths = {}

    with open(fai) as handle:
        for line in handle:
            x = line.rstrip('\n').split('\t')
            if len(x) >= 2:
                if x[0] in lengths:
                    raise SystemExit(f'ERROR: duplicate FASTA contig ID: {x[0]}')
                lengths[x[0]] = int(x[1])

    contig_to_mag = {}

    all_mags = set()

    if mag_map_file not in {'', 'NA', 'None'}:
        if not os.path.isfile(mag_map_file) or os.path.getsize(mag_map_file) == 0:
            raise SystemExit(f'ERROR: configured MAG_CONTIG_MAP is missing: {mag_map_file}')
        with open(mag_map_file) as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            fields = reader.fieldnames or []
            lower = {name.lower(): name for name in fields}
            mag_col = lower.get('mag_id') or lower.get('id')
            contig_col = lower.get('contig') or lower.get('original_contig')
            if not mag_col or not contig_col:
                raise SystemExit('ERROR: MAG_CONTIG_MAP needs MAG_ID and Contig columns')
            seen_map_rows = set()
            for row in reader:
                mag, contig = row[mag_col].strip(), row[contig_col].strip()
                if not mag or not contig or (mag, contig) in seen_map_rows:
                    raise SystemExit(f'ERROR: blank/duplicate MAG map row: {mag}/{contig}')
                seen_map_rows.add((mag, contig))
                if contig not in lengths:
                    raise SystemExit(f'ERROR: MAG map contig absent from assembly: {contig}')
                if contig in contig_to_mag and contig_to_mag[contig] != mag:
                    raise SystemExit(f'ERROR: contig assigned to multiple MAGs: {contig}')
                contig_to_mag[contig] = mag
                all_mags.add(mag)

    clipping_counts = defaultdict(int)

    if module_available['clipping']:
        with open(clipping_file) as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            required = {
                'contig', 'length', 'pos_1based', 'coverage', 'clipping_count',
                'clipping_ratio',
            }
            fields = set(reader.fieldnames or [])
            missing = required - fields
            if missing:
                raise SystemExit(
                    f'ERROR: clipping table lacks columns {sorted(missing)}: '
                    f'{clipping_file}'
                )
            if 'Clipping' in fields:
                call_column = 'Clipping'
            elif 'validated_clipping_breakpoint' in fields:
                call_column = 'validated_clipping_breakpoint'
            else:
                call_column = 'paper_high_confidence_clipping'
            if call_column not in fields:
                raise SystemExit('ERROR: clipping table lacks a Clipping call column')
            seen_sites = set()
            for row in reader:
                contig = row['contig']
                if contig not in lengths or int(row['length']) != lengths[contig]:
                    raise SystemExit(f'ERROR: stale clipping contig/length: {contig}')
                pos = int(row['pos_1based'])
                if pos < 1 or pos > lengths[contig] or (contig, pos) in seen_sites:
                    raise SystemExit(f'ERROR: invalid/duplicate clipping site: {contig}:{pos}')
                seen_sites.add((contig, pos))
                coverage = int(row['coverage'])
                clipped = int(row['clipping_count'])
                ratio = float(row['clipping_ratio'])
                expected = clipped / coverage if coverage else 0.0
                if coverage < 0 or clipped < 0 or clipped > coverage:
                    raise SystemExit(f'ERROR: invalid clipping counts: {contig}:{pos}')
                if abs(expected - ratio) > 1e-6:
                    raise SystemExit(f'ERROR: inconsistent clipping ratio: {contig}:{pos}')
                high = row[call_column].strip().upper()
                if high not in {'Y', 'N', 'YES', 'NO'}:
                    raise SystemExit(f'ERROR: invalid clipping call: {contig}:{pos}')
                if high in {'Y', 'YES'}:
                    clipping_counts[contig] += 1

    def read_call_table(path, key_col, call_col, allowed_calls):
        result = {}
        with open(path) as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            fields = set(reader.fieldnames or [])
            required = {key_col, call_col, 'length'}
            if not required.issubset(fields):
                raise SystemExit(f'ERROR: {path} lacks columns {sorted(required - fields)}')
            for row in reader:
                key = row[key_col]
                call = row[call_col]
                if key in result:
                    raise SystemExit(f'ERROR: duplicate contig result in {path}: {key}')
                if key not in lengths or int(row['length']) != lengths[key]:
                    raise SystemExit(f'ERROR: stale name/length in {path}: {key}')
                if call not in allowed_calls:
                    raise SystemExit(f'ERROR: unrecognized call in {path}: {key}/{call}')
                result[key] = call
        return result

    def read_gtdb_chimera(path):
        result = {}
        with open(path) as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            fields = set(reader.fieldnames or [])
            required = {
                'taxonomy_source', 'taxonomy_release', 'contig', 'length',
                'n_validated_clipping_sites', 'taxonomic_shift_ranks',
                'Taxonomic_shift', 'Taxonomic_chimera_with_clipping_support',
                'Clipping_independent_taxonomic_shift', 'GTDB_assessment_status',
                'Structural_chimera', 'High_confidence_structural_chimera',
                'Taxonomic_chimera', 'CHIMERA_EVIDENCE',
            }
            missing = required - fields
            if missing:
                raise SystemExit(
                    f'ERROR: GTDB taxonomic-chimera table lacks current columns '
                    f'{sorted(missing)}: {path}. Rerun 07_chimera.sh.'
                )
            release = None
            for row in reader:
                contig = row['contig']
                if contig in result:
                    raise SystemExit(f'ERROR: duplicate GTDB contig result: {contig}')
                if contig not in lengths or int(row['length']) != lengths[contig]:
                    raise SystemExit(f'ERROR: stale GTDB name/length: {contig}')
                if row['taxonomy_source'].strip().upper() != 'GTDB':
                    raise SystemExit(f'ERROR: non-GTDB taxonomy: {contig}')
                current_release = row['taxonomy_release'].strip()
                if not current_release:
                    raise SystemExit(f'ERROR: blank GTDB release: {contig}')
                if release is None:
                    release = current_release
                elif current_release != release:
                    raise SystemExit('ERROR: mixed GTDB releases in chimera results')
                row = dict(row)
                for flag in (
                    'Taxonomic_shift', 'Taxonomic_chimera_with_clipping_support',
                    'Clipping_independent_taxonomic_shift', 'Structural_chimera',
                    'High_confidence_structural_chimera', 'Taxonomic_chimera',
                ):
                    value = row[flag].strip().upper()
                    if value not in {'Y', 'N', 'NA'}:
                        raise SystemExit(f'ERROR: invalid GTDB {flag}: {contig}/{value}')
                    row[flag] = value
                if row['Taxonomic_chimera_with_clipping_support'] != row['Taxonomic_chimera']:
                    raise SystemExit(
                        f'ERROR: final Taxonomic_chimera must equal the validated-clipping-supported call: {contig}'
                    )
                allowed_evidence = {
                    'TAXONOMIC_SHIFT_WITH_STRINGENT_VALIDATED_CLIPPING',
                    'TAXONOMIC_SHIFT_WITH_VALIDATED_CLIPPING',
                    'TAXONOMIC_SHIFT_WITHOUT_VALIDATED_CLIPPING',
                    'NO_TAXONOMIC_SHIFT', 'NOT_ASSESSABLE',
                }
                evidence = row['CHIMERA_EVIDENCE'].strip()
                if evidence not in allowed_evidence:
                    raise SystemExit(
                        f'ERROR: invalid GTDB CHIMERA_EVIDENCE: {contig}/{evidence}'
                    )
                final_call = row['Taxonomic_chimera']
                high_call = row['High_confidence_taxonomic_chimera'].strip().upper() if 'High_confidence_taxonomic_chimera' in row else 'NA'
                if final_call == 'Y' and evidence not in {
                    'TAXONOMIC_SHIFT_WITH_STRINGENT_VALIDATED_CLIPPING',
                    'TAXONOMIC_SHIFT_WITH_VALIDATED_CLIPPING',
                }:
                    raise SystemExit(f'ERROR: positive taxonomic chimera has inconsistent CHIMERA_EVIDENCE: {contig}')
                if final_call == 'N' and row['Taxonomic_shift'] == 'Y' and evidence != 'TAXONOMIC_SHIFT_WITHOUT_VALIDATED_CLIPPING':
                    raise SystemExit(f'ERROR: taxonomic-shift-only contig has inconsistent CHIMERA_EVIDENCE: {contig}')
                if row['Taxonomic_shift'] == 'N' and evidence != 'NO_TAXONOMIC_SHIFT':
                    raise SystemExit(f'ERROR: negative GTDB shift has inconsistent CHIMERA_EVIDENCE: {contig}')
                if row['Taxonomic_shift'] == 'NA' and evidence != 'NOT_ASSESSABLE':
                    raise SystemExit(f'ERROR: non-assessable GTDB shift has inconsistent CHIMERA_EVIDENCE: {contig}')
                row['CHIMERA_EVIDENCE'] = evidence
                row['validated_clipping_breakpoints'] = row['n_validated_clipping_sites']
                result[contig] = row
        expected = set(lengths)
        if set(result) != expected:
            missing = sorted(expected - set(result))[:5]
            extra = sorted(set(result) - expected)[:5]
            raise SystemExit(
                f'ERROR: GTDB taxonomic-chimera table does not match all assembly contigs; '
                f'missing={missing}, extra={extra}'
            )
        return result

    if module_available['gtdb_chimera']:
        gtdb_chimera = read_gtdb_chimera(gtdb_chimera_file)
        gtdb_chimera_default = 'Missing_contig_result'
    else:
        gtdb_chimera = {}
        gtdb_chimera_default = unavailable_label(
            'gtdb_chimera', run_taxonomic_chimera
        )

    allele_calls = {
        'Unsupported_assembly_variant', 'Minority_assembly_variant',
        'No_assembly_allele_disagreement',
    }

    if module_available['assembly_variants']:
        allele = read_call_table(
            allele_file, 'contig', 'assembly_allele_support_call', allele_calls
        )
        expected = set(contig_to_mag)
        if set(allele) != expected:
            missing = sorted(expected - set(allele))[:5]
            extra = sorted(set(allele) - expected)[:5]
            raise SystemExit(
                'ERROR: variant table does not match MAG contigs; '
                f'missing={missing}, extra={extra}'
            )
        allele_default = 'Not_applicable_non_MAG_contig'
    else:
        allele = {}
        allele_default = unavailable_label('assembly_variants', run_allele)

    repeat = {}

    repeat_frac = {}

    if module_available['repeat']:
        with open(repeat_file) as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            fields = set(reader.fieldnames or [])
            base_required = {'contig', 'length'}
            if not base_required.issubset(fields):
                raise SystemExit(f'ERROR: malformed repeat table: {repeat_file}')
            new_schema = {'intra_repeat_fraction', 'EXCESSIVE_REPEAT_CONTIG'}.issubset(fields)
            legacy_schema = {'repeat_fraction', 'excessive_repeat_call'}.issubset(fields)
            if not (new_schema or legacy_schema):
                raise SystemExit(
                    'ERROR: repeat table lacks either current or legacy repeat-call columns: '
                    f'{repeat_file}'
                )
            for row in reader:
                contig = row['contig']
                if contig in repeat:
                    raise SystemExit(f'ERROR: duplicate repeat result: {contig}')
                if contig not in lengths or int(row['length']) != lengths[contig]:
                    raise SystemExit(f'ERROR: stale repeat name/length: {contig}')
                if new_schema:
                    fraction_text = row['intra_repeat_fraction']
                    analyzable = row.get('REPEAT_ANALYZABLE', 'Y')
                    flag = row['EXCESSIVE_REPEAT_CONTIG']
                    if analyzable == 'N':
                        call = 'Repeat_analysis_skipped'
                    elif analyzable != 'Y':
                        raise SystemExit(f'ERROR: invalid REPEAT_ANALYZABLE value: {contig}: {analyzable}')
                    elif flag == 'Y':
                        call = 'Excessive_repeat_candidate'
                    elif flag == 'N':
                        call = 'No_excessive_repeat_signal'
                    else:
                        raise SystemExit(f'ERROR: invalid EXCESSIVE_REPEAT_CONTIG value: {contig}: {flag}')
                else:
                    fraction_text = row['repeat_fraction']
                    call = row['excessive_repeat_call']
                    if call not in {
                        'Excessive_repeat_candidate', 'No_excessive_repeat_signal',
                        'Repeat_analysis_skipped',
                    }:
                        raise SystemExit(f'ERROR: unrecognized repeat call: {contig}')
                fraction = float(fraction_text)
                if not 0 <= fraction <= 1:
                    raise SystemExit(f'ERROR: invalid repeat fraction: {contig}')
                repeat[contig] = call
                repeat_frac[contig] = fraction_text
        if set(repeat) != set(lengths):
            raise SystemExit('ERROR: repeat table does not contain every assembly contig')

        repeat_default = 'Missing_contig_result'
    else:
        repeat_default = unavailable_label('repeat', run_repeat)

    mag_repeat = {}

    if module_available['repeat']:
        with open(repeat_mag_file) as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            fields = set(reader.fieldnames or [])
            base_required = {'MAG_ID', 'n_contigs', 'total_length_bp'}
            missing = base_required - fields
            if missing:
                raise SystemExit(
                    f'ERROR: MAG repeat table lacks columns {sorted(missing)}: '
                    f'{repeat_mag_file}'
                )
            new_schema = {'total_repeat_fraction', 'EXCESSIVE_REPEAT_MAG'}.issubset(fields)
            legacy_schema = {'repeat_fraction', 'excessive_repeat_call'}.issubset(fields)
            if not (new_schema or legacy_schema):
                raise SystemExit(
                    'ERROR: MAG repeat table lacks either current or legacy repeat-call columns: '
                    f'{repeat_mag_file}'
                )
            for row in reader:
                mag = row['MAG_ID']
                if not mag or mag in mag_repeat:
                    raise SystemExit(f'ERROR: blank or duplicate MAG repeat ID: {mag!r}')
                if new_schema:
                    fraction_text = row['total_repeat_fraction']
                    analyzable = row.get('REPEAT_ANALYZABLE', 'Y')
                    flag = row['EXCESSIVE_REPEAT_MAG']
                    if analyzable != 'Y':
                        raise SystemExit(f'ERROR: invalid MAG REPEAT_ANALYZABLE value: {mag}: {analyzable}')
                    if flag == 'Y':
                        call = 'Excessive_repeat_candidate'
                    elif flag == 'N':
                        call = 'No_excessive_repeat_signal'
                    else:
                        raise SystemExit(f'ERROR: invalid EXCESSIVE_REPEAT_MAG value: {mag}: {flag}')
                else:
                    fraction_text = row['repeat_fraction']
                    call = row['excessive_repeat_call']
                    if call not in {'Excessive_repeat_candidate', 'No_excessive_repeat_signal'}:
                        raise SystemExit(f'ERROR: invalid MAG repeat call: {mag}: {call}')
                fraction = float(fraction_text)
                if not 0 <= fraction <= 1:
                    raise SystemExit(f'ERROR: invalid MAG repeat fraction: {mag}')
                mag_repeat[mag] = {
                    'fraction': fraction_text,
                    'call': call,
                    'n_contigs': int(row['n_contigs']),
                }

    false_dup = {}

    false_dup_n_pairs = {}

    false_dup_n_fully = {}

    false_dup_n_positive = {}

    false_dup_n_insufficient = {}

    false_dup_detection = {}

    false_dup_assessment_call = {}

    false_dup_assessment_status = {}

    false_dup_assessment_complete = {}

    if module_available['false_duplication']:
        with open(false_dup_file) as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            fields = set(reader.fieldnames or [])
            required = {
                'contig', 'length', 'n_duplicated_pairs',
                'n_fully_evaluable_pairs',
                'n_false_duplication_candidates',
                'n_no_false_duplication_signal_pairs',
                'n_insufficient_evidence_pairs', 'structural_assessment_status',
                'structural_assessment_complete',
                'false_duplication_candidate_detected',
                'false_duplication_assessment_call', 'false_duplication_call',
            }
            if not required.issubset(fields):
                raise SystemExit(
                    'ERROR: false-duplication table uses an old or malformed schema; '
                    'rerun 08_false_duplication.sh'
                )
            allowed_calls = {
                'False_duplication_candidate',
                'No_duplicated_structure_detected',
                'No_false_duplication_signal',
                'No_false_duplication_candidate_detected_assessment_incomplete',
            }
            allowed_status = {
                'Not_applicable_no_candidate_pairs', 'Complete', 'Partial',
                'Insufficient_for_all_candidate_pairs',
            }
            for row in reader:
                contig = row['contig']
                if contig in false_dup:
                    raise SystemExit(f'ERROR: duplicate false-duplication result: {contig}')
                if contig not in lengths or int(row['length']) != lengths[contig]:
                    raise SystemExit(f'ERROR: stale false-duplication name/length: {contig}')
                if contig not in contig_to_mag:
                    raise SystemExit(
                        f'ERROR: false-duplication result is outside final-MAG scope: {contig}'
                    )
                n_pairs = int(row['n_duplicated_pairs'])
                n_fully = int(row['n_fully_evaluable_pairs'])
                n_positive = int(row['n_false_duplication_candidates'])
                n_negative = int(row['n_no_false_duplication_signal_pairs'])
                n_insufficient = int(row['n_insufficient_evidence_pairs'])
                if min(n_pairs, n_fully, n_positive, n_negative, n_insufficient) < 0:
                    raise SystemExit(f'ERROR: negative false-duplication count: {contig}')
                if n_pairs != n_positive + n_negative + n_insufficient:
                    raise SystemExit(f'ERROR: inconsistent false-duplication pair counts: {contig}')
                if n_fully != n_positive + n_negative:
                    raise SystemExit(f'ERROR: inconsistent evaluable-pair count: {contig}')
                call = row['false_duplication_call']
                status = row['structural_assessment_status']
                complete = row['structural_assessment_complete']
                detected = row['false_duplication_candidate_detected']
                assessment_call = row['false_duplication_assessment_call']
                if call not in allowed_calls:
                    raise SystemExit(f'ERROR: invalid false-duplication call: {contig}/{call}')
                if status not in allowed_status:
                    raise SystemExit(f'ERROR: invalid false-duplication status: {contig}/{status}')
                if complete not in {'Y', 'N'} or detected not in {'Y', 'N'}:
                    raise SystemExit(f'ERROR: invalid false-duplication Y/N field: {contig}')
                if assessment_call not in {'Y', 'N', 'U'}:
                    raise SystemExit(
                        f'ERROR: invalid false-duplication assessment call: '
                        f'{contig}/{assessment_call}'
                    )
                expected_complete = status in {
                    'Complete', 'Not_applicable_no_candidate_pairs'
                }
                if complete != ('Y' if expected_complete else 'N'):
                    raise SystemExit(
                        f'ERROR: inconsistent false-duplication completeness: {contig}'
                    )
                if detected != ('Y' if n_positive > 0 else 'N'):
                    raise SystemExit(
                        f'ERROR: inconsistent false-duplication detection field: {contig}'
                    )
                if n_positive > 0:
                    expected_assessment = 'Y'
                elif n_insufficient > 0:
                    expected_assessment = 'U'
                else:
                    expected_assessment = 'N'
                if assessment_call != expected_assessment:
                    raise SystemExit(
                        f'ERROR: inconsistent false-duplication Y/N/U call: {contig}'
                    )
                false_dup[contig] = call
                false_dup_n_pairs[contig] = n_pairs
                false_dup_n_fully[contig] = n_fully
                false_dup_n_positive[contig] = n_positive
                false_dup_n_insufficient[contig] = n_insufficient
                false_dup_detection[contig] = detected
                false_dup_assessment_call[contig] = assessment_call
                false_dup_assessment_status[contig] = status
                false_dup_assessment_complete[contig] = complete
        expected_false_dup_contigs = set(contig_to_mag)
        if set(false_dup) != expected_false_dup_contigs:
            missing = sorted(expected_false_dup_contigs - set(false_dup))[:5]
            extra = sorted(set(false_dup) - expected_false_dup_contigs)[:5]
            raise SystemExit(
                'ERROR: false-duplication table does not match final-MAG contigs; '
                f'missing={missing}, extra={extra}'
            )
        false_dup_default = 'Missing_contig_result'
    else:
        false_dup_default = unavailable_label('false_duplication', run_false_dup)

    false_dup_clipping_contig = {}

    false_dup_clipping_mag = {}

    if module_available['false_duplication_clipping']:
        required = {
            'MAG_ID', 'assembler', 'contig', 'length',
            'CLIPPING_DATA_AVAILABLE_Y_N', 'n_false_duplication_pairs',
            'n_false_duplication_pairs_with_any_validated_clipping',
            'n_false_duplication_defining_boundaries',
            'n_false_duplication_defining_boundaries_with_validated_clipping',
            'ANY_FALSE_DUPLICATION_PAIR_WITH_VALIDATED_CLIPPING_Y_N_NA',
            'FALSE_DUPLICATION_CLIPPING_SUMMARY',
        }
        allowed_summaries = {
            'CLIPPING_DATA_NOT_AVAILABLE',
            'FALSE_DUPLICATION_WITH_VALIDATED_CLIPPING_CO_OCCURRENCE',
            'FALSE_DUPLICATION_WITHOUT_VALIDATED_CLIPPING_CO_OCCURRENCE',
            'NO_FALSE_DUPLICATION_PAIR',
        }
        with open(false_dup_clipping_contig_file) as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            missing = required - set(reader.fieldnames or [])
            if missing:
                raise SystemExit(
                    'ERROR: false-duplication clipping contig table lacks columns '
                    f'{sorted(missing)}'
                )
            for row in reader:
                contig = row['contig']
                if contig in false_dup_clipping_contig:
                    raise SystemExit(
                        f'ERROR: duplicate false-duplication clipping contig: {contig}'
                    )
                if contig not in lengths or int(row['length']) != lengths[contig]:
                    raise SystemExit(
                        f'ERROR: stale false-duplication clipping contig: {contig}'
                    )
                available_call = row['CLIPPING_DATA_AVAILABLE_Y_N']
                any_call = row[
                    'ANY_FALSE_DUPLICATION_PAIR_WITH_VALIDATED_CLIPPING_Y_N_NA'
                ]
                if available_call not in {'Y', 'N'} or any_call not in {'Y', 'N', 'NA'}:
                    raise SystemExit(
                        f'ERROR: invalid false-duplication clipping call: {contig}'
                    )
                counts = [
                    int(row['n_false_duplication_pairs']),
                    int(row['n_false_duplication_pairs_with_any_validated_clipping']),
                    int(row['n_false_duplication_defining_boundaries']),
                    int(row['n_false_duplication_defining_boundaries_with_validated_clipping']),
                ]
                if min(counts) < 0 or counts[1] > counts[0] or counts[3] > counts[2]:
                    raise SystemExit(
                        f'ERROR: inconsistent false-duplication clipping counts: {contig}'
                    )
                expected_any = 'NA' if available_call == 'N' else ('Y' if counts[1] else 'N')
                if any_call != expected_any or row['FALSE_DUPLICATION_CLIPPING_SUMMARY'] not in allowed_summaries:
                    raise SystemExit(
                        f'ERROR: inconsistent false-duplication clipping summary: {contig}'
                    )
                false_dup_clipping_contig[contig] = row
        if set(false_dup_clipping_contig) != set(lengths):
            missing = sorted(set(lengths) - set(false_dup_clipping_contig))[:5]
            extra = sorted(set(false_dup_clipping_contig) - set(lengths))[:5]
            raise SystemExit(
                'ERROR: false-duplication clipping contig scope differs from assembly; '
                f'missing={missing}, extra={extra}'
            )

        mag_required = required - {'assembler', 'contig', 'length'} | {'n_contigs'}
        with open(false_dup_clipping_mag_file) as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            missing = mag_required - set(reader.fieldnames or [])
            if missing:
                raise SystemExit(
                    'ERROR: false-duplication clipping MAG table lacks columns '
                    f'{sorted(missing)}'
                )
            for row in reader:
                mag = row['MAG_ID']
                if not mag or mag in false_dup_clipping_mag:
                    raise SystemExit(
                        f'ERROR: blank or duplicate false-duplication clipping MAG: {mag!r}'
                    )
                if mag not in all_mags or int(row['n_contigs']) < 1:
                    raise SystemExit(
                        f'ERROR: invalid false-duplication clipping MAG: {mag}'
                    )
                available_call = row['CLIPPING_DATA_AVAILABLE_Y_N']
                any_call = row[
                    'ANY_FALSE_DUPLICATION_PAIR_WITH_VALIDATED_CLIPPING_Y_N_NA'
                ]
                counts = [
                    int(row['n_false_duplication_pairs']),
                    int(row['n_false_duplication_pairs_with_any_validated_clipping']),
                    int(row['n_false_duplication_defining_boundaries']),
                    int(row['n_false_duplication_defining_boundaries_with_validated_clipping']),
                ]
                expected_any = 'NA' if available_call == 'N' else ('Y' if counts[1] else 'N')
                if (
                    available_call not in {'Y', 'N'}
                    or any_call != expected_any
                    or min(counts) < 0
                    or counts[1] > counts[0]
                    or counts[3] > counts[2]
                    or row['FALSE_DUPLICATION_CLIPPING_SUMMARY'] not in allowed_summaries
                ):
                    raise SystemExit(
                        f'ERROR: inconsistent false-duplication clipping MAG summary: {mag}'
                    )
                false_dup_clipping_mag[mag] = row
        if set(false_dup_clipping_mag) != all_mags:
            missing = sorted(all_mags - set(false_dup_clipping_mag))[:5]
            extra = sorted(set(false_dup_clipping_mag) - all_mags)[:5]
            raise SystemExit(
                'ERROR: false-duplication clipping MAG scope differs from MAG map; '
                f'missing={missing}, extra={extra}'
            )

    phantom_counts = defaultdict(Counter)

    if module_available['phantom']:
        with open(phantom_file) as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            required = {
                'region_id', 'contig', 'start_0based', 'end_0based', 'length_bp',
                'raw_read_blast_hsps', 'phantom_sequence_call',
            }
            if not required.issubset(reader.fieldnames or []):
                raise SystemExit(f'ERROR: malformed phantom table: {phantom_file}')
            seen_regions = set()
            for row in reader:
                region = row['region_id']
                contig = row['contig']
                call = row['phantom_sequence_call']
                if region in seen_regions:
                    raise SystemExit(f'ERROR: duplicate phantom region: {region}')
                seen_regions.add(region)
                if contig not in lengths:
                    raise SystemExit(f'ERROR: phantom contig absent from assembly: {contig}')
                start, end, region_length = (
                    int(row['start_0based']), int(row['end_0based']),
                    int(row['length_bp']),
                )
                if start < 0 or end > lengths[contig] or end - start != region_length:
                    raise SystemExit(f'ERROR: invalid phantom coordinates: {region}')
                if int(row['raw_read_blast_hsps']) < 0 or call not in {
                    'No_raw_read_BLAST_hit_phantom_candidate',
                    'Raw_read_BLAST_hit_detected',
                }:
                    raise SystemExit(f'ERROR: invalid phantom result: {region}')
                phantom_counts[contig][call] += 1

    def phantom_call(contig):
        if not module_available['phantom']:
            return unavailable_label('phantom', run_phantom)
        counts = phantom_counts[contig]
        if not counts:
            return 'No_zero_coverage_region_ge_threshold'
        phantom_labels = {
            'No_raw_read_BLAST_hit_phantom_candidate',
        }
        ambiguity_labels = {
            'Raw_read_BLAST_hit_detected',
        }
        if any(counts[label] for label in phantom_labels):
            return 'Phantom_sequence_candidate'
        if any(counts[label] for label in ambiguity_labels):
            return 'Zero_region_has_raw_read_BLAST_hit'
        return 'Unrecognized_phantom_module_call'

    manual_review_phantom = {'Zero_region_has_raw_read_BLAST_hit'}

    positive_allele_disagreement = {
        'Unsupported_assembly_variant', 'Minority_assembly_variant',
    }

    positive_phantom = {'Phantom_sequence_candidate'}

    positive_repeat = {'Excessive_repeat_candidate'}

    positive_false_dup = {'False_duplication_candidate'}

    manual_review_false_dup = set()

    incomplete_labels = {
        'Disabled', 'Not_run', 'Missing_contig_result', 'Repeat_analysis_skipped',
        'False_duplication_assessment_incomplete',
    }

    with open(out_contig, 'w', newline='') as out:
        out.write(
            'MAG_ID\tcontig\tlength\tvalidated_clipping_breakpoints\t'
            'Taxonomic_chimera_GTDB\t'
            'GTDB_assessment_status\tGTDB_taxonomic_chimera_ranks\t'
            'GTDB_taxonomic_chimera_with_clipping_support\t'
            'GTDB_clipping_independent_taxonomic_shift\t'
            'GTDB_CHIMERA_EVIDENCE\t'
            'Taxonomic_chimera\t'
            'assembly_allele_support_call\tphantom_sequence_call\texcessive_repeat_call\t'
            'repeat_fraction\tfalse_duplication_call\t'
            'false_duplication_candidate_detected\tfalse_duplication_assessment_call\t'
            'n_false_duplication_candidate_pairs\t'
            'n_false_duplication_fully_evaluable_pairs\t'
            'n_false_duplication_candidates\t'
            'n_false_duplication_insufficient_pairs\t'
            'false_duplication_assessment_status\t'
            'false_duplication_assessment_complete\t'
            'false_duplication_with_validated_clipping_Y_N_NA\t'
            'n_false_duplication_pairs_with_validated_clipping\t'
            'n_false_duplication_defining_boundaries\t'
            'n_false_duplication_defining_boundaries_with_validated_clipping\t'
            'false_duplication_clipping_summary\t'
            'completed_contig_modules\tn_flagged_modules\t'
            'n_manual_review_modules\tintegrated_interpretation\n'
        )
        for contig in sorted(lengths):
            if module_available['clipping']:
                c_clip_count = str(clipping_counts[contig])
                c_clip_status = (
                    'Validated_clipping_breakpoint_detected'
                    if clipping_counts[contig] > 0
                    else 'No_validated_clipping_breakpoint'
                )
            else:
                c_clip_count = 'NA'
                c_clip_status = unavailable_label('clipping', run_clipping)
            gtdb_row = gtdb_chimera.get(contig)
            if gtdb_row is None:
                c_tax_gtdb = gtdb_chimera_default
                c_gtdb_status = gtdb_chimera_default
                c_gtdb_ranks = 'NA'
                c_gtdb_clip_support = 'NA'
                c_gtdb_clip_independent = 'NA'
                c_gtdb_evidence = gtdb_chimera_default
            else:
                c_tax_gtdb = gtdb_row['Taxonomic_chimera']
                c_gtdb_status = gtdb_row['GTDB_assessment_status']
                c_gtdb_ranks = gtdb_row['taxonomic_shift_ranks']
                c_gtdb_clip_support = gtdb_row['Taxonomic_chimera_with_clipping_support']
                c_gtdb_clip_independent = gtdb_row['Clipping_independent_taxonomic_shift']
                c_gtdb_evidence = gtdb_row['CHIMERA_EVIDENCE']
            c_tax_chimera = c_tax_gtdb
            c_chimera_status = c_tax_chimera
            c_allele = allele.get(contig, allele_default)
            c_pha = phantom_call(contig)
            c_rep = repeat.get(contig, repeat_default)

            if run_false_dup != 'Y':
                c_fdup = 'Disabled'
                c_fdup_detected = 'NA'
                c_fdup_assessment_call = 'NA'
                c_fdup_pairs = 'NA'
                c_fdup_fully = 'NA'
                c_fdup_positive = 'NA'
                c_fdup_insufficient = 'NA'
                c_fdup_assessment = 'Disabled'
                c_fdup_complete = 'N'
                c_fdup_completion_call = 'Disabled'
                c_fdup_clip_call = 'NA'
                c_fdup_clip_pairs = 'NA'
                c_fdup_clip_boundaries = 'NA'
                c_fdup_clip_boundaries_validated = 'NA'
                c_fdup_clip_summary = 'Disabled'
            elif contig in contig_to_mag:
                c_fdup = false_dup.get(contig, false_dup_default)
                c_fdup_detected = false_dup_detection.get(contig, 'NA')
                c_fdup_assessment_call = false_dup_assessment_call.get(contig, 'NA')
                c_fdup_pairs = false_dup_n_pairs.get(contig, 'NA')
                c_fdup_fully = false_dup_n_fully.get(contig, 'NA')
                c_fdup_positive = false_dup_n_positive.get(contig, 'NA')
                c_fdup_insufficient = false_dup_n_insufficient.get(contig, 'NA')
                c_fdup_assessment = false_dup_assessment_status.get(
                    contig, 'Missing_contig_result'
                )
                c_fdup_complete = false_dup_assessment_complete.get(contig, 'N')
                c_fdup_completion_call = (
                    c_fdup if c_fdup_complete == 'Y'
                    else 'False_duplication_assessment_incomplete'
                )
                clipping_row = false_dup_clipping_contig[contig]
                c_fdup_clip_call = clipping_row[
                    'ANY_FALSE_DUPLICATION_PAIR_WITH_VALIDATED_CLIPPING_Y_N_NA'
                ]
                c_fdup_clip_pairs = clipping_row[
                    'n_false_duplication_pairs_with_any_validated_clipping'
                ]
                c_fdup_clip_boundaries = clipping_row[
                    'n_false_duplication_defining_boundaries'
                ]
                c_fdup_clip_boundaries_validated = clipping_row[
                    'n_false_duplication_defining_boundaries_with_validated_clipping'
                ]
                c_fdup_clip_summary = clipping_row['FALSE_DUPLICATION_CLIPPING_SUMMARY']
            else:
                c_fdup = 'Not_applicable_non_MAG_contig'
                c_fdup_detected = 'N'
                c_fdup_assessment_call = 'N'
                c_fdup_pairs = 0
                c_fdup_fully = 0
                c_fdup_positive = 0
                c_fdup_insufficient = 0
                c_fdup_assessment = 'Not_applicable_non_MAG_contig'
                c_fdup_complete = 'Y'
                c_fdup_completion_call = c_fdup
                c_fdup_clip_call = 'NA'
                c_fdup_clip_pairs = 0
                c_fdup_clip_boundaries = 0
                c_fdup_clip_boundaries_validated = 0
                c_fdup_clip_summary = 'Not_applicable_non_MAG_contig'

            calls = [
                c_clip_status, c_chimera_status, c_allele, c_pha, c_rep,
                c_fdup_completion_call,
            ]
            completed = sum(call not in incomplete_labels for call in calls)
            flags = sum([
                clipping_counts[contig] > 0 if module_available['clipping'] else False,
                c_chimera_status == 'Y',
                c_allele in positive_allele_disagreement,
                c_pha in positive_phantom,
                c_rep in positive_repeat,
                c_fdup in positive_false_dup,
            ])
            manual_reviews = sum([
                c_pha in manual_review_phantom,
                c_fdup in manual_review_false_dup,
            ])
            if flags:
                interpretation = 'QC_signal_detected'
            elif manual_reviews:
                interpretation = 'Manual_review_required'
            elif completed < expected_contig_modules:
                interpretation = 'Incomplete_QC_modules'
            else:
                interpretation = 'No_major_signal_in_completed_modules'
            out.write(
                f'{contig_to_mag.get(contig,"NA")}\t{contig}\t{lengths[contig]}\t'
                f'{c_clip_count}\t'
                f'{c_tax_gtdb}\t{c_gtdb_status}\t'
                f'{c_gtdb_ranks}\t{c_gtdb_clip_support}\t{c_gtdb_clip_independent}\t'
                f'{c_gtdb_evidence}\t'
                f'{c_tax_chimera}\t'
                f'{c_allele}\t{c_pha}\t{c_rep}\t'
                f'{repeat_frac.get(contig,"NA")}\t{c_fdup}\t'
                f'{c_fdup_detected}\t{c_fdup_assessment_call}\t'
                f'{c_fdup_pairs}\t{c_fdup_fully}\t{c_fdup_positive}\t'
                f'{c_fdup_insufficient}\t{c_fdup_assessment}\t{c_fdup_complete}\t'
                f'{c_fdup_clip_call}\t{c_fdup_clip_pairs}\t'
                f'{c_fdup_clip_boundaries}\t{c_fdup_clip_boundaries_validated}\t'
                f'{c_fdup_clip_summary}\t'
                f'{completed}\t{flags}\t'
                f'{manual_reviews}\t{interpretation}\n'
            )

    circularity = {}

    circularity_contig_counts = {}

    circularity_eligibility = {}

    if module_available['circularity']:
        with open(circularity_file) as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            fields = set(reader.fieldnames or [])
            required = {
                'ID', 'Assembler_reported_circularity',
                'Circularity_evaluable', 'Circularity_eligibility',
                'Read_supported_circular_topology', 'High_confidence_circular',
                'Additional_circular_candidate', 'Circular',
            }
            missing = required - fields
            if missing:
                raise SystemExit(
                    f'ERROR: circularity-only table lacks columns '
                    f'{sorted(missing)}: {circularity_file}'
                )
            for row in reader:
                mag = row['ID']
                if mag in circularity:
                    raise SystemExit(f'ERROR: duplicate circularity ID: {mag}')
                assembler = row['Assembler_reported_circularity']
                evaluable = row['Circularity_evaluable']
                topology = row['Read_supported_circular_topology']
                high_conf = row['High_confidence_circular']
                additional_candidate = row['Additional_circular_candidate']
                circ = row['Circular']
                eligibility = row['Circularity_eligibility']
                if assembler not in {'Y', 'N', 'NA'}:
                    raise SystemExit(
                        f'ERROR: invalid assembler circularity for {mag}: {assembler}'
                    )
                if evaluable not in {'Y', 'N'}:
                    raise SystemExit(f'ERROR: invalid Circularity_evaluable for {mag}: {evaluable}')
                if topology not in {'Y', 'N', 'NA'}:
                    raise SystemExit(
                        f'ERROR: invalid Read_supported_circular_topology for {mag}: {topology}'
                    )
                derived_values = {'Y', 'N', 'NA'} if assembler == 'NA' else {'Y', 'N'}
                if (high_conf not in derived_values or additional_candidate not in derived_values
                        or circ not in {'Y', 'N', 'NA'}):
                    raise SystemExit(f'ERROR: invalid derived circularity field for {mag}')
                if eligibility not in {
                    'Y_single_contig', 'N_not_single_contig',
                    'N_contig_shorter_than_two_anchors',
                }:
                    raise SystemExit(
                        f'ERROR: invalid circularity eligibility for {mag}: {eligibility}'
                    )
                if eligibility != 'Y_single_contig':
                    if evaluable != 'N' or topology != 'NA':
                        raise SystemExit(
                            f'ERROR: non-single-contig circularity row must be evaluable=N and topology=NA for {mag}'
                        )
                elif evaluable == 'Y':
                    if topology not in {'Y', 'N'}:
                        raise SystemExit(
                            f'ERROR: evaluable single-contig MAG must have Y/N topology for {mag}'
                        )
                elif topology != 'NA':
                    raise SystemExit(
                        f'ERROR: non-evaluable single-contig MAG must have topology=NA for {mag}'
                    )
                if assembler == 'NA':
                    expected_high_conf = expected_additional = 'NA'
                    expected_circ = topology
                else:
                    expected_high_conf = 'Y' if assembler == 'Y' and topology == 'Y' else 'N'
                    expected_additional = 'Y' if assembler == 'N' and topology == 'Y' else 'N'
                    expected_circ = expected_high_conf
                if high_conf != expected_high_conf:
                    raise SystemExit(
                        f'ERROR: inconsistent high-confidence circularity for {mag}'
                    )
                if additional_candidate != expected_additional:
                    raise SystemExit(
                        f'ERROR: inconsistent additional circular candidate call for {mag}'
                    )
                if circ != expected_circ:
                    raise SystemExit(f'ERROR: inconsistent Circular call for {mag}')
                circularity[mag] = {'Circular': circ, 'Assembler': assembler}
                circularity_eligibility[mag] = eligibility
        all_mags.update(circularity)

        with open(circularity_audit_file) as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            required = {'ID', 'Contig_count', 'Eligibility'}
            missing = required - set(reader.fieldnames or [])
            if missing:
                raise SystemExit(
                    f'ERROR: circularity audit lacks columns {sorted(missing)}: '
                    f'{circularity_audit_file}'
                )
            for row in reader:
                mag = row['ID']
                if mag in circularity_contig_counts:
                    raise SystemExit(f'ERROR: duplicate circularity audit ID: {mag}')
                circularity_contig_counts[mag] = int(row['Contig_count'])
                if row['Eligibility'] != circularity_eligibility.get(mag):
                    raise SystemExit(
                        f'ERROR: circularity final/audit eligibility differs for {mag}'
                    )
        if set(circularity_contig_counts) != set(circularity):
            raise SystemExit('ERROR: circularity final/audit MAG IDs differ')

    false_closure = {}

    if module_available['false_closure']:
        with open(false_closure_file) as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            fields = set(reader.fieldnames or [])
            required = {
                'ID', 'MAG_CONTIG_COUNT', 'FALSE_CLOSURE_EVIDENCE',
                'Paper_small_genome_proxy',
                'Premature_circularization_candidate',
                'Premature_circularization_screen',
            }
            missing = required - fields
            if missing:
                raise SystemExit(
                    f'ERROR: false-closure table lacks columns '
                    f'{sorted(missing)}: {false_closure_file}'
                )
            for row in reader:
                mag = row['ID']
                if mag in false_closure:
                    raise SystemExit(f'ERROR: duplicate false-closure ID: {mag}')
                assembler_circular = row.get(
                    'ASSEMBLER_REPORTED_CIRCULARITY',
                    circularity.get(mag, {}).get('Assembler', 'NA')
                )
                if assembler_circular not in {'Y', 'N', 'NA'}:
                    raise SystemExit(
                        f'ERROR: assembler circularity unavailable for '
                        f'false-closure validation of {mag}'
                    )
                false_evidence = row['FALSE_CLOSURE_EVIDENCE']
                proxy = row['Paper_small_genome_proxy']
                premature = row['Premature_circularization_candidate']
                screen = row['Premature_circularization_screen']
                values = {
                    'FALSE_CLOSURE_EVIDENCE': false_evidence,
                    'Paper_small_genome_proxy': proxy,
                    'Premature_circularization_candidate': premature,
                    'Premature_circularization_screen': screen,
                }
                invalid = {k: v for k, v in values.items() if v not in {'Y', 'N'}}
                if invalid:
                    raise SystemExit(
                        f'ERROR: non-binary false-closure call for {mag}: {invalid}'
                    )
                try:
                    false_closure_contig_count = int(row['MAG_CONTIG_COUNT'])
                except ValueError as exc:
                    raise SystemExit(
                        f'ERROR: invalid MAG_CONTIG_COUNT for false-closure row {mag}: {exc}'
                    )
                expected_premature = (
                    'Y' if false_closure_contig_count == 1 and assembler_circular == 'Y'
                    and false_evidence == 'Y' else 'N'
                )
                expected_screen = (
                    'Y' if expected_premature == 'Y' or proxy == 'Y' else 'N'
                )
                if premature != expected_premature or screen != expected_screen:
                    raise SystemExit(
                        f'ERROR: inconsistent false-closure calls for {mag}: '
                        f'expected premature={expected_premature}, screen={expected_screen}; '
                        f'observed premature={premature}, screen={screen}'
                    )
                false_closure[mag] = {
                    'Premature': premature,
                    'Premature_screen': screen,
                    'FALSE_CLOSURE_EVIDENCE': false_evidence,
                    'Candidate_additional_contigs': row.get(
                        'Candidate_additional_target_contigs', 'NA'
                    ),
                    'Paper_small_genome_proxy': proxy,
                }
        all_mags.update(false_closure)

        with open(false_closure_audit_file) as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            fields = set(reader.fieldnames or [])
            if 'ID' not in fields:
                raise SystemExit(
                    f'ERROR: false-closure audit lacks ID: {false_closure_audit_file}'
                )
            audit_ids = [row['ID'] for row in reader]
            if len(audit_ids) != len(set(audit_ids)):
                raise SystemExit('ERROR: duplicate false-closure audit ID')
            if set(audit_ids) != set(false_closure):
                raise SystemExit('ERROR: false-closure final/audit MAG IDs differ')

    mag_stats = defaultdict(lambda: defaultdict(int))

    gtdb_mag_contig_rows = defaultdict(list)

    with open(out_contig) as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for row in reader:
            mag = row['MAG_ID']
            if mag == 'NA':
                continue
            all_mags.add(mag)
            stats = mag_stats[mag]
            gtdb_mag_contig_rows[mag].append(row)
            stats['n_contigs'] += 1
            if row['validated_clipping_breakpoints'] != 'NA':
                stats['clip_sites'] += int(row['validated_clipping_breakpoints'])
            stats['taxonomic_chimera_gtdb'] += (
                row['Taxonomic_chimera_GTDB'] == 'Y'
            )
            stats['taxonomic_chimera'] += row['Taxonomic_chimera'] == 'Y'
            stats['allele_disagreement'] += (
                row['assembly_allele_support_call'] in positive_allele_disagreement
            )
            stats['phantom'] += row['phantom_sequence_call'] in positive_phantom
            stats['repeat'] += row['excessive_repeat_call'] in positive_repeat
            stats['false_dup'] += row['false_duplication_call'] in positive_false_dup
            if row['n_false_duplication_candidate_pairs'] != 'NA':
                stats['false_dup_pairs'] += int(
                    row['n_false_duplication_candidate_pairs']
                )
                stats['false_dup_fully'] += int(
                    row['n_false_duplication_fully_evaluable_pairs']
                )
                stats['false_dup_insufficient_pairs'] += int(
                    row['n_false_duplication_insufficient_pairs']
                )
            stats['manual_review'] += int(row['n_manual_review_modules'])
            stats['incomplete_contigs'] += (
                int(row['completed_contig_modules']) < expected_contig_modules
            )

    gtdb_mag_evidence = {}

    for mag, rows in gtdb_mag_contig_rows.items():
        labels = {row['GTDB_CHIMERA_EVIDENCE'] for row in rows}
        for label in (
            'HIGH_CONFIDENCE_HIFI_STRUCTURAL_CHIMERA',
            'HIFI_STRUCTURAL_CHIMERA',
            'UNRESOLVED_HIFI_STRUCTURAL_EVIDENCE',
            'HIFI_CONTINUITY_SUPPORTED',
        ):
            if label in labels:
                gtdb_mag_evidence[mag] = label
                break
        else:
            if labels == {'NO_TAXONOMIC_SHIFT'}:
                gtdb_mag_evidence[mag] = 'NO_TAXONOMIC_SHIFT'
            else:
                gtdb_mag_evidence[mag] = 'NOT_ASSESSABLE'

    with open(out_mag, 'w', newline='') as out:
        out.write(
            'MAG_ID\tn_contigs\t'
            'Clipping_Y_N\tTaxonomic_chimera_Y_N\t'
            'Assembly_allele_disagreement_Y_N\tPhantom_sequence_Y_N\t'
            'Excessive_repeat_Y_N\tFalse_duplication_Y_N\t'
            'False_closure_Y_N\tPremature_circularization_Y_N\t'
            'Any_detected_QC_signal_Y_N\t'
            'Circularity_eligibility\tCircular\t'
            'False_closure_evidence\t'
            'False_closure_candidate_additional_contigs\t'
            'Paper_small_genome_proxy\t'
            'Premature_circularization_candidate\tPremature_circularization\t'
            'Premature_circularization_screen\t'
            'n_validated_clipping_breakpoints\t'
            'n_taxonomic_chimera_GTDB_flags\t'
            'n_taxonomic_chimera_flags\t'
            'GTDB_CHIMERA_EVIDENCE\t'
            'n_assembly_allele_disagreement_flags\t'
            'n_phantom_sequence_flags\texcessive_repeat_call\trepeat_fraction\t'
            'n_excessive_repeat_contig_flags\t'
            'n_false_duplication_flags\t'
            'false_duplication_candidate_detected\t'
            'false_duplication_assessment_call\t'
            'n_false_duplication_candidate_pairs\t'
            'n_false_duplication_fully_evaluable_pairs\t'
            'n_false_duplication_insufficient_pairs\t'
            'false_duplication_assessment_status\t'
            'false_duplication_assessment_complete\t'
            'false_duplication_with_validated_clipping_Y_N_NA\t'
            'n_false_duplication_pairs_with_validated_clipping\t'
            'n_false_duplication_defining_boundaries\t'
            'n_false_duplication_defining_boundaries_with_validated_clipping\t'
            'false_duplication_clipping_summary\t'
            'n_manual_review_modules\tn_contigs_with_incomplete_modules\t'
            'MAG_integrated_QC_call\n'
        )
        for mag in sorted(all_mags):
            stats = mag_stats[mag]
            if module_available['circularity'] and mag in circularity_contig_counts:
                audit_n = circularity_contig_counts[mag]
                if stats['n_contigs'] and stats['n_contigs'] != audit_n:
                    raise SystemExit(
                        f'ERROR: MAG contig-map/audit count mismatch for {mag}: '
                        f'{stats["n_contigs"]} versus {audit_n}'
                    )
                if not stats['n_contigs']:
                    stats['n_contigs'] = audit_n
            if run_circularity != 'Y':
                eligibility, circ = 'Disabled', 'Disabled'
            elif not module_available['circularity']:
                eligibility, circ = 'Not_run', 'Not_run'
            elif mag in circularity:
                eligibility = circularity_eligibility[mag]
                circ = circularity[mag]['Circular']
            else:
                eligibility, circ = 'Not_target', 'Not_target'

            if run_false_closure != 'Y':
                premature = false_evidence = false_targets = proxy = screen = 'Disabled'
            elif not module_available['false_closure']:
                premature = false_evidence = false_targets = proxy = screen = 'Not_run'
            elif mag in false_closure:
                fdata = false_closure[mag]
                premature = fdata['Premature']
                false_evidence = fdata['FALSE_CLOSURE_EVIDENCE']
                false_targets = fdata['Candidate_additional_contigs']
                proxy = fdata['Paper_small_genome_proxy']
                screen = fdata['Premature_screen']
            else:
                premature = false_evidence = false_targets = proxy = screen = 'Not_target'

            if run_false_dup != 'Y':
                fd_detected = 'NA'
                fd_assessment_call = 'NA'
                fd_status = 'Disabled'
                fd_complete = 'N'
            elif stats['false_dup_pairs'] == 0:
                fd_detected = 'N'
                fd_assessment_call = 'N'
                fd_status = 'Not_applicable_no_candidate_pairs'
                fd_complete = 'Y'
            else:
                fd_detected = 'Y' if stats['false_dup'] > 0 else 'N'
                if stats['false_dup_insufficient_pairs'] == 0:
                    fd_status = 'Complete'
                    fd_assessment_call = 'Y' if stats['false_dup'] > 0 else 'N'
                    fd_complete = 'Y'
                else:
                    fd_status = (
                        'Partial' if stats['false_dup_fully'] > 0
                        else 'Insufficient_for_all_candidate_pairs'
                    )
                    fd_assessment_call = 'Y' if stats['false_dup'] > 0 else 'U'
                    fd_complete = 'N'

            if run_false_dup != 'Y':
                fd_clip_call = 'NA'
                fd_clip_pairs = 'NA'
                fd_clip_boundaries = 'NA'
                fd_clip_boundaries_validated = 'NA'
                fd_clip_summary = 'Disabled'
            else:
                clipping_row = false_dup_clipping_mag[mag]
                fd_clip_call = clipping_row[
                    'ANY_FALSE_DUPLICATION_PAIR_WITH_VALIDATED_CLIPPING_Y_N_NA'
                ]
                fd_clip_pairs = clipping_row[
                    'n_false_duplication_pairs_with_any_validated_clipping'
                ]
                fd_clip_boundaries = clipping_row[
                    'n_false_duplication_defining_boundaries'
                ]
                fd_clip_boundaries_validated = clipping_row[
                    'n_false_duplication_defining_boundaries_with_validated_clipping'
                ]
                fd_clip_summary = clipping_row['FALSE_DUPLICATION_CLIPPING_SUMMARY']

            if run_repeat != 'Y':
                repeat_call, mag_repeat_fraction = 'Disabled', 'NA'
            elif not module_available['repeat']:
                repeat_call, mag_repeat_fraction = 'Not_run', 'NA'
            elif mag in mag_repeat:
                repeat_call = mag_repeat[mag]['call']
                mag_repeat_fraction = mag_repeat[mag]['fraction']
                if stats['n_contigs'] and stats['n_contigs'] != mag_repeat[mag]['n_contigs']:
                    raise SystemExit(
                        f'ERROR: MAG repeat/contig count mismatch for {mag}: '
                        f'{mag_repeat[mag]["n_contigs"]} versus {stats["n_contigs"]}'
                    )
            else:
                repeat_call, mag_repeat_fraction = 'Not_target', 'NA'

            if run_taxonomic_chimera != 'Y':
                gtdb_evidence = 'Disabled'
            elif not module_available['gtdb_chimera']:
                gtdb_evidence = 'Not_run'
            else:
                gtdb_evidence = gtdb_mag_evidence.get(mag, 'Not_target')

            clip_yn = 'Y' if stats['clip_sites'] > 0 else 'N'
            chimera_yn = 'Y' if stats['taxonomic_chimera'] > 0 else 'N'
            allele_yn = 'Y' if stats['allele_disagreement'] > 0 else 'N'
            phantom_yn = 'Y' if stats['phantom'] > 0 else 'N'
            repeat_yn = 'Y' if repeat_call == 'Excessive_repeat_candidate' else 'N'
            false_dup_yn = 'Y' if stats['false_dup'] > 0 else 'N'
            false_closure_yn = 'Y' if false_evidence == 'Y' else 'N'
            premature_yn = 'Y' if premature == 'Y' else 'N'

            any_flag = any(
                x == 'Y' for x in [
                    clip_yn, chimera_yn, allele_yn, phantom_yn, repeat_yn,
                    false_dup_yn, false_closure_yn, premature_yn,
                ]
            )
            any_flag_yn = 'Y' if any_flag else 'N'
            incomplete = (
                stats['incomplete_contigs'] > 0
                or circ in {'Not_run', 'Disabled'}
                or (
                    run_false_closure == 'Y'
                    and not module_available['false_closure']
                )
            )
            circularity_not_evaluable = (
                run_circularity == 'Y'
                and module_available['circularity']
                and eligibility in {
                    'N_not_single_contig',
                    'N_contig_shorter_than_two_anchors',
                }
            )
            if any_flag:
                call = 'QC_signal_detected'
            elif stats['manual_review'] > 0:
                call = 'Manual_review_required'
            elif circularity_not_evaluable:
                call = 'Circularity_not_evaluable'
            elif incomplete:
                call = 'Incomplete_QC_modules'
            else:
                call = 'No_major_signal_in_completed_modules'
            out.write(
                f'{mag}\t{stats["n_contigs"]}\t'
                f'{clip_yn}\t{chimera_yn}\t{allele_yn}\t{phantom_yn}\t'
                f'{repeat_yn}\t{false_dup_yn}\t{false_closure_yn}\t'
                f'{premature_yn}\t{any_flag_yn}\t'
                f'{eligibility}\t{circ}\t'
                f'{false_evidence}\t{false_targets}\t{proxy}\t'
                f'{premature}\t{premature}\t{screen}\t'
                f'{stats["clip_sites"]}\t'
                f'{stats["taxonomic_chimera_gtdb"]}\t'
                f'{stats["taxonomic_chimera"]}\t'
                f'{gtdb_evidence}\t'
                f'{stats["allele_disagreement"]}\t'
                f'{stats["phantom"]}\t{repeat_call}\t{mag_repeat_fraction}\t'
                f'{stats["repeat"]}\t{stats["false_dup"]}\t'
                f'{fd_detected}\t{fd_assessment_call}\t'
                f'{stats["false_dup_pairs"]}\t{stats["false_dup_fully"]}\t'
                f'{stats["false_dup_insufficient_pairs"]}\t{fd_status}\t{fd_complete}\t'
                f'{fd_clip_call}\t{fd_clip_pairs}\t{fd_clip_boundaries}\t'
                f'{fd_clip_boundaries_validated}\t{fd_clip_summary}\t'
                f'{stats["manual_review"]}\t'
                f'{stats["incomplete_contigs"]}\t{call}\n'
            )

    print('Module availability: ' + ', '.join(f'{k}={v}' for k, v in module_available.items()))

    print(f'Contigs summarized: {len(lengths)}')

    print(f'MAGs summarized: {len(all_mags)}')


LEAD = ('MAG_ID', 'n_contigs', 'MAG_integrated_QC_call')

DROP = ('Premature_circularization_candidate',)


def arrange_mag(path):
    """Put the overall call next to the MAG ID and drop duplicated columns, in place."""
    with open(path, newline='', encoding='utf-8') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = list(reader.fieldnames or [])
        records = list(reader)
    missing = [field for field in LEAD if field not in fields]
    if missing:
        raise SystemExit(f'ERROR: MAG summary lacks columns: {missing}')
    for field in DROP:
        if field in fields and 'Premature_circularization' in fields:
            if any(row[field] != row['Premature_circularization'] for row in records):
                raise SystemExit(f'ERROR: {field} differs from Premature_circularization')
    order = [*LEAD, *(f for f in fields if f not in LEAD and f not in DROP)]
    temp = Path(f'{path}.tmp')
    with open(temp, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=order, delimiter='\t', lineterminator='\n',
                                extrasaction='ignore')
        writer.writeheader()
        writer.writerows(records)
    temp.replace(path)


def rows(path, required):
    with open(path, newline='', encoding='utf-8') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        header = reader.fieldnames or []
        if len(header) != len(set(header)) or not set(required).issubset(header):
            raise ValueError(f"Invalid TSV header: {path}")
        result = list(reader)
    if not result:
        raise ValueError(f"Empty QC table: {path}")
    return result


def validate(contig_path, mag_path):
    contigs = rows(contig_path, ('MAG_ID', 'contig', 'length'))
    mags = rows(mag_path, ('MAG_ID', 'n_contigs', 'MAG_integrated_QC_call', 'Any_detected_QC_signal_Y_N'))
    ids = [row['MAG_ID'] for row in mags]
    if any(not mag or mag == 'NA' for mag in ids) or len(set(ids)) != len(ids):
        raise ValueError('Invalid or duplicated MAG_ID in MAG summary')
    counts = {}
    keys = set()
    for row in contigs:
        key = row['contig']
        if not key or key in keys:
            raise ValueError(f'Duplicated or blank contig in contig summary: {key}')
        keys.add(key)
        if int(row['length']) < 1:
            raise ValueError(f'Invalid contig length: {key}')
        mag = row['MAG_ID']
        if mag != 'NA':
            counts[mag] = counts.get(mag, 0) + 1
    if set(counts) != set(ids):
        raise ValueError('MAG membership differs between MAG and contig summaries')
    for row in mags:
        if int(row['n_contigs']) != counts[row['MAG_ID']]:
            raise ValueError(f'Contig count mismatch for {row["MAG_ID"]}')
    print(f'Validated QC summary: {len(mags)} MAGs, {len(contigs)} contigs '
          f'({sum(counts.values())} in MAGs)')


def _arrange_main():
    if len(sys.argv) != 2:
        raise SystemExit('Usage: summary_tables.py arrange MAG_TSV')
    arrange_mag(sys.argv[1])


def cmd_arrange(argv):
    """Put the overall call next to the MAG ID and drop duplicated columns"""
    sys.argv = [sys.argv[0] + " arrange", *argv]
    return _arrange_main()


def _validate_main():
    if len(sys.argv) != 3:
        raise SystemExit('Usage: summary_tables.py validate CONTIG_TSV MAG_TSV')
    try:
        validate(*(Path(arg) for arg in sys.argv[1:]))
    except (ValueError, OSError, KeyError) as exc:
        raise SystemExit(f'ERROR: {exc}') from None


def cmd_validate(argv):
    """Check that the MAG and contig summaries agree"""
    sys.argv = [sys.argv[0] + " validate", *argv]
    return _validate_main()


COMMANDS = {
    'build': cmd_build,
    'arrange': cmd_arrange,
    'validate': cmd_validate,
}

if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] not in COMMANDS:
        sys.exit(__doc__)
    result = COMMANDS[sys.argv[1]](sys.argv[2:])
    sys.exit(result if isinstance(result, int) else 0)
