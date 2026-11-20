#!/usr/bin/env python3
"""Create MAG-level false-closure final and audit tables."""
from __future__ import annotations

import argparse
import csv
import os
import tempfile
from typing import Dict, Iterable, List, Mapping, NoReturn, Sequence

VERSION = "summarize_false_closure"


def fail(message: str) -> NoReturn:
    raise SystemExit(f"ERROR: {message}")


def read_table(path: str) -> tuple[list[str], Dict[str, Dict[str, str]]]:
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        if "ID" not in fields:
            fail(f"table lacks ID column: {path}")
        rows: Dict[str, Dict[str, str]] = {}
        for row in reader:
            mag = (row.get("ID") or "").strip()
            if not mag or mag in rows:
                fail(f"blank/duplicate ID in {path}: {mag!r}")
            rows[mag] = row
    return fields, rows


def atomic_dict_table(path: str, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=f".{os.path.basename(path)}.", suffix=".tmp",
                               dir=os.path.dirname(os.path.abspath(path)), text=True)
    os.close(fd)
    try:
        with open(tmp, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in rows:
                writer.writerow({f: row.get(f, "NA") for f in fields})
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def call(value: str, allowed: set[str], label: str, mag: str) -> str:
    v = (value or "NA").strip().upper()
    if v not in allowed:
        fail(f"invalid {label} for {mag}: {v!r}")
    return v


def intval(row: Mapping[str, str], field: str, mag: str, minimum: int = 0) -> int:
    try:
        value = int(row.get(field, "0") or 0)
    except ValueError as exc:
        fail(f"invalid integer {field} for {mag}: {exc}")
    if value < minimum:
        fail(f"{field} below {minimum} for {mag}: {value}")
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--false-closure-tsv", required=True)
    parser.add_argument("--false-final", required=True)
    parser.add_argument("--false-audit", required=True)
    parser.add_argument("--metadata", required=True)
    args = parser.parse_args()

    source_fields, source = read_table(args.false_closure_tsv)
    required = {
        "ID", "Contig_count", "Contig", "Genome_length_bp",
        "Assembler_reported_circularity", "MAG_level_input_assembler_circularity",
        "SOURCE_CIRCULAR_CONTIG_COUNT", "SOURCE_CIRCULAR_CONTIGS",
        "SOURCE_CONTIG_CIRCULARITY_RESOLVED_COUNT",
        "SOURCE_CONTIG_CIRCULARITY_UNRESOLVED_COUNT",
        "SOURCE_CONTIG_CIRCULARITY_COMPLETE_Y_N",
        "False_closure_analyzable",
        "Paper_ribosomal_protein_family_count", "Paper_ribosomal_protein_families",
        "Paper_small_genome_length_criterion_Y_N",
        "Paper_small_genome_ribosomal_criterion_Y_N", "Paper_small_genome_proxy",
        "Paper_small_genome_proxy_contigs", "Premature_circularization_candidate",
        "Premature_circularization_screen",
        "Large_insertion_sensitivity_analyzable", "Large_insertion_sensitivity_support",
        "Large_insertion_sensitivity_site_count", "Large_insertion_sensitivity_sites",
        "Large_insertion_sensitivity_min_bp", "Large_insertion_sensitivity_min_QNAMEs",
        "All_MAG_validated_clipping_breakpoint_count",
        "Validated_clipping_breakpoint_count",
        "Validated_internal_clipping_breakpoint_count",
        "Validated_terminal_clipping_breakpoint_count",
        "False_closure_breakpoint_count", "Internal_false_closure_breakpoint_count",
        "Terminal_false_closure_breakpoint_count", "False_closure_breakpoints",
        "False_closure_regions", "Candidate_additional_target_contigs",
        "Candidate_additional_target_MAGs", "Internal_false_closure_evidence",
        "Terminal_false_closure_evidence", "False_closure_evidence", "False_closure",
        "Min_connected_QNAMEs_per_side", "False_closure_reason",
    }
    missing = required - set(source_fields)
    if missing:
        fail(f"false-closure table lacks columns {sorted(missing)}")

    final_fields = [
        "ID", "MAG_CONTIG_COUNT", "FALSE_CLOSURE_ANALYZABLE", "ASSEMBLER_REPORTED_CIRCULARITY",
        "MAG_LEVEL_INPUT_ASSEMBLER_CIRCULARITY",
        "SOURCE_CIRCULAR_CONTIG_COUNT", "SOURCE_CIRCULAR_CONTIGS",
        "SOURCE_CONTIG_CIRCULARITY_COMPLETE_Y_N",
        "Paper_ribosomal_protein_family_count", "Paper_ribosomal_protein_families",
        "Paper_small_genome_length_criterion_Y_N", "Paper_small_genome_ribosomal_criterion_Y_N",
        "Paper_small_genome_proxy", "Paper_small_genome_proxy_contigs",
        "Premature_circularization_candidate", "Premature_circularization_screen",
        "LARGE_INSERTION_SENSITIVITY_ANALYZABLE", "LARGE_INSERTION_SENSITIVITY_SUPPORT",
        "LARGE_INSERTION_SENSITIVITY_SITE_COUNT", "LARGE_INSERTION_SENSITIVITY_SITES",
        "Candidate_additional_target_contigs", "Candidate_additional_target_MAGs",
        "VALIDATED_CLIPPING_PRESENT", "INTERNAL_VALIDATED_CLIPPING_PRESENT",
        "TERMINAL_VALIDATED_CLIPPING_PRESENT", "INTERNAL_FALSE_CLOSURE_EVIDENCE",
        "TERMINAL_FALSE_CLOSURE_EVIDENCE", "FALSE_CLOSURE_EVIDENCE", "FALSE_CLOSURE",
        "FALSE_CLOSURE_BREAKPOINTS", "FALSE_CLOSURE_REGIONS",
    ]
    audit_fields = [
        "ID", "Contig_count", "Contig", "Genome_length_bp",
        "ASSEMBLER_REPORTED_CIRCULARITY", "MAG_LEVEL_INPUT_ASSEMBLER_CIRCULARITY",
        "FALSE_CLOSURE_ANALYZABLE",
        "SOURCE_CIRCULAR_CONTIG_COUNT", "SOURCE_CIRCULAR_CONTIGS",
        "SOURCE_CONTIG_CIRCULARITY_RESOLVED_COUNT", "SOURCE_CONTIG_CIRCULARITY_UNRESOLVED_COUNT",
        "SOURCE_CONTIG_CIRCULARITY_COMPLETE_Y_N",
        "Paper_ribosomal_protein_family_count", "Paper_ribosomal_protein_families",
        "Paper_small_genome_length_criterion_Y_N", "Paper_small_genome_ribosomal_criterion_Y_N",
        "Paper_small_genome_proxy", "Paper_small_genome_proxy_contigs",
        "Premature_circularization_candidate", "Premature_circularization_screen",
        "LARGE_INSERTION_SENSITIVITY_ANALYZABLE", "LARGE_INSERTION_SENSITIVITY_SUPPORT",
        "LARGE_INSERTION_SENSITIVITY_SITE_COUNT", "LARGE_INSERTION_SENSITIVITY_SITES",
        "Large_insertion_sensitivity_min_bp", "Large_insertion_sensitivity_min_QNAMEs",
        "All_MAG_validated_clipping_breakpoint_count", "Validated_clipping_breakpoint_count",
        "Validated_internal_clipping_breakpoint_count", "Validated_terminal_clipping_breakpoint_count",
        "False_closure_breakpoint_count", "Internal_false_closure_breakpoint_count",
        "Terminal_false_closure_breakpoint_count", "FALSE_CLOSURE_BREAKPOINTS",
        "FALSE_CLOSURE_REGIONS", "Candidate_additional_target_contigs",
        "Candidate_additional_target_MAGs", "INTERNAL_FALSE_CLOSURE_EVIDENCE",
        "TERMINAL_FALSE_CLOSURE_EVIDENCE", "FALSE_CLOSURE_EVIDENCE", "FALSE_CLOSURE",
        "Min_connected_QNAMEs_per_side", "False_closure_reason",
    ]

    final_rows: List[Dict[str, object]] = []
    audit_rows: List[Dict[str, object]] = []
    for mag, row in source.items():
        contig_count = intval(row, "Contig_count", mag, 1)
        circular_count = intval(row, "SOURCE_CIRCULAR_CONTIG_COUNT", mag)
        resolved_count = intval(row, "SOURCE_CONTIG_CIRCULARITY_RESOLVED_COUNT", mag)
        unresolved_count = intval(row, "SOURCE_CONTIG_CIRCULARITY_UNRESOLVED_COUNT", mag)
        if resolved_count + unresolved_count != contig_count:
            fail(f"assembler circularity counts do not sum to contig count for {mag}")
        if circular_count > resolved_count:
            fail(f"circular contig count exceeds resolved count for {mag}")

        assembler = call(row["Assembler_reported_circularity"], {"Y", "N", "NA"}, "source-contig aggregate assembler circularity", mag)
        mag_level_input = call(row["MAG_level_input_assembler_circularity"], {"Y", "N", "NA"}, "MAG-level input assembler circularity", mag)
        complete = call(row["SOURCE_CONTIG_CIRCULARITY_COMPLETE_Y_N"], {"Y", "N"}, "source circularity completeness", mag)
        expected_aggregate = "Y" if circular_count > 0 else ("N" if complete == "Y" else "NA")
        if assembler != expected_aggregate:
            fail(f"contig aggregate circularity inconsistent for {mag}: {assembler} vs {expected_aggregate}")
        if complete != ("Y" if unresolved_count == 0 else "N"):
            fail(f"assembler circularity completeness inconsistent for {mag}")
        analyzable = call(row["False_closure_analyzable"], {"Y", "N"}, "false-closure analyzability", mag)
        false_call = call(row["False_closure"], {"Y", "N", "NA"}, "false-closure call", mag)
        if (analyzable == "Y") != (false_call in {"Y", "N"}):
            fail(f"false-closure analyzability/call mismatch for {mag}: {analyzable}/{false_call}")

        evidence = call(row["False_closure_evidence"], {"Y", "N"}, "false-closure evidence", mag)
        if false_call == "Y" and evidence != "Y":
            fail(f"FALSE_CLOSURE=Y without structural evidence for {mag}")
        premature = call(row["Premature_circularization_candidate"], {"Y", "N"}, "premature candidate", mag)
        proxy = call(row["Paper_small_genome_proxy"], {"Y", "N"}, "paper proxy", mag)
        screen = call(row["Premature_circularization_screen"], {"Y", "N"}, "premature screen", mag)
        expected_premature = "Y" if contig_count == 1 and assembler == "Y" and evidence == "Y" else "N"
        expected_screen = "Y" if expected_premature == "Y" or proxy == "Y" else "N"
        if premature != expected_premature or screen != expected_screen:
            fail(f"premature-circularization logic inconsistent for {mag}")

        validated = intval(row, "Validated_clipping_breakpoint_count", mag)
        internal_validated = intval(row, "Validated_internal_clipping_breakpoint_count", mag)
        terminal_validated = intval(row, "Validated_terminal_clipping_breakpoint_count", mag)
        false_n = intval(row, "False_closure_breakpoint_count", mag)
        internal_false_n = intval(row, "Internal_false_closure_breakpoint_count", mag)
        terminal_false_n = intval(row, "Terminal_false_closure_breakpoint_count", mag)
        if internal_validated + terminal_validated != validated:
            fail(f"validated clipping regional counts inconsistent for {mag}")
        if internal_false_n + terminal_false_n != false_n or false_n > validated:
            fail(f"false-closure breakpoint counts inconsistent for {mag}")
        if (false_n > 0) != (evidence == "Y"):
            fail(f"false-closure evidence/count mismatch for {mag}")

        large_an = call(row["Large_insertion_sensitivity_analyzable"], {"Y", "N"}, "large insertion analyzability", mag)
        large_call = call(row["Large_insertion_sensitivity_support"], {"Y", "N", "NA"}, "large insertion support", mag)
        if (large_an == "Y") != (large_call in {"Y", "N"}):
            fail(f"large-insertion analyzability/call mismatch for {mag}")

        final = {
            "ID": mag,
            "MAG_CONTIG_COUNT": contig_count,
            "FALSE_CLOSURE_ANALYZABLE": analyzable,
            "ASSEMBLER_REPORTED_CIRCULARITY": assembler,
            "MAG_LEVEL_INPUT_ASSEMBLER_CIRCULARITY": mag_level_input,
            "SOURCE_CIRCULAR_CONTIG_COUNT": circular_count,
            "SOURCE_CIRCULAR_CONTIGS": row["SOURCE_CIRCULAR_CONTIGS"],
            "SOURCE_CONTIG_CIRCULARITY_COMPLETE_Y_N": complete,
            "Paper_ribosomal_protein_family_count": row["Paper_ribosomal_protein_family_count"],
            "Paper_ribosomal_protein_families": row["Paper_ribosomal_protein_families"],
            "Paper_small_genome_length_criterion_Y_N": row["Paper_small_genome_length_criterion_Y_N"],
            "Paper_small_genome_ribosomal_criterion_Y_N": row["Paper_small_genome_ribosomal_criterion_Y_N"],
            "Paper_small_genome_proxy": proxy,
            "Paper_small_genome_proxy_contigs": row["Paper_small_genome_proxy_contigs"],
            "Premature_circularization_candidate": premature,
            "Premature_circularization_screen": screen,
            "LARGE_INSERTION_SENSITIVITY_ANALYZABLE": large_an,
            "LARGE_INSERTION_SENSITIVITY_SUPPORT": large_call,
            "LARGE_INSERTION_SENSITIVITY_SITE_COUNT": row["Large_insertion_sensitivity_site_count"],
            "LARGE_INSERTION_SENSITIVITY_SITES": row["Large_insertion_sensitivity_sites"],
            "Candidate_additional_target_contigs": row["Candidate_additional_target_contigs"],
            "Candidate_additional_target_MAGs": row["Candidate_additional_target_MAGs"],
            "VALIDATED_CLIPPING_PRESENT": "Y" if validated else "N",
            "INTERNAL_VALIDATED_CLIPPING_PRESENT": "Y" if internal_validated else "N",
            "TERMINAL_VALIDATED_CLIPPING_PRESENT": "Y" if terminal_validated else "N",
            "INTERNAL_FALSE_CLOSURE_EVIDENCE": row["Internal_false_closure_evidence"],
            "TERMINAL_FALSE_CLOSURE_EVIDENCE": row["Terminal_false_closure_evidence"],
            "FALSE_CLOSURE_EVIDENCE": evidence,
            "FALSE_CLOSURE": false_call,
            "FALSE_CLOSURE_BREAKPOINTS": row["False_closure_breakpoints"],
            "FALSE_CLOSURE_REGIONS": row["False_closure_regions"],
        }
        audit = {
            **final,
            "Contig_count": row["Contig_count"], "Contig": row["Contig"],
            "Genome_length_bp": row["Genome_length_bp"],
            "SOURCE_CONTIG_CIRCULARITY_RESOLVED_COUNT": resolved_count,
            "SOURCE_CONTIG_CIRCULARITY_UNRESOLVED_COUNT": unresolved_count,
            "Large_insertion_sensitivity_min_bp": row["Large_insertion_sensitivity_min_bp"],
            "Large_insertion_sensitivity_min_QNAMEs": row["Large_insertion_sensitivity_min_QNAMEs"],
            "All_MAG_validated_clipping_breakpoint_count": row["All_MAG_validated_clipping_breakpoint_count"],
            "Validated_clipping_breakpoint_count": validated,
            "Validated_internal_clipping_breakpoint_count": internal_validated,
            "Validated_terminal_clipping_breakpoint_count": terminal_validated,
            "False_closure_breakpoint_count": false_n,
            "Internal_false_closure_breakpoint_count": internal_false_n,
            "Terminal_false_closure_breakpoint_count": terminal_false_n,
            "Min_connected_QNAMEs_per_side": row["Min_connected_QNAMEs_per_side"],
            "False_closure_reason": row["False_closure_reason"],
        }
        final_rows.append(final)
        audit_rows.append(audit)

    atomic_dict_table(args.false_final, final_fields, final_rows)
    atomic_dict_table(args.false_audit, audit_fields, audit_rows)
    with open(args.metadata, "w", encoding="utf-8") as handle:
        handle.write(f"summarizer_version\t{VERSION}\n")
        handle.write("false_closure_main_reference\tsource_assembler_original_total_assembly\n")
        handle.write("circularity_gate\tsource_contig_assembler_reported_circularity\n")
        handle.write("final_assembler_reported_circularity\tsource_contig_aggregate_Y_if_any_source_contig_Y_N_if_all_resolved_N_else_NA\n")
        handle.write("MAG_level_input_circularity\taudit_and_single_contig_fallback_only\n")
        handle.write("multi_contig_MAG_policy\tFALSE_CLOSURE_Y_if_any_source_circular_contig_supports_false_closure; premature-circularization candidate limited to single-contig MAGs\n")
        handle.write("large_insertion_evidence\tsensitivity_only_no_change_to_false_closure_call\n")
        handle.write("paper_small_genome_proxy_changes_false_closure_call\tN\n")
    print(f"{VERSION}: MAGs={len(source)}")
    print(f"Final: {args.false_final}")
    print(f"Audit: {args.false_audit}")


if __name__ == "__main__":
    main()
