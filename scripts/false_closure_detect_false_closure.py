#!/usr/bin/env python3
"""Detect HiFi-supported false-closure evidence from validated clipping sites.

The false-closure module consumes the ALL-region clipping table produced by
03_clipping.sh and evaluates every final HiFi Clipping=Y site. Internal and
terminal sites are analyzed identically; region_type is retained only as an
annotation in the outputs. Exact physical reference boundaries cannot be
Clipping=Y upstream because a two-sided continuity test is impossible there.

For each validated site, clipped QNAMEs are recovered from the SAME source
assembler's HiFi BAM used for clipping analysis. False-closure structural
evidence uses primary + supplementary alignments only; secondary alignments are
excluded. A cross-contig structural connection requires a split alignment
involving at least one supplementary segment. QNAMEs are counted once per
breakpoint side. The two sides may connect to different contigs.

"""
from __future__ import annotations

import argparse
import csv
import os
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, NoReturn, Sequence, Set, Tuple

from false_closure_scg_utils import load_ribosomal_families
from false_closure_contig_connection_utils import (
    Alignment,
    atomic_tsv,
    cigar_operations,
    load_alignments_for_qnames,
    load_contig_to_mag,
    parse_alignment,
    representative,
    stream_sam,
    summarize_external_connection,
    target_count_text,
)

VERSION = "detect_false_closure"


@dataclass(frozen=True)
class MagInfo:
    mag_id: str
    contigs: Tuple[str, ...]
    genome_length: int

    @property
    def contig_count(self) -> int:
        return len(self.contigs)

    @property
    def contig(self) -> str:
        return self.contigs[0] if len(self.contigs) == 1 else "NA"

    @property
    def eligibility(self) -> str:
        return "Y" if len(self.contigs) == 1 else "N"


@dataclass(frozen=True)
class SourceMapping:
    assembler: str
    assembly_fasta: str
    hifi_all_bam: str
    hifi_primary_supplementary_bam: str
    lengths: Dict[str, int]


@dataclass(frozen=True)
class Breakpoint:
    mag_id: str
    source_assembler: str
    contig: str
    contig_length: int
    position_1based: int
    region_type: str
    distance_to_nearest_end_bp: int
    coverage: int
    clipping_count: int
    clipping_fraction: float
    n_all_clipped_qnames_from_validation: int
    n_exact_spanning_qnames_from_validation: int


@dataclass(frozen=True)
class LargeInsertionSignal:
    mag_id: str
    source_assembler: str
    contig: str
    position_1based: int
    region_type: str
    n_spanning_insertion_qnames: int
    insertion_length_bp_min: int
    insertion_length_bp_median: str
    insertion_length_bp_max: int


def fail(message: str) -> NoReturn:
    raise SystemExit(f"ERROR: {message}")


def yes(value: str) -> bool:
    return (value or "").strip().upper() in {"Y", "YES", "TRUE", "1"}


def load_resolved_contig_circularity(path: str) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, str]]:
    """Load assembly-wide source-contig circularity calls.

    The table may include contigs outside the final MAG catalogue. MAG-level
    circularity is retained only as compatibility metadata for MAG-assigned
    contigs and is never used as a structural false-closure gate.
    """
    calls: Dict[str, str] = {}
    methods: Dict[str, str] = {}
    mag_calls: Dict[str, str] = {}
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="	")
        required = {
            "MAG_ID", "Contig", "SOURCE_CONTIG_ASSEMBLER_CIRCULARITY",
            "Circularity_resolution_method", "MAG_level_assembler_circularity",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            fail(f"resolved contig circularity table lacks columns {sorted(missing)}: {path}")
        for line_no, row in enumerate(reader, start=2):
            mag = (row.get("MAG_ID") or "NA").strip() or "NA"
            contig = (row.get("Contig") or "").strip()
            call = (row.get("SOURCE_CONTIG_ASSEMBLER_CIRCULARITY") or "NA").strip().upper()
            method = (row.get("Circularity_resolution_method") or "").strip()
            mag_call = (row.get("MAG_level_assembler_circularity") or "NA").strip().upper()
            if not contig or contig in calls:
                fail(f"blank/duplicate contig in resolved circularity table at line {line_no}: {contig}")
            if call not in {"Y", "N", "NA"} or mag_call not in {"Y", "N", "NA"}:
                fail(f"invalid circularity call at {path}:{line_no}: {call}/{mag_call}")
            calls[contig] = call
            methods[contig] = method or "NA"
            if mag != "NA":
                previous = mag_calls.get(mag)
                if previous is not None and previous != mag_call:
                    fail(f"MAG-level circularity differs among contig rows for {mag}")
                mag_calls[mag] = mag_call
    if not calls:
        fail(f"resolved circularity table contains no contigs: {path}")
    return calls, methods, mag_calls


def load_fai_lengths(path: str) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    with open(path, encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) < 2:
                fail(f"malformed FAI row {line_number}: {path}")
            name = fields[0]
            try:
                length = int(fields[1])
            except ValueError as exc:
                fail(f"invalid FAI length at {path}:{line_number}: {exc}")
            if not name or length < 1 or name in lengths:
                fail(f"blank/duplicate contig in FAI: {name}")
            lengths[name] = length
    if not lengths:
        fail(f"FAI contains no sequences: {path}")
    return lengths


def load_mag_info(map_path: str, assemblies_fai: str) -> List[MagInfo]:
    lengths = load_fai_lengths(assemblies_fai)
    grouped: Dict[str, List[str]] = defaultdict(list)
    seen_contigs: Dict[str, str] = {}
    with open(map_path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        lower = {name.lower(): name for name in fields}
        mag_col = lower.get("mag_id") or lower.get("id")
        contig_col = lower.get("contig") or lower.get("original_contig")
        if not mag_col or not contig_col:
            fail("MAG contig map requires MAG_ID/ID and Contig columns")
        for line_number, row in enumerate(reader, start=2):
            mag_id = (row.get(mag_col) or "").strip()
            contig = (row.get(contig_col) or "").strip()
            if not mag_id or not contig:
                fail(f"blank MAG/contig in MAG contig map line {line_number}")
            if contig not in lengths:
                fail(f"MAG contig absent from all-assemblies FAI: {contig}")
            previous = seen_contigs.get(contig)
            if previous is not None and previous != mag_id:
                fail(f"contig assigned to multiple MAGs: {contig}")
            if previous is None:
                grouped[mag_id].append(contig)
                seen_contigs[contig] = mag_id
    if not grouped:
        fail("MAG contig map contains no MAGs")
    return [
        MagInfo(
            mag_id=mag_id,
            contigs=tuple(contigs),
            genome_length=sum(lengths[contig] for contig in contigs),
        )
        for mag_id, contigs in sorted(grouped.items())
    ]


def load_mapping_manifest(path: str) -> Dict[str, SourceMapping]:
    mappings: Dict[str, SourceMapping] = {}
    seen_contigs: Dict[str, str] = {}
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "Assembler", "Assembly_FASTA", "HiFi_all_BAM",
            "HiFi_primary_supplementary_BAM",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            fail(f"mapping manifest lacks columns {sorted(missing)}: {path}")
        for line_number, row in enumerate(reader, start=2):
            assembler = (row.get("Assembler") or "").strip()
            assembly = (row.get("Assembly_FASTA") or "").strip()
            bam = (row.get("HiFi_all_BAM") or "").strip()
            ps_bam = (row.get("HiFi_primary_supplementary_BAM") or "").strip()
            if not assembler or not assembly or not bam or not ps_bam or assembler in mappings:
                fail(f"invalid/duplicate mapping manifest row {line_number}")
            fai = assembly + ".fai"
            for required_path, label in (
                (assembly, "assembly"),
                (fai, "FAI"),
                (bam, "HiFi all-mapped BAM"),
                (ps_bam, "HiFi primary+supplementary BAM"),
            ):
                if not os.path.isfile(required_path) or os.path.getsize(required_path) == 0:
                    fail(f"{label} missing for {assembler}: {required_path}")
            lengths = load_fai_lengths(fai)
            for contig in lengths:
                previous = seen_contigs.get(contig)
                if previous is not None:
                    fail(
                        "contig ID occurs in multiple source assemblers: "
                        f"{contig} ({previous}, {assembler})"
                    )
                seen_contigs[contig] = assembler
            mappings[assembler] = SourceMapping(assembler, assembly, bam, ps_bam, lengths)
    if not mappings:
        fail("mapping manifest contains no assembler rows")
    return mappings


def load_breakpoints(
    path: str,
    mappings: Dict[str, SourceMapping],
    contig_to_mag: Dict[str, str],
) -> Tuple[List[Breakpoint], Set[str]]:
    """Load assembly-wide final HiFi Clipping=Y sites.

    All source-assembly contigs are eligible. MAG membership is carried only as
    an annotation so the contig-level false-closure call is independent of
    binning. A circular contig with any clipping candidate for which
    Clipping_analyzable=N cannot receive a definitive negative call.
    """
    source_of = {contig: assembler for assembler, src in mappings.items() for contig in src.lengths}
    out: List[Breakpoint] = []
    unresolved_contigs: Set[str] = set()
    seen: Set[Tuple[str, int]] = set()
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="	")
        fields = set(reader.fieldnames or [])
        required = {
            "source_assembler", "contig", "length", "pos_1based",
            "region_type", "distance_to_nearest_end_bp", "reference_boundary",
            "coverage", "clipping_count", "clipping_ratio",
            "Clipping_analyzable", "Clipping",
            "validated_clipping", "n_unique_all_clipped_qnames",
        }
        missing = required - fields
        if missing:
            fail(f"all-region clipping table lacks columns {sorted(missing)}: {path}")
        spanning_column = (
            "n_unique_exact_spanning_qnames"
            if "n_unique_exact_spanning_qnames" in fields
            else ("n_unique_clean_spanning_qnames" if "n_unique_clean_spanning_qnames" in fields else None)
        )
        if spanning_column is None:
            fail(
                "all-region clipping table lacks a spanning-QNAME column; "
                f"expected n_unique_exact_spanning_qnames: {path}"
            )
        for line_number, row in enumerate(reader, start=2):
            contig = (row.get("contig") or "").strip()
            assembler = (row.get("source_assembler") or "").strip()
            if contig not in source_of:
                fail(f"clipping contig absent from assemblies at line {line_number}: {contig}")
            if source_of[contig] != assembler:
                fail(f"assembler mismatch for clipping contig {contig}: {assembler} vs {source_of[contig]}")
            clipping_analyzable = (row.get("Clipping_analyzable") or "").strip().upper()
            clipping_call = (row.get("Clipping") or "NA").strip().upper()
            if clipping_analyzable not in {"Y", "N"}:
                fail(f"invalid Clipping_analyzable at line {line_number}: {clipping_analyzable}")
            if clipping_analyzable == "N":
                if clipping_call != "NA":
                    fail(f"Clipping_analyzable=N requires Clipping=NA at line {line_number}")
                unresolved_contigs.add(contig)
                continue
            if clipping_call not in {"Y", "N"}:
                fail(f"Clipping_analyzable=Y requires Clipping=Y/N at line {line_number}")
            if clipping_call != "Y":
                continue
            if not yes(row.get("validated_clipping", "")):
                fail(
                    "Clipping=Y is inconsistent with validated_clipping=N at "
                    f"{contig}:{row.get('pos_1based', '?')}"
                )
            if yes(row.get("reference_boundary", "")):
                fail(
                    "validated_clipping=Y is inconsistent with reference_boundary=Y at "
                    f"{contig}:{row.get('pos_1based', '?')}"
                )
            source = mappings[assembler]
            region = (row.get("region_type") or "").strip().upper()
            if region not in {"INTERNAL", "TERMINAL"}:
                fail(f"invalid region_type at line {line_number}: {region}")
            try:
                length = int(row["length"]); position = int(row["pos_1based"]); distance = int(row["distance_to_nearest_end_bp"])
                coverage = int(row["coverage"]); clipping_count = int(row["clipping_count"]); clipping_fraction = float(row["clipping_ratio"])
                n_clipped = int(row["n_unique_all_clipped_qnames"]); n_spanning = int(row[spanning_column])
            except (TypeError, ValueError) as exc:
                fail(f"invalid clipping row at line {line_number}: {exc}")
            if length != source.lengths[contig]:
                fail(f"clipping length disagrees with the assembler FAI for {contig}")
            if not 1 <= position <= length or distance < 0:
                fail(f"invalid clipping coordinate for {contig}: {position}")
            if n_clipped < 2 or n_spanning != 0:
                fail(
                    "Clipping=Y must mean >=2 clipped QNAMEs and zero exact "
                    f"spanning QNAMEs: {contig}:{position}"
                )
            key = (contig, position)
            if key in seen:
                fail(f"duplicate validated clipping breakpoint: {contig}:{position}")
            seen.add(key)
            out.append(Breakpoint(
                mag_id=contig_to_mag.get(contig, "NA"), source_assembler=assembler, contig=contig,
                contig_length=length, position_1based=position, region_type=region,
                distance_to_nearest_end_bp=distance, coverage=coverage,
                clipping_count=clipping_count, clipping_fraction=clipping_fraction,
                n_all_clipped_qnames_from_validation=n_clipped,
                n_exact_spanning_qnames_from_validation=n_spanning,
            ))
    return (
        sorted(out, key=lambda b: (b.source_assembler, b.contig, b.position_1based)),
        unresolved_contigs,
    )


def load_large_insertion_signals(
    path: str,
    contig_to_mag: Dict[str, str],
    min_insertion_bp: int,
    min_qnames: int,
) -> List[LargeInsertionSignal]:
    """Load the assembly-wide large read-side insertion sensitivity screen."""
    out: List[LargeInsertionSignal] = []
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="	")
        fields = set(reader.fieldnames or [])
        required = {
            "source_assembler", "contig", "pos_1based", "region_type",
            "Clipping_analyzable", "n_unique_spanning_with_I_qnames",
            "spanning_I_length_bp_min", "spanning_I_length_bp_median",
            "spanning_I_length_bp_max",
        }
        missing = required - fields
        if missing:
            fail(
                "all-region clipping table lacks large-insertion audit columns "
                f"{sorted(missing)}: {path}"
            )
        for line_number, row in enumerate(reader, start=2):
            contig = (row.get("contig") or "").strip()
            if (row.get("Clipping_analyzable") or "").strip().upper() != "Y":
                continue
            try:
                n_qnames = int(row.get("n_unique_spanning_with_I_qnames") or 0)
            except ValueError as exc:
                fail(f"invalid insertion-QNAME count at line {line_number}: {exc}")
            if n_qnames < min_qnames:
                continue
            min_text = (row.get("spanning_I_length_bp_min") or "NA").strip()
            median_text = (row.get("spanning_I_length_bp_median") or "NA").strip()
            max_text = (row.get("spanning_I_length_bp_max") or "NA").strip()
            if min_text.upper() == "NA" or max_text.upper() == "NA":
                fail(f"insertion-bearing spanning QNAMEs lack insertion-length audit at line {line_number}")
            try:
                min_len = int(float(min_text)); max_len = int(float(max_text)); position = int(row["pos_1based"])
            except (TypeError, ValueError) as exc:
                fail(f"invalid large-insertion audit row at line {line_number}: {exc}")
            if min_len < min_insertion_bp:
                continue
            region = (row.get("region_type") or "").strip().upper()
            if region not in {"INTERNAL", "TERMINAL"}:
                fail(f"invalid region_type at line {line_number}: {region}")
            out.append(LargeInsertionSignal(
                mag_id=contig_to_mag.get(contig, "NA"),
                source_assembler=(row.get("source_assembler") or "").strip(),
                contig=contig, position_1based=position, region_type=region,
                n_spanning_insertion_qnames=n_qnames, insertion_length_bp_min=min_len,
                insertion_length_bp_median=median_text, insertion_length_bp_max=max_len,
            ))
    return sorted(out, key=lambda x: (x.source_assembler, x.contig, x.position_1based))


def load_breakpoint_sources(
    bam: str, breakpoints: Sequence[Breakpoint]
) -> Dict[Tuple[str, int], Dict[str, Dict[str, Alignment]]]:
    """Recover exact-position clipped QNAMEs from one source-assembler BAM."""
    if not breakpoints:
        return {}
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".bed") as bed:
        bed_path = bed.name
        for bp in breakpoints:
            pos0 = bp.position_1based - 1
            bed.write(f"{bp.contig}\t{pos0}\t{pos0 + 1}\n")

    by_contig: Dict[str, Dict[str, List[Alignment]]] = defaultdict(
        lambda: defaultdict(list)
    )
    try:
        for line in stream_sam(
            ["samtools", "view", "-F", "4", "-q", "0", "-L", bed_path, bam]
        ):
            aln = parse_alignment(line)
            if aln is not None:
                by_contig[aln.rname][aln.qname].append(aln)
    finally:
        os.unlink(bed_path)

    result: Dict[Tuple[str, int], Dict[str, Dict[str, Alignment]]] = {}
    for bp in breakpoints:
        pos0 = bp.position_1based - 1
        sides: Dict[str, Dict[str, Alignment]] = {"left": {}, "right": {}}
        for qname, alignments in by_contig.get(bp.contig, {}).items():
            left_candidates = [
                aln
                for aln in alignments
                if aln.right_clip > 0 and (aln.rend - 1) == pos0
            ]
            right_candidates = [
                aln
                for aln in alignments
                if aln.left_clip > 0 and aln.rstart == pos0
            ]
            rep = representative(left_candidates)
            if rep is not None:
                sides["left"][qname] = rep
            rep = representative(right_candidates)
            if rep is not None:
                sides["right"][qname] = rep
        recovered = set(sides["left"]) | set(sides["right"])
        if len(recovered) != bp.n_all_clipped_qnames_from_validation:
            fail(
                "exact clipped-QNAME recovery disagrees with 03_clipping.sh validation at "
                f"{bp.contig}:{bp.position_1based}; "
                f"validation={bp.n_all_clipped_qnames_from_validation}, recovered={len(recovered)}"
            )
        result[(bp.contig, bp.position_1based)] = sides
    return result


def false_closure_connection_rule(
    left_connected_qnames: int,
    right_connected_qnames: int,
    min_connected_qnames_per_side: int,
) -> bool:
    return (
        left_connected_qnames >= min_connected_qnames_per_side
        and right_connected_qnames >= min_connected_qnames_per_side
    )


def region_call(records: Iterable[Tuple[Breakpoint, bool, Dict[str, Set[str]], Dict[str, Set[str]]]], region: str) -> bool:
    return any(record[1] for record in records if record[0].region_type == region)


def run(args: argparse.Namespace) -> None:
    mag_info = load_mag_info(args.mag_contig_map, args.assemblies_fai)
    mag_by_id = {row.mag_id: row for row in mag_info}
    combined_lengths = load_fai_lengths(args.assemblies_fai)
    ribosomal_families = load_ribosomal_families(args.scg_tsv)
    contig_to_mag = load_contig_to_mag(args.mag_contig_map)
    mappings = load_mapping_manifest(args.mapping_manifest)
    contig_circularity, contig_circularity_method, assembler_calls = load_resolved_contig_circularity(
        args.contig_circularity
    )
    source_contigs = {contig for source in mappings.values() for contig in source.lengths}
    if set(contig_circularity) != source_contigs:
        missing = sorted(source_contigs - set(contig_circularity))[:10]
        extra = sorted(set(contig_circularity) - source_contigs)[:10]
        fail(f"resolved circularity/assembly contig mismatch: missing={missing}, extra={extra}")
    breakpoints, unresolved_clipping_contigs = load_breakpoints(
        args.clipping_all_final_tsv, mappings, contig_to_mag
    )
    large_insertion_signals = load_large_insertion_signals(
        args.clipping_all_final_tsv, contig_to_mag,
        args.large_insertion_min_bp, args.large_insertion_min_qnames,
    )
    large_insertion_by_mag: Dict[str, List[LargeInsertionSignal]] = defaultdict(list)
    for signal in large_insertion_signals:
        large_insertion_by_mag[signal.mag_id].append(signal)

    if args.large_insertion_details:
        sensitivity_rows: List[List[object]] = []
        for signal in large_insertion_signals:
            mag = mag_by_id.get(signal.mag_id)
            mag_call = assembler_calls.get(signal.mag_id, "NA") if signal.mag_id != "NA" else "NA"
            contig_call = contig_circularity.get(signal.contig, "NA")
            sensitivity_analyzable = contig_call == "Y"
            sensitivity_rows.append([
                signal.mag_id, signal.source_assembler, signal.contig,
                signal.position_1based, signal.region_type,
                signal.n_spanning_insertion_qnames,
                signal.insertion_length_bp_min, signal.insertion_length_bp_median,
                signal.insertion_length_bp_max, mag.contig_count if mag is not None else 0, contig_call, mag_call,
                contig_circularity_method.get(signal.contig, "NA"),
                "Y" if sensitivity_analyzable else "N",
                "Y" if sensitivity_analyzable else "NA",
                ("SOURCE_CIRCULAR_CONTIG_LARGE_INSERTION_SENSITIVITY"
                 if sensitivity_analyzable
                 else ("SOURCE_CONTIG_CIRCULARITY_UNRESOLVED" if contig_call == "NA"
                       else "SOURCE_CONTIG_NOT_REPORTED_CIRCULAR")),
            ])
        atomic_tsv(
            args.large_insertion_details,
            [
                "ID", "Source_assembler", "Contig", "Breakpoint_1based",
                "Region_type", "Spanning_insertion_QNAMEs",
                "Insertion_length_bp_min", "Insertion_length_bp_median",
                "Insertion_length_bp_max", "Contig_count",
                "SOURCE_CONTIG_ASSEMBLER_CIRCULARITY", "MAG_level_assembler_circularity",
                "Circularity_resolution_method",
                "LARGE_INSERTION_SENSITIVITY_ANALYZABLE_Y_N",
                "LARGE_INSERTION_SENSITIVITY_SUPPORT_Y_N_NA", "Sensitivity_reason",
            ],
            sensitivity_rows,
        )

    sources: Dict[Tuple[str, int], Dict[str, Dict[str, Alignment]]] = {}
    primary_supplementary_by_source_read: Dict[str, Dict[str, List[Alignment]]] = {}
    by_assembler: Dict[str, List[Breakpoint]] = defaultdict(list)
    for bp in breakpoints:
        by_assembler[bp.source_assembler].append(bp)

    for assembler, assembler_breakpoints in sorted(by_assembler.items()):
        bam = mappings[assembler].hifi_all_bam
        source_rows = load_breakpoint_sources(bam, assembler_breakpoints)
        sources.update(source_rows)
        candidate_qnames: Set[str] = set()
        for sides in source_rows.values():
            candidate_qnames.update(sides.get("left", {}))
            candidate_qnames.update(sides.get("right", {}))
        primary_supplementary_by_source_read[assembler] = load_alignments_for_qnames(
            mappings[assembler].hifi_primary_supplementary_bam, candidate_qnames
        )


    if args.candidate_read_output:
        candidate_rows: List[List[object]] = []
        for bp in breakpoints:
            sides = sources.get((bp.contig, bp.position_1based), {})
            for side in ("left", "right"):
                for qname in sorted(sides.get(side, {})):
                    candidate_rows.append(
                        [
                            bp.mag_id,
                            bp.source_assembler,
                            bp.contig,
                            bp.position_1based,
                            bp.region_type,
                            side.upper(),
                            qname,
                        ]
                    )
        atomic_tsv(
            args.candidate_read_output,
            [
                "ID", "Source_assembler", "Contig", "Breakpoint_1based",
                "Region_type", "Breakpoint_side", "Read_ID",
            ],
            candidate_rows,
        )

    detail_header = [
        "ID", "Source_assembler", "Contig", "Breakpoint_1based",
        "Region_type", "Distance_to_nearest_end_bp",
        "Validated_clipping", "Coverage", "Clipping_count", "Clipping_fraction",
        "Validated_all_clipped_QNAMEs", "Validated_clean_spanning_QNAMEs",
        "HIFI_PRIMARY_SUPPLEMENTARY_LEFT_CLIPPED_UNIQUE_READS",
        "HIFI_PRIMARY_SUPPLEMENTARY_LEFT_CROSS_CONTIG_CONNECTED_UNIQUE_READS",
        "HIFI_PRIMARY_SUPPLEMENTARY_LEFT_TARGET_CONTIGS",
        "HIFI_PRIMARY_SUPPLEMENTARY_RIGHT_CLIPPED_UNIQUE_READS",
        "HIFI_PRIMARY_SUPPLEMENTARY_RIGHT_CROSS_CONTIG_CONNECTED_UNIQUE_READS",
        "HIFI_PRIMARY_SUPPLEMENTARY_RIGHT_TARGET_CONTIGS",
        "HIFI_PRIMARY_SUPPLEMENTARY_BILATERAL_CROSS_CONTIG_CONNECTION_SUPPORTED_Y_N",
        "Left_clipped_QNAMEs", "Left_cross_contig_connected_QNAMEs",
        "Left_target_contigs", "Right_clipped_QNAMEs",
        "Right_cross_contig_connected_QNAMEs", "Right_target_contigs",
        "Bilateral_cross_contig_connection_evidence", "False_closure_evidence",
        "Assembler_reported_circularity", "MAG_level_assembler_reported_circularity",
        "Circularity_resolution_method", "False_closure_analyzable",
        "HIFI_PRIMARY_SUPPLEMENTARY_FALSE_CLOSURE_Y_N",
        "False_closure", "False_closure_reason",
    ]

    detail_rows: List[List[object]] = []
    records: List[Tuple] = []
    for bp in breakpoints:
        sides = sources.get(
            (bp.contig, bp.position_1based), {"left": {}, "right": {}}
        )
        left = sides.get("left", {})
        right = sides.get("right", {})
        ps_left = {q: a for q, a in left.items() if not a.secondary}
        ps_right = {q: a for q, a in right.items() if not a.secondary}
        primary_supplementary_by_read = primary_supplementary_by_source_read.get(
            bp.source_assembler, {}
        )

        left_connected, left_targets = summarize_external_connection(
            ps_left, primary_supplementary_by_read, bp.contig
        )
        right_connected, right_targets = summarize_external_connection(
            ps_right, primary_supplementary_by_read, bp.contig
        )
        evidence = false_closure_connection_rule(
            len(left_connected),
            len(right_connected),
            args.min_linked_reads_per_side,
        )

        mag_call = assembler_calls.get(bp.mag_id, "NA") if bp.mag_id != "NA" else "NA"
        contig_call = contig_circularity.get(bp.contig, "NA")
        circularity_method = contig_circularity_method.get(bp.contig, "NA")
        false_closure_analyzable = (
            contig_call == "Y"
            and bp.contig not in unresolved_clipping_contigs
        )
        false_closure = (
            "Y" if false_closure_analyzable and evidence
            else ("N" if false_closure_analyzable else "NA")
        )

        if contig_call == "NA":
            reason = "SOURCE_CONTIG_CIRCULARITY_UNRESOLVED"
        elif contig_call != "Y":
            reason = "SOURCE_CONTIG_NOT_REPORTED_CIRCULAR"
        elif bp.contig in unresolved_clipping_contigs:
            reason = "UNRESOLVED_HIFI_CLIPPING_CANDIDATE_ON_CIRCULAR_CONTIG"
        elif len(left_connected) < args.min_linked_reads_per_side:
            reason = "INSUFFICIENT_LEFT_PRIMARY_SUPPLEMENTARY_CROSS_CONTIG_CONNECTED_READS"
        elif len(right_connected) < args.min_linked_reads_per_side:
            reason = "INSUFFICIENT_RIGHT_PRIMARY_SUPPLEMENTARY_CROSS_CONTIG_CONNECTED_READS"
        else:
            reason = "VALIDATED_CLIPPING_WITH_BILATERAL_PRIMARY_SUPPLEMENTARY_CROSS_CONTIG_SPLIT_CONNECTIONS"

        detail_rows.append([
            bp.mag_id,
            bp.source_assembler,
            bp.contig,
            bp.position_1based,
            bp.region_type,
            bp.distance_to_nearest_end_bp,
            "Y",
            bp.coverage,
            bp.clipping_count,
            f"{bp.clipping_fraction:.6f}",
            bp.n_all_clipped_qnames_from_validation,
            bp.n_exact_spanning_qnames_from_validation,
            len(ps_left), len(left_connected), target_count_text(left_targets),
            len(ps_right), len(right_connected), target_count_text(right_targets),
            "Y" if evidence else "N",
            len(ps_left), len(left_connected), target_count_text(left_targets),
            len(ps_right), len(right_connected), target_count_text(right_targets),
            "Y" if evidence else "N", "Y" if evidence else "N",
            contig_call, mag_call, circularity_method,
            "Y" if false_closure_analyzable else "N",
            false_closure, false_closure, reason,
        ])
        records.append((bp, evidence, left_targets, right_targets))

    atomic_tsv(args.details, detail_header, detail_rows)

    by_mag: Dict[str, List[Tuple]] = defaultdict(list)
    for record in records:
        by_mag[record[0].mag_id].append(record)

    summary_header = [
        "ID", "Contig_count", "Contig", "Genome_length_bp",
        "Assembler_reported_circularity",
        "MAG_level_input_assembler_circularity",
        "SOURCE_CIRCULAR_CONTIG_COUNT", "SOURCE_CIRCULAR_CONTIGS",
        "SOURCE_CONTIG_CIRCULARITY_RESOLVED_COUNT",
        "SOURCE_CONTIG_CIRCULARITY_UNRESOLVED_COUNT",
        "SOURCE_CONTIG_CIRCULARITY_COMPLETE_Y_N",
        "False_closure_analyzable",
        "Paper_ribosomal_protein_family_count",
        "Paper_ribosomal_protein_families",
        "Paper_small_genome_length_criterion_Y_N",
        "Paper_small_genome_ribosomal_criterion_Y_N",
        "Paper_small_genome_proxy", "Paper_small_genome_proxy_contigs",
        "Premature_circularization_candidate",
        "Premature_circularization_screen",
        "Large_insertion_sensitivity_analyzable",
        "Large_insertion_sensitivity_support",
        "Large_insertion_sensitivity_site_count",
        "Large_insertion_sensitivity_sites",
        "Large_insertion_sensitivity_min_bp",
        "Large_insertion_sensitivity_min_QNAMEs",
        "All_MAG_validated_clipping_breakpoint_count",
        "Validated_clipping_breakpoint_count",
        "Validated_internal_clipping_breakpoint_count",
        "Validated_terminal_clipping_breakpoint_count",
        "HIFI_PRIMARY_SUPPLEMENTARY_FALSE_CLOSURE_BREAKPOINT_COUNT",
        "HIFI_PRIMARY_SUPPLEMENTARY_FALSE_CLOSURE_Y_N",
        "False_closure_breakpoint_count",
        "Internal_false_closure_breakpoint_count",
        "Terminal_false_closure_breakpoint_count",
        "False_closure_breakpoints", "False_closure_regions",
        "Candidate_additional_target_contigs", "Candidate_additional_target_MAGs",
        "Internal_false_closure_evidence", "Terminal_false_closure_evidence",
        "False_closure_evidence", "False_closure",
        "Min_connected_QNAMEs_per_side", "False_closure_reason",
    ]
    summary_rows: List[List[object]] = []
    for mag in mag_info:
        recs = by_mag.get(mag.mag_id, [])
        circular_contigs = [c for c in mag.contigs if contig_circularity.get(c, "NA") == "Y"]
        unresolved_circularity_contigs = [c for c in mag.contigs if contig_circularity.get(c, "NA") == "NA"]
        resolved_circularity_count = mag.contig_count - len(unresolved_circularity_contigs)
        circularity_complete = not unresolved_circularity_contigs
        circular_recs = [r for r in recs if contig_circularity.get(r[0].contig, "NA") == "Y"]
        positives = [r for r in circular_recs if r[1]]
        internal = [r for r in circular_recs if r[0].region_type == "INTERNAL"]
        terminal = [r for r in circular_recs if r[0].region_type == "TERMINAL"]
        internal_pos = [r for r in positives if r[0].region_type == "INTERNAL"]
        terminal_pos = [r for r in positives if r[0].region_type == "TERMINAL"]

        circular_clipping_unresolved = any(c in unresolved_clipping_contigs for c in circular_contigs)
        definitive_negative_possible = (
            bool(circular_contigs)
            and circularity_complete
            and not circular_clipping_unresolved
        )
        evidence = bool(positives)
        false_closure = "Y" if evidence else ("N" if definitive_negative_possible else "NA")
        analyzable = false_closure in {"Y", "N"}

        positions = (
            ";".join(
                f"{r[0].contig}|{r[0].region_type}:{r[0].position_1based}"
                for r in positives
            ) or "NA"
        )
        regions = ";".join(sorted({r[0].region_type for r in positives})) or "NA"
        targets: Set[str] = set()
        for r in positives:
            targets.update(r[2].keys())
            targets.update(r[3].keys())
        target_contigs = ";".join(sorted(targets)) if targets else "NA"
        target_mags = (
            ";".join(sorted({contig_to_mag[t] for t in targets if t in contig_to_mag})) or "NA"
        )

        large_sites = [
            site for site in large_insertion_by_mag.get(mag.mag_id, [])
            if contig_circularity.get(site.contig, "NA") == "Y"
        ]
        if large_sites:
            large_insertion_support = "Y"
        elif definitive_negative_possible:
            large_insertion_support = "N"
        else:
            large_insertion_support = "NA"
        large_insertion_analyzable = large_insertion_support in {"Y", "N"}
        large_insertion_sites_text = (
            ";".join(
                f"{site.contig}|{site.region_type}:{site.position_1based}"
                f"[n={site.n_spanning_insertion_qnames},minI={site.insertion_length_bp_min}]"
                for site in large_sites
            ) if large_sites else "NA"
        )


        circular_family_map = {c: ribosomal_families.get(c, ()) for c in circular_contigs}
        paper_family_count = max((len(v) for v in circular_family_map.values()), default=0)
        max_family_contigs = [c for c, v in circular_family_map.items() if len(v) == paper_family_count and paper_family_count > 0]
        paper_families = tuple(sorted({f for c in max_family_contigs for f in circular_family_map[c]}))
        paper_length_ok = any(combined_lengths[c] < args.premature_max_length_bp for c in circular_contigs)
        paper_ribosomal_ok = any(
            len(circular_family_map[c]) >= args.premature_min_ribosomal_families for c in circular_contigs
        )
        paper_proxy_contigs = [
            c for c in circular_contigs
            if combined_lengths[c] < args.premature_max_length_bp
            and len(circular_family_map[c]) >= args.premature_min_ribosomal_families
        ]
        paper_proxy = bool(paper_proxy_contigs)
        source_aggregate_call = (
            "Y" if circular_contigs else ("N" if circularity_complete else "NA")
        )
        mag_level_input_call = assembler_calls.get(mag.mag_id, "NA")
        premature_candidate = mag.contig_count == 1 and source_aggregate_call == "Y" and evidence
        premature_screen = premature_candidate or paper_proxy

        if evidence:
            reason = "SOURCE_CIRCULAR_CONTIG_PRIMARY_SUPPLEMENTARY_FALSE_CLOSURE_STRUCTURE_SUPPORTED"
        elif not circular_contigs and unresolved_circularity_contigs:
            reason = "SOURCE_CONTIG_CIRCULARITY_UNRESOLVED"
        elif not circular_contigs:
            reason = "NOT_APPLICABLE_NO_SOURCE_CIRCULAR_CONTIG"
        elif not circularity_complete:
            reason = "SOURCE_CONTIG_CIRCULARITY_INCOMPLETE_FOR_DEFINITIVE_NEGATIVE"
        elif circular_clipping_unresolved:
            reason = "UNRESOLVED_HIFI_CLIPPING_CANDIDATE_ON_SOURCE_CIRCULAR_CONTIG"
        elif not circular_recs:
            reason = "ANALYZABLE_SOURCE_CIRCULAR_CONTIG_WITH_NO_HIFI_CLIPPING_Y_BREAKPOINT"
        else:
            reason = "ANALYZABLE_CIRCULAR_CONTIG_CLIPPING_LACKED_BILATERAL_PRIMARY_SUPPLEMENTARY_CROSS_CONTIG_SPLIT_CONNECTIONS"

        summary_rows.append([
            mag.mag_id, mag.contig_count, mag.contig, mag.genome_length,
            source_aggregate_call, mag_level_input_call,
            len(circular_contigs), ";".join(circular_contigs) if circular_contigs else "NA",
            resolved_circularity_count, len(unresolved_circularity_contigs),
            "Y" if circularity_complete else "N",
            "Y" if analyzable else "N",
            paper_family_count, ";".join(paper_families) if paper_families else "NA",
            "Y" if paper_length_ok else "N", "Y" if paper_ribosomal_ok else "N",
            "Y" if paper_proxy else "N",
            ";".join(paper_proxy_contigs) if paper_proxy_contigs else "NA",
            "Y" if premature_candidate else "N", "Y" if premature_screen else "N",
            "Y" if large_insertion_analyzable else "N", large_insertion_support,
            len(large_sites), large_insertion_sites_text,
            args.large_insertion_min_bp, args.large_insertion_min_qnames,
            len(recs), len(circular_recs), len(internal), len(terminal),
            len(positives), false_closure,
            len(positives), len(internal_pos), len(terminal_pos),
            positions, regions, target_contigs, target_mags,
            "Y" if internal_pos else "N", "Y" if terminal_pos else "N",
            "Y" if evidence else "N", false_closure,
            args.min_linked_reads_per_side, reason,
        ])

    atomic_tsv(args.output, summary_header, summary_rows)
    print(
        f"{VERSION}: validated_breakpoints={len(detail_rows)}; "
        f"false_closure_structure_Y={sum(record[1] for record in records)}"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assemblies-fai", required=True)
    parser.add_argument("--mapping-manifest", required=True)
    parser.add_argument("--clipping-all-final-tsv", required=True)
    parser.add_argument("--contig-circularity", required=True)
    parser.add_argument("--mag-contig-map", required=True)
    parser.add_argument("--scg-tsv", required=True)
    parser.add_argument("--premature-max-length-bp", type=int, default=500000)
    parser.add_argument("--premature-min-ribosomal-families", type=int, default=3)
    parser.add_argument("--min-linked-reads-per-side", type=int, default=2)
    parser.add_argument("--large-insertion-min-bp", type=int, default=100)
    parser.add_argument("--large-insertion-min-qnames", type=int, default=2)
    parser.add_argument("--large-insertion-details", default="")
    parser.add_argument("--candidate-read-output", default="")
    parser.add_argument("--output", required=True)
    parser.add_argument("--details", required=True)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.premature_max_length_bp < 1:
        fail("premature-circularization maximum length must be at least 1 bp")
    if args.premature_min_ribosomal_families < 1:
        fail("minimum ribosomal-protein family count must be at least 1")
    if args.min_linked_reads_per_side < 1:
        fail("minimum linked reads per side must be at least 1")
    if args.large_insertion_min_bp < 1:
        fail("large-insertion minimum length must be at least 1 bp")
    if args.large_insertion_min_qnames < 1:
        fail("large-insertion minimum QNAME count must be at least 1")
    run(args)


if __name__ == "__main__":
    main()
