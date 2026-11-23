#!/usr/bin/env python3
"""Taxonomic chimera calls: transition candidates, HiFi validation, MAG summary and clipping-based final calls.

Subcommands:
  evaluate           Find gene-taxonomy transitions along each contig
  validate           Validate transitions with HiFi alignments and SCG redundancy
  mag                Summarize transitions per MAG
  reclassify         Final calls: taxonomic transition plus validated clipping
"""
from __future__ import annotations

import argparse
import bisect
import csv
import math
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, Sequence, Set, Tuple


CALL_RANKS = ("domain", "phylum", "class", "order")

ALL_RANKS = CALL_RANKS

RANK_INDEX = {rank: i for i, rank in enumerate(ALL_RANKS)}

RANK_ALIASES = {
    "domain": ("domain", "superkingdom"),
    "phylum": ("phylum",),
    "class": ("class", "tax_class"),
    "order": ("order", "tax_order"),
}

UNCLASSIFIED_EXACT = {
    "", "0", "1", "na", "n/a", "none", "null", "unknown", "unclassified",
    "no rank", "root", "cellular organisms",
}

AMBIGUOUS_PATTERNS = (
    "uncultured", "environmental sample", "environmental samples", "metagenome",
    "unidentified", "incertae sedis",
)


@dataclass(frozen=True, slots=True)
class Gene:
    gene_id: str
    start: int
    end: int
    taxa: tuple[str | None, ...]

    def taxon(self, rank: str) -> str | None:
        return self.taxa[RANK_INDEX[rank]]


@dataclass(frozen=True, slots=True)
class SideBlock:
    genes: tuple[Gene, ...]
    searched_genes: tuple[Gene, ...]
    dominant_taxon: str
    dominant_count: int
    dominance: float

    @property
    def start(self) -> int:
        return min(gene.start for gene in self.genes)

    @property
    def end(self) -> int:
        return max(gene.end for gene in self.genes)

    @property
    def gene_ids(self) -> str:
        return ",".join(gene.gene_id for gene in self.genes)

    @property
    def searched_gene_ids(self) -> str:
        return ",".join(gene.gene_id for gene in self.searched_genes)

    def classified_count(self, rank: str) -> int:
        return len(self.genes)

    def unclassified_count(self, rank: str) -> int:
        return sum(gene.taxon(rank) is None for gene in self.searched_genes)


@dataclass(frozen=True, slots=True)
class GeneIndex:
    by_start: tuple[Gene, ...]
    starts: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class PositionEvaluation:
    position: int
    transition_rank: str | None
    deepest_evaluable_rank: str | None
    left: SideBlock | None
    right: SideBlock | None
    ambiguous: bool
    reason: str

    @property
    def is_shift(self) -> bool:
        return self.transition_rank is not None

    @property
    def shift_rank(self) -> str | None:
        return self.transition_rank


@dataclass(frozen=True, slots=True)
class ClipRecord:
    position: int
    coverage: int
    clipping_count: int
    clipping_ratio: float
    validated: bool


@dataclass(slots=True)
class RawShift:
    contig: str
    mag: str
    evaluation: PositionEvaluation
    from_default_clipping: bool
    clip_record: ClipRecord | None


@dataclass(slots=True)
class ShiftCluster:
    contig: str
    mag: str
    rank: str
    left_taxon: str
    right_taxon: str
    members: list[RawShift]


EVENT_FIELDS = [
    "event_id", "taxonomy_source", "taxonomy_release", "MAG_ID", "contig",
    "assessment_class", "taxonomic_transition", "shift_rank",
    "analysis_pos_1based", "transition_interval_start_1based",
    "transition_interval_end_1based", "default_clipping_support",
    "n_supporting_default_clipping_positions", "supporting_default_clipping_positions_1based",
    "max_default_clipping_coverage", "max_default_clipping_count",
    "max_default_clipping_ratio", "validated_clipping_support",
    "n_supporting_validated_clipping_positions",
    "supporting_validated_clipping_positions_1based",
    "CHIMERA_EVIDENCE",
    "MAG_GTDB_taxon_at_shift_rank", "left_matches_MAG_at_shift_rank",
    "right_matches_MAG_at_shift_rank", "MAG_discordant_side_at_shift_rank",
    "left_dominant_taxon", "right_dominant_taxon",
    "left_block_genes", "right_block_genes",
    "left_classified_genes_at_rank", "right_classified_genes_at_rank",
    "left_searched_predicted_genes", "right_searched_predicted_genes",
    "left_unclassified_genes_at_rank", "right_unclassified_genes_at_rank",
    "left_gene_ids", "right_gene_ids",
    "left_searched_gene_ids", "right_searched_gene_ids",
    "left_dominant_count", "right_dominant_count", "left_dominance", "right_dominance",
    "left_block_start_1based", "left_block_end_1based", "left_block_span_bp",
    "right_block_start_1based", "right_block_end_1based", "right_block_span_bp",
    "left_GC_percent", "right_GC_percent", "delta_GC_percentage_points",
    "left_GC_callable_bases", "right_GC_callable_bases", "depth_bam",
    "depth_mapq_filter", "left_mean_depth", "right_mean_depth", "depth_right_over_left",
    "abs_log2_depth_ratio", "maximum_search_genes_per_side",
    "target_classified_genes_per_side", "minimum_classified_genes_per_side",
    "minimum_block_dominance", "unclassified_gene_policy", "final_call_ranks",
]

CONTIG_FIELDS = [
    "taxonomy_source", "taxonomy_release", "MAG_ID", "contig", "length",
    "n_predicted_genes", "n_taxonomy_classified_genes",
    "n_gene_count_eligible_boundaries", "n_taxonomically_evaluable_boundaries",
    "Taxonomic_chimera_assessable", "n_default_clipping_candidates",
    "n_validated_clipping_sites", "n_taxonomic_shift_events",
    "n_taxonomic_chimera_with_clipping_support_events",
    "n_clipping_independent_taxonomic_shift_events",
    "taxonomic_shift_ranks", "taxonomic_shift_positions_1based",
    "taxonomic_chimera_with_clipping_support_ranks",
    "taxonomic_chimera_with_clipping_support_positions_1based",
    "clipping_independent_taxonomic_shift_positions_1based",
    "max_delta_GC_percentage_points", "max_abs_log2_depth_ratio",
    "Taxonomic_shift", "Taxonomic_chimera_with_clipping_support",
    "Clipping_independent_taxonomic_shift",
    "CHIMERA_EVIDENCE",
    "High_confidence_taxonomic_chimera",
    "Clipping_independent_taxonomic_chimera", "GTDB_assessment_status", "Taxonomic_chimera",
]


class FastaIndex:
    def __init__(self, fasta: str, fai: str) -> None:
        self.fasta_path = fasta
        self.handle = open(fasta, "rb")
        self.entries: dict[str, tuple[int, int, int, int]] = {}
        with open(fai, encoding="utf-8") as fh:
            for line_number, raw in enumerate(fh, start=1):
                fields = raw.rstrip("\n\r").split("\t")
                if len(fields) < 5:
                    raise SystemExit(f"ERROR: malformed FAI row {line_number}")
                name = fields[0]
                try:
                    length, offset, line_bases, line_width = map(int, fields[1:5])
                except ValueError as exc:
                    raise SystemExit(f"ERROR: malformed FAI numeric field at row {line_number}") from exc
                if not name or length < 1 or line_bases < 1 or line_width < line_bases:
                    raise SystemExit(f"ERROR: invalid FAI row {line_number}")
                self.entries[name] = (length, offset, line_bases, line_width)

    def close(self) -> None:
        self.handle.close()

    def fetch(self, contig: str, start_1based: int, end_1based: int) -> bytes:
        if contig not in self.entries:
            raise SystemExit(f"ERROR: FASTA contig absent from FAI: {contig}")
        length, offset, line_bases, line_width = self.entries[contig]
        if not 1 <= start_1based <= end_1based <= length:
            raise SystemExit(
                f"ERROR: invalid FASTA interval {contig}:{start_1based}-{end_1based}/{length}"
            )
        start0 = start_1based - 1
        remaining = end_1based - start_1based + 1
        parts: list[bytes] = []
        position0 = start0
        while remaining:
            line_index, column = divmod(position0, line_bases)
            take = min(remaining, line_bases - column)
            byte_offset = offset + line_index * line_width + column
            self.handle.seek(byte_offset)
            chunk = self.handle.read(take)
            if len(chunk) != take:
                raise SystemExit(f"ERROR: truncated FASTA while reading {contig}")
            parts.append(chunk)
            position0 += take
            remaining -= take
        return b"".join(parts).upper()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--fai", required=True)
    parser.add_argument("--depth-manifest", required=True)
    parser.add_argument("--clipping", default=None)
    parser.add_argument("--gene-tax-tsv", required=True)
    parser.add_argument("--mag-map", required=True)
    parser.add_argument("--mag-taxonomy-tsv", default=None)
    parser.add_argument("--events-output", required=True)
    parser.add_argument("--supported-events-output", required=True)
    parser.add_argument("--independent-events-output", required=True)
    parser.add_argument("--contig-output", required=True)
    parser.add_argument("--supported-contig-output", required=True)
    parser.add_argument("--independent-contig-output", required=True)
    parser.add_argument("--max-search-genes", required=True, type=int)
    parser.add_argument("--target-classified-genes", required=True, type=int)
    parser.add_argument("--min-classified-genes", required=True, type=int)
    parser.add_argument("--min-dominance", required=True, type=float)
    parser.add_argument("--taxonomy-source", default="GTDB")
    parser.add_argument("--taxonomy-release", required=True)
    parser.add_argument("--parallel-bam-jobs", type=int, default=3)
    return parser.parse_args()


def require_file(path: str, label: str) -> None:
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise SystemExit(f"ERROR: {label} is missing or empty: {path}")


def atomic_write(path: str, fieldnames: Sequence[str], rows: Iterable[Mapping[str, str]]) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=destination.name + ".tmp.", dir=str(destination.parent), text=True)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="raise")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, destination)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def find_field(fieldnames: Sequence[str], aliases: Sequence[str]) -> str | None:
    normalized = {name.strip().casefold(): name for name in fieldnames}
    for alias in aliases:
        field = normalized.get(alias.casefold())
        if field is not None:
            return field
    return None


def normalize_taxon(value: str) -> str | None:
    value = value.strip()
    low = value.casefold()
    if low in UNCLASSIFIED_EXACT:
        return None
    for prefix in (
        "d__", "p__", "c__", "o__", "d_", "p_", "c_", "o_",
        "superkingdom:", "domain:", "phylum:", "class:", "order:",
    ):
        if low.startswith(prefix):
            value = value[len(prefix):].strip()
            low = value.casefold()
            break
    if low in UNCLASSIFIED_EXACT or low.startswith("uc_"):
        return None
    if value.isdigit() or low.startswith("taxid:"):
        return None
    if any(pattern in low for pattern in AMBIGUOUS_PATTERNS):
        return None
    value = re.sub(r"\s+", " ", value).strip()
    return sys.intern(value) if value else None


def enforce_hierarchy(values: Sequence[str | None]) -> tuple[str | None, ...]:
    out: list[str | None] = []
    resolved = True
    for value in values:
        if not resolved:
            out.append(None)
            continue
        out.append(value)
        if value is None:
            resolved = False
    return tuple(out)


def taxonomic_vote(values: Sequence[str | None]) -> tuple[str | None, int, int, float | None]:
    """Return the dominant taxon and its unweighted gene-count support."""
    resolved = [value for value in values if value is not None]
    if not resolved:
        return None, 0, 0, None
    dominant_taxon, dominant_count = Counter(resolved).most_common(1)[0]
    total = len(resolved)
    return str(dominant_taxon), dominant_count, total, dominant_count / total


def read_fai_lengths(path: str) -> dict[str, int]:
    lengths: dict[str, int] = {}
    with open(path, encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            fields = raw.rstrip("\n\r").split("\t")
            if len(fields) < 2:
                raise SystemExit(f"ERROR: malformed FAI row {line_number}")
            try:
                length = int(fields[1])
            except ValueError as exc:
                raise SystemExit(f"ERROR: invalid FAI length at row {line_number}") from exc
            if not fields[0] or length < 1 or fields[0] in lengths:
                raise SystemExit(f"ERROR: invalid or duplicate FAI sequence at row {line_number}")
            lengths[fields[0]] = length
    return lengths


def read_mag_map(path: str, lengths: Mapping[str, int]) -> tuple[dict[str, str], dict[str, list[str]]]:
    contig_to_mag: dict[str, str] = {}
    mag_to_contigs: dict[str, list[str]] = defaultdict(list)
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        lower = {name.strip().lower(): name for name in fields}
        mag_col = lower.get("mag_id") or lower.get("id")
        contig_col = lower.get("contig") or lower.get("original_contig")
        if not mag_col or not contig_col:
            raise SystemExit("ERROR: MAG map requires MAG_ID and Contig columns")
        for line_number, row in enumerate(reader, start=2):
            mag = (row.get(mag_col) or "").strip()
            contig = (row.get(contig_col) or "").strip()
            if not mag or not contig:
                raise SystemExit(f"ERROR: blank MAG/contig in MAG map line {line_number}")
            if contig not in lengths:
                raise SystemExit(f"ERROR: MAG-map contig absent from assembly: {contig}")
            if contig in contig_to_mag:
                raise SystemExit(f"ERROR: duplicate/multiply assigned MAG contig: {contig}")
            contig_to_mag[contig] = mag
            mag_to_contigs[mag].append(contig)
    if not contig_to_mag:
        raise SystemExit("ERROR: MAG map contains no rows")
    return contig_to_mag, mag_to_contigs


def read_clipping(path: str, lengths: Mapping[str, int], target_contigs: set[str]) -> dict[str, list[ClipRecord]]:
    """Read the validated clipping table while retaining default-candidate audit context."""
    result: dict[str, list[ClipRecord]] = defaultdict(list)
    seen: set[tuple[str, int]] = set()
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        required = {
            "contig", "pos_1based", "coverage", "clipping_count",
            "clipping_ratio", "default_candidate", "Clipping",
        }
        missing = required - fields
        if missing:
            raise SystemExit(f"ERROR: validated clipping table lacks columns: {sorted(missing)}")
        for line_number, row in enumerate(reader, start=2):
            contig = (row.get("contig") or "").strip()
            if contig not in lengths:
                raise SystemExit(f"ERROR: clipping contig absent from assembly: {contig}")
            try:
                position = int(row["pos_1based"])
                coverage = int(row["coverage"])
                clipping_count = int(row["clipping_count"])
                ratio = float(row["clipping_ratio"])
            except (TypeError, ValueError) as exc:
                raise SystemExit(f"ERROR: invalid clipping row {line_number}") from exc
            if not 1 <= position <= lengths[contig] or (contig, position) in seen:
                raise SystemExit(f"ERROR: invalid or duplicate clipping site: {contig}:{position}")
            seen.add((contig, position))
            default_flag = (row.get("default_candidate") or "").strip().upper()
            clipping_flag = (row.get("Clipping") or "").strip().upper()
            if default_flag not in {"Y", "N"} or clipping_flag not in {"Y", "N"}:
                raise SystemExit(f"ERROR: invalid clipping flags at line {line_number}")
            if default_flag != "Y":
                raise SystemExit(
                    f"ERROR: validated clipping table contains a non-default candidate at line {line_number}"
                )
            if clipping_flag == "Y" and default_flag != "Y":
                raise SystemExit(f"ERROR: Clipping=Y without default_candidate=Y at line {line_number}")
            if coverage < 1 or clipping_count < 1 or clipping_count > coverage or not 0 < ratio <= 1:
                raise SystemExit(f"ERROR: invalid clipping support at {contig}:{position}")
            if contig in target_contigs:
                result[contig].append(
                    ClipRecord(position, coverage, clipping_count, ratio, clipping_flag == "Y")
                )
    for contig in result:
        result[contig].sort(key=lambda item: item.position)
    return result


def read_mag_taxonomy(path: str, required_mags: set[str]) -> dict[str, tuple[str | None, ...]]:
    result: dict[str, tuple[str | None, ...]] = {}
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"ID", "GTDB_Domain", "GTDB_Phylum", "GTDB_Class", "GTDB_Order"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"ERROR: MAG GTDB-Tk table lacks columns: {sorted(missing)}")
        for line_number, row in enumerate(reader, start=2):
            mag = (row.get("ID") or "").strip()
            if not mag or mag in result:
                raise SystemExit(f"ERROR: blank/duplicate MAG ID in GTDB-Tk table at line {line_number}: {mag}")
            raw = tuple(normalize_taxon(row.get(col) or "") for col in (
                "GTDB_Domain", "GTDB_Phylum", "GTDB_Class", "GTDB_Order"
            ))
            result[mag] = enforce_hierarchy(raw)
    missing_mags = sorted(required_mags - set(result))
    if missing_mags:
        raise SystemExit(
            f"ERROR: {len(missing_mags)} MAG(s) are absent from GTDB-Tk table; first: "
            + ", ".join(missing_mags[:10])
        )
    return result


def simple_chimera_evidence(validated_hifi_clipping: bool) -> str:
    """Return one conservative candidate-stage evidence label for an event."""
    return "HIFI_SUPPORTED" if validated_hifi_clipping else "TAXONOMIC_SHIFT_ONLY"


def aggregate_chimera_evidence(event_labels: Sequence[str], taxonomic_shift: str) -> str:
    """Collapse candidate-stage event labels to one contig-level label."""
    if taxonomic_shift == "NA":
        return "NOT_ASSESSABLE"
    if taxonomic_shift == "N":
        return "NO_TAXONOMIC_SHIFT"
    labels = set(event_labels)
    if "HIFI_SUPPORTED" in labels:
        return "HIFI_SUPPORTED"
    if "TAXONOMIC_SHIFT_ONLY" in labels:
        return "TAXONOMIC_SHIFT_ONLY"
    raise AssertionError("positive taxonomic shift lacks an evidence label")


def gene_groups(path: str, lengths: Mapping[str, int], target_contigs: set[str]) -> Iterator[tuple[str, list[Gene]]]:
    """Yield one contig at a time; current taxonomy output is contig-block sorted."""
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        contig_col = find_field(fields, ("contig", "contig_id", "seqid", "sequence"))
        start_col = find_field(fields, ("start", "start_1based", "gene_start"))
        end_col = find_field(fields, ("end", "end_1based", "gene_end"))
        gene_col = find_field(fields, ("gene_id", "gene", "query", "protein_id"))
        rank_columns = {rank: find_field(fields, RANK_ALIASES[rank]) for rank in ALL_RANKS}
        missing_ranks = [rank for rank, col in rank_columns.items() if col is None]
        if contig_col is None or start_col is None or end_col is None or gene_col is None:
            raise SystemExit("ERROR: GTDB gene-taxonomy table requires contig,start,end,gene_id")
        if missing_ranks:
            raise SystemExit("ERROR: GTDB gene-taxonomy table lacks ranks: " + ",".join(missing_ranks))

        current_contig: str | None = None
        genes: list[Gene] = []
        closed: set[str] = set()
        seen_gene_ids: set[str] = set()
        last_start = -1

        def flush() -> tuple[str, list[Gene]] | None:
            nonlocal current_contig, genes, seen_gene_ids, last_start
            if current_contig is None:
                return None
            out = (current_contig, genes)
            closed.add(current_contig)
            current_contig = None
            genes = []
            seen_gene_ids = set()
            last_start = -1
            return out

        for line_number, row in enumerate(reader, start=2):
            contig = (row.get(contig_col) or "").strip()
            if not contig:
                raise SystemExit(f"ERROR: blank contig at taxonomy line {line_number}")
            if contig not in lengths:
                raise SystemExit(f"ERROR: taxonomy contig absent from assembly: {contig}")
            if current_contig is None:
                if contig in closed:
                    raise SystemExit(f"ERROR: taxonomy table is not contig-block sorted: {contig}")
                current_contig = contig
            elif contig != current_contig:
                emitted = flush()
                if emitted is not None and emitted[0] in target_contigs:
                    yield emitted
                if contig in closed:
                    raise SystemExit(f"ERROR: taxonomy table is not contig-block sorted: {contig}")
                current_contig = contig

            if contig not in target_contigs:
                continue
            try:
                start = int(float(row[start_col]))
                end = int(float(row[end_col]))
            except (TypeError, ValueError, OverflowError) as exc:
                raise SystemExit(f"ERROR: invalid gene coordinates at taxonomy line {line_number}") from exc
            if end < start:
                start, end = end, start
            if start < 1 or end > lengths[contig]:
                raise SystemExit(f"ERROR: gene outside contig at line {line_number}: {contig}:{start}-{end}")
            if start < last_start:
                raise SystemExit(f"ERROR: genes are not start-sorted within contig {contig}")
            last_start = start
            gene_id = (row.get(gene_col) or "").strip()
            if not gene_id or gene_id in seen_gene_ids:
                raise SystemExit(f"ERROR: blank/duplicate gene_id within {contig} at line {line_number}: {gene_id}")
            seen_gene_ids.add(gene_id)
            raw = tuple(normalize_taxon(row.get(rank_columns[rank]) or "") for rank in ALL_RANKS)  # type: ignore[index]
            taxa = enforce_hierarchy(raw)
            genes.append(Gene(gene_id, start, end, taxa))

        emitted = flush()
        if emitted is not None and emitted[0] in target_contigs:
            yield emitted


def build_gene_index(genes: Sequence[Gene]) -> GeneIndex:
    by_start = tuple(sorted(genes, key=lambda gene: (gene.start, gene.end, gene.gene_id)))
    return GeneIndex(
        by_start=by_start,
        starts=tuple(gene.start for gene in by_start),
    )


def build_side_block(
    index: GeneIndex,
    position: int,
    max_search_genes: int,
    target_classified_genes: int,
    min_classified_genes: int,
    left: bool,
    rank: str,
) -> SideBlock | None:
    """Build a bounded adaptive taxonomic block on one side of a boundary."""
    split = bisect.bisect_left(index.starts, position)
    if left:
        candidates = reversed(index.by_start[max(0, split - max_search_genes):split])
    else:
        candidates = iter(index.by_start[split:split + max_search_genes])

    searched_near_to_far: list[Gene] = []
    selected_near_to_far: list[Gene] = []
    for gene in candidates:
        searched_near_to_far.append(gene)
        if gene.taxon(rank) is not None:
            selected_near_to_far.append(gene)
            if len(selected_near_to_far) >= target_classified_genes:
                break

    if len(selected_near_to_far) < min_classified_genes:
        return None

    selected = tuple(sorted(selected_near_to_far, key=lambda gene: (gene.start, gene.end, gene.gene_id)))
    searched = tuple(sorted(searched_near_to_far, key=lambda gene: (gene.start, gene.end, gene.gene_id)))
    taxa = [gene.taxon(rank) for gene in selected]
    if any(value is None for value in taxa):
        raise AssertionError("selected taxonomic block contains an unclassified gene")
    dominant_taxon, dominant_count, n_resolved, dominance = taxonomic_vote(taxa)
    if dominant_taxon is None or dominance is None or n_resolved != len(selected):
        raise AssertionError("failed to summarize a fully classified local taxonomic block")
    return SideBlock(selected, searched, dominant_taxon, dominant_count, dominance)


def qualifying_transition(left: SideBlock | None, right: SideBlock | None, min_dominance: float) -> bool:
    return bool(
        left is not None and right is not None
        and left.dominance >= min_dominance and right.dominance >= min_dominance
        and left.dominant_taxon != right.dominant_taxon
    )


def evaluate_position(
    index: GeneIndex,
    position: int,
    max_search_genes: int,
    target_classified_genes: int,
    min_classified_genes: int,
    min_dominance: float,
) -> PositionEvaluation:
    blocks: dict[str, tuple[SideBlock, SideBlock]] = {}
    for rank in ALL_RANKS:
        left = build_side_block(
            index, position, max_search_genes, target_classified_genes, min_classified_genes, True, rank
        )
        right = build_side_block(
            index, position, max_search_genes, target_classified_genes, min_classified_genes, False, rank
        )
        if left is not None and right is not None:
            blocks[rank] = (left, right)
    evaluable = [rank for rank in CALL_RANKS if rank in blocks]
    if not evaluable:
        return PositionEvaluation(
            position, None, None, None, None, False,
            "insufficient_classified_genes_within_bounded_search",
        )
    deepest = evaluable[-1]
    for rank in CALL_RANKS:
        pair = blocks.get(rank)
        if pair and qualifying_transition(pair[0], pair[1], min_dominance):
            return PositionEvaluation(position, rank, deepest, pair[0], pair[1], False, "NA")
    left, right = blocks[deepest]
    ambiguous = left.dominance < min_dominance or right.dominance < min_dominance
    reason = "low_dominance_among_selected_classified_genes" if ambiguous else "NA"
    return PositionEvaluation(position, None, deepest, left, right, ambiguous, reason)


def candidate_boundaries(index: GeneIndex) -> set[int]:
    return set(index.starts[1:])


def _block_gene_ids(item: RawShift) -> set[str]:
    ev = item.evaluation
    ids: set[str] = set()
    if ev.left is not None:
        ids.update(g.gene_id for g in ev.left.genes)
    if ev.right is not None:
        ids.update(g.gene_id for g in ev.right.genes)
    return ids


def cluster_raw_shifts(raw: Sequence[RawShift]) -> list[ShiftCluster]:
    grouped: dict[tuple[str, str, str, str], list[RawShift]] = defaultdict(list)
    for item in raw:
        ev = item.evaluation
        if not ev.is_shift or ev.left is None or ev.right is None or ev.shift_rank is None:
            continue
        grouped[(item.contig, ev.shift_rank, ev.left.dominant_taxon, ev.right.dominant_taxon)].append(item)

    clusters: list[ShiftCluster] = []
    for (_contig, rank, left_taxon, right_taxon), items in grouped.items():
        items.sort(key=lambda item: item.evaluation.position)
        current: list[RawShift] = []
        current_gene_ids: set[str] = set()
        for item in items:
            item_ids = _block_gene_ids(item)
            if current and current_gene_ids.isdisjoint(item_ids):
                clusters.append(ShiftCluster(current[0].contig, current[0].mag, rank, left_taxon, right_taxon, current))
                current = []
                current_gene_ids = set()
            current.append(item)
            current_gene_ids.update(item_ids)
        if current:
            clusters.append(ShiftCluster(current[0].contig, current[0].mag, rank, left_taxon, right_taxon, current))
    clusters.sort(key=lambda c: (c.mag, c.contig, min(m.evaluation.position for m in c.members), RANK_INDEX[c.rank]))
    return clusters


def choose_representative(cluster: ShiftCluster) -> RawShift:
    clip_members = [member for member in cluster.members if member.from_default_clipping and member.clip_record is not None]
    validated_members = [member for member in clip_members if member.clip_record and member.clip_record.validated]
    preferred = validated_members or clip_members
    if preferred:
        return max(
            preferred,
            key=lambda member: (
                member.clip_record.coverage if member.clip_record else -1,
                member.clip_record.clipping_ratio if member.clip_record else -1.0,
                member.clip_record.clipping_count if member.clip_record else -1,
                -member.evaluation.position,
            ),
        )
    return max(
        cluster.members,
        key=lambda member: (
            min(member.evaluation.left.dominance, member.evaluation.right.dominance)  # type: ignore[union-attr]
            if member.evaluation.left and member.evaluation.right else -1.0,
            -member.evaluation.position,
        ),
    )


def gc_metrics(fasta: FastaIndex, contig: str, start: int, end: int) -> tuple[float | None, int]:
    seq = fasta.fetch(contig, start, end)
    a = seq.count(b"A")
    c = seq.count(b"C")
    g = seq.count(b"G")
    t = seq.count(b"T")
    callable_bases = a + c + g + t
    if callable_bases == 0:
        return None, 0
    return 100.0 * (g + c) / callable_bases, callable_bases


def read_depth_routes(manifest_path: str, target_contigs: set[str]) -> dict[str, tuple[str, str]]:
    """Route each target contig to its assembler-specific primary HiFi BAM."""
    routes: dict[str, tuple[str, str]] = {}
    with open(manifest_path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Assembler", "Assembly_FASTA", "HiFi_primary_BAM"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"ERROR: depth manifest lacks columns: {sorted(missing)}")
        rows = list(reader)
    if not rows:
        raise SystemExit(f"ERROR: depth manifest contains no assembler rows: {manifest_path}")
    for raw in rows:
        assembler = (raw.get("Assembler") or "").strip()
        assembly = (raw.get("Assembly_FASTA") or "").strip()
        bam = (raw.get("HiFi_primary_BAM") or "").strip()
        if not assembler or not assembly or not bam or bam == "NA":
            raise SystemExit("ERROR: blank assembler-specific depth manifest field")
        require_file(assembly, f"{assembler} assembly FASTA")
        fai = assembly + ".fai"
        require_file(fai, f"{assembler} assembly FAI")
        require_file(bam, f"{assembler} primary-only HiFi BAM")
        with open(fai, encoding="utf-8") as fai_handle:
            for line_number, raw_fai in enumerate(fai_handle, start=1):
                fields = raw_fai.rstrip("\n").split("\t")
                if len(fields) < 2:
                    raise SystemExit(f"ERROR: malformed FAI line {line_number}: {fai}")
                contig = fields[0]
                if contig not in target_contigs:
                    continue
                if contig in routes:
                    prev_assembler, _ = routes[contig]
                    raise SystemExit(
                        f"ERROR: target contig occurs in multiple assembler references: "
                        f"{contig} ({prev_assembler}, {assembler})"
                    )
                routes[contig] = (assembler, bam)
    missing_contigs = sorted(target_contigs - set(routes))
    if missing_contigs:
        preview = ", ".join(missing_contigs[:10])
        raise SystemExit(
            f"ERROR: {len(missing_contigs)} target contig(s) lack a unique assembler-specific "
            f"primary HiFi BAM route; first: {preview}"
        )
    return routes


def bedcov_means(
    depth_routes: Mapping[str, tuple[str, str]],
    intervals: Sequence[tuple[str, int, int, str]],
    parallel_bam_jobs: int = 3,
) -> dict[str, float]:
    if not intervals:
        return {}
    grouped: dict[str, list[tuple[str, int, int, str]]] = defaultdict(list)
    for interval in intervals:
        contig = interval[0]
        try:
            _assembler, bam = depth_routes[contig]
        except KeyError as exc:
            raise SystemExit(f"ERROR: no assembler-specific depth route for contig: {contig}") from exc
        grouped[bam].append(interval)

    if parallel_bam_jobs < 1:
        raise SystemExit("ERROR: --parallel-bam-jobs must be >=1")

    def run_group(group_number: int, depth_bam: str, bam_intervals, tmpdir: str):
        bed = Path(tmpdir) / f"windows.{group_number}.bed"
        with bed.open("w", encoding="utf-8", newline="") as handle:
            for contig, start0, end0, name in bam_intervals:
                handle.write(f"{contig}\t{start0}\t{end0}\t{name}\n")
        command = ["samtools", "bedcov", str(bed), depth_bam]
        local: dict[str, float] = {}
        with tempfile.TemporaryFile(mode="w+t", encoding="utf-8") as error_output:
            with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=error_output, text=True) as proc:
                assert proc.stdout is not None
                for line_number, raw in enumerate(proc.stdout, start=1):
                    fields = raw.rstrip("\n\r").split("\t")
                    if len(fields) < 5:
                        raise SystemExit(f"ERROR: malformed samtools bedcov output row {line_number}")
                    try:
                        start0, end0 = int(fields[1]), int(fields[2])
                        coverage_sum = float(fields[-1])
                    except ValueError as exc:
                        raise SystemExit(f"ERROR: malformed samtools bedcov numeric output row {line_number}") from exc
                    length = end0 - start0
                    if length <= 0:
                        raise SystemExit(f"ERROR: non-positive context window in bedcov row {line_number}")
                    local[fields[3]] = coverage_sum / length
            if proc.returncode:
                error_output.seek(0)
                raise subprocess.CalledProcessError(proc.returncode, command, stderr=error_output.read())
        return local

    result: dict[str, float] = {}
    items = list(sorted(grouped.items()))
    jobs = min(parallel_bam_jobs, len(items))
    with tempfile.TemporaryDirectory(prefix="chimera_bedcov_") as tmpdir:
        if jobs <= 1:
            for group_number, (depth_bam, bam_intervals) in enumerate(items, start=1):
                result.update(run_group(group_number, depth_bam, bam_intervals, tmpdir))
        else:
            with ThreadPoolExecutor(max_workers=jobs) as executor:
                futures = {
                    executor.submit(run_group, n, bam, ivals, tmpdir): bam
                    for n, (bam, ivals) in enumerate(items, start=1)
                }
                for future in as_completed(futures):
                    result.update(future.result())
    if len(result) != len(intervals):
        raise SystemExit("ERROR: samtools bedcov did not return every requested context window")
    return result


def fmt_float(value: float | None, digits: int = 6) -> str:
    if value is None or not math.isfinite(value):
        return "NA"
    return f"{value:.{digits}f}"


def _evaluate_main() -> None:
    args = parse_args()
    if args.taxonomy_source.strip().upper() != "GTDB":
        raise SystemExit("ERROR: this evaluator accepts GTDB taxonomy only")
    if args.max_search_genes < 1:
        raise SystemExit("ERROR: --max-search-genes must be positive")
    if args.target_classified_genes < 1 or args.target_classified_genes > args.max_search_genes:
        raise SystemExit("ERROR: --target-classified-genes must be in [1, --max-search-genes]")
    if args.min_classified_genes < 1 or args.min_classified_genes > args.target_classified_genes:
        raise SystemExit("ERROR: --min-classified-genes must be in [1, --target-classified-genes]")
    if not 0.5 <= args.min_dominance <= 1.0:
        raise SystemExit("ERROR: --min-dominance must be between 0.5 and 1.0")
    if args.parallel_bam_jobs < 1:
        raise SystemExit("ERROR: --parallel-bam-jobs must be >=1")
    for path, label in (
        (args.assembly, "assembly FASTA"), (args.fai, "assembly FAI"),
        (args.depth_manifest, "assembler mapping/depth manifest"),
        (args.gene_tax_tsv, "GTDB gene-taxonomy table"), (args.mag_map, "MAG-contig map"),
    ):
        require_file(path, label)
    if args.clipping is not None:
        require_file(args.clipping, "validated clipping table")
    if not shutil_which("samtools"):
        raise SystemExit("ERROR: samtools is not available in PATH")

    lengths = read_fai_lengths(args.fai)
    contig_to_mag, mag_to_contigs = read_mag_map(args.mag_map, lengths)
    target_contigs = set(lengths)
    mag_taxonomy = read_mag_taxonomy(args.mag_taxonomy_tsv, set(mag_to_contigs)) if args.mag_taxonomy_tsv else {}
    depth_routes = read_depth_routes(args.depth_manifest, target_contigs)
    clipping = read_clipping(args.clipping, lengths, target_contigs) if args.clipping else defaultdict(list)
    clip_lookup = {contig: {item.position: item for item in records} for contig, records in clipping.items()}

    raw_shifts: list[RawShift] = []
    predicted_gene_counts: dict[str, int] = {contig: 0 for contig in target_contigs}
    classified_counts: dict[str, int] = {contig: 0 for contig in target_contigs}
    eligible_boundary_counts: dict[str, int] = {contig: 0 for contig in target_contigs}
    evaluable_boundary_counts: dict[str, int] = {contig: 0 for contig in target_contigs}
    for contig, genes in gene_groups(args.gene_tax_tsv, lengths, target_contigs):
        predicted_gene_counts[contig] = len(genes)
        classified_counts[contig] = sum(gene.taxon("domain") is not None for gene in genes)

        if len(genes) < 2 * args.min_classified_genes:
            continue
        if classified_counts[contig] < 2 * args.min_classified_genes:
            continue

        index = build_gene_index(genes)

        evaluations: dict[int, PositionEvaluation] = {}
        for position in sorted(candidate_boundaries(index)):
            split = bisect.bisect_left(index.starts, position)
            if split >= args.min_classified_genes and len(index.by_start) - split >= args.min_classified_genes:
                eligible_boundary_counts[contig] += 1
            evaluation = evaluate_position(
                index, position, args.max_search_genes, args.target_classified_genes,
                args.min_classified_genes, args.min_dominance
            )
            evaluations[position] = evaluation
            if evaluation.deepest_evaluable_rank is not None:
                evaluable_boundary_counts[contig] += 1

        for record in clipping.get(contig, []):
            if record.position not in evaluations:
                evaluations[record.position] = evaluate_position(
                    index, record.position, args.max_search_genes, args.target_classified_genes,
                    args.min_classified_genes, args.min_dominance
                )

        for position, evaluation in sorted(evaluations.items()):
            if position < 1 or position > lengths[contig] or not evaluation.is_shift:
                continue
            clip_record = clip_lookup.get(contig, {}).get(position)
            raw_shifts.append(
                RawShift(
                    contig=contig,
                    mag=contig_to_mag.get(contig, "NA"),
                    evaluation=evaluation,
                    from_default_clipping=clip_record is not None,
                    clip_record=clip_record,
                )
            )

    clusters = cluster_raw_shifts(raw_shifts)
    fasta = FastaIndex(args.assembly, args.fai)
    event_rows: list[dict[str, str]] = []
    bed_intervals: list[tuple[str, int, int, str]] = []
    event_meta: list[tuple[str, str]] = []

    try:
        for event_number, cluster in enumerate(clusters, start=1):
            representative = choose_representative(cluster)
            ev = representative.evaluation
            if ev.left is None or ev.right is None:
                raise AssertionError("shift representative lacks taxonomic side blocks")
            supporting = sorted(
                (member.clip_record for member in cluster.members if member.from_default_clipping and member.clip_record is not None),
                key=lambda record: record.position,
            )
            validated_supporting = [record for record in supporting if record.validated]
            evidence_summary = simple_chimera_evidence(bool(validated_supporting))
            expected_taxon = mag_taxonomy.get(cluster.mag, (None,) * len(ALL_RANKS))[RANK_INDEX[cluster.rank]]
            left_matches_mag = "NA" if expected_taxon is None else ("Y" if ev.left.dominant_taxon == expected_taxon else "N")
            right_matches_mag = "NA" if expected_taxon is None else ("Y" if ev.right.dominant_taxon == expected_taxon else "N")
            if expected_taxon is None:
                mag_discordant_side = "NA"
            elif left_matches_mag == "Y" and right_matches_mag == "N":
                mag_discordant_side = "RIGHT"
            elif left_matches_mag == "N" and right_matches_mag == "Y":
                mag_discordant_side = "LEFT"
            elif left_matches_mag == "N" and right_matches_mag == "N":
                mag_discordant_side = "BOTH"
            else:
                mag_discordant_side = "NA"
            if validated_supporting:
                assessment_class = "HIGH_CONFIDENCE_TAXONOMIC_CHIMERA"
            else:
                assessment_class = "CLIPPING_INDEPENDENT_TAXONOMIC_CHIMERA"

            analysis_pos = ev.position
            left_start, left_end = ev.left.start, ev.left.end
            right_start, right_end = ev.right.start, ev.right.end
            left_gc, left_callable = gc_metrics(fasta, cluster.contig, left_start, left_end)
            right_gc, right_callable = gc_metrics(fasta, cluster.contig, right_start, right_end)
            delta_gc = abs(left_gc - right_gc) if left_gc is not None and right_gc is not None else None

            event_id = f"gtdb_shift_{event_number:08d}"
            left_window_id = event_id + ":L"
            right_window_id = event_id + ":R"
            bed_intervals.append((cluster.contig, left_start - 1, left_end, left_window_id))
            bed_intervals.append((cluster.contig, right_start - 1, right_end, right_window_id))
            event_meta.append((left_window_id, right_window_id))

            interval_start = min(member.evaluation.position for member in cluster.members)
            interval_end = max(member.evaluation.position for member in cluster.members)
            max_clip = max(supporting, key=lambda record: (record.coverage, record.clipping_ratio, record.clipping_count), default=None)
            event_rows.append({
                "event_id": event_id,
                "taxonomy_source": "GTDB",
                "taxonomy_release": args.taxonomy_release,
                "MAG_ID": cluster.mag,
                "contig": cluster.contig,
                "assessment_class": assessment_class,
                "taxonomic_transition": "Y",
                "shift_rank": cluster.rank,
                "analysis_pos_1based": str(analysis_pos),
                "transition_interval_start_1based": str(interval_start),
                "transition_interval_end_1based": str(interval_end),
                "default_clipping_support": "Y" if supporting else "N",
                "n_supporting_default_clipping_positions": str(len(supporting)),
                "supporting_default_clipping_positions_1based": ";".join(str(record.position) for record in supporting) or "NA",
                "max_default_clipping_coverage": str(max_clip.coverage) if max_clip else "NA",
                "max_default_clipping_count": str(max_clip.clipping_count) if max_clip else "NA",
                "max_default_clipping_ratio": f"{max_clip.clipping_ratio:.6f}" if max_clip else "NA",
                "validated_clipping_support": "Y" if validated_supporting else "N",
                "n_supporting_validated_clipping_positions": str(len(validated_supporting)),
                "supporting_validated_clipping_positions_1based": ";".join(
                    str(record.position) for record in validated_supporting
                ) or "NA",
                "CHIMERA_EVIDENCE": evidence_summary,
                "MAG_GTDB_taxon_at_shift_rank": expected_taxon or "NA",
                "left_matches_MAG_at_shift_rank": left_matches_mag,
                "right_matches_MAG_at_shift_rank": right_matches_mag,
                "MAG_discordant_side_at_shift_rank": mag_discordant_side,
                "left_dominant_taxon": ev.left.dominant_taxon,
                "right_dominant_taxon": ev.right.dominant_taxon,
                "left_block_genes": str(len(ev.left.genes)),
                "right_block_genes": str(len(ev.right.genes)),
                "left_classified_genes_at_rank": str(ev.left.classified_count(cluster.rank)),
                "right_classified_genes_at_rank": str(ev.right.classified_count(cluster.rank)),
                "left_searched_predicted_genes": str(len(ev.left.searched_genes)),
                "right_searched_predicted_genes": str(len(ev.right.searched_genes)),
                "left_unclassified_genes_at_rank": str(ev.left.unclassified_count(cluster.rank)),
                "right_unclassified_genes_at_rank": str(ev.right.unclassified_count(cluster.rank)),
                "left_gene_ids": ev.left.gene_ids,
                "right_gene_ids": ev.right.gene_ids,
                "left_searched_gene_ids": ev.left.searched_gene_ids,
                "right_searched_gene_ids": ev.right.searched_gene_ids,
                "left_dominant_count": str(ev.left.dominant_count),
                "right_dominant_count": str(ev.right.dominant_count),
                "left_dominance": f"{ev.left.dominance:.4f}",
                "right_dominance": f"{ev.right.dominance:.4f}",
                "left_block_start_1based": str(ev.left.start),
                "left_block_end_1based": str(ev.left.end),
                "left_block_span_bp": str(ev.left.end - ev.left.start + 1),
                "right_block_start_1based": str(ev.right.start),
                "right_block_end_1based": str(ev.right.end),
                "right_block_span_bp": str(ev.right.end - ev.right.start + 1),
                "left_GC_percent": fmt_float(left_gc),
                "right_GC_percent": fmt_float(right_gc),
                "delta_GC_percentage_points": fmt_float(delta_gc),
                "left_GC_callable_bases": str(left_callable),
                "right_GC_callable_bases": str(right_callable),
                "depth_bam": depth_routes[cluster.contig][1],
                "depth_mapq_filter": "none",
                "left_mean_depth": "PENDING",
                "right_mean_depth": "PENDING",
                "depth_right_over_left": "PENDING",
                "abs_log2_depth_ratio": "PENDING",
                "maximum_search_genes_per_side": str(args.max_search_genes),
                "target_classified_genes_per_side": str(args.target_classified_genes),
                "minimum_classified_genes_per_side": str(args.min_classified_genes),
                "minimum_block_dominance": f"{args.min_dominance:.4f}",
                "unclassified_gene_policy": "skipped_for_voting_within_bounded_predicted_gene_search",
                "final_call_ranks": ",".join(CALL_RANKS),
            })
    finally:
        fasta.close()

    depth = bedcov_means(depth_routes, bed_intervals, args.parallel_bam_jobs)
    for row, (left_id, right_id) in zip(event_rows, event_meta):
        left_depth = depth[left_id]
        right_depth = depth[right_id]
        row["left_mean_depth"] = f"{left_depth:.6f}"
        row["right_mean_depth"] = f"{right_depth:.6f}"
        if left_depth > 0:
            row["depth_right_over_left"] = f"{right_depth / left_depth:.6f}"
        else:
            row["depth_right_over_left"] = "NA"
        if left_depth > 0 and right_depth > 0:
            row["abs_log2_depth_ratio"] = f"{abs(math.log2(right_depth / left_depth)):.6f}"
        else:
            row["abs_log2_depth_ratio"] = "NA"

    event_rows.sort(key=lambda row: (row["MAG_ID"], row["contig"], int(row["analysis_pos_1based"]), RANK_INDEX[row["shift_rank"]]))
    supported_events = [row for row in event_rows if row["assessment_class"] == "HIGH_CONFIDENCE_TAXONOMIC_CHIMERA"]
    independent_events = [row for row in event_rows if row["assessment_class"] == "CLIPPING_INDEPENDENT_TAXONOMIC_CHIMERA"]

    by_contig: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in event_rows:
        by_contig[row["contig"]].append(row)
    contig_rows: list[dict[str, str]] = []
    for contig in sorted(target_contigs, key=lambda c: (contig_to_mag.get(c, "NA"), c)):
        rows = by_contig.get(contig, [])
        call_rows = [row for row in rows if row["taxonomic_transition"] == "Y"]
        supported_rows = [row for row in rows if row["assessment_class"] == "HIGH_CONFIDENCE_TAXONOMIC_CHIMERA"]
        independent_rows = [row for row in rows if row["assessment_class"] == "CLIPPING_INDEPENDENT_TAXONOMIC_CHIMERA"]

        gc_values = [float(row["delta_GC_percentage_points"]) for row in call_rows if row["delta_GC_percentage_points"] != "NA"]
        depth_values = [float(row["abs_log2_depth_ratio"]) for row in call_rows if row["abs_log2_depth_ratio"] != "NA"]
        assessable = evaluable_boundary_counts.get(contig, 0) > 0
        if supported_rows and independent_rows:
            status = "HIGH_CONFIDENCE_AND_CLIPPING_INDEPENDENT_TAXONOMIC_CHIMERAS"
        elif supported_rows:
            status = "HIGH_CONFIDENCE_TAXONOMIC_CHIMERA"
        elif independent_rows:
            status = "CLIPPING_INDEPENDENT_TAXONOMIC_CHIMERA"
        elif eligible_boundary_counts.get(contig, 0) == 0:
            status = "NOT_ASSESSABLE_INSUFFICIENT_LOCAL_GENE_COUNT"
        elif not assessable:
            status = "NOT_ASSESSABLE_INSUFFICIENT_TAXONOMIC_RESOLUTION"
        else:
            status = "NO_TAXONOMIC_CHIMERA"
        tax_call = "Y" if call_rows else ("N" if assessable else "NA")
        supported_call = "Y" if supported_rows else ("N" if assessable else "NA")
        independent_call = "Y" if independent_rows else ("N" if assessable else "NA")
        chimera_evidence = aggregate_chimera_evidence(
            [row["CHIMERA_EVIDENCE"] for row in call_rows], tax_call
        )

        contig_rows.append({
            "taxonomy_source": "GTDB",
            "taxonomy_release": args.taxonomy_release,
            "MAG_ID": contig_to_mag.get(contig, "NA"),
            "contig": contig,
            "length": str(lengths[contig]),
            "n_predicted_genes": str(predicted_gene_counts.get(contig, 0)),
            "n_taxonomy_classified_genes": str(classified_counts.get(contig, 0)),
            "n_gene_count_eligible_boundaries": str(eligible_boundary_counts.get(contig, 0)),
            "n_taxonomically_evaluable_boundaries": str(evaluable_boundary_counts.get(contig, 0)),
            "Taxonomic_chimera_assessable": "Y" if assessable else "N",
            "n_default_clipping_candidates": str(len(clipping.get(contig, []))),
            "n_validated_clipping_sites": str(sum(record.validated for record in clipping.get(contig, []))),
            "n_taxonomic_shift_events": str(len(call_rows)),
            "n_taxonomic_chimera_with_clipping_support_events": str(len(supported_rows)),
            "n_clipping_independent_taxonomic_shift_events": str(len(independent_rows)),
            "taxonomic_shift_ranks": ";".join(row["shift_rank"] for row in call_rows) or "NA",
            "taxonomic_shift_positions_1based": ";".join(row["analysis_pos_1based"] for row in call_rows) or "NA",
            "taxonomic_chimera_with_clipping_support_ranks": ";".join(row["shift_rank"] for row in supported_rows) or "NA",
            "taxonomic_chimera_with_clipping_support_positions_1based": ";".join(row["analysis_pos_1based"] for row in supported_rows) or "NA",
            "clipping_independent_taxonomic_shift_positions_1based": ";".join(row["analysis_pos_1based"] for row in independent_rows) or "NA",
            "max_delta_GC_percentage_points": f"{max(gc_values):.6f}" if gc_values else "NA",
            "max_abs_log2_depth_ratio": f"{max(depth_values):.6f}" if depth_values else "NA",
            "Taxonomic_shift": tax_call,
            "Taxonomic_chimera_with_clipping_support": supported_call,
            "Clipping_independent_taxonomic_shift": independent_call,
            "CHIMERA_EVIDENCE": chimera_evidence,
            "High_confidence_taxonomic_chimera": supported_call,
            "Clipping_independent_taxonomic_chimera": independent_call,
            "GTDB_assessment_status": status,
            "Taxonomic_chimera": tax_call,
        })

    supported_contigs = [row for row in contig_rows if row["Taxonomic_chimera_with_clipping_support"] == "Y"]
    independent_contigs = [row for row in contig_rows if row["Clipping_independent_taxonomic_shift"] == "Y"]

    atomic_write(args.events_output, EVENT_FIELDS, event_rows)
    atomic_write(args.supported_events_output, EVENT_FIELDS, supported_events)
    atomic_write(args.independent_events_output, EVENT_FIELDS, independent_events)
    atomic_write(args.contig_output, CONTIG_FIELDS, contig_rows)
    atomic_write(args.supported_contig_output, CONTIG_FIELDS, supported_contigs)
    atomic_write(args.independent_contig_output, CONTIG_FIELDS, independent_contigs)

    print(f"GTDB all-assembly contigs scanned: {len(contig_rows)}")
    print(f"MAG-associated contigs: {sum(row['MAG_ID'] != 'NA' for row in contig_rows)}")
    print(f"Taxonomically assessable contigs: {sum(row['Taxonomic_chimera_assessable'] == 'Y' for row in contig_rows)}")
    print(f"Whole-contig domain-to-order taxonomic-transition candidates: {sum(int(row['n_taxonomic_shift_events']) for row in contig_rows)}")
    print(f"Taxonomic-transition candidates overlapping validated clipping: {len(supported_events)}")
    print(f"Taxonomic-transition candidates without validated HiFi clipping: {len(independent_events)}")
    print("This evaluator is candidate-generation only; final structural-chimera calls require downstream HiFi continuity validation.")


def shutil_which(command: str) -> str | None:
    for directory in os.environ.get("PATH", "").split(os.pathsep):
        candidate = os.path.join(directory or ".", command)
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate
    return None


def cmd_evaluate(argv):
    """Find gene-taxonomy transitions along each contig"""
    sys.argv = [sys.argv[0] + " evaluate", *argv]
    return _evaluate_main()


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")

MATCH_OPS = {"M", "=", "X"}

REF_OPS = {"M", "D", "N", "=", "X"}

CLIP_OPS = {"S", "H"}


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def _validate_parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--candidate-events", required=True)
    p.add_argument("--candidate-contigs", required=True)
    p.add_argument("--mapping-manifest", required=True)
    p.add_argument("--mag-map", required=True)
    p.add_argument("--mag-taxonomy-tsv", default=None)
    p.add_argument("--clipping", default=None)
    p.add_argument("--scg-tsv", default=None)
    p.add_argument("--events-output", required=True)
    p.add_argument("--structural-events-output", required=True)
    p.add_argument("--high-confidence-events-output", required=True)
    p.add_argument("--unresolved-events-output", required=True)
    p.add_argument("--contig-output", required=True)
    p.add_argument("--structural-contig-output", required=True)
    p.add_argument("--high-confidence-contig-output", required=True)
    p.add_argument("--unresolved-contig-output", required=True)
    p.add_argument("--scg-audit-output", required=True)
    p.add_argument("--anchor-bp", type=int, default=500)
    p.add_argument("--min-support-qnames", type=int, default=2)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--parallel-bam-jobs", type=int, default=3)
    return p.parse_args()


def _validate_require_file(path: str | None, label: str, optional: bool = False) -> None:
    if path is None:
        if optional:
            return
        fail(f"{label} was not provided")
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        fail(f"{label} is missing or empty: {path}")


def read_tsv(path: str) -> tuple[list[str], list[dict[str, str]]]:
    with open(path, encoding="utf-8-sig", newline="") as h:
        r = csv.DictReader(h, delimiter="\t")
        fields = list(r.fieldnames or [])
        return fields, [dict(row) for row in r]


def _validate_atomic_write(path: str, fields: Sequence[str], rows: Iterable[Mapping[str, str]]) -> None:
    dst = Path(path)
    dst.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=dst.name + ".tmp.", dir=str(dst.parent), text=True)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as h:
            w = csv.DictWriter(h, fieldnames=list(fields), delimiter="\t", lineterminator="\n", extrasaction="raise")
            w.writeheader(); w.writerows(rows)
        os.replace(tmp, dst)
    except BaseException:
        try: os.unlink(tmp)
        except FileNotFoundError: pass
        raise


def cigar_ops(cigar: str) -> list[tuple[int, str]]:
    ops = [(int(n), op) for n, op in CIGAR_RE.findall(cigar)]
    if not ops or "".join(f"{n}{op}" for n, op in ops) != cigar:
        fail(f"invalid CIGAR: {cigar}")
    return ops


def terminal_clip(ops: Sequence[tuple[int, str]], left: bool) -> int:
    seq = ops if left else reversed(ops)
    total = 0
    for n, op in seq:
        if op in CLIP_OPS: total += n
        else: break
    return total


def matched_runs(ref_start0: int, ops: Sequence[tuple[int, str]]) -> list[tuple[int, int]]:
    runs: list[tuple[int, int]] = []
    ref = ref_start0
    start = end = None
    for n, op in ops:
        if op in MATCH_OPS:
            a, b = ref, ref + n
            if start is None:
                start, end = a, b
            elif end == a:
                end = b
            else:
                runs.append((start, end)); start, end = a, b
            ref = b
        elif op in {"I", "P"}:
            continue
        elif op in {"D", "N"}:
            if start is not None:
                runs.append((start, end)); start = end = None
            ref += n
        elif op in {"S", "H"}:
            if start is not None:
                runs.append((start, end)); start = end = None
        else:
            fail(f"unsupported CIGAR operation: {op}")
    if start is not None:
        runs.append((start, end))
    return runs


def matched_bp(runs: Sequence[tuple[int, int]], start0: int, end0: int) -> int:
    total = 0
    for a, b in runs:
        total += max(0, min(b, end0) - max(a, start0))
    return total


def _validate_read_mag_map(path: str) -> dict[str, str]:
    fields, rows = read_tsv(path)
    lowered = {x.lower(): x for x in fields}
    mag_col = lowered.get("mag_id") or lowered.get("mag") or lowered.get("id")
    contig_col = lowered.get("contig") or lowered.get("contig_id")
    if not mag_col or not contig_col:
        fail("MAG-contig map needs MAG_ID/ID and contig columns")
    out = {}
    for row in rows:
        mag, contig = row[mag_col].strip(), row[contig_col].strip()
        if not mag or not contig or contig in out:
            fail(f"invalid/duplicate MAG-contig map row: {mag}/{contig}")
        out[contig] = mag
    return out


def read_fai_names(path: str) -> set[str]:
    out = set()
    with open(path, encoding="utf-8") as h:
        for raw in h:
            if raw.strip(): out.add(raw.split("\t", 1)[0])
    return out


def load_routes(path: str, target_contigs: set[str]) -> dict[str, dict[str, str]]:
    fields, rows = read_tsv(path)
    required = {"Assembler", "Assembly_FASTA", "HiFi_primary_supplementary_BAM"}
    missing = required - set(fields)
    if missing: fail(f"mapping manifest lacks columns: {sorted(missing)}")
    routes: dict[str, dict[str, str]] = {}
    for row in rows:
        assembler = row["Assembler"].strip(); assembly = row["Assembly_FASTA"].strip()
        ps_bam = row["HiFi_primary_supplementary_BAM"].strip()
        if not assembler or not assembly or not ps_bam:
            fail("blank structural-validation mapping manifest field")
        for p, label in ((assembly + ".fai", f"{assembler} assembly FAI"), (ps_bam, f"{assembler} primary+supplementary HiFi BAM")):
            _validate_require_file(p, label)
        refs = read_fai_names(assembly + ".fai")
        for contig in target_contigs & refs:
            if contig in routes:
                fail(f"contig occurs in multiple assembler references: {contig}")
            routes[contig] = {"assembler": assembler, "primary_supp": ps_bam}
    missing = sorted(target_contigs - set(routes))
    if missing: fail(f"{len(missing)} target contigs lack a unique assembler route; first: {', '.join(missing[:10])}")
    return routes


def load_clipping(path: str | None) -> dict[str, list[dict[str, str]]]:
    result: dict[str, list[dict[str, str]]] = defaultdict(list)
    if not path: return result
    fields, rows = read_tsv(path)
    required = {"contig", "pos_1based", "Clipping", "Stringent_clipping"}
    missing = required - set(fields)
    if missing: fail(f"clipping table lacks columns: {sorted(missing)}")
    for row in rows:
        if row["Clipping"].strip().upper() == "Y":
            result[row["contig"].strip()].append(row)
    return result


def normalize_domain(value: str) -> str | None:
    value = value.strip()
    if not value or value.upper() in {"NA", "N/A", "NONE", "UNCLASSIFIED"}:
        return None
    token = value.split(";", 1)[0].strip()
    if token.startswith("d__"):
        token = token[3:]
    low = token.casefold()
    if low in {"bacteria", "bacterial"}:
        return "Bacteria"
    if low in {"archaea", "archaeal"}:
        return "Archaea"
    return None


def read_mag_domains(path: str | None, required_mags: set[str]) -> dict[str, str]:
    if not path:
        return {}
    fields, rows = read_tsv(path)
    low = {x.lower(): x for x in fields}
    id_col = low.get("id") or low.get("mag_id") or low.get("mag")
    domain_col = low.get("gtdb_domain") or low.get("domain")
    if not id_col or not domain_col:
        fail("MAG GTDB-Tk table needs ID/MAG_ID and GTDB_Domain columns")
    out: dict[str, str] = {}
    for row in rows:
        mag = row[id_col].strip()
        if not mag or mag in out:
            fail(f"blank/duplicate MAG ID in GTDB-Tk table: {mag}")
        domain = normalize_domain(row[domain_col])
        if domain:
            out[mag] = domain
    missing = sorted(required_mags - set(out))
    if missing:
        fail(f"{len(missing)} MAG(s) lack a resolved Bacteria/Archaea GTDB domain; first: {', '.join(missing[:10])}")
    return out


def marker_source_and_family(row: Mapping[str, str]) -> tuple[str, str]:
    marker = (row.get("hmm_name") or row.get("marker") or "").strip()
    source = (row.get("HMM_collection") or row.get("source") or row.get("hmm_source") or "").strip()
    if "|" in marker:
        prefix, rest = marker.split("|", 1)
        if prefix in {"Bacteria_71", "Archaea_76"}:
            source, marker = prefix, rest
    return source, marker


def empty_scg_annotation() -> dict[str, str]:
    return {
        "SCG_domain": "NA",
        "SCG_HMM_source": "NA",
        "n_SCG_hits": "NA",
        "n_distinct_SCG_families": "NA",
        "n_duplicated_SCG_families": "NA",
        "SCG_redundancy_fraction": "NA",
        "SCG_redundancy_observed": "NA",
        "duplicated_SCG_families": "NA",
    }


def load_scg(
    path: str | None,
    target_contigs: set[str],
    contig_to_mag: Mapping[str, str],
    mag_domains: Mapping[str, str],
) -> dict[str, dict[str, str]]:
    """Annotate SCG redundancy only for MAG-associated contigs."""
    if not path:
        return {contig: empty_scg_annotation() for contig in target_contigs}

    counts: dict[str, Counter[str]] = {c: Counter() for c in target_contigs if c in contig_to_mag}
    expected_source = {"Bacteria": "Bacteria_71", "Archaea": "Archaea_76"}
    fields, rows = read_tsv(path)
    if "contig" not in fields or not ({"hmm_name", "marker"} & set(fields)):
        fail("SCG table needs contig and hmm_name/marker columns")
    if counts and not mag_domains:
        fail("MAG GTDB-Tk taxonomy is required when SCG redundancy annotation is enabled")
    for row in rows:
        contig = row["contig"].strip()
        if contig not in counts:
            continue
        mag = contig_to_mag[contig]
        source, family = marker_source_and_family(row)
        wanted = expected_source[mag_domains[mag]]
        if source != wanted or not family:
            continue
        counts[contig][family] += 1

    out: dict[str, dict[str, str]] = {}
    for contig in target_contigs:
        if contig not in contig_to_mag:
            out[contig] = empty_scg_annotation()
            continue
        c = counts.get(contig, Counter())
        duplicated = sorted((fam, n) for fam, n in c.items() if n > 1)
        n_distinct = len(c)
        n_dup = len(duplicated)
        mag = contig_to_mag[contig]
        domain = mag_domains.get(mag)
        source = expected_source.get(domain, "NA")
        out[contig] = {
            "SCG_domain": domain or "NA",
            "SCG_HMM_source": source,
            "n_SCG_hits": str(sum(c.values())),
            "n_distinct_SCG_families": str(n_distinct),
            "n_duplicated_SCG_families": str(n_dup),
            "SCG_redundancy_fraction": f"{(n_dup / n_distinct):.6f}" if n_distinct else "NA",
            "SCG_redundancy_observed": "Y" if n_dup else "N",
            "duplicated_SCG_families": ";".join(f"{fam}x{n}" for fam, n in duplicated) or "NA",
        }
    return out


def event_bounds(row: Mapping[str, str], contig_length: int) -> tuple[int, int, int, int]:
    """Return the localized uncertain transition interval."""
    try:
        left = int(row.get("transition_interval_start_1based") or row["analysis_pos_1based"])
        right = int(row.get("transition_interval_end_1based") or row["analysis_pos_1based"])
    except (KeyError, ValueError):
        fail(f"candidate event lacks valid transition coordinates: {row.get('event_id','NA')}")
    if left > right: left, right = right, left
    left = max(1, min(left, contig_length)); right = max(1, min(right, contig_length))
    return left, right, left, right


def direct_call(stats: Mapping[str, Set[str]], min_q: int, validated_clip: bool = False) -> str:
    if len(stats["spanning"]) >= min_q:
        return "N"
    if len(stats["left"]) >= min_q and len(stats["right"]) >= min_q and (len(stats["clipped"]) >= min_q or validated_clip):
        return "Y"
    return "NA"


@dataclass(frozen=True, slots=True)
class ValidationWindow:
    """Precomputed search and anchor coordinates for one candidate event."""
    index: int
    contig: str
    start0: int
    end0: int
    left_start0: int
    left_end0: int
    right_start0: int
    right_end0: int
    transition_start1: int
    transition_end1: int


def overlapping_validation_windows(
    windows: Sequence[ValidationWindow], starts: Sequence[int],
    max_ends: Sequence[int], start0: int, end0: int,
) -> Iterator[ValidationWindow]:
    """Locate overlapping windows without scanning all events on a contig."""
    limit = bisect.bisect_left(starts, end0)
    for offset in range(limit - 1, -1, -1):
        if max_ends[offset] <= start0:
            break
        window = windows[offset]
        if window.end0 > start0:
            yield window


def validate_mode(
    bam: str, events: list[dict[str, str]], indices: list[int], lengths: Mapping[str, int],
    anchor: int, threads: int,
) -> dict[int, dict[str, Set[str]]]:
    stats: dict[int, dict[str, Set[str]]] = defaultdict(lambda: {"left": set(), "right": set(), "spanning": set(), "clipped": set()})
    if not indices:
        return stats

    windows_by_contig: dict[str, list[ValidationWindow]] = defaultdict(list)
    bed_windows: list[ValidationWindow] = []
    for idx in indices:
        row = events[idx]
        contig = row["contig"]
        length = lengths[contig]
        left, right, lo, hi = event_bounds(row, length)
        window = ValidationWindow(
            index=idx, contig=contig,
            start0=max(0, left - anchor - 1), end0=min(length, right + anchor),
            left_start0=max(0, left - anchor), left_end0=left,
            right_start0=right - 1, right_end0=min(length, right - 1 + anchor),
            transition_start1=lo, transition_end1=hi,
        )
        windows_by_contig[contig].append(window)
        bed_windows.append(window)

    window_index: dict[str, tuple[list[ValidationWindow], list[int], list[int]]] = {}
    for contig, windows in windows_by_contig.items():
        windows.sort(key=lambda w: (w.start0, w.index))
        starts = [w.start0 for w in windows]
        max_ends = []
        current_max = 0
        for window in windows:
            current_max = max(current_max, window.end0)
            max_ends.append(current_max)
        window_index[contig] = (windows, starts, max_ends)

    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", suffix=".bed", delete=False) as bed:
        bed_path = bed.name
        for window in bed_windows:
            bed.write(f"{window.contig}\t{window.start0}\t{window.end0}\n")
    try:
        cmd = ["samtools", "view", "-@", str(threads), "-F", "4", "-L", bed_path, bam]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, encoding="utf-8", errors="replace")
        assert proc.stdout is not None and proc.stderr is not None
        for raw in proc.stdout:
            f = raw.rstrip("\n\r").split("\t")
            if len(f) < 11 or f[5] == "*":
                continue
            qname, contig, cigar = f[0], f[2], f[5]
            if contig not in window_index:
                continue
            try:
                ref_start0 = int(f[3]) - 1
            except ValueError:
                continue
            ops = cigar_ops(cigar)
            ref_end0 = ref_start0 + sum(n for n, op in ops if op in REF_OPS)
            runs = matched_runs(ref_start0, ops)
            left_clip = terminal_clip(ops, True) > 0
            right_clip = terminal_clip(ops, False) > 0
            windows, starts, max_ends = window_index[contig]
            for window in overlapping_validation_windows(windows, starts, max_ends, ref_start0, ref_end0):
                lbp = matched_bp(runs, window.left_start0, window.left_end0)
                rbp = matched_bp(runs, window.right_start0, window.right_end0)
                if lbp >= anchor:
                    stats[window.index]["left"].add(qname)
                if rbp >= anchor:
                    stats[window.index]["right"].add(qname)
                if lbp >= anchor and rbp >= anchor:
                    stats[window.index]["spanning"].add(qname)
                if ((left_clip and window.transition_start1 - 1 <= ref_start0 + 1 <= window.transition_end1 + 1)
                        or (right_clip and window.transition_start1 - 1 <= ref_end0 <= window.transition_end1 + 1)):
                    stats[window.index]["clipped"].add(qname)
        proc.stdout.close()
        stderr = proc.stderr.read()
        code = proc.wait()
        if code != 0:
            fail(f"samtools view failed for {bam}: {stderr}")
    finally:
        try:
            os.unlink(bed_path)
        except FileNotFoundError:
            pass
    return stats


def evidence_label(call: str, high_conf: bool) -> str:
    if call == "Y": return "HIGH_CONFIDENCE_HIFI_STRUCTURAL_CHIMERA" if high_conf else "HIFI_STRUCTURAL_CHIMERA"
    if call == "N": return "HIFI_CONTINUITY_SUPPORTED"
    return "UNRESOLVED_HIFI_STRUCTURAL_EVIDENCE"


def _validate_main() -> None:
    args = _validate_parse_args()
    if args.anchor_bp < 1 or args.min_support_qnames < 1 or args.threads < 1 or args.parallel_bam_jobs < 1:
        fail("anchor/min-support/threads/parallel-bam-jobs must be positive")
    for p, label in ((args.candidate_events, "candidate events"), (args.candidate_contigs, "candidate contigs"), (args.mapping_manifest, "mapping manifest"), (args.mag_map, "MAG-contig map")):
        _validate_require_file(p, label)
    _validate_require_file(args.clipping, "validated clipping", optional=True)
    _validate_require_file(args.scg_tsv, "SCG table", optional=True)
    _validate_require_file(args.mag_taxonomy_tsv, "MAG GTDB-Tk taxonomy table", optional=True)

    event_fields, events = read_tsv(args.candidate_events)
    contig_fields, contig_rows = read_tsv(args.candidate_contigs)
    required_event = {"event_id", "MAG_ID", "contig", "left_block_end_1based", "right_block_start_1based", "taxonomic_transition", "shift_rank"}
    if required_event - set(event_fields): fail(f"candidate event table lacks columns: {sorted(required_event-set(event_fields))}")
    required_contig = {"MAG_ID", "contig", "length", "Taxonomic_shift"}
    if required_contig - set(contig_fields): fail(f"candidate contig table lacks columns: {sorted(required_contig-set(contig_fields))}")

    lengths = {r["contig"]: int(r["length"]) for r in contig_rows}
    target_contigs = set(lengths)
    contig_to_mag = _validate_read_mag_map(args.mag_map)
    routes = load_routes(args.mapping_manifest, target_contigs)
    clipping = load_clipping(args.clipping)
    target_mags = {contig_to_mag[c] for c in target_contigs if c in contig_to_mag}
    mag_domains = read_mag_domains(args.mag_taxonomy_tsv, target_mags) if args.scg_tsv and target_mags else {}
    scg = load_scg(args.scg_tsv, target_contigs, contig_to_mag, mag_domains)

    overlap_clip: dict[int, list[dict[str, str]]] = defaultdict(list)
    for idx, row in enumerate(events):
        contig = row["contig"]; left, right, lo, hi = event_bounds(row, lengths[contig])
        for cr in clipping.get(contig, []):
            try: pos = int(cr["pos_1based"])
            except ValueError: continue
            if lo - 1 <= pos <= hi + 1:
                overlap_clip[idx].append(cr)

    by_bam_ps: dict[str, list[int]] = defaultdict(list)
    for idx, row in enumerate(events):
        route = routes[row["contig"]]
        by_bam_ps[route["primary_supp"]].append(idx)

    ps_stats: dict[int, dict[str, Set[str]]] = defaultdict(lambda: {"left": set(), "right": set(), "spanning": set(), "clipped": set()})
    bam_items = list(by_bam_ps.items())
    jobs = min(args.parallel_bam_jobs, len(bam_items)) if bam_items else 1
    threads_per_bam = max(1, args.threads // jobs)
    if jobs <= 1:
        for bam, idxs in bam_items:
            ps_stats.update(validate_mode(bam, events, idxs, lengths, args.anchor_bp, threads_per_bam))
    else:
        with ThreadPoolExecutor(max_workers=jobs) as executor:
            futures = {
                executor.submit(
                    validate_mode, bam, events, idxs, lengths, args.anchor_bp, threads_per_bam
                ): bam
                for bam, idxs in bam_items
            }
            for future in as_completed(futures):
                ps_stats.update(future.result())

    appended_event = [
        "source_assembler", "junction_left_boundary_1based", "junction_right_boundary_1based",
        "hifi_anchor_bp", "hifi_min_support_qnames",
        "hifi_primary_supp_left_anchor_qnames", "hifi_primary_supp_right_anchor_qnames", "hifi_primary_supp_spanning_qnames", "hifi_primary_supp_breakpoint_clipped_qnames",
        "HIFI_PRIMARY_SUPP_STRUCTURAL_CALL",
        "Stringent_clipping_support", "Structural_chimera", "High_confidence_structural_chimera",
        "SCG_domain", "SCG_HMM_source", "n_SCG_hits", "n_distinct_SCG_families", "n_duplicated_SCG_families", "SCG_redundancy_fraction", "SCG_redundancy_observed", "duplicated_SCG_families",
    ]
    final_event_fields = list(event_fields)
    for f in appended_event:
        if f not in final_event_fields: final_event_fields.append(f)

    final_events = []
    for idx, row0 in enumerate(events):
        row = dict(row0); contig = row["contig"]; route = routes[contig]
        left, right, _, _ = event_bounds(row, lengths[contig])
        clips = overlap_clip.get(idx, [])
        validated_clip = bool(clips)
        stringent = any((c.get("Stringent_clipping") or "").strip().upper() == "Y" for c in clips)
        p = ps_stats[idx]
        ps_call = direct_call(p, args.min_support_qnames, validated_clip=validated_clip)
        high_conf = ps_call == "Y" and stringent

        row["assessment_class"] = (
            "HIGH_CONFIDENCE_STRUCTURAL_CHIMERA" if high_conf else
            "STRUCTURAL_CHIMERA" if ps_call == "Y" else
            "HIFI_CONTINUITY_SUPPORTED" if ps_call == "N" else
            "UNRESOLVED_STRUCTURAL_CHIMERA"
        )
        row["default_clipping_support"] = "Y" if clips else "N"
        row["validated_clipping_support"] = "Y" if clips else "N"
        row["n_supporting_validated_clipping_positions"] = str(len(clips))
        row["supporting_validated_clipping_positions_1based"] = ";".join(c["pos_1based"] for c in clips) or "NA"
        row["Stringent_clipping_support"] = "Y" if stringent else "N"
        row["CHIMERA_EVIDENCE"] = evidence_label(ps_call, high_conf)
        row["source_assembler"] = route["assembler"]
        row["junction_left_boundary_1based"] = str(left); row["junction_right_boundary_1based"] = str(right)
        row["hifi_anchor_bp"] = str(args.anchor_bp); row["hifi_min_support_qnames"] = str(args.min_support_qnames)
        row["hifi_primary_supp_left_anchor_qnames"] = str(len(p["left"])); row["hifi_primary_supp_right_anchor_qnames"] = str(len(p["right"]))
        row["hifi_primary_supp_spanning_qnames"] = str(len(p["spanning"])); row["hifi_primary_supp_breakpoint_clipped_qnames"] = str(len(p["clipped"]))
        row["HIFI_PRIMARY_SUPP_STRUCTURAL_CALL"] = ps_call
        row["Structural_chimera"] = ps_call; row["High_confidence_structural_chimera"] = "Y" if high_conf else ("N" if ps_call in {"Y", "N"} else "NA")
        row.update(scg[contig])
        final_events.append(row)

    structural_events = [r for r in final_events if r["Structural_chimera"] == "Y"]
    high_conf_events = [r for r in final_events if r["High_confidence_structural_chimera"] == "Y"]
    unresolved_events = [r for r in final_events if r["Structural_chimera"] == "NA"]

    appended_contig = [
        "source_assembler",
        "n_structural_chimera_events", "n_high_confidence_structural_chimera_events", "n_hifi_continuity_supported_events", "n_unresolved_structural_events",
        "structural_chimera_positions_1based", "high_confidence_structural_chimera_positions_1based", "unresolved_structural_positions_1based",
        "Structural_chimera", "High_confidence_structural_chimera",
        "SCG_domain", "SCG_HMM_source", "n_SCG_hits", "n_distinct_SCG_families", "n_duplicated_SCG_families", "SCG_redundancy_fraction", "SCG_redundancy_observed", "duplicated_SCG_families",
    ]
    final_contig_fields = list(contig_fields)
    for f in appended_contig:
        if f not in final_contig_fields: final_contig_fields.append(f)

    events_by_contig: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in final_events: events_by_contig[row["contig"]].append(row)
    final_contigs = []
    for base in contig_rows:
        row = dict(base); contig = row["contig"]; ers = events_by_contig.get(contig, [])
        y = [e for e in ers if e["Structural_chimera"] == "Y"]
        n = [e for e in ers if e["Structural_chimera"] == "N"]
        na = [e for e in ers if e["Structural_chimera"] == "NA"]
        hc = [e for e in ers if e["High_confidence_structural_chimera"] == "Y"]
        clip_y = [e for e in y if e["validated_clipping_support"] == "Y"]
        independent_y = [e for e in y if e["validated_clipping_support"] != "Y"]
        independent_shift = [e for e in ers if e["validated_clipping_support"] != "Y"]
        tax_shift = (row.get("Taxonomic_shift") or "NA").strip().upper()
        if y: structural = "Y"
        elif na: structural = "NA"
        elif ers: structural = "N"
        elif tax_shift == "N": structural = "N"
        else: structural = "NA"
        high_call = "Y" if hc else ("N" if structural in {"Y", "N"} else "NA")
        if structural == "Y": evidence = "HIGH_CONFIDENCE_HIFI_STRUCTURAL_CHIMERA" if hc else "HIFI_STRUCTURAL_CHIMERA"
        elif structural == "N" and ers: evidence = "HIFI_CONTINUITY_SUPPORTED"
        elif structural == "NA" and ers: evidence = "UNRESOLVED_HIFI_STRUCTURAL_EVIDENCE"
        elif tax_shift == "N": evidence = "NO_TAXONOMIC_SHIFT"
        else: evidence = "NOT_ASSESSABLE"
        if structural == "Y": status = "STRUCTURAL_CHIMERA_DETECTED"
        elif structural == "N" and ers: status = "HIFI_CONTINUITY_SUPPORTED_FOR_TAXONOMIC_TRANSITIONS"
        elif structural == "NA" and ers: status = "UNRESOLVED_TAXONOMIC_TRANSITION_STRUCTURE"
        elif tax_shift == "N": status = "NO_TAXONOMIC_SHIFT"
        else: status = "NOT_ASSESSABLE_INSUFFICIENT_TAXONOMIC_RESOLUTION"
        row["n_structural_chimera_events"] = str(len(y)); row["n_high_confidence_structural_chimera_events"] = str(len(hc))
        row["n_hifi_continuity_supported_events"] = str(len(n)); row["n_unresolved_structural_events"] = str(len(na))
        row["structural_chimera_positions_1based"] = ";".join(e["analysis_pos_1based"] for e in y) or "NA"
        row["high_confidence_structural_chimera_positions_1based"] = ";".join(e["analysis_pos_1based"] for e in hc) or "NA"
        row["unresolved_structural_positions_1based"] = ";".join(e["analysis_pos_1based"] for e in na) or "NA"
        row["source_assembler"] = routes[contig]["assembler"]
        row["Structural_chimera"] = structural; row["High_confidence_structural_chimera"] = high_call
        row.update(scg[contig])
        row["n_taxonomic_chimera_with_clipping_support_events"] = str(len(clip_y))
        row["n_clipping_independent_taxonomic_shift_events"] = str(len(independent_shift))
        row["taxonomic_chimera_with_clipping_support_positions_1based"] = ";".join(e["analysis_pos_1based"] for e in clip_y) or "NA"
        row["clipping_independent_taxonomic_shift_positions_1based"] = ";".join(e["analysis_pos_1based"] for e in independent_shift) or "NA"
        row["Taxonomic_chimera_with_clipping_support"] = "Y" if clip_y else ("N" if structural in {"Y", "N"} else "NA")
        row["Clipping_independent_taxonomic_shift"] = "Y" if independent_shift else ("N" if tax_shift in {"Y", "N"} else "NA")
        row["CHIMERA_EVIDENCE"] = evidence
        row["High_confidence_taxonomic_chimera"] = high_call
        row["Clipping_independent_taxonomic_chimera"] = "Y" if independent_y else ("N" if structural in {"Y", "N"} else "NA")
        row["GTDB_assessment_status"] = status
        row["Taxonomic_chimera"] = structural
        final_contigs.append(row)

    structural_contigs = [r for r in final_contigs if r["Structural_chimera"] == "Y"]
    high_conf_contigs = [r for r in final_contigs if r["High_confidence_structural_chimera"] == "Y"]
    unresolved_contigs = [r for r in final_contigs if r["Structural_chimera"] == "NA" and r.get("Taxonomic_shift") == "Y"]
    scg_audit = [r for r in final_contigs if r["SCG_redundancy_observed"] == "Y"]

    _validate_atomic_write(args.events_output, final_event_fields, final_events)
    _validate_atomic_write(args.structural_events_output, final_event_fields, structural_events)
    _validate_atomic_write(args.high_confidence_events_output, final_event_fields, high_conf_events)
    _validate_atomic_write(args.unresolved_events_output, final_event_fields, unresolved_events)
    _validate_atomic_write(args.contig_output, final_contig_fields, final_contigs)
    _validate_atomic_write(args.structural_contig_output, final_contig_fields, structural_contigs)
    _validate_atomic_write(args.high_confidence_contig_output, final_contig_fields, high_conf_contigs)
    _validate_atomic_write(args.unresolved_contig_output, final_contig_fields, unresolved_contigs)
    _validate_atomic_write(args.scg_audit_output, final_contig_fields, scg_audit)

    print(f"Taxonomic-transition candidates: {len(final_events)}")
    print(f"HiFi-supported structural chimeras: {len(structural_events)}")
    print(f"High-confidence structural chimeras (Stringent_clipping=Y): {len(high_conf_events)}")
    print(f"HiFi continuity-supported transitions: {sum(r['Structural_chimera']=='N' for r in final_events)}")
    print(f"Unresolved transitions: {len(unresolved_events)}")
    print(f"Contigs with any SCG family redundancy (annotation only): {len(scg_audit)}")


def cmd_validate(argv):
    """Validate transitions with HiFi alignments and SCG redundancy"""
    sys.argv = [sys.argv[0] + " validate", *argv]
    return _validate_main()


def _mag_fail(msg: str) -> None: raise SystemExit(f"ERROR: {msg}")


def _mag_parse_args():
    p=argparse.ArgumentParser()
    p.add_argument('--contig-table', required=True)
    p.add_argument('--mag-map', required=True)
    p.add_argument('--mag-contig-out', required=True)
    p.add_argument('--mag-summary-out', required=True)
    return p.parse_args()


def _mag_read_tsv(path):
    with open(path, encoding='utf-8-sig', newline='') as h:
        r=csv.DictReader(h, delimiter='\t'); return list(r.fieldnames or []), [dict(x) for x in r]


def _mag_atomic_write(path: str, fields: Sequence[str], rows: Iterable[Mapping[str,str]]):
    dst=Path(path); dst.parent.mkdir(parents=True, exist_ok=True)
    fd,tmp=tempfile.mkstemp(prefix=dst.name+'.tmp.', dir=str(dst.parent), text=True)
    try:
        with os.fdopen(fd,'w',encoding='utf-8',newline='') as h:
            w=csv.DictWriter(h,fieldnames=list(fields),delimiter='\t',lineterminator='\n',extrasaction='raise'); w.writeheader(); w.writerows(rows)
        os.replace(tmp,dst)
    except BaseException:
        try: os.unlink(tmp)
        except FileNotFoundError: pass
        raise


def _mag_read_mag_map(path):
    fields,rows=_mag_read_tsv(path); low={x.lower():x for x in fields}
    mcol=low.get('mag_id') or low.get('mag') or low.get('id'); ccol=low.get('contig') or low.get('contig_id')
    if not mcol or not ccol: _mag_fail('MAG-contig map needs MAG_ID/ID and contig columns')
    c2m={}; m2c=defaultdict(list)
    for r in rows:
        m=r[mcol].strip(); c=r[ccol].strip()
        if not m or not c or c in c2m: _mag_fail(f'invalid/duplicate MAG-contig row: {m}/{c}')
        c2m[c]=m; m2c[m].append(c)
    return c2m,m2c


def intv(row,key):
    try:return int(row.get(key,'0') or 0)
    except ValueError: _mag_fail(f'invalid integer {key} for {row.get("contig","NA")}')


def aggregate_call(rows):
    calls=[r['Structural_chimera'].strip().upper() for r in rows]
    if 'Y' in calls:return 'Y'
    if 'NA' in calls:return 'NA'
    return 'N'


def evidence(rows, call):
    if call=='Y':
        return 'HIGH_CONFIDENCE_HIFI_STRUCTURAL_CHIMERA' if any(r.get('High_confidence_structural_chimera')=='Y' for r in rows) else 'HIFI_STRUCTURAL_CHIMERA'
    if call=='NA': return 'UNRESOLVED_HIFI_STRUCTURAL_EVIDENCE'
    if any(r.get('Taxonomic_shift')=='Y' for r in rows): return 'HIFI_CONTINUITY_SUPPORTED'
    return 'NO_TAXONOMIC_SHIFT'


def _mag_main():
    a=_mag_parse_args()
    fields,rows=_mag_read_tsv(a.contig_table)
    required={'taxonomy_source','taxonomy_release','MAG_ID','contig','length','Taxonomic_shift','Taxonomic_chimera','Structural_chimera','High_confidence_structural_chimera','CHIMERA_EVIDENCE'}
    if required-set(fields): _mag_fail(f'contig table lacks columns: {sorted(required-set(fields))}')
    c2m,m2c=_mag_read_mag_map(a.mag_map)
    all_by={r['contig']:r for r in rows}
    missing_mag_contigs = sorted(set(c2m) - set(all_by))
    if missing_mag_contigs:
        _mag_fail(f'contig table is missing MAG contigs: {missing_mag_contigs[:5]}')
    by={c:all_by[c] for c in c2m}
    for c,r in by.items():
        if r['MAG_ID']!=c2m[c]: _mag_fail(f'MAG mismatch for {c}: {r["MAG_ID"]} vs {c2m[c]}')
    sorted_rows=sorted(by.values(),key=lambda r:(r['MAG_ID'],r['contig']))
    _mag_atomic_write(a.mag_contig_out, fields, sorted_rows)

    out_fields=[
        'taxonomy_source','taxonomy_release','MAG_ID','n_contigs','total_length_bp',
        'n_taxonomic_shift_contigs','n_structural_chimera_contigs','n_high_confidence_structural_chimera_contigs','n_unresolved_structural_contigs',
        'n_taxonomic_shift_events','n_structural_chimera_events','n_high_confidence_structural_chimera_events','n_hifi_continuity_supported_events','n_unresolved_structural_events',
        'n_SCG_redundancy_contigs',
        'structural_chimera_contig_ids','high_confidence_structural_chimera_contig_ids','unresolved_structural_contig_ids','SCG_redundancy_contig_ids',
        'Taxonomic_shift','Structural_chimera','Taxonomic_chimera','High_confidence_structural_chimera','CHIMERA_EVIDENCE','GTDB_assessment_status',
    ]
    mag_rows=[]
    for mag in sorted(m2c):
        rs=[by[c] for c in sorted(m2c[mag])]
        release={r['taxonomy_release'] for r in rs}; source={r['taxonomy_source'] for r in rs}
        if len(release)!=1 or len(source)!=1: _mag_fail(f'mixed taxonomy metadata within {mag}')
        structural=aggregate_call(rs)
        taxshift='Y' if any(r['Taxonomic_shift']=='Y' for r in rs) else ('NA' if any(r['Taxonomic_shift']=='NA' for r in rs) else 'N')
        high='Y' if any(r['High_confidence_structural_chimera']=='Y' for r in rs) else ('N' if structural in {'Y','N'} else 'NA')
        y=[r['contig'] for r in rs if r['Structural_chimera']=='Y']
        hc=[r['contig'] for r in rs if r['High_confidence_structural_chimera']=='Y']
        unresolved=[r['contig'] for r in rs if r['Structural_chimera']=='NA' and r['Taxonomic_shift']=='Y']
        scg=[r['contig'] for r in rs if r.get('SCG_redundancy_observed')=='Y']
        ev=evidence(rs,structural)
        status=('STRUCTURAL_CHIMERA_DETECTED' if structural=='Y' else 'UNRESOLVED_TAXONOMIC_TRANSITION_STRUCTURE' if structural=='NA' and taxshift=='Y' else 'HIFI_CONTINUITY_SUPPORTED_FOR_TAXONOMIC_TRANSITIONS' if taxshift=='Y' else 'NO_TAXONOMIC_SHIFT')
        mag_rows.append({
            'taxonomy_source':next(iter(source)),'taxonomy_release':next(iter(release)),'MAG_ID':mag,'n_contigs':str(len(rs)),
            'total_length_bp':str(sum(int(r['length']) for r in rs)),
            'n_taxonomic_shift_contigs':str(sum(r['Taxonomic_shift']=='Y' for r in rs)),
            'n_structural_chimera_contigs':str(len(y)),'n_high_confidence_structural_chimera_contigs':str(len(hc)),'n_unresolved_structural_contigs':str(len(unresolved)),
            'n_taxonomic_shift_events':str(sum(intv(r,'n_taxonomic_shift_events') for r in rs)),
            'n_structural_chimera_events':str(sum(intv(r,'n_structural_chimera_events') for r in rs)),
            'n_high_confidence_structural_chimera_events':str(sum(intv(r,'n_high_confidence_structural_chimera_events') for r in rs)),
            'n_hifi_continuity_supported_events':str(sum(intv(r,'n_hifi_continuity_supported_events') for r in rs)),
            'n_unresolved_structural_events':str(sum(intv(r,'n_unresolved_structural_events') for r in rs)),
            'n_SCG_redundancy_contigs':str(len(scg)),
            'structural_chimera_contig_ids':';'.join(y) or 'NA','high_confidence_structural_chimera_contig_ids':';'.join(hc) or 'NA',
            'unresolved_structural_contig_ids':';'.join(unresolved) or 'NA','SCG_redundancy_contig_ids':';'.join(scg) or 'NA',
            'Taxonomic_shift':taxshift,'Structural_chimera':structural,'Taxonomic_chimera':structural,'High_confidence_structural_chimera':high,
            'CHIMERA_EVIDENCE':ev,'GTDB_assessment_status':status,
        })
    _mag_atomic_write(a.mag_summary_out,out_fields,mag_rows)
    print(f'MAGs summarized: {len(mag_rows)}')
    print(f'MAGs with structurally validated taxonomic chimera: {sum(r["Structural_chimera"]=="Y" for r in mag_rows)}')
    print(f'MAGs unresolved at one or more taxonomic transitions: {sum(r["Structural_chimera"]=="NA" and r["Taxonomic_shift"]=="Y" for r in mag_rows)}')


def cmd_mag(argv):
    """Summarize transitions per MAG"""
    sys.argv = [sys.argv[0] + " mag", *argv]
    return _mag_main()


NA_VALUES = {"", "NA", "N/A", "NONE", "."}


def _reclassify_fail(msg: str) -> None:
    raise SystemExit(f"ERROR: {msg}")


def _reclassify_parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--regions-in", required=True,
                   help="Existing taxonomic_transition_structural_validation.*.events.tsv")
    p.add_argument("--clipping", required=True,
                   help="Current 03_clipping.sh final clipping table; only Clipping=Y sites are used")
    p.add_argument("--contigs-in", required=True,
                   help="Existing contig_structural_chimera.*.final.tsv")
    p.add_argument("--mag-contigs-in", required=True,
                   help="Existing MAG_contig_structural_chimera.*.final.tsv")
    p.add_argument("--mags-in", required=True,
                   help="Existing MAG_structural_chimera.*.final.tsv")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--aux-dir", help="directory for chimera.MAG_contig.tsv and run_metadata.tsv (default: --out-dir)")
    p.add_argument("--release", default="GTDB_R232")
    return p.parse_args()


def _reclassify_read_tsv(path: str) -> Tuple[List[str], List[Dict[str, str]]]:
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        _reclassify_fail(f"missing or empty input: {path}")
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = [dict(row) for row in reader]
    if not fields:
        _reclassify_fail(f"missing header: {path}")
    return fields, rows


def _reclassify_atomic_write(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({f: row.get(f, "NA") for f in fields})
    os.replace(tmp, path)


def ensure_fields(fields: List[str], additions: Sequence[str]) -> List[str]:
    out = list(fields)
    for f in additions:
        if f not in out:
            out.append(f)
    return out


def pos_sort_key(value: str) -> Tuple[int, str]:
    try:
        return (int(value), value)
    except (TypeError, ValueError):
        return (10**30, value or "")


def join_values(rows: Sequence[Mapping[str, str]], key: str) -> str:
    vals = sorted({r.get(key, "") for r in rows if r.get(key, "") not in NA_VALUES}, key=pos_sort_key)
    return ";".join(vals) or "NA"


def join_contigs(values: Iterable[str]) -> str:
    vals = sorted({v for v in values if v not in NA_VALUES})
    return ";".join(vals) or "NA"


def load_validated_clipping(path: str) -> Dict[str, List[Dict[str, str]]]:
    fields, rows = _reclassify_read_tsv(path)
    required = {"contig", "pos_1based", "Clipping", "Stringent_clipping"}
    missing = required - set(fields)
    if missing:
        _reclassify_fail(f"clipping input lacks columns: {sorted(missing)}")
    out: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    seen = set()
    for row in rows:
        clip = (row.get("Clipping") or "").strip().upper()
        stringent = (row.get("Stringent_clipping") or "").strip().upper()
        if clip not in {"Y", "N"}:
            _reclassify_fail(f"invalid Clipping flag at {row.get('contig','NA')}:{row.get('pos_1based','NA')}")
        if stringent not in {"Y", "N", "NA"}:
            _reclassify_fail(f"invalid Stringent_clipping flag at {row.get('contig','NA')}:{row.get('pos_1based','NA')}")
        if stringent == "Y" and clip != "Y":
            _reclassify_fail(f"Stringent_clipping=Y without Clipping=Y at {row.get('contig','NA')}:{row.get('pos_1based','NA')}")
        if clip != "Y":
            continue
        contig = (row.get("contig") or "").strip()
        try:
            pos = int(row.get("pos_1based") or "")
        except ValueError:
            _reclassify_fail(f"invalid clipping position: {contig}:{row.get('pos_1based','NA')}")
        key = (contig, pos)
        if key in seen:
            _reclassify_fail(f"duplicate validated clipping site: {contig}:{pos}")
        seen.add(key)
        out[contig].append(dict(row))
    for contig in out:
        out[contig].sort(key=lambda r: int(r["pos_1based"]))
    return out


def annotate_region_with_clipping(row: Dict[str, str], clipping: Mapping[str, Sequence[Mapping[str, str]]]) -> None:
    contig = row["contig"]
    try:
        lo = int(row.get("transition_interval_start_1based") or row["analysis_pos_1based"])
        hi = int(row.get("transition_interval_end_1based") or row["analysis_pos_1based"])
    except (KeyError, ValueError):
        _reclassify_fail(f"invalid transition interval: {row.get('event_id','NA')}")
    if lo > hi:
        lo, hi = hi, lo
    supporting = []
    for clip in clipping.get(contig, []):
        pos = int(clip["pos_1based"])
        if lo - 1 <= pos <= hi + 1:
            supporting.append(clip)
    row["validated_clipping_support"] = "Y" if supporting else "N"
    row["n_supporting_validated_clipping_positions"] = str(len(supporting))
    row["supporting_validated_clipping_positions_1based"] = ";".join(c["pos_1based"] for c in supporting) or "NA"
    row["Stringent_clipping_support"] = "Y" if any((c.get("Stringent_clipping") or "").strip().upper() == "Y" for c in supporting) else "N"


def final_region_call(row: Mapping[str, str]) -> Tuple[str, str, str, str]:
    if (row.get("taxonomic_transition") or "").strip().upper() != "Y":
        _reclassify_fail(f"non-taxonomic-shift row in region table: {row.get('event_id','NA')}")
    validated = (row.get("validated_clipping_support") or "").strip().upper()
    stringent = (row.get("Stringent_clipping_support") or "").strip().upper()
    if validated not in {"Y", "N"}:
        _reclassify_fail(f"invalid validated_clipping_support: {row.get('event_id','NA')}/{validated}")
    if stringent not in {"Y", "N"}:
        _reclassify_fail(f"invalid Stringent_clipping_support: {row.get('event_id','NA')}/{stringent}")
    if stringent == "Y" and validated != "Y":
        _reclassify_fail(f"Stringent_clipping_support=Y without validated clipping: {row.get('event_id','NA')}")

    chimera = "Y" if validated == "Y" else "N"
    high = "Y" if stringent == "Y" else "N"
    if high == "Y":
        evidence = "TAXONOMIC_SHIFT_WITH_STRINGENT_VALIDATED_CLIPPING"
        status = "HIGH_CONFIDENCE_TAXONOMIC_CHIMERA"
    elif chimera == "Y":
        evidence = "TAXONOMIC_SHIFT_WITH_VALIDATED_CLIPPING"
        status = "TAXONOMIC_CHIMERA"
    else:
        evidence = "TAXONOMIC_SHIFT_WITHOUT_VALIDATED_CLIPPING"
        status = "TAXONOMIC_SHIFT_ONLY"
    return chimera, high, evidence, status


def classify_contig(base: Dict[str, str], regions: Sequence[Mapping[str, str]]) -> Dict[str, str]:
    row = dict(base)
    tax_shift = (row.get("Taxonomic_shift") or "NA").strip().upper()
    if tax_shift not in {"Y", "N", "NA"}:
        _reclassify_fail(f"invalid Taxonomic_shift on contig {row.get('contig','NA')}: {tax_shift}")

    chim = [r for r in regions if r["Taxonomic_chimera"] == "Y"]
    high = [r for r in regions if r["High_confidence_taxonomic_chimera"] == "Y"]
    no_clip = [r for r in regions if r["Taxonomic_chimera"] == "N"]

    if tax_shift == "Y":
        call = "Y" if chim else "N"
        high_call = "Y" if high else "N"
        if high:
            evidence = "TAXONOMIC_SHIFT_WITH_STRINGENT_VALIDATED_CLIPPING"
            status = "HIGH_CONFIDENCE_TAXONOMIC_CHIMERA"
        elif chim:
            evidence = "TAXONOMIC_SHIFT_WITH_VALIDATED_CLIPPING"
            status = "TAXONOMIC_CHIMERA"
        else:
            evidence = "TAXONOMIC_SHIFT_WITHOUT_VALIDATED_CLIPPING"
            status = "TAXONOMIC_SHIFT_ONLY"
    elif tax_shift == "N":
        call = high_call = "N"
        evidence = status = "NO_TAXONOMIC_SHIFT"
    else:
        call = high_call = "NA"
        evidence = status = "NOT_ASSESSABLE"

    row["n_taxonomic_shift_regions"] = str(len(regions))
    row["n_taxonomic_chimera_regions"] = str(len(chim))
    row["n_high_confidence_taxonomic_chimera_regions"] = str(len(high))
    row["n_taxonomic_shift_only_regions"] = str(len(no_clip))
    row["taxonomic_chimera_positions_1based"] = join_values(chim, "analysis_pos_1based")
    row["high_confidence_taxonomic_chimera_positions_1based"] = join_values(high, "analysis_pos_1based")
    row["taxonomic_shift_only_positions_1based"] = join_values(no_clip, "analysis_pos_1based")
    row["Taxonomic_chimera_with_clipping_support"] = call
    row["Clipping_independent_taxonomic_shift"] = "Y" if no_clip else ("N" if tax_shift in {"Y", "N"} else "NA")
    row["High_confidence_taxonomic_chimera"] = high_call
    row["Taxonomic_chimera"] = call
    row["CHIMERA_EVIDENCE"] = evidence
    row["GTDB_assessment_status"] = status
    row["final_taxonomic_chimera_policy"] = "taxonomic_shift_plus_validated_Clipping_Y"
    row["high_confidence_policy"] = "taxonomic_shift_plus_Stringent_clipping_Y"
    return row


def classify_mag(base: Dict[str, str], contigs: Sequence[Mapping[str, str]], regions: Sequence[Mapping[str, str]]) -> Dict[str, str]:
    row = dict(base)
    tax_shift_contigs = [c for c in contigs if c.get("Taxonomic_shift") == "Y"]
    chim_contigs = [c for c in contigs if c.get("Taxonomic_chimera") == "Y"]
    high_contigs = [c for c in contigs if c.get("High_confidence_taxonomic_chimera") == "Y"]
    shift_only_contigs = [c for c in contigs if c.get("Taxonomic_shift") == "Y" and c.get("Taxonomic_chimera") != "Y"]
    chim_regions = [r for r in regions if r["Taxonomic_chimera"] == "Y"]
    high_regions = [r for r in regions if r["High_confidence_taxonomic_chimera"] == "Y"]
    shift_only_regions = [r for r in regions if r["Taxonomic_chimera"] == "N"]

    any_assessable = any(c.get("Taxonomic_chimera_assessable") == "Y" for c in contigs)
    tax_shift = "Y" if tax_shift_contigs else ("N" if any_assessable else "NA")
    chim_call = "Y" if chim_contigs else ("N" if tax_shift in {"Y", "N"} else "NA")
    high_call = "Y" if high_contigs else ("N" if chim_call in {"Y", "N"} else "NA")

    if high_call == "Y":
        evidence = "TAXONOMIC_SHIFT_WITH_STRINGENT_VALIDATED_CLIPPING"
        status = "HIGH_CONFIDENCE_TAXONOMIC_CHIMERA"
    elif chim_call == "Y":
        evidence = "TAXONOMIC_SHIFT_WITH_VALIDATED_CLIPPING"
        status = "TAXONOMIC_CHIMERA"
    elif tax_shift == "Y":
        evidence = "TAXONOMIC_SHIFT_WITHOUT_VALIDATED_CLIPPING"
        status = "TAXONOMIC_SHIFT_ONLY"
    elif tax_shift == "N":
        evidence = status = "NO_TAXONOMIC_SHIFT"
    else:
        evidence = status = "NOT_ASSESSABLE"

    row["n_taxonomic_shift_contigs"] = str(len(tax_shift_contigs))
    row["n_taxonomic_chimera_contigs"] = str(len(chim_contigs))
    row["n_high_confidence_taxonomic_chimera_contigs"] = str(len(high_contigs))
    row["n_taxonomic_shift_only_contigs"] = str(len(shift_only_contigs))
    row["n_taxonomic_shift_regions"] = str(len(regions))
    row["n_taxonomic_chimera_regions"] = str(len(chim_regions))
    row["n_high_confidence_taxonomic_chimera_regions"] = str(len(high_regions))
    row["n_taxonomic_shift_only_regions"] = str(len(shift_only_regions))
    row["taxonomic_chimera_contig_ids"] = join_contigs(c.get("contig", "") for c in chim_contigs)
    row["high_confidence_taxonomic_chimera_contig_ids"] = join_contigs(c.get("contig", "") for c in high_contigs)
    row["taxonomic_shift_only_contig_ids"] = join_contigs(c.get("contig", "") for c in shift_only_contigs)
    row["Taxonomic_shift"] = tax_shift
    row["Taxonomic_chimera"] = chim_call
    row["High_confidence_taxonomic_chimera"] = high_call
    row["CHIMERA_EVIDENCE"] = evidence
    row["GTDB_assessment_status"] = status
    row["final_taxonomic_chimera_policy"] = "taxonomic_shift_plus_validated_Clipping_Y"
    row["high_confidence_policy"] = "taxonomic_shift_plus_Stringent_clipping_Y"
    return row


def _reclassify_main() -> None:
    args = _reclassify_parse_args()
    region_fields, region_rows = _reclassify_read_tsv(args.regions_in)
    clipping = load_validated_clipping(args.clipping)
    contig_fields, contig_rows = _reclassify_read_tsv(args.contigs_in)
    mag_contig_fields, mag_contig_rows = _reclassify_read_tsv(args.mag_contigs_in)
    mag_fields, mag_rows = _reclassify_read_tsv(args.mags_in)

    required_regions = {
        "event_id", "MAG_ID", "contig", "taxonomic_transition", "analysis_pos_1based",
        "transition_interval_start_1based", "transition_interval_end_1based",
    }
    missing = required_regions - set(region_fields)
    if missing:
        _reclassify_fail(f"region input lacks columns: {sorted(missing)}")
    required_contigs = {"MAG_ID", "contig", "length", "Taxonomic_shift", "Taxonomic_chimera_assessable"}
    missing = required_contigs - set(contig_fields)
    if missing:
        _reclassify_fail(f"contig input lacks columns: {sorted(missing)}")
    if set(mag_contig_fields) != set(contig_fields):
        missing_mc = required_contigs - set(mag_contig_fields)
        if missing_mc:
            _reclassify_fail(f"MAG-contig input lacks columns: {sorted(missing_mc)}")
    if "MAG_ID" not in mag_fields:
        _reclassify_fail("MAG input lacks MAG_ID")

    new_region_fields = ensure_fields(region_fields, [
        "Taxonomic_chimera", "High_confidence_taxonomic_chimera",
        "final_taxonomic_chimera_policy", "high_confidence_policy",
        "final_taxonomic_chimera_evidence", "final_taxonomic_chimera_status",
    ])

    classified_regions: List[Dict[str, str]] = []
    regions_by_contig: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    regions_by_mag: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for base in region_rows:
        row = dict(base)
        annotate_region_with_clipping(row, clipping)
        call, high, evidence, status = final_region_call(row)
        row["Taxonomic_chimera"] = call
        row["High_confidence_taxonomic_chimera"] = high
        row["final_taxonomic_chimera_policy"] = "taxonomic_shift_plus_validated_Clipping_Y"
        row["high_confidence_policy"] = "taxonomic_shift_plus_Stringent_clipping_Y"
        row["final_taxonomic_chimera_evidence"] = evidence
        row["final_taxonomic_chimera_status"] = status
        classified_regions.append(row)
        regions_by_contig[row["contig"]].append(row)
        mag = (row.get("MAG_ID") or "").strip()
        if mag not in NA_VALUES:
            regions_by_mag[mag].append(row)

    contig_add = [
        "n_taxonomic_shift_regions", "n_taxonomic_chimera_regions",
        "n_high_confidence_taxonomic_chimera_regions", "n_taxonomic_shift_only_regions",
        "taxonomic_chimera_positions_1based", "high_confidence_taxonomic_chimera_positions_1based",
        "taxonomic_shift_only_positions_1based", "final_taxonomic_chimera_policy", "high_confidence_policy",
    ]
    new_contig_fields = ensure_fields(contig_fields, contig_add)
    classified_contigs: List[Dict[str, str]] = []
    contig_index: Dict[str, Dict[str, str]] = {}
    for base in contig_rows:
        c = base["contig"]
        row = classify_contig(base, regions_by_contig.get(c, []))
        classified_contigs.append(row)
        contig_index[c] = row

    classified_mag_contigs: List[Dict[str, str]] = []
    mag_contigs_by_mag: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for base in mag_contig_rows:
        c = base["contig"]
        if c not in contig_index:
            _reclassify_fail(f"MAG contig absent from assembly-wide contig table: {c}")
        row = dict(contig_index[c])
        classified_mag_contigs.append(row)
        mag = row["MAG_ID"]
        if mag in NA_VALUES:
            _reclassify_fail(f"MAG-contig table has blank MAG_ID: {c}")
        mag_contigs_by_mag[mag].append(row)

    mag_add = [
        "n_taxonomic_chimera_contigs", "n_high_confidence_taxonomic_chimera_contigs",
        "n_taxonomic_shift_only_contigs", "n_taxonomic_shift_regions", "n_taxonomic_chimera_regions",
        "n_high_confidence_taxonomic_chimera_regions", "n_taxonomic_shift_only_regions",
        "taxonomic_chimera_contig_ids", "high_confidence_taxonomic_chimera_contig_ids",
        "taxonomic_shift_only_contig_ids", "High_confidence_taxonomic_chimera",
        "final_taxonomic_chimera_policy", "high_confidence_policy",
    ]
    new_mag_fields = ensure_fields(mag_fields, mag_add)
    classified_mags: List[Dict[str, str]] = []
    for base in mag_rows:
        mag = base["MAG_ID"]
        classified_mags.append(classify_mag(base, mag_contigs_by_mag.get(mag, []), regions_by_mag.get(mag, [])))

    out_dir = Path(args.out_dir)
    aux_dir = Path(args.aux_dir) if args.aux_dir else out_dir
    aux_dir.mkdir(parents=True, exist_ok=True)
    all_regions_path = out_dir / "chimera.region.tsv"
    contig_path = out_dir / "chimera.contig.tsv"
    mag_contig_path = aux_dir / "chimera.MAG_contig.tsv"
    mag_path = out_dir / "chimera.MAG.tsv"
    meta_path = aux_dir / "run_metadata.final.tsv"

    chim_regions = [r for r in classified_regions if r["Taxonomic_chimera"] == "Y"]
    high_regions = [r for r in classified_regions if r["High_confidence_taxonomic_chimera"] == "Y"]

    _reclassify_atomic_write(all_regions_path, new_region_fields, classified_regions)
    _reclassify_atomic_write(contig_path, new_contig_fields, classified_contigs)
    _reclassify_atomic_write(mag_contig_path, new_contig_fields, classified_mag_contigs)
    _reclassify_atomic_write(mag_path, new_mag_fields, classified_mags)

    mag_positive = [m for m in classified_mags if m["Taxonomic_chimera"] == "Y"]
    mag_high = [m for m in classified_mags if m["High_confidence_taxonomic_chimera"] == "Y"]
    with open(meta_path, "w", encoding="utf-8") as handle:
        handle.write("key\tvalue\n")
        handle.write("reclassification_policy\ttaxonomic_shift_plus_validated_Clipping_Y\n")
        handle.write("validated_clipping_definition\tClipping_Y_from_03_clipping_exact_position_validation\n")
        handle.write(f"clipping_source\t{args.clipping}\n")
        handle.write("high_confidence_policy\ttaxonomic_shift_plus_Stringent_clipping_Y\n")
        handle.write("primary_supp_structural_call_policy\taudit_only_not_final_classifier\n")
        handle.write(f"n_taxonomic_shift_regions\t{len(classified_regions)}\n")
        handle.write(f"n_taxonomic_chimera_regions\t{len(chim_regions)}\n")
        handle.write(f"n_high_confidence_taxonomic_chimera_regions\t{len(high_regions)}\n")
        handle.write(f"n_taxonomic_chimera_contigs\t{sum(r['Taxonomic_chimera']=='Y' for r in classified_contigs)}\n")
        handle.write(f"n_high_confidence_taxonomic_chimera_contigs\t{sum(r['High_confidence_taxonomic_chimera']=='Y' for r in classified_contigs)}\n")
        handle.write(f"n_taxonomic_chimera_MAGs\t{len(mag_positive)}\n")
        handle.write(f"n_high_confidence_taxonomic_chimera_MAGs\t{len(mag_high)}\n")

    print(f"Taxonomic-shift regions: {len(classified_regions)}")
    print(f"Taxonomic chimeras (validated Clipping=Y): {len(chim_regions)}")
    print(f"High-confidence taxonomic chimeras (Stringent_clipping=Y): {len(high_regions)}")
    print(f"Assembly contigs with taxonomic chimera: {sum(r['Taxonomic_chimera']=='Y' for r in classified_contigs)}")
    print(f"MAG-associated regions with taxonomic chimera: {sum(r['Taxonomic_chimera']=='Y' and r.get('MAG_ID','') not in NA_VALUES for r in classified_regions)}")
    print(f"MAGs with taxonomic chimera: {len(mag_positive)}")
    print(f"MAGs with high-confidence taxonomic chimera: {len(mag_high)}")
    print(f"Outputs written to: {out_dir}")


def cmd_reclassify(argv):
    """Final calls: taxonomic transition plus validated clipping"""
    sys.argv = [sys.argv[0] + " reclassify", *argv]
    return _reclassify_main()


COMMANDS = {
    'evaluate': cmd_evaluate,
    'validate': cmd_validate,
    'mag': cmd_mag,
    'reclassify': cmd_reclassify,
}

if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] not in COMMANDS:
        sys.exit(__doc__)
    result = COMMANDS[sys.argv[1]](sys.argv[2:])
    sys.exit(result if isinstance(result, int) else 0)
