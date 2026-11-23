#!/usr/bin/env python3
"""Assembly variant tables: per-assembler classification, merge across assemblers and MAG summaries.

Subcommands:
  classify           Classify anvi'o variable sites for one assembler
  merge              Merge the assembler results
  mag                Write the MAG-level variant tables
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sqlite3
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple


CLASS_UNSUPPORTED = "Unsupported_assembly_variant"

CLASS_MINORITY = "Minority_assembly_variant"

CLASS_NO_DISAGREEMENT = "No_assembly_allele_disagreement"

BASES = ("A", "C", "G", "T")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fai", required=True)
    split_group = parser.add_mutually_exclusive_group(required=True)
    split_group.add_argument("--contigs-db")
    split_group.add_argument(
        "--split-map",
        help="Two-column TSV with anvi'o split name and renamed parent contig; used for batch-local contigs DBs",
    )
    parser.add_argument("--snv-table", required=True)
    parser.add_argument("--indel-table", required=True)
    parser.add_argument("--detail-out", required=True)
    parser.add_argument("--contig-out", required=True)
    parser.add_argument("--summary-out", required=True)
    parser.add_argument("--mag-map", default="NA")
    parser.add_argument(
        "--name-map",
        default="NA",
        help="Two-column anvi'o contig rename report",
    )
    parser.add_argument(
        "--renamed-fai",
        default="NA",
        help="FAI for the anvi'o-compatible renamed FASTA",
    )
    return parser.parse_args()


def normalize(name: str) -> str:
    return name.strip().lstrip("#").lower().replace("-", "_").replace(" ", "_")


def normalized_row(row: Mapping[str, str]) -> Dict[str, str]:
    return {normalize(k): (v.strip() if v is not None else "") for k, v in row.items()}


def make_event_id(source_table: str, row: Mapping[str, str]) -> str:
    """Return a deterministic ID for one exported anvi'o event row."""
    payload = {"source_table": source_table, "row": dict(sorted(row.items()))}
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    return f"{source_table}:{digest}"


def require(row: Mapping[str, str], key: str, table: str) -> str:
    if key not in row:
        raise SystemExit(f"ERROR: {table} lacks required column: {key}")
    return row[key]


def as_int(value: str, label: str) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError) as exc:
        raise SystemExit(f"ERROR: invalid integer for {label}: {value!r}") from exc


def as_float(value: str, label: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise SystemExit(f"ERROR: invalid numeric value for {label}: {value!r}") from exc


def load_lengths(fai: str) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    with open(fai, encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            contig = fields[0]
            if contig in lengths:
                raise SystemExit(f"ERROR: duplicate contig in FAI: {contig}")
            lengths[contig] = int(fields[1])
    if not lengths:
        raise SystemExit(f"ERROR: no contigs found in FAI: {fai}")
    return lengths


def load_name_map(
    path: str,
    original_lengths: Mapping[str, int],
    renamed_lengths: Mapping[str, int],
) -> Dict[str, str]:
    """Return anvi'o-compatible contig name -> original contig name."""
    if path in {"", "NA", "None", None}:  # type: ignore[arg-type]
        return {}
    if not renamed_lengths:
        raise SystemExit(
            "ERROR: --renamed-fai is required when --name-map is provided"
        )

    raw_pairs: list[tuple[int, str, str]] = []
    with open(path, encoding="utf-8-sig") as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise SystemExit(
                    f"ERROR: malformed anvi'o name-map line {line_no}: {line.rstrip()!r}"
                )
            first, second = fields[0].strip(), fields[1].strip()
            if not first or not second:
                raise SystemExit(
                    f"ERROR: blank contig name in anvi'o name map line {line_no}"
                )
            raw_pairs.append((line_no, first, second))

    if not raw_pairs:
        raise SystemExit(f"ERROR: no contig renaming records found: {path}")

    def normalize_original(value: str) -> str:
        value = value.strip()
        if value.startswith(">"):
            value = value[1:]
        return value.split(None, 1)[0] if value else value

    new_old_valid = all(
        first in renamed_lengths and normalize_original(second) in original_lengths
        for _, first, second in raw_pairs
    )
    old_new_valid = all(
        normalize_original(first) in original_lengths and second in renamed_lengths
        for _, first, second in raw_pairs
    )

    if new_old_valid and not old_new_valid:
        orientation = "renamed_original"
    elif old_new_valid and not new_old_valid:
        orientation = "original_renamed"
    elif new_old_valid and old_new_valid:
        if all(first == second for _, first, second in raw_pairs):
            orientation = "renamed_original"
        else:
            raise SystemExit(
                "ERROR: anvi'o name-map orientation is ambiguous against the original "
                "and renamed FASTA indices"
            )
    else:
        line_no, first, second = raw_pairs[0]
        raise SystemExit(
            "ERROR: anvi'o name map does not match the original and renamed FASTA "
            f"indices; first record at line {line_no}: {first!r} -> {second!r}"
        )

    result: Dict[str, str] = {}
    originals_seen: set[str] = set()
    renamed_seen: set[str] = set()
    for line_no, first, second in raw_pairs:
        if orientation == "renamed_original":
            renamed, original = first, normalize_original(second)
        else:
            original, renamed = normalize_original(first), second

        if original_lengths[original] != renamed_lengths[renamed]:
            raise SystemExit(
                "ERROR: contig length differs between original and renamed FASTA for "
                f"{original!r} -> {renamed!r}: "
                f"{original_lengths[original]} != {renamed_lengths[renamed]}"
            )
        if renamed in result and result[renamed] != original:
            raise SystemExit(
                f"ERROR: duplicate renamed contig in anvi'o name map: {renamed}"
            )
        if original in originals_seen:
            raise SystemExit(
                f"ERROR: duplicate original contig in anvi'o name map: {original}"
            )
        result[renamed] = original
        originals_seen.add(original)
        renamed_seen.add(renamed)

    missing_original = set(original_lengths) - originals_seen
    if missing_original:
        examples = ", ".join(sorted(missing_original)[:5])
        raise SystemExit(
            "ERROR: anvi'o name map does not cover all original assembly contigs; "
            f"examples: {examples}"
        )
    missing_renamed = set(renamed_lengths) - renamed_seen
    if missing_renamed:
        examples = ", ".join(sorted(missing_renamed)[:5])
        raise SystemExit(
            "ERROR: anvi'o name map does not cover all renamed assembly contigs; "
            f"examples: {examples}"
        )
    return result


def load_split_map(
    contigs_db: str,
    lengths: Mapping[str, int],
    renamed_to_original: Mapping[str, str],
) -> Dict[str, str]:
    con = sqlite3.connect(contigs_db)
    try:
        rows = con.execute("SELECT split, parent FROM splits_basic_info").fetchall()
    except sqlite3.Error as exc:
        raise SystemExit(
            f"ERROR: could not read splits_basic_info from anvi'o contigs DB: {contigs_db}: {exc}"
        ) from exc
    finally:
        con.close()
    split_to_contig: Dict[str, str] = {}
    for split, parent in rows:
        parent = str(parent)
        original_parent = renamed_to_original.get(parent, parent)
        if original_parent not in lengths:
            raise SystemExit(
                f"ERROR: contigs DB parent absent from assembly FAI/name map: {parent}"
            )
        if split in split_to_contig and split_to_contig[split] != original_parent:
            raise SystemExit(f"ERROR: split maps to multiple contigs: {split}")
        split_to_contig[str(split)] = original_parent
    if not split_to_contig:
        raise SystemExit("ERROR: anvi'o contigs DB contains no split-to-contig records")
    return split_to_contig


def load_split_map_tsv(
    path: str,
    lengths: Mapping[str, int],
    renamed_to_original: Mapping[str, str],
) -> Dict[str, str]:
    split_to_contig: Dict[str, str] = {}
    with open(path, encoding="utf-8-sig") as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line_no == 1 and len(fields) >= 2 and fields[0].lower() in {"split", "split_name"}:
                continue
            if len(fields) < 2:
                raise SystemExit(f"ERROR: malformed split-map line {line_no}: {line.rstrip()!r}")
            split, parent = fields[0].strip(), fields[1].strip()
            if not split or not parent:
                raise SystemExit(f"ERROR: blank split/parent in split-map line {line_no}")
            original_parent = renamed_to_original.get(parent, parent)
            if original_parent not in lengths:
                raise SystemExit(
                    f"ERROR: split-map parent absent from assembly FAI/name map: {parent}"
                )
            previous = split_to_contig.get(split)
            if previous is not None and previous != original_parent:
                raise SystemExit(f"ERROR: split maps to multiple contigs: {split}")
            split_to_contig[split] = original_parent
    if not split_to_contig:
        raise SystemExit(f"ERROR: split-map contains no records: {path}")
    return split_to_contig


def load_mag_map(path: str, lengths: Mapping[str, int]) -> Dict[str, str]:
    if path in {"", "NA", "None", None}:  # type: ignore[arg-type]
        return {}
    result: Dict[str, str] = {}
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        aliases = {normalize(name): name for name in fields}
        mag_field = aliases.get("mag_id") or aliases.get("id")
        contig_field = aliases.get("contig") or aliases.get("original_contig")
        if not mag_field or not contig_field:
            raise SystemExit("ERROR: MAG_CONTIG_MAP needs MAG_ID and Contig columns")
        for row in reader:
            contig = row[contig_field].strip()
            mag = row[mag_field].strip()
            if not contig or not mag:
                raise SystemExit("ERROR: blank MAG_ID/Contig value in MAG_CONTIG_MAP")
            if contig not in lengths:
                continue
            if contig in result:
                raise SystemExit(f"ERROR: duplicate contig in MAG_CONTIG_MAP: {contig}")
            result[contig] = mag
    return result


def read_table(path: str) -> Iterable[Dict[str, str]]:
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise SystemExit(f"ERROR: missing header in table: {path}")
        for raw in reader:
            if raw is None:
                continue
            yield normalized_row(raw)


def resolve_contig_and_position(
    row: Mapping[str, str],
    split_to_contig: Mapping[str, str],
    lengths: Mapping[str, int],
    table: str,
    *,
    allow_pre_contig_start_indel: bool = False,
) -> tuple[str, int, bool]:
    """Resolve an anvi'o coordinate to a usable 1-based output position."""
    split = require(row, "split_name", table)
    if split not in split_to_contig:
        raise SystemExit(f"ERROR: {table} split absent from contigs DB: {split}")
    contig = split_to_contig[split]
    pos0 = as_int(require(row, "pos_in_contig", table), f"{table}.pos_in_contig")

    if allow_pre_contig_start_indel and pos0 == -1:
        return contig, 1, True

    pos1 = pos0 + 1
    if pos1 < 1 or pos1 > lengths[contig]:
        raise SystemExit(
            f"ERROR: {table} position outside contig: {contig}:{pos1}/{lengths[contig]} "
            f"(raw pos_in_contig={pos0})"
        )
    return contig, pos1, False


def clean_gene(value: str) -> str:
    value = value.strip()
    return value if value and value not in {"-1", "None", "NA"} else "NA"


def add_disagreement(
    records: list[dict[str, object]],
    counts: MutableMapping[str, Counter],
    genes: MutableMapping[str, MutableMapping[str, set[str]]],
    *,
    contig: str,
    pos1: int,
    variant_type: str,
    assembly_allele: str,
    alternative_allele: str,
    assembly_count: object,
    alternative_count: object,
    coverage: object,
    variant_class: str,
    gene: str,
    source_table: str,
    event_id: str,
) -> None:
    numeric_coverage = float(coverage) if isinstance(coverage, (int, float)) else 0.0
    alternative_fraction = (
        float(alternative_count) / numeric_coverage
        if numeric_coverage and isinstance(alternative_count, (int, float))
        else None
    )
    assembly_fraction = (
        float(assembly_count) / numeric_coverage
        if numeric_coverage and isinstance(assembly_count, (int, float))
        else None
    )
    records.append(
        {
            "contig": contig,
            "pos_1based": pos1,
            "variant_type": variant_type,
            "assembly_allele": assembly_allele,
            "dominant_alternative_allele": alternative_allele,
            "assembly_support_count": assembly_count,
            "alternative_support_count": alternative_count,
            "coverage": coverage,
            "assembly_support_fraction": "NA" if assembly_fraction is None else f"{assembly_fraction:.6f}",
            "alternative_support_fraction": "NA" if alternative_fraction is None else f"{alternative_fraction:.6f}",
            "assembly_support_class": variant_class,
            "corresponding_gene_call": gene,
            "source_table": source_table,
            "event_id": event_id,
        }
    )
    counts[contig][(variant_class, variant_type)] += 1
    if gene != "NA":
        genes[contig][variant_class].add(gene)


def analyze_snvs(
    path: str,
    split_to_contig: Mapping[str, str],
    lengths: Mapping[str, int],
    records: list[dict[str, object]],
    counts: MutableMapping[str, Counter],
    genes: MutableMapping[str, MutableMapping[str, set[str]]],
    seen_event_ids: set[str],
    examined_by_contig: MutableMapping[str, int],
    qualify_gene_ids: bool = False,
) -> tuple[int, int]:
    examined = 0
    exact_duplicates_skipped = 0
    for row in read_table(path):
        examined += 1
        contig, pos1, _ = resolve_contig_and_position(
            row, split_to_contig, lengths, "variable_nucleotides"
        )
        examined_by_contig[contig] += 1
        reference = require(row, "reference", "variable_nucleotides").upper()
        if reference not in BASES:
            continue
        base_counts = {
            base: as_int(require(row, base.lower(), "variable_nucleotides"), f"SNV.{base}")
            for base in BASES
        }
        coverage = sum(base_counts.values())
        if coverage <= 0:
            continue
        assembly_count = base_counts[reference]
        alt_base, alt_count = max(
            ((base, count) for base, count in base_counts.items() if base != reference),
            key=lambda item: (item[1], item[0]),
        )
        if assembly_count == 0 and alt_count > 0:
            variant_class = CLASS_UNSUPPORTED
        elif assembly_count > 0 and alt_count / coverage > 0.5:
            variant_class = CLASS_MINORITY
        else:
            continue
        gene = clean_gene(require(row, "corresponding_gene_call", "variable_nucleotides"))
        if qualify_gene_ids and gene != "NA":
            gene = f"{contig}:{gene}"
        event_id = make_event_id("variable_nucleotides", row)
        if event_id in seen_event_ids:
            exact_duplicates_skipped += 1
            continue
        seen_event_ids.add(event_id)
        add_disagreement(
            records,
            counts,
            genes,
            contig=contig,
            pos1=pos1,
            variant_type="SNV",
            assembly_allele=reference,
            alternative_allele=alt_base,
            assembly_count=assembly_count,
            alternative_count=alt_count,
            coverage=coverage,
            variant_class=variant_class,
            gene=gene,
            source_table="variable_nucleotides",
            event_id=event_id,
        )
    return examined, exact_duplicates_skipped


def analyze_indels(
    path: str,
    split_to_contig: Mapping[str, str],
    lengths: Mapping[str, int],
    records: list[dict[str, object]],
    counts: MutableMapping[str, Counter],
    genes: MutableMapping[str, MutableMapping[str, set[str]]],
    seen_event_ids: set[str],
    examined_by_contig: MutableMapping[str, int],
    qualify_gene_ids: bool = False,
) -> tuple[int, int]:
    examined = 0
    exact_duplicates_skipped = 0
    for row in read_table(path):
        examined += 1
        contig, pos1, _pre_contig_start = resolve_contig_and_position(
            row, split_to_contig, lengths, "indels",
            allow_pre_contig_start_indel=True,
        )
        examined_by_contig[contig] += 1
        event_count = as_int(require(row, "count", "indels"), "INDEL.count")
        coverage = as_float(require(row, "coverage", "indels"), "INDEL.coverage")
        if coverage <= 0 or event_count <= 0:
            continue

        if event_count == coverage:
            variant_class = CLASS_UNSUPPORTED
            assembly_count: object = 0
        elif event_count / coverage > 0.5:
            variant_class = CLASS_MINORITY
            assembly_count = coverage - event_count if event_count < coverage else "NA"
        else:
            continue
        indel_type = require(row, "type", "indels") or "INDEL"
        sequence = require(row, "sequence", "indels")
        reference = require(row, "reference", "indels")
        alt_label = f"{indel_type}:{sequence}" if sequence else indel_type
        gene = clean_gene(require(row, "corresponding_gene_call", "indels"))
        if qualify_gene_ids and gene != "NA":
            gene = f"{contig}:{gene}"
        event_id = make_event_id("indels", row)
        if event_id in seen_event_ids:
            exact_duplicates_skipped += 1
            continue
        seen_event_ids.add(event_id)
        add_disagreement(
            records,
            counts,
            genes,
            contig=contig,
            pos1=pos1,
            variant_type="INDEL",
            assembly_allele=reference or "reference_no_indel",
            alternative_allele=alt_label,
            assembly_count=assembly_count,
            alternative_count=event_count,
            coverage=coverage,
            variant_class=variant_class,
            gene=gene,
            source_table="indels",
            event_id=event_id,
        )
    return examined, exact_duplicates_skipped


def write_detail(path: str, records: list[dict[str, object]]) -> None:
    fields = [
        "contig",
        "pos_1based",
        "variant_type",
        "assembly_allele",
        "dominant_alternative_allele",
        "assembly_support_count",
        "alternative_support_count",
        "coverage",
        "assembly_support_fraction",
        "alternative_support_fraction",
        "assembly_support_class",
        "corresponding_gene_call",
        "source_table",
        "event_id",
    ]
    records.sort(
        key=lambda r: (
            str(r["contig"]),
            int(r["pos_1based"]),
            str(r["variant_type"]),
            str(r["event_id"]),
        )
    )
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def write_contig(
    path: str,
    lengths: Mapping[str, int],
    mag_map: Mapping[str, str],
    counts: Mapping[str, Counter],
    genes: Mapping[str, Mapping[str, set[str]]],
    examined_snvs_by_contig: Mapping[str, int],
    examined_indels_by_contig: Mapping[str, int],
) -> None:
    fields = [
        "MAG_ID",
        "contig",
        "length",
        "Assembly_variant",
        "anvio_variable_nucleotide_rows_examined",
        "anvio_indel_rows_examined",
        "minority_assembly_snv_count",
        "unsupported_assembly_snv_count",
        "minority_assembly_indel_count",
        "unsupported_assembly_indel_count",
        "minority_impacted_gene_count",
        "unsupported_impacted_gene_count",
        "assembly_variant_count",
        "assembly_allele_support_call",
    ]
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for contig in sorted(lengths):
            c = counts.get(contig, Counter())
            partial_snv = c[(CLASS_MINORITY, "SNV")]
            full_snv = c[(CLASS_UNSUPPORTED, "SNV")]
            partial_indel = c[(CLASS_MINORITY, "INDEL")]
            full_indel = c[(CLASS_UNSUPPORTED, "INDEL")]
            total = partial_snv + full_snv + partial_indel + full_indel
            if full_snv + full_indel:
                call = CLASS_UNSUPPORTED
            elif partial_snv + partial_indel:
                call = CLASS_MINORITY
            else:
                call = CLASS_NO_DISAGREEMENT
            writer.writerow(
                {
                    "MAG_ID": mag_map.get(contig, "NA"),
                    "contig": contig,
                    "length": lengths[contig],
                    "Assembly_variant": "Y" if total > 0 else "N",
                    "anvio_variable_nucleotide_rows_examined": examined_snvs_by_contig.get(contig, 0),
                    "anvio_indel_rows_examined": examined_indels_by_contig.get(contig, 0),
                    "minority_assembly_snv_count": partial_snv,
                    "unsupported_assembly_snv_count": full_snv,
                    "minority_assembly_indel_count": partial_indel,
                    "unsupported_assembly_indel_count": full_indel,
                    "minority_impacted_gene_count": len(genes.get(contig, {}).get(CLASS_MINORITY, set())),
                    "unsupported_impacted_gene_count": len(genes.get(contig, {}).get(CLASS_UNSUPPORTED, set())),
                    "assembly_variant_count": total,
                    "assembly_allele_support_call": call,
                }
            )


def write_summary(
    path: str,
    lengths: Mapping[str, int],
    counts: Mapping[str, Counter],
    genes: Mapping[str, Mapping[str, set[str]]],
    examined_snvs: int,
    examined_indels: int,
) -> None:
    total_bp = sum(lengths.values())
    total_mb = total_bp / 1_000_000.0
    aggregate = Counter()
    partial_genes: set[str] = set()
    full_genes: set[str] = set()
    for contig in lengths:
        aggregate.update(counts.get(contig, Counter()))
        partial_genes.update(genes.get(contig, {}).get(CLASS_MINORITY, set()))
        full_genes.update(genes.get(contig, {}).get(CLASS_UNSUPPORTED, set()))

    def value(variant_class: str, variant_type: str) -> int:
        return aggregate[(variant_class, variant_type)]

    fields = [
        "analyzed_contig_count",
        "assembly_length_bp",
        "anvio_variable_nucleotide_rows_examined",
        "anvio_indel_rows_examined",
        "minority_assembly_snv_count",
        "minority_assembly_snv_per_Mbp",
        "unsupported_assembly_snv_count",
        "unsupported_assembly_snv_per_Mbp",
        "minority_assembly_indel_count",
        "minority_assembly_indel_per_Mbp",
        "unsupported_assembly_indel_count",
        "unsupported_assembly_indel_per_Mbp",
        "minority_impacted_gene_count",
        "unsupported_impacted_gene_count",
        "contigs_with_assembly_allele_disagreement",
        "assembly_variant_count",
        "assembly_variant_per_Mbp",
    ]
    partial_snv = value(CLASS_MINORITY, "SNV")
    full_snv = value(CLASS_UNSUPPORTED, "SNV")
    partial_indel = value(CLASS_MINORITY, "INDEL")
    full_indel = value(CLASS_UNSUPPORTED, "INDEL")
    denom = total_mb if total_mb > 0 else 1.0
    total_variants = partial_snv + full_snv + partial_indel + full_indel
    contigs_with_variants = sum(
        1
        for contig in lengths
        if sum(counts.get(contig, Counter()).values()) > 0
    )
    row = {
        "analyzed_contig_count": len(lengths),
        "assembly_length_bp": total_bp,
        "anvio_variable_nucleotide_rows_examined": examined_snvs,
        "anvio_indel_rows_examined": examined_indels,
        "minority_assembly_snv_count": partial_snv,
        "minority_assembly_snv_per_Mbp": f"{partial_snv / denom:.6f}",
        "unsupported_assembly_snv_count": full_snv,
        "unsupported_assembly_snv_per_Mbp": f"{full_snv / denom:.6f}",
        "minority_assembly_indel_count": partial_indel,
        "minority_assembly_indel_per_Mbp": f"{partial_indel / denom:.6f}",
        "unsupported_assembly_indel_count": full_indel,
        "unsupported_assembly_indel_per_Mbp": f"{full_indel / denom:.6f}",
        "minority_impacted_gene_count": len(partial_genes),
        "unsupported_impacted_gene_count": len(full_genes),
        "contigs_with_assembly_allele_disagreement": contigs_with_variants,
        "assembly_variant_count": total_variants,
        "assembly_variant_per_Mbp": f"{total_variants / denom:.6f}",
    }
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)


def _classify_main() -> None:
    args = parse_args()
    required_paths = [args.fai, args.snv_table, args.indel_table]
    required_paths.append(args.split_map if args.split_map else args.contigs_db)
    if args.name_map not in {"", "NA", "None", None}:
        required_paths.append(args.renamed_fai)
    for path in required_paths:
        if not Path(path).is_file():
            raise SystemExit(f"ERROR: required input is missing: {path}")

    lengths = load_lengths(args.fai)
    renamed_lengths = (
        load_lengths(args.renamed_fai)
        if args.name_map not in {"", "NA", "None", None}
        else {}
    )
    renamed_to_original = load_name_map(
        args.name_map, lengths, renamed_lengths
    )
    if args.split_map:
        split_to_contig = load_split_map_tsv(args.split_map, lengths, renamed_to_original)
    else:
        split_to_contig = load_split_map(args.contigs_db, lengths, renamed_to_original)
    mag_map = load_mag_map(args.mag_map, lengths)
    records: list[dict[str, object]] = []
    counts: MutableMapping[str, Counter] = defaultdict(Counter)
    examined_snvs_by_contig: MutableMapping[str, int] = defaultdict(int)
    examined_indels_by_contig: MutableMapping[str, int] = defaultdict(int)
    genes: MutableMapping[str, MutableMapping[str, set[str]]] = defaultdict(
        lambda: defaultdict(set)
    )
    seen_event_ids: set[str] = set()
    examined_snvs, duplicate_snvs = analyze_snvs(
        args.snv_table, split_to_contig, lengths, records, counts, genes, seen_event_ids,
        examined_snvs_by_contig, qualify_gene_ids=bool(args.split_map),
    )
    examined_indels, duplicate_indels = analyze_indels(
        args.indel_table, split_to_contig, lengths, records, counts, genes, seen_event_ids,
        examined_indels_by_contig, qualify_gene_ids=bool(args.split_map),
    )
    write_detail(args.detail_out, records)
    write_contig(
        args.contig_out, lengths, mag_map, counts, genes,
        examined_snvs_by_contig, examined_indels_by_contig,
    )
    write_summary(
        args.summary_out, lengths, counts, genes, examined_snvs, examined_indels
    )
    print(f"anvi'o variable-nucleotide rows examined: {examined_snvs}")
    print(f"anvi'o INDEL rows examined: {examined_indels}")
    print(
        "exact duplicate exported event rows skipped: "
        f"{duplicate_snvs + duplicate_indels}"
    )
    print(f"assembly-allele disagreement rows written: {len(records)}")


def cmd_classify(argv):
    """Classify anvi'o variable sites for one assembler"""
    sys.argv = [sys.argv[0] + " classify", *argv]
    return _classify_main()


POSITIVE = {"Minority_assembly_variant", "Unsupported_assembly_variant"}

COUNT_FIELDS = (
    "minority_assembly_snv_count",
    "unsupported_assembly_snv_count",
    "minority_assembly_indel_count",
    "unsupported_assembly_indel_count",
)

EXAMINED_FIELDS = (
    "anvio_variable_nucleotide_rows_examined",
    "anvio_indel_rows_examined",
)


def _merge_parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", required=True)
    p.add_argument("--detail-out", required=True)
    p.add_argument("--contig-out", required=True)
    p.add_argument("--summary-out", required=True)
    p.add_argument("--assembler-summary-out", required=True)
    return p.parse_args()


def read_tsv(path: str) -> tuple[List[str], List[Dict[str, str]]]:
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise SystemExit(f"ERROR: TSV lacks a header: {path}")
        return list(reader.fieldnames), [dict(row) for row in reader]


def parse_int(value: str, label: str) -> int:
    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise SystemExit(f"ERROR: invalid integer for {label}: {value!r}") from exc
    if result < 0:
        raise SystemExit(f"ERROR: negative integer for {label}: {value!r}")
    return result


def atomic_write(path: str, fields: Iterable[str], rows: Iterable[Mapping[str, object]]) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    tmp = f"{path}.tmp.{os.getpid()}"
    field_list = list(fields)
    with open(tmp, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=field_list,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


def per_mbp(value: int, bp: int) -> str:
    return f"{value / (bp / 1_000_000.0):.6f}" if bp > 0 else "0.000000"


def summarize_contigs(
    rows: List[Dict[str, str]],
    detail_rows: List[Dict[str, str]],
) -> Dict[str, object]:
    total_bp = sum(parse_int(row["length"], "length") for row in rows)
    totals = {field: sum(parse_int(row[field], field) for row in rows) for field in COUNT_FIELDS}
    examined = {
        field: sum(parse_int(row[field], field) for row in rows) for field in EXAMINED_FIELDS
    }
    total_variants = sum(totals.values())
    assembly_variant_result = "Y" if total_variants > 0 else "N"
    contigs_with = sum((row.get("assembly_allele_support_call") or "") in POSITIVE for row in rows)
    contigs_minority = sum(
        (row.get("assembly_allele_support_call") or "") == "Minority_assembly_variant"
        for row in rows
    )
    contigs_unsupported = sum(
        (row.get("assembly_allele_support_call") or "") == "Unsupported_assembly_variant"
        for row in rows
    )

    minority_genes = set()
    unsupported_genes = set()
    for row in detail_rows:
        gene = (row.get("corresponding_gene_call") or "").strip()
        if gene in {"", "NA", "None", "-1"}:
            continue
        gene_key = ((row.get("Source_assembler") or "").strip(), gene)
        cls = (row.get("assembly_support_class") or "").strip()
        if cls == "Minority_assembly_variant":
            minority_genes.add(gene_key)
        elif cls == "Unsupported_assembly_variant":
            unsupported_genes.add(gene_key)

    return {
        "analyzed_contig_count": len(rows),
        "Assembly_variant": assembly_variant_result,
        "assembly_length_bp": total_bp,
        **examined,
        "minority_assembly_snv_count": totals["minority_assembly_snv_count"],
        "minority_assembly_snv_per_Mbp": per_mbp(totals["minority_assembly_snv_count"], total_bp),
        "unsupported_assembly_snv_count": totals["unsupported_assembly_snv_count"],
        "unsupported_assembly_snv_per_Mbp": per_mbp(totals["unsupported_assembly_snv_count"], total_bp),
        "minority_assembly_indel_count": totals["minority_assembly_indel_count"],
        "minority_assembly_indel_per_Mbp": per_mbp(totals["minority_assembly_indel_count"], total_bp),
        "unsupported_assembly_indel_count": totals["unsupported_assembly_indel_count"],
        "unsupported_assembly_indel_per_Mbp": per_mbp(totals["unsupported_assembly_indel_count"], total_bp),
        "minority_impacted_gene_count": len(minority_genes),
        "unsupported_impacted_gene_count": len(unsupported_genes),
        "contigs_with_assembly_allele_disagreement": contigs_with,
        "contigs_classified_as_minority_assembly_variant": contigs_minority,
        "contigs_classified_as_unsupported_assembly_variant": contigs_unsupported,
        "assembly_variant_count": total_variants,
        "assembly_variant_per_Mbp": per_mbp(total_variants, total_bp),
    }


def _merge_main() -> None:
    args = _merge_parse_args()
    manifest_fields, manifest_rows = read_tsv(args.manifest)
    required_manifest = {
        "Assembler",
        "Assembly_FASTA",
        "Detail_TSV",
        "Contig_TSV",
        "Summary_TSV",
        "SNV_raw_TSV",
        "INDEL_raw_TSV",
    }
    missing = sorted(required_manifest - set(manifest_fields))
    if missing:
        raise SystemExit("ERROR: per-assembler result manifest lacks columns: " + ", ".join(missing))
    if not manifest_rows:
        raise SystemExit("ERROR: per-assembler result manifest contains no rows")

    all_detail: List[Dict[str, str]] = []
    all_contigs: List[Dict[str, str]] = []
    detail_fields: List[str] | None = None
    contig_fields: List[str] | None = None
    seen_assemblers = set()
    seen_contigs: Dict[str, str] = {}
    seen_events = set()
    per_assembler: List[Dict[str, object]] = []

    for manifest_row in manifest_rows:
        assembler = (manifest_row.get("Assembler") or "").strip()
        if not assembler or assembler in seen_assemblers:
            raise SystemExit(f"ERROR: blank or duplicate assembler in manifest: {assembler!r}")
        seen_assemblers.add(assembler)
        detail_path = manifest_row["Detail_TSV"]
        contig_path = manifest_row["Contig_TSV"]
        summary_path = manifest_row["Summary_TSV"]
        for path in (detail_path, contig_path, summary_path):
            if not os.path.isfile(path) or os.path.getsize(path) == 0:
                raise SystemExit(f"ERROR: missing/empty per-assembler result for {assembler}: {path}")

        d_fields, d_rows = read_tsv(detail_path)
        c_fields, c_rows = read_tsv(contig_path)
        if detail_fields is None:
            detail_fields = d_fields
        elif d_fields != detail_fields:
            raise SystemExit(f"ERROR: detail schema differs for assembler {assembler}")
        if contig_fields is None:
            contig_fields = c_fields
        elif c_fields != contig_fields:
            raise SystemExit(f"ERROR: contig schema differs for assembler {assembler}")

        assembler_detail = []
        for row in d_rows:
            raw_event = (row.get("event_id") or "").strip()
            if not raw_event:
                raise SystemExit(f"ERROR: blank event_id in {detail_path}")
            event_id = f"{assembler}::{raw_event}"
            if event_id in seen_events:
                raise SystemExit(f"ERROR: duplicate merged event_id: {event_id}")
            seen_events.add(event_id)
            merged = {"Source_assembler": assembler, **row, "event_id": event_id}
            all_detail.append(merged)
            assembler_detail.append(merged)

        for row in c_rows:
            contig = (row.get("contig") or "").strip()
            if not contig:
                raise SystemExit(f"ERROR: blank contig in {contig_path}")
            previous = seen_contigs.get(contig)
            if previous is not None:
                raise SystemExit(
                    "ERROR: a contig ID occurs in more than one assembly and cannot be "
                    f"merged unambiguously: {contig} ({previous}, {assembler})"
                )
            seen_contigs[contig] = assembler
            all_contigs.append({"Source_assembler": assembler, **row})

        stats = summarize_contigs(
            [{"Source_assembler": assembler, **row} for row in c_rows], assembler_detail
        )
        mag_rows = [row for row in c_rows if (row.get("MAG_ID") or "").strip() not in {"", "NA", "None"}]
        mag_bp = sum(parse_int(row["length"], "length") for row in mag_rows)
        mag_variants = sum(
            sum(parse_int(row[field], field) for field in COUNT_FIELDS) for row in mag_rows
        )
        mag_contigs_with = sum(
            (row.get("assembly_allele_support_call") or "") in POSITIVE for row in mag_rows
        )
        per_assembler.append(
            {
                "Assembler": assembler,
                "Assembly_FASTA": os.path.realpath(manifest_row["Assembly_FASTA"]),
                **stats,
                "MAG_assigned_contig_count": len(mag_rows),
                "MAG_assigned_bp": mag_bp,
                "MAG_assigned_contigs_with_assembly_allele_disagreement": mag_contigs_with,
                "MAG_assigned_assembly_variant_count": mag_variants,
                "MAG_assigned_assembly_variant_per_Mbp": per_mbp(mag_variants, mag_bp),
            }
        )

    assert detail_fields is not None and contig_fields is not None
    all_detail.sort(
        key=lambda row: (
            row["Source_assembler"],
            row["contig"],
            int(row["pos_1based"]),
            row["variant_type"],
            row["event_id"],
        )
    )
    all_contigs.sort(key=lambda row: (row["Source_assembler"], row["contig"]))
    per_assembler.sort(key=lambda row: str(row["Assembler"]))

    atomic_write(args.detail_out, ["Source_assembler", *detail_fields], all_detail)
    atomic_write(args.contig_out, ["Source_assembler", *contig_fields], all_contigs)

    overall = summarize_contigs(all_contigs, all_detail)
    summary_fields = list(overall)
    atomic_write(args.summary_out, summary_fields, [overall])

    assembler_fields = [
        "Assembler",
        "Assembly_FASTA",
        *summary_fields,
        "MAG_assigned_contig_count",
        "MAG_assigned_bp",
        "MAG_assigned_contigs_with_assembly_allele_disagreement",
        "MAG_assigned_assembly_variant_count",
        "MAG_assigned_assembly_variant_per_Mbp",
    ]
    atomic_write(args.assembler_summary_out, assembler_fields, per_assembler)

    print(f"Assemblers merged: {len(per_assembler)}")
    print(f"Assembly-wide contigs: {len(all_contigs)}")
    print(f"Assembly-wide disagreement events: {len(all_detail)}")
    print(f"Assembly-wide sequence length: {overall['assembly_length_bp']} bp")


def cmd_merge(argv):
    """Merge the assembler results"""
    sys.argv = [sys.argv[0] + " merge", *argv]
    return _merge_main()


VALID_CLASSES = {
    CLASS_UNSUPPORTED,
    CLASS_MINORITY,
    CLASS_NO_DISAGREEMENT,
}

POSITIVE_CLASSES = {CLASS_UNSUPPORTED, CLASS_MINORITY}

VARIANT_TYPES = {"SNV", "INDEL"}

CONTIG_REQUIRED = (
    "MAG_ID",
    "contig",
    "length",
    "Assembly_variant",
    "anvio_variable_nucleotide_rows_examined",
    "anvio_indel_rows_examined",
    "minority_assembly_snv_count",
    "unsupported_assembly_snv_count",
    "minority_assembly_indel_count",
    "unsupported_assembly_indel_count",
    "minority_impacted_gene_count",
    "unsupported_impacted_gene_count",
    "assembly_variant_count",
    "assembly_allele_support_call",
)

DETAIL_REQUIRED = (
    "contig",
    "pos_1based",
    "variant_type",
    "assembly_allele",
    "dominant_alternative_allele",
    "assembly_support_count",
    "alternative_support_count",
    "coverage",
    "assembly_support_fraction",
    "alternative_support_fraction",
    "assembly_support_class",
    "corresponding_gene_call",
    "source_table",
)

MAG_SUMMARY_FIELDS = (
    "MAG_ID",
    "Assembly_variant",
    "Minority_assembly_variant",
    "Unsupported_assembly_variant",
    "n_contigs",
    "MAG_length_bp",
    "n_contigs_with_assembly_allele_disagreement",
    "n_contigs_classified_as_minority_assembly_variant",
    "n_contigs_classified_as_unsupported_assembly_variant",
    "minority_assembly_snv_count",
    "minority_assembly_snv_per_Mbp",
    "unsupported_assembly_snv_count",
    "unsupported_assembly_snv_per_Mbp",
    "minority_assembly_indel_count",
    "minority_assembly_indel_per_Mbp",
    "unsupported_assembly_indel_count",
    "unsupported_assembly_indel_per_Mbp",
    "minority_impacted_gene_count",
    "unsupported_impacted_gene_count",
    "assembly_variant_count",
    "assembly_variant_per_Mbp",
    "MAG_assembly_allele_support_call",
)


def _mag_parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create MAG-associated and MAG-level assembly-allele support tables "
            "from existing assembly-variant outputs."
        )
    )
    parser.add_argument("--contig-summary", required=True)
    parser.add_argument("--detail", required=True)
    parser.add_argument("--mag-contig-map", required=True)
    parser.add_argument("--site-output", required=True)
    parser.add_argument("--summary-output", required=True)
    parser.add_argument(
        "--contig-output",
        default="NA",
        help="Optional MAG-only contig table retaining the input contig-summary columns",
    )
    parser.add_argument(
        "--overall-summary-output",
        default="NA",
        help="Optional one-row summary restricted to final-MAG contigs",
    )
    return parser.parse_args()


def derive_event_id(row: Mapping[str, str]) -> str:
    """Create a stable fallback event ID for legacy disagreement tables."""
    payload = {key: value for key, value in sorted(row.items()) if key != "event_id"}
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    return f"legacy:{digest}"


def require_nonempty(path: str, label: str) -> None:
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise SystemExit(f"ERROR: {label} not found or empty: {path}")


def resolve_column(
    fieldnames: Sequence[str], candidates: Iterable[str], label: str
) -> str:
    normalized = {normalize(name): name for name in fieldnames}
    for candidate in candidates:
        if candidate in normalized:
            return normalized[candidate]
    raise SystemExit(
        f"ERROR: MAG_CONTIG_MAP needs a {label} column; observed columns: {fieldnames}"
    )


def parse_nonnegative_int(value: str, field: str, line_number: int) -> int:
    try:
        result = int(value)
    except (TypeError, ValueError) as error:
        raise SystemExit(
            f"ERROR: invalid integer in {field} at line {line_number}: {value!r}"
        ) from error
    if result < 0:
        raise SystemExit(
            f"ERROR: negative integer in {field} at line {line_number}: {value!r}"
        )
    return result


def _mag_load_mag_map(
    path: str,
) -> Tuple[Dict[str, str], Dict[str, List[str]]]:
    contig_to_mag: Dict[str, str] = {}
    mag_to_contigs: MutableMapping[str, List[str]] = defaultdict(list)

    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        mag_column = resolve_column(
            fields,
            ("mag_id", "entity_id", "genome_id", "bin_id", "genome", "bin", "id"),
            "MAG ID",
        )
        contig_column = resolve_column(fields, ("contig", "original_contig"), "Contig")
        for line_number, row in enumerate(reader, start=2):
            mag_id = (row.get(mag_column) or "").strip()
            contig = (row.get(contig_column) or "").strip()
            if not mag_id or not contig:
                raise SystemExit(
                    f"ERROR: blank MAG/contig in MAG_CONTIG_MAP at line {line_number}"
                )
            previous = contig_to_mag.get(contig)
            if previous is not None:
                if previous == mag_id:
                    raise SystemExit(
                        f"ERROR: duplicate MAG-map assignment at line {line_number}: "
                        f"{mag_id}/{contig}"
                    )
                raise SystemExit(
                    f"ERROR: contig assigned to multiple MAGs: {contig} "
                    f"({previous}, {mag_id})"
                )
            contig_to_mag[contig] = mag_id
            mag_to_contigs[mag_id].append(contig)

    if not mag_to_contigs:
        raise SystemExit(f"ERROR: MAG_CONTIG_MAP contains no assignments: {path}")
    for contigs in mag_to_contigs.values():
        contigs.sort()
    return contig_to_mag, dict(mag_to_contigs)


def expected_call(
    minority_snv: int,
    unsupported_snv: int,
    minority_indel: int,
    unsupported_indel: int,
) -> str:
    if unsupported_snv + unsupported_indel > 0:
        return CLASS_UNSUPPORTED
    if minority_snv + minority_indel > 0:
        return CLASS_MINORITY
    return CLASS_NO_DISAGREEMENT


def load_contig_summary(
    path: str,
    contig_to_mag: Mapping[str, str],
) -> Dict[str, Dict[str, object]]:
    rows: Dict[str, Dict[str, object]] = {}
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        missing = [field for field in CONTIG_REQUIRED if field not in fields]
        if missing:
            raise SystemExit(f"ERROR: contig summary lacks columns: {missing}")

        for line_number, raw in enumerate(reader, start=2):
            contig = (raw.get("contig") or "").strip()
            if not contig:
                raise SystemExit(f"ERROR: blank contig in contig summary at line {line_number}")
            if contig in rows:
                raise SystemExit(f"ERROR: duplicate contig in contig summary: {contig}")
            if contig not in contig_to_mag:
                continue
            mag_id = (raw.get("MAG_ID") or "").strip()
            expected_mag = contig_to_mag[contig]
            if mag_id != expected_mag:
                raise SystemExit(
                    f"ERROR: MAG_ID mismatch for {contig}: contig summary={mag_id!r}, "
                    f"MAG_CONTIG_MAP={expected_mag!r}"
                )

            numeric_fields = (
                "length",
                "anvio_variable_nucleotide_rows_examined",
                "anvio_indel_rows_examined",
                "minority_assembly_snv_count",
                "unsupported_assembly_snv_count",
                "minority_assembly_indel_count",
                "unsupported_assembly_indel_count",
                "minority_impacted_gene_count",
                "unsupported_impacted_gene_count",
                "assembly_variant_count",
            )
            parsed = {
                field: parse_nonnegative_int(raw.get(field, ""), field, line_number)
                for field in numeric_fields
            }
            if parsed["length"] < 1:
                raise SystemExit(f"ERROR: non-positive contig length for {contig}")
            count_sum = (
                parsed["minority_assembly_snv_count"]
                + parsed["unsupported_assembly_snv_count"]
                + parsed["minority_assembly_indel_count"]
                + parsed["unsupported_assembly_indel_count"]
            )
            if parsed["assembly_variant_count"] != count_sum:
                raise SystemExit(
                    f"ERROR: assembly_variant_count does not equal component counts for {contig}: "
                    f"{parsed['assembly_variant_count']} != {count_sum}"
                )
            call = (raw.get("assembly_allele_support_call") or "").strip()
            if call not in VALID_CLASSES:
                raise SystemExit(
                    f"ERROR: unrecognized assembly_allele_support_call for {contig}: {call}"
                )
            calculated_call = expected_call(
                parsed["minority_assembly_snv_count"],
                parsed["unsupported_assembly_snv_count"],
                parsed["minority_assembly_indel_count"],
                parsed["unsupported_assembly_indel_count"],
            )
            if call != calculated_call:
                raise SystemExit(
                    f"ERROR: contig call/count mismatch for {contig}: {call} != {calculated_call}"
                )
            result = (raw.get("Assembly_variant") or "N").strip().upper()
            expected_result = "Y" if parsed["assembly_variant_count"] > 0 else "N"
            if result != expected_result:
                raise SystemExit(
                    f"ERROR: disagreement/result mismatch for {contig}: {result}; "
                    f"expected result={expected_result}"
                )
            rows[contig] = {
                "MAG_ID": mag_id,
                "call": call,
                "result": result,
                "raw": {field: (raw.get(field) or "").strip() for field in fields},
                **parsed,
            }

    missing_contigs = set(contig_to_mag) - set(rows)
    if missing_contigs:
        examples = ", ".join(sorted(missing_contigs)[:5])
        raise SystemExit(
            "ERROR: contig summary does not contain every MAG-map contig; "
            f"examples: {examples}"
        )
    return rows


def load_detail(
    path: str,
    contig_to_mag: Mapping[str, str],
) -> Tuple[
    List[Dict[str, str]],
    List[str],
    Mapping[str, Counter],
    Mapping[str, Mapping[str, set[str]]],
]:
    output_rows: List[Dict[str, str]] = []
    counts: MutableMapping[str, Counter] = defaultdict(Counter)
    genes: MutableMapping[str, MutableMapping[str, set[str]]] = defaultdict(
        lambda: defaultdict(set)
    )
    repeated_locus_keys = 0
    seen_locus_keys: set[tuple[str, str, str, str]] = set()
    seen_event_ids: set[str] = set()

    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        missing = [field for field in DETAIL_REQUIRED if field not in fields]
        if missing:
            raise SystemExit(f"ERROR: disagreement detail table lacks columns: {missing}")
        if "MAG_ID" in fields:
            raise SystemExit(
                "ERROR: input detail table already contains MAG_ID; provide the original "
                "assembly_allele_disagreements.tsv"
            )
        output_fields = list(fields)
        if "event_id" not in output_fields:
            output_fields.append("event_id")

        for line_number, raw in enumerate(reader, start=2):
            clean = {field: (raw.get(field) or "").strip() for field in fields}
            contig = clean["contig"]
            if not contig:
                raise SystemExit(f"ERROR: blank contig in detail table at line {line_number}")
            if contig not in contig_to_mag:
                continue
            position = parse_nonnegative_int(clean["pos_1based"], "pos_1based", line_number)
            if position < 1:
                raise SystemExit(
                    f"ERROR: non-positive position in detail table at line {line_number}"
                )
            variant_type = clean["variant_type"]
            variant_class = clean["assembly_support_class"]
            source_table = clean["source_table"]
            if variant_type not in VARIANT_TYPES:
                raise SystemExit(
                    f"ERROR: unrecognized variant_type at line {line_number}: {variant_type}"
                )
            if variant_class not in POSITIVE_CLASSES:
                raise SystemExit(
                    f"ERROR: detail row has a non-disagreement class at line {line_number}: "
                    f"{variant_class}"
                )
            event_id = clean.get("event_id", "") or derive_event_id(clean)
            clean["event_id"] = event_id
            if event_id in seen_event_ids:
                raise SystemExit(
                    "ERROR: exact duplicate disagreement event_id at line "
                    f"{line_number}: {event_id}. Regenerate the disagreement table "
                    "with the updated variants_summarize_assembly_variants.py parser."
                )
            seen_event_ids.add(event_id)

            locus_key = (contig, clean["pos_1based"], variant_type, source_table)
            if locus_key in seen_locus_keys:
                repeated_locus_keys += 1
            else:
                seen_locus_keys.add(locus_key)

            counts[contig][(variant_class, variant_type)] += 1
            gene = clean["corresponding_gene_call"]
            if gene and gene not in {"NA", "None", "-1"}:
                source_assembler = clean.get("Source_assembler", "")
                gene_key = f"{source_assembler}::{gene}" if source_assembler else gene
                genes[contig][variant_class].add(gene_key)
            output_rows.append({"MAG_ID": contig_to_mag[contig], **clean})

    if repeated_locus_keys:
        print(
            "Same-position distinct disagreement events retained: "
            f"{repeated_locus_keys}"
        )
    return output_rows, output_fields, counts, genes


def validate_detail_against_contigs(
    contig_rows: Mapping[str, Mapping[str, object]],
    detail_counts: Mapping[str, Counter],
    detail_genes: Mapping[str, Mapping[str, set[str]]],
) -> None:
    count_fields = {
        (CLASS_MINORITY, "SNV"): "minority_assembly_snv_count",
        (CLASS_UNSUPPORTED, "SNV"): "unsupported_assembly_snv_count",
        (CLASS_MINORITY, "INDEL"): "minority_assembly_indel_count",
        (CLASS_UNSUPPORTED, "INDEL"): "unsupported_assembly_indel_count",
    }
    for contig, row in contig_rows.items():
        for key, field in count_fields.items():
            observed = detail_counts.get(contig, Counter())[key]
            expected = int(row[field])
            if observed != expected:
                raise SystemExit(
                    f"ERROR: detail/contig count mismatch for {contig}/{field}: "
                    f"{observed} != {expected}"
                )
        for variant_class, field in (
            (CLASS_MINORITY, "minority_impacted_gene_count"),
            (CLASS_UNSUPPORTED, "unsupported_impacted_gene_count"),
        ):
            observed = len(detail_genes.get(contig, {}).get(variant_class, set()))
            expected = int(row[field])
            if observed != expected:
                raise SystemExit(
                    f"ERROR: detail/contig gene-count mismatch for {contig}/{field}: "
                    f"{observed} != {expected}"
                )


def write_site_output(path: str, rows: List[Dict[str, str]], input_fields: List[str]) -> None:
    fields = ["MAG_ID", *input_fields]
    rows.sort(
        key=lambda row: (
            row["MAG_ID"],
            row["contig"],
            int(row["pos_1based"]),
            row["variant_type"],
            row["source_table"],
            row["dominant_alternative_allele"],
            row["assembly_support_class"],
            row["event_id"],
        )
    )
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_mag_summary(
    path: str,
    mag_to_contigs: Mapping[str, List[str]],
    contig_rows: Mapping[str, Mapping[str, object]],
    detail_genes: Mapping[str, Mapping[str, set[str]]],
) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=MAG_SUMMARY_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()

        for mag_id in sorted(mag_to_contigs):
            contigs = mag_to_contigs[mag_id]
            mag_length = sum(int(contig_rows[contig]["length"]) for contig in contigs)
            denom_mbp = mag_length / 1_000_000.0

            aggregate = Counter()
            minority_genes: set[str] = set()
            unsupported_genes: set[str] = set()
            calls = Counter()
            for contig in contigs:
                row = contig_rows[contig]
                calls[str(row["call"])] += 1
                aggregate["minority_snv"] += int(row["minority_assembly_snv_count"])
                aggregate["unsupported_snv"] += int(row["unsupported_assembly_snv_count"])
                aggregate["minority_indel"] += int(row["minority_assembly_indel_count"])
                aggregate["unsupported_indel"] += int(row["unsupported_assembly_indel_count"])
                minority_genes.update(
                    detail_genes.get(contig, {}).get(CLASS_MINORITY, set())
                )
                unsupported_genes.update(
                    detail_genes.get(contig, {}).get(CLASS_UNSUPPORTED, set())
                )

            total = sum(aggregate.values())
            mag_result = "Y" if total > 0 else "N"
            mag_call = expected_call(
                aggregate["minority_snv"],
                aggregate["unsupported_snv"],
                aggregate["minority_indel"],
                aggregate["unsupported_indel"],
            )

            def per_mbp(value: int) -> str:
                return f"{value / denom_mbp:.6f}" if denom_mbp > 0 else "0.000000"

            writer.writerow(
                {
                    "MAG_ID": mag_id,
                    "Assembly_variant": mag_result,
                    "Minority_assembly_variant": (
                        "Y" if aggregate["minority_snv"] + aggregate["minority_indel"] > 0 else "N"
                    ),
                    "Unsupported_assembly_variant": (
                        "Y" if aggregate["unsupported_snv"] + aggregate["unsupported_indel"] > 0 else "N"
                    ),
                    "n_contigs": len(contigs),
                    "MAG_length_bp": mag_length,
                    "n_contigs_with_assembly_allele_disagreement": (
                        calls[CLASS_MINORITY] + calls[CLASS_UNSUPPORTED]
                    ),
                    "n_contigs_classified_as_minority_assembly_variant": calls[CLASS_MINORITY],
                    "n_contigs_classified_as_unsupported_assembly_variant": calls[CLASS_UNSUPPORTED],
                    "minority_assembly_snv_count": aggregate["minority_snv"],
                    "minority_assembly_snv_per_Mbp": per_mbp(aggregate["minority_snv"]),
                    "unsupported_assembly_snv_count": aggregate["unsupported_snv"],
                    "unsupported_assembly_snv_per_Mbp": per_mbp(
                        aggregate["unsupported_snv"]
                    ),
                    "minority_assembly_indel_count": aggregate["minority_indel"],
                    "minority_assembly_indel_per_Mbp": per_mbp(
                        aggregate["minority_indel"]
                    ),
                    "unsupported_assembly_indel_count": aggregate["unsupported_indel"],
                    "unsupported_assembly_indel_per_Mbp": per_mbp(
                        aggregate["unsupported_indel"]
                    ),
                    "minority_impacted_gene_count": len(minority_genes),
                    "unsupported_impacted_gene_count": len(unsupported_genes),
                    "assembly_variant_count": total,
                    "assembly_variant_per_Mbp": per_mbp(total),
                    "MAG_assembly_allele_support_call": mag_call,
                }
            )


def write_mag_contig_subset(
    path: str,
    contig_rows: Mapping[str, Mapping[str, object]],
) -> None:
    rows = [dict(contig_rows[contig]["raw"]) for contig in sorted(contig_rows)]
    if not rows:
        raise SystemExit("ERROR: no MAG contigs are available for MAG-only contig output")
    fields = list(rows[0])
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_mag_overall_summary(
    path: str,
    contig_rows: Mapping[str, Mapping[str, object]],
    detail_genes: Mapping[str, Mapping[str, set[str]]],
) -> None:
    total_bp = sum(int(row["length"]) for row in contig_rows.values())
    denom_mbp = total_bp / 1_000_000.0
    minority_snv = sum(int(row["minority_assembly_snv_count"]) for row in contig_rows.values())
    unsupported_snv = sum(int(row["unsupported_assembly_snv_count"]) for row in contig_rows.values())
    minority_indel = sum(int(row["minority_assembly_indel_count"]) for row in contig_rows.values())
    unsupported_indel = sum(int(row["unsupported_assembly_indel_count"]) for row in contig_rows.values())
    examined_snvs = sum(int(row["anvio_variable_nucleotide_rows_examined"]) for row in contig_rows.values())
    examined_indels = sum(int(row["anvio_indel_rows_examined"]) for row in contig_rows.values())
    total_variants = minority_snv + unsupported_snv + minority_indel + unsupported_indel
    minority_genes: set[str] = set()
    unsupported_genes: set[str] = set()
    for contig in contig_rows:
        minority_genes.update(detail_genes.get(contig, {}).get(CLASS_MINORITY, set()))
        unsupported_genes.update(detail_genes.get(contig, {}).get(CLASS_UNSUPPORTED, set()))
    contigs_with_variants = sum(
        str(row["call"]) in POSITIVE_CLASSES for row in contig_rows.values()
    )

    def per_mbp(value: int) -> str:
        return f"{value / denom_mbp:.6f}" if denom_mbp > 0 else "0.000000"

    fields = [
        "analyzed_contig_count",
        "assembly_length_bp",
        "anvio_variable_nucleotide_rows_examined",
        "anvio_indel_rows_examined",
        "minority_assembly_snv_count",
        "minority_assembly_snv_per_Mbp",
        "unsupported_assembly_snv_count",
        "unsupported_assembly_snv_per_Mbp",
        "minority_assembly_indel_count",
        "minority_assembly_indel_per_Mbp",
        "unsupported_assembly_indel_count",
        "unsupported_assembly_indel_per_Mbp",
        "minority_impacted_gene_count",
        "unsupported_impacted_gene_count",
        "contigs_with_assembly_allele_disagreement",
        "assembly_variant_count",
        "assembly_variant_per_Mbp",
    ]
    row = {
        "analyzed_contig_count": len(contig_rows),
        "assembly_length_bp": total_bp,
        "anvio_variable_nucleotide_rows_examined": examined_snvs,
        "anvio_indel_rows_examined": examined_indels,
        "minority_assembly_snv_count": minority_snv,
        "minority_assembly_snv_per_Mbp": per_mbp(minority_snv),
        "unsupported_assembly_snv_count": unsupported_snv,
        "unsupported_assembly_snv_per_Mbp": per_mbp(unsupported_snv),
        "minority_assembly_indel_count": minority_indel,
        "minority_assembly_indel_per_Mbp": per_mbp(minority_indel),
        "unsupported_assembly_indel_count": unsupported_indel,
        "unsupported_assembly_indel_per_Mbp": per_mbp(unsupported_indel),
        "minority_impacted_gene_count": len(minority_genes),
        "unsupported_impacted_gene_count": len(unsupported_genes),
        "contigs_with_assembly_allele_disagreement": contigs_with_variants,
        "assembly_variant_count": total_variants,
        "assembly_variant_per_Mbp": per_mbp(total_variants),
    }
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)


def _mag_main() -> None:
    args = _mag_parse_args()
    for path, label in (
        (args.contig_summary, "contig assembly-allele summary"),
        (args.detail, "assembly-allele disagreement detail"),
        (args.mag_contig_map, "MAG_CONTIG_MAP"),
    ):
        require_nonempty(path, label)

    outputs = [args.site_output, args.summary_output]
    if args.contig_output not in {"", "NA", "None", None}:
        outputs.append(args.contig_output)
    if args.overall_summary_output not in {"", "NA", "None", None}:
        outputs.append(args.overall_summary_output)
    for output in outputs:
        parent = Path(output).resolve().parent
        parent.mkdir(parents=True, exist_ok=True)

    contig_to_mag, mag_to_contigs = _mag_load_mag_map(args.mag_contig_map)
    contig_rows = load_contig_summary(args.contig_summary, contig_to_mag)
    site_rows, detail_fields, detail_counts, detail_genes = load_detail(
        args.detail, contig_to_mag
    )
    validate_detail_against_contigs(contig_rows, detail_counts, detail_genes)
    write_site_output(args.site_output, site_rows, detail_fields)
    write_mag_summary(args.summary_output, mag_to_contigs, contig_rows, detail_genes)
    if args.contig_output not in {"", "NA", "None", None}:
        write_mag_contig_subset(args.contig_output, contig_rows)
    if args.overall_summary_output not in {"", "NA", "None", None}:
        write_mag_overall_summary(args.overall_summary_output, contig_rows, detail_genes)

    n_positive_mags = 0
    with open(args.summary_output, encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["MAG_assembly_allele_support_call"] in POSITIVE_CLASSES:
                n_positive_mags += 1
    print(f"MAG-associated disagreement rows written: {len(site_rows)}")
    print(f"MAG summary rows written: {len(mag_to_contigs)}")
    print(f"MAGs with assembly-allele disagreement: {n_positive_mags}")


def cmd_mag(argv):
    """Write the MAG-level variant tables"""
    sys.argv = [sys.argv[0] + " mag", *argv]
    return _mag_main()


COMMANDS = {
    'classify': cmd_classify,
    'merge': cmd_merge,
    'mag': cmd_mag,
}

if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] not in COMMANDS:
        sys.exit(__doc__)
    result = COMMANDS[sys.argv[1]](sys.argv[2:])
    sys.exit(result if isinstance(result, int) else 0)
