#!/usr/bin/env python3
"""Summarize site-wise assembly-allele support from anvi'o tables.

The Trigodet et al. workflow exports anvi'o `variable_nucleotides` and
`indels` profile tables. For SNVs, an assembly allele is classified as a minority
assembly variant when the most abundant non-assembly nucleotide accounts for
more than half of nucleotide coverage while the assembly nucleotide retains
support, and as an unsupported assembly variant when the assembly nucleotide
has no alignment support. For indels, the event is classified as a minority
assembly variant when it accounts for more than half but less than all coverage,
and as unsupported when it accounts for all coverage.

These thresholds follow the >50% and all-coverage criteria reported in Fig. 1
of Trigodet et al. Ties and non-majority alternative nucleotides are not classified.
Classification is site-level. No separate contig-level evaluability state is inferred
from the presence of one or more mapped HiFi alignments.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, Mapping, MutableMapping

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
    """Return a deterministic ID for one exported anvi'o event row.

    anvi'o may report multiple INDEL events at the same contig coordinate.
    Consequently, contig/position/type is not a unique event key.  Hashing the
    complete normalized source row preserves distinct same-position events,
    including deletion rows whose exported ``sequence`` field is empty.
    """
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
    """Return anvi'o-compatible contig name -> original contig name.

    Rename-report orientation is detected from the original and renamed FASTA
    indices rather than assumed. Original entries may be full FASTA deflines,
    so they are normalized to the first whitespace-delimited token before
    comparison with the source FAI.
    """
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
    """Resolve an anvi'o coordinate to a usable 1-based output position.

    anvi'o INDEL tables can represent an event anchored immediately before the
    first reference base with ``pos_in_contig=-1``.  That boundary event is
    valid for assembly-variant counting but has no ordinary 1-based reference
    base coordinate.  Retain it at display/join coordinate 1 and return a flag
    so downstream validators can keep the event non-evaluable rather than
    interrogating base 1 as though it were the same boundary.
    """
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


def main() -> None:
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


if __name__ == "__main__":
    main()
