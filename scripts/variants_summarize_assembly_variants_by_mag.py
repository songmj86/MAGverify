#!/usr/bin/env python3
"""Generate MAG-associated and MAG-level assembly-variant outputs.

This script does not rerun mapping, anvi'o profiling, or variant classification.
It joins the existing site-wise and contig-wise assembly-allele support tables to
the MAG-contig ownership map, validates their internal consistency, and writes:

1. a site table with MAG_ID added to every disagreement locus; and
2. one aggregate row for every MAG, including MAGs with no disagreement loci.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple

CLASS_UNSUPPORTED = "Unsupported_assembly_variant"
CLASS_MINORITY = "Minority_assembly_variant"
CLASS_NO_DISAGREEMENT = "No_assembly_allele_disagreement"
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


def parse_args() -> argparse.Namespace:
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
    """Create a stable fallback event ID for legacy disagreement tables.

    New tables contain ``event_id`` from the complete exported anvi'o row.
    Legacy tables do not, so all available disagreement fields are hashed. This
    distinguishes same-position INDEL rows with different alleles or support
    counts without assuming that a genomic coordinate identifies one event.
    """
    payload = {key: value for key, value in sorted(row.items()) if key != "event_id"}
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    return f"legacy:{digest}"


def require_nonempty(path: str, label: str) -> None:
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise SystemExit(f"ERROR: {label} not found or empty: {path}")


def normalize(name: str) -> str:
    return name.strip().lstrip("#").lower().replace("-", "_").replace(" ", "_")


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


def load_mag_map(
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

def main() -> None:
    args = parse_args()
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

    contig_to_mag, mag_to_contigs = load_mag_map(args.mag_contig_map)
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


if __name__ == "__main__":
    main()
