#!/usr/bin/env python3
"""Annotate false-duplication boundaries with validated HiFi clipping.

Clipping is reported as a co-occurring mapping signal and never changes the
HiFi boundary-based false-duplication call. A false-duplication boundary is
matched only to clipping positions at the last base on its left or the first
base on its right. In 1-based coordinates these are boundary-1 and boundary.
This two-position rule accommodates the way left- and right-clipped alignments
record the same junction without introducing an arbitrary distance window.
"""

from __future__ import annotations

import argparse
import csv
import os
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, Sequence, Set, Tuple


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", required=True)
    parser.add_argument("--boundaries", required=True)
    parser.add_argument("--contigs", required=True)
    parser.add_argument(
        "--clipping",
        required=True,
        help="Validated clipping TSV, or NA when clipping was not run.",
    )
    parser.add_argument("--out-boundaries", required=True)
    parser.add_argument("--out-pairs", required=True)
    parser.add_argument("--out-contigs", required=True)
    parser.add_argument("--out-mags", required=True)
    return parser.parse_args()


def read_rows(path: str, required: Set[str], label: str) -> Tuple[List[str], List[Dict[str, str]]]:
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        fail(f"{label} missing or empty: {path}")
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        missing = required - set(fields)
        if missing:
            fail(f"{label} lacks columns: {sorted(missing)}")
        rows = [dict(row) for row in reader]
    return fields, rows


def parse_positive_int(value: str, label: str) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        fail(f"invalid integer for {label}: {value!r}")
    if parsed < 1:
        fail(f"{label} must be positive: {parsed}")
    return parsed


@dataclass(frozen=True)
class ClippingSite:
    position: int
    coverage: int
    clipping_count: int
    clipping_ratio: float
    validated: bool


def read_contigs(path: str) -> Tuple[List[Dict[str, str]], Dict[str, int]]:
    _fields, rows = read_rows(
        path,
        {"MAG_ID", "assembler", "contig", "length"},
        "false-duplication contig table",
    )
    lengths: Dict[str, int] = {}
    for line_number, row in enumerate(rows, start=2):
        contig = row["contig"].strip()
        if not contig or contig in lengths:
            fail(f"blank or duplicate contig in contig table at line {line_number}: {contig!r}")
        lengths[contig] = parse_positive_int(row["length"], f"contig length at line {line_number}")
    if not rows:
        fail("false-duplication contig table contains no rows")
    return rows, lengths


def read_clipping(
    path: str,
    lengths: Mapping[str, int],
) -> Tuple[bool, Dict[str, Dict[int, ClippingSite]]]:
    if not path or path.upper() == "NA":
        return False, {}
    _fields, rows = read_rows(
        path,
        {
            "contig", "length", "pos_1based", "coverage",
            "clipping_count", "clipping_ratio", "Clipping",
        },
        "validated clipping table",
    )
    result: Dict[str, Dict[int, ClippingSite]] = defaultdict(dict)
    for line_number, row in enumerate(rows, start=2):
        contig = row["contig"].strip()
        if contig not in lengths:
            fail(f"clipping contig absent from false-duplication scope at line {line_number}: {contig}")
        length = parse_positive_int(row["length"], f"clipping contig length at line {line_number}")
        if length != lengths[contig]:
            fail(f"stale clipping contig length at line {line_number}: {contig}")
        position = parse_positive_int(row["pos_1based"], f"clipping position at line {line_number}")
        if position > length or position in result[contig]:
            fail(f"invalid or duplicate clipping position at line {line_number}: {contig}:{position}")
        coverage = parse_positive_int(row["coverage"], f"clipping coverage at line {line_number}")
        clipping_count = parse_positive_int(
            row["clipping_count"], f"clipping count at line {line_number}"
        )
        if clipping_count > coverage:
            fail(f"clipping count exceeds coverage at line {line_number}")
        try:
            ratio = float(row["clipping_ratio"])
        except (TypeError, ValueError) as exc:
            fail(f"invalid clipping ratio at line {line_number}: {row['clipping_ratio']!r}")
        expected = clipping_count / coverage
        if not 0 < ratio <= 1 or abs(ratio - expected) > 5.0e-7:
            fail(f"inconsistent clipping ratio at line {line_number}: {contig}:{position}")
        call = row["Clipping"].strip().upper()
        if call not in {"Y", "N", "NA"}:
            fail(f"invalid Clipping call at line {line_number}: {call!r}")
        result[contig][position] = ClippingSite(
            position=position,
            coverage=coverage,
            clipping_count=clipping_count,
            clipping_ratio=ratio,
            validated=(call == "Y"),
        )
    return True, result


def defining_boundaries(pair: Mapping[str, str]) -> Set[str]:
    """Return the LR-unsupported boundary names that establish a positive call."""
    if pair["LR_FALSE_DUPLICATION"] != "Y":
        return set()
    duplication_type = pair["DUPLICATION_TYPE"]
    if duplication_type == "TANDEM":
        if pair["LR_COPY1_RIGHT_STATUS"] != "UNSUPPORTED":
            fail(f"positive tandem pair lacks an unsupported internal boundary: {pair['pair_id']}")
        return {"copy1_right"}

    copy1_unsupported = (
        pair["LR_COPY1_LEFT_STATUS"] == "UNSUPPORTED"
        and pair["LR_COPY1_RIGHT_STATUS"] == "UNSUPPORTED"
    )
    copy2_unsupported = (
        pair["LR_COPY2_LEFT_STATUS"] == "UNSUPPORTED"
        and pair["LR_COPY2_RIGHT_STATUS"] == "UNSUPPORTED"
    )
    copy1_supported = (
        pair["LR_COPY1_LEFT_STATUS"] == "SUPPORTED"
        and pair["LR_COPY1_RIGHT_STATUS"] == "SUPPORTED"
    )
    copy2_supported = (
        pair["LR_COPY2_LEFT_STATUS"] == "SUPPORTED"
        and pair["LR_COPY2_RIGHT_STATUS"] == "SUPPORTED"
    )
    if copy1_unsupported and copy2_supported:
        return {"copy1_left", "copy1_right"}
    if copy2_unsupported and copy1_supported:
        return {"copy2_left", "copy2_right"}
    fail(f"positive non-tandem pair lacks the expected supported/unsupported pattern: {pair['pair_id']}")


def matched_sites(
    contig: str,
    boundary_position: int,
    clipping: Mapping[str, Mapping[int, ClippingSite]],
) -> List[ClippingSite]:
    positions = {boundary_position - 1, boundary_position}
    return [
        clipping[contig][position]
        for position in sorted(positions)
        if position >= 1 and position in clipping.get(contig, {})
    ]


def format_positions(sites: Sequence[ClippingSite]) -> str:
    return ";".join(str(site.position) for site in sites) or "NA"


def format_offsets(sites: Sequence[ClippingSite], boundary_position: int) -> str:
    return ";".join(str(site.position - boundary_position) for site in sites) or "NA"


def max_or_na(values: Iterable[float | int], decimal: bool = False) -> str:
    materialized = list(values)
    if not materialized:
        return "NA"
    maximum = max(materialized)
    return f"{float(maximum):.6f}" if decimal else str(maximum)


def write_tsv(path: str, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    tmp = f"{path}.tmp.{os.getpid()}"
    with open(tmp, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    os.replace(tmp, path)


def main() -> None:
    args = parse_args()
    pair_required = {
        "pair_id", "MAG_ID", "assembler", "PAIR_SCOPE", "DUPLICATION_TYPE", "ORIENTATION",
        "copy1_contig", "copy2_contig", "LR_COPY1_LEFT_STATUS", "LR_COPY1_RIGHT_STATUS",
        "LR_COPY2_LEFT_STATUS", "LR_COPY2_RIGHT_STATUS", "LR_FALSE_DUPLICATION",
    }
    _pair_fields, pairs = read_rows(args.pairs, pair_required, "false-duplication pair table")
    _boundary_fields, boundaries = read_rows(
        args.boundaries,
        {
            "pair_id", "MAG_ID", "assembler", "PAIR_SCOPE", "DUPLICATION_TYPE",
            "ORIENTATION", "boundary", "contig", "boundary_position_1based",
        },
        "false-duplication boundary table",
    )
    contigs, lengths = read_contigs(args.contigs)
    clipping_available, clipping = read_clipping(args.clipping, lengths)

    pair_by_id: Dict[str, Dict[str, str]] = {}
    for pair in pairs:
        pair_id = pair["pair_id"]
        if not pair_id or pair_id in pair_by_id:
            fail(f"blank or duplicate pair_id: {pair_id!r}")
        if pair["LR_FALSE_DUPLICATION"] not in {"Y", "N", "NA"}:
            fail(f"invalid LR false-duplication call for {pair_id}")
        pair_by_id[pair_id] = pair

    boundary_by_pair: Dict[str, Dict[str, Dict[str, str]]] = defaultdict(dict)
    for row in boundaries:
        pair_id = row["pair_id"]
        boundary_name = row["boundary"]
        if pair_id not in pair_by_id:
            fail(f"boundary row references an unknown pair: {pair_id}")
        if boundary_name not in {"copy1_left", "copy1_right", "copy2_left", "copy2_right"}:
            fail(f"invalid boundary name for {pair_id}: {boundary_name}")
        if boundary_name in boundary_by_pair[pair_id]:
            fail(f"duplicate boundary row for {pair_id}/{boundary_name}")
        if row["contig"] not in lengths:
            fail(f"boundary contig absent from contig table: {row['contig']}")
        boundary_by_pair[pair_id][boundary_name] = row
    for pair_id in pair_by_id:
        if set(boundary_by_pair[pair_id]) != {
            "copy1_left", "copy1_right", "copy2_left", "copy2_right"
        }:
            fail(f"pair does not have exactly four boundary rows: {pair_id}")

    boundary_output: List[Dict[str, object]] = []
    boundary_annotations: Dict[Tuple[str, str], Dict[str, object]] = {}
    for pair in pairs:
        pair_id = pair["pair_id"]
        defining = defining_boundaries(pair)
        for boundary_name in ("copy1_left", "copy1_right", "copy2_left", "copy2_right"):
            row = boundary_by_pair[pair_id][boundary_name]
            position_text = row["boundary_position_1based"].strip()
            is_defining = boundary_name in defining
            if position_text.upper() == "NA":
                sites: List[ClippingSite] = []
                position: int | None = None
            else:
                position = parse_positive_int(
                    position_text, f"boundary position for {pair_id}/{boundary_name}"
                )
                if position > lengths[row["contig"]]:
                    fail(f"boundary position exceeds contig length: {pair_id}/{boundary_name}")
                sites = matched_sites(row["contig"], position, clipping) if clipping_available else []
            validated = [site for site in sites if site.validated]

            if not is_defining:
                interpretation = "NOT_A_FALSE_DUPLICATION_DEFINING_BOUNDARY"
            elif not clipping_available:
                interpretation = "CLIPPING_DATA_NOT_AVAILABLE"
            elif validated:
                interpretation = "FALSE_DUPLICATION_BOUNDARY_WITH_VALIDATED_CLIPPING"
            elif sites:
                interpretation = "CLIPPING_CANDIDATE_NOT_VALIDATED_AT_FALSE_DUPLICATION_BOUNDARY"
            else:
                interpretation = "NO_CLIPPING_CANDIDATE_AT_FALSE_DUPLICATION_BOUNDARY"

            annotated: Dict[str, object] = {
                "pair_id": pair_id,
                "MAG_ID": row["MAG_ID"],
                "assembler": row["assembler"],
                "PAIR_SCOPE": row["PAIR_SCOPE"],
                "DUPLICATION_TYPE": row["DUPLICATION_TYPE"],
                "ORIENTATION": row["ORIENTATION"],
                "boundary": boundary_name,
                "contig": row["contig"],
                "boundary_position_1based": position_text,
                "FALSE_DUPLICATION_PAIR_Y_N_NA": pair["LR_FALSE_DUPLICATION"],
                "FALSE_DUPLICATION_DEFINING_BOUNDARY_Y_N": "Y" if is_defining else "N",
                "CLIPPING_DATA_AVAILABLE_Y_N": "Y" if clipping_available else "N",
                "CLIPPING_CANDIDATE_AT_BOUNDARY_Y_N_NA": (
                    "Y" if sites else "N" if clipping_available and position is not None else "NA"
                ),
                "VALIDATED_CLIPPING_AT_BOUNDARY_Y_N_NA": (
                    "Y" if validated else "N" if clipping_available and position is not None else "NA"
                ),
                "matched_clipping_positions_1based": format_positions(sites),
                "matched_clipping_offsets_from_boundary_bp": (
                    format_offsets(sites, position) if position is not None else "NA"
                ),
                "matched_validated_clipping_positions_1based": format_positions(validated),
                "max_matched_clipping_coverage": max_or_na(site.coverage for site in sites),
                "max_matched_clipping_count": max_or_na(site.clipping_count for site in sites),
                "max_matched_clipping_ratio": max_or_na(
                    (site.clipping_ratio for site in sites), decimal=True
                ),
                "max_validated_clipping_ratio": max_or_na(
                    (site.clipping_ratio for site in validated), decimal=True
                ),
                "FALSE_DUPLICATION_CLIPPING_INTERPRETATION": interpretation,
            }
            boundary_output.append(annotated)
            boundary_annotations[(pair_id, boundary_name)] = annotated

    boundary_output_fields = [
        "pair_id", "MAG_ID", "assembler", "PAIR_SCOPE", "DUPLICATION_TYPE", "ORIENTATION",
        "boundary", "contig", "boundary_position_1based", "FALSE_DUPLICATION_PAIR_Y_N_NA",
        "FALSE_DUPLICATION_DEFINING_BOUNDARY_Y_N", "CLIPPING_DATA_AVAILABLE_Y_N",
        "CLIPPING_CANDIDATE_AT_BOUNDARY_Y_N_NA", "VALIDATED_CLIPPING_AT_BOUNDARY_Y_N_NA",
        "matched_clipping_positions_1based", "matched_clipping_offsets_from_boundary_bp",
        "matched_validated_clipping_positions_1based", "max_matched_clipping_coverage",
        "max_matched_clipping_count", "max_matched_clipping_ratio",
        "max_validated_clipping_ratio", "FALSE_DUPLICATION_CLIPPING_INTERPRETATION",
    ]
    write_tsv(args.out_boundaries, boundary_output_fields, boundary_output)

    pair_output: List[Dict[str, object]] = []
    for pair in pairs:
        pair_id = pair["pair_id"]
        defining = defining_boundaries(pair)
        defining_rows = [boundary_annotations[(pair_id, name)] for name in sorted(defining)]
        candidate_count = sum(
            row["CLIPPING_CANDIDATE_AT_BOUNDARY_Y_N_NA"] == "Y" for row in defining_rows
        )
        validated_count = sum(
            row["VALIDATED_CLIPPING_AT_BOUNDARY_Y_N_NA"] == "Y" for row in defining_rows
        )
        positive = pair["LR_FALSE_DUPLICATION"] == "Y"
        if not positive:
            any_call = all_call = "NA"
            interpretation = "NOT_APPLICABLE_FALSE_DUPLICATION_NOT_DETECTED"
        elif not clipping_available:
            any_call = all_call = "NA"
            interpretation = "CLIPPING_DATA_NOT_AVAILABLE"
        elif validated_count == len(defining):
            any_call = all_call = "Y"
            interpretation = "ALL_FALSE_DUPLICATION_BOUNDARIES_HAVE_VALIDATED_CLIPPING"
        elif validated_count > 0:
            any_call, all_call = "Y", "N"
            interpretation = "SOME_FALSE_DUPLICATION_BOUNDARIES_HAVE_VALIDATED_CLIPPING"
        else:
            any_call = all_call = "N"
            interpretation = "NO_VALIDATED_CLIPPING_AT_FALSE_DUPLICATION_BOUNDARIES"
        pair_output.append({
            "pair_id": pair_id,
            "MAG_ID": pair["MAG_ID"],
            "assembler": pair["assembler"],
            "PAIR_SCOPE": pair["PAIR_SCOPE"],
            "DUPLICATION_TYPE": pair["DUPLICATION_TYPE"],
            "ORIENTATION": pair["ORIENTATION"],
            "FALSE_DUPLICATION_PAIR_Y_N_NA": pair["LR_FALSE_DUPLICATION"],
            "CLIPPING_DATA_AVAILABLE_Y_N": "Y" if clipping_available else "N",
            "n_false_duplication_defining_boundaries": len(defining),
            "n_defining_boundaries_with_clipping_candidate": candidate_count,
            "n_defining_boundaries_with_validated_clipping": validated_count,
            "ANY_VALIDATED_CLIPPING_AT_FALSE_DUPLICATION_BOUNDARY_Y_N_NA": any_call,
            "ALL_FALSE_DUPLICATION_BOUNDARIES_HAVE_VALIDATED_CLIPPING_Y_N_NA": all_call,
            "FALSE_DUPLICATION_CLIPPING_INTERPRETATION": interpretation,
        })

    pair_output_fields = [
        "pair_id", "MAG_ID", "assembler", "PAIR_SCOPE", "DUPLICATION_TYPE", "ORIENTATION",
        "FALSE_DUPLICATION_PAIR_Y_N_NA", "CLIPPING_DATA_AVAILABLE_Y_N",
        "n_false_duplication_defining_boundaries",
        "n_defining_boundaries_with_clipping_candidate",
        "n_defining_boundaries_with_validated_clipping",
        "ANY_VALIDATED_CLIPPING_AT_FALSE_DUPLICATION_BOUNDARY_Y_N_NA",
        "ALL_FALSE_DUPLICATION_BOUNDARIES_HAVE_VALIDATED_CLIPPING_Y_N_NA",
        "FALSE_DUPLICATION_CLIPPING_INTERPRETATION",
    ]
    write_tsv(args.out_pairs, pair_output_fields, pair_output)

    # A contig-level co-occurrence is assigned only to the contig carrying the
    # LR-unsupported boundary. This avoids assigning clipping from a false-copy
    # contig to the supported copy in a same-MAG inter-contig pair.
    defining_rows_by_contig: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for row in boundary_output:
        if (
            row["FALSE_DUPLICATION_PAIR_Y_N_NA"] == "Y"
            and row["FALSE_DUPLICATION_DEFINING_BOUNDARY_Y_N"] == "Y"
        ):
            defining_rows_by_contig[str(row["contig"])].append(row)

    summary_fields = [
        "MAG_ID", "assembler", "contig", "length", "CLIPPING_DATA_AVAILABLE_Y_N",
        "n_false_duplication_pairs", "n_false_duplication_pairs_with_any_validated_clipping",
        "n_false_duplication_defining_boundaries",
        "n_false_duplication_defining_boundaries_with_validated_clipping",
        "ANY_FALSE_DUPLICATION_PAIR_WITH_VALIDATED_CLIPPING_Y_N_NA",
        "FALSE_DUPLICATION_CLIPPING_SUMMARY",
    ]
    contig_output: List[Dict[str, object]] = []
    for contig_row in contigs:
        defining_rows = defining_rows_by_contig.get(contig_row["contig"], [])
        positive_pair_ids = {str(row["pair_id"]) for row in defining_rows}
        pairs_with_clipping = {
            str(row["pair_id"])
            for row in defining_rows
            if row["VALIDATED_CLIPPING_AT_BOUNDARY_Y_N_NA"] == "Y"
        }
        pair_with_clipping = len(pairs_with_clipping)
        if not clipping_available:
            any_call = "NA"
            summary = "CLIPPING_DATA_NOT_AVAILABLE"
        elif pair_with_clipping:
            any_call = "Y"
            summary = "FALSE_DUPLICATION_WITH_VALIDATED_CLIPPING_CO_OCCURRENCE"
        elif positive_pair_ids:
            any_call = "N"
            summary = "FALSE_DUPLICATION_WITHOUT_VALIDATED_CLIPPING_CO_OCCURRENCE"
        else:
            any_call = "N"
            summary = "NO_FALSE_DUPLICATION_PAIR"
        contig_output.append({
            "MAG_ID": contig_row["MAG_ID"],
            "assembler": contig_row["assembler"],
            "contig": contig_row["contig"],
            "length": contig_row["length"],
            "CLIPPING_DATA_AVAILABLE_Y_N": "Y" if clipping_available else "N",
            "n_false_duplication_pairs": len(positive_pair_ids),
            "n_false_duplication_pairs_with_any_validated_clipping": pair_with_clipping,
            "n_false_duplication_defining_boundaries": len(defining_rows),
            "n_false_duplication_defining_boundaries_with_validated_clipping": sum(
                row["VALIDATED_CLIPPING_AT_BOUNDARY_Y_N_NA"] == "Y" for row in defining_rows
            ),
            "ANY_FALSE_DUPLICATION_PAIR_WITH_VALIDATED_CLIPPING_Y_N_NA": any_call,
            "FALSE_DUPLICATION_CLIPPING_SUMMARY": summary,
        })
    write_tsv(args.out_contigs, summary_fields, contig_output)

    mag_pair_rows: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for row in pair_output:
        mag = str(row["MAG_ID"])
        if mag not in {"", "NA"} and row["FALSE_DUPLICATION_PAIR_Y_N_NA"] == "Y":
            mag_pair_rows[mag].append(row)
    mag_contig_count: Dict[str, int] = defaultdict(int)
    for row in contigs:
        mag = row["MAG_ID"]
        if mag not in {"", "NA"}:
            mag_contig_count[mag] += 1

    mag_fields = [
        "MAG_ID", "n_contigs", "CLIPPING_DATA_AVAILABLE_Y_N", "n_false_duplication_pairs",
        "n_false_duplication_pairs_with_any_validated_clipping",
        "n_false_duplication_defining_boundaries",
        "n_false_duplication_defining_boundaries_with_validated_clipping",
        "ANY_FALSE_DUPLICATION_PAIR_WITH_VALIDATED_CLIPPING_Y_N_NA",
        "FALSE_DUPLICATION_CLIPPING_SUMMARY",
    ]
    mag_output: List[Dict[str, object]] = []
    for mag in sorted(mag_contig_count):
        positive_rows = mag_pair_rows.get(mag, [])
        pair_with_clipping = sum(
            row["ANY_VALIDATED_CLIPPING_AT_FALSE_DUPLICATION_BOUNDARY_Y_N_NA"] == "Y"
            for row in positive_rows
        )
        if not clipping_available:
            any_call = "NA"
            summary = "CLIPPING_DATA_NOT_AVAILABLE"
        elif pair_with_clipping:
            any_call = "Y"
            summary = "FALSE_DUPLICATION_WITH_VALIDATED_CLIPPING_CO_OCCURRENCE"
        elif positive_rows:
            any_call = "N"
            summary = "FALSE_DUPLICATION_WITHOUT_VALIDATED_CLIPPING_CO_OCCURRENCE"
        else:
            any_call = "N"
            summary = "NO_FALSE_DUPLICATION_PAIR"
        mag_output.append({
            "MAG_ID": mag,
            "n_contigs": mag_contig_count[mag],
            "CLIPPING_DATA_AVAILABLE_Y_N": "Y" if clipping_available else "N",
            "n_false_duplication_pairs": len(positive_rows),
            "n_false_duplication_pairs_with_any_validated_clipping": pair_with_clipping,
            "n_false_duplication_defining_boundaries": sum(
                int(row["n_false_duplication_defining_boundaries"]) for row in positive_rows
            ),
            "n_false_duplication_defining_boundaries_with_validated_clipping": sum(
                int(row["n_defining_boundaries_with_validated_clipping"]) for row in positive_rows
            ),
            "ANY_FALSE_DUPLICATION_PAIR_WITH_VALIDATED_CLIPPING_Y_N_NA": any_call,
            "FALSE_DUPLICATION_CLIPPING_SUMMARY": summary,
        })
    write_tsv(args.out_mags, mag_fields, mag_output)

    print(f"Clipping data available: {'Y' if clipping_available else 'N'}")
    print(f"False-duplication pairs annotated: {len(pair_output)}")
    print(
        "Positive false-duplication pairs with any validated clipping: "
        f"{sum(row['ANY_VALIDATED_CLIPPING_AT_FALSE_DUPLICATION_BOUNDARY_Y_N_NA'] == 'Y' for row in pair_output)}"
    )


if __name__ == "__main__":
    main()
