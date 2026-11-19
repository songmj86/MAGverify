#!/usr/bin/env python3
"""Create intuitive MAG-level clipping summaries from all-region site calls."""

from __future__ import annotations

import argparse
import csv
import math
import os
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple

# Put the biological Y/N calls first so the final table can be interpreted at a glance.
SUMMARY_HEADER = (
    "MAG_ID",
    "n_contigs",
    "total_length_bp",
    "Internal_clipping_ge1",
    "Terminal_clipping_ge1",
    "Clipping_ge1",
    "Internal_clipping",
    "Terminal_clipping",
    "Clipping_analyzable",
    "Clipping",
    "Stringent_clipping",
    "n_internal_clipping_ge1_sites",
    "n_terminal_clipping_ge1_sites",
    "n_clipping_ge1_sites",
    "n_internal_clipping_sites",
    "n_terminal_clipping_sites",
    "n_clipping_sites",
    "n_stringent_clipping_sites",
    "n_internal_candidate_sites",
    "n_terminal_candidate_sites",
    "n_candidate_sites",
    "n_stringent_candidate_sites",
    "n_candidates_with_spanning",
    "n_candidates_with_spanning_I",
    "n_candidates_with_spanning_D",
    "n_candidates_with_spanning_N",
    "n_candidates_with_spanning_any_IDN",
    "n_contigs_with_candidate",
    "n_contigs_with_stringent_candidate",
    "n_contigs_with_clipping",
    "n_contigs_with_stringent_clipping",
    "candidate_contigs",
    "stringent_candidate_contigs",
    "clipping_contigs",
    "stringent_clipping_contigs",
    "max_anvio_coverage",
    "max_anvio_clipping_count",
    "max_unique_clipped_hifi_reads",
    "max_anvio_clipping_ratio",
    "Candidate",
    "Stringent_candidate",
)

# Contig-level table: same calls and counts as the MAG summary, computed per
# contig. MAG_ID is NA for contigs that are not assigned to any MAG.
CONTIG_DROPPED = {
    "MAG_ID", "n_contigs", "total_length_bp",
    "n_contigs_with_candidate", "n_contigs_with_stringent_candidate",
    "n_contigs_with_clipping", "n_contigs_with_stringent_clipping",
    "candidate_contigs", "stringent_candidate_contigs",
    "clipping_contigs", "stringent_clipping_contigs",
}
CONTIG_HEADER = ("contig", "MAG_ID", "length_bp",
                 *(field for field in SUMMARY_HEADER if field not in CONTIG_DROPPED))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--clipping", required=True)
    parser.add_argument("--mag-contig-map", required=True)
    parser.add_argument("--fai", required=True)
    parser.add_argument("--positive-site-output", required=True)
    parser.add_argument("--mag-site-output", required=True)
    parser.add_argument("--mag-positive-site-output", required=True)
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--positive-summary-output", required=True)
    parser.add_argument("--contig-summary-output", required=True,
                        help="One row per assembly contig (all contigs in the FAI)")
    return parser.parse_args()


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def require_nonempty(path: str, label: str) -> None:
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        fail(f"{label} not found or empty: {path}")


def atomic_write(path: str, header: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    destination = Path(path).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="", delete=False,
        dir=str(destination.parent), prefix=destination.name + ".", suffix=".tmp",
    )
    try:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
        handle.close()
        os.replace(handle.name, destination)
    finally:
        try:
            handle.close()
        except Exception:
            pass
        if os.path.exists(handle.name):
            os.unlink(handle.name)


def load_fai(path: str) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    with open(path, encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            fields = raw.rstrip("\n\r").split("\t")
            if len(fields) < 2:
                fail(f"malformed FAI row {line_number}")
            try:
                length = int(fields[1])
            except ValueError as error:
                fail(f"invalid FAI length at row {line_number}: {error}")
            if not fields[0] or length < 1 or fields[0] in lengths:
                fail(f"invalid or duplicate FAI contig: {fields[0]}")
            lengths[fields[0]] = length
    return lengths


def resolve_column(fields: Sequence[str], choices: Sequence[str], label: str) -> str:
    normalized = {field.strip().lower(): field for field in fields}
    for choice in choices:
        if choice in normalized:
            return normalized[choice]
    fail(f"MAG-contig map lacks a {label} column")


def load_mag_map(path: str, lengths: Mapping[str, int]) -> Tuple[Dict[str, str], Dict[str, List[str]]]:
    contig_to_mag: Dict[str, str] = {}
    mag_to_contigs: Dict[str, List[str]] = defaultdict(list)
    seen = set()
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        mag_col = resolve_column(reader.fieldnames or [], ("mag_id", "id"), "MAG")
        contig_col = resolve_column(reader.fieldnames or [], ("contig", "original_contig"), "contig")
        for line_number, row in enumerate(reader, start=2):
            mag = (row.get(mag_col) or "").strip()
            contig = (row.get(contig_col) or "").strip()
            if not mag or not contig or (mag, contig) in seen:
                fail(f"blank or duplicate MAG map row {line_number}")
            seen.add((mag, contig))
            if contig not in lengths:
                fail(f"MAG contig absent from FAI: {contig}")
            if contig in contig_to_mag and contig_to_mag[contig] != mag:
                fail(f"contig assigned to multiple MAGs: {contig}")
            contig_to_mag[contig] = mag
            mag_to_contigs[mag].append(contig)
    for mag in mag_to_contigs:
        mag_to_contigs[mag].sort()
    return contig_to_mag, dict(mag_to_contigs)


def normalize_call(value: str, line_number: int, label: str) -> str:
    value = value.strip().upper()
    if value in {"Y", "YES"}:
        return "Y"
    if value in {"N", "NO"}:
        return "N"
    fail(f"invalid {label} at row {line_number}: {value}")


def load_clipping(path: str, lengths: Mapping[str, int]) -> Tuple[List[str], List[Dict[str, str]]]:
    records: List[Dict[str, str]] = []
    seen = set()
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        required = {
            "contig", "length", "pos_1based", "coverage", "clipping_count",
            "clipping_ratio", "default_candidate", "stringent_candidate", "Clipping_ge1",
            "Clipping_analyzable", "Clipping", "Stringent_clipping", "region_type", "n_unique_all_clipped_qnames",
            "n_unique_exact_spanning_qnames",
            "exact_spanning_min_flank_bp_min", "exact_spanning_min_flank_bp_median",
            "exact_spanning_min_flank_bp_max",
            "n_unique_spanning_with_I_qnames",
            "spanning_I_length_bp_min", "spanning_I_length_bp_median",
            "spanning_I_length_bp_max",
            "n_unique_spanning_with_D_qnames",
            "spanning_D_length_bp_min", "spanning_D_length_bp_median",
            "spanning_D_length_bp_max",
            "n_unique_spanning_with_N_qnames",
            "spanning_N_length_bp_min", "spanning_N_length_bp_median",
            "spanning_N_length_bp_max",
            "n_unique_spanning_with_any_IDN_qnames",
        }
        missing = required - set(fields)
        if missing:
            fail(f"clipping table lacks columns: {sorted(missing)}")

        for line_number, row in enumerate(reader, start=2):
            contig = (row.get("contig") or "").strip()
            try:
                length = int(row["length"])
                position = int(row["pos_1based"])
                coverage = int(row["coverage"])
                clipped = int(row["clipping_count"])
                ratio = float(row["clipping_ratio"])
            except (TypeError, ValueError) as error:
                fail(f"invalid clipping row {line_number}: {error}")
            if contig not in lengths or length != lengths[contig] or not 1 <= position <= length:
                fail(f"stale clipping coordinate: {contig}:{position}")
            if (contig, position) in seen:
                fail(f"duplicate clipping site: {contig}:{position}")
            seen.add((contig, position))
            if coverage < 1 or not 0 <= clipped <= coverage:
                fail(f"invalid clipping counts: {contig}:{position}")
            if not math.isclose(ratio, clipped / coverage, rel_tol=0, abs_tol=5.0e-7):
                fail(f"inconsistent clipping ratio: {contig}:{position}")
            for field in ("default_candidate", "stringent_candidate", "Clipping_ge1"):
                row[field] = normalize_call(row[field], line_number, field)
            analyzable = (row.get("Clipping_analyzable") or "").strip().upper()
            if analyzable not in {"Y", "N"}:
                fail(f"invalid Clipping_analyzable at row {line_number}: {analyzable}")
            row["Clipping_analyzable"] = analyzable
            for field in ("Clipping", "Stringent_clipping"):
                value = (row.get(field) or "").strip().upper()
                if value not in {"Y", "N", "NA"}:
                    fail(f"invalid {field} at row {line_number}: {value}")
                row[field] = value
            if analyzable == "Y" and row["Clipping"] not in {"Y", "N"}:
                fail(f"Clipping_analyzable=Y requires Clipping=Y/N at row {line_number}")
            if analyzable == "N" and row["Clipping"] != "NA":
                fail(f"Clipping_analyzable=N requires Clipping=NA at row {line_number}")
            if row["Clipping"] == "Y" and row["Clipping_ge1"] != "Y":
                fail(f"Clipping=Y requires Clipping_ge1=Y at row {line_number}")
            if row["Stringent_clipping"] == "Y" and row["Clipping"] != "Y":
                fail(f"Stringent_clipping=Y requires Clipping=Y at row {line_number}")
            if analyzable == "N" and row["Stringent_clipping"] != "NA":
                fail(f"Clipping_analyzable=N requires Stringent_clipping=NA at row {line_number}")

            hifi_span_fields = (
                "n_unique_exact_spanning_qnames", "n_unique_spanning_with_I_qnames",
                "n_unique_spanning_with_D_qnames", "n_unique_spanning_with_N_qnames",
                "n_unique_spanning_with_any_IDN_qnames",
            )
            hifi_span_counts = {}
            for field in hifi_span_fields:
                try:
                    value = int((row.get(field) or "0").strip())
                except ValueError:
                    fail(f"invalid {field} at row {line_number}: {row.get(field)!r}")
                if value < 0:
                    fail(f"negative {field} at row {line_number}: {value}")
                row[field] = str(value)
                hifi_span_counts[field] = value
            hifi_span_total = hifi_span_counts["n_unique_exact_spanning_qnames"]
            hifi_span_any = hifi_span_counts["n_unique_spanning_with_any_IDN_qnames"]
            hifi_specific = [
                hifi_span_counts["n_unique_spanning_with_I_qnames"],
                hifi_span_counts["n_unique_spanning_with_D_qnames"],
                hifi_span_counts["n_unique_spanning_with_N_qnames"],
            ]
            if hifi_span_any > hifi_span_total or any(v > hifi_span_any for v in hifi_specific):
                fail(f"inconsistent HiFi spanning I/D/N audit counts at row {line_number}")

            # Length annotations are audit-only. A missing event is represented
            # by NA; observed lengths must be non-negative numeric values.
            length_audit_fields = (
                "exact_spanning_min_flank_bp_min", "exact_spanning_min_flank_bp_median",
                "exact_spanning_min_flank_bp_max",
                "spanning_I_length_bp_min", "spanning_I_length_bp_median", "spanning_I_length_bp_max",
                "spanning_D_length_bp_min", "spanning_D_length_bp_median", "spanning_D_length_bp_max",
                "spanning_N_length_bp_min", "spanning_N_length_bp_median", "spanning_N_length_bp_max",
            )
            for field in length_audit_fields:
                raw = (row.get(field) or "NA").strip()
                if raw.upper() == "NA":
                    row[field] = "NA"
                    continue
                try:
                    value = float(raw)
                except ValueError:
                    fail(f"invalid {field} at row {line_number}: {raw}")
                if value < 0:
                    fail(f"negative {field} at row {line_number}: {raw}")
                row[field] = raw

            region = (row.get("region_type") or "").strip().upper()
            if region not in {"INTERNAL", "TERMINAL"}:
                fail(f"invalid region_type at row {line_number}: {region}")
            row["region_type"] = region
            records.append(dict(row))
    preferred = [
        "contig", "length", "pos_1based", "region_type",
        "Clipping_analyzable", "Clipping", "Clipping_ge1", "Stringent_clipping",
    ]
    fields = [*preferred, *(field for field in fields if field not in preferred)]
    records.sort(key=lambda row: (row["contig"], int(row["pos_1based"])))
    return fields, records


def joined_contigs(rows: Sequence[Mapping[str, str]]) -> str:
    values = sorted({row["contig"] for row in rows})
    return ";".join(values) if values else "NA"


def summary_rows(
    lengths: Mapping[str, int], mag_to_contigs: Mapping[str, Sequence[str]],
    by_mag: Mapping[str, Sequence[Mapping[str, str]]],
) -> List[List[object]]:
    result: List[List[object]] = []
    for mag in sorted(mag_to_contigs):
        members = list(mag_to_contigs[mag])
        records = list(by_mag.get(mag, []))
        candidates = [row for row in records if row["default_candidate"] == "Y"]
        stringent = [row for row in records if row["stringent_candidate"] == "Y"]
        analyzable_candidates = [row for row in candidates if row.get("Clipping_analyzable") == "Y"]
        not_analyzable_candidates = [row for row in candidates if row.get("Clipping_analyzable") == "N"]
        permissive_positives = [row for row in records if row["Clipping_ge1"] == "Y"]
        positives = [row for row in records if row["Clipping"] == "Y"]
        # Presence is decidable as soon as one analyzable positive site exists.
        # A negative MAG call requires every candidate site to be analyzable.
        if positives:
            clipping_analyzable_call, clipping_call = "Y", "Y"
        elif not_analyzable_candidates:
            clipping_analyzable_call, clipping_call = "N", "NA"
        else:
            clipping_analyzable_call, clipping_call = "Y", "N"
        stringent_positives = [row for row in records if row["Stringent_clipping"] == "Y"]
        internal_candidates = [row for row in candidates if row["region_type"] == "INTERNAL"]
        terminal_candidates = [row for row in candidates if row["region_type"] == "TERMINAL"]
        internal_permissive = [row for row in permissive_positives if row["region_type"] == "INTERNAL"]
        terminal_permissive = [row for row in permissive_positives if row["region_type"] == "TERMINAL"]
        internal_positive = [row for row in positives if row["region_type"] == "INTERNAL"]
        terminal_positive = [row for row in positives if row["region_type"] == "TERMINAL"]
        spanned = [row for row in candidates if int(row.get("n_unique_exact_spanning_qnames") or 0) > 0]
        spanned_i = [row for row in candidates if int(row.get("n_unique_spanning_with_I_qnames") or 0) > 0]
        spanned_d = [row for row in candidates if int(row.get("n_unique_spanning_with_D_qnames") or 0) > 0]
        spanned_n = [row for row in candidates if int(row.get("n_unique_spanning_with_N_qnames") or 0) > 0]
        spanned_any_idn = [row for row in candidates if int(row.get("n_unique_spanning_with_any_IDN_qnames") or 0) > 0]
        result.append(
            [
                mag,
                len(members),
                sum(lengths[contig] for contig in members),
                "Y" if internal_permissive else "N",
                "Y" if terminal_permissive else "N",
                "Y" if permissive_positives else "N",
                "Y" if internal_positive else "N",
                "Y" if terminal_positive else "N",
                clipping_analyzable_call,
                clipping_call,
                ("Y" if stringent_positives else "N") if clipping_analyzable_call == "Y" else "NA",
                len(internal_permissive),
                len(terminal_permissive),
                len(permissive_positives),
                len(internal_positive),
                len(terminal_positive),
                len(positives),
                len(stringent_positives),
                len(internal_candidates),
                len(terminal_candidates),
                len(candidates),
                len(stringent),
                len(spanned),
                len(spanned_i),
                len(spanned_d),
                len(spanned_n),
                len(spanned_any_idn),
                len({row["contig"] for row in candidates}),
                len({row["contig"] for row in stringent}),
                len({row["contig"] for row in positives}),
                len({row["contig"] for row in stringent_positives}),
                joined_contigs(candidates),
                joined_contigs(stringent),
                joined_contigs(positives),
                joined_contigs(stringent_positives),
                max((int(row["coverage"]) for row in records), default=0),
                max((int(row["clipping_count"]) for row in records), default=0),
                max((int(row.get("n_unique_all_clipped_qnames") or 0) for row in records), default=0),
                f"{max((float(row['clipping_ratio']) for row in records), default=0):.6f}",
                "Y" if candidates else "N",
                "Y" if stringent else "N",
            ]
        )
    return result


def contig_rows(
    lengths: Mapping[str, int], contig_to_mag: Mapping[str, str],
    records: Sequence[Mapping[str, str]],
) -> List[List[object]]:
    """Reuse the MAG summary logic with each contig as its own group."""
    by_contig: Dict[str, List[Mapping[str, str]]] = defaultdict(list)
    for row in records:
        by_contig[row["contig"]].append(row)
    singletons = {contig: [contig] for contig in lengths}
    keep = [index for index, field in enumerate(SUMMARY_HEADER) if field not in CONTIG_DROPPED]
    result: List[List[object]] = []
    for row in summary_rows(lengths, singletons, by_contig):
        contig = row[0]
        result.append([contig, contig_to_mag.get(contig, "NA"), lengths[contig],
                       *(row[index] for index in keep)])
    return result


def main() -> None:
    args = parse_args()
    for path, label in (
        (args.clipping, "clipping table"),
        (args.mag_contig_map, "MAG-contig map"),
        (args.fai, "assembly FAI"),
    ):
        require_nonempty(path, label)
    lengths = load_fai(args.fai)
    contig_to_mag, mag_to_contigs = load_mag_map(args.mag_contig_map, lengths)
    fields, records = load_clipping(args.clipping, lengths)

    positives = [row for row in records if row["Clipping"] == "Y"]
    mapped = [
        (contig_to_mag[row["contig"]], row)
        for row in records
        if row["contig"] in contig_to_mag
    ]
    mapped.sort(key=lambda item: (item[0], item[1]["contig"], int(item[1]["pos_1based"])))
    mapped_positive = [item for item in mapped if item[1]["Clipping"] == "Y"]
    by_mag: Dict[str, List[Mapping[str, str]]] = defaultdict(list)
    for mag, row in mapped:
        by_mag[mag].append(row)
    summaries = summary_rows(lengths, mag_to_contigs, by_mag)

    atomic_write(args.positive_site_output, fields, ([row[field] for field in fields] for row in positives))
    atomic_write(args.mag_site_output, ("MAG_ID", *fields), ([mag, *(row[field] for field in fields)] for mag, row in mapped))
    atomic_write(args.mag_positive_site_output, ("MAG_ID", *fields), ([mag, *(row[field] for field in fields)] for mag, row in mapped_positive))
    atomic_write(args.summary_output, SUMMARY_HEADER, summaries)
    clipping_index = SUMMARY_HEADER.index("Clipping")
    atomic_write(args.positive_summary_output, SUMMARY_HEADER, (row for row in summaries if row[clipping_index] == "Y"))
    contigs = contig_rows(lengths, contig_to_mag, records)
    atomic_write(args.contig_summary_output, CONTIG_HEADER, contigs)

    print(f"MAG summary rows: {len(summaries)}")
    print(f"Contig summary rows: {len(contigs)}")
    print(f"All-assembly Clipping=Y sites: {len(positives)}")
    print(f"MAG-associated Clipping=Y sites: {len(mapped_positive)}")
    print(f"MAGs with Clipping_ge1=Y: {sum(row[SUMMARY_HEADER.index('Clipping_ge1')] == 'Y' for row in summaries)}")
    print(f"MAGs with Internal_clipping=Y: {sum(row[SUMMARY_HEADER.index('Internal_clipping')] == 'Y' for row in summaries)}")
    print(f"MAGs with Terminal_clipping=Y: {sum(row[SUMMARY_HEADER.index('Terminal_clipping')] == 'Y' for row in summaries)}")
    print(f"MAGs with Clipping=Y: {sum(row[clipping_index] == 'Y' for row in summaries)}")
    print(f"MAGs with Stringent_clipping=Y: {sum(row[SUMMARY_HEADER.index('Stringent_clipping')] == 'Y' for row in summaries)}")


if __name__ == "__main__":
    main()
