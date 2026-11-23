#!/usr/bin/env python3
"""Clipping tables for all assemblers: site table, per-assembler density and contig/MAG summaries.

Subcommands:
  sites              Merge the assembler site tables and classify each site as INTERNAL or TERMINAL
  assembler          Terminal versus internal site density per assembler
  mag                Add MAG_ID to the site table and write the contig and MAG summaries
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mapping-manifest", required=True)
    parser.add_argument("--assemblies-fai", required=True)
    parser.add_argument("--terminal-window", type=int, default=100)
    parser.add_argument("--input", action="append", required=True,
                        help="assembler=TSV; may be supplied repeatedly")
    parser.add_argument("--all-output", required=True)
    parser.add_argument("--internal-output", required=True)
    return parser.parse_args()


def load_fai(path: str) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    with open(path, encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            fields = raw.rstrip("\n\r").split("\t")
            if len(fields) < 2:
                fail(f"malformed FAI row {line_number}: {path}")
            name = fields[0]
            try:
                length = int(fields[1])
            except ValueError as error:
                fail(f"invalid FAI length at row {line_number}: {error}")
            if not name or length < 1 or name in lengths:
                fail(f"blank/duplicate FAI sequence ID: {name}")
            lengths[name] = length
    if not lengths:
        fail(f"FAI contains no sequences: {path}")
    return lengths


def load_assemblers(manifest: str) -> Tuple[Dict[str, str], Dict[str, Dict[str, int]]]:
    assembler_to_fasta: Dict[str, str] = {}
    assembler_lengths: Dict[str, Dict[str, int]] = {}
    with open(manifest, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Assembler", "Assembly_FASTA"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            fail(f"mapping manifest lacks columns: {sorted(missing)}")
        for line_number, row in enumerate(reader, start=2):
            assembler = (row.get("Assembler") or "").strip()
            fasta = (row.get("Assembly_FASTA") or "").strip()
            if not assembler or not fasta or assembler in assembler_to_fasta:
                fail(f"invalid/duplicate mapping manifest row {line_number}")
            fai = fasta + ".fai"
            if not os.path.isfile(fai) or os.path.getsize(fai) == 0:
                fail(f"assembler FAI missing for {assembler}: {fai}")
            assembler_to_fasta[assembler] = fasta
            assembler_lengths[assembler] = load_fai(fai)
    if not assembler_to_fasta:
        fail("mapping manifest contains no assembler rows")
    return assembler_to_fasta, assembler_lengths


def parse_inputs(values: Sequence[str]) -> Dict[str, str]:
    result: Dict[str, str] = {}
    for value in values:
        if "=" not in value:
            fail(f"--input must be assembler=TSV: {value}")
        assembler, path = value.split("=", 1)
        assembler, path = assembler.strip(), path.strip()
        if not assembler or not path or assembler in result:
            fail(f"invalid/duplicate --input: {value}")
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            fail(f"input table missing or empty: {path}")
        result[assembler] = path
    return result


def atomic_write(path: str, fields: Sequence[str], rows: Iterable[Mapping[str, str]]) -> None:
    destination = Path(path).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="", delete=False,
        dir=str(destination.parent), prefix=destination.name + ".", suffix=".tmp",
    )
    try:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
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


def _sites_main() -> None:
    args = parse_args()
    if args.terminal_window < 1:
        fail("terminal-window must be >=1")
    if not os.path.isfile(args.mapping_manifest) or os.path.getsize(args.mapping_manifest) == 0:
        fail(f"mapping manifest missing: {args.mapping_manifest}")
    if not os.path.isfile(args.assemblies_fai) or os.path.getsize(args.assemblies_fai) == 0:
        fail(f"all-assemblies FAI missing: {args.assemblies_fai}")

    assembler_to_fasta, assembler_lengths = load_assemblers(args.mapping_manifest)
    inputs = parse_inputs(args.input)
    if set(inputs) != set(assembler_to_fasta):
        missing = sorted(set(assembler_to_fasta) - set(inputs))
        extra = sorted(set(inputs) - set(assembler_to_fasta))
        fail(f"input assemblers do not match mapping manifest; missing={missing}, extra={extra}")

    combined = load_fai(args.assemblies_fai)
    union: Dict[str, Tuple[str, int]] = {}
    for assembler, lengths in assembler_lengths.items():
        for contig, length in lengths.items():
            if contig in union:
                fail(f"contig ID occurs in multiple assemblers: {contig}")
            union[contig] = (assembler, length)
    if set(union) != set(combined):
        missing = sorted(set(combined) - set(union))[:10]
        extra = sorted(set(union) - set(combined))[:10]
        fail(
            "assembler contig union differs from the all-assemblies FAI; "
            f"missing_from_assemblers={missing}, extra_in_assemblers={extra}"
        )
    for contig, length in combined.items():
        if union[contig][1] != length:
            fail(f"assembler/all-assemblies contig length mismatch: {contig}")

    base_fields: List[str] | None = None
    rows: List[Dict[str, str]] = []
    seen = set()
    for assembler in sorted(inputs):
        path = inputs[assembler]
        with open(path, encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = list(reader.fieldnames or [])
            required = {"contig", "length", "anvio_pos_0based", "pos_1based"}
            missing = required - set(fields)
            if missing:
                fail(f"{path} lacks columns: {sorted(missing)}")
            if base_fields is None:
                base_fields = fields
            elif fields != base_fields:
                fail(f"input clipping schemas differ: {path}")
            for line_number, raw in enumerate(reader, start=2):
                contig = (raw.get("contig") or "").strip()
                try:
                    length = int(raw["length"])
                    pos0 = int(raw["anvio_pos_0based"])
                    pos1 = int(raw["pos_1based"])
                except (TypeError, ValueError) as error:
                    fail(f"invalid clipping coordinate at {path}:{line_number}: {error}")
                if contig not in assembler_lengths[assembler]:
                    fail(f"clipping contig is not in the {assembler} assembly: {contig}")
                if assembler_lengths[assembler][contig] != length or not 0 <= pos0 < length or pos1 != pos0 + 1:
                    fail(f"stale clipping coordinate: {assembler}:{contig}:{pos1}")
                key = (contig, pos0)
                if key in seen:
                    fail(f"duplicate clipping position across assembler tables: {contig}:{pos1}")
                seen.add(key)
                distance = min(pos0, length - 1 - pos0)
                region = "TERMINAL" if distance < args.terminal_window else "INTERNAL"
                row = {"Assembler": assembler}
                row.update(raw)
                row["distance_to_nearest_end_bp"] = str(distance)
                row["region_type"] = region
                row["reference_boundary"] = "Y" if distance == 0 else "N"
                rows.append(row)

    if base_fields is None:
        fail("no clipping table schema was loaded")
    preferred = [
        "Assembler", "contig", "length", "pos_1based",
        "region_type", "distance_to_nearest_end_bp", "reference_boundary",
        "Clipping_analyzable", "Clipping_analysis_reason",
        "Clipping_ge1", "Clipping", "Stringent_clipping",
        "default_candidate", "stringent_candidate",
        "validated_clipping", "clipping_confidence",
    ]
    all_fields = ["Assembler", *base_fields,
                  "distance_to_nearest_end_bp", "region_type", "reference_boundary"]
    output_fields = []
    for field in preferred + all_fields:
        if field in all_fields and field not in output_fields:
            output_fields.append(field)
    rows.sort(key=lambda row: (row["Assembler"], row["contig"], int(row["pos_1based"])))
    internal = [row for row in rows if row["region_type"] == "INTERNAL"]

    atomic_write(args.all_output, output_fields, rows)
    atomic_write(args.internal_output, output_fields, internal)

    print(f"Merged clipping rows: {len(rows)}")
    print(f"Internal rows: {len(internal)}")
    print(f"Terminal rows: {len(rows) - len(internal)}")


def cmd_sites(argv):
    """Merge the assembler site tables and classify each site as INTERNAL or TERMINAL"""
    sys.argv = [sys.argv[0] + " sites", *argv]
    return _sites_main()


def _assembler_parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mapping-manifest", required=True)
    parser.add_argument("--terminal-window", type=int, required=True)
    parser.add_argument("--sites", required=True, help="clipping site table (all regions)")
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def assembler_lengths(manifest: str) -> Dict[str, Dict[str, int]]:
    result: Dict[str, Dict[str, int]] = {}
    with open(manifest, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            assembler = (row.get("Assembler") or "").strip()
            fasta = (row.get("Assembly_FASTA") or "").strip()
            if not assembler or not fasta or assembler in result:
                fail("invalid mapping manifest")
            result[assembler] = load_fai(fasta + ".fai")
    if not result:
        fail("mapping manifest contains no assemblies")
    return result


def load_rows(path: str) -> list[dict[str, str]]:
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Assembler", "contig", "region_type"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            fail(f"{path} lacks columns: {sorted(missing)}")
        return [dict(row) for row in reader]


def _assembler_atomic_write(path: str, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    destination = Path(path).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="", delete=False,
        dir=str(destination.parent), prefix=destination.name + ".", suffix=".tmp",
    )
    try:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
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


def density(count: int, bp: int) -> float:
    return count * 1_000_000.0 / bp if bp > 0 else math.nan


def fmt(value: float) -> str:
    return "NA" if math.isnan(value) else f"{value:.6f}"


def ratio(numerator: float, denominator: float) -> str:
    if math.isnan(numerator) or math.isnan(denominator) or denominator <= 0:
        return "NA"
    return f"{numerator / denominator:.6f}"


METRICS = (
    ("candidate_sites", "candidate"),
    ("stringent_sites", "stringent"),
    ("Clipping_ge1_Y_sites", "Clipping_ge1"),
    ("validated_clipping_sites", "validated"),
    ("Clipping_Y_sites", "Clipping"),
    ("candidate_clipped_alignment_records", "clipped_records"),
)


def _assembler_main() -> None:
    args = _assembler_parse_args()
    if args.terminal_window < 1:
        fail("terminal-window must be >=1")
    for path in (args.mapping_manifest, args.sites):
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            fail(f"required input missing or empty: {path}")

    lengths = assembler_lengths(args.mapping_manifest)
    analyzed_bp: Dict[Tuple[str, str], int] = defaultdict(int)
    for assembler, contigs in lengths.items():
        for length in contigs.values():
            terminal_bp = min(length, 2 * args.terminal_window)
            analyzed_bp[(assembler, "TERMINAL")] += terminal_bp
            analyzed_bp[(assembler, "INTERNAL")] += length - terminal_bp

    counts: Dict[Tuple[str, str], Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for row in load_rows(args.sites):
        key = (row["Assembler"], row["region_type"])
        c = counts[key]
        c["candidate"] += 1
        try:
            c["clipped_records"] += int(row.get("clipping_count") or 0)
        except ValueError:
            fail("invalid clipping_count in clipping site table")
        if (row.get("stringent_candidate") or "").strip().upper() == "Y":
            c["stringent"] += 1
        for field, metric in (("Clipping_ge1", "Clipping_ge1"), ("validated_clipping", "validated"),
                              ("Clipping", "Clipping")):
            if (row.get(field) or "").strip().upper() == "Y":
                c[metric] += 1

    assemblers = sorted(lengths)
    fields = ["Assembler", "INTERNAL_analyzed_bp", "TERMINAL_analyzed_bp"]
    for name, _ in METRICS:
        fields += [f"INTERNAL_{name}", f"INTERNAL_{name}_per_Mbp",
                   f"TERMINAL_{name}", f"TERMINAL_{name}_per_Mbp",
                   f"terminal_to_internal_{name}_density_ratio"]
    rows = []
    for assembler in [*assemblers, "ALL"]:
        group = assemblers if assembler == "ALL" else [assembler]
        out: Dict[str, object] = {"Assembler": assembler}
        bp = {region: sum(analyzed_bp[(a, region)] for a in group) for region in ("INTERNAL", "TERMINAL")}
        out["INTERNAL_analyzed_bp"], out["TERMINAL_analyzed_bp"] = bp["INTERNAL"], bp["TERMINAL"]
        for name, metric in METRICS:
            dens = {}
            for region in ("INTERNAL", "TERMINAL"):
                value = sum(counts[(a, region)][metric] for a in group)
                dens[region] = density(value, bp[region])
                out[f"{region}_{name}"] = value
                out[f"{region}_{name}_per_Mbp"] = fmt(dens[region])
            out[f"terminal_to_internal_{name}_density_ratio"] = ratio(dens["TERMINAL"], dens["INTERNAL"])
        rows.append(out)
    _assembler_atomic_write(args.output, fields, rows)

    all_row = rows[-1]
    print("Terminal/internal density ratios (ALL assemblers):")
    for name, _ in METRICS:
        print(f"  {name}: {all_row[f'terminal_to_internal_{name}_density_ratio']}")


def cmd_assembler(argv):
    """Terminal versus internal site density per assembler"""
    sys.argv = [sys.argv[0] + " assembler", *argv]
    return _assembler_main()


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

CONTIG_DROPPED = {
    "MAG_ID", "n_contigs", "total_length_bp",
    "n_contigs_with_candidate", "n_contigs_with_stringent_candidate",
    "n_contigs_with_clipping", "n_contigs_with_stringent_clipping",
    "candidate_contigs", "stringent_candidate_contigs",
    "clipping_contigs", "stringent_clipping_contigs",
}

CONTIG_HEADER = ("contig", "MAG_ID", "length_bp",
                 *(field for field in SUMMARY_HEADER if field not in CONTIG_DROPPED))


def _mag_parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--clipping", required=True)
    parser.add_argument("--mag-contig-map", required=True)
    parser.add_argument("--fai", required=True)
    parser.add_argument("--site-output", required=True,
                        help="site table rewritten with a MAG_ID column (NA for unbinned contigs)")
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--contig-summary-output", required=True,
                        help="One row per assembly contig (all contigs in the FAI)")
    return parser.parse_args()


def require_nonempty(path: str, label: str) -> None:
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        fail(f"{label} not found or empty: {path}")


def _mag_atomic_write(path: str, header: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
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
        fields = [field for field in (reader.fieldnames or []) if field != "MAG_ID"]
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
        "Assembler", "contig", "length", "pos_1based", "region_type",
        "Clipping_analyzable", "Clipping", "Clipping_ge1", "Stringent_clipping",
    ]
    preferred = [field for field in preferred if field in fields]
    fields = [*preferred, *(field for field in fields if field not in preferred)]
    records.sort(key=lambda row: (row.get("Assembler", ""), row["contig"], int(row["pos_1based"])))
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


def _mag_main() -> None:
    args = _mag_parse_args()
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
    mapped_positive = [row for row in positives if row["contig"] in contig_to_mag]
    by_mag: Dict[str, List[Mapping[str, str]]] = defaultdict(list)
    for row in records:
        if row["contig"] in contig_to_mag:
            by_mag[contig_to_mag[row["contig"]]].append(row)
    summaries = summary_rows(lengths, mag_to_contigs, by_mag)
    clipping_index = SUMMARY_HEADER.index("Clipping")

    site_fields = list(fields)
    site_fields.insert(site_fields.index("contig"), "MAG_ID")
    _mag_atomic_write(args.site_output, site_fields, (
        [contig_to_mag.get(row["contig"], "NA") if field == "MAG_ID" else row[field] for field in site_fields]
        for row in records))
    _mag_atomic_write(args.summary_output, SUMMARY_HEADER, summaries)
    contigs = contig_rows(lengths, contig_to_mag, records)
    _mag_atomic_write(args.contig_summary_output, CONTIG_HEADER, contigs)

    print(f"MAG summary rows: {len(summaries)}")
    print(f"Contig summary rows: {len(contigs)}")
    print(f"All-assembly Clipping=Y sites: {len(positives)}")
    print(f"MAG-associated Clipping=Y sites: {len(mapped_positive)}")
    print(f"MAGs with Clipping_ge1=Y: {sum(row[SUMMARY_HEADER.index('Clipping_ge1')] == 'Y' for row in summaries)}")
    print(f"MAGs with Internal_clipping=Y: {sum(row[SUMMARY_HEADER.index('Internal_clipping')] == 'Y' for row in summaries)}")
    print(f"MAGs with Terminal_clipping=Y: {sum(row[SUMMARY_HEADER.index('Terminal_clipping')] == 'Y' for row in summaries)}")
    print(f"MAGs with Clipping=Y: {sum(row[clipping_index] == 'Y' for row in summaries)}")
    print(f"MAGs with Stringent_clipping=Y: {sum(row[SUMMARY_HEADER.index('Stringent_clipping')] == 'Y' for row in summaries)}")


def cmd_mag(argv):
    """Add MAG_ID to the site table and write the contig and MAG summaries"""
    sys.argv = [sys.argv[0] + " mag", *argv]
    return _mag_main()


COMMANDS = {
    'sites': cmd_sites,
    'assembler': cmd_assembler,
    'mag': cmd_mag,
}

if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] not in COMMANDS:
        sys.exit(__doc__)
    result = COMMANDS[sys.argv[1]](sys.argv[2:])
    sys.exit(result if isinstance(result, int) else 0)
