#!/usr/bin/env python3
"""Summarize terminal versus internal clipping density per analyzed base."""

from __future__ import annotations

import argparse
import csv
import math
import os
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, Mapping, Sequence, Tuple


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mapping-manifest", required=True)
    parser.add_argument("--terminal-window", type=int, required=True)
    parser.add_argument("--default", required=True)
    parser.add_argument("--stringent", required=True)
    parser.add_argument("--final", required=True)
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--ratio-output", required=True)
    return parser.parse_args()


def load_fai(path: str) -> Dict[str, int]:
    result: Dict[str, int] = {}
    with open(path, encoding="utf-8") as handle:
        for raw in handle:
            fields = raw.rstrip("\n\r").split("\t")
            if len(fields) >= 2:
                name, length = fields[0], int(fields[1])
                if name in result:
                    fail(f"duplicate contig in FAI: {name}")
                result[name] = length
    return result


def source_lengths(manifest: str) -> Dict[str, Dict[str, int]]:
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
        required = {"source_assembler", "contig", "region_type"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            fail(f"{path} lacks columns: {sorted(missing)}")
        return [dict(row) for row in reader]


def atomic_write(path: str, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
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


def main() -> None:
    args = parse_args()
    if args.terminal_window < 1:
        fail("terminal-window must be >=1")
    for path in (args.mapping_manifest, args.default, args.stringent, args.final):
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            fail(f"required input missing or empty: {path}")

    lengths = source_lengths(args.mapping_manifest)
    analyzed_bp: Dict[Tuple[str, str], int] = defaultdict(int)
    for assembler, contigs in lengths.items():
        for length in contigs.values():
            terminal_bp = min(length, 2 * args.terminal_window)
            internal_bp = length - terminal_bp
            analyzed_bp[(assembler, "TERMINAL")] += terminal_bp
            analyzed_bp[(assembler, "INTERNAL")] += internal_bp

    default_rows = load_rows(args.default)
    stringent_rows = load_rows(args.stringent)
    final_rows = load_rows(args.final)

    counts: Dict[Tuple[str, str], Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for row in default_rows:
        key = (row["source_assembler"], row["region_type"])
        counts[key]["candidate"] += 1
        try:
            counts[key]["candidate_clipped_alignment_records"] += int(row.get("clipping_count") or 0)
        except ValueError:
            fail("invalid clipping_count in default table")
    for row in stringent_rows:
        counts[(row["source_assembler"], row["region_type"])]["stringent"] += 1
    for row in final_rows:
        key = (row["source_assembler"], row["region_type"])
        if (row.get("Clipping_ge1") or "").strip().upper() == "Y":
            counts[key]["clipping_ge1_y"] += 1
        if (row.get("validated_clipping") or "").strip().upper() == "Y":
            counts[key]["validated"] += 1
        if (row.get("Clipping") or "").strip().upper() == "Y":
            counts[key]["clipping_y"] += 1

    assemblers = sorted(lengths)
    summary_rows = []
    metrics_by_group: Dict[str, Dict[str, Dict[str, float]]] = defaultdict(dict)
    for assembler in [*assemblers, "ALL"]:
        for region in ("INTERNAL", "TERMINAL"):
            if assembler == "ALL":
                bp = sum(analyzed_bp[(a, region)] for a in assemblers)
                c = defaultdict(int)
                for a in assemblers:
                    for metric, value in counts[(a, region)].items():
                        c[metric] += value
            else:
                bp = analyzed_bp[(assembler, region)]
                c = counts[(assembler, region)]
            candidate = c["candidate"]
            stringent = c["stringent"]
            clipping_ge1_y = c["clipping_ge1_y"]
            validated = c["validated"]
            clipping_y = c["clipping_y"]
            clipped_records = c["candidate_clipped_alignment_records"]
            metrics = {
                "candidate": density(candidate, bp),
                "stringent": density(stringent, bp),
                "clipping_ge1_y": density(clipping_ge1_y, bp),
                "validated": density(validated, bp),
                "clipping_y": density(clipping_y, bp),
                "candidate_clipped_alignment_records": density(clipped_records, bp),
            }
            metrics_by_group[assembler][region] = metrics
            summary_rows.append({
                "source_assembler": assembler,
                "region_type": region,
                "analyzed_bp": bp,
                "candidate_sites": candidate,
                "candidate_sites_per_Mbp": fmt(metrics["candidate"]),
                "stringent_sites": stringent,
                "stringent_sites_per_Mbp": fmt(metrics["stringent"]),
                "Clipping_ge1_Y_sites": clipping_ge1_y,
                "Clipping_ge1_Y_sites_per_Mbp": fmt(metrics["clipping_ge1_y"]),
                "validated_clipping_sites": validated,
                "validated_clipping_sites_per_Mbp": fmt(metrics["validated"]),
                "Clipping_Y_sites": clipping_y,
                "Clipping_Y_sites_per_Mbp": fmt(metrics["clipping_y"]),
                "candidate_clipped_alignment_records": clipped_records,
                "candidate_clipped_alignment_records_per_Mbp": fmt(metrics["candidate_clipped_alignment_records"]),
            })

    summary_fields = [
        "source_assembler", "region_type", "analyzed_bp",
        "candidate_sites", "candidate_sites_per_Mbp",
        "stringent_sites", "stringent_sites_per_Mbp",
        "Clipping_ge1_Y_sites", "Clipping_ge1_Y_sites_per_Mbp",
        "validated_clipping_sites", "validated_clipping_sites_per_Mbp",
        "Clipping_Y_sites", "Clipping_Y_sites_per_Mbp",
        "candidate_clipped_alignment_records", "candidate_clipped_alignment_records_per_Mbp",
    ]
    atomic_write(args.summary_output, summary_fields, summary_rows)

    ratio_rows = []
    for assembler in [*assemblers, "ALL"]:
        internal = metrics_by_group[assembler]["INTERNAL"]
        terminal = metrics_by_group[assembler]["TERMINAL"]
        ratio_rows.append({
            "source_assembler": assembler,
            "terminal_to_internal_candidate_density_ratio": ratio(terminal["candidate"], internal["candidate"]),
            "terminal_to_internal_stringent_density_ratio": ratio(terminal["stringent"], internal["stringent"]),
            "terminal_to_internal_Clipping_ge1_Y_density_ratio": ratio(terminal["clipping_ge1_y"], internal["clipping_ge1_y"]),
            "terminal_to_internal_validated_density_ratio": ratio(terminal["validated"], internal["validated"]),
            "terminal_to_internal_Clipping_Y_density_ratio": ratio(terminal["clipping_y"], internal["clipping_y"]),
            "terminal_to_internal_candidate_clipped_alignment_record_density_ratio": ratio(
                terminal["candidate_clipped_alignment_records"], internal["candidate_clipped_alignment_records"]
            ),
        })
    ratio_fields = list(ratio_rows[0])
    atomic_write(args.ratio_output, ratio_fields, ratio_rows)

    all_ratio = next(row for row in ratio_rows if row["source_assembler"] == "ALL")
    print("Terminal/internal density ratios (ALL assemblers):")
    for key, value in all_ratio.items():
        if key != "source_assembler":
            print(f"  {key}: {value}")


if __name__ == "__main__":
    main()
