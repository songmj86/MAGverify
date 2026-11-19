#!/usr/bin/env python3
"""Normalize exact clipping positions reported by anvi-script-find-misassemblies.

The anvi'o `pos` column is a zero-based reference-base coordinate. This script
preserves it as `anvio_pos_0based` and additionally reports `pos_1based` for
human-readable tables. Every input row is a `default_candidate`. A
`stringent_candidate` is the exact-position subset with the configured minimum
anvi'o alignment-record depth and clipping ratio.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence

OUTPUT_FIELDS = (
    "contig",
    "length",
    "anvio_pos_0based",
    "pos_1based",
    "relative_pos",
    "coverage",
    "clipping_count",
    "clipping_ratio",
    "default_candidate",
    "stringent_candidate",
    "candidate_clipping",
)


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--anvio-clipping", required=True)
    parser.add_argument("--fai", required=True)
    parser.add_argument("--default-output", required=True)
    parser.add_argument("--stringent-output", required=True)
    parser.add_argument("--stringent-min-depth", type=int, default=10)
    parser.add_argument("--stringent-ratio", type=float, default=1.0)
    return parser.parse_args()


def atomic_write(path: str, rows: Iterable[Mapping[str, object]]) -> None:
    destination = Path(path).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="",
        delete=False,
        dir=str(destination.parent),
        prefix=destination.name + ".",
        suffix=".tmp",
    )
    try:
        writer = csv.DictWriter(
            handle, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\n"
        )
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
                fail(f"invalid or duplicate FAI sequence: {fields[0]}")
            lengths[fields[0]] = length
    if not lengths:
        fail("FAI contains no sequences")
    return lengths


def resolve_columns(fields: Sequence[str]) -> Dict[str, str]:
    normalized = {field.strip().lower(): field for field in fields}
    aliases = {
        "contig": ("contig",),
        "length": ("length",),
        "pos": ("pos", "position"),
        "relative_pos": ("relative_pos", "relative_position"),
        "cov": ("cov", "coverage"),
        "clipping": ("clipping", "clipping_coverage", "clipping_count"),
        "clipping_ratio": ("clipping_ratio",),
    }
    result: Dict[str, str] = {}
    for label, choices in aliases.items():
        for choice in choices:
            if choice in normalized:
                result[label] = normalized[choice]
                break
        if label not in result:
            fail(f"anvi'o clipping table lacks column for {label}")
    return result


def main() -> None:
    args = parse_args()
    if args.stringent_min_depth < 1:
        fail("stringent minimum depth must be at least 1")
    if not 0 < args.stringent_ratio <= 1:
        fail("stringent clipping ratio must be in (0,1]")
    for path, label in ((args.anvio_clipping, "anvi'o clipping table"), (args.fai, "FAI")):
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            fail(f"{label} is missing or empty: {path}")

    lengths = load_fai(args.fai)
    rows: List[Dict[str, object]] = []
    seen = set()
    with open(args.anvio_clipping, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = resolve_columns(reader.fieldnames or [])
        for line_number, raw in enumerate(reader, start=2):
            contig = (raw.get(columns["contig"]) or "").strip()
            try:
                length = int(raw[columns["length"]])
                pos0 = int(raw[columns["pos"]])
                relative_pos = float(raw[columns["relative_pos"]])
                coverage = int(raw[columns["cov"]])
                clipping = int(raw[columns["clipping"]])
                ratio = float(raw[columns["clipping_ratio"]])
            except (TypeError, ValueError) as error:
                fail(f"invalid numeric value in anvi'o row {line_number}: {error}")
            if contig not in lengths or lengths[contig] != length:
                fail(f"anvi'o contig/length is stale at row {line_number}: {contig}")
            if not 0 <= pos0 < length or (contig, pos0) in seen:
                fail(f"invalid or duplicate anvi'o position: {contig}:{pos0}")
            seen.add((contig, pos0))
            if coverage < 1 or clipping < 1 or clipping > coverage:
                fail(f"invalid clipping counts at {contig}:{pos0}")
            expected_ratio = clipping / coverage
            if not math.isclose(ratio, expected_ratio, rel_tol=0, abs_tol=5.0e-7):
                fail(f"inconsistent clipping ratio at {contig}:{pos0}")
            expected_relative = pos0 / length
            if not math.isclose(relative_pos, expected_relative, rel_tol=0, abs_tol=5.0e-7):
                fail(f"inconsistent relative position at {contig}:{pos0}")
            stringent = coverage >= args.stringent_min_depth and math.isclose(
                ratio, args.stringent_ratio, rel_tol=0, abs_tol=5.0e-7
            )
            rows.append(
                {
                    "contig": contig,
                    "length": length,
                    "anvio_pos_0based": pos0,
                    "pos_1based": pos0 + 1,
                    "relative_pos": f"{relative_pos:.12g}",
                    "coverage": coverage,
                    "clipping_count": clipping,
                    "clipping_ratio": f"{ratio:.6f}",
                    "default_candidate": "Y",
                    "stringent_candidate": "Y" if stringent else "N",
                    "candidate_clipping": "Y",
                }
            )

    rows.sort(key=lambda row: (str(row["contig"]), int(row["anvio_pos_0based"])))
    stringent_rows = [row for row in rows if row["stringent_candidate"] == "Y"]
    atomic_write(args.default_output, rows)
    atomic_write(args.stringent_output, stringent_rows)
    print(f"anvi'o default candidates: {len(rows)}")
    print(f"Stringent candidates: {len(stringent_rows)}")


if __name__ == "__main__":
    main()
