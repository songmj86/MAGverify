#!/usr/bin/env python3
"""Reclassify taxonomic-shift regions using validated clipping as the final chimera criterion.

Final policy
------------
Taxonomic_chimera = Y only when a taxonomic-shift region contains at least one
validated clipping site from 03_clipping.sh (Clipping=Y).

High_confidence_taxonomic_chimera = Y when the same region contains at least one
Stringent_clipping=Y site (Clipping=Y, anvi'o coverage >=10x, clipping ratio 1.0).

The pre-existing primary+supplementary HiFi structural fields are retained only
as audit annotations and do not determine the final Taxonomic_chimera call.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple

NA_VALUES = {"", "NA", "N/A", "NONE", "."}


def fail(msg: str) -> None:
    raise SystemExit(f"ERROR: {msg}")


def parse_args() -> argparse.Namespace:
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
    p.add_argument("--release", default="GTDB_R232")
    return p.parse_args()


def read_tsv(path: str) -> Tuple[List[str], List[Dict[str, str]]]:
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        fail(f"missing or empty input: {path}")
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = [dict(row) for row in reader]
    if not fields:
        fail(f"missing header: {path}")
    return fields, rows


def atomic_write(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, str]]) -> None:
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
    fields, rows = read_tsv(path)
    required = {"contig", "pos_1based", "Clipping", "Stringent_clipping"}
    missing = required - set(fields)
    if missing:
        fail(f"clipping input lacks columns: {sorted(missing)}")
    out: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    seen = set()
    for row in rows:
        clip = (row.get("Clipping") or "").strip().upper()
        stringent = (row.get("Stringent_clipping") or "").strip().upper()
        if clip not in {"Y", "N"}:
            fail(f"invalid Clipping flag at {row.get('contig','NA')}:{row.get('pos_1based','NA')}")
        if stringent not in {"Y", "N", "NA"}:
            fail(f"invalid Stringent_clipping flag at {row.get('contig','NA')}:{row.get('pos_1based','NA')}")
        if stringent == "Y" and clip != "Y":
            fail(f"Stringent_clipping=Y without Clipping=Y at {row.get('contig','NA')}:{row.get('pos_1based','NA')}")
        if clip != "Y":
            continue
        contig = (row.get("contig") or "").strip()
        try:
            pos = int(row.get("pos_1based") or "")
        except ValueError:
            fail(f"invalid clipping position: {contig}:{row.get('pos_1based','NA')}")
        key = (contig, pos)
        if key in seen:
            fail(f"duplicate validated clipping site: {contig}:{pos}")
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
        fail(f"invalid transition interval: {row.get('event_id','NA')}")
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
        fail(f"non-taxonomic-shift row in region table: {row.get('event_id','NA')}")
    validated = (row.get("validated_clipping_support") or "").strip().upper()
    stringent = (row.get("Stringent_clipping_support") or "").strip().upper()
    if validated not in {"Y", "N"}:
        fail(f"invalid validated_clipping_support: {row.get('event_id','NA')}/{validated}")
    if stringent not in {"Y", "N"}:
        fail(f"invalid Stringent_clipping_support: {row.get('event_id','NA')}/{stringent}")
    if stringent == "Y" and validated != "Y":
        fail(f"Stringent_clipping_support=Y without validated clipping: {row.get('event_id','NA')}")

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
        fail(f"invalid Taxonomic_shift on contig {row.get('contig','NA')}: {tax_shift}")

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


def main() -> None:
    args = parse_args()
    region_fields, region_rows = read_tsv(args.regions_in)
    clipping = load_validated_clipping(args.clipping)
    contig_fields, contig_rows = read_tsv(args.contigs_in)
    mag_contig_fields, mag_contig_rows = read_tsv(args.mag_contigs_in)
    mag_fields, mag_rows = read_tsv(args.mags_in)

    required_regions = {
        "event_id", "MAG_ID", "contig", "taxonomic_transition", "analysis_pos_1based",
        "transition_interval_start_1based", "transition_interval_end_1based",
    }
    missing = required_regions - set(region_fields)
    if missing:
        fail(f"region input lacks columns: {sorted(missing)}")
    required_contigs = {"MAG_ID", "contig", "length", "Taxonomic_shift", "Taxonomic_chimera_assessable"}
    missing = required_contigs - set(contig_fields)
    if missing:
        fail(f"contig input lacks columns: {sorted(missing)}")
    if set(mag_contig_fields) != set(contig_fields):
        missing_mc = required_contigs - set(mag_contig_fields)
        if missing_mc:
            fail(f"MAG-contig input lacks columns: {sorted(missing_mc)}")
    if "MAG_ID" not in mag_fields:
        fail("MAG input lacks MAG_ID")

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
            fail(f"MAG contig absent from assembly-wide contig table: {c}")
        row = dict(contig_index[c])
        classified_mag_contigs.append(row)
        mag = row["MAG_ID"]
        if mag in NA_VALUES:
            fail(f"MAG-contig table has blank MAG_ID: {c}")
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
    rel = args.release
    all_regions_path = out_dir / f"taxonomic_transition_clipping_validation.{rel}.regions.tsv"
    chim_regions_path = out_dir / f"taxonomic_chimera.{rel}.regions.tsv"
    high_regions_path = out_dir / f"high_confidence_taxonomic_chimera.{rel}.regions.tsv"
    contig_path = out_dir / "chimera.contig.tsv"
    mag_contig_path = out_dir / "chimera.MAG_contig.tsv"
    mag_path = out_dir / "chimera.MAG.tsv"
    meta_path = out_dir / "run_metadata.tsv"

    chim_regions = [r for r in classified_regions if r["Taxonomic_chimera"] == "Y"]
    high_regions = [r for r in classified_regions if r["High_confidence_taxonomic_chimera"] == "Y"]

    atomic_write(all_regions_path, new_region_fields, classified_regions)
    atomic_write(chim_regions_path, new_region_fields, chim_regions)
    atomic_write(high_regions_path, new_region_fields, high_regions)
    atomic_write(contig_path, new_contig_fields, classified_contigs)
    atomic_write(mag_contig_path, new_contig_fields, classified_mag_contigs)
    atomic_write(mag_path, new_mag_fields, classified_mags)

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


if __name__ == "__main__":
    main()
