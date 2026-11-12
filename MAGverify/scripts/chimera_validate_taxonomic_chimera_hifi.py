#!/usr/bin/env python3
"""Convert taxonomic-transition candidates into HiFi-supported structural-chimera calls.

The upstream GTDB evaluator is used only to identify local taxonomic transitions across
all assembly contigs and provide GC/depth context. This script performs the structural decision using primary and
supplementary HiFi alignments, deduplicated by QNAME. Secondary alignments are not used
for structural validation.

A taxonomic transition is called STRUCTURAL_CHIMERA only when both sides have sufficient
HiFi support, fewer than the minimum number of anchored spanning reads cross the full
uncertain transition interval, and systematic breakpoint evidence is present. Existing
validated clipping from 03_clipping.sh is reused and is never recomputed. Taxonomic
transitions with >= minimum anchored spanning QNAMEs are HIFI_CONTINUITY_SUPPORTED.
Insufficient or mixed evidence is UNRESOLVED.

SCG redundancy is annotation only. It is not a structural call criterion because SCG
copy-number anomalies do not localize a breakpoint by themselves.
"""
from __future__ import annotations

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import re
import subprocess
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, Mapping, Sequence, Set, Tuple

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
MATCH_OPS = {"M", "=", "X"}
REF_OPS = {"M", "D", "N", "=", "X"}
CLIP_OPS = {"S", "H"}


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def parse_args() -> argparse.Namespace:
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


def require_file(path: str | None, label: str, optional: bool = False) -> None:
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


def atomic_write(path: str, fields: Sequence[str], rows: Iterable[Mapping[str, str]]) -> None:
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


def read_mag_map(path: str) -> dict[str, str]:
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
            require_file(p, label)
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
    source = (row.get("source") or row.get("hmm_source") or "").strip()
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
    """Annotate SCG redundancy only for MAG-associated contigs.

    Structural validation is assembly-wide. Non-MAG contigs therefore receive NA
    for MAG-domain-dependent SCG fields rather than being excluded from analysis.
    """
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
    """Return the localized uncertain transition interval.

    The upstream detector clusters adjacent gene-boundary calls that share the same
    local taxonomic transition. Those clustered boundary coordinates are a safer
    structural interval than the full span between the nearest classified genes,
    which can become wide when unresolved genes are skipped for voting.
    """
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


def validate_mode(
    bam: str, events: list[dict[str, str]], indices: list[int], lengths: Mapping[str, int],
    anchor: int, threads: int,
) -> dict[int, dict[str, Set[str]]]:
    stats: dict[int, dict[str, Set[str]]] = defaultdict(lambda: {"left": set(), "right": set(), "spanning": set(), "clipped": set()})
    if not indices: return stats
    by_contig: dict[str, list[int]] = defaultdict(list)
    for idx in indices: by_contig[events[idx]["contig"]].append(idx)
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", suffix=".bed", delete=False) as bed:
        bed_path = bed.name
        for idx in indices:
            row = events[idx]; contig = row["contig"]; length = lengths[contig]
            left, right, _, _ = event_bounds(row, length)
            start0 = max(0, left - anchor - 1)
            end0 = min(length, right + anchor)
            bed.write(f"{contig}\t{start0}\t{end0}\n")
    try:
        cmd = ["samtools", "view", "-@", str(threads), "-F", "4", "-L", bed_path, bam]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, encoding="utf-8", errors="replace")
        assert proc.stdout is not None and proc.stderr is not None
        for raw in proc.stdout:
            f = raw.rstrip("\n\r").split("\t")
            if len(f) < 11 or f[5] == "*": continue
            qname, contig, cigar = f[0], f[2], f[5]
            if contig not in by_contig: continue
            try: ref_start0 = int(f[3]) - 1
            except ValueError: continue
            ops = cigar_ops(cigar)
            ref_end0 = ref_start0 + sum(n for n, op in ops if op in REF_OPS)
            runs = matched_runs(ref_start0, ops)
            for idx in by_contig[contig]:
                row = events[idx]; length = lengths[contig]
                left1, right1, lo1, hi1 = event_bounds(row, length)
                search_start0 = max(0, left1 - anchor - 1); search_end0 = min(length, right1 + anchor)
                if ref_end0 <= search_start0 or ref_start0 >= search_end0: continue
                left_start0, left_end0 = max(0, left1 - anchor), left1
                right_start0, right_end0 = right1 - 1, min(length, right1 - 1 + anchor)
                lbp = matched_bp(runs, left_start0, left_end0)
                rbp = matched_bp(runs, right_start0, right_end0)
                if lbp >= anchor: stats[idx]["left"].add(qname)
                if rbp >= anchor: stats[idx]["right"].add(qname)
                if lbp >= anchor and rbp >= anchor: stats[idx]["spanning"].add(qname)
                left_clip = terminal_clip(ops, True) > 0 and (lo1 - 1) <= (ref_start0 + 1) <= (hi1 + 1)
                right_clip = terminal_clip(ops, False) > 0 and (lo1 - 1) <= ref_end0 <= (hi1 + 1)
                if left_clip or right_clip: stats[idx]["clipped"].add(qname)
        proc.stdout.close(); stderr = proc.stderr.read(); code = proc.wait()
        if code != 0: fail(f"samtools view failed for {bam}: {stderr}")
    finally:
        try: os.unlink(bed_path)
        except FileNotFoundError: pass
    return stats


def evidence_label(call: str, high_conf: bool) -> str:
    if call == "Y": return "HIGH_CONFIDENCE_HIFI_STRUCTURAL_CHIMERA" if high_conf else "HIFI_STRUCTURAL_CHIMERA"
    if call == "N": return "HIFI_CONTINUITY_SUPPORTED"
    return "UNRESOLVED_HIFI_STRUCTURAL_EVIDENCE"


def main() -> None:
    args = parse_args()
    if args.anchor_bp < 1 or args.min_support_qnames < 1 or args.threads < 1 or args.parallel_bam_jobs < 1:
        fail("anchor/min-support/threads/parallel-bam-jobs must be positive")
    for p, label in ((args.candidate_events, "candidate events"), (args.candidate_contigs, "candidate contigs"), (args.mapping_manifest, "mapping manifest"), (args.mag_map, "MAG-contig map")):
        require_file(p, label)
    require_file(args.clipping, "validated clipping", optional=True)
    require_file(args.scg_tsv, "SCG table", optional=True)
    require_file(args.mag_taxonomy_tsv, "MAG GTDB-Tk taxonomy table", optional=True)

    event_fields, events = read_tsv(args.candidate_events)
    contig_fields, contig_rows = read_tsv(args.candidate_contigs)
    required_event = {"event_id", "MAG_ID", "contig", "left_block_end_1based", "right_block_start_1based", "taxonomic_transition", "shift_rank"}
    if required_event - set(event_fields): fail(f"candidate event table lacks columns: {sorted(required_event-set(event_fields))}")
    required_contig = {"MAG_ID", "contig", "length", "Taxonomic_shift"}
    if required_contig - set(contig_fields): fail(f"candidate contig table lacks columns: {sorted(required_contig-set(contig_fields))}")

    lengths = {r["contig"]: int(r["length"]) for r in contig_rows}
    target_contigs = set(lengths)
    contig_to_mag = read_mag_map(args.mag_map)
    routes = load_routes(args.mapping_manifest, target_contigs)
    clipping = load_clipping(args.clipping)
    target_mags = {contig_to_mag[c] for c in target_contigs if c in contig_to_mag}
    mag_domains = read_mag_domains(args.mag_taxonomy_tsv, target_mags) if args.scg_tsv and target_mags else {}
    scg = load_scg(args.scg_tsv, target_contigs, contig_to_mag, mag_domains)

    # Identify existing validated clipping that lies anywhere in the uncertain interval
    # between the nearest classified left and right taxonomic blocks.
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
        # Backward-compatible names now reflect structurally validated calls, not raw taxonomy shifts.
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

    atomic_write(args.events_output, final_event_fields, final_events)
    atomic_write(args.structural_events_output, final_event_fields, structural_events)
    atomic_write(args.high_confidence_events_output, final_event_fields, high_conf_events)
    atomic_write(args.unresolved_events_output, final_event_fields, unresolved_events)
    atomic_write(args.contig_output, final_contig_fields, final_contigs)
    atomic_write(args.structural_contig_output, final_contig_fields, structural_contigs)
    atomic_write(args.high_confidence_contig_output, final_contig_fields, high_conf_contigs)
    atomic_write(args.unresolved_contig_output, final_contig_fields, unresolved_contigs)
    atomic_write(args.scg_audit_output, final_contig_fields, scg_audit)

    print(f"Taxonomic-transition candidates: {len(final_events)}")
    print(f"HiFi-supported structural chimeras: {len(structural_events)}")
    print(f"High-confidence structural chimeras (Stringent_clipping=Y): {len(high_conf_events)}")
    print(f"HiFi continuity-supported transitions: {sum(r['Structural_chimera']=='N' for r in final_events)}")
    print(f"Unresolved transitions: {len(unresolved_events)}")
    print(f"Contigs with any SCG family redundancy (annotation only): {len(scg_audit)}")

if __name__ == "__main__":
    main()
