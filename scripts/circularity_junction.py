#!/usr/bin/env python3
"""Circularity END|START junction evidence for one assembler group.

Subcommands:
  summarize          Count HiFi reads reaching and spanning each END|START junction
  integrate          Join candidate selection, rotation and junction counts into the group table
"""

import collections
import csv
import re
import subprocess
import sys


def load_manifest(path):
    data = {}
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        index = {name: position for position, name in enumerate(header)}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            contig = fields[index["Contig_ID"]]
            data[contig] = {
                "mag": fields[index["MAG_ID"]],
                "length": int(fields[index["Contig_length_bp"]]),
                "junction0": int(fields[index["Rotated_END_START_left_base"]]),
            }
    return data


def reference_end(start0, operations):
    position = start0
    for length, operation in operations:
        if operation in "MDN=X":
            position += length
    return position


def aligned_reference_bases_each_side(start0, operations, junction0):
    """Count actual M/=X-aligned reference bases on each side of the junction."""
    position = start0
    left_aligned = 0
    right_aligned = 0
    for length, operation in operations:
        if operation in "M=X":
            block_start = position
            block_end = position + length
            if block_start < junction0:
                left_aligned += max(0, min(block_end, junction0) - block_start)
            if block_end > junction0:
                right_aligned += max(0, block_end - max(block_start, junction0))
            position = block_end
        elif operation in "DN":
            position += length
    return left_aligned, right_aligned


def junction_crossing_events(start0, operations, junction0):
    """Return I/D/N events that cross the internalized END|START boundary."""
    position = start0
    events = []
    for length, operation in operations:
        if operation in "M=X":
            position += length
        elif operation == "I":
            if position == junction0:
                events.append((operation, length))
        elif operation in "DN":
            next_position = position + length
            if position < junction0 < next_position:
                events.append((operation, length))
            position = next_position
    return events


def cmd_summarize(argv):
    """Count HiFi reads reaching and spanning each END|START junction"""
    (
        bam_path,
        junction_bed,
        manifest_path,
        anchor_text,
        minimum_hifi_reads_text,
        span_details_path,
        summary_path,
    ) = argv

    anchor = int(anchor_text)

    minimum_hifi_reads = int(minimum_hifi_reads_text)

    cigar_pattern = re.compile(r"(\d+)([MIDNSHP=X])")

    def cigar_operations(cigar):
        if cigar == "*":
            return []
        operations = [(int(length), operation) for length, operation in cigar_pattern.findall(cigar)]
        if "".join(f"{length}{operation}" for length, operation in operations) != cigar:
            raise ValueError(f"Unsupported CIGAR string: {cigar}")
        return operations

    manifest = load_manifest(manifest_path)

    ps_spans = collections.defaultdict(set)

    ps_end_reaching = collections.defaultdict(set)

    ps_start_reaching = collections.defaultdict(set)

    all_spans = collections.defaultdict(set)

    all_end_reaching = collections.defaultdict(set)

    all_start_reaching = collections.defaultdict(set)

    primary_spans = collections.defaultdict(set)

    supplementary_spans = collections.defaultdict(set)

    secondary_spans = collections.defaultdict(set)

    ps_clean_spans = collections.defaultdict(set)

    all_clean_spans = collections.defaultdict(set)

    ps_hifi_i_qnames = collections.defaultdict(set)

    ps_hifi_d_qnames = collections.defaultdict(set)

    ps_hifi_n_qnames = collections.defaultdict(set)

    ps_hifi_i_lengths = collections.defaultdict(list)

    ps_hifi_d_lengths = collections.defaultdict(list)

    ps_hifi_n_lengths = collections.defaultdict(list)

    hifi_i_qnames = collections.defaultdict(set)

    hifi_d_qnames = collections.defaultdict(set)

    hifi_n_qnames = collections.defaultdict(set)

    hifi_i_lengths = collections.defaultdict(list)

    hifi_d_lengths = collections.defaultdict(list)

    hifi_n_lengths = collections.defaultdict(list)

    span_rows = {}

    command = ["samtools", "view", "-L", junction_bed, bam_path]

    process = subprocess.Popen(command, stdout=subprocess.PIPE, text=True)

    assert process.stdout is not None

    for raw in process.stdout:
        fields = raw.rstrip("\n").split("\t")
        if len(fields) < 11:
            continue
        qname = fields[0]
        flag = int(fields[1])
        contig = fields[2]
        if contig not in manifest or flag & 0x4:
            continue
        mapq = int(fields[4])
        operations = cigar_operations(fields[5])
        if not operations:
            continue
        start0 = int(fields[3]) - 1
        end0 = reference_end(start0, operations)
        junction0 = manifest[contig]["junction0"]
        left_aligned, right_aligned = aligned_reference_bases_each_side(start0, operations, junction0)
        has_left_anchor = left_aligned >= anchor
        has_right_anchor = right_aligned >= anchor
        crossing_events = junction_crossing_events(start0, operations, junction0)
        has_crossing_dn = any(event_type in {"D", "N"} for event_type, _ in crossing_events)
        junction_span = has_left_anchor and has_right_anchor and not has_crossing_dn
        aln_type = "secondary" if flag & 0x100 else ("supplementary" if flag & 0x800 else "primary")
        is_ps = aln_type in {"primary", "supplementary"}

        if has_left_anchor:
            all_end_reaching[contig].add(qname)
            if is_ps:
                ps_end_reaching[contig].add(qname)
        if has_right_anchor:
            all_start_reaching[contig].add(qname)
            if is_ps:
                ps_start_reaching[contig].add(qname)

        for event_type, event_length in crossing_events:
            if event_type == "I":
                hifi_i_qnames[contig].add(qname)
                hifi_i_lengths[contig].append(event_length)
                if is_ps:
                    ps_hifi_i_qnames[contig].add(qname)
                    ps_hifi_i_lengths[contig].append(event_length)
            elif event_type == "D":
                hifi_d_qnames[contig].add(qname)
                hifi_d_lengths[contig].append(event_length)
                if is_ps:
                    ps_hifi_d_qnames[contig].add(qname)
                    ps_hifi_d_lengths[contig].append(event_length)
            elif event_type == "N":
                hifi_n_qnames[contig].add(qname)
                hifi_n_lengths[contig].append(event_length)
                if is_ps:
                    ps_hifi_n_qnames[contig].add(qname)
                    ps_hifi_n_lengths[contig].append(event_length)

        if junction_span:
            all_spans[contig].add(qname)
            if is_ps:
                ps_spans[contig].add(qname)
            if aln_type == "primary":
                primary_spans[contig].add(qname)
            elif aln_type == "supplementary":
                supplementary_spans[contig].add(qname)
            else:
                secondary_spans[contig].add(qname)

            clean_alignment = not crossing_events
            if clean_alignment:
                all_clean_spans[contig].add(qname)
                if is_ps:
                    ps_clean_spans[contig].add(qname)

            key = (contig, qname, aln_type)
            row = (
                manifest[contig]["mag"], contig, qname, mapq, start0 + 1, end0,
                left_aligned, right_aligned, "Y" if clean_alignment else "N",
                sum(1 for t, _ in crossing_events if t == "I"),
                sum(1 for t, _ in crossing_events if t == "D"),
                sum(1 for t, _ in crossing_events if t == "N"),
                sum(n for t, n in crossing_events if t == "I"),
                sum(n for t, n in crossing_events if t == "D"),
                sum(n for t, n in crossing_events if t == "N"),
                max((n for t, n in crossing_events if t == "I"), default=0),
                max((n for t, n in crossing_events if t == "D"), default=0),
                max((n for t, n in crossing_events if t == "N"), default=0),
                aln_type, fields[5],
            )
            previous = span_rows.get(key)
            if previous is None or (row[8] == "Y" and previous[8] != "Y"):
                span_rows[key] = row

    return_code = process.wait()

    if return_code != 0:
        raise SystemExit(f"ERROR: samtools view failed with exit code {return_code}")

    with open(span_details_path, "w") as out:
        out.write(
            "MAG_ID\tContig_ID\tRead_name\tMAPQ\tAlignment_start\tAlignment_end\t"
            "Actual_aligned_bp_original_END_side\tActual_aligned_bp_original_START_side\t"
            "CLEAN_JUNCTION_ALIGNMENT_Y_N\tJUNCTION_I_EVENT_COUNT\tJUNCTION_D_EVENT_COUNT\tJUNCTION_N_EVENT_COUNT\t"
            "TOTAL_JUNCTION_I_BP\tTOTAL_JUNCTION_D_BP\tTOTAL_JUNCTION_N_BP\t"
            "MAX_JUNCTION_I_BP\tMAX_JUNCTION_D_BP\tMAX_JUNCTION_N_BP\tAlignment_type\tCIGAR\n"
        )
        for row in sorted(span_rows.values()):
            out.write("\t".join(map(str, row)) + "\n")

    with open(summary_path, "w") as out:
        out.write(
            "MAG_ID\tContig_ID\t"
            "HIFI_PRIMARY_SUPPLEMENTARY_END_REACHING_UNIQUE_READS\t"
            "HIFI_PRIMARY_SUPPLEMENTARY_START_REACHING_UNIQUE_READS\t"
            "HIFI_PRIMARY_SUPPLEMENTARY_ANALYZABLE_Y_N\t"
            "HIFI_PRIMARY_SUPPLEMENTARY_END_START_SPANNING_UNIQUE_READS\t"
            "HIFI_PRIMARY_SUPPLEMENTARY_END_START_SUPPORTED_Y_N\t"
            "HIFI_ALL_MAPPED_END_REACHING_UNIQUE_READS\t"
            "HIFI_ALL_MAPPED_START_REACHING_UNIQUE_READS\t"
            "HIFI_ALL_MAPPED_ANALYZABLE_Y_N\t"
            "HIFI_ALL_MAPPED_END_START_SPANNING_UNIQUE_READS\t"
            "HIFI_ALL_MAPPED_END_START_SUPPORTED_Y_N\t"
            "HIFI_PRIMARY_END_START_SPANNING_UNIQUE_READS\t"
            "HIFI_SUPPLEMENTARY_END_START_SPANNING_UNIQUE_READS\t"
            "HIFI_SECONDARY_END_START_SPANNING_UNIQUE_READS\t"
            "HIFI_SECONDARY_ADDED_END_START_SUPPORT_Y_N\t"
            "HIFI_PRIMARY_SUPPLEMENTARY_CLEAN_END_START_SPANNING_UNIQUE_READS\t"
            "HIFI_PRIMARY_SUPPLEMENTARY_CLEAN_JUNCTION_Y_N\t"
            "HIFI_ALL_MAPPED_CLEAN_END_START_SPANNING_UNIQUE_READS\t"
            "HIFI_ALL_MAPPED_CLEAN_JUNCTION_Y_N\t"
            "HIFI_PRIMARY_SUPPLEMENTARY_JUNCTION_I_UNIQUE_READS\t"
            "HIFI_PRIMARY_SUPPLEMENTARY_JUNCTION_I_EVENT_COUNT\t"
            "HIFI_PRIMARY_SUPPLEMENTARY_JUNCTION_I_TOTAL_BP\t"
            "HIFI_PRIMARY_SUPPLEMENTARY_MAX_JUNCTION_I_BP\t"
            "HIFI_PRIMARY_SUPPLEMENTARY_JUNCTION_D_UNIQUE_READS\t"
            "HIFI_PRIMARY_SUPPLEMENTARY_JUNCTION_D_EVENT_COUNT\t"
            "HIFI_PRIMARY_SUPPLEMENTARY_JUNCTION_D_TOTAL_BP\t"
            "HIFI_PRIMARY_SUPPLEMENTARY_MAX_JUNCTION_D_BP\t"
            "HIFI_PRIMARY_SUPPLEMENTARY_JUNCTION_N_UNIQUE_READS\t"
            "HIFI_PRIMARY_SUPPLEMENTARY_JUNCTION_N_EVENT_COUNT\t"
            "HIFI_PRIMARY_SUPPLEMENTARY_JUNCTION_N_TOTAL_BP\t"
            "HIFI_PRIMARY_SUPPLEMENTARY_MAX_JUNCTION_N_BP\t"
            "HIFI_JUNCTION_I_UNIQUE_READS\tHIFI_JUNCTION_D_UNIQUE_READS\tHIFI_JUNCTION_N_UNIQUE_READS\t"
            "HIFI_MAX_JUNCTION_I_BP\tHIFI_MAX_JUNCTION_D_BP\tHIFI_MAX_JUNCTION_N_BP\t"
            "HIFI_END_START_SUPPORTED\n"
        )
        for contig in sorted(manifest, key=lambda item: (manifest[item]["mag"], item)):
            ps_end = len(ps_end_reaching.get(contig, set()))
            ps_start = len(ps_start_reaching.get(contig, set()))
            all_end = len(all_end_reaching.get(contig, set()))
            all_start = len(all_start_reaching.get(contig, set()))
            ps_an = ps_end >= minimum_hifi_reads and ps_start >= minimum_hifi_reads
            all_an = all_end >= minimum_hifi_reads and all_start >= minimum_hifi_reads
            ps_count = len(ps_spans.get(contig, set()))
            all_count = len(all_spans.get(contig, set()))
            ps_call = "Y" if ps_count >= minimum_hifi_reads else ("N" if ps_an else "NA")
            all_call = "Y" if all_count >= minimum_hifi_reads else ("N" if all_an else "NA")
            secondary_added = (
                "Y" if all_call == "Y" and ps_call != "Y"
                else ("N" if ps_call in {"Y", "N"} and all_call in {"Y", "N"} else "NA")
            )
            ps_clean_count = len(ps_clean_spans.get(contig, set()))
            all_clean_count = len(all_clean_spans.get(contig, set()))
            ps_clean_call = "Y" if ps_clean_count >= minimum_hifi_reads else ("N" if ps_an else "NA")
            all_clean_call = "Y" if all_clean_count >= minimum_hifi_reads else ("N" if all_an else "NA")
            out.write(
                f"{manifest[contig]['mag']}\t{contig}\t"
                f"{ps_end}\t{ps_start}\t{'Y' if ps_an else 'N'}\t{ps_count}\t{ps_call}\t"
                f"{all_end}\t{all_start}\t{'Y' if all_an else 'N'}\t{all_count}\t{all_call}\t"
                f"{len(primary_spans.get(contig, set()))}\t"
                f"{len(supplementary_spans.get(contig, set()))}\t"
                f"{len(secondary_spans.get(contig, set()))}\t{secondary_added}\t"
                f"{ps_clean_count}\t{ps_clean_call}\t{all_clean_count}\t{all_clean_call}\t"
                f"{len(ps_hifi_i_qnames.get(contig, set()))}\t{len(ps_hifi_i_lengths.get(contig, []))}\t{sum(ps_hifi_i_lengths.get(contig, []))}\t{max(ps_hifi_i_lengths.get(contig, [0]))}\t"
                f"{len(ps_hifi_d_qnames.get(contig, set()))}\t{len(ps_hifi_d_lengths.get(contig, []))}\t{sum(ps_hifi_d_lengths.get(contig, []))}\t{max(ps_hifi_d_lengths.get(contig, [0]))}\t"
                f"{len(ps_hifi_n_qnames.get(contig, set()))}\t{len(ps_hifi_n_lengths.get(contig, []))}\t{sum(ps_hifi_n_lengths.get(contig, []))}\t{max(ps_hifi_n_lengths.get(contig, [0]))}\t"
                f"{len(hifi_i_qnames.get(contig, set()))}\t{len(hifi_d_qnames.get(contig, set()))}\t{len(hifi_n_qnames.get(contig, set()))}\t"
                f"{max(hifi_i_lengths.get(contig, [0]))}\t{max(hifi_d_lengths.get(contig, [0]))}\t{max(hifi_n_lengths.get(contig, [0]))}\t"
                f"{ps_call}\n"
            )


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def index_by_contig(rows, label):
    out = {}
    for row in rows:
        contig = row.get("Contig_ID", "")
        if not contig:
            continue
        if contig in out:
            raise SystemExit(f"ERROR: duplicate {label} row for {contig}")
        out[contig] = row
    return out


def cmd_integrate(argv):
    """Join candidate selection, rotation and junction counts into the group table"""
    selection_path, manifest_path, hifi_summary_path, output_path = argv

    selection = read_tsv(selection_path)

    manifest = index_by_contig(read_tsv(manifest_path), "rotation manifest")

    hifi = index_by_contig(read_tsv(hifi_summary_path), "HiFi summary")

    hifi_fields = list(next(iter(hifi.values())).keys()) if hifi else []

    hifi_fields = [f for f in hifi_fields if f not in {"MAG_ID", "Contig_ID"}]

    output_header = [
        "MAG_ID", "Number_of_contigs", "Contig_ID", "Contig_length_bp",
        "SINGLE_CONTIG", "CIRCULARITY_ANALYZABLE", "ROTATION_VALIDATED",
        "Longest_exact_terminal_overlap_bp_up_to_100000",
    ] + hifi_fields

    if "HIFI_END_START_SUPPORTED" not in output_header:
        output_header.append("HIFI_END_START_SUPPORTED")

    output_header.append("CIRCULARITY_SUPPORTED")

    with open(output_path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=output_header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for selected in selection:
            mag_id = selected["MAG_ID"]
            n = selected["Number_of_contigs"]
            if selected["SINGLE_CONTIG_MAG"] != "Y":
                row = {field: "NA" for field in output_header}
                row.update({"MAG_ID": mag_id, "Number_of_contigs": n, "Contig_ID": "NA",
                            "SINGLE_CONTIG": "N", "CIRCULARITY_ANALYZABLE": "N",
                            "HIFI_PRIMARY_SUPPLEMENTARY_ANALYZABLE_Y_N": "N",
                            "HIFI_ALL_MAPPED_ANALYZABLE_Y_N": "N",
                            "CIRCULARITY_SUPPORTED": "NA"})
                writer.writerow(row); continue
            contig = selected["Contig_ID"]
            if contig not in manifest or contig not in hifi:
                raise SystemExit(f"ERROR: missing HiFi circularity evidence for {mag_id}/{contig}")
            meta = manifest[contig]; h = hifi[contig]
            status = h["HIFI_PRIMARY_SUPPLEMENTARY_END_START_SUPPORTED_Y_N"]
            analyzable = h["HIFI_PRIMARY_SUPPLEMENTARY_ANALYZABLE_Y_N"]
            row = {field: "NA" for field in output_header}
            row.update({
                "MAG_ID": mag_id, "Number_of_contigs": n, "Contig_ID": contig,
                "Contig_length_bp": meta["Contig_length_bp"], "SINGLE_CONTIG": "Y",
                "CIRCULARITY_ANALYZABLE": analyzable,
                "ROTATION_VALIDATED": meta["ROTATION_VALIDATED"],
                "Longest_exact_terminal_overlap_bp_up_to_100000": meta["Longest_exact_terminal_overlap_bp_up_to_100000"],
            })
            for field in hifi_fields:
                row[field] = h.get(field, "NA")
            row["HIFI_END_START_SUPPORTED"] = status
            row["CIRCULARITY_SUPPORTED"] = status
            writer.writerow(row)


COMMANDS = {
    'summarize': cmd_summarize,
    'integrate': cmd_integrate,
}

if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] not in COMMANDS:
        sys.exit(__doc__)
    result = COMMANDS[sys.argv[1]](sys.argv[2:])
    sys.exit(result if isinstance(result, int) else 0)
