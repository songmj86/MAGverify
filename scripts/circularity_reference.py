#!/usr/bin/env python3
"""Circularity reference: match MAG contigs to their assembler and build the rotated reference.

Subcommands:
  match              Match single-contig MAGs to the assembler that produced them
  rotate             Build the assembler reference with rotated single-contig MAG contigs
"""

import gzip
import hashlib
import os
import sys
from collections import defaultdict
from pathlib import Path


def open_text(path):
    return gzip.open(path, "rt") if str(path).lower().endswith(".gz") else open(path, "rt")


def fasta_records(path):
    with open_text(path) as handle:
        header = None
        sequence_parts = []
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence_parts).upper()
                header = line[1:]
                sequence_parts = []
            else:
                if header is None:
                    raise SystemExit(f"ERROR: Sequence before FASTA header in {path}")
                sequence_parts.append(line)
        if header is not None:
            yield header, "".join(sequence_parts).upper()


def fasta_id(header):
    return header.split()[0]


def md5_sequence(sequence):
    return hashlib.md5(sequence.encode("ascii")).hexdigest()


def cmd_match(argv):
    """Match single-contig MAGs to the assembler that produced them"""
    (
        assembly_dir,
        assembly_fasta_names,
        mag_dir,
        grouped_mag_root,
        plan_path,
        groups_path,
        validation_path,
        hifi_anchor_text,
    ) = argv

    hifi_anchor = int(hifi_anchor_text)

    allowed_suffixes = (
        ".fa", ".fna", ".fasta", ".fas",
        ".fa.gz", ".fna.gz", ".fasta.gz", ".fas.gz",
    )

    def mag_id_from_path(path):
        name = Path(path).name
        lower = name.lower()
        for suffix in sorted(allowed_suffixes, key=len, reverse=True):
            if lower.endswith(suffix):
                return name[:-len(suffix)]
        return name

    assembly_names = assembly_fasta_names.split()

    if not assembly_names:
        raise SystemExit("ERROR: No assembler FASTA names were configured")

    assemblies = []

    seen_labels = set()

    for name in assembly_names:
        path = os.path.realpath(os.path.join(assembly_dir, name))
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            raise SystemExit(f"ERROR: Assembler FASTA missing or empty: {path}")
        label = Path(name).name
        if label.lower().endswith(".gz"):
            label = label[:-3]
        label = Path(label).stem
        if label in seen_labels:
            raise SystemExit(f"ERROR: Duplicate assembler FASTA label: {label}")
        seen_labels.add(label)
        assemblies.append((label, path))

    mag_files = []

    for path in sorted(Path(mag_dir).iterdir()):
        if path.is_file() and path.name.lower().endswith(allowed_suffixes):
            mag_files.append(path)

    if not mag_files:
        raise SystemExit(f"ERROR: No MAG FASTA files found in {mag_dir}")

    eligible = {}

    candidate_by_contig = {}

    validation_rows = []

    for mag_path in mag_files:
        mag_id = mag_id_from_path(mag_path)
        if mag_id in eligible:
            raise SystemExit(f"ERROR: Duplicate MAG_ID derived from filenames: {mag_id}")
        records = list(fasta_records(mag_path))
        if not records:
            raise SystemExit(f"ERROR: Empty MAG FASTA: {mag_path}")
        if len(records) != 1:
            validation_rows.append((mag_id, len(records), "NA", "NA", "NA", "N"))
            continue
        header, sequence = records[0]
        contig_id = fasta_id(header)
        if len(sequence) < 2 * hifi_anchor + 2:
            raise SystemExit(
                f"ERROR: {mag_id} / {contig_id} is too short for {hifi_anchor}-bp anchors"
            )
        if contig_id in candidate_by_contig:
            raise SystemExit(
                f"ERROR: Contig ID occurs in two single-contig MAGs: {contig_id}"
            )
        eligible[mag_id] = {
            "path": str(mag_path.resolve()),
            "contig": contig_id,
            "length": len(sequence),
            "md5": md5_sequence(sequence),
        }
        candidate_by_contig[contig_id] = mag_id

    if not eligible:
        raise SystemExit("ERROR: No single-contig MAGs were found")

    matches = defaultdict(list)

    for label, assembly_path in assemblies:
        seen_ids = set()
        for header, sequence in fasta_records(assembly_path):
            contig_id = fasta_id(header)
            if contig_id in seen_ids:
                raise SystemExit(
                    f"ERROR: Duplicate contig ID in assembly {assembly_path}: {contig_id}"
                )
            seen_ids.add(contig_id)
            mag_id = candidate_by_contig.get(contig_id)
            if mag_id is None:
                continue
            candidate = eligible[mag_id]
            if len(sequence) == candidate["length"] and md5_sequence(sequence) == candidate["md5"]:
                matches[mag_id].append((label, assembly_path))

    for mag_id, item in sorted(eligible.items()):
        source_matches = matches.get(mag_id, [])
        if not source_matches:
            raise SystemExit(
                "ERROR: Single-contig MAG contig was not found by exact ID and sequence in any "
                f"configured assembly: {mag_id} / {item['contig']}"
            )
        if len(source_matches) > 1:
            sources = ", ".join(label for label, _ in source_matches)
            raise SystemExit(
                f"ERROR: MAG contig matches multiple assemblies: {mag_id} / {sources}"
            )
        item["assembler"], item["assembly"] = source_matches[0]

    Path(grouped_mag_root).mkdir(parents=True, exist_ok=True)

    group_counts = defaultdict(int)

    plan_rows = []

    for mag_id, item in sorted(eligible.items()):
        group_id = item["assembler"]
        group_dir = Path(grouped_mag_root) / group_id
        group_dir.mkdir(parents=True, exist_ok=True)
        link_path = group_dir / Path(item["path"]).name
        if link_path.exists() or link_path.is_symlink():
            raise SystemExit(f"ERROR: Duplicate grouped MAG path: {link_path}")
        os.symlink(item["path"], link_path)
        group_counts[(group_id, item["assembly"])] += 1
        plan_rows.append(
            (
                group_id,
                item["assembler"],
                item["assembly"],
                mag_id,
                item["path"],
                item["contig"],
                item["length"],
                item["md5"],
                "Y",
            )
        )
        validation_rows.append(
            (
                mag_id,
                1,
                item["contig"],
                item["assembler"],
                item["assembly"],
                "Y",
            )
        )

    with open(plan_path, "w") as out:
        out.write(
            "Group_ID\tAssembler\tAssembly_FASTA\tMAG_ID\tMAG_FASTA\tContig_ID\t"
            "Contig_length_bp\tSequence_MD5\tASSEMBLER_EXACT_MATCH\n"
        )
        for row in plan_rows:
            out.write("\t".join(map(str, row)) + "\n")

    with open(groups_path, "w") as out:
        out.write("Group_ID\tAssembler\tAssembly_FASTA\tMAG_count\tGrouped_MAG_directory\n")
        for (group_id, assembly_path), count in sorted(group_counts.items()):
            out.write(
                f"{group_id}\t{group_id}\t{assembly_path}\t{count}\t"
                f"{(Path(grouped_mag_root) / group_id).resolve()}\n"
            )

    with open(validation_path, "w") as out:
        out.write(
            "MAG_ID\tNumber_of_contigs\tContig_ID\tAssembler\tAssembly_FASTA\t"
            "ASSEMBLER_EXACT_MATCH\n"
        )
        for row in sorted(validation_rows, key=lambda value: value[0]):
            out.write("\t".join(map(str, row)) + "\n")

    print(f"Single-contig MAGs matched to an assembler FASTA: {len(plan_rows)}")

    print(f"Assembler groups: {len(group_counts)}")


def mag_id_from_path(path):
    name = Path(path).name
    if name.endswith(".gz"):
        name = name[:-3]
    for suffix in (".fasta", ".fna", ".fa", ".fas"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return name


def longest_exact_terminal_overlap(sequence, maximum_checked_bp=100000):
    """Return the longest exact prefix/suffix match within a bounded audit window."""
    checked_bp = min(maximum_checked_bp, len(sequence) // 2)
    if checked_bp == 0:
        return 0
    combined = sequence[:checked_bp] + "\0" + sequence[-checked_bp:]
    prefix_lengths = [0] * len(combined)
    for position in range(1, len(combined)):
        matched = prefix_lengths[position - 1]
        while matched and combined[position] != combined[matched]:
            matched = prefix_lengths[matched - 1]
        if combined[position] == combined[matched]:
            matched += 1
        prefix_lengths[position] = matched
    return prefix_lengths[-1]


def write_fasta_record(handle, header, sequence, width=80):
    handle.write(f">{header}\n")
    for start in range(0, len(sequence), width):
        handle.write(sequence[start : start + width] + "\n")


def cmd_rotate(argv):
    """Build the assembler reference with rotated single-contig MAG contigs"""
    (
        assembly_path,
        mag_dir,
        output_reference,
        manifest_path,
        selection_path,
        candidate_bed_path,
        junction_bed_path,
        hifi_anchor,
    ) = argv

    hifi_anchor = int(hifi_anchor)

    allowed_suffixes = (".fa", ".fna", ".fasta", ".fas", ".fa.gz", ".fna.gz", ".fasta.gz", ".fas.gz")

    mag_files = sorted(
        str(path)
        for path in Path(mag_dir).iterdir()
        if path.is_file() and path.name.lower().endswith(allowed_suffixes)
    )

    if not mag_files:
        raise SystemExit(f"ERROR: No MAG FASTA files found in {mag_dir}")

    candidates = {}

    selection_rows = []

    for mag_file in mag_files:
        mag_id = mag_id_from_path(mag_file)
        records = list(fasta_records(mag_file))
        if len(records) != 1:
            selection_rows.append((mag_id, len(records), "NA", "N"))
            continue
        header, sequence = records[0]
        contig_id = fasta_id(header)
        if not sequence:
            raise SystemExit(f"ERROR: Empty sequence in {mag_file}")
        if contig_id in candidates:
            other = candidates[contig_id]["mag_id"]
            raise SystemExit(
                f"ERROR: Contig ID {contig_id} occurs in two MAG files: {other} and {mag_id}"
            )
        candidates[contig_id] = {
            "mag_id": mag_id,
            "length": len(sequence),
            "md5": md5_sequence(sequence),
            "found": False,
        }
        selection_rows.append((mag_id, 1, contig_id, "Y"))

    if not candidates:
        raise SystemExit("ERROR: No single-contig MAGs were found")

    with open(selection_path, "w") as out:
        out.write("MAG_ID\tNumber_of_contigs\tContig_ID\tSINGLE_CONTIG_MAG\n")
        for row in selection_rows:
            out.write("\t".join(map(str, row)) + "\n")

    manifest_rows = []

    seen_assembly_ids = set()

    with open(output_reference, "w") as reference_out:
        for header, sequence in fasta_records(assembly_path):
            contig_id = fasta_id(header)
            if contig_id in seen_assembly_ids:
                raise SystemExit(f"ERROR: Duplicate contig ID in the assembler assembly: {contig_id}")
            seen_assembly_ids.add(contig_id)

            if contig_id not in candidates:
                write_fasta_record(reference_out, header, sequence)
                continue

            candidate = candidates[contig_id]
            if len(sequence) != candidate["length"] or md5_sequence(sequence) != candidate["md5"]:
                raise SystemExit(
                    "ERROR: The single-contig MAG sequence does not exactly match the assembler assembly "
                    f"for {candidate['mag_id']} / {contig_id}"
                )
            if len(sequence) < 2 * hifi_anchor + 2:
                raise SystemExit(
                    f"ERROR: {contig_id} is too short for a {hifi_anchor}-bp HiFi anchor on both sides"
                )

            original_split_after = len(sequence) // 2
            original_first_half = sequence[:original_split_after]
            original_second_half = sequence[original_split_after:]
            rotated_sequence = original_second_half + original_first_half
            end_start_left_base = len(original_second_half)
            end_start_right_base = end_start_left_base + 1

            if len(rotated_sequence) != len(sequence):
                raise SystemExit(f"ERROR: Rotation changed sequence length for {contig_id}")
            restored = rotated_sequence[len(original_second_half):] + rotated_sequence[:len(original_second_half)]
            if restored != sequence:
                raise SystemExit(f"ERROR: Rotation was not exactly reversible for {contig_id}")

            exact_terminal_overlap = longest_exact_terminal_overlap(sequence)

            write_fasta_record(reference_out, header, rotated_sequence)
            candidate["found"] = True
            candidate["junction0"] = end_start_left_base
            candidate["rotated_md5"] = md5_sequence(rotated_sequence)
            manifest_rows.append(
                (
                    candidate["mag_id"],
                    contig_id,
                    len(sequence),
                    original_split_after,
                    original_split_after + 1,
                    end_start_left_base,
                    end_start_right_base,
                    exact_terminal_overlap,
                    candidate["md5"],
                    candidate["rotated_md5"],
                    "Y",
                )
            )

    missing = [contig_id for contig_id, data in candidates.items() if not data["found"]]

    if missing:
        shown = ", ".join(missing[:10])
        raise SystemExit(f"ERROR: {len(missing)} single-contig MAG contigs were absent from the assembler assembly: {shown}")

    manifest_rows.sort(key=lambda row: (row[0], row[1]))

    with open(manifest_path, "w") as out:
        out.write(
            "MAG_ID\tContig_ID\tContig_length_bp\t"
            "New_reference_END_original_position\tNew_reference_START_original_position\t"
            "Rotated_END_START_left_base\tRotated_END_START_right_base\t"
            "Longest_exact_terminal_overlap_bp_up_to_100000\t"
            "Original_sequence_MD5\tRotated_sequence_MD5\tROTATION_VALIDATED\n"
        )
        for row in manifest_rows:
            out.write("\t".join(map(str, row)) + "\n")

    with open(candidate_bed_path, "w") as bed_out, open(junction_bed_path, "w") as junction_out:
        for row in manifest_rows:
            _, contig_id, length, _, _, left_base, _, *_ = row
            junction0 = int(left_base)
            bed_out.write(f"{contig_id}\t0\t{length}\n")
            flank = hifi_anchor + 1
            region_start = max(0, junction0 - flank)
            region_end = min(int(length), junction0 + flank)
            junction_out.write(f"{contig_id}\t{region_start}\t{region_end}\n")

    print(f"Single-contig MAGs: {len(manifest_rows)}")

    print(f"Rotated reference: {output_reference}")


COMMANDS = {
    'match': cmd_match,
    'rotate': cmd_rotate,
}

if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] not in COMMANDS:
        sys.exit(__doc__)
    result = COMMANDS[sys.argv[1]](sys.argv[2:])
    sys.exit(result if isinstance(result, int) else 0)
