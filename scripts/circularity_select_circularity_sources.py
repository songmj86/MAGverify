import gzip
import hashlib
import os
import sys
from collections import defaultdict
from pathlib import Path

(
    assembly_dir,
    assembly_fasta_names,
    mag_dir,
    grouped_mag_root,
    plan_path,
    groups_path,
    validation_path,
    hifi_anchor_text,
) = sys.argv[1:]
hifi_anchor = int(hifi_anchor_text)

allowed_suffixes = (
    ".fa", ".fna", ".fasta", ".fas",
    ".fa.gz", ".fna.gz", ".fasta.gz", ".fas.gz",
)


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


def mag_id_from_path(path):
    name = Path(path).name
    lower = name.lower()
    for suffix in sorted(allowed_suffixes, key=len, reverse=True):
        if lower.endswith(suffix):
            return name[:-len(suffix)]
    return name


def md5_sequence(sequence):
    return hashlib.md5(sequence.encode("ascii")).hexdigest()


assembly_names = assembly_fasta_names.split()
if not assembly_names:
    raise SystemExit("ERROR: No source assembly FASTA names were configured")

assemblies = []
seen_labels = set()
for name in assembly_names:
    path = os.path.realpath(os.path.join(assembly_dir, name))
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise SystemExit(f"ERROR: Source assembly missing or empty: {path}")
    label = Path(name).name
    if label.lower().endswith(".gz"):
        label = label[:-3]
    label = Path(label).stem
    if label in seen_labels:
        raise SystemExit(f"ERROR: Duplicate source assembly label: {label}")
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
    raise SystemExit("ERROR: No eligible single-contig MAGs were found")

matches = defaultdict(list)
for label, assembly_path in assemblies:
    seen_ids = set()
    for header, sequence in fasta_records(assembly_path):
        contig_id = fasta_id(header)
        if contig_id in seen_ids:
            raise SystemExit(
                f"ERROR: Duplicate contig ID in source assembly {assembly_path}: {contig_id}"
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
            "ERROR: Eligible MAG contig was not found by exact ID and sequence in any "
            f"configured source assembly: {mag_id} / {item['contig']}"
        )
    if len(source_matches) > 1:
        sources = ", ".join(label for label, _ in source_matches)
        raise SystemExit(
            f"ERROR: MAG contig matches multiple source assemblies: {mag_id} / {sources}"
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
        "Contig_length_bp\tSequence_MD5\tSOURCE_ASSEMBLY_EXACT_MATCH\n"
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
        "SOURCE_ASSEMBLY_EXACT_MATCH\n"
    )
    for row in sorted(validation_rows, key=lambda value: value[0]):
        out.write("\t".join(map(str, row)) + "\n")

print(f"Eligible source-matched MAGs: {len(plan_rows)}")
print(f"Source groups: {len(group_counts)}")
