import csv
import os
import sys

plan_path, groups_path, grouped_root, validation_path, output_path, group_outputs_path = sys.argv[1:]

def read_tsv(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)

_, plan_rows = read_tsv(plan_path)
expected = {
    row["MAG_ID"]: {
        "group": row["Group_ID"], "assembler": row["Assembler"],
        "assembly": row["Assembly_FASTA"], "contig": row["Contig_ID"],
    }
    for row in plan_rows
}
_, group_rows = read_tsv(groups_path)
merged_rows = []
output_fields = None
seen = set()
group_output_rows = []
for group in group_rows:
    group_id = group["Group_ID"]
    final_path = os.path.join(grouped_root, group_id, "MAG_END_START_closure.final.tsv")
    if not os.path.isfile(final_path) or os.path.getsize(final_path) == 0:
        raise SystemExit(f"ERROR: Missing source-group result: {final_path}")
    fields, rows = read_tsv(final_path)
    if output_fields is None:
        output_fields = fields
    elif fields != output_fields:
        raise SystemExit(f"ERROR: Source-group final schema mismatch: {final_path}")
    group_output_rows.append((group_id, group["Assembler"], group["Assembly_FASTA"], group["MAG_count"], final_path))
    for row in rows:
        mag_id = row.get("MAG_ID", "")
        contig_id = row.get("Contig_ID", "")
        if mag_id not in expected:
            raise SystemExit(f"ERROR: Unexpected MAG_ID in {final_path}: {mag_id}")
        meta = expected[mag_id]
        if meta["group"] != group_id or meta["contig"] != contig_id:
            raise SystemExit(f"ERROR: Source provenance mismatch: {mag_id} / {contig_id}")
        if mag_id in seen:
            raise SystemExit(f"ERROR: MAG_ID occurs in multiple groups: {mag_id}")
        seen.add(mag_id)
        row["Source_assembler"] = meta["assembler"]
        row["Source_assembly_FASTA"] = meta["assembly"]
        merged_rows.append(row)

if output_fields is None:
    raise SystemExit("ERROR: no source-group circularity outputs")
for field in ("Source_assembler", "Source_assembly_FASTA"):
    if field not in output_fields:
        output_fields.append(field)

_, validation_rows = read_tsv(validation_path)
for item in validation_rows:
    mag_id = item["MAG_ID"]
    n = int(item["Number_of_contigs"])
    if n == 1:
        if mag_id not in seen:
            raise SystemExit(f"ERROR: Single-contig MAG missing from merged evidence: {mag_id}")
        continue
    row = {field: "NA" for field in output_fields}
    row.update({
        "MAG_ID": mag_id,
        "Number_of_contigs": str(n),
        "Contig_ID": "NA",
        "Contig_length_bp": "NA",
        "SINGLE_CONTIG": "N",
        "CIRCULARITY_ANALYZABLE": "N",
        "HIFI_PRIMARY_SUPPLEMENTARY_ANALYZABLE_Y_N": "N",
        "HIFI_ALL_MAPPED_ANALYZABLE_Y_N": "N",
        "HIFI_END_START_SUPPORTED": "NA",
        "CIRCULARITY_SUPPORTED": "NA",
        "Source_assembler": "NA",
        "Source_assembly_FASTA": "NA",
    })
    merged_rows.append(row)

merged_rows.sort(key=lambda r: r["MAG_ID"])
with open(output_path, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=output_fields, delimiter="\t", lineterminator="\n")
    writer.writeheader(); writer.writerows(merged_rows)

with open(group_outputs_path, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(("Group_ID", "Assembler", "Assembly_FASTA", "MAG_count", "HiFi_final_result"))
    writer.writerows(group_output_rows)
