import csv
import sys
selection_path, manifest_path, hifi_summary_path, output_path = sys.argv[1:]

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
