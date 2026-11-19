# Database preparation

Activate the `MAGverify` Conda environment and run preparation commands from the repository root. Chimera analysis uses an MMseqs2 GTDB protein taxonomy database (`--database`) or an existing compatible gene taxonomy table (`--taxonomy`). MAG-level GTDB classifications (`--mag-taxonomy`) are a separate input used for domain-matched SCG annotation. Prepare these inputs once; the main analysis modules do not download GTDB data.

## 1. GTDB R232 MMseqs2 protein taxonomy database

Skip this step if an existing compatible gene taxonomy table will be supplied with `--taxonomy`.

Download the three files from the fixed **GTDB R232 (232.0)** release:

| Download | Content |
| --- | --- |
| [`gtdb_proteins_aa_reps_r232.tar.gz`](https://data.gtdb.ecogenomic.org/releases/release232/232.0/genomic_files_reps/gtdb_proteins_aa_reps_r232.tar.gz) | Predicted amino-acid sequences from GTDB representative genomes (~123 GB compressed) |
| [`bac120_taxonomy_r232.tsv`](https://data.gtdb.ecogenomic.org/releases/release232/232.0/bac120_taxonomy_r232.tsv) | Bacterial genome accession → full GTDB lineage |
| [`ar53_taxonomy_r232.tsv`](https://data.gtdb.ecogenomic.org/releases/release232/232.0/ar53_taxonomy_r232.tsv) | Archaeal genome accession → full GTDB lineage |

**Download with aria2c** (the commands only download the files):

```bash
GTDB_DIR=/path/to/GTDB_R232_files
BASE=https://data.gtdb.ecogenomic.org/releases/release232/232.0
mkdir -p "$GTDB_DIR"

aria2c -c -x 4 -s 4 -d "$GTDB_DIR" "$BASE/genomic_files_reps/gtdb_proteins_aa_reps_r232.tar.gz"
aria2c -c -d "$GTDB_DIR" "$BASE/bac120_taxonomy_r232.tsv"
aria2c -c -d "$GTDB_DIR" "$BASE/ar53_taxonomy_r232.tsv"
```

**Build the MMseqs2 taxonomy database** from the downloaded files:

```bash
DB=/path/to/databases/GTDB_R232
mkdir -p "$(dirname "$DB")"

MAGverify prep-chimera \
  --gtdb-dir "$GTDB_DIR" \
  --database-prefix "$DB" \
  -o /path/to/MAGverify_output \
  -t 100

ls -lh "${DB}.dbtype" "${DB}_mapping" "${DB}_taxonomy"
```

`--gtdb-dir` points to the three downloaded files. `--database-prefix` specifies the full MMseqs2 database prefix, **not a directory**; omit it to use `/path/to/MAGverify_output/0_prep/GTDB_R232`. Temporary MMseqs2 files are written to `OUTPUT_DIR/tmp_prep_chimera/`, which needs substantially more disk space than the compressed archive alone and can be deleted after the database is built. Supply the chosen prefix to `MAGverify chimera --database`.

## 2. MAG taxonomy TSV

For domain-matched SCG annotation, supply a tab-separated MAG taxonomy table with these columns:

```text
ID	GTDB_Domain	GTDB_Phylum	GTDB_Class	GTDB_Order
ExampleMAG	d__Bacteria	p__Nitrospirota	c__Nitrospiria	o__Nitrospirales
```

The example classification illustrates the format only. Each `ID` must match a MAG FASTA filename stem exactly (for example, `ExampleMAG.fa` → `ExampleMAG`), with one row per MAG and no duplicates. GTDB rank prefixes are `d__`, `p__`, `c__` and `o__`; unassigned ranks can use an empty prefix (for example, `o__`).

If a compatible table already exists, use it directly as `--mag-taxonomy FILE`. If existing **GTDB-Tk** output is available, convert its `gtdbtk.bac120.summary.tsv` and `gtdbtk.ar53.summary.tsv` files without rerunning classification. Run the converter from the MAGverify repository root:

```bash
PREP=/path/to/MAGverify_output/0_prep
mkdir -p "$PREP"
python -m scripts.prep_gtdbtk_summary_to_mag_taxonomy \
  --gtdbtk-dir /path/to/existing_GTDBTk_results \
  -m /path/to/MAGs \
  --output "$PREP/MAG_taxonomy.tsv"
```

If no GTDB-Tk result is available, classify the MAGs first with `gtdbtk classify_wf`, then convert its summary files as above. GTDB-Tk is not installed with MAGverify; see [GTDB-Tk](https://github.com/ecogenomics/gtdbtk) for installation, reference data and options.

## Using the prepared database

Supply the MMseqs2 DB prefix to `MAGverify chimera --database` and the MAG-level taxonomy table to `--mag-taxonomy`, when available. If a compatible gene taxonomy TSV already exists, use `--taxonomy FILE` in place of `--database PREFIX`. The user-prepared `0_prep/` is separate from MAGverify's internal `00_prep/`.
