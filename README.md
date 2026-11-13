# MAGverify

MAGverify assesses the integrity of metagenome assemblies and metagenome-assembled genomes (MAGs) using PacBio HiFi reads. Its ten modules cover read mapping, circularity, clipping, false closure, phantom sequences, repeats, taxonomic chimeras, false duplication, assembly variants and integrated QC.

## 1. Installation

MAGverify supports Linux with Conda or Mamba. Create the environment from the repository root:

```bash
git clone https://github.com/songmj86/MAGverify.git
cd MAGverify
mamba env create -f environment.yml
conda activate MAGverify
MAGverify doctor
MAGverify list
```

The environment includes anvi'o v9 and the command-line tools used by MAGverify. `MAGverify doctor` checks executable availability; it is not a full pipeline test. See [Installation](docs/INSTALL.md) for troubleshooting.

## 2. Database preparation

Chimera analysis uses an MMseqs2 GTDB protein taxonomy database (`--database`) or an existing compatible gene taxonomy table (`--taxonomy`). MAG-level GTDB classifications (`--mag-taxonomy`) are a separate input used for domain-matched SCG annotation. Prepare these inputs once; the main analysis modules do not download GTDB data.

### 3.1. GTDB R232 MMseqs2 protein taxonomy database

Skip this step if an existing compatible gene taxonomy table will be supplied with `--taxonomy`.

Download the three files from the fixed **GTDB R232 (232.0)** release:

| Download | Content |
| --- | --- |
| [`gtdb_proteins_aa_reps_r232.tar.gz`](https://data.gtdb.ecogenomic.org/releases/release232/232.0/genomic_files_reps/gtdb_proteins_aa_reps_r232.tar.gz) | Predicted amino-acid sequences from GTDB representative genomes (~123 GB compressed) |
| [`bac120_taxonomy_r232.tsv`](https://data.gtdb.ecogenomic.org/releases/release232/232.0/bac120_taxonomy_r232.tsv) | Bacterial genome accession → full GTDB lineage |
| [`ar53_taxonomy_r232.tsv`](https://data.gtdb.ecogenomic.org/releases/release232/232.0/ar53_taxonomy_r232.tsv) | Archaeal genome accession → full GTDB lineage |

**Download with aria2c** (the commands only download the files):

```bash
GTDB_DIR=/path/to/GTDB_R232
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
  --output-dir /path/to/MAGverify_output \
  --threads 100

ls -lh "${DB}.dbtype" "${DB}_mapping" "${DB}_taxonomy"
```

`--gtdb-dir` points to the three downloaded files. `--database-prefix` specifies the full MMseqs2 database prefix, **not a directory**; omit it to use `/path/to/MAGverify_output/0_prep/GTDB_R232`. The build uses `OUTPUT_DIR/0_prep/` for intermediate files and needs substantially more disk space than the compressed archive alone. Supply the chosen prefix to `MAGverify chimera --database`.

### 3.2. MAG taxonomy TSV

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
  --mag-dir /path/to/MAGs \
  --output "$PREP/MAG_taxonomy.tsv"
```

### 3.3. GTDB-Tk classification (only if MAG taxonomy is unavailable)

GTDB-Tk uses a **separate** R232 reference package, not the MMseqs2 protein database. Skip this section if a compatible MAG taxonomy table or existing GTDB-Tk classification is available. The [GTDB-Tk R232 package](https://data.gtdb.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz) is approximately 56.6 GB compressed and 100 GB extracted.

```bash
PREP=/path/to/MAGverify_output/0_prep
mkdir -p "$PREP/GTDBTk_R232"
aria2c -c -x 4 -s 4 -d "$PREP/GTDBTk_R232" \
  https://data.gtdb.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz

tar -xzf "$PREP/GTDBTk_R232/gtdbtk_r232_data.tar.gz" -C "$PREP/GTDBTk_R232"
META_DIR=$(find "$PREP/GTDBTk_R232" -maxdepth 3 -type d -name metadata -print -quit)
test -n "$META_DIR" || { echo 'GTDB-Tk metadata directory not found' >&2; exit 1; }
export GTDBTK_DATA_PATH="$(dirname "$META_DIR")"
gtdbtk check_install

gtdbtk classify_wf \
  --genome_dir /path/to/MAGs \
  --out_dir "$PREP/GTDBTk_results" \
  --extension fa \
  --cpus 100

python -m scripts.prep_gtdbtk_summary_to_mag_taxonomy \
  --gtdbtk-dir "$PREP/GTDBTk_results" \
  --mag-dir /path/to/MAGs \
  --output "$PREP/MAG_taxonomy.tsv"
```

Change `--extension fa` for `.fna` or `.fasta` inputs. GTDB-Tk bacterial classification may need approximately 140 GB RAM. Set `GTDBTK_DATA_PATH` in any new shell session before using GTDB-Tk. The user-prepared `0_prep/` is distinct from MAGverify's internally generated `00_prep/`.


## 3. Input data

- **Assemblies** — Original, unbinned assembly FASTA files (`.fa`, `.fasta` or `.fna`). Supply a directory with `--assemblies-dir`, or use repeated `--source-assembly /path/to/file.fa` options when files are in different locations.
- **MAGs** — One FASTA per MAG in `--mag-dir`. MAG contig IDs must match those in the original assembly FASTAs, and assembly contig IDs must be unique across input files.
- **Reads** — PacBio HiFi reads (`.fastq` or `.fastq.gz`) supplied with `--reads` when needed.
- **Output** — Use the same `--output-dir` for every module so that completed results can be reused.

Example assembly directory:

```text
/path/to/assemblies/
├── hifiasm_meta.fa
├── metaFlye.fa
└── metaMDBG.fa
```

**Keep assembler FASTAs separate.** MAGverify maps reads to each original assembly separately to avoid competing alignments between similar contigs from different assemblers. Do not merge these files into one mapping reference.


## 4. Run the analysis modules

Run each command separately in the following order, using the **same `--output-dir`** throughout. The examples use 100 threads; adjust this to the CPU allocation. The mapping manifest created in Step 1 is found and reused automatically by later modules. No `00_config.sh`, separate gene-taxonomy module or separate SCG module is needed.

### 1. Read mapping

```bash
MAGverify map \
  --assemblies-dir /path/to/assemblies \
  --reads /path/to/hifi.fastq.gz \
  --output-dir /path/to/MAGverify_output \
  --threads 100
```

Reads are mapped separately to each assembler FASTA. The module writes `01_mapping/assembler_mapping_outputs.tsv` and links its manifest under `results/mapping_manifest.tsv`. Rerun mapping if the input reads or assemblies change.

### 2. Circularity

```bash
MAGverify circularity \
  --assemblies-dir /path/to/assemblies \
  --mag-dir /path/to/MAGs \
  --reads /path/to/hifi.fastq.gz \
  --output-dir /path/to/MAGverify_output \
  --threads 100
```

### 3. Clipping

```bash
MAGverify clipping \
  --assemblies-dir /path/to/assemblies \
  --mag-dir /path/to/MAGs \
  --output-dir /path/to/MAGverify_output \
  --threads 100
```

This and later modules reuse Step 1 mapping from the same output directory. If no mapping exists, supply `--reads` alongside the assemblies so mapping can be prepared.

### 4. False closure

```bash
MAGverify false-closure \
  --assemblies-dir /path/to/assemblies \
  --mag-dir /path/to/MAGs \
  --assembler-circularity /path/to/assembler_circularity.tsv \
  --output-dir /path/to/MAGverify_output \
  --threads 100
```

`--assembler-circularity` takes the assembler-reported circularity table, not the read-supported circularity result from Step 2.

### 5. Phantom sequences

```bash
MAGverify phantom \
  --assemblies-dir /path/to/assemblies \
  --mag-dir /path/to/MAGs \
  --reads /path/to/hifi.fastq.gz \
  --output-dir /path/to/MAGverify_output \
  --threads 100
```

### 6. Repeats

```bash
MAGverify repeats \
  --assemblies-dir /path/to/assemblies \
  --mag-dir /path/to/MAGs \
  --output-dir /path/to/MAGverify_output \
  --threads 100
```

Repeat detection and false-duplication testing are separate assessments.

### 7. Taxonomic chimera

```bash
MAGverify chimera \
  --assemblies-dir /path/to/assemblies \
  --mag-dir /path/to/MAGs \
  --database /path/to/databases/GTDB_R232 \
  --mag-taxonomy /path/to/MAGverify_output/0_prep/MAG_taxonomy.tsv \
  --output-dir /path/to/MAGverify_output \
  --threads 100
```

Use the MMseqs2 database prefix from Section 3.1, or replace `--database PREFIX` with `--taxonomy /path/to/gene_taxonomy.tsv` when a compatible gene taxonomy table already exists. `--mag-taxonomy` is a separate MAG-level classification input. Gene taxonomy and SCG processing occur within this module as needed; chimera assessment also uses validated clipping results.

### 8. False duplication

```bash
MAGverify false-duplication \
  --assemblies-dir /path/to/assemblies \
  --mag-dir /path/to/MAGs \
  --output-dir /path/to/MAGverify_output \
  --threads 100
```

### 9. Assembly variants

```bash
MAGverify variants \
  --assemblies-dir /path/to/assemblies \
  --mag-dir /path/to/MAGs \
  --reads /path/to/hifi.fastq.gz \
  --output-dir /path/to/MAGverify_output \
  --threads 100
```

This step may need substantial RAM; set `--memory-gb` when appropriate.

### 10. Integrated QC summary

```bash
MAGverify summary \
  --assemblies-dir /path/to/assemblies \
  --mag-dir /path/to/MAGs \
  --output-dir /path/to/MAGverify_output \
  --threads 100
```

The summary combines completed assessments and validates MAG IDs, contig memberships, counts and table consistency. Validation failures return a nonzero exit status. For a partial run, use `--include clipping,chimera` (or the modules completed). The summary does not rerun analysis modules.

## 5. Results

Final tables are linked under `OUTPUT_DIR/results/`; module directories retain the intermediate files.

| Module | Main output files |
| --- | --- |
| Read mapping | `mapping_manifest.tsv` |
| Circularity | `circularity_mag.tsv` |
| Clipping | `clipping_sites.tsv`, `clipping_mag.tsv` |
| False closure | `false_closure_contig.tsv`, `false_closure_mag.tsv` |
| Phantom sequences | `phantom_contig.tsv`, `phantom_mag.tsv` |
| Repeats | `repeats_contig.tsv`, `repeats_mag.tsv` |
| Taxonomic chimera | `chimera_contig.tsv`, `chimera_mag.tsv` (plus gene taxonomy and SCG outputs when generated) |
| False duplication | `false_duplication_contig.tsv`, `false_duplication_mag.tsv` |
| Assembly variants | `variants_contig.tsv`, `variants_mag.tsv` |
| Integrated QC summary | `summary_contig.tsv`, `summary_mag.tsv` |

Keep the results directory and its module subdirectories together so the symbolic links remain valid.

## 6. Additional information

- [Installation and environment checks](docs/INSTALL.md)
- [Module inputs and outputs](docs/MODULES.md)
- [Input formats](docs/INPUT_FORMATS.md)
- [Recorded tool versions](docs/TOOL_VERSIONS.md)

MAGverify uses command-line options to specify paths and threads; no editable `00_config.sh` is needed. The repository's automated tests cover command routing and selected synthetic fixtures. A complete end-to-end run on real HiFi data is not yet verified.
