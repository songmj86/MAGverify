# MAGverify analysis protocol

## Installation (once)

Create one environment with MAGverify, anvi'o and its required bioinformatics programs. Run from the repository root:

```bash
conda env create -f environment.yml
conda activate MAGverify
MAGverify doctor
```

Use the active `MAGverify` environment for **every** module below. Do not pass `--anvio-env`. See [docs/INSTALL.md](docs/INSTALL.md) for setup and verification.

## Database preparation (once, before analysis)

Prepare external databases and MAG taxonomy under `/path/to/MAGverify_output/0_prep/`. This **user-prepared `0_prep/`** is not the internal `00_prep/` created by MAGverify. All ten analysis commands are unchanged; the Chimera module still uses `--database PREFIX` and `--mag-taxonomy FILE`.

### A. GTDB R232 protein taxonomy database for MMseqs2

The MMseqs2 input is **the protein FASTA archive for GTDB R232 species representatives** plus **two genome-to-GTDB taxonomy tables**. Download the **fixed R232 release (232.0)**, not an unpinned `latest` release. Official files:

| Download | Content | Purpose |
| --- | --- | --- |
| [`gtdb_proteins_aa_reps_r232.tar.gz`](https://data.gtdb.ecogenomic.org/releases/release232/232.0/genomic_files_reps/gtdb_proteins_aa_reps_r232.tar.gz) | Predicted amino-acid sequences from GTDB representative genomes (~123 GB compressed) | Input protein sequences to `mmseqs tar2db` |
| [`bac120_taxonomy_r232.tsv`](https://data.gtdb.ecogenomic.org/releases/release232/232.0/bac120_taxonomy_r232.tsv) | Bacterial genome accession → full GTDB lineage | Bacterial taxon mapping |
| [`ar53_taxonomy_r232.tsv`](https://data.gtdb.ecogenomic.org/releases/release232/232.0/ar53_taxonomy_r232.tsv) | Archaeal genome accession → full GTDB lineage | Archaeal taxon mapping |

Official [R232 protein directory](https://data.gtdb.ecogenomic.org/releases/release232/232.0/genomic_files_reps/) and [R232 taxonomy directory](https://data.gtdb.ecogenomic.org/releases/release232/232.0/). **Do not use** `bac120_msa_reps_r232.faa.gz` or `ar53_msa_reps_r232.faa.gz`; these are concatenated marker alignments rather than the complete representative predicted protein sets. `gtdb_proteins_nt_reps_r232.tar.gz` is nucleotides, not proteins. The protein archive is large, and decompression, temporary databases and taxonomy indexing require substantially more disk space than the compressed download alone.

**Run with the single MAGverify Conda environment active.** The `prep-chimera` command invokes `modules/00_prep_chimera.sh` internally. It downloads exactly these three files, builds the protein sequence DB and assigns GTDB lineage IDs to sequences using the MMseqs2 GTDB taxonomy workflow:

```bash
conda activate MAGverify
PREP=/path/to/MAGverify_output/0_prep
MAGverify prep-chimera \
  --output-dir /path/to/MAGverify_output \
  --threads 100
```

The `prep-chimera` command performs the following workflow steps (view it for full, executable code):

1. Downloads R232 proteins to `$PREP/gtdb_r232_downloads/gtdb_proteins_aa_reps_r232.tar.gz` and bacterial/archaeal taxonomy TSVs to the same folder.
2. Reads nested `*_protein.faa.gz` using `mmseqs tar2db` and creates the final MMseqs2 sequence database with `mmseqs createdb`.
3. Uses **both** GTDB rank tables to create a GTDB taxonomic tree and accession-to-taxon mapping (`nodes.dmp`, `names.dmp`, `mapping_genomes`).
4. Runs `mmseqs createtaxdb` with `--ncbi-tax-dump` (the locally generated GTDB-derived taxdump), `--tax-mapping-file` and `--tax-mapping-mode 1`. Do not use a downloaded NCBI taxonomy mapping for GTDB protein IDs.

Verify the completed database:

```bash
ls -lh "$PREP/GTDB_R232.dbtype" "$PREP/GTDB_R232_mapping" "$PREP/GTDB_R232_taxonomy"
cat "$PREP/GTDB_R232.version"
```

The command automatically creates the `0_prep/` directory inside the specified `--output-dir` (do not append `/0_prep` to that option). The resulting **MMseqs2 database prefix** is `/path/to/MAGverify_output/0_prep/GTDB_R232`, which is the exact `--database` argument passed to `MAGverify chimera`. The files in `gtdb_r232_downloads/` are input downloads, not MAGverify analysis results. The `mmseqs databases GTDB ...` convenience command fetches the upstream **latest** release and is not recommended when R232 must remain fixed. A compatible, previously created gene taxonomy table may be provided using `--taxonomy FILE` instead of `--database`.

### B. MAG taxonomy TSV (`--mag-taxonomy`)

`MAGverify chimera` expects a prepared tab-separated file and **does not read GTDB-Tk summary files directly**. Its supported header includes:

```text
ID	GTDB_Domain	GTDB_Phylum	GTDB_Class	GTDB_Order
Bin13L	d__Bacteria	p__Nitrospirota	c__Nitrospiria	o__Nitrospirales
```

The sample row shows **format only**, not a verified classification of Bin13L. Each `ID` must exactly match a MAG FASTA filename stem (e.g. `Bin13L.fa` → `Bin13L`). GTDB ranks use the prefixes `d__`, `p__`, `c__` and `o__`; an unassigned taxon rank may be an empty prefix such as `o__`. Supply one row per input MAG, with no duplicate IDs. If an existing MAG taxonomy table already follows this format, copy it to `$PREP/MAG_taxonomy.tsv` and **skip GTDB-Tk and the conversion step**.

For **existing GTDB-Tk output**, keep the original `gtdbtk.bac120.summary.tsv` and `gtdbtk.ar53.summary.tsv` files unchanged and use the provided **preparation-only** converter. It reads the original `user_genome` and `classification` columns and generates the TSV expected by the existing main-analysis code:

```bash
PREP=/path/to/MAGverify_output/0_prep
python -m scripts.prep_gtdbtk_summary_to_mag_taxonomy \
  --gtdbtk-dir /path/to/existing_GTDBTk_results \
  --mag-dir /path/to/MAGs \
  --output "$PREP/MAG_taxonomy.tsv"
head -n 3 "$PREP/MAG_taxonomy.tsv"
```

The converter validates duplicate IDs, GTDB rank strings and MAG FASTA membership. You do **not** need to rerun GTDB-Tk when an existing `classify_wf` output covers the MAGs at the intended GTDB release.

### C. GTDB-Tk R232 classification — only if MAG taxonomy is unavailable

The optional `gtdbtk classify_wf` step requires a **separate GTDB-Tk reference-data package**. It is **not** the MMseqs2 protein DB above. GTDB-Tk is installed in the same `MAGverify` Conda environment, but external reference data must be downloaded. Download the official pinned R232 package (~56.6 GB compressed, ~100 GB extracted):

[GTDB-Tk R232 full package](https://data.gtdb.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz) · [Package directory](https://data.gtdb.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/)

```bash
conda activate MAGverify
PREP=/path/to/MAGverify_output/0_prep
mkdir -p "$PREP/GTDBTk_R232"
wget -c -O "$PREP/GTDBTk_R232/gtdbtk_r232_data.tar.gz" \
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

Adapt `--extension fa` if MAG FASTAs are `.fna` or `.fasta`. GTDB-Tk R232 classification may require at least ~140 GB RAM for bacterial placements. `GTDBTK_DATA_PATH` is only needed to run **GTDB-Tk**; it is not an argument to MAGverify modules. An existing installation with the R232 GTDB-Tk reference package can be reused rather than downloaded again. These large data downloads have not been executed as part of repository packaging.

### Prepared inputs and their use

```text
/path/to/MAGverify_output/
├── 0_prep/                                # user-created (before analysis)
│   ├── gtdb_r232_downloads/               # R232 protein archive, bac120/ar53 TSVs
│   ├── GTDB_R232                          # MMseqs2 protein database prefix
│   ├── GTDB_R232.dbtype
│   ├── GTDB_R232_mapping
│   ├── GTDB_R232_taxonomy
│   ├── MAG_taxonomy.tsv                   # input to --mag-taxonomy
│   └── GTDBTk_results/                   # optional (only if classification run)
└── 00_prep/                               # automatically generated by MAGverify
```

The relevant analysis command is unchanged:

```bash
MAGverify chimera \
  --assemblies-dir /path/to/assemblies/ \
  --mag-dir /path/to/MAGs \
  --database /path/to/MAGverify_output/0_prep/GTDB_R232 \
  --mag-taxonomy /path/to/MAGverify_output/0_prep/MAG_taxonomy.tsv \
  --output-dir /path/to/MAGverify_output \
  --threads 100
```

## Usage

Run each analysis separately and reuse the **same `--output-dir`** across all modules. The `--assemblies-dir` input must be a **directory** of original, unbinned assembly FASTA files (`.fa`, `.fasta` or `.fna`), not a single combined `assemblies.fa` file or a MAG directory. Each original assembler FASTA is mapped separately. If the FASTAs are in different locations, use repeated `--source-assembly /path/to/file.fa` options instead of `--assemblies-dir`. MAGverify automatically creates `00_prep/combined_assembly.fasta` when a downstream analysis needs it. Examples request 100 threads and assume at least 100 CPU threads are allocated.

Expected input example:

```text
/path/to/assemblies/
├── hifiasm_meta.fa
├── metaFlye.fa
└── metaMDBG.fa
```

## 1. Read mapping

```bash
MAGverify map \
  --assemblies-dir /path/to/assemblies/ \
  --reads /path/to/hifi.fastq.gz \
  --output-dir /path/to/MAGverify_output \
  --threads 100
```

Maps HiFi reads **separately** against each original assembler FASTA. Step 1 creates `/path/to/MAGverify_output/01_mapping/assembler_mapping_outputs.tsv` and exposes `/path/to/MAGverify_output/results/mapping_manifest.tsv` as a symbolic link (`/path/to/MAGverify_output` is the chosen `--output-dir`). It records each assembler FASTA and the generated BAM paths. Subsequent modules discover the existing manifest automatically from the same output directory when `--mapping-manifest` is omitted. The manifest is reused automatically even when downstream modules also receive `--reads`. You do not need to supply `--mapping-manifest` in the usual workflow. If no mapping exists, a module receiving both `--reads` and original assemblies can initiate mapping. To update existing mapping after changing the reads or assemblies, rerun `MAGverify map` first.

## 2. Circularity

```bash
MAGverify circularity \
  --assemblies-dir /path/to/assemblies/ \
  --mag-dir /path/to/MAGs \
  --reads /path/to/hifi.fastq.gz \
  --output-dir /path/to/MAGverify_output \
  --threads 100
```

Evaluates read support for circular junctions. Main output: `circularity_mag.tsv`.

## 3. Clipping

```bash
MAGverify clipping \
  --assemblies-dir /path/to/assemblies/ \
  --mag-dir /path/to/MAGs \
  --output-dir /path/to/MAGverify_output \
  --threads 100
```

Detects clipping candidates and evaluates HiFi read support. The main outputs are `clipping_sites.tsv` and `clipping_mag.tsv`. Step 1 mapping is reused automatically when the same `--output-dir` is specified. If no mapping manifest exists, provide `--reads /path/to/hifi.fastq.gz` instead.

## 4. False closure

```bash
MAGverify false-closure \
  --assemblies-dir /path/to/assemblies/ \
  --mag-dir /path/to/MAGs \
  --assembler-circularity /path/to/assembler_circularity.tsv \
  --output-dir /path/to/MAGverify_output \
  --threads 100
```

This assessment uses the assembler-reported circularity calls, validated clipping and HiFi mapping. SCG annotations are generated internally if missing; they need not be run as a separate module. Main outputs: `false_closure_contig.tsv`, `false_closure_mag.tsv`. The assembler circularity table is distinct from the MAGverify read-supported circularity result.

## 5. Phantom sequences

```bash
MAGverify phantom \
  --assemblies-dir /path/to/assemblies/ \
  --mag-dir /path/to/MAGs \
  --reads /path/to/hifi.fastq.gz \
  --output-dir /path/to/MAGverify_output \
  --threads 100
```

Assesses zero-coverage regions and raw-read k-mer evidence. Main outputs: `phantom_contig.tsv`, `phantom_mag.tsv`.

## 6. Repeats

```bash
MAGverify repeats \
  --assemblies-dir /path/to/assemblies/ \
  --mag-dir /path/to/MAGs \
  --output-dir /path/to/MAGverify_output \
  --threads 100
```

Quantifies repeat regions. Main outputs: `repeats_contig.tsv`, `repeats_mag.tsv`. This is a separate analysis from false duplication: a repeat is not itself evidence of false duplication.

## 7. Taxonomic chimera (gene taxonomy and SCG included)

```bash
MAGverify chimera \
  --assemblies-dir /path/to/assemblies/ \
  --mag-dir /path/to/MAGs \
  --database /path/to/MAGverify_output/0_prep/GTDB_R232 \
  --mag-taxonomy /path/to/MAGverify_output/0_prep/MAG_taxonomy.tsv \
  --output-dir /path/to/MAGverify_output \
  --threads 100
```

The Chimera module generates the all-contig gene taxonomy table via MMseqs2 only when the table is missing. `--database` is a taxonomy database **prefix** with a `.dbtype` file, not simply a directory. An existing gene taxonomy table can be supplied with `--taxonomy /path/to/gene_taxonomy.tsv` in place of `--database`. SCG annotation is generated internally, when needed and when anvi'o HMM programs are available. `--mag-taxonomy` is used to assign domain-matched SCG annotations; without it, that supplementary annotation is unavailable. Chimera uses validated clipping results from Step 3. Main outputs: `chimera_contig.tsv` and `chimera_mag.tsv` (plus gene taxonomy and SCG intermediates).

## 8. False duplication

```bash
MAGverify false-duplication \
  --assemblies-dir /path/to/assemblies/ \
  --mag-dir /path/to/MAGs \
  --output-dir /path/to/MAGverify_output \
  --threads 100
```

Evaluates candidate duplication artifacts using mapping and validated clipping. Main outputs: `false_duplication_contig.tsv`, `false_duplication_mag.tsv`.

## 9. Assembly variants

```bash
MAGverify variants \
  --assemblies-dir /path/to/assemblies/ \
  --mag-dir /path/to/MAGs \
  --reads /path/to/hifi.fastq.gz \
  --output-dir /path/to/MAGverify_output \
  --threads 100
```

Analyze assembly variation and read-level allele support. A machine with sufficient available memory is necessary; use `--memory-gb` when appropriate. Main outputs: `variants_contig.tsv`, `variants_mag.tsv`.

## 10. Integrated QC summary

```bash
MAGverify summary \
  --assemblies-dir /path/to/assemblies/ \
  --mag-dir /path/to/MAGs \
  --output-dir /path/to/MAGverify_output \
  --threads 100
```

Integrates completed assessments into `summary_mag.tsv` and `summary_contig.tsv`. The module then validates MAG IDs, contig memberships, contig counts and the consistency of the other summary tables. Validation failure terminates the command with a nonzero exit status. Previously completed analysis modules are not re-executed during summary validation. For partial runs, pass `--include clipping,chimera` (or the other completed modules) so the summary does not expect absent analysis results.


## Input and mapping notes

- `--assemblies-dir DIR` accepts a folder of original assembler FASTAs. For a single original assembly use `--source-assembly /path/to/assembly.fa`, not `--assemblies-dir /path/to/assembly.fa`.
- Original FASTA files should be separated by assembler; do not merge competing assemblers into one mapping reference.
- The mapping manifest is **created in Step 1** and reused across later steps. It is not an input file to prepare manually.
- Modules automatically reuse the manifest from the same `--output-dir`, including those that need `--reads`. The optional `--mapping-manifest FILE` flag is only for using a manifest stored elsewhere. Rerun Step 1 if the reads or assemblies change.
- `--output-dir` is the documented output option. The older `--output-directory` spelling remains accepted as a backwards-compatible alias.

## Data prerequisites

- Original assembly FASTAs must use unique contig IDs across assembler files.
- MAG FASTAs must use the same contig identifiers as their assembly FASTAs.
- Prepare the GTDB R232 MMseqs2 gene taxonomy database in `0_prep/` (or supply a compatible `--taxonomy` TSV).
- `--mag-taxonomy` (MAG-level GTDB classification) and `--database` (gene-level MMseqs2 database) are different inputs.
- Install anvi'o in the same `MAGverify` Conda environment. Activate it once before running commands; see [INSTALL.md](docs/INSTALL.md).
- There is no separate `00_config.sh`, `taxonomy` or `scg` command in the user-facing protocol.
