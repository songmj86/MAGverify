# MAGverify

MAGverify assesses the integrity of metagenome assemblies and metagenome-assembled genomes (MAGs) using PacBio HiFi reads. Its modules cover read mapping, circularity, clipping, false closure, phantom sequences, repeats, taxonomic chimeras, false duplication, assembly variants and an integrated QC summary.

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

MAGverify doctor checks the installation path and availability of required programs.
MAGverify list displays the available commands and analysis modules with brief descriptions.

## 2. Preparation

Run the four preparation commands before analysis using the same -o (not included in end-to-end module).

| Command | Required | Output in `00_prep/` | Used by |
| --- | --- | --- | --- |
| [`prep-chimera`](#21-prep-chimera-gtdb-r232-database) | `GTDB_R232*` | `chimera` |
| [`prep-scg`](#22-prep-scg-single-copy-genes) | `scg.contig.tsv`, `genes.*` | `false-closure`, `chimera` |
| [`prep-circularity`](#23-prep-circularity-assembler-circularity) | `assembler_circularity.MAG.tsv`, `assembler_circularity.contig.tsv` | `false-closure`, `circularity` |
| [`prep-taxonomy`](#24-prep-taxonomy-mag-taxonomy) | `gtdbtk_taxonomy.MAG.tsv` | `chimera` |

Download commands and the format of each table are in [Input preparation](docs/INPUT_PREPARATION.md).

### 2.1 `prep-chimera`: GTDB R232 database

Builds the MMseqs2 protein taxonomy database that `chimera` uses to classify the predicted genes. Build it once; other projects can use it through `end-to-end --database`.

```bash
MAGverify prep-chimera -d /path/to/GTDB_R232_downloads -o /path/to/MAGverify_output -t 100
```

- `-d DIR`: Directory of the three GTDB R232 downloads below.
- Output: `00_prep/GTDB_R232*` (several files with this prefix). The temporary directory `00_prep/.cache/prep_chimera_work/` needs much more disk space than the archive and can be deleted afterwards.

| Download | Content |
| --- | --- |
| [`gtdb_proteins_aa_reps_r232.tar.gz`](https://data.gtdb.ecogenomic.org/releases/release232/232.0/genomic_files_reps/gtdb_proteins_aa_reps_r232.tar.gz) | Predicted amino-acid sequences from GTDB representative genomes (~123 GB compressed) |
| [`bac120_taxonomy_r232.tsv`](https://data.gtdb.ecogenomic.org/releases/release232/232.0/bac120_taxonomy_r232.tsv) | Bacterial genome accession → full GTDB lineage |
| [`ar53_taxonomy_r232.tsv`](https://data.gtdb.ecogenomic.org/releases/release232/232.0/ar53_taxonomy_r232.tsv) | Archaeal genome accession → full GTDB lineage |

### 2.2 `prep-scg`: single-copy genes

Predicts genes on all assemblies (Pyrodigal) and annotates single-copy genes with anvi'o (Bacteria_71, Archaea_76). `false-closure` uses the ribosomal proteins to screen small circular genomes, and `chimera` uses SCG redundancy.

```bash
MAGverify prep-scg -a /path/to/assemblies -o /path/to/MAGverify_output -t 100
```

- `-a DIR`: Input directory of contig files from each assembler (same `-a` as in the analysis).
- Output: `00_prep/scg.contig.tsv` (SCG hits per contig) and `00_prep/genes.faa`, `genes.gff`, `genes.coordinates.tsv` (reused by `chimera` and `variants`).

### 2.3 `prep-circularity`: assembler circularity

Collects the circularity that each assembler reported for its contigs. `false-closure` tests the contigs reported circular; `circularity` combines this call with the read evidence.

```bash
MAGverify prep-circularity \
  -c /path/to/hifiasm_meta/Read.p_ctg.gfa \
  -c /path/to/metaMDBG/contigs.fasta.gz \
  -c /path/to/metaflye/assembly_info.txt \
  -g /path/to/MAGs \
  -o /path/to/MAGverify_output
```

- `-c FILE`: one per assembler: hifiasm-meta `*.p_ctg.gfa` (or `*.p_ctg.noseq.gfa`), metaMDBG `contigs.fasta.gz`, metaFlye `assembly_info.txt`. The assembler is detected from the file content.
- `-g DIR`: Input directory of MAG files (one file per MAG; same `-g` as in the analysis).
- Output: `00_prep/assembler_circularity.MAG.tsv` and `assembler_circularity.contig.tsv`.

### 2.4 `prep-taxonomy`: MAG taxonomy

Converts GTDB-Tk results into a MAG taxonomy table. `chimera` uses the MAG domain to select the SCG collection and compares the taxa on both sides of a transition with the MAG taxon. Optional: without it, `chimera` runs without SCG redundancy.

```bash
MAGverify prep-taxonomy -d /path/to/GTDBTk_results -g /path/to/MAGs -o /path/to/MAGverify_output
```

- `-d DIR`: GTDB-Tk `classify_wf` output directory (`gtdbtk.bac120.summary.tsv`, `gtdbtk.ar53.summary.tsv`). Any GTDB-Tk release can be used.
- `-g DIR`: Input directory of MAG files (one file per MAG; same `-g` as in the analysis).
- Output: `00_prep/gtdbtk_taxonomy.MAG.tsv`.

## 3. Analysis (end-to-end)

After the preparation, one command runs all analysis modules in order. Inputs and output:

- **Assemblies** — An input directory of contig files (`.fa`, `.fasta`, or `.fna`), specified using `-a`. Contig IDs must be unique across assembly files.
- **MAGs** — An input directory of MAG files (one file per MAG), specified using `-g`. MAG contig IDs must match those in the assembly files.
- **Reads** — An input directory of PacBio HiFi read files (`.fastq`, `.fq`, `.fastq.gz`, or `.fq.gz`), specified using `-r` for read mapping. Several files are joined once into `00_prep/hifi_reads.fastq[.gz]`. The mapping manifest records the reads, so later modules reuse them without `-r`.
- **Output** — An output directory specified using `-o`. Use the same directory across all modules to reuse completed results.

```bash
MAGverify end-to-end \
  -a /path/to/assemblies \
  -r /path/to/hifi_reads \
  -g /path/to/MAGs \
  -o /path/to/MAGverify_output \
  -t 100 -m 900
```

Steps, in order: `map` → `circularity` → `clipping` → `false-closure` → `phantom` → `repeats` → `false-duplication` → `variants` → `chimera` → `summary`.

- **Before anything starts**, the preparation results are checked; if one is missing, the run stops at once and names the `prep-` command to run (or `--skip` to leave the module out).
- **Finished steps are reused.** Each step keeps fingerprints of its inputs, so an interrupted run is restarted with the same command and continues where it stopped.
- The preparation outputs are read from `-o DIR/00_prep/` automatically, including the GTDB database of `prep-chimera`. Only a database built under another `-o` needs `--database PREFIX` (for example `--database /path/to/other_output/00_prep/GTDB_R232`).
- **`--skip MODULES`** leaves out modules, for example `--skip chimera,variants`; `summary` then reports only the modules that ran.
- `-t` sets the threads and `-m` (GB) the memory limit of `map`, `circularity`, `repeats` and `variants`.
- **Progress and remaining time** are printed at the start of each step and every 5 minutes (`PROGRESS_INTERVAL=SECONDS` changes the interval): the progress (%), elapsed and remaining time of the current step and of the whole run. The estimate comes from the earlier run time of each step in the same `-o` (kept in `.cache/run_times.tsv`), so it is available from the second run on; read mapping (`map`, `circularity`) and the BLAST steps of `repeats` also report their own progress from the first run. Each step ends with its run time (`Elapsed`).
- `MAGVERIFY_DRY_RUN=1 MAGverify end-to-end ...` only prints the step commands.

Each module can also be run on its own, for example to rerun one step with other settings; see [Modules](docs/MODULES.md) for every command, what each module tests and its options.

## 4. Results

The summary is written in two versions at the top of the output directory:

- **`summary.MAG.tsv` (Y/N version) — start here:** one row per MAG, the number of checks that flagged it (`n_flagged_checks`) and one `Y`/`N`/`NA` column per check (`NA`: not assessable or not run; `Checks_not_run` names the checks that were not run). `summary.contig.tsv` has the same checks for every contig (`MAG_ID` is `NA` for unbinned contigs).
- **`summary_full.MAG.tsv` (full version):** the same columns, then key counts, then every column of every module MAG table in one row per MAG, named `MODULE.COLUMN` (for example `clipping.n_candidate_sites`). `summary_full.contig.tsv` does the same for every MAG contig.

Each module writes its tables to its own directory as `MODULE.LEVEL.tsv`, where `LEVEL` is the unit of one row. The modules themselves read the complete internal tables in each hidden `.cache/` directory, which also holds intermediate files and fingerprints that later runs reuse; do not delete it. Keep `01_mapping/` (BAMs, `.csi` indexes and `.cache/`) together when moving results.

| Module | Directory | Tables |
| --- | --- | --- |
| Preparation | `00_prep/` | `all_assemblies.fa`, `MAG_contig_map.tsv`, genes, `scg.contig.tsv`, `assembler_circularity.*.tsv`, `gtdbtk_taxonomy.MAG.tsv`, `taxonomy.gene.tsv` |
| Read mapping | `01_mapping/` | `mapping.assembler.tsv` (BAMs are in the same directory) |
| Circularity | `02_circularity/` | `circularity.MAG.tsv`, `circularity.contig.tsv` |
| Clipping | `03_clipping/` | `clipping.site.tsv`, `clipping.contig.tsv`, `clipping.MAG.tsv`, `clipping.assembler.tsv` |
| False closure | `04_false_closure/` | `false_closure.contig.tsv`, `false_closure.MAG.tsv` |
| Phantom sequences | `05_phantom/` | `phantom.region.tsv`, `phantom.contig.tsv`, `phantom.MAG.tsv` |
| Repeats | `06_repeats/` | `repeats.contig.tsv`, `repeats.MAG.tsv`, `repeats.region.tsv`, `repeats.HSP.tsv` |
| Taxonomic chimera | `07_chimera/` | `chimera.region.tsv`, `chimera.contig.tsv`, `chimera.MAG.tsv` |
| False duplication | `08_false_duplication/` | `false_duplication.pair.tsv`, `false_duplication.contig.tsv`, `false_duplication.MAG.tsv` |
| Assembly variants | `09_assembly_variants/` | `variants.site.tsv`, `variants.contig.tsv`, `variants.MAG.tsv` |
| Integrated QC summary | top of the output directory | `summary.MAG.tsv`, `summary.contig.tsv` (Y/N); `summary_full.MAG.tsv`, `summary_full.contig.tsv` (full) |

Every column of every table is defined in [Result tables](docs/OUTPUTS.md).

## Citation

If you use MAGverify, please cite:

Song MJ, Yun T, Yoon S, Lee JY, Pinto AJ, Yoon S. HiFi long-read metagenomics reshapes genome-resolved interpretations of activated sludge microbiome. *bioRxiv* (2026). https://doi.org/10.64898/2026.11.18.767492
