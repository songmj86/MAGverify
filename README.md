# MAGverify

MAGverify assesses the integrity of metagenome assemblies and metagenome-assembled genomes (MAGs) using PacBio HiFi reads. Its numbered modules follow the assembly-validation methods: read mapping, circularity, clipping, false closure, phantom sequences, repeats, taxonomic chimeras, false duplication, assembly variants and integrated QC.

## Installation

Linux and Conda/Mamba are supported. One Conda environment installs MAGverify, anvi'o v9 and the command-line tools needed by its modules.

```bash
git clone https://github.com/songmj86/MAGverify.git
cd MAGverify
mamba env create -f environment.yml -n MAGverify
conda activate MAGverify
MAGverify doctor
MAGverify list
```

**Run the creation command from the repository root**, since `environment.yml` installs the local MAGverify Python package via its `pip:` section. See [INSTALL.md](docs/INSTALL.md) for package checks and platform limitations.

## Database preparation

Prepare external inputs once under `/path/to/MAGverify_output/0_prep/`. Full fixed-release download links, required source files, preparation scripts, taxonomy TSV format and GTDB-Tk workflow appear in [PROTOCOL.md](PROTOCOL.md#database-preparation-once-before-analysis).

**GTDB R232 MMseqs2 protein taxonomy DB** (required for chimera unless `--taxonomy FILE` is supplied):

```bash
conda activate MAGverify
MAGverify prep-chimera \
  --output-dir /path/to/MAGverify_output \
  --threads 100
```

This downloads the fixed R232 representative **protein** archive (`gtdb_proteins_aa_reps_r232.tar.gz`) plus `bac120_taxonomy_r232.tsv` and `ar53_taxonomy_r232.tsv` from the official [GTDB R232 release](https://data.gtdb.ecogenomic.org/releases/release232/232.0/), and creates `$PREP/GTDB_R232` (MMseqs2 taxonomy DB prefix). Plan for large storage needs; the protein download alone is ~123 GB compressed.

**MAG taxonomy TSV:** Prepare `$PREP/MAG_taxonomy.tsv` with `ID`, `GTDB_Domain`, `GTDB_Phylum`, `GTDB_Class`, `GTDB_Order`. Use an existing table or convert existing GTDB-Tk `classify_wf` summary TSVs without rerunning classification:

```bash
python -m scripts.prep_gtdbtk_summary_to_mag_taxonomy \
  --gtdbtk-dir /path/to/GTDBTk_results \
  --mag-dir /path/to/MAGs \
  --output "$PREP/MAG_taxonomy.tsv"
```

If MAG classification is not yet available, [download the separate GTDB-Tk R232 reference package](https://data.gtdb.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz) and run `gtdbtk classify_wf` as shown in [PROTOCOL.md](PROTOCOL.md#c-gtdb-tk-r232-classification--only-if-mag-taxonomy-is-unavailable). No change to MAGverify analysis options or internal taxonomy parser is needed.

## Usage

### Input data

Use paths to existing data. Do not create a special directory tree or prepare `combined.fasta` yourself.

- `--assemblies-dir` — directory with original source-assembly FASTA files (unique contig IDs across files)
- `--mag-dir` — directory with one `.fa`, `.fasta` or `.fna` file per MAG
- `--reads` — PacBio HiFi reads (`.fastq` or `.fastq.gz`)
- `--output-dir` — common results directory, reused across separately executed modules (`--output-directory` remains accepted for compatibility)

Gene-level taxonomy for `chimera` also needs a GTDB MMseqs2 taxonomy database (`--database`) or an existing gene taxonomy table (`--taxonomy`). False closure additionally needs the assembler-reported circularity table (`--assembler-circularity`).

### Run modules individually

Run the ten steps in the order shown in [PROTOCOL.md](PROTOCOL.md). `MAGverify summary` includes final-result validation; no separate `check` module is necessary. Each command is independent and can reuse prior results saved in the same output directory.

**Clipping**

```bash
MAGverify clipping \
  --assemblies-dir /path/to/assemblies \
  --mag-dir /path/to/MAGs \
  --reads /path/to/hifi.fastq.gz \
  --output-dir /path/to/MAGverify_output \
  --threads 100
```

**Taxonomic chimera** — generates gene taxonomy internally if a compatible result does not already exist; SCG annotation is also prepared when necessary and when anvi'o is available.

```bash
MAGverify chimera \
  --assemblies-dir /path/to/assemblies \
  --mag-dir /path/to/MAGs \
  --database /path/to/MAGverify_output/0_prep/GTDB_R232 \
  --mag-taxonomy /path/to/MAGverify_output/0_prep/MAG_taxonomy.tsv \
  --output-dir /path/to/MAGverify_output \
  --threads 100
```

Chimera classification also uses the validated clipping table from the preceding clipping step. The MMseqs2 taxonomy database is a **database prefix** for which a `.dbtype` file exists. `--mag-taxonomy` supplies genome-level GTDB classification for SCG annotation, not the gene taxonomy database.

**Show available modules**

```bash
MAGverify list
MAGverify chimera --help
MAGverify doctor
```

`MAGverify prepare`, `MAGverify structure`, `MAGverify run` and similar grouped commands are not part of this interface. The preprocessing that is specific to chimera (gene taxonomy, SCG annotation) belongs to that analysis, rather than being separate user-run commands.

The `map` module creates `01_mapping/assembler_mapping_outputs.tsv` and publishes `results/mapping_manifest.tsv` as a symbolic link in `--output-dir`. Later modules automatically reuse the mapping, including when `--reads` is also supplied. `--assemblies-dir` refers to a directory of original FASTA files, not a single merged FASTA. See [PROTOCOL.md](PROTOCOL.md) for command examples.

## Results

Final tables are linked under `OUTPUT/results/`. Detailed intermediate files remain in their original module directories for diagnostics and reproducibility.

| Analysis | Main output |
| --- | --- |
| Mapping | `mapping_manifest.tsv` |
| Clipping | `clipping_sites.tsv`, `clipping_mag.tsv` |
| Chimera (including taxonomy and SCG) | `chimera_contig.tsv`, `chimera_mag.tsv`, `gene_taxonomy.tsv`, `scg.tsv` (when generated) |
| Circularity | `circularity_mag.tsv` |
| False closure | `false_closure_mag.tsv`, `false_closure_contig.tsv` |
| Assembly variants | `variants_contig.tsv`, `variants_mag.tsv` |
| Phantom sequences | `phantom_contig.tsv`, `phantom_mag.tsv` |
| Repeats | `repeats_contig.tsv`, `repeats_mag.tsv` |
| False duplication | `false_duplication_contig.tsv`, `false_duplication_mag.tsv` |
| Integrated QC | `summary_contig.tsv`, `summary_mag.tsv` |

The short filenames are symbolic links to the generated tables, not separate data copies. Keep the results directory and its analysis subdirectories together when moving the output.

To inspect the commands without installing, run `./bin/MAGverify --help` from the repository root.

## Configuration

**No `00_config.sh` or `config/modules/` setup is needed.** Data locations, output paths and the number of threads are supplied as command-line options. Internal Shell defaults are in `modules/lib_config.sh` and `modules/stage_paths.sh`; users do not edit these files. Summary-table validation is integrated into `MAGverify summary`. It checks the resulting tables without repeating the analysis-stage tool checks. Gene taxonomy preprocessing is integrated into `modules/07_chimera.sh`; `modules/lib_scg.sh` supplies SCG annotations as an internal helper.

All protocol examples request `--threads 100`. Use a smaller allocation when fewer than 100 CPU threads are available.

## Documentation

- [Installation](docs/INSTALL.md)
- [Step-by-step protocol](PROTOCOL.md)
- [Module inputs and outputs](docs/MODULES.md)
- [Project architecture](docs/ARCHITECTURE.md)
- [Input formats](docs/INPUT_FORMATS.md)
- [Recorded tool versions](docs/TOOL_VERSIONS.md)

The command-line routing, syntax and selected fixture tests were checked; the full pipeline has **not** been rerun with real HiFi data. Compare results against the original analyses before publishing a scientific benchmark.

## Development and reproducibility

The repository separates the installable `modules/` and `scripts/` directories, validation tests, documentation, and environment definitions. Run `python -m unittest discover -s tests -p "test_*.py"` and the three shell tests before a release. The bundled tests exercise input handling, CLI behavior, and summary consistency with synthetic fixtures. A complete end-to-end verification against real PacBio HiFi datasets is not yet established.

The user-facing command and project name are both **MAGverify**. Python helpers live in the top-level `scripts/` directory and Bash modules in `modules/`. No dataset-specific paths belong in the GitHub repository.

A software license and author-approved citation should be added before describing this repository as openly licensed or publication-ready.

Python helper scripts are grouped by module prefix under `scripts/`; see [Script index](docs/SCRIPT_INDEX.md).
