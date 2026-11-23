# MAGverify

MAGverify assesses the integrity of metagenome assemblies and metagenome-assembled genomes (MAGs) using PacBio HiFi reads. Its modules cover read mapping, circularity, clipping, false closure, phantom sequences, repeats, taxonomic chimeras, false duplication, assembly variants and an integrated summary.

## 1. Installation

Install MAGverify on Linux using Conda or Mamba from the repository root:

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

Run the preparation commands before analysis using the same `-o` to save their outputs in `00_prep/`. They are not included in `end-to-end`.

| Command | Output in `00_prep/` | Used by |
| --- | --- | --- |
| [`prep-chimera`](docs/USAGE.md#gtdb-r232-database) | `GTDB_R232*` | `chimera` |
| [`prep-scg`](docs/USAGE.md#single-copy-gene-annotation) | `scg.contig.tsv`, `genes.*` | `false-closure`, `chimera` |
| [`prep-circularity`](docs/USAGE.md#assembler-circularity) | `assembler_circularity.MAG.tsv`, `assembler_circularity.contig.tsv` | `circularity`, `false-closure` |
| [`prep-taxonomy`](docs/USAGE.md#mag-taxonomy) | `gtdbtk_taxonomy.MAG.tsv` | `chimera` |

See [Usage](docs/USAGE.md#preparation) for commands, download links and input formats.

## 3. Quick start

After preparation, run all analysis modules using end-to-end with the following options:

- `-a DIR`: Input directory of contig files (`.fa`, `.fasta`, or `.fna`). Contig IDs must be unique across assembly files.
- `-g DIR`: Input directory of MAG files (one file per MAG). MAG contig IDs must match those in the assembly files.
- `-r DIR`: Input directory of PacBio HiFi read files (`.fastq`, `.fq`, `.fastq.gz`, or `.fq.gz`). Multiple files are merged into `00_prep/hifi_reads.fastq[.gz]` and reused by subsequent modules without `-r`.
- `-o DIR`: Output directory. Use the same directory across all modules to reuse completed results.
- `-t N`: Number of threads.
- `-m N`: Memory limit (GB) for memory-managed modules.
- `--skip MODULES`: Comma-separated modules to exclude (e.g., `--skip chimera,variants`).


```bash
MAGverify end-to-end \
  -a /path/to/assemblies \
  -r /path/to/hifi_reads \
  -g /path/to/MAGs \
  -o /path/to/MAGverify_output \
  -t 100 -m 900
```

**Note:** `end-to-end` runs the modules in the following order: `map` → `circularity` → `clipping` → `false-closure` → `phantom` → `repeats` → `false-duplication` → `variants` → `chimera` → `summary`. Before analysis, MAGverify checks for missing preparation outputs in `00_prep/` and identifies the necessary `prep-` commands. Interrupted analyses can be resumed by rerunning the same command, automatically skipping completed modules with unchanged inputs.

Each module can also be run separately. See [Usage](docs/USAGE.md#analysis-modules) for detailed descriptions and options.


## 4. Results

MAGverify generates two summary formats: `summary.MAG.tsv` and `summary.contig.tsv` provide `Y`/`N`/`NA` flags for each analysis, whereas `summary_detail.MAG.tsv` and `summary_detail.contig.tsv` additionally include detailed results from each module.

Each module saves its results in a separate directory, as listed below. Intermediate files are stored in hidden `.cache/` directories for reuse in subsequent runs. Keep these directories to enable analysis resumption.

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
| Integrated summary | top of the output directory | `summary.MAG.tsv`, `summary.contig.tsv` (Y/N); `summary_detail.MAG.tsv`, `summary_detail.contig.tsv` (detailed) |

See [Result tables](docs/OUTPUTS.md) for detailed descriptions of all output columns.


## Citation

If you use MAGverify, please cite:

Song MJ, Yun T, Yoon S, Lee JY, Pinto AJ, Yoon S. HiFi long-read metagenomics reshapes genome-resolved interpretations of activated sludge microbiome. *bioRxiv* (2026). https://doi.org/10.64898/2026.11.18.767492
