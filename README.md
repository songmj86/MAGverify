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

The environment includes anvi'o 9 and the other programs used by MAGverify (versions pinned in `environment.yml`) and registers the `MAGverify` command from this directory. `MAGverify doctor` shows the code directory in use and checks every required program and its version; it ends with `All 22 required programs found.`

- **Updating:** update the code in this directory (`git pull`, or copy new files over it) instead of deleting and recreating it. After `environment.yml` changes, run `mamba env update -f environment.yml --prune`.
- **If the directory was moved**, or `MAGverify` reports `ModuleNotFoundError: No module named 'magverify_cli'`, register the current directory once with `pip install -e . --no-deps --no-build-isolation` (run inside the activated environment, from the repository root).

## 2. Input data

| Input | Content |
| --- | --- |
| Assemblies (`-a DIR`) | The original, unbinned assembly FASTA of each assembler (.fa, .fasta or .fna, uncompressed). Do not give MAG FASTAs here: reads are mapped to the whole assemblies so that alignments to unbinned contigs are kept. Contig IDs (the first word after `>`) must be unique across the files. |
| MAGs (`-g DIR`) | MAG FASTA files, one per MAG. The file name without `.fa`, `.fna` or `.fasta` is the MAG ID. MAG contig IDs and sequences must match the assembly files. |
| HiFi reads (`-r DIR`) | PacBio HiFi reads (.fastq, .fq, .fastq.gz or .fq.gz). Several files are joined once into `00_prep/hifi_reads.fastq[.gz]`. |
| Output (`-o DIR`) | Output directory. Use the **same `-o`** for the preparation and analysis commands: each command finds the results of the earlier ones there. |

The preparation step also uses the circularity output of each assembler, GTDB-Tk results of the MAGs (optional) and three GTDB R232 downloads.

## 3. Preparation (run first, separately)

Run these commands before the analysis, with the same `-o`. `prep-chimera` builds the GTDB database once and can be reused by other projects through `--database`.

```bash
# MMseqs2 database of GTDB R232 (once; for chimera)
MAGverify prep-chimera -d /path/to/GTDB_R232_downloads -o /path/to/MAGverify_output -t 100

# Single-copy genes with anvi'o (for false-closure and chimera)
MAGverify prep-scg -a /path/to/assemblies -o /path/to/MAGverify_output -t 100

# Circularity reported by each assembler, one -c per assembler (for circularity and false-closure)
MAGverify prep-circularity \
  -c /path/to/hifiasm_meta/Read.p_ctg.gfa \
  -c /path/to/metaMDBG/contigs.fasta.gz \
  -c /path/to/metaflye/assembly_info.txt \
  -g /path/to/MAGs -o /path/to/MAGverify_output

# MAG taxonomy from GTDB-Tk (optional; adds SCG redundancy to chimera)
MAGverify prep-taxonomy -d /path/to/GTDBTk_results -g /path/to/MAGs -o /path/to/MAGverify_output
```

| Command | Output in `00_prep/` | Used by |
| --- | --- | --- |
| `prep-chimera` | `GTDB_R232*` (MMseqs2 database) | `chimera` |
| `prep-scg` | `scg.contig.tsv`, `genes.*` | `false-closure`, `chimera` |
| `prep-circularity` | `assembler_circularity.MAG.tsv`, `assembler_circularity.contig.tsv` | `circularity` (optional), `false-closure` |
| `prep-taxonomy` | `gtdbtk_taxonomy.MAG.tsv` | `chimera` (optional) |

`prep-chimera` reads three files from **GTDB R232**:

| Download | Content |
| --- | --- |
| [`gtdb_proteins_aa_reps_r232.tar.gz`](https://data.gtdb.ecogenomic.org/releases/release232/232.0/genomic_files_reps/gtdb_proteins_aa_reps_r232.tar.gz) | Predicted amino-acid sequences from GTDB representative genomes (~123 GB compressed) |
| [`bac120_taxonomy_r232.tsv`](https://data.gtdb.ecogenomic.org/releases/release232/232.0/bac120_taxonomy_r232.tsv) | Bacterial genome accession → full GTDB lineage |
| [`ar53_taxonomy_r232.tsv`](https://data.gtdb.ecogenomic.org/releases/release232/232.0/ar53_taxonomy_r232.tsv) | Archaeal genome accession → full GTDB lineage |

Download commands, supported assembler files and the format of each table are in [Input preparation](docs/INPUT_PREPARATION.md).

## 4. Analysis (end-to-end)

After the preparation, one command runs all analysis modules in order:

```bash
MAGverify end-to-end \
  -a /path/to/assemblies \
  -r /path/to/hifi_reads \
  -g /path/to/MAGs \
  -o /path/to/MAGverify_output \
  --database /path/to/MAGverify_output/00_prep/GTDB_R232 \
  -t 100 -m 900
```

Steps, in order: `map` → `circularity` → `clipping` → `false-closure` → `phantom` → `repeats` → `false-duplication` → `variants` → `chimera` → `summary`.

- **Before anything starts**, the preparation results are checked; if one is missing, the run stops at once and names the `prep-` command to run (or `--skip` to leave the module out).
- **Finished steps are reused.** Each step keeps fingerprints of its inputs, so an interrupted run is restarted with the same command and continues where it stopped.
- `--database` can be left out when `prep-chimera` used the same `-o`.
- **`--skip MODULES`** leaves out modules, for example `--skip chimera,variants`; `summary` then reports only the modules that ran.
- `-t` sets the threads and `-m` (GB) the memory limit of `map`, `circularity`, `repeats` and `variants`.
- Each step prints its commands and run time (`Elapsed`). `MAGVERIFY_DRY_RUN=1 MAGverify end-to-end ...` only prints the step commands.

Each module can also be run on its own, for example to rerun one step with other settings; see [Modules](docs/MODULES.md) for every command, what each module tests and its options.

## 5. Results

**Start with `summary.MAG.tsv`** at the top of the output folder: one row per MAG, the number of checks that flagged it (`n_flagged_checks`), then one `Y`/`N`/`NA` column per check (`NA`: not assessable or not run; `Checks_not_run` names the checks that were not run), then counts. `summary.contig.tsv` has the same checks for every contig (`MAG_ID` is `NA` for unbinned contigs).

Each module writes its tables to its own folder as `MODULE.LEVEL.tsv`, where `LEVEL` is the unit of one row. The modules themselves read the full internal tables in each hidden `.cache/` folder, which also holds intermediate files and fingerprints that later runs reuse; do not delete it. Keep `01_mapping/` (BAMs, `.csi` indexes and `.cache/`) together when moving results.

| Module | Folder | Tables |
| --- | --- | --- |
| Preparation | `00_prep/` | `all_assemblies.fa`, `MAG_contig_map.tsv`, genes, `scg.contig.tsv`, `assembler_circularity.*.tsv`, `gtdbtk_taxonomy.MAG.tsv`, `taxonomy.gene.tsv` |
| Read mapping | `01_mapping/` | `mapping.assembler.tsv` (BAMs are in the same folder) |
| Circularity | `02_circularity/` | `circularity.MAG.tsv`, `circularity.contig.tsv` |
| Clipping | `03_clipping/` | `clipping.site.tsv`, `clipping.contig.tsv`, `clipping.MAG.tsv`, `clipping.assembler.tsv` |
| False closure | `04_false_closure/` | `false_closure.contig.tsv`, `false_closure.MAG.tsv` |
| Phantom sequences | `05_phantom/` | `phantom.region.tsv`, `phantom.contig.tsv`, `phantom.MAG.tsv` |
| Repeats | `06_repeats/` | `repeats.contig.tsv`, `repeats.MAG.tsv`, `repeats.region.tsv`, `repeats.HSP.tsv` |
| Taxonomic chimera | `07_chimera/` | `chimera.region.tsv`, `chimera.contig.tsv`, `chimera.MAG.tsv` |
| False duplication | `08_false_duplication/` | `false_duplication.pair.tsv`, `false_duplication.contig.tsv`, `false_duplication.MAG.tsv` |
| Assembly variants | `09_assembly_variants/` | `variants.site.tsv`, `variants.contig.tsv`, `variants.MAG.tsv` |
| Integrated QC summary | top of the output folder | `summary.MAG.tsv`, `summary.contig.tsv` |

Every column of every table is defined in [Result tables](docs/OUTPUTS.md).

## Citation

If you use MAGverify, please cite:

Song MJ, Yun T, Yoon S, Lee JY, Pinto AJ, Yoon S. HiFi long-read metagenomics reshapes genome-resolved interpretations of activated sludge microbiome. *bioRxiv* (2026). https://doi.org/10.64898/2026.11.18.767492
