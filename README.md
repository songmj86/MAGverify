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

| Option | Input |
| --- | --- |
| `-a DIR` | The original, unbinned assembly FASTA of each assembler (.fa, .fasta or .fna, uncompressed). Do not give MAG FASTAs here: reads are mapped to the whole assemblies so that alignments to unbinned contigs are kept. Contig IDs (the first word after `>`) must be unique across the files. |
| `-g DIR` | MAG FASTA files, one per MAG. The file name without `.fa`, `.fna` or `.fasta` is the MAG ID. MAG contig IDs and sequences must match the assembly files. |
| `-r DIR` | PacBio HiFi reads (.fastq, .fq, .fastq.gz or .fq.gz). Several files are joined once into `00_prep/hifi_reads.fastq[.gz]`. |
| `-c FILE` | The circularity output of each assembler, one `-c` per assembler (hifiasm-meta `p_ctg.gfa`, metaMDBG `contigs.fasta.gz`, metaFlye `assembly_info.txt`). |
| `-d DIR` | GTDB-Tk `classify_wf` results of the MAGs (optional; adds SCG redundancy to `chimera`). |
| `-o DIR` | Output directory. Every command uses the same `-o`, so finished results are reused. |

The GTDB R232 MMseqs2 database for `chimera` is built once with `MAGverify prep-chimera` from three GTDB R232 downloads:

| Download | Content |
| --- | --- |
| [`gtdb_proteins_aa_reps_r232.tar.gz`](https://data.gtdb.ecogenomic.org/releases/release232/232.0/genomic_files_reps/gtdb_proteins_aa_reps_r232.tar.gz) | Predicted amino-acid sequences from GTDB representative genomes (~123 GB compressed) |
| [`bac120_taxonomy_r232.tsv`](https://data.gtdb.ecogenomic.org/releases/release232/232.0/bac120_taxonomy_r232.tsv) | Bacterial genome accession → full GTDB lineage |
| [`ar53_taxonomy_r232.tsv`](https://data.gtdb.ecogenomic.org/releases/release232/232.0/ar53_taxonomy_r232.tsv) | Archaeal genome accession → full GTDB lineage |

```bash
MAGverify prep-chimera -d /path/to/GTDB_R232_downloads -o /path/to/database_dir -t 100
```

Download commands and the format of each preparation table are in [Input preparation](docs/INPUT_PREPARATION.md).

## 3. Run

One command runs every step in order with the same inputs:

```bash
conda activate MAGverify

MAGverify end-to-end \
  -a /path/to/assemblies \
  -r /path/to/hifi_reads \
  -g /path/to/MAGs \
  -c /path/to/hifiasm_meta/Read.p_ctg.gfa \
  -c /path/to/metaMDBG/contigs.fasta.gz \
  -c /path/to/metaflye/assembly_info.txt \
  -d /path/to/GTDBTk_results \
  --database /path/to/database_dir/00_prep/GTDB_R232 \
  -o /path/to/MAGverify_output \
  -t 100 -m 900
```

Steps, in order: `prep-circularity` → `prep-taxonomy` → `prep-scg` → `map` → `circularity` → `clipping` → `false-closure` → `phantom` → `repeats` → `false-duplication` → `variants` → `chimera` → `summary`.

- **Before anything starts**, the inputs are checked: the run stops at once if a step would fail later (for example, no `-c` and no assembler circularity table for `false-closure`, or no MMseqs2 database for `chimera`).
- **Finished steps are reused.** Each step keeps fingerprints of its inputs, so an interrupted run is restarted with the same command and continues where it stopped.
- **`-c` and `-d` are needed only once:** later runs with the same `-o` reuse `00_prep/assembler_circularity.MAG.tsv` and `00_prep/gtdbtk_taxonomy.MAG.tsv`. `--database` can be left out when the database was built with the same `-o`.
- **`--skip MODULES`** leaves out modules, for example `--skip chimera,variants`; `summary` then reports only the modules that ran.
- `-t` sets the threads and `-m` (GB) the memory limit of `map`, `circularity`, `repeats` and `variants`.
- Each step prints its commands and run time (`Elapsed`). `MAGVERIFY_DRY_RUN=1 MAGverify end-to-end ...` only prints the step commands.

Each module can also be run on its own, for example to rerun one step with other settings; see [Modules](docs/MODULES.md) for every command, what each module tests and its options.

## 4. Results

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
