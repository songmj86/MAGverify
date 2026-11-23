# Input preparation

MAGverify uses four preparation commands to build or collect the inputs needed by its analysis modules. Run the applicable commands from the activated MAGverify environment **before** `end-to-end`. Use the same `-o` for preparation and analysis so that all outputs are saved in `00_prep/` and can be reused.

| Command | Main output in `00_prep/` | Used by |
| --- | --- | --- |
| [`prep-chimera`](#1-gtdb-r232-database) | `GTDB_R232*` | `chimera` |
| [`prep-scg`](#2-single-copy-gene-annotation) | `scg.contig.tsv`, `genes.*` | `false-closure`, `chimera` |
| [`prep-circularity`](#3-assembler-circularity) | `assembler_circularity.MAG.tsv`, `assembler_circularity.contig.tsv` | `circularity`, `false-closure` |
| [`prep-taxonomy`](#4-mag-taxonomy) | `gtdbtk_taxonomy.MAG.tsv` | `chimera` (optional annotation) |

The `end-to-end` command does not run these preparation steps. See the [README](../README.md#2-preparation) for the overall workflow and [Modules](MODULES.md) for the analysis commands.

## 1. GTDB R232 database

`prep-chimera` builds an MMseqs2 taxonomy database from three files in the [GTDB R232 (232.0) release](https://data.gtdb.ecogenomic.org/releases/release232/232.0/).

| Input file | Purpose |
| --- | --- |
| [`gtdb_proteins_aa_reps_r232.tar.gz`](https://data.gtdb.ecogenomic.org/releases/release232/232.0/genomic_files_reps/gtdb_proteins_aa_reps_r232.tar.gz) | Protein sequences (~123 GB compressed) |
| [`bac120_taxonomy_r232.tsv`](https://data.gtdb.ecogenomic.org/releases/release232/232.0/bac120_taxonomy_r232.tsv) | Bacterial taxonomy |
| [`ar53_taxonomy_r232.tsv`](https://data.gtdb.ecogenomic.org/releases/release232/232.0/ar53_taxonomy_r232.tsv) | Archaeal taxonomy |

Download all three files into the same directory with `aria2c`:

```bash
GTDB_DIR=/path/to/GTDB_R232_downloads
BASE=https://data.gtdb.ecogenomic.org/releases/release232/232.0
mkdir -p "$GTDB_DIR"

aria2c -c -x 4 -s 4 -d "$GTDB_DIR" "$BASE/genomic_files_reps/gtdb_proteins_aa_reps_r232.tar.gz"
aria2c -c -d "$GTDB_DIR" "$BASE/bac120_taxonomy_r232.tsv"
aria2c -c -d "$GTDB_DIR" "$BASE/ar53_taxonomy_r232.tsv"
```

Build the database:

```bash
MAGverify prep-chimera -d "$GTDB_DIR" -o /path/to/MAGverify_output -t 100
```

**Output:** `00_prep/GTDB_R232*` (MMseqs2 database files).

**Disk space:** The build needs substantial temporary storage in `00_prep/.cache/prep_chimera_work/`. This temporary working directory can be deleted after the build completes. To reuse an existing database from another output directory, pass its prefix using `--database` with `chimera` or `end-to-end`.

## 2. Single-copy gene annotation

`prep-scg` predicts genes with Pyrodigal and identifies single-copy genes (SCGs) using anvi'o's `Bacteria_71` and `Archaea_76` collections.

```bash
MAGverify prep-scg -a /path/to/assemblies -o /path/to/MAGverify_output -t 100
```

- `-a DIR`: Directory of assembly contig files (`.fa`, `.fasta`, or `.fna`), one per assembler.
- **Output:** `00_prep/scg.contig.tsv` and reusable gene predictions (`genes.faa`, `genes.gff`, `genes.coordinates.tsv`).

The SCG table is used by `false-closure` and for SCG annotation in `chimera`. The gene predictions are reused by `chimera` and `variants`.

## 3. Assembler circularity

`prep-circularity` collects assembler-reported circularity calls. **This step is mandatory for `circularity` and `false-closure`.** Provide one `-c` file per assembler.

```bash
MAGverify prep-circularity \
  -c /path/to/hifiasm_meta/Read.p_ctg.gfa \
  -c /path/to/metaMDBG/contigs.fasta.gz \
  -c /path/to/metaFlye/assembly_info.txt \
  -g /path/to/MAGs \
  -o /path/to/MAGverify_output
```

| Assembler | Input (`-c`) | Circularity information |
| --- | --- | --- |
| hifiasm-meta | `*.p_ctg.gfa` or `*.p_ctg.noseq.gfa` | Contig-name suffix `c` (circular) or `l` (linear) |
| metaMDBG | `contigs.fasta.gz` | Header `circular=yes` or `circular=no` |
| metaFlye | `assembly_info.txt` | `circ.` column (`Y` / `N`) |

- `-g DIR`: Directory of MAG FASTA files, one per MAG.
- **Output:** `00_prep/assembler_circularity.MAG.tsv` and `00_prep/assembler_circularity.contig.tsv`.

Provide complete assembler outputs so every MAG contig can be matched. Missing contigs are reported by the command. For hifiasm-meta, `*.p_ctg.noseq.gfa` is faster to read because only contig IDs are needed. For individual module execution, a compatible existing MAG circularity table can instead be supplied with `--assembler-circularity FILE`.

## 4. MAG taxonomy

`prep-taxonomy` converts GTDB-Tk `classify_wf` results to a MAG taxonomy table for `chimera`. **This step is optional:** without MAG taxonomy, `chimera` runs without SCG redundancy annotation.

```bash
MAGverify prep-taxonomy \
  -d /path/to/GTDBTk_results \
  -g /path/to/MAGs \
  -o /path/to/MAGverify_output
```

- `-d DIR`: GTDB-Tk output directory with `gtdbtk.bac120.summary.tsv` and/or `gtdbtk.ar53.summary.tsv`.
- `-g DIR`: Directory of MAG FASTA files, one per MAG.
- **Output:** `00_prep/gtdbtk_taxonomy.MAG.tsv`.

MAG IDs must match the FASTA filenames without their extensions. Alternatively, pass an existing table to `chimera` using `--mag-taxonomy TSV`:

```text
ID	GTDB_Domain	GTDB_Phylum	GTDB_Class	GTDB_Order
ExampleMAG	d__Bacteria	p__Nitrospirota	c__Nitrospiria	o__Nitrospirales
```

Use a bare rank prefix for an unassigned rank (for example, `o__`). If GTDB-Tk results are unavailable, generate them first with [GTDB-Tk](https://github.com/ecogenomics/gtdtk), which is installed separately from MAGverify.
