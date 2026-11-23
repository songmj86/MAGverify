# Input preparation

Run these commands before analysis from the activated MAGverify environment. Use the same `-o` for all commands so that subsequent modules can find their outputs in `00_prep/`. The `end-to-end` command does not run preparation steps.

For a summary of the preparation outputs and downstream modules, see the [README](../README.md#2-preparation).

## 1. GTDB R232 database

`prep-chimera` builds the MMseqs2 taxonomy database from three files in the [GTDB R232 (232.0) release](https://data.gtdb.ecogenomic.org/releases/release232/232.0/).

| Input file | Purpose |
| --- | --- |
| [`gtdb_proteins_aa_reps_r232.tar.gz`](https://data.gtdb.ecogenomic.org/releases/release232/232.0/genomic_files_reps/gtdb_proteins_aa_reps_r232.tar.gz) | Protein sequences (~123 GB compressed) |
| [`bac120_taxonomy_r232.tsv`](https://data.gtdb.ecogenomic.org/releases/release232/232.0/bac120_taxonomy_r232.tsv) | Bacterial taxonomy |
| [`ar53_taxonomy_r232.tsv`](https://data.gtdb.ecogenomic.org/releases/release232/232.0/ar53_taxonomy_r232.tsv) | Archaeal taxonomy |

Download all three files into one directory with `aria2c`:

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

The database is written as `00_prep/GTDB_R232*`. Building it needs substantial temporary disk space in `00_prep/.cache/prep_chimera_work/`. This temporary directory can be deleted after the build completes. To reuse a database built under another output directory, specify its prefix with `--database` when running `chimera` or `end-to-end`.

## 2. Single-copy gene annotation

`prep-scg` predicts genes with Pyrodigal and identifies single-copy genes (SCGs) using anvi'o (Bacteria_71 and Archaea_76).

```bash
MAGverify prep-scg -a /path/to/assemblies -o /path/to/MAGverify_output -t 100
```

- `-a DIR`: Input directory of contig files (`.fa`, `.fasta`, or `.fna`), one per assembler.
- Output: `00_prep/scg.contig.tsv` and reusable gene predictions (`genes.faa`, `genes.gff`, `genes.coordinates.tsv`).

The SCG table is used by `false-closure` and by the SCG annotation in `chimera`. The gene predictions are reused by `chimera` and `variants`.

## 3. Assembler circularity

`prep-circularity` collects the circularity calls reported by the assemblers. **This preparation step is mandatory for both `circularity` and `false-closure`.** Supply one `-c` file for each assembler:

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
| hifiasm-meta | `*.p_ctg.gfa` or `*.p_ctg.noseq.gfa` | Contig name suffix `c` (circular) or `l` (linear) |
| metaMDBG | `contigs.fasta.gz` | Header `circular=yes` or `circular=no` |
| metaFlye | `assembly_info.txt` | `circ.` column (`Y`/`N`) |

- `-g DIR`: Input directory of MAG files (one file per MAG).
- Output: `00_prep/assembler_circularity.MAG.tsv` and `assembler_circularity.contig.tsv`.

Provide the complete assembler outputs. Each MAG contig must be present, or the command reports the missing contigs. For hifiasm-meta, `*.p_ctg.noseq.gfa` is faster to read because only contig identifiers are needed. A compatible pre-existing MAG circularity table can also be supplied using `--assembler-circularity FILE` for individual module execution.

## 4. MAG taxonomy

`prep-taxonomy` converts existing GTDB-Tk `classify_wf` results into a MAG taxonomy table for `chimera`. It is optional: without MAG taxonomy, `chimera` runs without SCG redundancy annotation.

```bash
MAGverify prep-taxonomy \
  -d /path/to/GTDBTk_results \
  -g /path/to/MAGs \
  -o /path/to/MAGverify_output
```

- `-d DIR`: GTDB-Tk output directory with `gtdbtk.bac120.summary.tsv` and/or `gtdbtk.ar53.summary.tsv`.
- `-g DIR`: Input directory of MAG files (one file per MAG).
- Output: `00_prep/gtdbtk_taxonomy.MAG.tsv`.

MAG IDs must match the FASTA filenames without their extensions. An equivalent taxonomy table can be supplied to `chimera` with `--mag-taxonomy TSV`:

```text
ID	GTDB_Domain	GTDB_Phylum	GTDB_Class	GTDB_Order
ExampleMAG	d__Bacteria	p__Nitrospirota	c__Nitrospiria	o__Nitrospirales
```

Unassigned ranks use their bare prefix (for example, `o__`). If GTDB-Tk results are not available, generate them first with [GTDB-Tk](https://github.com/ecogenomics/gtdbtk), which is installed separately from MAGverify.
