# Input preparation

Activate the `MAGverify` Conda environment before running the preparation steps. These steps build the reference data and annotations that the analysis modules read; the analysis modules do not download data or run these steps themselves.

| Step | Command | Output | Used by |
| --- | --- | --- | --- |
| 1 | `MAGverify prep-chimera` | `DB_DIR/GTDB_R232*` | `chimera --database` |
| 2 | `prep_gtdbtk_summary_to_mag_taxonomy.py` | `DB_DIR/MAG_taxonomy.tsv` | `chimera --mag-taxonomy` |
| 3 | `MAGverify prep-scg` | `OUTPUT/00_prep/SCG/` | `false-closure`, `chimera` (found automatically) |
| 4 | `MAGverify prep-circularity` | `OUTPUT/00_prep/assembler_circularity.MAG.tsv` | `circularity` (optional), `false-closure` (found automatically) |

Steps 1 and 2 write to a database directory (`DB_DIR`) that can be reused across projects. Steps 3 and 4 write to the analysis output directory (`OUTPUT`, the same `-o` as the analysis modules).

## 1. GTDB R232 MMseqs2 protein taxonomy database

Skip this step if an existing compatible gene taxonomy table will be supplied with `--taxonomy`.

Download the three files from the fixed **GTDB R232 (232.0)** release:

| Download | Content |
| --- | --- |
| [`gtdb_proteins_aa_reps_r232.tar.gz`](https://data.gtdb.ecogenomic.org/releases/release232/232.0/genomic_files_reps/gtdb_proteins_aa_reps_r232.tar.gz) | Predicted amino-acid sequences from GTDB representative genomes (~123 GB compressed) |
| [`bac120_taxonomy_r232.tsv`](https://data.gtdb.ecogenomic.org/releases/release232/232.0/bac120_taxonomy_r232.tsv) | Bacterial genome accession → full GTDB lineage |
| [`ar53_taxonomy_r232.tsv`](https://data.gtdb.ecogenomic.org/releases/release232/232.0/ar53_taxonomy_r232.tsv) | Archaeal genome accession → full GTDB lineage |

**Download with aria2c** into one input directory (the commands only download the files):

```bash
GTDB_DIR=/path/to/GTDB_R232_files     # input: the three GTDB R232 files
BASE=https://data.gtdb.ecogenomic.org/releases/release232/232.0
mkdir -p "$GTDB_DIR"

aria2c -c -x 4 -s 4 -d "$GTDB_DIR" "$BASE/genomic_files_reps/gtdb_proteins_aa_reps_r232.tar.gz"
aria2c -c -d "$GTDB_DIR" "$BASE/bac120_taxonomy_r232.tsv"
aria2c -c -d "$GTDB_DIR" "$BASE/ar53_taxonomy_r232.tsv"
```

**Build the MMseqs2 taxonomy database.** Steps 1 and 2 write to the same database directory (`-o`):

```bash
DB_OUT=/path/to/MAGverify_db          # database directory: shared by steps 1 and 2

MAGverify prep-chimera \
  --gtdb-dir "$GTDB_DIR" \
  -o "$DB_OUT" \
  -t 100

ls -lh "$DB_OUT"/GTDB_R232.dbtype "$DB_OUT"/GTDB_R232_mapping "$DB_OUT"/GTDB_R232_taxonomy
```

| Option | Input or output |
| --- | --- |
| `--gtdb-dir DIR` | Directory containing `gtdb_proteins_aa_reps_r232.tar.gz`, `bac120_taxonomy_r232.tsv` and `ar53_taxonomy_r232.tsv` |
| `-o DIR` | Output directory. The database is written as `DIR/GTDB_R232*`; temporary files go to `DIR/tmp_prep_chimera/` |

Use `DIR/GTDB_R232` as `MAGverify chimera --database`. `DIR/tmp_prep_chimera/` needs substantially more disk space than the compressed archive and can be deleted once the database is built.

## 2. MAG taxonomy TSV

For domain-matched SCG annotation, `MAGverify chimera --mag-taxonomy` takes a tab-separated table with one row per MAG:

```text
ID	GTDB_Domain	GTDB_Phylum	GTDB_Class	GTDB_Order
ExampleMAG	d__Bacteria	p__Nitrospirota	c__Nitrospiria	o__Nitrospirales
```

`ID` must match the MAG FASTA file name without its extension (`ExampleMAG.fa` → `ExampleMAG`). Unassigned ranks keep the bare prefix (for example, `o__`). A compatible table can be used directly.

To build it from existing GTDB-Tk results, write it into the same output directory as step 1:

```bash
python /path/to/MAGverify/scripts/prep_gtdbtk_summary_to_mag_taxonomy.py \
  --gtdbtk-dir /path/to/GTDBTk_results \
  -g /path/to/MAGs \
  -o "$DB_OUT"
```

| Option | Input or output |
| --- | --- |
| `--gtdbtk-dir DIR` | GTDB-Tk `classify_wf` output directory containing `gtdbtk.bac120.summary.tsv` and/or `gtdbtk.ar53.summary.tsv` |
| `-g DIR` | Directory of MAG FASTA files (`.fa`, `.fna`, `.fasta`; file name = MAG ID). Only these MAGs are written; other genomes in the GTDB-Tk results are ignored. Every MAG must have a Bacteria or Archaea classification |
| `-o DIR` | Output directory. The table is written as `DIR/MAG_taxonomy.tsv` |

If no GTDB-Tk result is available, classify the MAGs first with `gtdbtk classify_wf`, then convert its summary files as above. GTDB-Tk is not installed with MAGverify; see [GTDB-Tk](https://github.com/ecogenomics/gtdbtk) for installation, reference data and options.

After steps 1 and 2, `$DB_OUT` contains both chimera inputs:

```text
$DB_OUT/GTDB_R232*            -> MAGverify chimera --database $DB_OUT/GTDB_R232
$DB_OUT/MAG_taxonomy.tsv      -> MAGverify chimera --mag-taxonomy $DB_OUT/MAG_taxonomy.tsv
```

## 3. Single-copy gene annotation (anvi'o)

`false-closure` and `chimera` use single-copy gene (SCG) hits on all assemblies. Build them once per project with anvi'o, using the same `-o` as the analysis modules:

```bash
MAGverify prep-scg \
  -a /path/to/assemblies \
  -o /path/to/MAGverify_output \
  -t 100
```

| Option | Input or output |
| --- | --- |
| `-a DIR` | Directory of the original assembler FASTA files (`.fa`, `.fna`, `.fasta`) |
| `-o DIR` | Analysis output directory. All assemblies are written to `DIR/00_prep/all_assemblies.fasta` and the SCG table to `DIR/00_prep/SCG/SCG_from_anvio.original_contig.tsv` |

The command builds an anvi'o contigs database and runs `anvi-run-hmms` with the Bacteria_71 and Archaea_76 collections. `false-closure` stops if the table is missing; `chimera` runs without SCG annotation and prints a warning. A completed table is reused when the assemblies have not changed.

## 4. Assembler circularity table

`false-closure` compares the circularity reported by each assembler with read evidence, and `circularity` uses it to separate high-confidence circular MAGs. Build the table from the assemblers' own outputs, giving one `-c` per assembler:

```bash
MAGverify prep-circularity \
  -c /path/to/hifiasm_meta/Read.p_ctg.gfa \
  -c /path/to/metaMDBG/contigs.fasta.gz \
  -c /path/to/metaflye/assembly_info.txt \
  -g /path/to/MAGs \
  -o /path/to/MAGverify_output
```

| Assembler | File for `-c` | Circularity field |
| --- | --- | --- |
| hifiasm-meta | `*.p_ctg.gfa` (or its contig FASTA) | Contig name ends in `c` (circular) or `l` (linear) |
| metaMDBG | `contigs.fasta.gz` | Header `circular=yes` / `circular=no` |
| metaFlye | `assembly_info.txt` | `circ.` column (`Y`/`N`) |

The format and assembler are detected from the file content, gzip-compressed files are accepted, and the files are read in parallel. Give the complete files: every contig of every MAG must be listed, otherwise the command stops and names the missing contigs. Outputs in `DIR/00_prep/`:

| File | Content |
| --- | --- |
| `assembler_circularity.MAG.tsv` | `MAG_ID`, `Assembler`, `Assembler_circularity` (`Y` when the MAG is one contig that its assembler reported as circular), `Number_of_contigs`, `Contig` |
| `assembler_circularity.contig.tsv` | `Contig`, `Assembler`, `Assembler_circularity`, `Source_file` for every contig in the given files |

