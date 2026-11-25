# Usage

Run the applicable preparation steps, followed by `end-to-end` or individual QC modules. Use the same output directory (`-o`) throughout.

## 1. Input configuration

The examples use the directory structure of the [test dataset](../README.md#2-test-dataset).

```bash
Data=/path/to/input_dataset
Base=/path/to/MAGverify_output
Threads=100
Memory=900

mkdir -p "$Base/00_prep"
```

Adjust the example resource values to your system.

| Option | Description |
| --- | --- |
| `-a DIR` | Assembly FASTA directory, one file per assembler |
| `-g DIR` | MAG FASTA directory, one file per MAG |
| `-r DIR` | HiFi FASTQ directory |
| `-o DIR` | Shared output directory |
| `-t N` | Threads (default: 80% of available CPUs) |
| `-m N` | Memory budget in GB for `map`, `circularity`, `repeats` and `variants` (default: 80% of available memory) |

Assemblies and MAGs accept `.fa`, `.fasta` and `.fna`. Use unique contig IDs across assemblies and retain the original IDs in MAGs. FASTA filenames must be unique and must not include whitespace. Reads accept `.fastq`, `.fq`, `.fastq.gz` and `.fq.gz`. Supply the complete HiFi read set used to generate all input assemblies. Missing or unrelated reads can distort QC calls.

`PARALLEL_JOBS` controls concurrent assembler jobs. Mapping, circularity and anvi'o profiling in `variants` process assemblers sequentially. Clipping, phantom and false duplication can distribute threads across assemblers.

## 2. Preparation

Preparation is not part of `end-to-end`. Outputs are stored in `00_prep/`.

| Command | Main output | Used by |
| --- | --- | --- |
| `prep-chimera` | `GTDB_R232*` | `chimera` |
| `prep-scg` | `scg.contig.tsv`, `genes.*` | `false-closure`, `chimera`, `variants` |
| `prep-circularity` | `assembler_circularity.*.tsv` | `circularity`, `false-closure` |
| `prep-taxonomy` | `gtdbtk_taxonomy.MAG.tsv` | `chimera` |

### 2.1. GTDB R232 database

Choose one option to prepare the MMseqs2 taxonomy database.

#### Option A. Prebuilt database (recommended)

Download the prebuilt database (131 GB) from [Hugging Face](https://huggingface.co/datasets/ATCC19718/MAGverify-resources) to avoid downloading the original GTDB protein archive and building the database locally.

```bash
aria2c -c -x 8 -s 8 -k 8M \
  -d "$Base" -o GTDB_R232.tar.gz \
  "https://huggingface.co/datasets/ATCC19718/MAGverify-resources/resolve/main/GTDB_R232.tar.gz"

tar -xzf "$Base/GTDB_R232.tar.gz" -C "$Base/00_prep"
```

The database prefix should be `$Base/00_prep/GTDB_R232`. Keep all database companion files together. No `prep-chimera` run is needed.

#### Option B. Local database build

```bash
GTDB_DIR=/path/to/GTDB_R232_downloads
GTDB_URL=https://data.gtdb.ecogenomic.org/releases/release232/232.0
mkdir -p "$GTDB_DIR"

aria2c -c -x 8 -s 8 -d "$GTDB_DIR" "$GTDB_URL/genomic_files_reps/gtdb_proteins_aa_reps_r232.tar.gz"
aria2c -c -d "$GTDB_DIR" "$GTDB_URL/bac120_taxonomy_r232.tsv"
aria2c -c -d "$GTDB_DIR" "$GTDB_URL/ar53_taxonomy_r232.tsv"

MAGverify prep-chimera -d "$GTDB_DIR" -o "$Base" -t "$Threads"
```

Allow space for the compressed archive (~123 GB), database and temporary build files in `00_prep/.cache/prep_chimera_work/`.

### 2.2. Single-copy gene annotation

```bash
MAGverify prep-scg -a "$Data/contigs" -o "$Base" -t "$Threads"
```

### 2.3. Assembler circularity

```bash
MAGverify prep-circularity \
  -c "$Data/assembler_circularity/Read.p_ctg.noseq.gfa" \
  -c "$Data/assembler_circularity/contigs.fasta.gz" \
  -c "$Data/assembler_circularity/assembly_info.txt" \
  -g "$Data/MAGs" -o "$Base"
```

Supply one `-c` per assembler used, with complete contig information.

| Assembler | Input | Circularity field |
| --- | --- | --- |
| hifiasm-meta | `*.p_ctg.gfa` or `*.p_ctg.noseq.gfa` | Contig suffix `c` or `l` |
| metaMDBG | `contigs.fasta.gz` | `circular=yes` or `circular=no` |
| metaFlye | `assembly_info.txt` | `circ.` column |

### 2.4. MAG taxonomy

```bash
MAGverify prep-taxonomy -d "$Data/gtdbtk" -g "$Data/MAGs" -o "$Base"
```

Supply GTDB-Tk `classify_wf` results (`gtdbtk.bac120.summary.tsv` and/or `gtdbtk.ar53.summary.tsv`). MAG IDs must match FASTA basenames without extensions. GTDB-Tk is installed separately.

## 3. End-to-end analysis

### Run all modules

```bash
MAGverify end-to-end \
  -a "$Data/contigs" -r "$Data/reads" -g "$Data/MAGs" \
  -o "$Base" -t "$Threads" -m "$Memory"
```

### Skip selected modules

```bash
MAGverify end-to-end \
  -a "$Data/contigs" -r "$Data/reads" -g "$Data/MAGs" \
  -o "$Base" -t "$Threads" -m "$Memory" \
  --skip chimera,variants
```

Mapping cannot be skipped. Retain dependencies for all selected modules (e.g. clipping for false closure, false duplication and chimera).

### Reuse a database

Add `--database /path/to/GTDB_R232` to `end-to-end` or `chimera` to use an existing database. Specify the MMseqs2 prefix, not its directory or `.dbtype` file.

### Resume analysis

Rerun the original command with unchanged inputs. Preserve `.cache/` directories and intermediate files to reuse completed work.

## 4. Individual modules

Run commands in dependency order. Module outputs are retrieved from the shared output directory.

| Module | Preparation | Upstream modules |
| --- | --- | --- |
| `map` | — | — |
| `circularity` | `prep-circularity` | `map` |
| `clipping` | — | `map` |
| `false-closure` | `prep-scg`, `prep-circularity` | `clipping` |
| `phantom` | — | `map` |
| `repeats` | — | `map` |
| `false-duplication` | — | `clipping` |
| `variants` | — | `map` |
| `chimera` | GTDB database, `prep-scg`, `prep-taxonomy` | `clipping` |
| `summary` | — | Selected QC modules |

### 4.1. Read mapping

```bash
MAGverify map -a "$Data/contigs" -r "$Data/reads" -o "$Base" -t "$Threads" -m "$Memory"
```

### 4.2. Circularity

```bash
MAGverify circularity -g "$Data/MAGs" -o "$Base" -t "$Threads" -m "$Memory"
```

### 4.3. Clipping

```bash
MAGverify clipping -g "$Data/MAGs" -o "$Base" -t "$Threads"
```

### 4.4. False closure

```bash
MAGverify false-closure -g "$Data/MAGs" -o "$Base" -t "$Threads"
```

### 4.5. Phantom sequences

```bash
MAGverify phantom -g "$Data/MAGs" -o "$Base" -t "$Threads"
```

### 4.6. Repeats

```bash
MAGverify repeats -g "$Data/MAGs" -o "$Base" -t "$Threads" -m "$Memory"
```

### 4.7. False duplication

```bash
MAGverify false-duplication -g "$Data/MAGs" -o "$Base" -t "$Threads"
```

### 4.8. Assembly variants

```bash
MAGverify variants -g "$Data/MAGs" -o "$Base" -t "$Threads" -m "$Memory"
```

### 4.9. Taxonomic chimera

```bash
MAGverify chimera -g "$Data/MAGs" -o "$Base" -t "$Threads"
```

An existing gene taxonomy table can be supplied with `--taxonomy TSV` instead of building or downloading the database for individual chimera analysis.

### 4.10. Integrated QC summary

```bash
MAGverify summary -g "$Data/MAGs" -o "$Base"
```

For a partial analysis:

```bash
MAGverify summary -g "$Data/MAGs" -o "$Base" --include clipping,chimera
```

The summary reads completed results without rerunning QC. Use `MAGverify list` and `MAGverify <command> -h` for command-specific options.

## 5. Output files

| Module | Directory |
| --- | --- |
| Preparation | `00_prep/` |
| Mapping | `01_mapping/` |
| Circularity | `02_circularity/` |
| Clipping | `03_clipping/` |
| False closure | `04_false_closure/` |
| Phantom | `05_phantom/` |
| Repeats | `06_repeats/` |
| False duplication | `08_false_duplication/` |
| Variants | `09_assembly_variants/` |
| Chimera | `07_chimera/` |
| Summary | Output root |

Directory numbers are retained independently of execution order.

The output root includes `summary.MAG.tsv`, `summary.contig.tsv`, `detail.MAG.tsv` and `detail.contig.tsv`.

See [Analysis modules](MODULES.md) for detection criteria and [Result tables](OUTPUTS.md) for column definitions.
