# Usage

MAGverify runs four preparation commands and ten analysis modules. Run the applicable preparation steps first, then `MAGverify end-to-end` for all QC modules or the commands below for individual steps.

Use the same `-o` for all commands to reuse intermediate results. All outputs are organized under the specified base directory:

```
BASE/
├── 00_prep/             # Preparation files
├── 01_mapping/          # Read mapping
├── 02_circularity/      # Circularity
├── 03_clipping/         # Clipping
├── ...                  # Other analysis modules
├── summary.MAG.tsv      # Integrated MAG-level summary
├── summary.contig.tsv   # Integrated contig-level summary
├── detail.MAG.tsv       # MAG summary with every module column
└── detail.contig.tsv    # Contig summary with every module column
```

| Step | Command | Preparation | Upstream modules |
| --- | --- | --- | --- |
| 1 | [`map`](#1-read-mapping-map) | — | — |
| 2 | [`circularity`](#2-circularity-circularity) | `prep-circularity` | `map` |
| 3 | [`clipping`](#3-clipping-clipping) | — | `map` |
| 4 | [`false-closure`](#4-false-closure-false-closure) | `prep-scg`, `prep-circularity` | `clipping` |
| 5 | [`phantom`](#5-phantom-sequences-phantom) | — | `map` |
| 6 | [`repeats`](#6-repeats-repeats) | — | `map` |
| 7 | [`false-duplication`](#7-false-duplication-false-duplication) | — | `clipping` |
| 8 | [`variants`](#8-assembly-variants-variants) | — | `map` |
| 9 | [`chimera`](#9-taxonomic-chimera-chimera) | `prep-chimera`, `prep-scg`, `prep-taxonomy` | `clipping` |
| 10 | [`summary`](#10-integrated-summary-summary) | — | the modules to summarize |

## Preparation

Preparation is not part of `end-to-end`. Its outputs are written to `/path/to/MAGverify_output/00_prep/`.

| Command | Output in `00_prep/` | Used by |
| --- | --- | --- |
| [`prep-chimera`](#1-gtdb-r232-database-prep-chimera) | `GTDB_R232*` | `chimera` |
| [`prep-scg`](#2-single-copy-gene-annotation-prep-scg) | `scg.contig.tsv`, `genes.*` | `false-closure`, `chimera`, `variants` |
| [`prep-circularity`](#3-assembler-circularity-prep-circularity) | `assembler_circularity.MAG.tsv`, `assembler_circularity.contig.tsv` | `circularity`, `false-closure` |
| [`prep-taxonomy`](#4-mag-taxonomy-prep-taxonomy) | `gtdbtk_taxonomy.MAG.tsv` | `chimera` |

### 1. GTDB R232 database (`prep-chimera`)

`prep-chimera` builds an MMseqs2 taxonomy database from three [GTDB R232](https://data.gtdb.ecogenomic.org/releases/release232/232.0/) files. The download is 123 GB and the build takes hours, so a [prebuilt database](#alternative-prebuilt-database) is offered at the end of this section.

**Download input files** into the same directory with `aria2c`:

```bash
GTDB_DIR=/path/to/GTDB_R232_downloads
BASE=https://data.gtdb.ecogenomic.org/releases/release232/232.0
mkdir -p "$GTDB_DIR"

# Protein sequences (~123 GB compressed)
aria2c -c -x 8 -s 8 -d "$GTDB_DIR" "$BASE/genomic_files_reps/gtdb_proteins_aa_reps_r232.tar.gz"
# Bacterial taxonomy
aria2c -c -d "$GTDB_DIR" "$BASE/bac120_taxonomy_r232.tsv"
# Archaeal taxonomy
aria2c -c -d "$GTDB_DIR" "$BASE/ar53_taxonomy_r232.tsv"
```

**Prerequisites:** the three GTDB R232 files downloaded above, given with `-d`

```bash
MAGverify prep-chimera -d "$GTDB_DIR" -o /path/to/MAGverify_output -t 100
```

**Output:** `GTDB_R232*`

#### Alternative: prebuilt database

A prebuilt database is available on [Hugging Face](https://huggingface.co/datasets/ATCC19718/MAGverify-resources), which skips both the 123 GB download and the build. Extract `GTDB_R232.tar.gz` (131 GB) into `$Base/00_prep/`, where `chimera` then finds it under the same `-o` without running `prep-chimera`:

```bash
Base=/path/to/MAGverify_output
mkdir -p "$Base/00_prep"

aria2c -c -x 8 -s 8 -k 8M \
  -d "$Base" -o GTDB_R232.tar.gz \
  "https://huggingface.co/datasets/ATCC19718/MAGverify-resources/resolve/main/GTDB_R232.tar.gz"

tar -xzf "$Base/GTDB_R232.tar.gz" -C "$Base/00_prep"
```

### 2. Single-copy gene annotation (`prep-scg`)

`prep-scg` predicts genes with Pyrodigal and identifies single-copy genes (SCGs) using anvi'o's `Bacteria_71` and `Archaea_76` collections. The gene predictions are also reused by `variants`.

**Prerequisites:** none. The assemblies are given with `-a`

```bash
MAGverify prep-scg -a /path/to/input_dataset/contigs -o /path/to/MAGverify_output -t 100
```

- `-a DIR`: Input directory of contig files (`.fa`, `.fasta`, or `.fna`), one per assembler. Contig IDs must be unique across assembly files.

**Output:** `scg.contig.tsv`, `genes.faa`, `genes.gff`, `genes.coordinates.tsv`

### 3. Assembler circularity (`prep-circularity`)

`prep-circularity` collects assembler-reported circularity calls for `circularity` and `false-closure`. Provide one `-c` file per assembler: `*.p_ctg.gfa` or `*.p_ctg.noseq.gfa` (faster) for hifiasm-meta, `contigs.fasta.gz` for metaMDBG, and `assembly_info.txt` for metaFlye.

**Prerequisites:** none. The assembler output files are given with `-c` and the MAGs with `-g`

```bash
MAGverify prep-circularity \
  -c /path/to/input_dataset/assembler_circularity/Read.p_ctg.noseq.gfa \
  -c /path/to/input_dataset/assembler_circularity/contigs.fasta.gz \
  -c /path/to/input_dataset/assembler_circularity/assembly_info.txt \
  -g /path/to/input_dataset/MAGs \
  -o /path/to/MAGverify_output
```

- `-g DIR`: Input directory of MAG files (one file per MAG). MAG contig IDs must match those in the assembly files.

**Output:** `assembler_circularity.MAG.tsv` and `assembler_circularity.contig.tsv`

### 4. MAG taxonomy (`prep-taxonomy`)

`prep-taxonomy` converts [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) `classify_wf` results to the MAG taxonomy table that `chimera` needs for its SCG redundancy annotation.

**Prerequisites:** none. The GTDB-Tk results are given with `-d` and the MAGs with `-g`

```bash
MAGverify prep-taxonomy \
  -d /path/to/input_dataset/gtdbtk \
  -g /path/to/input_dataset/MAGs \
  -o /path/to/MAGverify_output
```

- `-d DIR`: GTDB-Tk output directory with `gtdbtk.bac120.summary.tsv` and/or `gtdbtk.ar53.summary.tsv`.

**Output:** `gtdbtk_taxonomy.MAG.tsv`

## Analysis modules

Each module reads its prerequisites' outputs automatically from the same `-o`.

### 1. Read mapping (`map`)

Maps HiFi reads separately to each assembler's complete assembly with minimap2, including unbinned contigs. Saves BAMs of all mapped alignments, primary alignments, and primary plus supplementary alignments for later modules. A mapping manifest (`01_mapping/.cache/mapping_manifest.tsv`) records the reads, assemblies and BAM paths so later modules do not need `-a` or `-r`.

**Prerequisites:** none. The assemblies are given with `-a` and the HiFi reads with `-r`

```bash
MAGverify map -a /path/to/input_dataset/contigs -r /path/to/input_dataset/reads -o /path/to/MAGverify_output -t 100 -m 900
```

**Output:** `01_mapping/mapping.assembler.tsv` (alignment counts per assembler)

### 2. Circularity (`circularity`)

Tests single-contig MAGs for read-supported circular topology. The contig is rotated to position its END|START junction in the middle. Circularity is read-supported when at least two HiFi reads span the junction with at least 2,000 bp aligned on each side.

**Prerequisites:** `map`, `prep-circularity`

```bash
MAGverify circularity -g /path/to/input_dataset/MAGs -o /path/to/MAGverify_output -t 100 -m 900
```

**Output:** `02_circularity/circularity.MAG.tsv`, `circularity.contig.tsv`

**High-confidence circularity** needs two independent lines of evidence to agree: `High_confidence_circularity` is `Y` only when the assembler reports the contig circular **and** the reads span the junction. An assembler claim with no junction-spanning reads stays `Read_supported_circularity` `N`, and the call is `NA` without the assembler circularity table from `prep-circularity`.

The module tables report `Assembler_reported_circularity` and `Read_supported_circular_topology` separately, and read-supported contigs the assembler did not report circular appear as `Additional_circular_candidate`. `summary.MAG.tsv` and `summary.contig.tsv` carry all three calls.

### 3. Clipping (`clipping`)

Finds candidate positions using `anvi-script-find-misassemblies` (default clipping ratio > 0.8, `--min-dist-to-end 0`) and validates them against HiFi alignments. Positions within 100 bp of a contig end are `TERMINAL`; all others are `INTERNAL`.

**Prerequisites:** `map`

```bash
MAGverify clipping -g /path/to/input_dataset/MAGs -o /path/to/MAGverify_output -t 100
```

**Output:** `03_clipping/clipping.site.tsv`, `clipping.contig.tsv`, `clipping.MAG.tsv`, `clipping.assembler.tsv`

The module reports three calls. `Clipping` is the read-validated call: `Y` when at least two distinct HiFi reads are clipped at the exact position and no read spans that position with aligned bases on both sides, and `NA` at the first or last base of a contig. `Stringent_candidate` applies the anvi'o thresholds alone, without read validation: `Y` when coverage is ≥ 10 and the clipping ratio is 1.0.

**High-confidence clipping** is `Stringent_clipping`, which needs both lines of evidence: it is `Y` only when `Clipping` and `Stringent_candidate` are both `Y`, following [Trigodet et al., 2026](https://doi.org/10.1038/s41587-025-02971-8). The thresholds alone are not sufficient, so a position with ≥ 10× coverage and a clipping ratio of 1.0 that is not validated against the reads is `Stringent_candidate` `Y` but `Stringent_clipping` `N`. `summary.MAG.tsv` and `summary.contig.tsv` carry this call as `High_confidence_clipping`, which is `NA` when clipping was not assessable.

### 4. False closure (`false-closure`)

Examines contigs reported circular by their assembler. A validated clipping site with at least two split reads on **each** side connecting the contig to another contig indicates a false-closure breakpoint. The module also reports a small-genome proxy for circular contigs < 500,000 bp with ≥ 3 ribosomal protein families (`prep-scg`). This screening field is distinct from the final `Premature_circularization` call.

**Prerequisites:** `clipping`, `prep-scg`, `prep-circularity`

```bash
MAGverify false-closure -g /path/to/input_dataset/MAGs -o /path/to/MAGverify_output -t 100
```

**Output:** `04_false_closure/false_closure.contig.tsv`, `false_closure.MAG.tsv`

### 5. Phantom sequences (`phantom`)

Finds sequence stretches longer than 500 bp with zero HiFi read depth, no BLASTN hit to the raw reads and no ambiguous bases. At the MAG level, Meryl also measures the fraction of distinct 21-mers in phantom regions that are absent from the reads.

**Prerequisites:** `map`

```bash
MAGverify phantom -g /path/to/input_dataset/MAGs -o /path/to/MAGverify_output -t 100
```

**Output:** `05_phantom/phantom.region.tsv`, `phantom.contig.tsv`, `phantom.MAG.tsv`

### 6. Repeats (`repeats`)

Finds BLASTN repeats (longer than 200 bp, ≥ 80% identity). For single-contig MAGs, the contig is scanned against itself, and the MAG call follows that contig. For multi-contig MAGs, each contig is scanned against itself and the contigs are additionally compared with each other.

A contig is flagged when repeats within it cover ≥ 0.70 of its bases. A MAG is flagged when repeats within and between its contigs together — overlapping bases counted once — cover ≥ 0.70 of its bases.

**Prerequisites:** `map`

```bash
MAGverify repeats -g /path/to/input_dataset/MAGs -o /path/to/MAGverify_output -t 100 -m 900
```

**Output:** `06_repeats/repeats.contig.tsv`, `repeats.MAG.tsv`, `repeats.region.tsv`, `repeats.HSP.tsv`

### 7. False duplication (`false-duplication`)

Identifies duplicated copies within a MAG from a minimap2 self-alignment of each assembly (`-x asm5 -DP`, ≥ 1,000 bp, ≥ 95% identity) and checks each copy boundary against HiFi alignments. A duplicated pair is flagged when one copy is spanned at both boundaries and the other at neither. For tandem copies, both outer boundaries must be spanned but the intervening junction must not be. Clipping evidence at the copy boundaries is also reported.

**Prerequisites:** `clipping`

```bash
MAGverify false-duplication -g /path/to/input_dataset/MAGs -o /path/to/MAGverify_output -t 100
```

**Output:** `08_false_duplication/false_duplication.pair.tsv`, `false_duplication.contig.tsv`, `false_duplication.MAG.tsv`

### 8. Assembly variants (`variants`)

Profiles the HiFi reads against each assembly with anvi'o and identifies assembly alleles that disagree with the reads. A **minority** variant occurs when another allele or an indel is supported by > 50% of reads. An **unsupported** variant occurs when no read supports the assembly allele.

**Prerequisites:** `map`

```bash
MAGverify variants -g /path/to/input_dataset/MAGs -o /path/to/MAGverify_output -t 100 -m 900
```

**Output:** `09_assembly_variants/variants.site.tsv`, `variants.contig.tsv`, `variants.MAG.tsv`

### 9. Taxonomic chimera (`chimera`)

Classifies genes against GTDB R232 with MMseqs2, then finds neighbouring gene blocks whose dominant taxa differ at domain, phylum, class or order.

**Prerequisites:** `clipping`, `prep-chimera`, `prep-scg`, `prep-taxonomy`

To reuse a GTDB database built under another output directory, pass its prefix with `--database PREFIX` — the MMseqs2 prefix itself, not its directory or its `.dbtype` file. With an existing gene taxonomy table, `--taxonomy TSV` runs this module without any database.

```bash
MAGverify chimera -g /path/to/input_dataset/MAGs -o /path/to/MAGverify_output -t 100
```

**Output:** `07_chimera/chimera.region.tsv`, `chimera.contig.tsv`, `chimera.MAG.tsv`

A transition is called `Taxonomic_chimera` when a clipping site with `Clipping` `Y` lies inside it, and `High_confidence_taxonomic_chimera` when that same site is also `Stringent_clipping` `Y`, the high-confidence clipping call. `summary.MAG.tsv` and `summary.contig.tsv` carry the high-confidence call as `Chimera`. The HiFi structural check across the boundary is reported separately as `Structural_chimera`.

### 10. Integrated summary (`summary`)

Combines completed module calls into MAG- and contig-level summaries and checks MAG IDs, contig membership, counts and consistency across tables. Validation failures return a nonzero exit status. For a partial run, select completed modules with `--include clipping,chimera`. The summary step does not rerun analyses.

**Prerequisites:** the modules to be summarized

```bash
MAGverify summary -g /path/to/input_dataset/MAGs -o /path/to/MAGverify_output
```

**Output:** `summary.MAG.tsv` and `summary.contig.tsv` (concise `Y` / `N` / `NA` flags); `detail.MAG.tsv` and `detail.contig.tsv` (full diagnostic tables)
