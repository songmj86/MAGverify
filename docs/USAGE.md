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

Downstream modules automatically reuse mapping BAMs, clipping sites and preparation tables from this directory without rerunning earlier steps.

| Step | Command | Needs first |
| --- | --- | --- |
| Preparation | [`prep-chimera`](#1-gtdb-r232-database), [`prep-scg`](#2-single-copy-gene-annotation), [`prep-circularity`](#3-assembler-circularity), [`prep-taxonomy`](#4-mag-taxonomy) | Run applicable commands separately before analysis |
| 1 | [`map`](#1-read-mapping) | |
| 2 | [`circularity`](#2-circularity) | `map`, `prep-circularity` |
| 3 | [`clipping`](#3-clipping) | `map` |
| 4 | [`false-closure`](#4-false-closure) | `clipping`, `prep-scg`, `prep-circularity` |
| 5 | [`phantom`](#5-phantom-sequences) | `map` |
| 6 | [`repeats`](#6-repeats) | `map` |
| 7 | [`false-duplication`](#7-false-duplication) | `clipping` |
| 8 | [`variants`](#8-assembly-variants) | `map` |
| 9 | [`chimera`](#9-taxonomic-chimera) | `clipping`, `prep-chimera`, `prep-scg`, `prep-taxonomy` |
| 10 | [`summary`](#10-integrated-summary) | Modules to be summarized |

## Preparation

Run the applicable preparation commands **before** `end-to-end`, which does not include preparation. Keep the same `-o` for preparation and analysis.

**Every preparation output is written to `/path/to/MAGverify_output/00_prep/`.**

### 1. GTDB R232 database

`prep-chimera` builds an MMseqs2 taxonomy database from three [GTDB R232](https://data.gtdb.ecogenomic.org/releases/release232/232.0/) files.

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

**Command**

```bash
MAGverify prep-chimera -d "$GTDB_DIR" -o /path/to/MAGverify_output -t 100
```

**Output:** `GTDB_R232*`.

---

### 2. Single-copy gene annotation

`prep-scg` predicts genes with Pyrodigal and identifies single-copy genes (SCGs) using anvi'o's `Bacteria_71` and `Archaea_76` collections.

**Command**

```bash
MAGverify prep-scg -a /path/to/assemblies -o /path/to/MAGverify_output -t 100
```

- `-a DIR`: Input directory of contig files (`.fa`, `.fasta`, or `.fna`), one per assembler. Contig IDs must be unique across assembly files.

The SCG table is used by `false-closure` and for SCG annotation in `chimera`. The gene predictions are reused by `chimera` and `variants`.

**Output:** `scg.contig.tsv`, `genes.faa`, `genes.gff`, `genes.coordinates.tsv`.

---

### 3. Assembler circularity

`prep-circularity` collects assembler-reported circularity calls for `circularity` and `false-closure`. Provide one `-c` file per assembler: `*.p_ctg.gfa` or `*.p_ctg.noseq.gfa` (faster) for hifiasm-meta, `contigs.fasta.gz` for metaMDBG, and `assembly_info.txt` for metaFlye.

**Command**

```bash
MAGverify prep-circularity \
  -c /path/to/hifiasm_meta/Read.p_ctg.gfa \
  -c /path/to/metaMDBG/contigs.fasta.gz \
  -c /path/to/metaFlye/assembly_info.txt \
  -g /path/to/MAGs \
  -o /path/to/MAGverify_output
```

- `-g DIR`: Input directory of MAG files (one file per MAG). MAG contig IDs must match those in the assembly files.

**Output:** `assembler_circularity.MAG.tsv` and `assembler_circularity.contig.tsv`.

---

### 4. MAG taxonomy

`prep-taxonomy` converts [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) `classify_wf` results to the MAG taxonomy table that `chimera` needs for its SCG redundancy annotation.

**Command**

```bash
MAGverify prep-taxonomy \
  -d /path/to/GTDBTk_results \
  -g /path/to/MAGs \
  -o /path/to/MAGverify_output
```

- `-d DIR`: GTDB-Tk output directory with `gtdbtk.bac120.summary.tsv` and/or `gtdbtk.ar53.summary.tsv`.
- `-g DIR`: Input directory of MAG files (one file per MAG). MAG contig IDs must match those in the assembly files.

**Output:** `gtdbtk_taxonomy.MAG.tsv`.

## Analysis modules

### 1. Read mapping

Maps HiFi reads separately to each assembler's complete assembly with minimap2, including unbinned contigs. Saves BAMs of all mapped alignments, primary alignments, and primary plus supplementary alignments for later modules. A mapping manifest (`01_mapping/.cache/mapping_manifest.tsv`) records the reads, assemblies and BAM paths so later modules do not need `-a` or `-r`.

**Command**

```bash
MAGverify map -a /path/to/assemblies -r /path/to/hifi_reads -o /path/to/MAGverify_output -t 100 -m 900
```

**Output:** `01_mapping/mapping.assembler.tsv` (alignment counts per assembler).

---

### 2. Circularity

Tests single-contig MAGs for read-supported circular topology. The contig is rotated to position its END|START junction in the middle. Circularity is read-supported when at least two HiFi reads span the junction with at least 2,000 bp aligned on each side. The `prep-circularity` assembler call is also needed.

**High-confidence circularity** needs two independent lines of evidence to agree: `High_confidence_circularity` is `Y` only when the assembler reports the contig circular **and** the reads span the junction. An assembler claim with no junction-spanning reads stays `Read_supported_circularity` `N`, and the call is `NA` without the assembler circularity table from `prep-circularity`.

The module tables report `Assembler_reported_circularity` and `Read_supported_circular_topology` separately, and read-supported contigs the assembler did not report circular appear as `Additional_circular_candidate`. `summary.MAG.tsv` and `summary.contig.tsv` carry all three calls.

**Command**

```bash
MAGverify circularity -g /path/to/MAGs -o /path/to/MAGverify_output -t 100 -m 900
```

**Output:** `02_circularity/circularity.MAG.tsv`, `circularity.contig.tsv`.

---

### 3. Clipping

Finds candidate positions using `anvi-script-find-misassemblies` (default clipping ratio > 0.8, `--min-dist-to-end 0`) and validates them against HiFi alignments. Positions within 100 bp of a contig end are `TERMINAL`; all others are `INTERNAL`.

**Clipping calls**

| Call | Criterion |
| --- | --- |
| `Clipping` | ≥ 2 distinct HiFi reads clipped at the exact position **and** no read spanning that position with aligned bases on both sides; `NA` at the first or last base of a contig |
| `Stringent_candidate` | anvi'o coverage ≥ 10 **and** anvi'o clipping ratio = 1.0 (thresholds only, without read validation) |
| `Stringent_clipping` | **High-confidence clipping:** `Clipping` is `Y` **and** `Stringent_candidate` is `Y`, following [Trigodet et al., 2026](https://doi.org/10.1038/s41587-025-02971-8) |

**High-confidence clipping** is `Stringent_clipping`, and it needs both lines of evidence: the thresholds alone are not sufficient, so a position with ≥ 10× coverage and a clipping ratio of 1.0 that is not validated against the reads is `Stringent_candidate` `Y` but `Stringent_clipping` `N`. `summary.MAG.tsv` and `summary.contig.tsv` carry this call as `High_confidence_clipping`, which is `NA` when clipping was not assessable. `chimera` builds on it: a taxonomic transition is a high-confidence chimera only when the clipping site it overlaps is itself high confidence.

Validated sites are reused by `false-closure`, `false-duplication`, `chimera` and `summary`.

**Command**

```bash
MAGverify clipping -g /path/to/MAGs -o /path/to/MAGverify_output -t 100
```

**Output:** `03_clipping/clipping.site.tsv`, `clipping.contig.tsv`, `clipping.MAG.tsv`, `clipping.assembler.tsv`.

---

### 4. False closure

Examines contigs reported circular by their assembler. A validated clipping site with at least two split reads on **each** side connecting the contig to another contig indicates a false-closure breakpoint. The module also reports a small-genome proxy for circular contigs < 500,000 bp with ≥ 3 ribosomal protein families (`prep-scg`). This screening field is distinct from the final `Premature_circularization` call.

Requires `prep-scg` and `prep-circularity` under the same output directory. To provide an existing assembler circularity table, use `--assembler-circularity TSV`.

**Command**

```bash
MAGverify false-closure -g /path/to/MAGs -o /path/to/MAGverify_output -t 100
```

**Output:** `04_false_closure/false_closure.contig.tsv`, `false_closure.MAG.tsv`.

---

### 5. Phantom sequences

Finds sequence stretches of at least 501 bp with zero HiFi read depth, no BLASTN hit to the raw reads and no ambiguous bases. At the MAG level, Meryl also measures the fraction of distinct 21-mers in phantom regions that are absent from the reads.

**Command**

```bash
MAGverify phantom -g /path/to/MAGs -o /path/to/MAGverify_output -t 100
```

**Output:** `05_phantom/phantom.region.tsv`, `phantom.contig.tsv`, `phantom.MAG.tsv`.

---

### 6. Repeats

Finds BLASTN repeats (≥ 200 bp, ≥ 80% identity) within contigs (`INTRA`) and between contigs of the same MAG (`INTER`). Flags a contig when its INTRA repeat fraction is ≥ 0.70 or a MAG when its total repeat fraction is ≥ 0.70. Repeat detection and the false-duplication check are separate assessments.

**Command**

```bash
MAGverify repeats -g /path/to/MAGs -o /path/to/MAGverify_output -t 100 -m 900
```

**Output:** `06_repeats/repeats.contig.tsv`, `repeats.MAG.tsv`, `repeats.region.tsv`, `repeats.HSP.tsv`.

---

### 7. False duplication

Identifies duplicated copies within a MAG from a minimap2 self-alignment of each assembly (`-x asm5 -DP`, ≥ 1,000 bp, ≥ 95% identity) and checks each copy boundary against HiFi alignments. A duplicated pair is flagged when one copy is spanned at both boundaries and the other at neither. For tandem copies, both outer boundaries must be spanned but the intervening junction must not be. Clipping evidence at the copy boundaries is also reported.

**Command**

```bash
MAGverify false-duplication -g /path/to/MAGs -o /path/to/MAGverify_output -t 100
```

**Output:** `08_false_duplication/false_duplication.pair.tsv`, `false_duplication.contig.tsv`, `false_duplication.MAG.tsv`.

---

### 8. Assembly variants

Profiles the HiFi reads against each assembly with anvi'o and identifies assembly alleles that disagree with the reads. A **minority** variant occurs when another allele or an indel is supported by > 50% of reads. An **unsupported** variant occurs when no read supports the assembly allele. This module can need substantial RAM; use `-m N` to set a memory limit in GB.

**Command**

```bash
MAGverify variants -g /path/to/MAGs -o /path/to/MAGverify_output -t 100 -m 900
```

**Output:** `09_assembly_variants/variants.site.tsv`, `variants.contig.tsv`, `variants.MAG.tsv`.

---

### 9. Taxonomic chimera

Classifies genes against GTDB R232 with MMseqs2, then finds neighbouring gene blocks whose dominant taxa differ at domain, phylum, class or order. A transition is called a taxonomic chimera when it overlaps a validated clipping site; high-confidence calls additionally need a high-confidence clipping site. HiFi structural evidence across the boundary is reported separately.

The database from `prep-chimera` is reused automatically under the same `-o`. Specify `--database PREFIX` to use one elsewhere or `--taxonomy TSV` to use an existing gene taxonomy table. SCG annotation uses `prep-scg` and MAG taxonomy uses `prep-taxonomy`.

**Command**

```bash
MAGverify chimera -g /path/to/MAGs -o /path/to/MAGverify_output -t 100
```

**Output:** `07_chimera/chimera.region.tsv`, `chimera.contig.tsv`, `chimera.MAG.tsv`.

---

### 10. Integrated summary

Combines completed module calls into MAG- and contig-level summaries and checks MAG IDs, contig membership, counts and consistency across tables. Validation failures return a nonzero exit status. For a partial run, select completed modules with `--include clipping,chimera`. The summary step does not rerun analyses.

**Command**

```bash
MAGverify summary -g /path/to/MAGs -o /path/to/MAGverify_output
```

**Output:** `summary.MAG.tsv` and `summary.contig.tsv` (concise `Y` / `N` / `NA` flags); `detail.MAG.tsv` and `detail.contig.tsv` (full diagnostic tables).

For all output column definitions, see [Result tables](OUTPUTS.md).
