# Modules

MAGverify runs read-based QC analyses in a defined order and combines their findings into MAG- and contig-level summaries. After completing [input preparation](INPUT_PREPARATION.md), run `MAGverify end-to-end` to execute all modules or use the commands below to run individual steps. Use the **same `-o`** throughout so modules can reuse mapping results, clipping tables and preparation files.

| Step | Command | Needs first |
| --- | --- | --- |
| Preparation | `prep-chimera`, `prep-scg`, `prep-circularity`, `prep-taxonomy` | Run separately before analysis |
| 1 | [`map`](#1-read-mapping) | — |
| 2 | [`circularity`](#2-circularity) | `map`, `prep-circularity` |
| 3 | [`clipping`](#3-clipping) | `map` |
| 4 | [`false-closure`](#4-false-closure) | `clipping`, `prep-scg`, `prep-circularity` |
| 5 | [`phantom`](#5-phantom-sequences) | `map` |
| 6 | [`repeats`](#6-repeats) | `map` |
| 7 | [`false-duplication`](#7-false-duplication) | `clipping` |
| 8 | [`variants`](#8-assembly-variants) | `map` |
| 9 | [`chimera`](#9-taxonomic-chimera) | `clipping`, `prep-chimera`; `prep-scg` and `prep-taxonomy` for SCG redundancy |
| 10 | [`summary`](#10-integrated-qc-summary) | Modules to be summarized |

The examples use these directory variables:

```bash
In_Contig=/path/to/assemblies
In_Read=/path/to/hifi_reads
In_MAG=/path/to/MAGs
Out=/path/to/MAGverify_output
```

**Resources:** `-t N` sets the thread count (default: all available CPUs, capped at 100). `map`, `circularity` and the anvi'o profiling in `variants` process one assembler at a time using all allocated threads. `clipping`, `phantom` and `false-duplication` process assemblers in parallel, dividing threads among jobs. Set `PARALLEL_JOBS=N` to limit concurrent assemblers. `-m N` sets a memory limit in GB for `map`, `circularity`, `repeats` and `variants`.

**Progress and resumption:** Modules print progress, elapsed time and estimated remaining time every five minutes. Set `PROGRESS_INTERVAL=SECONDS` to change the interval. Estimates use earlier run times saved in `$Out/.cache/run_times.tsv`; without run history, only elapsed time is shown. Read mapping and the BLAST steps of `repeats` also show progress on first runs. Completed results can be reused when inputs are unchanged.

## Preparation

Preparation is not part of `end-to-end`. Run the applicable steps first using the same `-o`:

```bash
MAGverify prep-circularity -c FILE [-c FILE ...] -g "$In_MAG" -o "$Out"
MAGverify prep-taxonomy    -d /path/to/GTDBTk_results -g "$In_MAG" -o "$Out"
MAGverify prep-scg         -a "$In_Contig" -o "$Out" -t 100
MAGverify prep-chimera     -d /path/to/GTDB_R232_downloads -o "$Out" -t 100
```

| Command | Output in `00_prep/` | Used by |
| --- | --- | --- |
| `prep-circularity` | `assembler_circularity.MAG.tsv`, `assembler_circularity.contig.tsv` | `circularity`, `false-closure` |
| `prep-taxonomy` | `gtdbtk_taxonomy.MAG.tsv` | `chimera` (optional) |
| `prep-scg` | `scg.contig.tsv`, `genes.*` | `false-closure`, `chimera` |
| `prep-chimera` | `GTDB_R232*` | `chimera` |

See [Input preparation](INPUT_PREPARATION.md) for the input formats, downloads and detailed options.

## 1. Read mapping

```bash
MAGverify map -a "$In_Contig" -r "$In_Read" -o "$Out" -t 100 -m 900
```

Maps HiFi reads separately to each assembler's complete assembly with minimap2, including unbinned contigs. Saves BAMs of all mapped alignments, primary alignments, and primary plus supplementary alignments for later modules. A mapping manifest (`01_mapping/.cache/mapping_manifest.tsv`) records the reads, assemblies and BAM paths so later modules do not need `-a` or `-r`.

**Output:** `01_mapping/mapping.assembler.tsv` (alignment counts per assembler).

## 2. Circularity

```bash
MAGverify circularity -g "$In_MAG" -o "$Out" -t 100 -m 900
```

Tests single-contig MAGs for read-supported circular topology. The contig is rotated to position its END|START junction in the middle. Circularity is read-supported when at least two HiFi reads span the junction with at least 2,000 bp aligned on each side. The `prep-circularity` assembler call is also needed.

The integrated MAG summary reports `Assembler_reported_circularity`, `Read_supported_circularity` and `High_confidence_circular` separately. High-confidence circularity is `Y` only when both the assembler and reads support the topology. Read-supported cases not reported circular by the assembler appear as `Additional_circular_candidate` in the module table.

**Output:** `02_circularity/circularity.MAG.tsv`, `circularity.contig.tsv`.

## 3. Clipping

```bash
MAGverify clipping -g "$In_MAG" -o "$Out" -t 100
```

Finds candidate positions using `anvi-script-find-misassemblies` (default clipping ratio > 0.8) and validates them against HiFi alignments. `Clipping` is `Y` when at least two distinct reads are clipped at the exact position and no read spans it. Stringent candidates additionally have anvi'o coverage ≥ 10 and clipping ratio = 1.0. Positions within 100 bp of a contig end are `TERMINAL`; all others are `INTERNAL`.

Validated sites are reused by `false-closure`, `false-duplication`, `chimera` and `summary`.

**Output:** `03_clipping/clipping.site.tsv`, `clipping.contig.tsv`, `clipping.MAG.tsv`, `clipping.assembler.tsv`.

## 4. False closure

```bash
MAGverify false-closure -g "$In_MAG" -o "$Out" -t 100
```

Examines contigs reported circular by their assembler. A validated clipping site with at least two split reads on **each** side connecting the contig to another contig indicates a false-closure breakpoint. The module also reports a small-genome proxy for circular contigs < 500,000 bp with ≥ 3 ribosomal protein families (`prep-scg`). This screening field is distinct from the final `Premature_circularization` call.

Requires `prep-scg` and `prep-circularity` under the same output directory. To provide an existing assembler circularity table, use `--assembler-circularity TSV`.

**Output:** `04_false_closure/false_closure.contig.tsv`, `false_closure.MAG.tsv`.

## 5. Phantom sequences

```bash
MAGverify phantom -g "$In_MAG" -o "$Out" -t 100
```

Finds sequence stretches of at least 501 bp with zero HiFi read depth, no BLASTN hit to the raw reads and no ambiguous bases. At the MAG level, Meryl also measures the fraction of distinct 21-mers in phantom regions that are absent from the reads.

**Output:** `05_phantom/phantom.region.tsv`, `phantom.contig.tsv`, `phantom.MAG.tsv`.

## 6. Repeats

```bash
MAGverify repeats -g "$In_MAG" -o "$Out" -t 100 -m 900
```

Finds BLASTN repeats (≥ 200 bp, ≥ 80% identity) within contigs (`INTRA`) and between contigs of the same MAG (`INTER`). Flags a contig when its INTRA repeat fraction is ≥ 0.70 or a MAG when its total repeat fraction is ≥ 0.70. Repeat detection and the false-duplication check are separate assessments.

**Output:** `06_repeats/repeats.contig.tsv`, `repeats.MAG.tsv`, `repeats.region.tsv`, `repeats.HSP.tsv`.

## 7. False duplication

```bash
MAGverify false-duplication -g "$In_MAG" -o "$Out" -t 100
```

Identifies duplicated copies within a MAG by BLASTN (≥ 1,000 bp, ≥ 95% identity) and checks each copy boundary against HiFi alignments. A duplicated pair is flagged when one copy is spanned at both boundaries and the other at neither. For tandem copies, both outer boundaries must be spanned but the intervening junction must not be. Clipping evidence at the copy boundaries is also reported.

**Output:** `08_false_duplication/false_duplication.pair.tsv`, `false_duplication.contig.tsv`, `false_duplication.MAG.tsv`.

## 8. Assembly variants

```bash
MAGverify variants -g "$In_MAG" -o "$Out" -t 100 -m 900
```

Profiles the HiFi reads against each assembly with anvi'o and identifies assembly alleles that disagree with the reads. A **minority** variant occurs when another allele or an indel is supported by > 50% of reads. An **unsupported** variant occurs when no read supports the assembly allele. This module can need substantial RAM; use `-m N` to set a memory limit in GB.

**Output:** `09_assembly_variants/variants.site.tsv`, `variants.contig.tsv`, `variants.MAG.tsv`.

## 9. Taxonomic chimera

```bash
MAGverify chimera -g "$In_MAG" -o "$Out" -t 100
```

Classifies genes against GTDB R232 with MMseqs2, then finds neighbouring gene blocks whose dominant taxa differ at domain, phylum, class or order. A transition is called a taxonomic chimera when it overlaps a validated clipping site; high-confidence calls additionally need a stringent clipping site. HiFi structural evidence across the boundary is reported separately.

The database from `prep-chimera` is reused automatically under the same `-o`. Specify `--database PREFIX` to use one elsewhere or `--taxonomy TSV` to use an existing gene taxonomy table. SCG annotation uses `prep-scg`; optional MAG taxonomy uses `prep-taxonomy`.

**Output:** `07_chimera/chimera.region.tsv`, `chimera.contig.tsv`, `chimera.MAG.tsv`.

## 10. Integrated QC summary

```bash
MAGverify summary -g "$In_MAG" -o "$Out"
```

Combines completed module calls into MAG- and contig-level summaries and checks MAG IDs, contig membership, counts and consistency across tables. Validation failures return a nonzero exit status. For a partial run, select completed modules with `--include clipping,chimera`. The summary step does not rerun analyses.

**Output:** `summary.MAG.tsv` and `summary.contig.tsv` (concise `Y` / `N` / `NA` flags); `summary_detail.MAG.tsv` and `summary_detail.contig.tsv` (full diagnostic tables).

For all output column definitions, see [Result tables](OUTPUTS.md).
