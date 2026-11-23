# Modules

Run the [preparation](#preparation) commands first. `MAGverify end-to-end` then runs analysis steps 1–10 below in this order (see the [README](../README.md#4-analysis-end-to-end)). Each step can also be run on its own with the commands here, for example to rerun one module or to run only some of them. Use the **same `-o`** for every command: each module finds the results of earlier steps (mapping manifest, clipping tables, preparation tables) there automatically. The column definitions of every output table are in [Result tables](OUTPUTS.md).

The examples use these variables:

```bash
In_Contig=/path/to/assemblies
In_Read=/path/to/hifi_reads
In_MAG=/path/to/MAGs
Out=/path/to/MAGverify_output
```

| Step | Command | Needs first |
| --- | --- | --- |
| Preparation (separately, first) | `prep-chimera` (once), `prep-scg`, `prep-circularity`, `prep-taxonomy` | — |
| 1 | `map` | — |
| 2 | `circularity` | `map` |
| 3 | `clipping` | `map` |
| 4 | `false-closure` | `clipping`, `prep-scg`, `prep-circularity` |
| 5 | `phantom` | `map` |
| 6 | `repeats` | `map` |
| 7 | `false-duplication` | `clipping` |
| 8 | `variants` | `map` |
| 9 | `chimera` | `clipping`, `prep-chimera` (and `prep-scg`, `prep-taxonomy` for SCG redundancy) |
| 10 | `summary` | the modules to summarize |

**Threads and memory.** `-t N` sets the threads (default: all available CPUs, at most 100). Read mapping (`map`, `circularity`) and anvi'o profiling (`variants`) process one assembler at a time with all threads; `clipping`, `phantom` and `false-duplication` process assemblers in parallel and divide `-t` between them (`PARALLEL_JOBS=N` before the command runs fewer at once). `-m N` (GB) keeps `map`, `circularity`, `repeats` and `variants` within a memory limit. Read mapping and the BLAST steps of `repeats` print progress with the remaining time every 5 minutes (`PROGRESS_INTERVAL=SECONDS` changes the interval).

## Preparation

Not part of `end-to-end`; run these first with the same `-o`.

```bash
MAGverify prep-circularity -c FILE [-c FILE ...] -g $In_MAG -o $Out
MAGverify prep-taxonomy    -d /path/to/GTDBTk_results -g $In_MAG -o $Out
MAGverify prep-scg         -a $In_Contig -o $Out -t 100
MAGverify prep-chimera     -d /path/to/GTDB_R232_downloads -o $Out -t 100
```

| Command | Output in `00_prep/` | Used by |
| --- | --- | --- |
| `prep-circularity` | `assembler_circularity.MAG.tsv`, `assembler_circularity.contig.tsv` | `circularity` (optional), `false-closure` |
| `prep-taxonomy` | `gtdbtk_taxonomy.MAG.tsv` | `chimera` (optional) |
| `prep-scg` | `scg.contig.tsv`, `genes.*` | `false-closure`, `chimera` |
| `prep-chimera` | `GTDB_R232*` (MMseqs2 database) | `chimera` |

Inputs, file formats and options of each step are in [Input preparation](INPUT_PREPARATION.md).

## 1. Read mapping

```bash
MAGverify map -a $In_Contig -r $In_Read -o $Out -t 100 -m 900
```

Maps the HiFi reads with minimap2 separately to each assembler FASTA (whole assemblies, so alignments to unbinned contigs are kept), and keeps one BAM of all mapped alignments plus BAMs of the primary and of the primary + supplementary alignments for the later modules. The reads, assemblies and BAMs are recorded in `01_mapping/.cache/mapping_manifest.tsv`, so later modules need neither `-a` nor `-r`.

Output: `01_mapping/mapping.assembler.tsv` (alignment counts per assembler).

## 2. Circularity

```bash
MAGverify circularity -g $In_MAG -o $Out -t 100 -m 900
```

Tests whether single-contig MAGs are circular with read evidence: the contig is rotated so that its END|START junction lies in the middle, the HiFi reads of its assembler are mapped to it, and the junction is supported when >= 2 reads span it with >= 2,000 bp aligned on each side. The assembler circularity table from `prep-circularity` is optional. Without it, `Circular` is the read-supported call. With it, `Circular` is Y only when the assembler and the reads agree (`High_confidence_circular`), and MAGs the reads support but the assembler did not call circular are listed as `Additional_circular_candidate`.

Output: `02_circularity/circularity.MAG.tsv`, `circularity.contig.tsv`.

## 3. Clipping

```bash
MAGverify clipping -g $In_MAG -o $Out -t 100
```

Finds positions where many reads are clipped with `anvi-script-find-misassemblies` (default candidates: clipping ratio > 0.8) and validates each position with the reads: `Clipping` is Y when >= 2 distinct HiFi reads are clipped at the exact position and no read spans it. Stringent sites also have anvi'o coverage >= 10 and clipping ratio 1.0. Positions within 100 bp of a contig end are TERMINAL, the others INTERNAL. `false-closure`, `false-duplication`, `chimera` and `summary` use these sites.

Output: `03_clipping/clipping.site.tsv`, `clipping.contig.tsv`, `clipping.MAG.tsv`, `clipping.assembler.tsv`.

## 4. False closure

```bash
MAGverify false-closure -g $In_MAG -o $Out -t 100
```

Tests contigs reported circular by their assembler: a validated clipping site where >= 2 reads on **each** side are split between this contig and another contig indicates that the circular contig joins sequences that do not belong together. Also screens small circular genomes (< 500,000 bp with >= 3 ribosomal protein families from `prep-scg`) as premature circularization candidates. Requires `prep-scg` and `prep-circularity` with the same `-o`; another assembler circularity table can be given with `--assembler-circularity TSV`.

Output: `04_false_closure/false_closure.contig.tsv`, `false_closure.MAG.tsv`.

## 5. Phantom sequences

```bash
MAGverify phantom -g $In_MAG -o $Out -t 100
```

Finds stretches >= 501 bp with zero HiFi read depth that also have no BLASTN hit to the raw reads (phantom regions), and measures, per MAG, the fraction of their 21-mers that are absent from the reads (Meryl).

Output: `05_phantom/phantom.region.tsv`, `phantom.contig.tsv`, `phantom.MAG.tsv`.

## 6. Repeats

```bash
MAGverify repeats -g $In_MAG -o $Out -t 100 -m 900
```

Finds repeats with BLASTN (>= 200 bp, >= 80% identity) within each contig (INTRA) and between contigs of the same MAG (INTER), and reports the repeat fractions. A contig with an INTRA repeat fraction >= 0.70, or a MAG with a total repeat fraction >= 0.70, is flagged as excessive repeat. Repeat detection and the false-duplication test are separate assessments.

Output: `06_repeats/repeats.contig.tsv`, `repeats.MAG.tsv`, `repeats.region.tsv`, `repeats.HSP.tsv`.

## 7. False duplication

```bash
MAGverify false-duplication -g $In_MAG -o $Out -t 100
```

Finds duplicated copies (BLASTN, >= 1,000 bp, >= 95% identity) within each MAG and tests each copy boundary with HiFi reads. A pair is a false duplication when one copy is spanned by reads at both boundaries and the other at neither (for tandem copies: both outer boundaries spanned and the junction between the copies not), i.e. the reads support only one copy. Clipping at the copy boundaries is reported alongside.

Output: `08_false_duplication/false_duplication.pair.tsv`, `false_duplication.contig.tsv`, `false_duplication.MAG.tsv`.

## 8. Assembly variants

```bash
MAGverify variants -g $In_MAG -o $Out -t 100 -m 900
```

Profiles the HiFi reads on each assembly with anvi'o and reports positions where the reads disagree with the assembly: **minority** assembly variants (another allele or the indel carried by > 50% of reads) and **unsupported** assembly variants (no read supports the assembly allele). This step may need substantial RAM; set `-m N` (GB) when appropriate.

Output: `09_assembly_variants/variants.site.tsv`, `variants.contig.tsv`, `variants.MAG.tsv`.

## 9. Taxonomic chimera

```bash
MAGverify chimera -g $In_MAG -o $Out -t 100 --database /path/to/database_dir/00_prep/GTDB_R232
```

Classifies the predicted genes of all assemblies against GTDB R232 (MMseqs2) and looks for gene boundaries where the dominant taxon of the neighbouring genes changes at domain, phylum, class or order. A transition is a taxonomic chimera when a validated clipping site lies in it, and high-confidence when that site is stringent; the HiFi structural evidence across the junction is reported alongside. The database (`prep-chimera`) and the MAG taxonomy table (`prep-taxonomy`) are read from `00_prep/` when they were built with the same `-o`; `--database PREFIX` gives a database elsewhere, and `--taxonomy TSV` a compatible gene taxonomy table that already exists. SCG annotation comes from `prep-scg`.

Output: `07_chimera/chimera.region.tsv`, `chimera.contig.tsv`, `chimera.MAG.tsv`.

## 10. Integrated QC summary

```bash
MAGverify summary -g $In_MAG -o $Out
```

Collects the calls of the completed modules into one row per MAG and per contig and validates MAG IDs, contig memberships, counts and table consistency; validation failures return a nonzero exit status. For a partial run, give the modules that were completed, for example `--include clipping,chimera`. The summary does not rerun analysis modules.

Output: `summary.MAG.tsv`, `summary.contig.tsv` at the top of the output folder.
