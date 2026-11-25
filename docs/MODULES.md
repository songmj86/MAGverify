# Modules

This document describes what each module measures and how it decides a call. See [Usage](USAGE.md) for the commands and [Result](OUTPUTS.md) for column definitions.

## 1. Read mapping

Maps HiFi reads separately to each assembler's complete assembly with minimap2, including unbinned contigs. Saves BAMs of all mapped alignments, primary alignments, and primary plus supplementary alignments for later modules. A mapping manifest (`01_mapping/.cache/mapping_manifest.tsv`) records the reads, assemblies and BAM paths so later modules do not need `-a` or `-r`.

## 2. Circularity

Tests single-contig MAGs for read-supported circular topology. The contig is rotated to position its END|START junction in the middle. Circularity is read-supported when at least two HiFi reads span the junction with at least 2,000 bp aligned on each side.

**High-confidence circularity** needs two independent lines of evidence to agree: `High_confidence_circularity` is `Y` only when the assembler reports the contig circular **and** the reads span the junction. An assembler claim with no junction-spanning reads stays `Read_supported_circularity` `N`, and the call is `NA` without the assembler circularity table from `prep-circularity`.

The module tables report `Assembler_reported_circularity` and `Read_supported_circular_topology` separately, and read-supported contigs the assembler did not report circular appear as `Additional_circular_candidate`. `summary.MAG.tsv` and `summary.contig.tsv` carry all three calls.

## 3. Clipping

Finds candidate positions using `anvi-script-find-misassemblies` (default clipping ratio > 0.8, `--min-dist-to-end 0`) and validates them against HiFi alignments. Positions within 100 bp of a contig end are `TERMINAL`; all others are `INTERNAL`.

The module reports three calls. `Clipping` is the read-validated call: `Y` when at least two distinct HiFi reads are clipped at the exact position and no read spans that position with aligned bases on both sides, and `NA` at the first or last base of a contig. `Stringent_candidate` applies the anvi'o thresholds alone, without read validation: `Y` when coverage is >= 10 and the clipping ratio is 1.0.

**High-confidence clipping** is `Stringent_clipping`, which needs both lines of evidence: it is `Y` only when `Clipping` and `Stringent_candidate` are both `Y`, following [Trigodet et al., 2026](https://doi.org/10.1038/s41587-025-02971-8). The thresholds alone are not sufficient, so a position with >= 10x coverage and a clipping ratio of 1.0 that is not validated against the reads is `Stringent_candidate` `Y` but `Stringent_clipping` `N`. `summary.MAG.tsv` and `summary.contig.tsv` carry this call as `High_confidence_clipping`, which is `NA` when clipping was not assessable.

## 4. False closure

Examines contigs reported circular by their assembler. A validated clipping site with at least two split reads on **each** side connecting the contig to another contig indicates a false-closure breakpoint. The module also reports a small-genome proxy for circular contigs < 500,000 bp with >= 3 ribosomal protein families (`prep-scg`). This screening field is distinct from the final `Premature_circularization` call.

## 5. Phantom sequences

Finds sequence stretches longer than 500 bp with zero HiFi read depth, no BLASTN hit to the raw reads and no ambiguous bases. At the MAG level, Meryl also measures the fraction of distinct 21-mers in phantom regions that are absent from the reads.

## 6. Repeats

Finds BLASTN repeats (longer than 200 bp, >= 80% identity). For single-contig MAGs, the contig is scanned against itself, and the MAG call follows that contig. For multi-contig MAGs, each contig is scanned against itself and the contigs are additionally compared with each other.

A contig is flagged when repeats within it cover >= 0.70 of its bases. A MAG is flagged when repeats within and between its contigs together — overlapping bases counted once — cover >= 0.70 of its bases.

## 7. False duplication

Identifies duplicated copies within a MAG from a minimap2 self-alignment of each assembly (`-x asm5 -DP`, >= 1,000 bp, >= 95% identity) and checks each copy boundary against HiFi alignments. A duplicated pair is flagged when one copy is spanned at both boundaries and the other at neither. For tandem copies, both outer boundaries must be spanned but the intervening junction must not be. Clipping evidence at the copy boundaries is also reported.

## 8. Assembly variants

Profiles the HiFi reads against each assembly with anvi'o and identifies assembly alleles that disagree with the reads. A **minority** variant occurs when another allele or an indel is supported by > 50% of reads. An **unsupported** variant occurs when no read supports the assembly allele.

## 9. Taxonomic chimera

Classifies genes against GTDB R232 with MMseqs2, then finds neighbouring gene blocks whose dominant taxa differ at domain, phylum, class or order.

A transition is called `Taxonomic_chimera` when a clipping site with `Clipping` `Y` lies inside it, and `High_confidence_taxonomic_chimera` when that same site is also `Stringent_clipping` `Y`, the high-confidence clipping call. `summary.MAG.tsv` and `summary.contig.tsv` carry the high-confidence call as `Chimera`. The HiFi structural check across the boundary is reported separately as `Structural_chimera`.

## 10. Integrated QC summary

Combines completed module calls into MAG- and contig-level summaries and checks MAG IDs, contig membership, counts and consistency across tables. Validation failures return a nonzero exit status. The summary step does not rerun analyses.

Circularity is reported as evidence, not as a defect. Read calls that are `NA` — a check that was not run, or not assessable for that MAG or contig — separately from `N`.
