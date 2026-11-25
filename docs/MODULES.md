# Modules

This document describes the analytical methods and detection criteria for each MAGverify module. See [Usage](USAGE.md) for execution commands and [Outputs](OUTPUTS.md) for output files and column definitions.

## 1. Read mapping (`map`)

HiFi reads are mapped to each metagenome assembly using minimap2, including unbinned contigs. Three sets of alignments are retained for downstream analyses: all mapped alignments (clipping, phantom sequences, assembly variants), primary plus supplementary alignments (read-spanning analyses), and primary alignments (depth).

## 2. Circularity (`circularity`)

Single-contig MAGs are tested for circular topology. The contig is rotated so that its original END|START junction lies in the middle, and HiFi alignments spanning that junction are counted.

| Criterion | Definition |
| --- | --- |
| Assembler-reported circularity | Circular topology reported by the assembler |
| Read-supported circularity | ≥ 2 HiFi reads spanning the circular junction with ≥ 2,000 bp aligned on each side |
| High-confidence circularity | Both assembler-reported and read-supported circularity |

Without the assembler circularity table from `prep-circularity`, assembler-reported circularity is unavailable and the high-confidence call is `NA`. Read-supported contigs that the assembler did not report circular are reported as `Additional_circular_candidate`.

## 3. Clipping (`clipping`)

Candidate positions are identified with `anvi-script-find-misassemblies` (clipping ratio > 0.8, `--min-dist-to-end 0`) and validated against the HiFi alignments. Positions within 100 bp of a contig end are terminal; all others are internal.

| Criterion | Definition |
| --- | --- |
| `Clipping` | ≥ 2 distinct HiFi reads clipped at the position and no read spanning it with aligned bases on both sides; `NA` at the first or last base of a contig |
| `Clipping_candidate` | Read depth ≥ 10× and clipping ratio = 1.0, without read validation |
| `High_confidence_clipping` | Both `Clipping` and `Clipping_candidate`, following [Trigodet et al., 2026](https://doi.org/10.1038/s41587-025-02971-8) |

The thresholds alone are not sufficient: a position with ≥ 10× coverage and a clipping ratio of 1.0 that the reads do not validate is `Clipping_candidate` `Y` but `High_confidence_clipping` `N`. The call is `NA` when clipping was not assessable.

## 4. False closure (`false-closure`)

Examines contigs reported circular by their assembler. A validated clipping site with at least two split reads on **each** side connecting the contig to another contig indicates a false-closure breakpoint. The module also reports a small-genome proxy for circular contigs < 500,000 bp with ≥ 3 ribosomal protein families (`prep-scg`). This screening field is distinct from the final `Premature_circularization` call.

## 5. Phantom sequences (`phantom`)

A region is called a phantom sequence when it is longer than 500 bp and has zero HiFi read depth, no BLASTN hit to the raw reads and no ambiguous bases. Meryl-based k-mer analysis provides additional evidence for evaluating phantom sequences: at the MAG level it measures the fraction of distinct 21-mers in phantom regions that are absent from the reads.

## 6. Repeats (`repeats`)

Repeats are BLASTN alignments longer than 200 bp with ≥ 80% identity.

| Level | Repeat assessment |
| --- | --- |
| Contig | Intra-contig repeats |
| MAG | Intra- and inter-contig repeats |

A contig is flagged when repeats cover ≥ 70% of its sequence, and a MAG when they cover ≥ 70% of its sequence. Overlapping repeat intervals are counted once. A single-contig MAG has no inter-contig repeats, so its MAG call follows that contig.

## 7. False duplication (`false-duplication`)

Identifies duplicated copies within a MAG from a minimap2 self-alignment of each assembly (`-x asm5 -DP`, ≥ 1,000 bp, ≥ 95% identity) and checks each copy boundary against HiFi alignments. A duplicated pair is flagged when one copy is spanned at both boundaries and the other at neither. For tandem copies, both outer boundaries must be spanned but the intervening junction must not be. Clipping evidence at the copy boundaries is also reported.

## 8. Assembly variants (`variants`)

The HiFi reads are profiled against each assembly with anvi'o to find positions where the assembly allele disagrees with the reads.

| Variant type | Definition |
| --- | --- |
| Minority variant | The assembly allele is supported by at least one read, but an alternative allele or indel is supported by > 50% of reads |
| Unsupported variant | No read supports the assembly allele |

The two types are mutually exclusive.

## 9. Taxonomic chimera (`chimera`)

Genes are classified against GTDB R232 with MMseqs2. At each gene boundary, up to ten genes on either side are searched and the nearest five classified genes form a block, with at least four required. A block is assigned a dominant taxon when ≥ 80% of its genes agree at a rank, and a boundary whose two blocks have different dominant taxa at domain, phylum, class or order is a taxonomic transition.

| Call | Evidence |
| --- | --- |
| `Clipping_supported_transition` | Taxonomic transition containing a `Clipping` `Y` site |
| `High_confidence_taxonomic_chimera` | The same site is also `High_confidence_clipping` `Y` |

`High_confidence_taxonomic_chimera` is the chimera call carried by `summary.MAG.tsv` and `summary.contig.tsv` as `Chimera`.

`Structural_discontinuity` is reported independently of both and does not use taxonomy. With 500 bp anchors on each side, it is `N` when ≥ 2 reads cross the boundary, `Y` when ≥ 2 reads anchor on each side but none cross it, and `NA` otherwise.

## 10. Integrated summary (`summary`)

Combines completed module calls into MAG- and contig-level summaries and checks MAG IDs, contig membership, counts and consistency across tables. Validation failures return a nonzero exit status. The summary step does not rerun analyses.

Circularity flags indicate evidence of circular topology rather than assembly defects. `NA` indicates that a check was not performed or was not assessable and should not be interpreted as `N`.
