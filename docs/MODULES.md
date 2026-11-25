# Modules

This document describes the analytical procedures and QC criteria of each MAGverify module. See [Usage](USAGE.md) for execution commands and [Outputs](OUTPUTS.md) for output files and column definitions.

## 1. Read mapping

The `MAGverify map` command maps HiFi reads to each complete metagenome assembly, including unbinned contigs, using minimap2. All mapped alignments, primary alignments, and primary plus supplementary alignments are stored for subsequent analyses. A mapping manifest records the input and alignment file paths.

## 2. Circularity

The `MAGverify circularity` command evaluates circular topology in single-contig MAGs. Each contig is rotated to place its end-to-start junction at the midpoint. Circular topology is supported when at least two HiFi reads span the junction with at least 2,000 bp aligned on each side. **High-confidence circularity** is defined by agreement between assembler-reported circularity and HiFi read-supported circular topology.

## 3. Clipping

The `MAGverify clipping` command identifies candidate clipping sites using `anvi-script-find-misassemblies` (default clipping ratio > 0.8 and `--min-dist-to-end 0`) and validates them using HiFi read alignments. **High-confidence clipping** is defined by a read depth of ≥ 10× and a clipping ratio of 1.0, as described by [Trigodet et al., 2026](https://doi.org/10.1038/s41587-025-02971-8), with two additional read-based validation criteria: at least two distinct HiFi reads clipped at the same position and no HiFi read spanning that position with aligned bases on both sides. The first and last bases of a contig are not assessable. Clipping sites within 100 bp of either contig end are classified as terminal, and all other sites as internal.

## 4. False closure

The `MAGverify false-closure` command evaluates premature circularization in contigs reported as circular by the assembler. A validated clipping site is identified as a false-closure breakpoint when at least two split HiFi reads on each side connect the circular contig to another contig.

## 5. Phantom sequences

The `MAGverify phantom` command identifies phantom sequences as regions longer than 500 bp with zero HiFi read depth, no BLASTN hit to the raw HiFi reads, and no ambiguous bases. Meryl-based k-mer analysis provides additional evidence by quantifying the MAG-level fraction of distinct 21-mers in phantom regions absent from the HiFi reads.

## 6. Repeats

The `MAGverify repeats` command identifies repeat regions using BLASTN alignments longer than 200 bp with ≥ 80% nucleotide identity. Contig-level assessment evaluates intra-contig repeats, whereas MAG-level assessment evaluates both intra- and inter-contig repeats. **Excessive repeats** are defined as repeat regions covering ≥ 70% of the total contig or MAG sequence at the respective assessment level. Overlapping repeat intervals are counted once.

## 7. False duplication

The `MAGverify false-duplication` command identifies duplicated regions within MAGs using minimap2 assembly self-alignments (`-x asm5 -DP`, alignment length ≥ 1,000 bp and nucleotide identity ≥ 95%) and evaluates HiFi read support across duplication boundaries. False duplication is identified when both boundaries of one duplicated copy are spanned by reads and neither boundary of the other copy is spanned. For tandem duplications, both outer boundaries must be spanned by reads, with no read spanning the intervening junction. Clipping evidence at duplication boundaries is also reported.

## 8. Assembly variants

The `MAGverify variants` command profiles HiFi read alignments using anvi'o to identify assembly alleles inconsistent with the reads. A **minority variant** is defined as a position where at least one read supports the assembly allele, but an alternative allele or indel is supported by > 50% of reads. An **unsupported variant** is defined as a position where no read supports the assembly allele.

## 9. Taxonomic chimera

The `MAGverify chimera` command classifies genes using MMseqs2 against GTDB R232 and identifies transitions between adjacent gene blocks with different dominant taxonomic assignments at the domain, phylum, class, or order level. A **high-confidence taxonomic chimera** is defined by a taxonomic transition containing a high-confidence clipping site.

## 10. Integrated summary

The `MAGverify summary` command integrates completed module results at the MAG and contig levels and validates MAG IDs, contig membership, counts, and consistency across output tables. It does not rerun the analytical modules. Validation failures return a nonzero exit status. Circularity is reported as supporting evidence rather than an assembly defect. `NA` indicates that a check was not performed or was not assessable, rather than a negative result (`N`).
