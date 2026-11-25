# Result

MAGverify writes four integrated tables to the output directory and module-specific tables to numbered subdirectories. This document defines every reported column, including diagnostic calls, supporting read evidence and measurements. See the [README](../README.md#5-main-results) for the output overview and [Usage](USAGE.md#4-individual-modules) for the corresponding commands.

**Table naming:** `MODULE.LEVEL.tsv`, where `LEVEL` identifies one row (for example, `MAG`, `contig`, `site`, `region`, `pair` or `assembler`). The four integrated `summary*.tsv` files are saved at the top of the output directory.

**Shared conventions**

- `Y` / `N` / `NA`: positive / negative / not assessable or not run. Columns with these values are calls; `n_*` fields are counts.
- Read counts (`*_reads`) count distinct read names. Unless noted otherwise, read evidence uses primary and supplementary HiFi alignments to the assembly of the contig. Fields beginning with `HiFi_all_*` also include secondary alignments.
- `I` / `D` / `N`: insertion, deletion or reference skip (CIGAR operation) in a read alignment at the tested position.
- Positions are 1-based unless the column name says `0based`. Lists within one cell are separated by `;`.
- `INTERNAL` / `TERMINAL`: more than 100 bp from a contig end / within 100 bp of an end.

**Keep the `.cache/` directories.** They store intermediate files, fingerprints and other data reused by the analysis modules. Deleting them can prevent normal resumption or reuse.

## Contents

<pre>
MAGverify_output
├── <a href="#summarymagtsv">summary.MAG.tsv</a>
├── <a href="#summarycontigtsv">summary.contig.tsv</a>
├── <a href="#detail-tables">detail.MAG.tsv</a>
├── <a href="#detail-tables">detail.contig.tsv</a>
├── 01_mapping
│   ╰── <a href="#mappingassemblertsv">mapping.assembler.tsv</a>
├── 02_circularity
│   ├── <a href="#circularitymagtsv">circularity.MAG.tsv</a>
│   ╰── <a href="#circularitycontigtsv">circularity.contig.tsv</a>
├── 03_clipping
│   ├── <a href="#clippingmagtsv">clipping.MAG.tsv</a>
│   ├── <a href="#clippingcontigtsv">clipping.contig.tsv</a>
│   ├── <a href="#clippingsitetsv">clipping.site.tsv</a>
│   ╰── <a href="#clippingassemblertsv">clipping.assembler.tsv</a>
├── 04_false_closure
│   ├── <a href="#false_closuremagtsv">false_closure.MAG.tsv</a>
│   ╰── <a href="#false_closurecontigtsv">false_closure.contig.tsv</a>
├── 05_phantom
│   ├── <a href="#phantommagtsv">phantom.MAG.tsv</a>
│   ├── <a href="#phantomcontigtsv">phantom.contig.tsv</a>
│   ╰── <a href="#phantomregiontsv">phantom.region.tsv</a>
├── 06_repeats
│   ├── <a href="#repeatsmagtsv">repeats.MAG.tsv</a>
│   ├── <a href="#repeatscontigtsv">repeats.contig.tsv</a>
│   ├── <a href="#repeatsregiontsv">repeats.region.tsv</a>
│   ╰── <a href="#repeatshsptsv">repeats.HSP.tsv</a>
├── 07_chimera
│   ├── <a href="#chimeramagtsv">chimera.MAG.tsv</a>
│   ├── <a href="#chimeracontigtsv">chimera.contig.tsv</a>
│   ╰── <a href="#chimeraregiontsv">chimera.region.tsv</a>
├── 08_false_duplication
│   ├── <a href="#false_duplicationmagtsv">false_duplication.MAG.tsv</a>
│   ├── <a href="#false_duplicationcontigtsv">false_duplication.contig.tsv</a>
│   ╰── <a href="#false_duplicationpairtsv">false_duplication.pair.tsv</a>
╰── 09_assembly_variants
    ├── <a href="#variantsmagtsv">variants.MAG.tsv</a>
    ├── <a href="#variantscontigtsv">variants.contig.tsv</a>
    ╰── <a href="#variantssitetsv">variants.site.tsv</a>
</pre>

## 1. Summary

Four integrated files are written at the top of the output directory:

- `summary.MAG.tsv` and `summary.contig.tsv`: concise Y/N/NA flags in end-to-end module order, with `Assembly_variant` last.
- `detail.MAG.tsv` and `detail.contig.tsv`: the same flags, diagnostic annotations, and every module column joined as `MODULE.COLUMN`.

`Y` indicates a positive call, `N` a negative call, and `NA` a check that was unavailable or not assessable. Circularity support fields describe positive circularity evidence rather than a defect. Counts and skipped-check names are kept in the detailed summaries.

### `summary.MAG.tsv`

One row per MAG; columns appear in the following order.

| Column | Definition |
|---|---|
| `MAG_ID` | MAG identifier |
| `Assembler_reported_circularity` | Y/N/NA call from `prep-circularity` |
| `Read_supported_circularity` | HiFi read evidence supporting circular topology, from `Read_supported_circular_topology` |
| `High_confidence_circularity` | High-confidence circularity: both assembler-reported and read-supported calls are Y; NA when the assembler call is unavailable |
| `High_confidence_clipping` | At least one high-confidence clipping site: `High_confidence_clipping` Y in the clipping module (validated clipping with anvi'o coverage >= 10 and clipping ratio = 1.0); NA when clipping was not assessable or not run |
| `Premature_circularization` | Single-contig MAG reported circular by its assembler with false-closure evidence |
| `Phantom_sequence` | At least one phantom region |
| `Excessive_repeat` | MAG-level excessive repeat flag |
| `False_duplication` | At least one false-duplication event |
| `Chimera` | High-confidence taxonomic chimera (`High_confidence_taxonomic_chimera`) from the chimera module |
| `Assembly_variant` | At least one assembly-allele disagreement |

### `summary.contig.tsv`

One row per contig across all assemblies. `MAG_ID` is `NA` for unbinned contigs. After the three identity columns the flags are the same, in the same order, as in `summary.MAG.tsv`.

| Column | Definition |
|---|---|
| `Contig` | Contig identifier |
| `MAG_ID` | Associated MAG identifier, or `NA` |
| `Length_bp` | Contig length in base pairs |
| `Assembler_reported_circularity` | Y/N/NA call from `prep-circularity` |
| `Read_supported_circularity` | HiFi read evidence supporting circular topology, from `Read_supported_circular_topology`; `NA` for contigs of multi-contig MAGs, which are not evaluated |
| `High_confidence_circularity` | High-confidence circularity: both assembler-reported and read-supported calls are Y; NA when the assembler call is unavailable or the contig was not evaluated |
| `High_confidence_clipping` | At least one high-confidence clipping site (`High_confidence_clipping` Y in the clipping module); NA when clipping was not assessable or not run |
| `Premature_circularization` | Contig reported circular by its assembler with false-closure evidence, from `Premature_circularization_candidate`; NA when false closure was not run |
| `Phantom_sequence` | Phantom sequence flag |
| `Excessive_repeat` | Intra-contig excessive repeat flag |
| `False_duplication` | False-duplication flag |
| `Chimera` | High-confidence taxonomic chimera (`High_confidence_taxonomic_chimera`) |
| `Assembly_variant` | Assembly-allele disagreement flag |

### Detail tables

The detail tables are grouped by module. `detail.MAG.tsv` opens with the identity columns `MAG_ID` and `n_contigs`, then carries one block per module in end-to-end order: circularity, clipping, false closure, phantom, repeats, false duplication, chimera, variants. Each block holds that module's summary flag and its diagnostic annotations, followed by every column of that module's MAG table as `MODULE.COLUMN`.

| Module block | Flags and annotations before the `MODULE.COLUMN` fields |
|---|---|
| circularity | `Assembler_reported_circularity`, `Read_supported_circularity`, `High_confidence_circularity` |
| clipping | `Clipping`, `High_confidence_clipping`, `n_clipping_sites` |
| false closure | `Premature_circularization`, `False_closure`, `Premature_circularization_screen`, `Small_genome_proxy`, `False_closure_partner_contigs` |
| phantom | `Phantom_sequence`, `n_phantom_contigs` |
| repeats | `Excessive_repeat`, `Repeat_fraction`, `n_excessive_repeat_contigs` |
| false duplication | `False_duplication`, `False_duplication_with_clipping`, `n_false_duplication_contigs`, `n_duplicated_pairs`, `n_unresolved_duplicated_pairs` |
| chimera | `Chimera`, `n_chimera_contigs` |
| variants | `Assembly_variant`, `n_variant_contigs` |

`detail.contig.tsv` is built the same way from the contig-level tables and includes only contigs assigned to MAGs. Its identity columns are `Contig`, `MAG_ID`, `Length_bp` and `Assembler`, and its module blocks carry the same names in the same order as the MAG table. Three MAG-level counts have no contig equivalent and are absent (`n_phantom_contigs`, `n_excessive_repeat_contigs`, `n_variant_contigs`), `n_false_duplication_contigs` is replaced by `n_false_duplication_pairs`, and `n_chimera_contigs` by `Chimera_ranks`. For the full meanings of `MODULE.COLUMN`, consult the relevant module table section below.

`Clipping` is the read-validated call (`Clipping` Y at >= 1 site), kept alongside `High_confidence_clipping` so both clipping tiers are visible. The detail tables carry no check counters and no list of skipped checks; a check that was not run shows `NA` in its own column.

They also retain the following summary-level fields that are omitted from the concise summaries:

| Column | Definition |
|---|---|
| `False_closure` | General false-closure check, distinct from `Premature_circularization` |
| `Premature_circularization_screen` | Combined screen for premature circularization or the small-genome proxy, not identical to the final premature-circularization flag |
| `Small_genome_proxy` | Small-genome indicator from false closure |
| `False_duplication_with_clipping` | Co-occurrence of false duplication and validated clipping; additional evidence, not a separate primary diagnostic flag |

Other annotations include supporting counts, affected contigs, and repeat fractions. The module analysis algorithms and original call criteria are unchanged.

## 2. Read mapping

**Directory:** `01_mapping/` · **Command:** `map`

### `mapping.assembler.tsv`

One row per assembler. `samtools flagstat` counts of the BAM with all mapped HiFi alignments (unmapped reads are not kept).

| Column | Definition |
|---|---|
| `Assembler` | Assembly file (assembler) the contig comes from |
| `Alignments_total` | All alignment records (primary + secondary + supplementary) |
| `Mapped_reads` | Mapped reads (primary alignments) |
| `Secondary_alignments` | Secondary alignment records |
| `Supplementary_alignments` | Supplementary alignment records (split parts of chimeric/split reads) |

## 3. Circularity

**Directory:** `02_circularity/` · **Command:** `circularity`

### `circularity.MAG.tsv`

One row per MAG. Only single-contig MAGs can be evaluated: the contig is rotated so that its END|START junction lies in the middle, HiFi reads of its assembler are mapped to it, and the reads crossing the junction are counted. Columns from `Terminal_overlap_bp` on are NA for multi-contig MAGs.

| Column | Definition |
|---|---|
| `MAG_ID` | MAG ID (MAG FASTA file name without extension); `NA` for unbinned contigs |
| `Assembler_reported_circularity` | Circularity reported by the assembler (`prep-circularity` table): Y, N, or NA when not available |
| `Circularity_evaluable` | Y when the MAG is a single contig and both contig ends are reached by >= 2 HiFi reads (primary + supplementary), so the END\|START junction can be tested |
| `Read_supported_circular_topology` | Y when >= 2 HiFi reads (primary + supplementary) span the END\|START junction with >= 2,000 bp aligned on each side and no deletion/skip across it; N when evaluable but not spanned; NA when not evaluable |
| `High_confidence_circularity` | Y when the assembler reports the contig circular **and** the reads support the junction; NA without the assembler circularity table |
| `Additional_circular_candidate` | Y when the reads support the junction but the assembler did not report the contig circular; NA without the assembler circularity table |
| `Circular` | Final call: `High_confidence_circularity` with the assembler table, otherwise `Read_supported_circular_topology` |
| `n_contigs` | Number of contigs in the MAG |
| `Contig` | Contig of a single-contig MAG; NA for multi-contig MAGs |
| `Terminal_overlap_bp` | Longest exact overlap between the contig end and start (checked up to 100,000 bp) |
| `HiFi_END_reaching_reads` | Reads (HiFi primary and supplementary alignments) with >= 2,000 bp aligned on the END side of the junction |
| `HiFi_START_reaching_reads` | Reads (HiFi primary and supplementary alignments) with >= 2,000 bp aligned on the START side of the junction |
| `HiFi_analyzable` | Y when both `HiFi_END_reaching_reads` and `HiFi_START_reaching_reads` are >= 2 |
| `HiFi_END_START_spanning_reads` | Reads (HiFi primary and supplementary alignments) spanning the junction with >= 2,000 bp on each side and no deletion/skip across it |
| `HiFi_END_START_supported` | Y when `HiFi_END_START_spanning_reads` >= 2; N when analyzable but fewer; NA otherwise |
| `HiFi_all_END_reaching_reads` | As `HiFi_END_reaching_reads`, counting all alignments (secondary included) |
| `HiFi_all_START_reaching_reads` | As `HiFi_START_reaching_reads`, counting all alignments |
| `HiFi_all_analyzable` | As `HiFi_analyzable`, counting all alignments |
| `HiFi_all_END_START_spanning_reads` | As `HiFi_END_START_spanning_reads`, counting all alignments |
| `HiFi_all_END_START_supported` | As `HiFi_END_START_supported`, counting all alignments |
| `HiFi_primary_END_START_spanning_reads` | Spanning reads with a primary alignment across the junction |
| `HiFi_supplementary_END_START_spanning_reads` | Spanning reads with a supplementary alignment across the junction |
| `HiFi_secondary_END_START_spanning_reads` | Spanning reads with a secondary alignment across the junction |
| `HiFi_secondary_added_END_START_support` | Y when the junction is supported only after secondary alignments are added |
| `HiFi_clean_END_START_spanning_reads` | Spanning reads (HiFi primary and supplementary alignments) with no insertion, deletion or skip at the junction |
| `HiFi_clean_junction` | Y when `HiFi_clean_END_START_spanning_reads` >= 2; N when analyzable but fewer; NA otherwise |
| `HiFi_all_clean_END_START_spanning_reads` | As `HiFi_clean_END_START_spanning_reads`, counting all alignments |
| `HiFi_all_clean_junction` | As `HiFi_clean_junction`, counting all alignments |
| `HiFi_junction_I_reads` | Reads (HiFi primary and supplementary alignments) with an insertion (I) at the junction |
| `HiFi_junction_I_event_count` | Number of I events at the junction |
| `HiFi_junction_I_total_bp` | Total length of the I events (bp) |
| `HiFi_max_junction_I_bp` | Longest I event (bp) |
| `HiFi_junction_D_reads` | Reads (HiFi primary and supplementary alignments) with a deletion (D) at the junction |
| `HiFi_junction_D_event_count` | Number of D events at the junction |
| `HiFi_junction_D_total_bp` | Total length of the D events (bp) |
| `HiFi_max_junction_D_bp` | Longest D event (bp) |
| `HiFi_junction_N_reads` | Reads (HiFi primary and supplementary alignments) with a reference skip (N) at the junction |
| `HiFi_junction_N_event_count` | Number of N events at the junction |
| `HiFi_junction_N_total_bp` | Total length of the N events (bp) |
| `HiFi_max_junction_N_bp` | Longest N event (bp) |
| `HiFi_all_junction_I_reads` | As `HiFi_junction_I_reads`, counting all alignments |
| `HiFi_all_junction_D_reads` | As `HiFi_junction_D_reads`, counting all alignments |
| `HiFi_all_junction_N_reads` | As `HiFi_junction_N_reads`, counting all alignments |
| `HiFi_all_max_junction_I_bp` | As `HiFi_max_junction_I_bp`, counting all alignments |
| `HiFi_all_max_junction_D_bp` | As `HiFi_max_junction_D_bp`, counting all alignments |
| `HiFi_all_max_junction_N_bp` | As `HiFi_max_junction_N_bp`, counting all alignments |

### `circularity.contig.tsv`

One row per contig of every assembly. The read evidence columns are filled for the contig of a single-contig MAG and NA for the other contigs.

| Column | Definition |
|---|---|
| `Contig` | Contig ID (first word of the FASTA header) |
| `MAG_ID` | MAG ID (MAG FASTA file name without extension); `NA` for unbinned contigs |
| `Assembler` | Assembly file (assembler) the contig comes from |
| `Length_bp` | Contig length (bp) |
| `n_contigs` | Number of contigs in the MAG of this contig (0 for unbinned contigs) |
| `Assembler_reported_circularity` | Circularity reported by the assembler (`prep-circularity` table): Y, N, or NA when not available |
| `Circularity_evaluable` | Y when the MAG is a single contig and both contig ends are reached by >= 2 HiFi reads (primary + supplementary), so the END\|START junction can be tested |
| `Read_supported_circular_topology` | Y when >= 2 HiFi reads (primary + supplementary) span the END\|START junction with >= 2,000 bp aligned on each side and no deletion/skip across it; N when evaluable but not spanned; NA when not evaluable |
| `High_confidence_circularity` | Y when the assembler reports the contig circular **and** the reads support the junction; NA without the assembler circularity table |
| `Additional_circular_candidate` | Y when the reads support the junction but the assembler did not report the contig circular; NA without the assembler circularity table |
| `Circular` | Final call: `High_confidence_circularity` with the assembler table, otherwise `Read_supported_circular_topology` |
| `Terminal_overlap_bp` | Longest exact overlap between the contig end and start (checked up to 100,000 bp) |
| `HiFi_END_reaching_reads` | Reads (HiFi primary and supplementary alignments) with >= 2,000 bp aligned on the END side of the junction |
| `HiFi_START_reaching_reads` | Reads (HiFi primary and supplementary alignments) with >= 2,000 bp aligned on the START side of the junction |
| `HiFi_analyzable` | Y when both `HiFi_END_reaching_reads` and `HiFi_START_reaching_reads` are >= 2 |
| `HiFi_END_START_spanning_reads` | Reads (HiFi primary and supplementary alignments) spanning the junction with >= 2,000 bp on each side and no deletion/skip across it |
| `HiFi_END_START_supported` | Y when `HiFi_END_START_spanning_reads` >= 2; N when analyzable but fewer; NA otherwise |
| `HiFi_all_END_reaching_reads` | As `HiFi_END_reaching_reads`, counting all alignments (secondary included) |
| `HiFi_all_START_reaching_reads` | As `HiFi_START_reaching_reads`, counting all alignments |
| `HiFi_all_analyzable` | As `HiFi_analyzable`, counting all alignments |
| `HiFi_all_END_START_spanning_reads` | As `HiFi_END_START_spanning_reads`, counting all alignments |
| `HiFi_all_END_START_supported` | As `HiFi_END_START_supported`, counting all alignments |
| `HiFi_primary_END_START_spanning_reads` | Spanning reads with a primary alignment across the junction |
| `HiFi_supplementary_END_START_spanning_reads` | Spanning reads with a supplementary alignment across the junction |
| `HiFi_secondary_END_START_spanning_reads` | Spanning reads with a secondary alignment across the junction |
| `HiFi_secondary_added_END_START_support` | Y when the junction is supported only after secondary alignments are added |
| `HiFi_clean_END_START_spanning_reads` | Spanning reads (HiFi primary and supplementary alignments) with no insertion, deletion or skip at the junction |
| `HiFi_clean_junction` | Y when `HiFi_clean_END_START_spanning_reads` >= 2; N when analyzable but fewer; NA otherwise |
| `HiFi_all_clean_END_START_spanning_reads` | As `HiFi_clean_END_START_spanning_reads`, counting all alignments |
| `HiFi_all_clean_junction` | As `HiFi_clean_junction`, counting all alignments |
| `HiFi_junction_I_reads` | Reads (HiFi primary and supplementary alignments) with an insertion (I) at the junction |
| `HiFi_junction_I_event_count` | Number of I events at the junction |
| `HiFi_junction_I_total_bp` | Total length of the I events (bp) |
| `HiFi_max_junction_I_bp` | Longest I event (bp) |
| `HiFi_junction_D_reads` | Reads (HiFi primary and supplementary alignments) with a deletion (D) at the junction |
| `HiFi_junction_D_event_count` | Number of D events at the junction |
| `HiFi_junction_D_total_bp` | Total length of the D events (bp) |
| `HiFi_max_junction_D_bp` | Longest D event (bp) |
| `HiFi_junction_N_reads` | Reads (HiFi primary and supplementary alignments) with a reference skip (N) at the junction |
| `HiFi_junction_N_event_count` | Number of N events at the junction |
| `HiFi_junction_N_total_bp` | Total length of the N events (bp) |
| `HiFi_max_junction_N_bp` | Longest N event (bp) |
| `HiFi_all_junction_I_reads` | As `HiFi_junction_I_reads`, counting all alignments |
| `HiFi_all_junction_D_reads` | As `HiFi_junction_D_reads`, counting all alignments |
| `HiFi_all_junction_N_reads` | As `HiFi_junction_N_reads`, counting all alignments |
| `HiFi_all_max_junction_I_bp` | As `HiFi_max_junction_I_bp`, counting all alignments |
| `HiFi_all_max_junction_D_bp` | As `HiFi_max_junction_D_bp`, counting all alignments |
| `HiFi_all_max_junction_N_bp` | As `HiFi_max_junction_N_bp`, counting all alignments |

## 4. Clipping

**Directory:** `03_clipping/` · **Command:** `clipping`

### `clipping.site.tsv`

One row per clipping candidate reported by `anvi-script-find-misassemblies` (default clipping ratio > 0.8, no minimum coverage, `--min-dist-to-end 0`) on the HiFi reads of each assembler, for all contigs and both region types.

| Column | Definition |
|---|---|
| `Contig` | Contig ID (first word of the FASTA header) |
| `MAG_ID` | MAG ID (MAG FASTA file name without extension); `NA` for unbinned contigs |
| `Assembler` | Assembly file (assembler) the contig comes from |
| `Length_bp` | Contig length (bp) |
| `Position` | Candidate position (1-based) |
| `Region_type` | TERMINAL within 100 bp of a contig end; INTERNAL otherwise |
| `Distance_to_nearest_end_bp` | Distance from the position to the nearer contig end (bp) |
| `Reference_boundary` | Y at the first or last base of the contig; reads cannot span such a position, so `Clipping` is NA |
| `Clipping_ge1` | Permissive call: >= 1 distinct clipped read at the exact position and no spanning read |
| `Clipping` | Final call: Y when >= 2 distinct clipped HiFi reads at the exact position and no read spans it; N otherwise; NA at a reference boundary |
| `High_confidence_clipping` | Y when `Clipping` is Y and the position is a clipping candidate |
| `Clipping_candidate` | anvi'o coverage >= 10 and clipping ratio = 1.0 |
| `Clipping_confidence` | HIGH_CONFIDENCE: `Clipping` Y and a clipping candidate; DEFAULT_ONLY: `Clipping` Y below the candidate thresholds; NOT_VALIDATED otherwise |
| `Anvio_position_0based` | Position as reported by anvi'o (0-based) |
| `Relative_position` | Position divided by the contig length (anvi'o) |
| `Coverage` | Read depth at the position (anvi'o) |
| `Clipping_count` | Clipped alignment records at the position (anvi'o) |
| `Clipping_ratio` | `Clipping_count` / `Coverage` (anvi'o) |
| `n_primary_alignment_records_covering_site` | Primary alignment records overlapping the position |
| `n_secondary_alignment_records_covering_site` | Secondary alignment records overlapping the position |
| `n_supplementary_alignment_records_covering_site` | Supplementary alignment records overlapping the position |
| `n_unique_all_overlapping_reads` | Distinct reads with any alignment overlapping the position |
| `n_unique_primary_clipped_reads` | Distinct reads whose primary alignment is clipped exactly at the position |
| `n_unique_secondary_clipped_reads` | Distinct reads whose secondary alignment is clipped exactly at the position |
| `n_unique_supplementary_clipped_reads` | Distinct reads whose supplementary alignment is clipped exactly at the position |
| `n_unique_all_clipped_reads` | Distinct reads clipped exactly at the position (any alignment class) |
| `Read_clipping_ratio` | `n_unique_all_clipped_reads` / `n_unique_all_overlapping_reads` |
| `n_unique_exact_spanning_reads` | Distinct reads with one alignment that has aligned bases on both sides of the position (counter-evidence) |
| `Exact_spanning_min_flank_bp_min` | Shorter aligned flank of the spanning reads: minimum (largest per read) |
| `Exact_spanning_min_flank_bp_median` | Shorter aligned flank of the spanning reads: median |
| `Exact_spanning_min_flank_bp_max` | Shorter aligned flank of the spanning reads: maximum |
| `n_unique_spanning_with_I_reads` | Spanning reads with an insertion next to the position |
| `Spanning_I_length_bp_min` | Min I length among those reads (largest per read) |
| `Spanning_I_length_bp_median` | Median I length among those reads (largest per read) |
| `Spanning_I_length_bp_max` | Max I length among those reads (largest per read) |
| `n_unique_spanning_with_D_reads` | Spanning reads with a deletion over the position |
| `Spanning_D_length_bp_min` | Min D length among those reads (largest per read) |
| `Spanning_D_length_bp_median` | Median D length among those reads (largest per read) |
| `Spanning_D_length_bp_max` | Max D length among those reads (largest per read) |
| `n_unique_spanning_with_N_reads` | Spanning reads with a reference skip over the position |
| `Spanning_N_length_bp_min` | Min N length among those reads (largest per read) |
| `Spanning_N_length_bp_median` | Median N length among those reads (largest per read) |
| `Spanning_N_length_bp_max` | Max N length among those reads (largest per read) |
| `n_unique_spanning_with_any_IDN_reads` | Spanning reads with any I, D or N event at the position |
| `n_unique_clip_and_span_reads` | Reads that are both clipped at and spanning the position (different alignments) |
| `Anvio_count_check` | PASS when the recounted coverage and clipping count equal the anvi'o values; MISMATCH otherwise; BOUNDARY_NOT_COMPARABLE at a reference boundary |
| `Mixed_clipping` | Y when the position has both clipped and spanning reads |
| `Exact_position_validation` | PASS: >= 2 clipped reads and no spanning read; INSUFFICIENT_CLIPPED_QNAMES: only one clipped read; FAIL: spanned or not clipped; BOUNDARY_UNRESOLVED at a reference boundary |

### `clipping.contig.tsv`

One row per contig of every assembly; counts and calls from `clipping.site.tsv`.

| Column | Definition |
|---|---|
| `Contig` | Contig ID (first word of the FASTA header) |
| `MAG_ID` | MAG ID (MAG FASTA file name without extension); `NA` for unbinned contigs |
| `Length_bp` | Contig length (bp) |
| `Internal_clipping_ge1` | Y when >= 1 INTERNAL site has `Clipping_ge1` Y |
| `Terminal_clipping_ge1` | Y when >= 1 TERMINAL site has `Clipping_ge1` Y |
| `Clipping_ge1` | Y when >= 1 site has `Clipping_ge1` Y |
| `Internal_clipping` | Y when >= 1 INTERNAL site has `Clipping` Y |
| `Terminal_clipping` | Y when >= 1 TERMINAL site has `Clipping` Y |
| `Clipping_analyzable` | N when there is no `Clipping` Y site and >= 1 candidate lies on a reference boundary; Y otherwise |
| `Clipping` | Y when >= 1 site has `Clipping` Y; N when none and every candidate was testable; NA otherwise |
| `High_confidence_clipping` | Y when >= 1 site has `High_confidence_clipping` Y; NA when `Clipping_analyzable` is N |
| `n_internal_clipping_ge1_sites` | INTERNAL sites with `Clipping_ge1` Y |
| `n_terminal_clipping_ge1_sites` | TERMINAL sites with `Clipping_ge1` Y |
| `n_clipping_ge1_sites` | Sites with `Clipping_ge1` Y |
| `n_internal_clipping_sites` | INTERNAL sites with `Clipping` Y |
| `n_terminal_clipping_sites` | TERMINAL sites with `Clipping` Y |
| `n_clipping_sites` | Sites with `Clipping` Y |
| `n_high_confidence_clipping_sites` | Sites with `High_confidence_clipping` Y |
| `n_internal_candidate_sites` | INTERNAL candidates |
| `n_terminal_candidate_sites` | TERMINAL candidates |
| `n_candidate_sites` | All candidates |
| `n_clipping_candidate_sites` | Clipping candidates |
| `n_candidates_with_spanning` | Candidates with >= 1 spanning read |
| `n_candidates_with_spanning_I` | Candidates with >= 1 spanning read carrying an insertion |
| `n_candidates_with_spanning_D` | Candidates with >= 1 spanning read carrying a deletion |
| `n_candidates_with_spanning_N` | Candidates with >= 1 spanning read carrying a reference skip |
| `n_candidates_with_spanning_any_IDN` | Candidates with >= 1 spanning read carrying I, D or N |
| `Max_anvio_coverage` | Highest anvi'o coverage among the candidates |
| `Max_anvio_clipping_count` | Highest anvi'o clipping count among the candidates |
| `Max_unique_clipped_hifi_reads` | Highest number of distinct clipped reads at one candidate |
| `Max_anvio_clipping_ratio` | Highest anvi'o clipping ratio among the candidates |
| `Candidate` | Y when >= 1 candidate |
| `Clipping_candidate` | Y when >= 1 clipping candidate |

### `clipping.MAG.tsv`

One row per MAG; counts and calls over the MAG contigs.

| Column | Definition |
|---|---|
| `MAG_ID` | MAG ID (MAG FASTA file name without extension); `NA` for unbinned contigs |
| `n_contigs` | Number of contigs in the MAG |
| `MAG_length_bp` | Sum of the MAG contig lengths (bp) |
| `Internal_clipping_ge1` | Y when >= 1 INTERNAL site has `Clipping_ge1` Y |
| `Terminal_clipping_ge1` | Y when >= 1 TERMINAL site has `Clipping_ge1` Y |
| `Clipping_ge1` | Y when >= 1 site has `Clipping_ge1` Y |
| `Internal_clipping` | Y when >= 1 INTERNAL site has `Clipping` Y |
| `Terminal_clipping` | Y when >= 1 TERMINAL site has `Clipping` Y |
| `Clipping_analyzable` | N when there is no `Clipping` Y site and >= 1 candidate lies on a reference boundary; Y otherwise |
| `Clipping` | Y when >= 1 site has `Clipping` Y; N when none and every candidate was testable; NA otherwise |
| `High_confidence_clipping` | Y when >= 1 site has `High_confidence_clipping` Y; NA when `Clipping_analyzable` is N |
| `n_internal_clipping_ge1_sites` | INTERNAL sites with `Clipping_ge1` Y |
| `n_terminal_clipping_ge1_sites` | TERMINAL sites with `Clipping_ge1` Y |
| `n_clipping_ge1_sites` | Sites with `Clipping_ge1` Y |
| `n_internal_clipping_sites` | INTERNAL sites with `Clipping` Y |
| `n_terminal_clipping_sites` | TERMINAL sites with `Clipping` Y |
| `n_clipping_sites` | Sites with `Clipping` Y |
| `n_high_confidence_clipping_sites` | Sites with `High_confidence_clipping` Y |
| `n_internal_candidate_sites` | INTERNAL candidates |
| `n_terminal_candidate_sites` | TERMINAL candidates |
| `n_candidate_sites` | All candidates |
| `n_clipping_candidate_sites` | Clipping candidates |
| `n_candidates_with_spanning` | Candidates with >= 1 spanning read |
| `n_candidates_with_spanning_I` | Candidates with >= 1 spanning read carrying an insertion |
| `n_candidates_with_spanning_D` | Candidates with >= 1 spanning read carrying a deletion |
| `n_candidates_with_spanning_N` | Candidates with >= 1 spanning read carrying a reference skip |
| `n_candidates_with_spanning_any_IDN` | Candidates with >= 1 spanning read carrying I, D or N |
| `n_contigs_with_candidate` | Contigs with >= 1 candidate |
| `n_contigs_with_clipping_candidate` | Contigs with >= 1 clipping candidate |
| `n_contigs_with_clipping` | Contigs with >= 1 `Clipping` Y site |
| `n_contigs_with_high_confidence_clipping` | Contigs with >= 1 `High_confidence_clipping` Y site |
| `Candidate_contigs` | Contigs with >= 1 candidate (`;`-separated) |
| `Clipping_candidate_contigs` | Contigs with >= 1 clipping candidate |
| `Clipping_contigs` | Contigs with >= 1 `Clipping` Y site |
| `High_confidence_clipping_contigs` | Contigs with >= 1 `High_confidence_clipping` Y site |
| `Max_anvio_coverage` | Highest anvi'o coverage among the candidates |
| `Max_anvio_clipping_count` | Highest anvi'o clipping count among the candidates |
| `Max_unique_clipped_hifi_reads` | Highest number of distinct clipped reads at one candidate |
| `Max_anvio_clipping_ratio` | Highest anvi'o clipping ratio among the candidates |
| `Candidate` | Y when >= 1 candidate |
| `Clipping_candidate` | Y when >= 1 clipping candidate |

### `clipping.assembler.tsv`

One row per assembler and one `ALL` row: site densities in INTERNAL and TERMINAL regions.

| Column | Definition |
|---|---|
| `Assembler` | Assembly file (assembler) the contig comes from |
| `Internal_analyzed_bp` | Bases in INTERNAL regions of all contigs |
| `Terminal_analyzed_bp` | Bases in TERMINAL regions (100 bp at each contig end) |
| `Internal_candidate_sites` | INTERNAL candidates |
| `Internal_candidate_sites_per_Mbp` | INTERNAL candidates per Mbp of INTERNAL bases |
| `Terminal_candidate_sites` | TERMINAL candidates |
| `Terminal_candidate_sites_per_Mbp` | TERMINAL candidates per Mbp of TERMINAL bases |
| `Terminal_to_internal_candidate_sites_density_ratio` | TERMINAL per-Mbp density / INTERNAL per-Mbp density (enrichment at contig ends) |
| `Internal_clipping_candidate_sites` | INTERNAL clipping candidates |
| `Internal_clipping_candidate_sites_per_Mbp` | INTERNAL clipping candidates per Mbp of INTERNAL bases |
| `Terminal_clipping_candidate_sites` | TERMINAL clipping candidates |
| `Terminal_clipping_candidate_sites_per_Mbp` | TERMINAL clipping candidates per Mbp of TERMINAL bases |
| `Terminal_to_internal_clipping_candidate_sites_density_ratio` | TERMINAL per-Mbp density / INTERNAL per-Mbp density (enrichment at contig ends) |
| `Internal_clipping_ge1_sites` | INTERNAL sites with `Clipping_ge1` Y |
| `Internal_clipping_ge1_sites_per_Mbp` | INTERNAL sites with `Clipping_ge1` Y per Mbp of INTERNAL bases |
| `Terminal_clipping_ge1_sites` | TERMINAL sites with `Clipping_ge1` Y |
| `Terminal_clipping_ge1_sites_per_Mbp` | TERMINAL sites with `Clipping_ge1` Y per Mbp of TERMINAL bases |
| `Terminal_to_internal_clipping_ge1_sites_density_ratio` | TERMINAL per-Mbp density / INTERNAL per-Mbp density (enrichment at contig ends) |
| `Internal_clipping_sites` | INTERNAL sites with `Clipping` Y |
| `Internal_clipping_sites_per_Mbp` | INTERNAL sites with `Clipping` Y per Mbp of INTERNAL bases |
| `Terminal_clipping_sites` | TERMINAL sites with `Clipping` Y |
| `Terminal_clipping_sites_per_Mbp` | TERMINAL sites with `Clipping` Y per Mbp of TERMINAL bases |
| `Terminal_to_internal_clipping_sites_density_ratio` | TERMINAL per-Mbp density / INTERNAL per-Mbp density (enrichment at contig ends) |
| `Internal_candidate_clipped_alignment_records` | INTERNAL clipped alignment records at candidates (anvi'o clipping counts) |
| `Internal_candidate_clipped_alignment_records_per_Mbp` | INTERNAL clipped alignment records at candidates (anvi'o clipping counts) per Mbp of INTERNAL bases |
| `Terminal_candidate_clipped_alignment_records` | TERMINAL clipped alignment records at candidates (anvi'o clipping counts) |
| `Terminal_candidate_clipped_alignment_records_per_Mbp` | TERMINAL clipped alignment records at candidates (anvi'o clipping counts) per Mbp of TERMINAL bases |
| `Terminal_to_internal_candidate_clipped_alignment_records_density_ratio` | TERMINAL per-Mbp density / INTERNAL per-Mbp density (enrichment at contig ends) |

## 5. False closure

**Directory:** `04_false_closure/` · **Command:** `false-closure`

### `false_closure.contig.tsv`

One row per contig reported circular by its assembler. A false closure is a contig reported circular whose sequence is not one continuous genome: at a validated clipping site (`Clipping` Y), >= 2 reads on **each** side are split between this contig and another contig (primary + supplementary alignments).

| Column | Definition |
|---|---|
| `Contig` | Contig ID (first word of the FASTA header) |
| `MAG_ID` | MAG ID (MAG FASTA file name without extension); `NA` for unbinned contigs |
| `Assembler` | Assembly file (assembler) the contig comes from |
| `Length_bp` | Contig length (bp) |
| `n_contigs` | Number of contigs in the MAG of this contig (0 for unbinned contigs) |
| `Assembler_reported_circularity` | Circularity reported by the assembler (Y for every row of the contig table) |
| `Circularity_resolution_method` | How the assembler call was read: the assembler circularity file, the FASTA header, the metaFlye assembly_info.txt, or the MAG-level call of a single-contig MAG |
| `False_closure_analyzable` | Y when `False_closure` is Y or N |
| `False_closure` | Y when >= 1 breakpoint is a false-closure breakpoint; N when none and every clipping candidate on the contig was testable; NA otherwise |
| `Validated_clipping_present` | Y when the contig has >= 1 `Clipping` Y site |
| `Unresolved_clipping_candidate_present` | Y when the contig has a clipping candidate on a reference boundary (cannot be tested) |
| `False_closure_breakpoints` | Number of false-closure breakpoints |
| `Internal_false_closure_breakpoints` | False-closure breakpoints in INTERNAL regions |
| `Terminal_false_closure_breakpoints` | False-closure breakpoints in TERMINAL regions |
| `False_closure_positions` | False-closure breakpoints as `REGION:position` |
| `Max_left_cross_contig_connected_reads` | Highest number of reads clipped on the left side of a breakpoint whose other part maps to another contig |
| `Max_right_cross_contig_connected_reads` | Same for the right side |
| `Partner_contigs` | Other contigs joined to the false-closure breakpoints by split reads |
| `Internal_false_closure_evidence` | Y when >= 1 INTERNAL false-closure breakpoint |
| `Terminal_false_closure_evidence` | Y when >= 1 TERMINAL false-closure breakpoint |
| `False_closure_evidence` | Y when >= 1 false-closure breakpoint |
| `Ribosomal_protein_family_count` | Ribosomal protein families annotated by `prep-scg` on the contig (MAG table: the highest among its circular contigs) |
| `Ribosomal_protein_families` | Those families |
| `Small_genome_length_criterion` | Y when the contig (MAG table: any circular contig) is < 500,000 bp |
| `Small_genome_ribosomal_criterion` | Y when the contig (MAG table: any circular contig) has >= 3 ribosomal protein families |
| `Small_genome_proxy` | Y when one circular contig meets both small-genome criteria (a small complete genome rather than a fragment) |
| `Premature_circularization_candidate` | Contig table: Y when `False_closure` is Y. MAG table: Y for a single-contig MAG reported circular with false-closure evidence |
| `Premature_circularization_screen` | Y when `Premature_circularization_candidate` or `Small_genome_proxy` is Y |
| `Validated_clipping_breakpoints` | Number of `Clipping` Y sites tested on the contig |
| `Breakpoint_details` | Every tested breakpoint, `;`-separated: `REGION@position(distance to end)\|L=connected/clipped reads[target contig:reads]\|R=...\|cov=,clip=,frac=\|spanning=\|evidence=Y/N\|reason` |

### `false_closure.MAG.tsv`

One row per MAG. A false closure is a contig reported circular whose sequence is not one continuous genome: at a validated clipping site (`Clipping` Y), >= 2 reads on **each** side are split between this contig and another contig (primary + supplementary alignments).

| Column | Definition |
|---|---|
| `MAG_ID` | MAG ID (MAG FASTA file name without extension); `NA` for unbinned contigs |
| `n_contigs` | Number of contigs in the MAG |
| `False_closure_analyzable` | Y when `False_closure` is Y or N |
| `Assembler_reported_circularity` | Assembler circularity of the MAG (`prep-circularity` table) |
| `n_circular_contigs` | MAG contigs reported circular by the assembler |
| `Circular_contigs` | Those contigs |
| `Contig_circularity_complete` | Y when the assembler circularity of every MAG contig is known |
| `Ribosomal_protein_family_count` | Ribosomal protein families annotated by `prep-scg` on the contig (MAG table: the highest among its circular contigs) |
| `Ribosomal_protein_families` | Those families |
| `Small_genome_length_criterion` | Y when the contig (MAG table: any circular contig) is < 500,000 bp |
| `Small_genome_ribosomal_criterion` | Y when the contig (MAG table: any circular contig) has >= 3 ribosomal protein families |
| `Small_genome_proxy` | Y when one circular contig meets both small-genome criteria (a small complete genome rather than a fragment) |
| `Small_genome_proxy_contigs` | Circular contigs that meet both small-genome criteria |
| `Premature_circularization_candidate` | Contig table: Y when `False_closure` is Y. MAG table: Y for a single-contig MAG reported circular with false-closure evidence |
| `Premature_circularization_screen` | Y when `Premature_circularization_candidate` or `Small_genome_proxy` is Y |
| `Large_insertion_sensitivity_analyzable` | Y when `Large_insertion_sensitivity_support` is Y or N |
| `Large_insertion_sensitivity_support` | Sensitivity screen: Y when a circular contig has a position spanned by >= 2 reads that each carry an insertion >= 100 bp |
| `Large_insertion_sensitivity_site_count` | Number of such positions |
| `Large_insertion_sensitivity_sites` | Those positions as `contig\|REGION:position[n=reads,minI=shortest insertion]` |
| `Partner_contigs` | Other contigs joined to the false-closure breakpoints by split reads |
| `Partner_MAGs` | MAGs of the partner contigs |
| `Validated_clipping_present` | Y when a circular MAG contig has >= 1 `Clipping` Y site |
| `Internal_validated_clipping_present` | Y when a circular contig has an INTERNAL `Clipping` Y site |
| `Terminal_validated_clipping_present` | Y when a circular contig has a TERMINAL `Clipping` Y site |
| `Internal_false_closure_evidence` | Y when >= 1 INTERNAL false-closure breakpoint |
| `Terminal_false_closure_evidence` | Y when >= 1 TERMINAL false-closure breakpoint |
| `False_closure_evidence` | Y when >= 1 false-closure breakpoint |
| `False_closure` | Y when a circular MAG contig has >= 1 false-closure breakpoint; N when none and the MAG has circular contigs whose candidates were all testable; NA otherwise |
| `False_closure_breakpoints` | False-closure breakpoints as `contig\|REGION:position` |
| `False_closure_regions` | Region types (INTERNAL, TERMINAL) of the false-closure breakpoints |

## 6. Phantom sequences

**Directory:** `05_phantom/` · **Command:** `phantom`

### `phantom.region.tsv`

One row per phantom region. A phantom region is a stretch longer than 500 bp with zero depth of all HiFi alignments (secondary included) of its own assembler, no BLASTN hit to the raw reads, and no ambiguous base.

| Column | Definition |
|---|---|
| `Contig` | Contig ID (first word of the FASTA header) |
| `Region_id` | Region ID |
| `Start_0based` | Region start (0-based, BED) |
| `End_0based` | Region end (0-based, exclusive, BED) |
| `Region_length_bp` | Region length (bp) |

### `phantom.contig.tsv`

One row per contig of every assembly. A phantom region is a stretch longer than 500 bp with zero depth of all HiFi alignments (secondary included) of its own assembler, no BLASTN hit to the raw reads, and no ambiguous base.

| Column | Definition |
|---|---|
| `Contig` | Contig ID (first word of the FASTA header) |
| `MAG_ID` | MAG ID (MAG FASTA file name without extension); `NA` for unbinned contigs |
| `Length_bp` | Contig length (bp) |
| `n_phantom_regions` | Phantom regions |
| `Phantom_bp` | Total phantom bases |
| `Phantom_fraction` | `Phantom_bp` / length |
| `Max_phantom_region_bp` | Longest phantom region (bp) |
| `Phantom_regions_0based` | Regions as `start-end` (0-based) |
| `Phantom_region_ids` | Region IDs |
| `Phantom` | Y when >= 1 phantom region |

### `phantom.MAG.tsv`

One row per MAG. A phantom region is a stretch longer than 500 bp with zero depth of all HiFi alignments (secondary included) of its own assembler, no BLASTN hit to the raw reads, and no ambiguous base.

| Column | Definition |
|---|---|
| `MAG_ID` | MAG ID (MAG FASTA file name without extension); `NA` for unbinned contigs |
| `n_contigs` | Number of contigs in the MAG |
| `MAG_length_bp` | Sum of the MAG contig lengths (bp) |
| `n_phantom_regions` | Phantom regions |
| `Phantom_bp` | Total phantom bases |
| `Phantom_fraction` | `Phantom_bp` / length |
| `Phantom_contigs` | MAG contigs with phantom regions |
| `Phantom_region_ids` | Region IDs |
| `Candidate_distinct_kmers` | Distinct 21-mers in the phantom regions of the MAG (Meryl) |
| `Absent_from_reads_distinct_kmers` | Of these, 21-mers absent from the HiFi reads |
| `Absent_kmer_fraction` | `Absent_from_reads_distinct_kmers` / `Candidate_distinct_kmers` (NA without phantom regions); a measurement, not a call |

## 7. Repeats

**Directory:** `06_repeats/` · **Command:** `repeats`

### `repeats.contig.tsv`

One row per contig of every assembly. Repeats are BLASTN alignments longer than 200 bp with >= 80% identity: within a contig (INTRA, self-alignment off the diagonal) or between contigs of the same MAG (INTER). Every contig is scanned for INTRA repeats; INTER repeats exist only for multi-contig MAGs. Contig-level repeat breadth uses INTRA only, while the MAG total is the union of INTRA and INTER with overlaps counted once.

| Column | Definition |
|---|---|
| `Contig` | Contig ID (first word of the FASTA header) |
| `MAG_ID` | MAG ID (MAG FASTA file name without extension); `NA` for unbinned contigs |
| `Assembler` | Assembly file (assembler) the contig comes from |
| `MAG_structure` | single_contig, multi_contig, or unbinned |
| `Length_bp` | Contig length (bp) |
| `Intra_repeat` | Y when the contig (MAG: any contig) has INTRA repeat bases |
| `Intra_repeat_bp` | Bases covered by INTRA repeats |
| `Intra_repeat_fraction` | `Intra_repeat_bp` / length |
| `n_intra_repeat_hsps` | INTRA repeat alignments (HSPs) |
| `Excessive_repeat_contig` | Y when `Intra_repeat_fraction` >= 0.70 |

### `repeats.MAG.tsv`

One row per MAG. Repeats are BLASTN alignments longer than 200 bp with >= 80% identity: within a contig (INTRA, self-alignment off the diagonal) or between contigs of the same MAG (INTER). Every contig is scanned for INTRA repeats; INTER repeats exist only for multi-contig MAGs. Contig-level repeat breadth uses INTRA only, while the MAG total is the union of INTRA and INTER with overlaps counted once.

| Column | Definition |
|---|---|
| `MAG_ID` | MAG ID (MAG FASTA file name without extension); `NA` for unbinned contigs |
| `Assembler` | Assembly file (assembler) the contig comes from |
| `MAG_structure` | single_contig, multi_contig, or unbinned |
| `n_contigs` | Number of contigs in the MAG |
| `MAG_length_bp` | Sum of the MAG contig lengths (bp) |
| `Intra_repeat` | Y when the contig (MAG: any contig) has INTRA repeat bases |
| `Intra_repeat_bp` | Bases covered by INTRA repeats |
| `Intra_repeat_fraction` | `Intra_repeat_bp` / length |
| `Inter_repeat` | Y when INTER repeat bases exist |
| `Inter_repeat_bp` | Bases covered by INTER repeats |
| `Inter_repeat_fraction` | `Inter_repeat_bp` / MAG length |
| `Total_repeat` | Y when INTRA or INTER repeat bases exist |
| `Total_repeat_bp` | Bases covered by INTRA or INTER repeats (union) |
| `Total_repeat_fraction` | `Total_repeat_bp` / MAG length |
| `n_excessive_repeat_contigs` | Contigs with `Excessive_repeat_contig` Y |
| `Has_excessive_intra_repeat_contig` | Y when `n_excessive_repeat_contigs` > 0 |
| `Excessive_repeat_MAG` | Y when `Total_repeat_fraction` >= 0.70 |

### `repeats.region.tsv`

One row per merged repeat region.

| Column | Definition |
|---|---|
| `Contig` | Contig ID (first word of the FASTA header) |
| `MAG_ID` | MAG ID (MAG FASTA file name without extension); `NA` for unbinned contigs |
| `MAG_structure` | single_contig, multi_contig, or unbinned |
| `Repeat_type` | INTRA or INTER |
| `Repeat_region_id` | Region ID |
| `Start_1based` | Region start (1-based) |
| `End_1based` | Region end (1-based, inclusive) |
| `Repeat_region_bp` | Region length (bp) |

### `repeats.HSP.tsv`

One row per repeat alignment (HSP).

| Column | Definition |
|---|---|
| `MAG_ID` | MAG ID (MAG FASTA file name without extension); `NA` for unbinned contigs |
| `Repeat_type` | INTRA or INTER |
| `Query_contig` | Contig of the query copy |
| `Subject_contig` | Contig of the subject copy |
| `Identity_percent` | Alignment identity (%) |
| `Alignment_length_bp` | Alignment length (bp) |
| `Query_start` | Query start (1-based) |
| `Query_end` | Query end |
| `Subject_start` | Subject start (1-based) |
| `Subject_end` | Subject end |
| `Query_length_bp` | Query contig length |
| `Subject_length_bp` | Subject contig length |
| `Orientation` | direct or inverted |

## 8. False duplication

**Directory:** `08_false_duplication/` · **Command:** `false-duplication`

### `false_duplication.pair.tsv`

One row per duplicated pair. Duplicated pairs come from a minimap2 self-alignment of each assembly (`-x asm5 -DP`), keeping copies >= 1,000 bp with >= 95% identity inside a MAG. Each copy boundary is tested with HiFi reads (primary + supplementary, 500-bp anchors): supported (Y) when >= 2 reads span it, N when >= 2 reads reach each side but none spans it, NA otherwise. A false duplication is a pair where one copy is spanned at both boundaries and the other at neither (for tandem copies: both outer boundaries spanned and the junction between the copies not).

| Column | Definition |
|---|---|
| `MAG_ID` | MAG ID (MAG FASTA file name without extension); `NA` for unbinned contigs |
| `Assembler` | Assembly file (assembler) the contig comes from |
| `Pair_id` | Pair ID |
| `Pair_scope` | INTRA_CONTIG (both copies on one contig) or INTER_CONTIG |
| `Duplication_type` | TANDEM (adjacent copies), INTERSPERSED, or OVERLAPPING |
| `Orientation` | DIRECT or INVERTED |
| `Copy1_contig` | Contig of copy 1 |
| `Copy1_contig_length` | Length of that contig |
| `Copy1_start` | Copy 1 start (1-based) |
| `Copy1_end` | Copy 1 end |
| `Copy2_contig` | Contig of copy 2 |
| `Copy2_contig_length` | Length of that contig |
| `Copy2_start` | Copy 2 start (1-based) |
| `Copy2_end` | Copy 2 end |
| `Identity_percent` | BLASTN identity (%) |
| `Alignment_length` | BLASTN alignment length (bp) |
| `Candidate_cluster_size` | Near-identical alignments merged into this pair |
| `Copy1_left_supported` | Left boundary of copy 1: Y/N/NA as above |
| `Copy1_right_supported` | Right boundary of copy 1 |
| `Copy2_left_supported` | Left boundary of copy 2 |
| `Copy2_right_supported` | Right boundary of copy 2 |
| `Assessment` | Call and reason: Y_false_duplication_pattern, N_no_false_duplication_pattern, N_tandem_junction_spanned, NA_insufficient_reads, NA_inverted_tandem, NA_overlapping_copies |
| `Assessment_complete` | Y when the call is Y or N |
| `False_duplication` | Pair call: Y, N or NA |

### `false_duplication.contig.tsv`

One row per MAG contig.

| Column | Definition |
|---|---|
| `Contig` | Contig ID (first word of the FASTA header) |
| `MAG_ID` | MAG ID (MAG FASTA file name without extension); `NA` for unbinned contigs |
| `Assembler` | Assembly file (assembler) the contig comes from |
| `Length_bp` | Contig length (bp) |
| `n_duplicated_pairs` | Duplicated pairs involving the contig (MAG: in the MAG) |
| `n_intra_contig_pairs` | INTRA_CONTIG pairs |
| `n_inter_contig_pairs` | INTER_CONTIG pairs |
| `n_tandem_pairs` | TANDEM pairs |
| `n_direct_tandem_pairs` | DIRECT TANDEM pairs |
| `n_inverted_tandem_pairs` | INVERTED TANDEM pairs |
| `n_interspersed_pairs` | INTERSPERSED pairs |
| `n_direct_pairs` | DIRECT pairs |
| `n_inverted_pairs` | INVERTED pairs |
| `Assessment_complete` | Y when every pair has a Y or N call (also Y without pairs) |
| `Analyzable` | Y when `False_duplication` is Y or N |
| `False_duplication` | Y when >= 1 pair is Y; NA when >= 1 pair is NA and none Y; N otherwise |
| `n_analyzable_pairs` | Pairs with a Y or N call |
| `n_false_duplication_pairs` | Pairs with call Y |
| `n_no_false_duplication_pairs` | Pairs with call N |
| `n_not_analyzable_pairs` | Pairs with call NA |

### `false_duplication.MAG.tsv`

One row per MAG.

| Column | Definition |
|---|---|
| `MAG_ID` | MAG ID (MAG FASTA file name without extension); `NA` for unbinned contigs |
| `n_contigs` | Number of contigs in the MAG |
| `n_duplicated_pairs` | Duplicated pairs involving the contig (MAG: in the MAG) |
| `n_intra_contig_pairs` | INTRA_CONTIG pairs |
| `n_inter_contig_pairs` | INTER_CONTIG pairs |
| `n_tandem_pairs` | TANDEM pairs |
| `n_direct_tandem_pairs` | DIRECT TANDEM pairs |
| `n_inverted_tandem_pairs` | INVERTED TANDEM pairs |
| `n_interspersed_pairs` | INTERSPERSED pairs |
| `n_false_duplication_pairs` | Pairs with call Y |
| `n_no_false_duplication_pairs` | Pairs with call N |
| `n_not_analyzable_pairs` | Pairs with call NA |
| `Assessment_complete` | Y when every pair has a Y or N call (also Y without pairs) |
| `Analyzable` | Y when `False_duplication` is Y or N |
| `False_duplication` | Y when >= 1 pair is Y; NA when >= 1 pair is NA and none Y; N otherwise |

## 9. Assembly variants

**Directory:** `09_assembly_variants/` · **Command:** `variants`

### `variants.site.tsv`

One row per minority or unsupported assembly variant. anvi'o profiles the HiFi reads on each assembly and reports variable positions and indels. A site is a **minority** assembly variant when the assembly allele is supported by some reads but another allele (SNV) or the indel is carried by > 50% of reads, and an **unsupported** assembly variant when no read supports the assembly allele.

| Column | Definition |
|---|---|
| `Contig` | Contig ID (first word of the FASTA header) |
| `Assembler` | Assembly file (assembler) the contig comes from |
| `Position` | Variant position (1-based) |
| `Variant_type` | SNV or INDEL |
| `Assembly_allele` | Base (or sequence) in the assembly |
| `Dominant_alternative_allele` | Most frequent other allele in the reads |
| `Assembly_support_count` | Reads supporting the assembly allele |
| `Alternative_support_count` | Reads supporting the alternative allele |
| `Coverage` | Read depth at the position |
| `Assembly_support_fraction` | `Assembly_support_count` / `Coverage` |
| `Alternative_support_fraction` | `Alternative_support_count` / `Coverage` |
| `Assembly_support_class` | Minority or Unsupported |
| `Corresponding_gene_call` | anvi'o gene call at the position (NA outside genes) |
| `Anvio_table` | anvi'o table the row comes from (variable nucleotides or indels) |

### `variants.contig.tsv`

One row per MAG contig. anvi'o profiles the HiFi reads on each assembly and reports variable positions and indels. A site is a **minority** assembly variant when the assembly allele is supported by some reads but another allele (SNV) or the indel is carried by > 50% of reads, and an **unsupported** assembly variant when no read supports the assembly allele.

| Column | Definition |
|---|---|
| `Contig` | Contig ID (first word of the FASTA header) |
| `MAG_ID` | MAG ID (MAG FASTA file name without extension); `NA` for unbinned contigs |
| `Length_bp` | Contig length (bp) |
| `Assembly_variant` | Y when >= 1 minority or unsupported variant |
| `Anvio_variable_nucleotide_rows_examined` | anvi'o variable-nucleotide rows read |
| `Anvio_indel_rows_examined` | anvi'o indel rows read |
| `Minority_assembly_snv_count` | Minority SNVs |
| `Unsupported_assembly_snv_count` | Unsupported SNVs |
| `Minority_assembly_indel_count` | Minority indels |
| `Unsupported_assembly_indel_count` | Unsupported indels |
| `Minority_impacted_gene_count` | Genes with >= 1 minority variant |
| `Unsupported_impacted_gene_count` | Genes with >= 1 unsupported variant |
| `Assembly_variant_count` | All minority and unsupported variants |
| `Assembly_variant_class` | Unsupported when >= 1 unsupported variant, otherwise Minority when >= 1 minority variant, otherwise None |

### `variants.MAG.tsv`

One row per MAG.

| Column | Definition |
|---|---|
| `MAG_ID` | MAG ID (MAG FASTA file name without extension); `NA` for unbinned contigs |
| `Assembly_variant` | Y when >= 1 minority or unsupported variant |
| `Minority_assembly_variant` | Y when >= 1 minority variant |
| `Unsupported_assembly_variant` | Y when >= 1 unsupported variant |
| `n_contigs` | Number of contigs in the MAG |
| `MAG_length_bp` | Sum of the MAG contig lengths (bp) |
| `n_contigs_with_assembly_allele_disagreement` | Contigs with >= 1 variant |
| `n_contigs_classified_as_minority_assembly_variant` | Contigs whose class is Minority |
| `n_contigs_classified_as_unsupported_assembly_variant` | Contigs whose class is Unsupported |
| `Minority_assembly_snv_count` | Minority SNVs |
| `Minority_assembly_snv_per_Mbp` | `Minority_assembly_snv_count` per Mbp of MAG |
| `Unsupported_assembly_snv_count` | Unsupported SNVs |
| `Unsupported_assembly_snv_per_Mbp` | `Unsupported_assembly_snv_count` per Mbp of MAG |
| `Minority_assembly_indel_count` | Minority indels |
| `Minority_assembly_indel_per_Mbp` | `Minority_assembly_indel_count` per Mbp of MAG |
| `Unsupported_assembly_indel_count` | Unsupported indels |
| `Unsupported_assembly_indel_per_Mbp` | `Unsupported_assembly_indel_count` per Mbp of MAG |
| `Minority_impacted_gene_count` | Genes with >= 1 minority variant |
| `Unsupported_impacted_gene_count` | Genes with >= 1 unsupported variant |
| `Assembly_variant_count` | All minority and unsupported variants |
| `Assembly_variant_per_Mbp` | `Assembly_variant_count` per Mbp of MAG |
| `Assembly_variant_class` | Unsupported when >= 1 unsupported variant, otherwise Minority when >= 1 minority variant, otherwise None |

## 10. Taxonomic chimera

**Directory:** `07_chimera/` · **Command:** `chimera`

### `chimera.region.tsv`

One row per taxonomic transition. Genes are classified against GTDB R232 (MMseqs2). At each gene boundary, up to 10 predicted genes are searched on each side and the nearest 5 classified genes form a block (>= 4 required). A taxonomic transition is a boundary where both blocks have >= 80% dominance at one rank (domain, phylum, class or order) and their dominant taxa differ. A transition becomes a taxonomic chimera when a validated clipping site (`Clipping` Y) lies in its interval (+/- 1 bp), and a high-confidence chimera when that site is also `High_confidence_clipping` Y. The HiFi structural check (primary + supplementary, 500-bp anchors, >= 2 reads) is reported alongside.

| Column | Definition |
|---|---|
| `Contig` | Contig ID (first word of the FASTA header) |
| `MAG_ID` | MAG ID (MAG FASTA file name without extension); `NA` for unbinned contigs |
| `Assembler` | Assembly file (assembler) the contig comes from |
| `Structural_assessment` | HiFi structural class: HIGH_CONFIDENCE_STRUCTURAL_CHIMERA, STRUCTURAL_CHIMERA, HIFI_CONTINUITY_SUPPORTED, UNRESOLVED_STRUCTURAL_CHIMERA |
| `Taxonomic_transition` | Y for every row |
| `Shift_rank` | Highest rank with a transition (domain > phylum > class > order) |
| `Analysis_position_1based` | Representative boundary position |
| `Transition_interval_start_1based` | Start of the interval of boundaries with the same transition |
| `Transition_interval_end_1based` | End of that interval |
| `n_supporting_default_clipping_positions` | anvi'o clipping candidates in the interval |
| `Supporting_default_clipping_positions_1based` | Their positions |
| `Max_default_clipping_coverage` | Highest anvi'o coverage among them |
| `Max_default_clipping_count` | Highest anvi'o clipping count |
| `Max_default_clipping_ratio` | Highest anvi'o clipping ratio |
| `Validated_clipping_support` | Y when a `Clipping` Y site lies in the interval (+/- 1 bp) |
| `n_supporting_validated_clipping_positions` | Number of such sites |
| `Supporting_validated_clipping_positions_1based` | Their positions |
| `MAG_GTDB_taxon_at_shift_rank` | GTDB-Tk taxon of the MAG at `Shift_rank` |
| `Left_matches_MAG_at_shift_rank` | Y when the left dominant taxon equals the MAG taxon |
| `Right_matches_MAG_at_shift_rank` | Y when the right dominant taxon equals the MAG taxon |
| `MAG_discordant_side_at_shift_rank` | Side whose taxon differs from the MAG: LEFT, RIGHT, BOTH or NA |
| `Left_dominant_taxon` | Most frequent taxon of the left block |
| `Right_dominant_taxon` | Most frequent taxon of the right block |
| `Left_block_genes` | Classified genes in the left block |
| `Right_block_genes` | Classified genes in the right block |
| `Left_classified_genes_at_rank` | Block genes classified at `Shift_rank` |
| `Right_classified_genes_at_rank` | Block genes classified at `Shift_rank` |
| `Left_searched_predicted_genes` | Predicted genes searched on the left |
| `Right_searched_predicted_genes` | Predicted genes searched on the right |
| `Left_unclassified_genes_at_rank` | Searched genes unclassified at `Shift_rank` |
| `Right_unclassified_genes_at_rank` | Searched genes unclassified at `Shift_rank` |
| `Left_gene_ids` | Gene IDs of the left block |
| `Right_gene_ids` | Gene IDs of the right block |
| `Left_searched_gene_ids` | Gene IDs searched on the left |
| `Right_searched_gene_ids` | Gene IDs searched on the right |
| `Left_dominant_count` | Block genes assigned to `Left_dominant_taxon` |
| `Right_dominant_count` | Block genes assigned to `Right_dominant_taxon` |
| `Left_dominance` | `Left_dominant_count` / `Left_block_genes` |
| `Right_dominance` | `Right_dominant_count` / `Right_block_genes` |
| `Left_block_start_1based` | Start of the left block |
| `Left_block_end_1based` | End of the left block |
| `Left_block_span_bp` | Length of the left block (bp) |
| `Right_block_start_1based` | Start of the right block |
| `Right_block_end_1based` | End of the right block |
| `Right_block_span_bp` | Length of the right block (bp) |
| `Left_GC_percent` | GC content of the left block (%) |
| `Right_GC_percent` | GC content of the right block (%) |
| `Delta_GC_percentage_points` | \|`Left_GC_percent` - `Right_GC_percent`\| |
| `Left_GC_callable_bases` | A/C/G/T bases used for `Left_GC_percent` |
| `Right_GC_callable_bases` | A/C/G/T bases used for `Right_GC_percent` |
| `Left_mean_depth` | Mean HiFi depth (primary alignments, `samtools bedcov`) over the left block |
| `Right_mean_depth` | Mean HiFi depth (primary alignments, `samtools bedcov`) over the right block |
| `Depth_right_over_left` | `Right_mean_depth` / `Left_mean_depth` |
| `Abs_log2_depth_ratio` | \|log2(`Depth_right_over_left`)\| |
| `Junction_left_boundary_1based` | Left boundary of the junction tested with HiFi reads |
| `Junction_right_boundary_1based` | Right boundary of that junction |
| `HiFi_left_anchor_reads` | Reads with >= 500 bp aligned left of the junction |
| `HiFi_right_anchor_reads` | Reads with >= 500 bp aligned right of the junction |
| `HiFi_spanning_reads` | Reads spanning the junction with 500-bp anchors on both sides |
| `HiFi_breakpoint_clipped_reads` | Reads clipped at the junction |
| `High_confidence_clipping_support` | Y when a supporting site is `High_confidence_clipping` Y |
| `Structural_discontinuity` | N when >= 2 reads span the junction; Y when >= 2 reads anchor on each side and >= 2 clipped reads (or a validated clipping site) exist; NA otherwise |
| `High_confidence_structural_discontinuity` | Y when `Structural_discontinuity` is Y and a supporting site is `High_confidence_clipping` Y |
| `SCG_domain` | Domain of the MAG (GTDB-Tk, `prep-taxonomy`) |
| `SCG_HMM_collection` | anvi'o SCG collection used (Bacteria_71 or Archaea_76) |
| `n_SCG_hits` | SCG hits on the contig |
| `n_distinct_SCG_families` | Distinct SCG families |
| `n_duplicated_SCG_families` | SCG families found more than once |
| `SCG_redundancy_fraction` | `n_duplicated_SCG_families` / `n_distinct_SCG_families` |
| `SCG_redundancy_observed` | Y when >= 1 SCG family is duplicated |
| `Duplicated_SCG_families` | Duplicated families as `family x count` |
| `Clipping_supported_transition` | Y when `Validated_clipping_support` is Y |
| `High_confidence_taxonomic_chimera` | Y when `High_confidence_clipping_support` is Y |
| `Taxonomic_assessment` | HIGH_CONFIDENCE_TAXONOMIC_CHIMERA, TAXONOMIC_CHIMERA or TAXONOMIC_SHIFT_ONLY (transition without validated clipping) |

### `chimera.contig.tsv`

One row per contig of every assembly. Genes are classified against GTDB R232 (MMseqs2). At each gene boundary, up to 10 predicted genes are searched on each side and the nearest 5 classified genes form a block (>= 4 required). A taxonomic transition is a boundary where both blocks have >= 80% dominance at one rank (domain, phylum, class or order) and their dominant taxa differ. A transition becomes a taxonomic chimera when a validated clipping site (`Clipping` Y) lies in its interval (+/- 1 bp), and a high-confidence chimera when that site is also `High_confidence_clipping` Y. The HiFi structural check (primary + supplementary, 500-bp anchors, >= 2 reads) is reported alongside.

| Column | Definition |
|---|---|
| `Contig` | Contig ID (first word of the FASTA header) |
| `MAG_ID` | MAG ID (MAG FASTA file name without extension); `NA` for unbinned contigs |
| `Assembler` | Assembly file (assembler) the contig comes from |
| `Length_bp` | Contig length (bp) |
| `n_predicted_genes` | Predicted genes (pyrodigal) |
| `n_taxonomy_classified_genes` | Genes with a GTDB classification |
| `n_gene_count_eligible_boundaries` | Gene boundaries with >= 4 genes on each side |
| `n_taxonomically_evaluable_boundaries` | Boundaries where both blocks could be built |
| `Clipping_supported_transition_assessable` | Y when >= 1 boundary is evaluable |
| `n_default_clipping_candidates` | anvi'o clipping candidates on the contig |
| `n_validated_clipping_sites` | `Clipping` Y sites on the contig |
| `Taxonomic_shift_ranks` | Ranks of the transitions |
| `Taxonomic_shift_positions_1based` | Positions of the transitions |
| `Max_delta_GC_percentage_points` | Largest GC difference across a transition |
| `Max_abs_log2_depth_ratio` | Largest \|log2 depth ratio\| across a transition |
| `Taxonomic_shift` | Y when >= 1 transition; N when assessable without one; NA when not assessable |
| `High_confidence_taxonomic_chimera` | Y when >= 1 region is `High_confidence_taxonomic_chimera` Y |
| `Taxonomic_assessment` | HIGH_CONFIDENCE_TAXONOMIC_CHIMERA, TAXONOMIC_CHIMERA, TAXONOMIC_SHIFT_ONLY, NO_TAXONOMIC_SHIFT, or NA |
| `Clipping_supported_transition` | Y when >= 1 region is `Clipping_supported_transition` Y; N when assessable without one; NA when not assessable |
| `n_structural_discontinuity_events` | Regions with `Structural_discontinuity` Y |
| `n_high_confidence_structural_discontinuity_events` | Regions with `High_confidence_structural_discontinuity` Y |
| `n_hifi_continuity_supported_events` | Regions with `Structural_discontinuity` N (spanned by reads) |
| `n_unresolved_structural_events` | Regions with `Structural_discontinuity` NA |
| `Structural_discontinuity_positions_1based` | Positions of `Structural_discontinuity` Y regions |
| `High_confidence_structural_discontinuity_positions_1based` | Positions of high-confidence structural regions |
| `Unresolved_structural_positions_1based` | Positions of unresolved regions |
| `Structural_discontinuity` | Y when >= 1 region is Y; NA when >= 1 is NA and none Y; N otherwise |
| `High_confidence_structural_discontinuity` | Y when >= 1 region is `High_confidence_structural_discontinuity` Y |
| `SCG_domain` | Domain of the MAG (GTDB-Tk, `prep-taxonomy`) |
| `SCG_HMM_collection` | anvi'o SCG collection used (Bacteria_71 or Archaea_76) |
| `n_SCG_hits` | SCG hits on the contig |
| `n_distinct_SCG_families` | Distinct SCG families |
| `n_duplicated_SCG_families` | SCG families found more than once |
| `SCG_redundancy_fraction` | `n_duplicated_SCG_families` / `n_distinct_SCG_families` |
| `SCG_redundancy_observed` | Y when >= 1 SCG family is duplicated |
| `Duplicated_SCG_families` | Duplicated families as `family x count` |
| `n_taxonomic_shift_regions` | Transition regions |
| `n_clipping_supported_transition_regions` | Regions with `Clipping_supported_transition` Y |
| `n_high_confidence_taxonomic_chimera_regions` | Regions with `High_confidence_taxonomic_chimera` Y |
| `n_taxonomic_shift_only_regions` | Regions without validated clipping |
| `Clipping_supported_transition_positions_1based` | Positions of `Clipping_supported_transition` Y regions |
| `High_confidence_taxonomic_chimera_positions_1based` | Positions of high-confidence regions |
| `Taxonomic_shift_only_positions_1based` | Positions of regions without validated clipping |

### `chimera.MAG.tsv`

One row per MAG; counts and calls over its contigs.

| Column | Definition |
|---|---|
| `MAG_ID` | MAG ID (MAG FASTA file name without extension); `NA` for unbinned contigs |
| `n_contigs` | Number of contigs in the MAG |
| `MAG_length_bp` | Sum of the MAG contig lengths (bp) |
| `n_taxonomic_shift_contigs` | Contigs with `Taxonomic_shift` Y |
| `n_structural_discontinuity_contigs` | Contigs with `Structural_discontinuity` Y |
| `n_high_confidence_structural_discontinuity_contigs` | Contigs with `High_confidence_structural_discontinuity` Y |
| `n_unresolved_structural_contigs` | Contigs with `Structural_discontinuity` NA |
| `n_structural_discontinuity_events` | Regions with `Structural_discontinuity` Y |
| `n_high_confidence_structural_discontinuity_events` | Regions with `High_confidence_structural_discontinuity` Y |
| `n_hifi_continuity_supported_events` | Regions with `Structural_discontinuity` N (spanned by reads) |
| `n_unresolved_structural_events` | Regions with `Structural_discontinuity` NA |
| `n_SCG_redundancy_contigs` | Contigs with `SCG_redundancy_observed` Y |
| `Structural_discontinuity_contig_ids` | Contigs with `Structural_discontinuity` Y |
| `High_confidence_structural_discontinuity_contig_ids` | Contigs with `High_confidence_structural_discontinuity` Y |
| `Unresolved_structural_contig_ids` | Contigs with `Structural_discontinuity` NA |
| `SCG_redundancy_contig_ids` | Contigs with duplicated SCG families |
| `Taxonomic_shift` | Y when >= 1 contig has `Taxonomic_shift` Y |
| `Structural_discontinuity` | Y when >= 1 contig has `Structural_discontinuity` Y |
| `Clipping_supported_transition` | Y when >= 1 contig has `Clipping_supported_transition` Y |
| `High_confidence_structural_discontinuity` | Y when >= 1 contig has `High_confidence_structural_discontinuity` Y |
| `Taxonomic_assessment` | HIGH_CONFIDENCE_TAXONOMIC_CHIMERA, TAXONOMIC_CHIMERA, TAXONOMIC_SHIFT_ONLY, NO_TAXONOMIC_SHIFT, or NA |
| `n_clipping_supported_transition_contigs` | Contigs with `Clipping_supported_transition` Y |
| `n_high_confidence_taxonomic_chimera_contigs` | Contigs with `High_confidence_taxonomic_chimera` Y |
| `n_taxonomic_shift_only_contigs` | Contigs with a transition but `Clipping_supported_transition` not Y |
| `n_taxonomic_shift_regions` | Transition regions |
| `n_clipping_supported_transition_regions` | Regions with `Clipping_supported_transition` Y |
| `n_high_confidence_taxonomic_chimera_regions` | Regions with `High_confidence_taxonomic_chimera` Y |
| `n_taxonomic_shift_only_regions` | Regions without validated clipping |
| `Clipping_supported_transition_contig_ids` | Contigs with `Clipping_supported_transition` Y |
| `High_confidence_taxonomic_chimera_contig_ids` | Contigs with `High_confidence_taxonomic_chimera` Y |
| `Taxonomic_shift_only_contig_ids` | Contigs with a transition only |
| `High_confidence_taxonomic_chimera` | Y when >= 1 contig has `High_confidence_taxonomic_chimera` Y |
