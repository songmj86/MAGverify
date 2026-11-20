# Script index

The CLI command `MAGverify prep-chimera` invokes `modules/00_prep_chimera.sh` to build the R232 MMseqs2 protein taxonomy database from existing GTDB files, without downloading data, for `MAGverify chimera`; the ten analysis shell modules are `modules/01_mapping.sh` to `modules/10_summary.sh`.
Python helper filenames begin with the analysis module responsible for their execution. Internal Python helper scripts are called by the numbered stage shell scripts and are not separate user-facing CLI commands.
| Prefix | Purpose | Files |
| --- | --- | ---: |
| `core_` | User-facing MAGverify CLI implementation | 1 |
| `prep_` | MAG-to-contig mapping and optional GTDB-Tk summary conversion | 2 |
| `circularity_` | Single-contig circularity and HiFi junction validation | 6 |
| `clipping_` | Clipped-read assessment and validation | 6 |
| `false_closure_` | False closure and contig connections | 5 |
| `phantom_` | Phantom sequence evidence | 1 |
| `repeats_` | Repeat regions and evidence summaries | 4 |
| `chimera_` | Gene taxonomy and taxonomic chimera assessment | 9 |
| `false_duplication_` | False duplication and clipping integration | 2 |
| `variants_` | Assembly allele support | 3 |
| `summary_` | Integration and QC validation | 2 |
| `common_` | Shared table helper used in multiple modules | 1 |

## Consolidated helpers

- `chimera_select_sequences.py` merges former FASTA and cached-gene selection scripts. Internal operations are `fasta` and `genes`.
- `summary_finalize.py` merges compact complete-genome QC table generation with summary validation. Internal operations are `complete-genome` and `validate`.

These are internal Python script interfaces. The ten `MAGverify MODULE` commands, output file names, result directories, and analysis criteria are unchanged. Other scripts remain separate when different input/output contracts or intermediate files need independent recomputation.

## Script inventory

### core

- `core_magverify_cli.py`

### prep

- `prep_assembler_circularity.py`
- `prep_build_mag_contig_map.py`
- `prep_genes_to_anvio_gene_calls.py`
- `prep_gtdbtk_summary_to_mag_taxonomy.py`

### circularity

- `circularity_build_rotated_reference.py`
- `circularity_integrate_circularity_results.py`
- `circularity_merge_assembler_groups.py`
- `circularity_refresh_circularity_compat.py`
- `circularity_match_mags_to_assemblers.py`
- `circularity_summarize_circularity_junctions.py`

### clipping

- `clipping_check_clipping_inputs.py`
- `clipping_classify_clipping_regions.py`
- `clipping_normalize_anvio_clipping.py`
- `clipping_summarize_clipping_by_mag.py`
- `clipping_summarize_clipping_regions.py`
- `clipping_validate_clipping_spanning.py`

### false_closure

- `false_closure_contig_connection_utils.py`
- `false_closure_detect_false_closure.py`
- `false_closure_resolve_false_closure_contig_circularity.py`
- `false_closure_summarize_false_closure.py`
- `false_closure_summarize_false_closure_by_contig.py`

### phantom

- `phantom_summarize_phantom_by_mag.py`

### repeats

- `repeats_consolidate_outputs.py`
- `repeats_normalize_repeat_tables.py`
- `repeats_summarize_repeat_by_assembler.py`
- `repeats_summarize_repeat_evidence.py`

### chimera

- `chimera_evaluate_taxonomic_chimera_whole_contig_GTDB.py`
- `chimera_filter_gene_predictions_by_min_contig_gene_count.py`
- `chimera_link_prodigal_accessions.py`
- `chimera_merge_gene_taxonomy_tables.py`
- `chimera_merge_mmseqs_taxonomy_with_coordinates.py`
- `chimera_reclassify_taxonomic_chimera_by_clipping.py`
- `chimera_select_sequences.py`
- `chimera_summarize_taxonomic_chimera_by_mag_GTDB_whole.py`
- `chimera_validate_taxonomic_chimera_hifi.py`

### false_duplication

- `false_duplication_annotate_false_duplication_clipping.py`
- `false_duplication_detect_false_duplication.py`

### variants

- `variants_merge_assembly_variant_assembler_results.py`
- `variants_summarize_assembly_variants.py`
- `variants_summarize_assembly_variants_by_mag.py`

### summary

- `summary_build_integrated_qc_summary.py`
- `summary_finalize.py`

### common

- `common_export_validation_view.py`

