# Script index

The `MAGverify` command runs the shell modules in `modules/` (`00_prep_chimera.sh`, `00_prep_taxonomy.sh`, `00_prep_scg.sh`, `00_prep_circularity.sh`, `01_mapping.sh` to `10_summary.sh`, and the shared `lib_*.sh` helpers). The Python scripts below are internal helpers that those modules call; they are not separate user commands.

Each script name starts with the module that uses it. A script with several steps takes the step as its first argument (for example `python3 scripts/clipping_summarize.py sites ...`); run it without arguments to list the steps.

| Script | Purpose | Steps |
| --- | --- | --- |
| `chimera_call.py` | Taxonomic chimera calls: transition candidates, HiFi validation, MAG summary and clipping-based final calls. | `evaluate`, `validate`, `mag`, `reclassify` |
| `chimera_gene_taxonomy.py` | Gene taxonomy for the chimera module: gene selection, prefilter and MMseqs2 result tables. | `select`, `filter`, `mmseqs`, `merge` |
| `circularity_finalize.py` | Circularity final tables for all assembler groups. | `merge`, `refresh`, `add-junction`, `contig` |
| `circularity_junction.py` | Circularity END\|START junction evidence for one assembler group. | `summarize`, `integrate` |
| `circularity_reference.py` | Circularity reference: match MAG contigs to their assembler and build the rotated reference. | `match`, `rotate` |
| `clipping_summarize.py` | Clipping tables for all assemblers: site table, per-assembler density and contig/MAG summaries. | `sites`, `assembler`, `mag` |
| `clipping_validate.py` | Clipping per assembler: input checks, anvi'o candidate normalization and HiFi read validation. | `check`, `normalize`, `validate` |
| `core_magverify_cli.py` | Command-line entry point for MAGverify. | - |
| `core_progress_reader.py` | Stream read files to stdout and report progress with the remaining time on stderr. | - |
| `core_public_tables.py` | Write the tables for people (plain column names, Y/N/NA, no paths) from the internal tables in `.cache/`. | `migrate`, `publish` |
| `false_closure.py` | False closure: contig circularity, breakpoint detection and contig/MAG summaries. | `resolve`, `detect`, `mag`, `contig` |
| `false_duplication.py` | False duplication: candidate pairs from self-alignment, HiFi validation and clipping co-occurrence. | `detect`, `annotate` |
| `phantom_mag.py` | Create the MAG-associated region, MAG-level and contig-level phantom tables. | - |
| `prep_assembler_circularity.py` | Build the assembler-reported circularity tables from native assembler outputs. | - |
| `prep_build_mag_contig_map.py` | Build the MAG-to-contig map (MAG_ID, Contig) from MAG FASTAs and the all-assemblies FAI. | - |
| `prep_genes.py` | Shared gene predictions: protein-to-coordinate table and anvi'o external gene calls. | `link`, `anvio` |
| `prep_gtdbtk_summary_to_mag_taxonomy.py` | Convert GTDB-Tk summaries to the MAGverify MAG taxonomy TSV. | - |
| `repeats_tables.py` | Repeat tables: evidence summary, normalization, per-assembler columns and consolidated region/HSP tables. | `evidence`, `normalize`, `assembler`, `consolidate` |
| `summary_tables.py` | Integrated QC summary: build the MAG and contig tables, put the overall call first and validate them. | `build`, `arrange`, `validate` |
| `variants_tables.py` | Assembly variant tables: per-assembler classification, merge across assemblers and MAG summaries. | `classify`, `merge`, `mag` |
