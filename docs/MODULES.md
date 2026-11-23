# Modules and dependencies

Before `MAGverify chimera`, download the GTDB R232 protein archive and both taxonomy TSVs manually, then build a local MMseqs2 database with `MAGverify prep-chimera --gtdb-dir /path/to/GTDB_R232_files -o OUTPUT -t 100` (database: `OUTPUT/00_prep/GTDB_R232`, read by `chimera` automatically). Build the MAG taxonomy table with `MAGverify prep-taxonomy -d DIR -g DIR -o OUTPUT` (`OUTPUT/00_prep/gtdbtk_taxonomy.MAG.tsv`). The CLI invokes `modules/00_prep_chimera.sh` to build the database but never downloads files. Build the SCG annotation used by `false-closure` and `chimera` with `MAGverify prep-scg -a DIR -o OUTPUT -t 100`. Preparation is separate from the ten analysis modules. See [Input preparation](INPUT_PREPARATION.md).

MAGverify modules are run individually. The preprocessing required by an analysis is performed internally, and all paths are supplied on the command line.

| Order | Command | Function | Main prerequisites | Result files |
| --- | --- | --- | --- | --- |
| 1 | `MAGverify map` | HiFi mapping per assembly | Original assemblies, reads | `mapping.assembler.tsv` |
| 2 | `MAGverify circularity` | Circular-junction evidence | MAGs, original assemblies, reads (from mapping) | `circularity.MAG.tsv`, `circularity.contig.tsv` |
| 3 | `MAGverify clipping` | Clipping discovery and read validation | Mapping, MAGs, anvi'o | `clipping.site.tsv`, `clipping.contig.tsv`, `clipping.MAG.tsv`, `clipping.assembler.tsv` |
| 4 | `MAGverify false-closure` | False closure assessment | Clipping, mapping, assembler circularity, SCG (`prep-scg`) | `false_closure.MAG.tsv`, `false_closure.contig.tsv` |
| 5 | `MAGverify phantom` | Phantom sequence assessment | Mapping (BAMs and reads) | `phantom.region.tsv`, `phantom.contig.tsv`, `phantom.MAG.tsv` |
| 6 | `MAGverify repeats` | Repeat-region summary | Original assemblies, MAGs | `repeats.contig.tsv`, `repeats.MAG.tsv`, `repeats.region.tsv` (INTRA/INTER/TOTAL), `repeats.HSP.tsv` (INTRA/INTER) |
| 7 | `MAGverify chimera` | Gene taxonomy and taxonomic chimera | Clipping, mapping, MMseqs2 GTDB DB, SCG (`prep-scg`) | `chimera.region.tsv`, `chimera.contig.tsv`, `chimera.MAG.tsv` |
| 8 | `MAGverify false-duplication` | False-duplication validation | Mapping, clipping | `false_duplication.pair.tsv`, `false_duplication.contig.tsv`, `false_duplication.MAG.tsv` |
| 9 | `MAGverify variants` | Assembly allele support | Mapping (BAMs and reads), anvi'o | `variants.site.tsv`, `variants.contig.tsv`, `variants.MAG.tsv` |
| 10 | `MAGverify summary` | Integrated QC and result validation | Completed module tables | `summary.MAG.tsv`, `summary.contig.tsv` (top of the output folder) |

For executable examples with `-t 100`, see the [README](../README.md). Use `MAGverify MODULE -h` for supported options such as `--mag-map`, `--taxonomy` and `--scg-table`.

Gene-level taxonomy is implemented within `07_chimera.sh`; `00_prep_scg.sh` is run by `MAGverify prep-scg` (and `00_prep_chimera.sh`, `00_prep_taxonomy.sh`, `00_prep_circularity.sh` by the other prep commands); `lib_genes.sh` predicts genes once with Pyrodigal (`-p meta`) into `00_prep/genes.*`, shared by prep-scg, chimera and variants; `modules/lib_config.sh` is an internal helper, not a standalone user command. `modules/lib_config.sh` retains common default paths and analysis parameters; no user edits are necessary.

The repeat and false-duplication modules remain separate since their biological and quality-control endpoints differ. Read mapping remains accessible as an independent module because other analyses reuse the BAM evidence. Per-assembler BAMs and BAM indexes are stored in `01_mapping/`; `01_mapping/.cache/mapping_manifest.tsv` identifies the alignments.

## Clipping outputs

`03_clipping/` holds four tables; the full internal tables and everything else are in `03_clipping/.cache/` (per-assembler anvi'o output and intermediate tables, kept so that reruns reuse them).

| File | One row per | Content |
|---|---|---|
| `clipping.site.tsv` | anvi'o clipping candidate | All candidates of all assemblers and both regions, with `Assembler`, `MAG_ID` (`NA` for unbinned contigs), `region_type`, read validation and calls. Filter it for subsets: `Stringent_candidate == Y`, `Region_type == INTERNAL/TERMINAL`, `Clipping == Y`. |
| `clipping.contig.tsv` | contig | Counts and calls per contig, with `MAG_ID` |
| `clipping.MAG.tsv` | MAG | Counts and calls per MAG; `Clipping == Y` gives the MAGs with clipping |
| `clipping.assembler.tsv` | assembler (+ `ALL`) | INTERNAL and TERMINAL site counts, densities per analyzed Mbp and terminal/internal density ratios |

| Column | Definition |
|---|---|
| (reference) | HiFi reads are mapped to all contigs of the assembler that produced the evaluated contig (all mapped reads, no MAPQ filter) |
| Default_candidate | Position reported by `anvi-script-find-misassemblies` with its default clipping-ratio threshold > 0.8; `--min-dist-to-end 0`; no minimum coverage |
| Region_type | TERMINAL for the first/last 100 bp of a contig (`TERMINAL_CLIPPING_WINDOW_BP`); INTERNAL otherwise |
| Reference_boundary | Y only at the first/last base of a contig; kept for terminal statistics, but two-sided spanning is impossible, so `Clipping` cannot be Y |
| Exact_position_validation | PASS: >= 2 distinct clipped HiFi reads at the exact position and no spanning read; INSUFFICIENT_CLIPPED_QNAMES: only the >= 1-read audit criterion is met; BOUNDARY_UNRESOLVED at contig edges; FAIL otherwise |
| Clipping_ge1 | Permissive audit call: >= 1 distinct clipped HiFi read and no spanning read |
| Exact_position_spanning | Counter-evidence: one alignment has aligned reference bases (M/=/X) on both sides of the position; I/D/N patterns do not change this call |
| Exact_spanning_min_flank_bp_min/median/max | Shorter aligned flank among spanning reads (largest per read); audit only |
| n_unique_spanning_with_I/D/N_qnames, Spanning_I/D/N_length_bp_min/median/max | Spanning reads with an insertion next to, or a deletion/skip over, the position, and their lengths (largest per read); audit only |
| n_unique_spanning_with_any_IDN_qnames | Spanning reads with any of the I/D/N events above |
| Stringent_candidate | anvi'o coverage >= 10 and clipping ratio = 1.0 |
| Clipping_confidence | STRINGENT: validated and stringent; DEFAULT_ONLY: validated, below the stringent threshold; NOT_VALIDATED otherwise |
| Clipping | Final call: Y for non-boundary positions with >= 2 distinct clipped HiFi reads and no spanning read |
| Stringent_clipping | Y when `Clipping` is Y and the position is a stringent candidate |
| Internal_clipping / Terminal_clipping (contig, MAG) | Y when >= 1 `Clipping == Y` site lies in an INTERNAL / TERMINAL region; `_ge1` columns use `Clipping_ge1` |
| *_per_Mbp, *_density_ratio (assembler) | Site counts divided by the analyzed bp of each region, per Mbp; terminal/internal ratios are enrichment ratios, not error-rate ratios |

## Other contig-level outputs

| File | Content |
|---|---|
| `01_mapping/mapping.assembler.tsv` | One row per assembler: the `samtools flagstat` counts of the mapped-reads BAM (`Alignments_total`, `Primary`, `Secondary`, `Supplementary`, `Duplicates`, `Mapped`, `Primary_mapped`, `Mapped_percent`, `Primary_mapped_percent`). The paths of the assemblies, BAMs and reads are in `01_mapping/.cache/mapping_manifest.tsv` |
| `02_circularity/circularity.contig.tsv` | Every contig: `MAG_ID` (`NA` if unbinned), assembler-reported circularity, and for the contig of a single-contig MAG the same calls and END\|START read evidence as `circularity.MAG.tsv`. `Circularity_eligibility` is `N_not_single_contig` or `N_unbinned` for contigs that are not tested |
| `04_false_closure/false_closure.contig.tsv` | Every assembler-reported circular contig. `Breakpoint_details` lists each validated clipping breakpoint as `REGION@position(distance to end)\|L=connected/clipped reads[target contig:reads]\|R=...\|cov=,clip=,frac=\|spanning=\|evidence=Y/N\|reason`, separated by `;` |
| `05_phantom/phantom.contig.tsv` | Every contig: number, total bp, fraction and 0-based coordinates of phantom regions; `LR_PHANTOM` Y/N |
| `08_false_duplication/false_duplication.pair.tsv` | `Copy1/2_left/right_supported`: Y = boundary spanned by >= 2 HiFi reads; N = reads reach both sides but none spans; NA = not analyzable (too few reads or a contig end). `Assessment` starts with the pair call (Y/N/NA) followed by the reason |

## Progress

Each command prints its total run time (`Elapsed`). Read mapping in `map` and `circularity` and the BLAST steps in `repeats` also print progress with the remaining time every 5 minutes, for example `[06_repeats][INTRA BLAST chunks] 812/1930 (42%) | elapsed 20m | remaining 27m`. Set `PROGRESS_INTERVAL=SECONDS` to change the interval.
