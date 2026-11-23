## Single-environment installation

- Renamed `Both_supported_circularity` to `High_confidence_circular` in `summary.MAG.tsv` (same criteria).

- Use `conda env create -f environment.yml` followed by `conda activate MAGverify` once for all modules.
- Install anvi'o v9 from Bioconda as `anvio-minimal=9` and install the MAGverify CLI through the YAML `pip:` section.
- Remove the separate anvi'o environment file; no module command needs an anvi'o environment argument.
- Keep `-a` (assembler FASTA directory, formerly `--contigs-dir`), `--output-dir` and automatic `mapping_manifest.tsv` reuse as the public workflow.

# Changelog

## Unreleased

- Clarified circularity evidence in `summary.MAG.tsv` using separate assembler-reported, read-supported, and both-supported flags. Replaced the ambiguous `Circular` summary alias; circularity algorithm unchanged.
- Simplified public summary flags in module order, renamed `Taxonomic_chimera` to `Chimera`, moved `n_flagged_checks` and `n_checks_not_run` to the end, and retained skipped-check names and complementary diagnostics in the detailed summary. The public MAG flagged-check count now covers only the seven displayed diagnostic checks, rather than the previously separate false-closure check.

- Public summary tables (`summary.MAG.tsv`, `summary.contig.tsv`, and `summary_detail.*.tsv`) list checks in end-to-end module order, with aggregate columns retained near the start/end. Detailed module columns follow the same order. Existing public summary tables are automatically republished from cached results when the previous column order is detected. No QC calls or values change.

<<<<<<< HEAD
- Documentation: `docs/` keeps `INPUT_PREPARATION.md` and `OUTPUTS.md` (result table columns). Installation is in the README; `INSTALL.md`, `MODULES.md`, `ARCHITECTURE.md`, `INPUT_FORMATS.md`, `TOOL_VERSIONS.md` and `SCRIPT_INDEX.md` were removed or merged.
=======
- Documentation examples use inline `/path/to/...` paths as in the README; the shell variable block (`In_Contig`, `In_Read`, `In_MAG`, `Out`) and the "Progress and resumption" paragraph are removed. Fixed the broken GTDB-Tk link (`ecogenomics/gtdtk` to `Ecogenomics/GTDBTk`).
- Documentation: `docs/` keeps `USAGE.md` (preparation and analysis commands) and `OUTPUTS.md` (result table columns). `INPUT_PREPARATION.md` and `MODULES.md` are merged into `USAGE.md`, which removes the duplicated preparation command block and output table and keeps the four preparation sections unnumbered (`#gtdb-r232-database`, `#single-copy-gene-annotation`, `#assembler-circularity`, `#mag-taxonomy`). Module anchors are unchanged. Installation is in the README; `INSTALL.md`, `ARCHITECTURE.md`, `INPUT_FORMATS.md`, `TOOL_VERSIONS.md` and `SCRIPT_INDEX.md` were removed or merged.
>>>>>>> 645efbd (Update MAGverify)
- `00_prep/all_assemblies.fasta` renamed `00_prep/all_assemblies.fa` (FASTA outputs use `.fa`; Pyrodigal outputs keep `.faa`/`.gff`). Existing files are renamed automatically.
- `-o` given as `OUTPUT/00_prep` or a module folder is read as `OUTPUT`. The old visible `00_prep/assembler_links/` is removed (the links are made in `00_prep/.cache/assembler_links/` when needed).
- `prep-taxonomy` writes `00_prep/gtdbtk_taxonomy.MAG.tsv` (earlier `taxonomy.MAG.tsv` and `gtdb_taxonomy.MAG.tsv` are renamed automatically). `-o` given as a file path (for example `.../x.tsv`) now stops with a message instead of creating a folder with that name.
- Every prep command has its module in `modules/`: `00_prep_chimera.sh`, `00_prep_taxonomy.sh`, `00_prep_scg.sh` (was `lib_scg.sh`), `00_prep_circularity.sh`.
- `chimera validate`: event windows are computed once and looked up through a per-contig interval index (same evidence, faster on contigs with many events); `chimera evaluate` reads `samtools bedcov` output line by line.
- `false-duplication detect`: `samtools view -F 260` skips secondary alignments, which were discarded afterwards anyway.
- Code comments and merge leftovers removed; scg.contig.tsv column `source` renamed `HMM_collection` (both are read).
- Readable result tables: each module folder now holds tables with plain column names (`Copy1_left_supported`, `HiFi_END_START_supported`, `Small_genome_proxy`, ...), `Y`/`N`/`NA` values and no file paths or repeated columns; the full tables the modules read are in `.cache/` (written by `scripts/core_public_tables.py`). Tables of earlier versions are moved to `.cache/` automatically and the result tables are written on the next run of any module.
- `summary.MAG.tsv` and `summary.contig.tsv` are at the top of the output folder (were in `10_summary/`): one `Y`/`N`/`NA` column per check, `n_flagged_checks`, counts and `Checks_not_run`; `MAG_integrated_QC_call`, `Circularity_eligibility`, `excessive_repeat_call` and `false_duplication_clipping_summary` are removed.
- `map`: `01_mapping/mapping.assembler.tsv` (read counts); the manifest with paths moved to `01_mapping/.cache/mapping_manifest.tsv`.
- `phantom.dataset.tsv` and `repeats.dataset.tsv` are kept in `.cache/` only.
- `00_prep/` has no subfolders: `genes.faa/gff/coordinates.tsv`, `scg.contig.tsv` (was `SCG/SCG_from_anvio.original_contig.tsv`), `taxonomy.MAG.tsv` (was `gtdb_taxonomy.MAG.tsv`), `taxonomy.gene.tsv` (was `gene_taxonomy/GTDB_R232/...all_contigs.tsv`); anvi'o files, gene-taxonomy work and assembler links are in `00_prep/.cache/`. Existing files are moved automatically.
- `map`: the `samtools flagstat` counts are columns of `mapping_manifest.tsv`; the `.flagstat.txt` files are in `01_mapping/.cache/`.
- `circularity`: new `circularity.contig.tsv` (every contig).
- `false-closure`: `false_closure.breakpoint.tsv` is folded into `false_closure.contig.tsv` (`Validated_clipping_breakpoints`, `Breakpoint_details`); the full table stays in `.cache/`.
- `phantom`: new `phantom.contig.tsv` (every contig).
- `false-duplication`: boundary values are Y/N/NA (were SUPPORTED/UNSUPPORTED/NOT_ANALYZABLE/NOT_EVALUABLE_CONTIG_END); columns `LR_COPY*_SUPPORTED`, `LR_ASSESSMENT`; duplicated `HIFI_PRIMARY_SUPPLEMENTARY_*` and pair-count columns removed from the pair and contig tables.
- Fixed: the merged `chimera_call.py`, `false_closure.py` and `false_duplication.py` had lost their `@dataclass` decorators (`ClipRecord() takes no arguments`).
- Python 3.10 or newer is required.
- `circularity` always runs minimap2 with `--secondary-seq`, as `map` does (`-ax map-hifi -p 1 --secondary-seq`); the version check and its "--secondary-seq is unavailable" note are removed.
- Messages: "eligible MAG" now reads "single-contig MAG"; "source assembly/assembler" wording removed from messages and documentation.
- Final tables in `results/` are named `MODULE.LEVEL.tsv` (`repeats.MAG.tsv`, `repeats.contig.tsv`, `clipping.site.tsv`, `phantom.region.tsv`, ...). MAG tables start with `MAG_ID`, `Assembler`; contig, site, region and gene tables start with `Contig`, `MAG_ID` (`NA` outside MAGs), `Assembler`. Contig and assembler ID columns are renamed to these names; module directories are unchanged. Links with the former names are removed.
- Help text: file types shown as `(*.fa)` and `(*.fastq.gz; *.fq.gz)`; `--assembler-circularity` described as `(prep-circularity output)`.
- Help output no longer lists `-h`; memory options show `N` (GB) as their value.
- New `MAGverify prep-circularity -c FILE ... -g DIR -o DIR` builds `00_prep/assembler_circularity.MAG.tsv` and `assembler_circularity.contig.tsv` (with assembler names) from hifiasm-meta GFA/FASTA, metaMDBG FASTA and metaFlye `assembly_info.txt`. `circularity` and `false-closure` use it automatically; `--assembler-circularity` still overrides it.
- `circularity` without `--assembler-circularity`: `Circular` now reports the read-supported call; `High_confidence_circular` and `Additional_circular_candidate` are `NA`. Previously every MAG was `N`.
- `00_prep/combined_assembly.fasta` renamed to `00_prep/all_assemblies.fasta` (the file concatenates the assembler FASTAs; it is not a co-assembly). Existing files are renamed automatically.
- `repeats` accepts `-m GB`: half goes to sorting, and parallel BLAST jobs are limited to about 2 GB each.
- `repeats`: INTRA BLAST compares each contig only with itself, and INTER BLAST compares each contig only with the other contigs of its MAG. Comparisons whose hits were discarded are no longer computed; results are unchanged. Longest contigs and largest MAGs start first, and BLAST runs as `-t` single-thread jobs.
- The MAG-contig map is rewritten only when its content changes, so completed module results are reused on reruns.
- MAG folder option renamed from `-m` to `-g` (also in `prep_gtdbtk_summary_to_mag_taxonomy.py`); `-m` is now the memory limit.
- New `MAGverify prep-scg -a DIR -o DIR -t N` builds the anvi'o SCG annotation (`00_prep/SCG/`) before the analysis modules. `false-closure` and `chimera` no longer run anvi'o for SCGs; `false-closure` stops if the table is missing and `chimera` continues without SCG annotation.
- "Database preparation" renamed to "Input preparation" (`docs/INPUT_PREPARATION.md`); README sections renumbered.
- Mapping outputs renamed from `LR_to_<assembler>.<class>.sorted.bam` to `<assembler>.<class>.bam` (flagstat: `<assembler>.flagstat.txt`). Cache fingerprints moved to `01_mapping/.cache/`. `map` renames existing files with the old names, so completed BAMs are still reused.
- `map` no longer builds secondary-only and supplementary-only BAMs, which no module used; they are removed from `mapping_manifest.tsv` and deleted on the next `map` run. Modules read the manifest by column name.
- `false-closure`: fixed the missing `classify_circularity` import. Ribosomal-protein family counts now come from `scripts/false_closure_scg_utils.py` (families whose HMM name starts with `Ribosom`, counted within Bacteria_71 or Archaea_76, whichever has more); the assembler circularity table is read with the same parser as `circularity`.
- `map`, `circularity` and `variants` process one assembler at a time with all threads; `map` builds the alignment-class BAMs of one assembler while the next is mapped. `clipping`, `phantom` (zero-depth scan) and `false-duplication` (minimap2 self-alignment) process assemblers in parallel and divide `-t` between them; `PARALLEL_JOBS=N` limits the number of parallel jobs.
- Memory limit option `-m N` (GB; formerly `--memory-gb`, `variants` only) for `map`, `circularity` and `variants`; in `map` and `circularity` half of the limit goes to samtools sort.
- Help text shortened: `-r` Folder of HiFi reads, `-a` Folder of assembled contigs, `-g` Folder of MAGs, `-t` Threads, `-o` Output folder. All previously accepted file extensions are still accepted.
- Removed explanatory comments from modules and scripts.
- `circularity`: `--assembler-circularity` is optional (assembler calls become `NA`; previously the final step failed without it). `--mag-map` and `--mapping-manifest` removed from `circularity`, where they had no effect.
- Help lists main inputs first (`-r`, `-a`, `-m`, then `-t`, `-o`) with accepted file extensions; help is `-h` only.
- `prep-chimera` always writes the database to `-o DIR` as `DIR/GTDB_R232` (`--database-prefix` removed). The GTDB-Tk converter takes `-o DIR` and writes `DIR/MAG_taxonomy.tsv`, so both chimera inputs share one directory.
- GTDB-Tk removed from `environment.yml`; MAGverify only converts existing GTDB-Tk summaries. The GTDB-Tk classification section was removed from `docs/INPUT_PREPARATION.md`.
- Short options only for common inputs: `-o` (output), `-t` (threads), `-r` (HiFi read directory), `-a`, `-m`. `--output-dir`, `--output-directory`, `--threads` and `--reads` were removed. `map` no longer lists `--mapping-manifest`.
- Mapping BAMs are written directly to `01_mapping/` (no `long_read/` subfolder).
- Assembler FASTA directory option renamed from `--contigs-dir` to `-a`.
- `--dry-run`, `--mag-fasta`, `-m` and `-i/--input` removed from the user interface; MAG input is `-m` only (formerly `--mag-dir`).
- `map` records the HiFi read file in `mapping_manifest.tsv` (`Reads`); circularity, phantom and variants reuse it when `--reads` is omitted.
- Clipping writes a contig-level table (`clipping_contig.tsv`).
- GTDB-Tk conversion writes only the MAGs in `-m` and ignores other GTDB-Tk genomes.
- Flattened repository structure: `modules/` for shell analysis stages and `scripts/` for Python helpers, without the old nested `magverify/` directory.
- Moved the GTDB R232 MMseqs2 preparation script to `modules/00_prep_chimera.sh`. Added `MAGverify prep-chimera --gtdb-dir DIR --output-dir DIR [--database-prefix PREFIX] --threads N` to build the taxonomy database from pre-downloaded GTDB files before chimera analysis. The ten analysis modules are unchanged.

- Public package layout with `MAGverify` command and a source-checkout launcher.
- Consolidated internal module path defaults; removed user-facing `config/modules/` files.
- Module-specific CLI, ordered README workflow, and automated syntax and unit tests.
- Integrated result consistency validation in the summary step.

The current codebase has not undergone complete end-to-end benchmarking with representative PacBio HiFi datasets.
