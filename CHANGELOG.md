## Single-environment installation

- Use `conda env create -f environment.yml` followed by `conda activate MAGverify` once for all modules.
- Install anvi'o v9 from Bioconda as `anvio-minimal=9` and install the MAGverify CLI through the YAML `pip:` section.
- Remove the separate anvi'o environment file; no module command needs an anvi'o environment argument.
- Keep `-a` (assembler FASTA directory, formerly `--contigs-dir`), `--output-dir` and automatic `mapping_manifest.tsv` reuse as the public workflow.

# Changelog

## Unreleased

- `repeats`: INTRA BLAST compares each contig only with itself, and INTER BLAST compares each contig only with the other contigs of its MAG. Comparisons whose hits were discarded are no longer computed; results are unchanged. Longest contigs and largest MAGs start first, and BLAST runs as `-t` single-thread jobs.
- The MAG-contig map is rewritten only when its content changes, so completed module results are reused on reruns.
- MAG folder option renamed from `-m` to `-g` (also in `prep_gtdbtk_summary_to_mag_taxonomy.py`); `-m` is now the memory limit.
- New `MAGverify prep-scg -a DIR -o DIR -t N` builds the anvi'o SCG annotation (`00_prep/SCG/`) before the analysis modules. `false-closure` and `chimera` no longer run anvi'o for SCGs; `false-closure` stops if the table is missing and `chimera` continues without SCG annotation.
- "Database preparation" renamed to "Input preparation" (`docs/INPUT_PREPARATION.md`); README sections renumbered.
- Mapping outputs renamed from `LR_to_<assembler>.<class>.sorted.bam` to `<assembler>.<class>.bam` (flagstat: `<assembler>.flagstat.txt`). Cache fingerprints moved to `01_mapping/.cache/`. `map` renames existing files with the old names, so completed BAMs are still reused.
- `map` no longer builds secondary-only and supplementary-only BAMs, which no module used; they are removed from `mapping_manifest.tsv` and deleted on the next `map` run. Modules read the manifest by column name.
- `false-closure`: fixed the missing `classify_circularity` import. Ribosomal-protein family counts now come from `scripts/false_closure_scg_utils.py` (families whose HMM name starts with `Ribosom`, counted within Bacteria_71 or Archaea_76, whichever has more); the assembler circularity table is read with the same parser as `circularity`.
- `map`, `circularity` and `variants` process one assembler at a time with all threads; `map` builds the alignment-class BAMs of one assembler while the next is mapped. `clipping`, `phantom` (zero-depth scan) and `false-duplication` (minimap2 self-alignment) process assemblers in parallel and divide `-t` between them; `PARALLEL_JOBS=N` limits the number of parallel jobs.
- Memory limit option `-m GB` (formerly `--memory-gb`, `variants` only) for `map`, `circularity` and `variants`; in `map` and `circularity` half of the limit goes to samtools sort.
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
