## Single-environment installation

- Use `conda env create -f environment.yml` followed by `conda activate MAGverify` once for all modules.
- Install anvi'o v9 from Bioconda as `anvio-minimal=9` and install the MAGverify CLI through the YAML `pip:` section.
- Remove the separate anvi'o environment file; no module command needs an anvi'o environment argument.
- Keep `--contigs-dir`, `--output-dir` and automatic `mapping_manifest.tsv` reuse as the public workflow.

# Changelog

## Unreleased

- Flattened repository structure: `modules/` for shell analysis stages and `scripts/` for Python helpers, without the old nested `magverify/` directory.
- Moved the GTDB R232 MMseqs2 preparation script to `modules/00_prep_chimera.sh`. Added `MAGverify prep-chimera --gtdb-dir DIR --output-dir DIR [--database-prefix PREFIX] --threads N` to build the taxonomy database from pre-downloaded GTDB files before chimera analysis. The ten analysis modules are unchanged.

- Public package layout with `MAGverify` command and a source-checkout launcher.
- Consolidated internal module path defaults; removed user-facing `config/modules/` files.
- Module-specific CLI, ordered README workflow, and automated syntax and unit tests.
- Integrated result consistency validation in the summary step.

The current codebase has not undergone complete end-to-end benchmarking with representative PacBio HiFi datasets.
