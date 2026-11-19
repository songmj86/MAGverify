# Project architecture

MAGverify separates the public command-line interface from the existing research-grade analyses.

## Repository layout

- `bin/MAGverify` runs the project from a checkout without installation.
- `magverify_cli.py` registers the Python entry point. `scripts/core_magverify_cli.py` implements input handling, stage selection, dependency checks, and result publishing.
- `modules/01_mapping.sh` through `modules/10_summary.sh` are the ten numbered analysis modules. Run them through `MAGverify MODULE` rather than editing the scripts.
- `scripts/` contains helpers named by their owning analysis (e.g. `circularity_*.py`, `chimera_*.py`, `summary_*.py`). Shared setup and utilities use `prep_*.py`, `core_*.py` or `common_*.py`. The redundant `tools/` wrapper scripts were removed after their calls were replaced by direct CLI/validator invocations.
- `modules/lib_config.sh`, `modules/lib_qc.sh`, and `modules/lib_scg.sh` are internal support code.
- `modules/defaults.sh` and `modules/stage_paths.sh` are internal defaults and computed output paths; they replace the earlier `config/project.sh` and ten `config/modules/*.sh` files. They are **not** user configuration files.
- `tests/` checks CLI behavior, input validation, summary consistency, and script compatibility.
- `docs/` explains installation, database preparation, accepted inputs, dependencies and tool versions. The README lists the ten analyses in execution order.

## Input and output contract

All inputs should be provided as command-line arguments. The same `-o` is reused across modules. Intermediate files stay in numbered analysis directories. Final public tables are placed (as links to their original results) in `OUTPUT/results/`; the summary produces `summary_mag.tsv` and `summary_contig.tsv`.

Output names and internal paths are managed by the CLI, not by a user-edited module configuration. The CLI may prepare the combined assembly and MAG-contig membership table as internal processing prerequisites. Chimera computes gene taxonomy internally; SCG annotation is built beforehand by `prep-scg`. Summary validates completed outputs without rerunning analyses.

## Workflow dependencies

Mapping provides shared read evidence. Clipping uses read mapping. False closure uses validated clipping, assembler circularity and SCGs. Chimera reuses mapping, clipping, gene taxonomy, and SCGs. False duplication uses mapping and clipping. Summary combines results from completed analyses only. See [Modules](MODULES.md) and the [README](../README.md).

## Testing scope

The automated test suite uses small fixtures and command stubs. Passing CI demonstrates packaging, syntax, command-line handling, and consistency checks, but **does not establish biological correctness, tool-version compatibility, or complete HiFi pipeline reproducibility**. Validate a representative real dataset and archive resolved environment versions before declaring a production release.

## Design patterns consulted

The layout adopts conventional separation of package code, user documentation, environmental dependencies, unit tests and continuous integration used by established bioinformatics repositories. Examples: [CheckM2](https://github.com/chklovski/CheckM2), [SemiBin](https://github.com/BigDataBiology/SemiBin), [geNomad](https://github.com/apcamargo/genomad), and [anvi'o](https://github.com/merenlab/anvio). These are organizational references, not software dependencies or imported code.

The preparation-only `modules/00_prep_chimera.sh` uses manually downloaded GTDB R232 inputs to build the MMseqs2 taxonomy database. It is invoked once by the `MAGverify prep-chimera` CLI command and is not one of the ten analysis modules.
