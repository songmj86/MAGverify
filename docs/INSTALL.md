# Installation

MAGverify uses **one Conda environment** for all ten analysis modules, including anvi'o. It is intended for Linux (and Linux under WSL where supported). The environment is named `MAGverify`. No separate anvi'o environment, activation script or `--anvio-env` argument is necessary.

## 1. Create and activate the environment

From the **MAGverify repository root**, run:

```bash
cd /path/to/MAGverify
mamba env create -f environment.yml      # or: conda env create -f environment.yml
conda activate MAGverify
MAGverify doctor
```

`environment.yml` installs anvi'o 9 (`anvio-minimal`), minimap2, samtools, bedtools, BLAST+, MMseqs2, Pyrodigal, Meryl and HMMER, and registers the `MAGverify` command from this directory (`pip: -e .`). No other installation command is needed.

`MAGverify doctor` shows the code directory in use, every required program with its version, and ends with `All 22 required programs found.` It exits with an error and lists what is missing otherwise. For a separate functional test of anvi'o, run `anvi-self-test --suite mini --no-interactive`.

## 2. Keep the code directory in place

The `MAGverify` command runs the code in the directory the environment was created from (`MAGverify path` in `MAGverify doctor`).

- **Update the code in this directory** rather than in a copy elsewhere, for example `git pull` or `rsync -rc --exclude .git NEW_CODE/ /path/to/MAGverify/`. Avoid deleting and recreating the directory: open shells then point to the deleted directory (`pip` reports "The folder you are executing pip from can no longer be found"; run `cd` into the directory again).
- **If the directory was moved**, or `MAGverify` reports `ModuleNotFoundError: No module named 'magverify_cli'`, register the current directory once:

```bash
conda activate MAGverify
cd /path/to/MAGverify
pip install -e . --no-deps --no-build-isolation
MAGverify doctor          # MAGverify path should now be this directory
```

`--no-deps --no-build-isolation` leaves the Conda-managed tools untouched and needs no internet access. Jobs already running are not affected.

- **To update an existing environment** after `environment.yml` changes, run `mamba env update -f environment.yml --prune` from the repository root instead of creating it again.

## 3. Run each module

A single activation is sufficient for all modules:

```bash
MAGverify map \
  -a /data/assemblies \
  -r /data/hifi_reads \
  -o /data/MAGverify_output \
  -t 100

MAGverify circularity \
  -a /data/assemblies \
  -g /data/MAGs \
  -o /data/MAGverify_output \
  -t 100
```

Continue the remaining modules in the order specified in the [README](../README.md). All commands accept `-a` for the directory of original assembler FASTA files and `-o` for the common output root. Neither `--anvio-env` nor an explicit `--mapping-manifest` is needed when using the same `-o` for every stage. The manifest is created by `MAGverify map` in `OUTPUT/01_mapping/.cache/mapping_manifest.tsv`.

## 4. Input preparation

The environment does not include the large GTDB datasets or GTDB-Tk. Follow [Input preparation](INPUT_PREPARATION.md) to download fixed GTDB R232 files, build the MMseqs2 protein taxonomy database, prepare the MAG taxonomy table and build the anvi'o SCG annotation.

## 5. Validation and reproducibility

```bash
conda list --explicit > MAGverify_conda.lock.txt
MAGverify map -h
```

**Scope of verification:** repository syntax checks and CLI fixture tests verify the command interface. An actual Conda dependency solve, anvi'o functional run and complete HiFi analysis have not been performed in the packaging environment. The Conda specification is therefore an install recipe rather than a tested lockfile. Compare key QC outputs with the original results before scientific use. See [TOOL_VERSIONS.md](TOOL_VERSIONS.md).
