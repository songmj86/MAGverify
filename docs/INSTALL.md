# Installation

MAGverify uses **one Conda environment** for all ten analysis modules, including anvi'o. It is intended for Linux (and Linux under WSL where supported). The environment is named `MAGverify`. No separate anvi'o environment, activation script or `--anvio-env` argument is necessary.

## 1. Create and activate the environment

From the **MAGverify repository root**, run:

```bash
conda env create -f environment.yml
conda activate MAGverify
```

This one creation command installs the tools declared in `environment.yml` and the local MAGverify package itself. The file includes `anvio-minimal=9` from Bioconda, together with minimap2, samtools, bedtools, BLAST+, MMseqs2, Prodigal, Meryl and HMMER. Anvi'o's Python dependencies are resolved through its Conda recipe. `pip:` inside `environment.yml` installs the MAGverify CLI automatically; there is **no additional manual installation command**. You may use `mamba env create -f environment.yml` as an alternative to Conda.

**Important:** change to the repository root (`cd MAGverify`) before invoking `conda env create`, since the YAML includes a local editable package (`-e .`). Keep this repository directory after installation, since the editable installation uses its Python files and shell modules.

## 2. Verify installation

With the environment activated, run:

```bash
which MAGverify
which anvi-profile
which anvi-script-find-misassemblies
which anvi-gen-contigs-database
which anvi-run-hmms
which anvi-script-reformat-bam
which anvi-export-table
MAGverify doctor
MAGverify list
python -c 'import anvio; print(anvio.__version__)'
```

The listed programs should be found under the same Conda environment's `bin` directory. `MAGverify doctor` lists available executables, but does not replace a functional test. Confirm that any prerequisite programs not found on PATH are installed before processing data. The `anvi-self-test --suite mini --no-interactive` command can be used for a separate anvi'o functional test.

`MAGverify doctor` also prints `MAGverify source`, the repository directory the command runs from. The installation is editable, so it keeps using the directory where the environment was created. If you update or move the repository, or `MAGverify <module> --help` lacks an option shown in this documentation, reinstall from the current repository root:

```bash
conda activate MAGverify
cd /path/to/MAGverify
pip install -e . --no-deps
MAGverify doctor          # MAGverify source should now be this directory
```

`--no-deps` leaves the Conda-managed tools untouched. Reinstalling does not modify files in the previous directory, so jobs already running from it are not affected.

Start a new shell session using only:

```bash
conda activate MAGverify
```

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
  -m /data/MAGs \
  -o /data/MAGverify_output \
  -t 100
```

Continue the remaining modules in the order specified in the [README](../README.md). All commands accept `-a` for the directory of original assembler FASTA files and `-o` for the common output root. Neither `--anvio-env` nor an explicit `--mapping-manifest` is needed when using the same `-o` for every stage. The manifest is created by `MAGverify map` in `OUTPUT/01_mapping/mapping_manifest.tsv` and linked under `OUTPUT/results/mapping_manifest.tsv`.

## 4. Database preparation

The environment does not include the large GTDB datasets or GTDB-Tk. Follow [Database preparation](DATABASE_PREPARATION.md) to download fixed GTDB R232 files, build the MMseqs2 protein taxonomy database and prepare the MAG taxonomy table.

## 5. Validation and reproducibility

```bash
conda list --explicit > MAGverify_conda.lock.txt
MAGverify map --help
```

**Scope of verification:** repository syntax checks and CLI fixture tests verify the command interface. An actual Conda dependency solve, anvi'o functional run and complete HiFi analysis have not been performed in the packaging environment. The Conda specification is therefore an install recipe rather than a tested lockfile. Compare key QC outputs with the original results before scientific use. See [TOOL_VERSIONS.md](TOOL_VERSIONS.md).
