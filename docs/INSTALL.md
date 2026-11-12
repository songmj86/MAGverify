# Installation

MAGverify uses **one Conda environment** for all ten analysis modules, including anvi'o. It is intended for Linux (and Linux under WSL where supported). The environment is named `MAGverify`. No separate anvi'o environment, activation script or `--anvio-env` argument is necessary.

## 1. Create and activate the environment

From the **MAGverify repository root**, run:

```bash
conda env create -f environment.yml
conda activate MAGverify
```

This one creation command installs the tools declared in `environment.yml` and the local MAGverify package itself. The file includes `anvio-minimal=9` and `gtdbtk=2.8.0` from Bioconda, together with minimap2, samtools, bedtools, BLAST+, MMseqs2, Prodigal, Meryl and HMMER. GTDB-Tk is used only when MAG-level taxonomy needs to be generated. Anvi'o's Python dependencies are resolved through its Conda recipe. `pip:` inside `environment.yml` installs the MAGverify CLI automatically; there is **no additional manual installation command**. You may use `mamba env create -f environment.yml` as an alternative to Conda.

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
gtdbtk --version
python -c 'import anvio; print(anvio.__version__)'
```

The listed programs should be found under the same Conda environment's `bin` directory. `MAGverify doctor` lists available executables, but does not replace a functional test. Confirm that any prerequisite programs not found on PATH are installed before processing data. The `anvi-self-test --suite mini --no-interactive` command can be used for a separate anvi'o functional test.

Start a new shell session using only:

```bash
conda activate MAGverify
```

## 3. Run each module

A single activation is sufficient for all modules:

```bash
MAGverify map \
  --assemblies-dir /data/assemblies \
  --reads /data/hifi.fastq.gz \
  --output-dir /data/MAGverify_output \
  --threads 100

MAGverify circularity \
  --assemblies-dir /data/assemblies \
  --mag-dir /data/MAGs \
  --reads /data/hifi.fastq.gz \
  --output-dir /data/MAGverify_output \
  --threads 100
```

Continue the remaining modules in the order specified in [PROTOCOL.md](../PROTOCOL.md). All commands accept `--assemblies-dir` for the directory of original assembler FASTA files and `--output-dir` for the common output root. Neither `--anvio-env` nor an explicit `--mapping-manifest` is needed when using the same `--output-dir` for every stage. The manifest is created by `MAGverify map` in `OUTPUT/01_mapping/assembler_mapping_outputs.tsv` and linked under `OUTPUT/results/mapping_manifest.tsv`.

## 4. Database preparation (not installed by Conda)

The single Conda environment installs the executables, not the large external GTDB data. The **Database preparation** section of [PROTOCOL.md](../PROTOCOL.md#database-preparation-once-before-analysis) gives the exact pinned R232 source URLs, the three files required for MMseqs2 (`gtdb_proteins_aa_reps_r232.tar.gz`, `bac120_taxonomy_r232.tsv`, `ar53_taxonomy_r232.tsv`), download/build code and verification steps.

```bash
PREP=/path/to/MAGverify_output/0_prep
MAGverify prep-chimera \
  --output-dir /path/to/MAGverify_output \
  --threads 100
```

The command creates the GTDB R232 MMseqs2 taxonomy database prefix at `$PREP/GTDB_R232`. To supply MAG taxonomy, use an existing `MAG_taxonomy.tsv` or run the bundled GTDB-Tk-summary converter on prior `classify_wf` results. Only if those results are absent, prepare the separate R232 GTDB-Tk reference package, run `gtdbtk classify_wf`, and convert its `summary.tsv` output as described in PROTOCOL.md. You do not need a second Conda environment and the main MAGverify code is unchanged.

## 5. Validation and reproducibility

```bash
conda list --explicit > MAGverify_conda.lock.txt
MAGverify map --help
```

**Scope of verification:** repository syntax checks and CLI fixture tests verify the command interface. An actual Conda dependency solve, anvi'o functional run and complete HiFi analysis have not been performed in the packaging environment. The Conda specification is therefore an install recipe rather than a tested lockfile. Compare key QC outputs with the original results before scientific use. See [TOOL_VERSIONS.md](TOOL_VERSIONS.md).
