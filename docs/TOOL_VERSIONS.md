# Tool versions and reproducibility

## Versions recorded in the original study

The following versions are documented in the original activated-sludge workflow. They describe the upstream dataset generation and genome characterization; they are **not** all dependencies of MAGverify.

| Program or resource | Recorded version | Role |
| --- | --- | --- |
| SMRT Link | 25.3 | HiFi read processing |
| metaMDBG | 1.1 | Upstream long-read assembly |
| hifiasm-meta | 0.1 | Upstream long-read assembly |
| metaFlye | 2.9.5 | Upstream long-read assembly |
| HiFi-MAG-Pipeline | 2.0 | Upstream MAG recovery |
| dRep | 3.5.0 | Upstream dereplication |
| CheckM2 | 1.1.0 | Genome quality estimates |
| GTDB-Tk | 2.4.0 | Upstream genome taxonomy (GTDB R226) |

MAGverify's gene-level taxonomic-chimera workflow is written for a **GTDB R232 MMseqs2 taxonomy database**. GTDB R226 in the GTDB-Tk analysis and R232 in this gene taxonomy module are distinct resources; do not substitute one for the other without reviewing the method and annotations.

## MAGverify QC dependencies

The QC code invokes minimap2, samtools, BLAST+, bedtools, MMseqs2, Prodigal, Meryl and anvi'o. Their **exact original run versions are not recoverable from the repository or the recorded version list**, so `environment.yml` declares them without unverified version pins. `environment.yml` installs MAGverify and external dependencies in **one Conda environment** using Bioconda `anvio-minimal=9`, but v9 is not claimed to be the original analysis version. For reproducibility, compare anvi'o script formats with the original analysis environment before publication.

After installing or activating the environment originally used for analysis, save its resolved packages:

```bash
conda list --explicit > MAGverify_conda.lock.txt
conda list --export > MAGverify_conda_versions.txt
minimap2 --version
samtools --version | head -n 1
mmseqs version
blastn -version | head -n 1
meryl --version
python --version
```

With the `MAGverify` environment activated, run `anvi-self-test --suite mini --no-interactive` and `python -c 'import anvio; print(anvio.__version__)'` to check anvi'o. A platform-specific `conda list --explicit` file can recreate *that exact resolved environment* on a compatible platform. Do not call an example YAML a locked/reproduced environment until it has been solved and tested with actual data.
