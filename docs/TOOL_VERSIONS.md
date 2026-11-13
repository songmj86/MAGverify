# Tool versions and reproducibility

## Versions recorded in the original study

The following versions are documented in the original activated-sludge workflow. They describe dataset generation and genome characterization; they are **not** all dependencies of MAGverify.

| Program or resource | Recorded version | Role |
| --- | --- | --- |
| SMRT Link | 25.3 | HiFi read processing |
| metaMDBG | 1.1 | Long-read assembly |
| hifiasm-meta | 0.1 | Long-read assembly |
| metaFlye | 2.9.5 | Long-read assembly |
| HiFi-MAG-Pipeline | 2.0 | MAG recovery |
| dRep | 3.5.0 | Dereplication |
| CheckM2 | 1.1.0 | Genome quality estimates |
| GTDB-Tk | 2.4.0 | Genome taxonomy (GTDB R226) |

MAGverify's gene-level taxonomic-chimera workflow is written for a **GTDB R232 MMseqs2 taxonomy database**. GTDB R226 in the GTDB-Tk analysis and R232 in this gene taxonomy module are distinct resources; do not substitute one for the other without reviewing the method and annotations.

