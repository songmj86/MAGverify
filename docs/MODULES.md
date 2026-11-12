# Modules and dependencies

Before using `MAGverify chimera`, build the external GTDB R232 taxonomy database once with `MAGverify prep-chimera --output-dir /path/to/MAGverify_output --threads 100`. The CLI invokes `modules/00_prep_chimera.sh` internally. Preparation is distinct from the ten analysis modules. See [Database preparation](../PROTOCOL.md#database-preparation-once-before-analysis).

MAGverify modules are run individually. The preprocessing required by an analysis is performed internally, and all paths are supplied on the command line.

| Order | Command | Function | Main prerequisites | Result files |
| --- | --- | --- | --- | --- |
| 1 | `MAGverify map` | HiFi mapping per assembly | Original assemblies, reads | `mapping_manifest.tsv` |
| 2 | `MAGverify circularity` | Circular-junction evidence | MAGs, reads, original assemblies | `circularity_mag.tsv` |
| 3 | `MAGverify clipping` | Clipping discovery and read validation | Mapping, MAGs, anvi'o | `clipping_sites.tsv`, `clipping_mag.tsv` |
| 4 | `MAGverify false-closure` | False closure assessment and SCG preparation | Clipping, mapping, assembler circularity | `false_closure_contig.tsv`, `false_closure_mag.tsv` |
| 5 | `MAGverify phantom` | Phantom sequence assessment | Mapping, reads | `phantom_contig.tsv`, `phantom_mag.tsv` |
| 6 | `MAGverify repeats` | Repeat-region summary | Original assemblies, MAGs | `repeats_contig.tsv`, `repeats_mag.tsv` |
| 7 | `MAGverify chimera` | Gene taxonomy, SCG annotation and taxonomic chimera | Clipping, mapping, MMseqs2 GTDB DB | `gene_taxonomy.tsv`, `chimera_contig.tsv`, `chimera_mag.tsv` |
| 8 | `MAGverify false-duplication` | False-duplication validation | Mapping, clipping | `false_duplication_contig.tsv`, `false_duplication_mag.tsv` |
| 9 | `MAGverify variants` | Assembly allele support | Mapping, reads, anvi'o | `variants_contig.tsv`, `variants_mag.tsv` |
| 10 | `MAGverify summary` | Integrated QC and result validation | Completed module tables | `summary_mag.tsv`, `summary_contig.tsv` |

For executable examples with all options and `--threads 100`, see [PROTOCOL.md](../PROTOCOL.md). Use `MAGverify MODULE --help` for alternative inputs such as `--source-assembly`, `--mag-map`, `--taxonomy` and `--scg-table`.

Gene-level taxonomy is implemented within `07_chimera.sh`; `lib_scg.sh` and `modules/lib_config.sh` are internal helpers, not standalone user commands. `modules/lib_config.sh` retains common default paths and analysis parameters; no user edits are necessary.

The repeat and false-duplication modules remain separate since their biological and quality-control endpoints differ. Read mapping remains accessible as an independent module because multiple other analyses reuse the BAM evidence.
