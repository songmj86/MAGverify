# Modules and dependencies

Before `MAGverify chimera`, download the GTDB R232 protein archive and both taxonomy TSVs manually, then build a local MMseqs2 database with `MAGverify prep-chimera --gtdb-dir /path/to/GTDB_R232_files -o /path/to/MAGverify_db -t 100` (database: `/path/to/MAGverify_db/GTDB_R232`). The CLI invokes `modules/00_prep_chimera.sh` to build the database but never downloads files. Build the SCG annotation used by `false-closure` and `chimera` with `MAGverify prep-scg -a DIR -o OUTPUT -t 100`. Preparation is separate from the ten analysis modules. See [Input preparation](INPUT_PREPARATION.md).

MAGverify modules are run individually. The preprocessing required by an analysis is performed internally, and all paths are supplied on the command line.

| Order | Command | Function | Main prerequisites | Result files |
| --- | --- | --- | --- | --- |
| 1 | `MAGverify map` | HiFi mapping per assembly | Original assemblies, reads | `mapping_manifest.tsv` |
| 2 | `MAGverify circularity` | Circular-junction evidence | MAGs, original assemblies, reads (from mapping) | `circularity.MAG.tsv` |
| 3 | `MAGverify clipping` | Clipping discovery and read validation | Mapping, MAGs, anvi'o | `clipping.site.tsv`, `clipping.contig.tsv`, `clipping.MAG.tsv` |
| 4 | `MAGverify false-closure` | False closure assessment | Clipping, mapping, assembler circularity, SCG (`prep-scg`) | `false_closure.contig.tsv`, `false_closure.MAG.tsv` |
| 5 | `MAGverify phantom` | Phantom sequence assessment | Mapping (BAMs and reads) | `phantom.region.tsv`, `phantom.MAG.tsv` |
| 6 | `MAGverify repeats` | Repeat-region summary | Original assemblies, MAGs | `repeats.contig.tsv`, `repeats.MAG.tsv` |
| 7 | `MAGverify chimera` | Gene taxonomy and taxonomic chimera | Clipping, mapping, MMseqs2 GTDB DB, SCG (`prep-scg`) | `gene_taxonomy.gene.tsv`, `chimera.contig.tsv`, `chimera.MAG.tsv` |
| 8 | `MAGverify false-duplication` | False-duplication validation | Mapping, clipping | `false_duplication.contig.tsv`, `false_duplication.MAG.tsv` |
| 9 | `MAGverify variants` | Assembly allele support | Mapping (BAMs and reads), anvi'o | `variants.contig.tsv`, `variants.MAG.tsv` |
| 10 | `MAGverify summary` | Integrated QC and result validation | Completed module tables | `summary.MAG.tsv`, `summary.contig.tsv` |

For executable examples with `-t 100`, see the [README](../README.md). Use `MAGverify MODULE -h` for supported options such as `--mag-map`, `--taxonomy` and `--scg-table`.

Gene-level taxonomy is implemented within `07_chimera.sh`; `lib_scg.sh` is run by `MAGverify prep-scg`; `modules/lib_config.sh` is an internal helper, not a standalone user command. `modules/lib_config.sh` retains common default paths and analysis parameters; no user edits are necessary.

The repeat and false-duplication modules remain separate since their biological and quality-control endpoints differ. Read mapping remains accessible as an independent module because other analyses reuse the BAM evidence. Per-assembler BAMs and BAM indexes are stored in `01_mapping/`; `results/mapping_manifest.tsv` links to the manifest identifying the alignments.
