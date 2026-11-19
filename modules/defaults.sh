#!/usr/bin/env bash
# Runtime defaults. The CLI supplies input paths; no user configuration is necessary.
CONFIG_DIR="${CONFIG_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
# ASSEMBLY_FA is the combined assembly FASTA with unique contig identifiers.
# READS is the original PacBio HiFi FASTQ or FASTQ.GZ (required by 05_phantom).
# MAG_DIR contains one .fa/.fasta/.fna file per MAG; filename = MAG ID.
# ASSEMBLER_ASSEMBLY_DIR contains the original co-assemblies listed in ASSEMBLER_FASTAS.
# Use full paths. Keep OUT as the common parent of individual module output folders.
ASSEMBLY_FA="${ASSEMBLY_FA:-}"
READS="${READS:-}"
MAG_DIR="${MAG_DIR:-}"
ASSEMBLER_ASSEMBLY_DIR="${ASSEMBLER_ASSEMBLY_DIR:-}"
ASSEMBLER_FASTAS="${ASSEMBLER_FASTAS:-hifiasm_meta.fa metaflye.fa metaMDBG.fa}"
MAG_CONTIG_MAP="${MAG_CONTIG_MAP:-}"
MMSEQS_TAX_DB="${MMSEQS_TAX_DB:-}"
TAXONOMY_MODE="${TAXONOMY_MODE:-contig}"
CIRCULARITY_TSV="${CIRCULARITY_TSV:-}"
MAG_GTDB_TAXONOMY_TSV="${MAG_GTDB_TAXONOMY_TSV:-}"
GTDB_GENE_TAX_TSV="${GTDB_GENE_TAX_TSV:-}"
OUT="${OUT:-$CONFIG_DIR/results}"

# Leave these empty to use commands from PATH or the active Conda environment.
MINIMAP2_ENV="${MINIMAP2_ENV:-}"
MERYL_ENV="${MERYL_ENV:-}"
ANVIO_ENV_PREFIX="${ANVIO_ENV_PREFIX:-${CONDA_PREFIX:-}}"
# The CLI uses the active Conda environment for anvi'o; no separate CLI option is needed.
