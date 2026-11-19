# Input files

## Original assembler FASTAs (`-a`)

Supply an existing directory with the original, unbinned assembler FASTA files in `.fa`, `.fasta` or `.fna` format. All FASTA files directly inside that directory are used. They should be original assemblies, **not** individual MAG FASTAs. Place the original assembler FASTAs together in one directory.

Clipping maps the original PacBio HiFi reads separately against these assemblies. This preserves the alignments to unbinned contigs, which would be lost if the reads were mapped only against the recovered MAGs.

Each FASTA record must use a distinct contig ID (the first non-whitespace string after `>`), including across separate assembler FASTAs. The tool rejects duplicate IDs rather than modifying sequence headers.

## MAG FASTAs (`--mag-dir`)

Supply an existing directory with one FASTA per MAG:

```text
/data/MAGs/Bin13L.fa
/data/MAGs/Bin23L.fa
/data/MAGs/Bin163L.fa
```

The filename without `.fa`, `.fna` or `.fasta` becomes the MAG ID. The MAG FASTA contig IDs and sequences must match the original assembler FASTAs. Circularity testing uses the appropriate single-contig MAGs.

## PacBio HiFi reads (`--reads`)

Supply PacBio HiFi reads in FASTQ or FASTQ.GZ format to `MAGverify map`. The mapping manifest records the read file, and circularity, phantom and variants reuse it when `--reads` is omitted. Modules assessing read support cannot infer it from assembly FASTAs alone.

## Automatically generated combined assembly

**Do not create or provide `combined.fasta` by hand.** When a module needs the union of source assemblies, the command-line interface builds `00_prep/combined_assembly.fasta` inside the output directory. It preserves the original contig identifiers and rejects duplicates. Source FASTA paths, sizes and modification times are recorded to determine whether this generated file can be reused.

## Mapping manifest (`--mapping-manifest`)

`MAGverify map` produces `01_mapping/mapping_manifest.tsv`, listing the assembler names, assembler FASTAs, HiFi mapping BAMs and the HiFi read file (`Reads`). Downstream modules can reuse this file without remapping. The manifest also identifies the source FASTAs from which the combined FASTA can be generated.

The tab-delimited manifest includes these columns:

```text
Assembler    Assembly_FASTA    HiFi_all_BAM    HiFi_primary_BAM    HiFi_secondary_BAM    HiFi_supplementary_BAM    HiFi_primary_supplementary_BAM    Output_directory    Reads
```

The manifest paths, associated BAMs and `.inputs.sha256` provenance files must remain accessible.

## Sequence identity and consistency

The combined FASTA represents the union of the original source assemblies. A MAG FASTA must match its source contig IDs and sequence lengths. Do not rename contigs between assembly, binning and QC. Do not substitute MAG FASTAs for the original assembler FASTAs, as this changes the mapping reference and the interpretation of read support.

## Assembler circularity (`--assembler-circularity`)

The false-closure module expects a tab-delimited table with a MAG identifier and assembler-reported circularity call. A minimal example is:

```text
MAG_ID\tcircularity
Bin13L\tY
Bin23L\tN
```

The accepted MAG identifier columns include `MAG_ID` and `ID`. The accepted circularity columns include `circularity`, `circular` and `assembler_circularity`. `Y` indicates an assembler-reported circular MAG and `N` indicates no such call. This input represents assembler metadata and is different from the HiFi junction-support result produced by `MAGverify circularity`.
