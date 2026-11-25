# MAGverify

MAGverify assesses the integrity of metagenome assemblies and metagenome-assembled genomes (MAGs) using PacBio HiFi reads. Its modules cover read mapping, circularity, clipping, false closure, phantom sequences, repeats, taxonomic chimeras, false duplication, assembly variants and an integrated summary.

## 1. Installation

Install MAGverify on Linux using Conda or Mamba from the repository root:

```bash
git clone https://github.com/songmj86/MAGverify.git
cd MAGverify
mamba env create -f environment.yml
conda activate MAGverify
MAGverify doctor
MAGverify list
```

MAGverify doctor checks the installation path and availability of required programs.
MAGverify list displays the available commands and analysis modules with brief descriptions.

## 2. Test dataset

Use the test dataset to check that MAGverify runs correctly before applying it to your own data. [input_dataset.tar.gz](https://huggingface.co/datasets/ATCC19718/MAGverify-resources/resolve/main/input_dataset.tar.gz) (1.1 GB) is available on [Hugging Face](https://huggingface.co/datasets/ATCC19718/MAGverify-resources). Download and extract it:

```bash
Data=/path/to/input_dataset
mkdir -p "$Data"

aria2c -c -x 8 -s 8 -k 8M \
  -d "$Data" -o input_dataset.tar.gz \
  "https://huggingface.co/datasets/ATCC19718/MAGverify-resources/resolve/main/input_dataset.tar.gz"

tar -xzf "$Data/input_dataset.tar.gz" -C "$Data" --strip-components=1
```

The dataset provides a set of inputs for MAGverify including 3 assemblies, 30 MAGs, 1 HiFi read file, 3 assembler circularity files and 2 GTDB-Tk results.

```
input_dataset
├── contigs
│   ├── metaMDBG.fa
│   ├── metaFlye.fa
│   ╰── hifiasm-meta.fa
├── MAGs
│   ├── Bin5L.fa
│   ├── Bin8L.fa
│   ├── ...
│   ╰── Bin444L.fa
├── reads
│   ╰── read.fastq.gz
├── assembler_circularity
│   ├── Read.p_ctg.noseq.gfa
│   ├── contigs.fasta.gz
│   ╰── assembly_info.txt
╰── gtdbtk
    ├── gtdbtk.bac120.summary.tsv
    ╰── gtdbtk.ar53.summary.tsv
```

## 3. Preparation

Complete the four preparation steps in [Usage](docs/USAGE.md#2-preparation). A prebuilt GTDB R232 database (131 GB) is available on Hugging Face as an alternative to a local database build. Use the same output directory for preparation and analysis.

## 4. Quick start

After preparation, run all analysis modules using end-to-end with the following options:

- `-a DIR`: Input directory of contig files (`.fa`, `.fasta`, or `.fna`), one per assembler. Contig IDs must be unique across assembly files, and filenames must be unique and free of whitespace.
- `-g DIR`: Input directory of MAG files (one file per MAG). MAG contig IDs must match those in the assembly files.
- `-r DIR`: Input directory of PacBio HiFi reads (`.fastq`, `.fq`, `.fastq.gz`, or `.fq.gz`); several files are pooled into `$Base/00_prep/hifi_reads.fastq[.gz]`. All assemblies specified with `-a` must be generated from the same HiFi read set.
- `-o DIR`: Output base directory.
- `-t N`: Number of threads (default: 80% of the available CPUs).
- `-m N`: Memory limit in GB (default: 80% of the available memory).
- `--skip MODULES`: Comma-separated modules to exclude (e.g., `--skip chimera,variants`).


With the [test dataset](#2-test-dataset):

```bash
Data=/path/to/input_dataset          # extracted test dataset
Base=/path/to/MAGverify_output

MAGverify end-to-end \
  -a "$Data/contigs" \
  -r "$Data/reads" \
  -g "$Data/MAGs" \
  -o "$Base" \
  -t 100 -m 900
```

**Note:** `end-to-end` runs the modules in the following order: `map` → `circularity` → `clipping` → `false-closure` → `phantom` → `repeats` → `false-duplication` → `variants` → `chimera` → `summary`. Before analysis, MAGverify checks for missing preparation outputs in `$Base/00_prep/` and identifies the necessary `prep-` commands. Interrupted analyses can be resumed by rerunning the same command, automatically skipping completed modules with unchanged inputs.

Each module can also be run separately. See [Usage](docs/USAGE.md#4-individual-modules) for the commands and [Module](docs/MODULES.md) for the detection criteria behind each call.


## 5. Main results

| File | Description |
| --- | --- |
| `summary.MAG.tsv` | MAG-level QC calls |
| `summary.contig.tsv` | Contig-level QC calls |
| `detail.MAG.tsv` | Detailed MAG results |
| `detail.contig.tsv` | Detailed contig results |

Module-specific tables are saved in numbered directories. Keep intermediate files and `.cache/` directories for resumption.

## 6. Documentation

| Document | Contents |
| --- | --- |
| [Usage](docs/USAGE.md) | Inputs, preparation and commands |
| [Module](docs/MODULES.md) | Detection methods and criteria |
| [Result](docs/OUTPUTS.md) | Output columns and interpretation |

## License

MIT. See [LICENSE](LICENSE).

## Citation

If you use MAGverify, please cite:

Song MJ, Yun T, Yoon S, Lee JY, Pinto AJ, Yoon S. HiFi long-read metagenomics reshapes genome-resolved interpretations of activated sludge microbiome. *bioRxiv* (2026). https://doi.org/10.64898/2026.11.18.767492
