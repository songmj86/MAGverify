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

Use the test dataset to check that MAGverify runs correctly before applying it to your own data. [test_inputs.tar.gz](https://huggingface.co/datasets/ATCC19718/MAGverify-resources/resolve/main/test_inputs.tar.gz) (1.1 GB) is available on [Hugging Face](https://huggingface.co/datasets/ATCC19718/MAGverify-resources). Download and extract it:

```bash
Data=/path/to/test_inputs
mkdir -p "$Data"

aria2c -c -x 8 -s 8 -k 8M \
  -d "$Data" -o test_inputs.tar.gz \
  "https://huggingface.co/datasets/ATCC19718/MAGverify-resources/resolve/main/test_inputs.tar.gz"

tar -xzf "$Data/test_inputs.tar.gz" -C "$Data"
```

The dataset provides a set of test inputs for MAGverify including 3 assemblies, 10 MAGs , 1 HiFi read, 3 assembler circularity files and 2 GTDB-Tk results.

```
test_inputs
├── contigs
│   ├── metaMDBG.fa
│   ├── metaFlye.fa
│   ╰── hifiasm-meta.fa
├── MAGs
│   ├── Bin37L.fa
│   ├── Bin65L.fa
│   ├── ...
│   ╰── Bin443L.fa
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

Complete the four preparation steps in [Usage](docs/USAGE.md#2-preparation). Use the same output directory for preparation and analysis.

## 4. Quick start

After preparation, run all analysis modules using end-to-end with the following options:

- `-a DIR`: Input directory of contig files (`.fa`, `.fasta`, or `.fna`), one per assembler. Contig IDs must be unique across assembly files, and filenames must be unique and free of whitespace.
- `-g DIR`: Input directory of MAG files (one file per MAG). MAG contig IDs must match those in the assembly files.
- `-r DIR`: Input directory of PacBio HiFi reads (`.fastq`, `.fq`, `.fastq.gz`, or `.fq.gz`). Several files are pooled into `$Base/00_prep/hifi_reads.fastq[.gz]`. All assemblies specified with `-a` must be generated from the same HiFi read set.
- `-o DIR`: Output base directory.
- `-t N`: Number of threads (default: 80% of the available CPUs).
- `-m N`: Memory limit in GB (default: 80% of the available memory).
- `--skip MODULES`: Comma-separated modules to exclude (e.g., `--skip chimera,variants`).


With the [test dataset](#2-test-dataset):

```bash
Data=/path/to/test_inputs          # extracted test dataset
Base=/path/to/MAGverify_output

MAGverify end-to-end \
  -a "$Data/contigs" \
  -r "$Data/reads" \
  -g "$Data/MAGs" \
  -o "$Base" \
  -t 100 -m 900
```

**Note:** `end-to-end` executes modules in the following order: `map` → `circularity` → `clipping` → `false-closure` → `phantom` → `repeats` → `false-duplication` → `variants` → `chimera` → `summary`. Before execution, MAGverify checks for missing preparation files in `$Base/00_prep/` and reports the necessary `prep-` commands. To resume an interrupted analysis, rerun the same command. MAGverify automatically skips completed modules if their inputs remain unchanged.

Each module can also be run individually. See [Usage](docs/USAGE.md#4-individual-modules) for execution commands and [Modules](docs/MODULES.md) for module-specific methods and detection criteria.


## 5. Main results

MAGverify generates four integrated tables at the genome and contig levels in the root of the output directory. The `summary` tables report a `Y`, `N`, or `NA` flag for each analysis, whereas the `detail` tables provide additional module-specific diagnostic information.

| File | Description |
| --- | --- |
| `summary.MAG.tsv` | MAG-level summary results |
| `summary.contig.tsv` | Contig-level summary results |
| `detail.MAG.tsv` | MAG-level detailed results |
| `detail.contig.tsv` | Contig-level detailed results |

Module-specific tables are saved in numbered directories, listed in [Usage](docs/USAGE.md#5-output-files).

## 6. Documentation

| Document | Contents |
| --- | --- |
| [Usage](docs/USAGE.md) | Inputs, preparation and commands |
| [Modules](docs/MODULES.md) | Detection methods and criteria |
| [Result](docs/OUTPUTS.md) | Output columns and interpretation |

## License

MIT. See [LICENSE](LICENSE).

## Citation

If you use MAGverify, please cite:

Song MJ, Yun T, Yoon S, Lee JY, Pinto AJ, Yoon S. HiFi long-read metagenomics reshapes genome-resolved interpretations of activated sludge microbiome. *bioRxiv* (2026). https://doi.org/10.64898/2026.11.18.767492
