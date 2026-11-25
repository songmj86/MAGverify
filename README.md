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

Use the test dataset to check that MAGverify runs correctly before applying it to your own data. [input_dataset.tar.gz](https://huggingface.co/datasets/ATCC19718/MAGverify-resources/resolve/main/input_dataset.tar.gz) (1.1 GB) is available on [Hugging Face](https://huggingface.co/datasets/ATCC19718/MAGverify-resources) and holds one complete set of MAGverify inputs:

```
input_dataset
├── contigs                        # -a: three assemblies of the same sample
│   ├── metaMDBG.fa
│   ├── metaFlye.fa
│   ╰── hifiasm-meta.fa
├── MAGs                           # -g: 30 MAGs binned from those assemblies
│   ├── Bin5L.fa
│   ├── Bin8L.fa
│   ├── ...
│   ╰── Bin444L.fa
├── reads
│   ╰── read.fastq.gz              # -r: the HiFi reads the assemblies were built from
├── assembler_circularity          # prep-circularity -c
│   ├── Read.p_ctg.noseq.gfa       #     hifiasm-meta
│   ├── contigs.fasta.gz           #     metaMDBG
│   ╰── assembly_info.txt          #     metaFlye
╰── gtdbtk                         # prep-taxonomy -d
    ├── gtdbtk.bac120.summary.tsv
    ╰── gtdbtk.ar53.summary.tsv
```

Download and extract it:

```bash
Data=/path/to/input_dataset
mkdir -p "$Data"

aria2c -c -x 8 -s 8 -k 8M \
  -d "$Data" -o input_dataset.tar.gz \
  "https://huggingface.co/datasets/ATCC19718/MAGverify-resources/resolve/main/input_dataset.tar.gz"

tar -xzf "$Data/input_dataset.tar.gz" -C "$Data" --strip-components=1
```

## 3. Preparation

Run the preparation commands before analysis using the same `-o` to save their outputs in `$Base/00_prep/`. They are not included in `end-to-end`.

With the [test dataset](#2-test-dataset), the four preparation commands are:

```bash
Data=/path/to/input_dataset            # extracted test dataset
Base=/path/to/MAGverify_output
GTDB_DIR=/path/to/GTDB_R232_downloads

MAGverify prep-chimera -d "$GTDB_DIR" -o "$Base" -t 100

MAGverify prep-scg -a "$Data/contigs" -o "$Base" -t 100

MAGverify prep-circularity \
  -c "$Data/assembler_circularity/Read.p_ctg.noseq.gfa" \
  -c "$Data/assembler_circularity/contigs.fasta.gz" \
  -c "$Data/assembler_circularity/assembly_info.txt" \
  -g "$Data/MAGs" \
  -o "$Base"

MAGverify prep-taxonomy -d "$Data/gtdbtk" -g "$Data/MAGs" -o "$Base"
```

See [Usage](docs/USAGE.md#preparation) for the input formats, the download links and the prebuilt GTDB database that replaces `prep-chimera`.

## 4. Quick start

After preparation, run all analysis modules using end-to-end with the following options:

- `-a DIR`: Input directory of contig files (`.fa`, `.fasta`, or `.fna`), one per assembler. Contig IDs must be unique across assembly files.
- `-g DIR`: Input directory of MAG files (one file per MAG). MAG contig IDs must match those in the assembly files.
- `-r DIR`: Input directory of PacBio HiFi reads (`.fastq`, `.fq`, `.fastq.gz`, or `.fq.gz`) the assemblies were built from; several files are pooled into `$Base/00_prep/hifi_reads.fastq[.gz]`. Every assembly given with `-a` must come from this same read set, because these reads are mapped to all of them: withholding reads the assembly used produces clipping and phantom calls that are not real, and adding unrelated reads hides real ones.
- `-o DIR`: Output base directory. Use the same `-o` for all commands to reuse intermediate results.
- `-t N`: Number of threads (default: 80% of the available CPUs).
- `-m N`: Memory limit in GB for `map`, `circularity`, `repeats` and `variants` (default: 80% of the detected memory).
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

Each module can also be run separately. See [Usage](docs/USAGE.md#analysis-modules) for detailed descriptions and options.


## 5. Results

MAGverify generates two summary formats: `summary.MAG.tsv` and `summary.contig.tsv` provide `Y`/`N`/`NA` flags for each analysis, whereas `detail.MAG.tsv` and `detail.contig.tsv` additionally include detailed results from each module.

Each module saves its results in a separate directory, as listed below. Intermediate files are stored in hidden `.cache/` directories for reuse in subsequent runs. Keep these directories to enable analysis resumption.

| Module | Directory | Tables |
| --- | --- | --- |
| Preparation | `00_prep/` | `all_assemblies.fa`, `MAG_contig_map.tsv`, genes, `scg.contig.tsv`, `assembler_circularity.*.tsv`, `gtdbtk_taxonomy.MAG.tsv`, `taxonomy.gene.tsv` |
| Read mapping | `01_mapping/` | `mapping.assembler.tsv` (BAMs are in the same directory) |
| Circularity | `02_circularity/` | `circularity.MAG.tsv`, `circularity.contig.tsv` |
| Clipping | `03_clipping/` | `clipping.site.tsv`, `clipping.contig.tsv`, `clipping.MAG.tsv`, `clipping.assembler.tsv` |
| False closure | `04_false_closure/` | `false_closure.contig.tsv`, `false_closure.MAG.tsv` |
| Phantom sequences | `05_phantom/` | `phantom.region.tsv`, `phantom.contig.tsv`, `phantom.MAG.tsv` |
| Repeats | `06_repeats/` | `repeats.contig.tsv`, `repeats.MAG.tsv`, `repeats.region.tsv`, `repeats.HSP.tsv` |
| Taxonomic chimera | `07_chimera/` | `chimera.region.tsv`, `chimera.contig.tsv`, `chimera.MAG.tsv` |
| False duplication | `08_false_duplication/` | `false_duplication.pair.tsv`, `false_duplication.contig.tsv`, `false_duplication.MAG.tsv` |
| Assembly variants | `09_assembly_variants/` | `variants.site.tsv`, `variants.contig.tsv`, `variants.MAG.tsv` |
| Integrated summary | top of the output directory | `summary.MAG.tsv`, `summary.contig.tsv` (Y/N); `detail.MAG.tsv`, `detail.contig.tsv` (detailed) |

See [Result tables](docs/OUTPUTS.md) for detailed descriptions of all output columns.


## License

MIT. See [LICENSE](LICENSE).

## Citation

If you use MAGverify, please cite:

Song MJ, Yun T, Yoon S, Lee JY, Pinto AJ, Yoon S. HiFi long-read metagenomics reshapes genome-resolved interpretations of activated sludge microbiome. *bioRxiv* (2026). https://doi.org/10.64898/2026.11.18.767492
