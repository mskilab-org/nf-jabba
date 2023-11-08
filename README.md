# NF-JaBbA (Nextflow - Junction Balance Analysis Pipeline)
```
▐▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▌
▐                                                                                                  ▌
▐                  ██████                   █████           ███████████  █████       █████████     ▌
▐                 ███░░███                 ░░███           ░░███░░░░░███░░███       ███░░░░░███    ▌
▐    ████████    ░███ ░░░                   ░███   ██████   ░███    ░███ ░███████  ░███    ░███    ▌
▐   ░░███░░███  ███████    ██████████       ░███  ░░░░░███  ░██████████  ░███░░███ ░███████████    ▌
▐    ░███ ░███ ░░░███░    ░░░░░░░░░░        ░███   ███████  ░███░░░░░███ ░███ ░███ ░███░░░░░███    ▌
▐    ░███ ░███   ░███                 ███   ░███  ███░░███  ░███    ░███ ░███ ░███ ░███    ░███    ▌
▐    ████ █████  █████               ░░████████  ░░████████ ███████████  ████████  █████   █████   ▌
▐   ░░░░ ░░░░░  ░░░░░                 ░░░░░░░░    ░░░░░░░░ ░░░░░░░░░░░  ░░░░░░░░  ░░░░░   ░░░░░    ▌
▐                                                                                                  ▌
▐▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▌
```

[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/mskilab-org/nf-jabba)

## Introduction

**mskilab-org/nf-JaBbA** is a new state-of-art bioinformatics pipeline from [`mskilab-org`](https://www.mskilab.org/) that is intended to run [`JaBbA`](https://github.com/mskilab-org/JaBbA/tree/master), an MIP based joint inference of copy number and rearrangement state in cancer whole genome sequence data. It runs all the pre-requisite modules necessary to run JaBbA and as followed in `mskilab-org`. This pipeline is built to handle only tumor-normal pairs as input (as of now) and is designed and tested to run on Human samples. 

We drew our inspiration and ideas from [`nf-core/Sarek`](https://github.com/nf-core/sarek), a workflow designed to detect variants on whole genome or targeted sequencing data. It is built using [`Nextflow`](https://www.nextflow.io/) and is implemented using `Nextflow DSL2`. All the modules uses [`Docker`](https://www.docker.com/) and [`Singularity`](https://sylabs.io/docs/) containers which makes the pipeline easily reproducible and maintain its dependencies. Some of the modules/processes are used from [`nf-core/modules`](https://github.com/nf-core/modules) that are available for the Nextflow Community.

This pipeline has been designed to start from scratch using **FASTQ** files or start directly from **BAM** files as input and should be supplied in a **CSV** file (*please refer to the documentation below for the input format of the .csv file*). We incorporated a modified version of the `Alignment` step of `nf-JaBbA` pipeline from `nf-core/Sarek`, many thanks to the Sarek community. 

## Workflow Summary:
1. Alignment to Reference Genome (currently support `BWA-MEM` & `BWA-MEM2`)
2. Quality Control (using [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
3. Perform trimming (must turn on using `--trim_fastq`) (using `fastp`)
4. Marking Duplicates (using `GATK MarkDuplicates`)
5. Perform baserecalibration (using `GATK BaseRecalibrator`)
6. Apply BQSR (using `GATK ApplyBQSR`)
7. Perform Structural Variants Calling (using [`SVABA`](https://github.com/walaj/svaba) and/or [`GRIDSS`](https://github.com/PapenfussLab/gridss); must mention using `--tools`)
8. Perform Pileups (using mskilab's custom `HetPileups`; must mention using `--tools`)
9. Generate raw coverages and corect for GC & Mappability bias (using [`fragCounter`](https://github.com/mskilab-org/fragCounter); must mention using `--tools`)
10. Remove biological and technical noise from coverage data. (using [`Dryclean`](https://github.com/mskilab-org/dryclean); must mention using `--tools`)
11. Perform Segmentation by using tumor/normal ratios of corrected read counts, (using `CBS` circular binary segmentation algorithm; must mention using `--tools`)
12. Get Purity & Ploidy separately to supply to JaBbA (currently support [`ASCAT`](https://www.crick.ac.uk/research/labs/peter-van-loo/software) to pass ploidy values to JaBbA; must mention using `--tools`)
13. Execute JaBbA (using inputs from `Dryclean`, `CBS`, `HetPileups` and/or `ASCAT`; must mention using `--tools`)


## Usage

> **Note**
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
> to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
> with `-profile test` before running the workflow on actual data.

### Setting up the ***samplesheet.csv*** file for input:

You need to create a samplesheet with information regarding the samples that you want to run the pipeline on. You need to specify the path of your **samplesheet** using the `--input` flag to specify the location. Make sure the input file is a *comma-separated* file and must contain headers with info discussed below. *It is highly recommended to provide the **absolute path** for inputs inside the samplesheet rather than relative paths.*

To mention a sample as paired tumor-normal, it has to be specified with the same `patient` ID, a different `sample`, and their respective `status`. For instance, a `tumor` sample should be mentioned **1** in `status` field for a sample, if it is normal mention **0**. If there are multiple `sample` IDs, `nf-jabba` will consider them as separate samples and output the results on separate folders based on `patient`, rest assured all the runs will be separate based on `patient`, so no need to be concerned with getting the outputs mixed.

You need to specify the desired output directory path using `--outdir` flag when you start a run so that the outputs get stored on your designated folder and separated by `tool` and `sample` names in folders.


<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

-->

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run mskilab-org/nf-jabba \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> **Warning:**
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
> provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

## Credits

`nf-jabba` was originally written by [`Tanubrata Dey`](https://github.com/tanubrata) and [`Shihab Dider`](https://github.com/shihabdider) at the Perlmutter Cancer Center and the New York Genome Center.

We thank the following people for their extensive guidance in the development of this pipeline:
- [Marcin Imielinski](https://github.com/imielinski)
- [Joel Rosiene](https://github.com/jrosiene)

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  mskilab-org/nf-jabba for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
