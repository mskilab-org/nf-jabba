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

This pipeline is built after being influenced by `nf-core/Sarek`, a workflow designed to detect variants on whole genome or targeted sequencing data. It is built using [`Nextflow`](https://www.nextflow.io/) and is implemented using `Nextflow DSL2`. All the modules uses [`Docker`](https://www.docker.com/) and [`Singularity`](https://sylabs.io/docs/) containers which makes the pipeline easily reproducible and maintain its dependencies. Some of the modules/processes are used from [`nf-core/modules`](https://github.com/nf-core/modules) that are available for the Nextflow Community.

This pipeline has been designed to start from scratch using **fastq** files or start using **BAM** files and should be supplied in a **csv** file as input (*please refer the documentation below for the input format of the .csv file*). 

<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

## Usage

> **Note**
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
> to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
> with `-profile test` before running the workflow on actual data.

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

mskilab-org/nf-jabba was originally written by Tanubrata Dey and Shihab Dider.

We thank the following people for their extensive assistance in the development of this pipeline:

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
