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

We drew our inspiration and ideas from [`nf-core/Sarek`](https://github.com/nf-core/sarek), a workflow designed to detect variants on whole genome or targeted sequencing data. **`nf-jabba`** is built using [`Nextflow`](https://www.nextflow.io/) and is implemented using `Nextflow DSL2`. All the modules uses [`Docker`](https://www.docker.com/) and [`Singularity`](https://sylabs.io/docs/) containers which makes the pipeline easily reproducible and maintain its dependencies. Some of the modules/processes are used from [`nf-core/modules`](https://github.com/nf-core/modules) that are available for the Nextflow Community.

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

To run the pipeline from the beginning, first create an `--input` `sampleSheet.csv` file with your file paths. A typical input should look like this:

```csv
patient,sex,status,sample,lane,fastq_1,fastq_2
TCXX49,XX,0,TCXX49_N,lane_1,/path/to/fastq_1.fq.gz,/path/to/fastq_2.gz
```
Each row represents a pair of fastq files (paired end) for each Sample.
After the input file is ready, you can run the pipeline using:

```bash
nextflow run mskilab-org/nf-jabba \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --tools <svaba,fragcounter,dryclean,cbs,hetpileups,ascat,jabba> \
   --genome <GATK.GRCh37/GATK.GRCh38>
```
> **Warning:**
> Please provide pipeline parameters via the CLI or Nextflow [`-params-file`](https://www.nextflow.io/blog/2020/cli-docs-release.html) option. Custom config files including those
> provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).


### Discussion of expected fields in input file and expected inputs for each `--step`

A typical sample sheet should populate with the column names as shown below:

|   Column Name   |                                               Description                                                                                                 |
|-----------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------|
|   patient       | Patient or Sample ID. This should differentiate each patient/sample. *Note*: Each patient can have multiple sample names.                                 |
|   sample        | Sample ID for each Patient. Should differentiate between tumor and normal. Sample IDs should be unique to Patient IDs                                     |
|   lane          | If starting with FASTQ files, if there are multiple lanes for each sample for each patient, mention lane name. **Required for `--step alignment`.         | 
|   sex           | If known, please provide the sex for the patient. For instance if **Male** type XY, else if **Female** type XX, else others and unknown should be NA.     |
|   status        | This should tell if your sample is **tumor** or **normal**. For **normal**, write 0, and for **tumor**, write 1.                                          |
|   fastq_1       | Full Path to FASTQ file read 1. The extension should be `.fastq.gz` or `.fq.gz`. **Required** for `--step alignment`.                                     |
|   fastq_2       | Full Path to FASTQ file read 2. The extension should be `.fastq.gz` or `.fq.gz`. **Required** for `--step alignment`.                                     |
|   bam           | Full Path to BAM file. The extension should be `.bam`. **Required** for `--step sv_calling`.                                                              |
|   bai           | Full Path to BAM index file. The extension should be `.bam.bai`. **Required** for `--step sv_calling`.                                                    |
|   cram          | Full Path to CRAM file. The extension should be `.cram`. **Required** for `--step sv_calling` if file is of type `CRAM`.                                  |
|   crai          | Full Path to CRAM index file. The extension should be `.cram.crai`. **Required** for `--step sv_calling` if file is of type `CRAM`.                       |
|   table         | Full path to Recalibration table file. **Required** for `--step recalibrate`.                                                                             |
|   vcf           | Full path to VCF file. **Required** for `--step jabba`.                                                                                                   |
|   hets          | Full path to HetPileups .txt file. **Required** for `--step jabba`.                                                                                       |


For more information  and further functionality regarding the pipeline usage and inputs necesaary for each step please follow the [Usage](docs//usage.md) documentation as suggested here.

### Helpful Core Nextflow Commands:

#### `-resume`
This is a life saving command which is part of Nextflow. If a Process of the pipeline fails at some point, Nextflow has the ability to start from that step where the job failed rather than starting all the way from the beginning. You must specify this in the `-CLI` or on the `command-line` when restarting a pipeline. You can also supply a run name to resume a specific run using: `-resume` [run-name]. Use the `nextflow log` command to show previous run names.

#### `-profile`
Use this parameter for choosing a configuration profile. Profiles can give configuration presets for different computing environments.

Several generic profiles have been provided with the pipeline which instruct the pipeline to use software packaged using different methods. You need to use this option to mention when using containers (singularity/Docker) which is highly recommended for running the pipeline. 

#### `-c`
You can mention custom configuration scripts to run the pipeline with using `-c` flag and providing the path to the `.config` file. This is advised when you want to submit processes into an executor like `slurm/LSF/..`.

#### `-bg`
The Nextflow `-bg` flag helps launching Nextflow pipeline in the background, and being detached from your terminal so that the curren run does not stop if you log out of your session and the log of the run are saved inside a file. Alternative ways include using `screen` or `tmux` sessions which you can easily detach and log back in at a later time.

## Debugging any step/process:

To debug any step or process that failed, please check your current `execution_trace*.txt` file inside the `<outdir>/pipeline_info/` folder and gather the `hash` number for that process. Then go inside the `work` folder and paste that `hash` number to locate thw working directory for that process. There should be multiple `.command.*` files inside that folder which corresponds to your run. This includes log, sh, trace, error files. One good thing is you can run `.command.sh` script locally to check where it is exactly breaking and replicat the issue (though you might need to edit the command a bit to run it successfully locally).

## Credits

`nf-jabba` was originally written by [`Tanubrata Dey`](https://github.com/tanubrata) and [`Shihab Dider`](https://github.com/shihabdider) at the Perlmutter Cancer Center and the New York Genome Center.

We thank the following people for their extensive guidance in the development of this pipeline:
- [Marcin Imielinski](https://github.com/imielinski)
- [Joel Rosiene](https://github.com/jrosiene)


## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
