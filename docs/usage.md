# mskilab-org/nf-jabba: Usage

# Introduction:

**mskilab-org/nf-JaBbA** is a new state-of-art bioinformatics pipeline from [`mskilab-org`](https://www.mskilab.org/) that is intended to run [`JaBbA`](https://github.com/mskilab-org/JaBbA/tree/master), an MIP based joint inference of copy number and rearrangement state in cancer whole genome sequence data. It runs all the pre-requisite modules necessary to run JaBbA and as followed in `mskilab-org`. This pipeline is built to handle only tumor-normal pairs as input (as of now) and is designed and tested to run on Human samples. 

We drew our inspiration and ideas from [`nf-core/Sarek`](https://github.com/nf-core/sarek), a workflow designed to detect variants on whole genome or targeted sequencing data. **`nf-jabba`** is built using [`Nextflow`](https://www.nextflow.io/) and is implemented using `Nextflow DSL2`. All the modules uses [`Docker`](https://www.docker.com/) and [`Singularity`](https://sylabs.io/docs/) containers which makes the pipeline easily reproducible and maintain its dependencies. Some of the modules/processes are used from [`nf-core/modules`](https://github.com/nf-core/modules) that are available for the Nextflow Community.

This pipeline has been designed to start from scratch using **FASTQ** files or start directly from **BAM** files as input and should be supplied in a **CSV** file (*please refer to the documentation below for the input format of the .csv file*). We incorporated a modified version of the `Alignment` step of `nf-JaBbA` pipeline from `nf-core/Sarek`, many thanks to the Sarek community. 

# Setting up a run: 

To run the pipeline from the beginning, first create an `--input` `samplesheet.csv` file with your file paths. A typical input should look like this:

```csv
patient,sex,status,sample,lane,fastq_1,fastq_2
TCXX49,XX,0,TCXX49_N,lane_1,/path/to/fastq_1.fq.gz,/path/to/fastq_2.gz
```
Each row represents a pair of fastq files (paired end) for each Sample.
A typical command for running the pipeline is as follows:

```bash
nextflow run mskilab-org/nf-jabba --input ./samplesheet.csv --outdir ./results --genome GATK.GRCh37 --tools <TOOLS> -profile singularity
```
This will launch the pipeline and run all the processes with the tools specified in `--tools` to `JaBbA`. It will create all the following files in the working directory from where the command is run:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> ⚠️ Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run mskilab-org/nf-jabba -profile singularity -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: './samplesheet.csv'
outdir: './results/'
genome: 'GATK.GRCh37'
tools:  'svaba,hetpileups,...,jabba'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

## Samplesheet input configurations (along with `--step`)

You need to create a samplesheet with information regarding the samples that you want to run the pipeline on. You need to specify the path of your **samplesheet** using the `--input` flag to specify the location. Make sure the input file is a *comma-separated* file and must contain headers with info discussed below. *It is highly recommended to provide the **absolute path** for inputs inside the samplesheet rather than relative paths.*

To mention a sample as paired tumor-normal, it has to be specified with the same `patient` ID, a different `sample`, and their respective `status`. For instance, a `tumor` sample should be mentioned **1** in `status` field for a sample, if it is normal mention **0**. If there are multiple `sample` IDs, `nf-jabba` will consider them as separate samples and output the results on separate folders based on `patient`, rest assured all the runs will be separate based on `patient`, so no need to be concerned with getting the outputs mixed.

```bash
--input '[path to samplesheet file]'
```

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


There are multiple `--steps` for `nf-jabba`. The main idea behind this was to make each of the tool separate so that one can only run their only desired tool rather than running the whole pipeline when it is provided with the required outputs. There are 2 primary `--steps` in the pipeline which can lead to `JaBbA` and run the module if mentioned with all the tools using a list of comma-separated names in `--tools`. A full input `csv` files for these 2 steps are shown below:
- **`--step alignment`**

```
patient,sex,status,sample,lane,fastq_1,fastq_2
TCXX49,XX,0,TCXX49_N,lane_1,/path/to/fastq_1.fq.gz,/path/to/fastq_2.gz
TCXX49,XX,0,TCXX49_N,lane_2,/path/to/fastq_1.fq.gz,/path/to/fastq_2.gz
TCXX49,XX,1,TCXX49_T,lane_1,/path/to/fastq_1.fq.gz,/path/to/fastq_2.gz
TCXX49,XX,1,TCXX49_T,lane_2,/path/to/fastq_1.fq.gz,/path/to/fastq_2.gz
TCXX52,NA,0,TCXX52_N,lane_1,/path/to/fastq_1.fq.gz,/path/to/fastq_2.gz
TCXX52,NA,1,TCXX52_T,lane_1,/path/to/fastq_1.fq.gz,/path/to/fastq_2.gz
```
- **`--step sv_calling`**

```
patient,sex,status,sample,bam,bai
TCXX49,XX,0,TCXX49_N,/path/to/alignment.bam,/path/to/alignment.bam.bai
TCXX49,XX,1,TCXX49_T,/path/to/alignment.bam,/path/to/alignment.bam.bai
TCXX52,NA,0,TCXX52_N,/path/to/alignment.bam,/path/to/alignment.bam.bai
TCXX52,NA,1,TCXX52_T,/path/to/alignment.bam,/path/to/alignment.bam.bai
```
> **Note**
> If you are using cram files, just replace *bam* and *bai* headers with *cram* and *crai* headers 
> and pass the *cram* and *crai* paths there


There are also many secondary `--steps` for this pipeline that are designed to only run a specific `--tools` for the module but is not adequate to run **`JaBbA`**. You can also run only **`JaBbA`** if you have all the inputs available and can be provided in the `--input` `csv` file. Below, we provide the other secondary steps and their desired input `csv` files for each step.

- **`--step markduplicates`**

```
patient,sex,status,sample,bam,bai
TCXX49,XX,0,TCXX49_N,/path/to/alignment.bam,/path/to/alignment.bam.bai
TCXX49,XX,1,TCXX49_T,/path/to/alignment.bam,/path/to/alignment.bam.bai
TCXX52,NA,0,TCXX52_N,/path/to/alignment.bam,/path/to/alignment.bam.bai
TCXX52,NA,1,TCXX52_T,/path/to/alignment.bam,/path/to/alignment.bam.bai
```
> **Note**
> If you are using cram files, just replace *bam* and *bai* headers with *cram* and *crai* headers 
> and pass the *cram* and *crai* paths there


- **`--step prepare_recalibration`**

```
patient,sex,status,sample,bam,bai
TCXX49,XX,0,TCXX49_N,/path/to/alignment.bam,/path/to/alignment.bam.bai
TCXX49,XX,1,TCXX49_T,/path/to/alignment.bam,/path/to/alignment.bam.bai
TCXX52,NA,0,TCXX52_N,/path/to/alignment.bam,/path/to/alignment.bam.bai
TCXX52,NA,1,TCXX52_T,/path/to/alignment.bam,/path/to/alignment.bam.bai
```
> **Note**
> If you are using cram files, just replace *bam* and *bai* headers with *cram* and *crai* headers 
> and pass the *cram* and *crai* paths there


- **`--step recalibrate`**

```
patient,sex,status,sample,bam,bai,table
TCXX49,XX,0,TCXX49_N,/path/to/alignment.bam,/path/to/alignment.bam.bai,TCXX49_N.table
TCXX49,XX,1,TCXX49_T,/path/to/alignment.bam,/path/to/alignment.bam.bai,TCXX49_T.table
TCXX52,NA,0,TCXX52_N,/path/to/alignment.bam,/path/to/alignment.bam.bai,TCXX52_N.table
TCXX52,NA,1,TCXX52_T,/path/to/alignment.bam,/path/to/alignment.bam.bai,TCXX52_T.table
```

- **`--step fragcounter`**

```
patient,sex,status,sample,bam,bai
TCXX49,XX,0,TCXX49_N,/path/to/alignment.bam,/path/to/alignment.bam.bai
TCXX49,XX,1,TCXX49_T,/path/to/alignment.bam,/path/to/alignment.bam.bai
TCXX52,NA,0,TCXX52_N,/path/to/alignment.bam,/path/to/alignment.bam.bai
TCXX52,NA,1,TCXX52_T,/path/to/alignment.bam,/path/to/alignment.bam.bai
```
> **Note**
> If you are using cram files, just replace *bam* and *bai* headers with *cram* and *crai* headers 
> and pass the *cram* and *crai* paths there


- **`--step dryclean`** (**Note**: you should also mention `--tools dryclean` to use Dryclean. This step also has `CBS`, if you want to perform both Dryclean and CBS, use `--tools dryclean,cbs`)

```
patient,sex,status,sample,cov
TCXX49,XX,0,TCXX49_N,/path/to/coverage.rds
TCXX49,XX,1,TCXX49_T,/path/to/coverage.rds
TCXX52,NA,0,TCXX52_N,/path/to/coverage.rds
TCXX52,NA,1,TCXX52_T,/path/to/coverage.rds
```

- **`--step hetpileups`** (**Note**: you should also mention `--tools hetpileups` to use HetPileups.)

```
patient,sex,status,sample,bam,bai
TCXX49,XX,0,TCXX49_N,/path/to/alignment.bam,/path/to/alignment.bam.bai
TCXX49,XX,1,TCXX49_T,/path/to/alignment.bam,/path/to/alignment.bam.bai
TCXX52,NA,0,TCXX52_N,/path/to/alignment.bam,/path/to/alignment.bam.bai
TCXX52,NA,1,TCXX52_T,/path/to/alignment.bam,/path/to/alignment.bam.bai
```
> **Note**
> HetPileups does not support for *`cram`* files, you must use *`bam`* files for this step.


- **`--step ascat`** (**Note**: you should also mention `--tools ascat` to use ASCAT.)

```
patient,sex,status,sample,cov,hets
TCXX49,XX,1,TCXX49_T,/path/to/coverage.rds,/path/to/hetpileups/sites.txt
TCXX52,NA,1,TCXX52_T,/path/to/coverage.rds,/path/to/hetpileups/sites.txt
```


- **`--step jabba`** (**Note**: you should also mention `--tools jabba` to use JaBbA.)

```
patient,sex,status,sample,cov,vcf
TCXX49,XX,1,TCXX49_T,/path/to/tumor/coverage.rds,/path/to/sv_caller/tumor/somatic.vcf
TCXX52,NA,1,TCXX52_T,/path/to/tumor/coverage.rds,/path/to/sv_caller/tumor/somatic.vcf
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull mskilab-org/nf-jabba
```

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

#### `-bg`
Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).


## Debugging any step/process:

To debug any step or process that failed, please check your current `execution_trace*.txt` file inside the `<outdir>/pipeline_info/` folder and gather the `hash` number for that process. Then go inside the `work` folder and paste that `hash` number to locate thw working directory for that process. There should be multiple `.command.*` files inside that folder which corresponds to your run. This includes log, sh, trace, error files. One good thing is you can run `.command.sh` script locally to check where it is exactly breaking and replicat the issue (though you might need to edit the command a bit to run it successfully locally).

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
