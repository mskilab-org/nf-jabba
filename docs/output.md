# mskilab-org/nf-jabba: Output

## Introduction

This document describes the output produced by the pipeline. 

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- **Directory Structure**
- **Alignment**
- **SV Calling**
- **Coverages**
- **CBS**
- **ASCAT**
- **JaBbA**
- **Pipeline information**

## Directory Structure

```
{outdir}
├── csv
├── CBS
├── pipeline_info
├── Alignment
│   ├── markduplicates
│       └── <sample>
│   ├── recal_table
│       └── <sample>
│   └── recalibrated
│       └── <sample>
├── Reference
├── SV_calling
│   ├── SVABA
│       └── <sample>
│   ├── GRIDSS
│       └── <sample>
├── HetPileups
├── JaBbA
├── Coverages
│   ├── fragCounter_normal
│       └── <sample>
│   ├── fragCounter_tumor
│       └── <sample>
│   └── Dryclean_normal
│       └── <sample>
│   └── Dryclean_tumor
│       └── <sample>
└── Reports
    ├── <tool1>
    └── <tool2>
work/
.nextflow.log
```

## Alignment
`nf-jabba` pre-processes raw FastQ files or unmapped BAM files, based on [GATK best practices](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows).

### Preparation of input files (FastQ or (u)BAM)

[FastP](https://github.com/OpenGene/fastp) is a tool designed to provide all-in-one preprocessing for FastQ files and as such is used for trimming and splitting. By default, these files are not published. However, if publishing is enabled, please be aware that these files are only published once, meaning if trimming and splitting is enabled, then the resulting files will be sharded FastQ files with trimmed reads. If only one of them is enabled then the files contain either trimmed or split reads, respectively.

#### Trim adapters

[FastP](https://github.com/OpenGene/fastp) supports global trimming, which means it trims all reads in the front or the tail. This function is useful since sometimes you want to drop some cycles of a sequencing run. In the current implementation in Sarek
`--detect_adapter_for_pe` is set by default which enables auto-detection of adapter sequences. For more information on how to fine-tune adapter trimming, take a look into the parameter docs.

The resulting files are intermediate and by default not kept in the final files delivered to users. Set `--save_trimmed` to enable publishing of the files in:

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/Alignment/fastp/<sample>`**

- `<sample>_<lane>_{1,2}.fastp.fastq.gz>`
  - Bgzipped FastQ file

</details>

#### Split FastQ files

[FastP](https://github.com/OpenGene/fastp) supports splitting of one FastQ file into multiple files allowing parallel alignment of sharded FastQ file. To enable splitting, the number of reads per output can be specified. For more information, take a look into the parameter `--split_fastq`in the parameter docs.

These files are intermediate and by default not placed in the output-folder kept in the final files delivered to users. Set `--save_split` to enable publishing of these files to:

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/Alignment/fastp/<sample>/`**

- `<sample_lane_{1,2}.fastp.fastq.gz>`
  - Bgzipped FastQ file

</details>

### Mapping to the Reference Genome

#### BWA

[BWA](https://github.com/lh3/bwa) is a software package for mapping low-divergent sequences against a large reference genome. The aligned reads are then coordinate-sorted (or name-sorted if [`GATK MarkDuplicatesSpark`](https://gatk.broadinstitute.org/hc/en-us/articles/5358833264411-MarkDuplicatesSpark) is used for duplicate marking) with [samtools](https://www.htslib.org/doc/samtools.html).

#### BWA-mem2

[BWA-mem2](https://github.com/bwa-mem2/bwa-mem2) is a software package for mapping low-divergent sequences against a large reference genome.The aligned reads are then coordinate-sorted (or name-sorted if [`GATK MarkDuplicatesSpark`](https://gatk.broadinstitute.org/hc/en-us/articles/5358833264411-MarkDuplicatesSpark) is used for duplicate marking) with [samtools](https://www.htslib.org/doc/samtools.html).

<details markdown="1">
<summary>Output files for all mappers and samples</summary>

The alignment files (BAM or CRAM) produced by the chosen aligner are not published by default. CRAM output files will not be saved in the output-folder (`outdir`), unless the flag `--save_mapped` is used. BAM output can be selected by setting the flag `--save_output_as_bam`.

**Output directory: `{outdir}/Alignment/mapped/<sample>/`**

- if `--save_mapped`: `<sample>.sorted.cram` and `<sample>.sorted.cram.crai`

  - CRAM file and index

- if `--save_mapped --save_output_as_bam`: `<sample>.sorted.bam` and `<sample>.sorted.bam.bai`
  - BAM file and index
  </details>

### Mark Duplicates

During duplicate marking, read pairs that are likely to have originated from duplicates of the same original DNA fragments through some artificial processes are identified. These are considered to be non-independent observations, so all but a single read pair within each set of duplicates are marked, causing the marked pairs to be ignored by default during the variant discovery process.

For further reading and documentation see the [data pre-processing for variant discovery from the GATK best practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery).

The resulting CRAM files are delivered to the users.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/preprocessing/markduplicates/<sample>/`**

- `<sample>.md.cram` and `<sample>.md.cram.crai`
  - CRAM file and index
- if `--save_output_as_bam`:
  - `<sample>.md.bam` and `<sample>.md.bam.bai`

</details>

### Base Quality Score Recalibration

During Base Quality Score Recalibration, systematic errors in the base quality scores are corrected by applying machine learning to detect and correct for them. This is important for evaluating the correct call of a variant during the variant discovery process. However, this is not needed for all combinations of tools in Sarek. Notably, this should be turned off when having UMI tagged reads or using DragMap (see [here](https://gatk.broadinstitute.org/hc/en-us/articles/4407897446939--How-to-Run-germline-single-sample-short-variant-discovery-in-DRAGEN-mode)) as mapper.

For further reading and documentation see the [technical documentation by GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-).

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/Alignment/recal_table/<sample>/`**

- `<sample>.recal.table`
  - Recalibration table associated to the duplicates-marked CRAM file.

</details>

### GATK ApplyBQSR 

[GATK ApplyBQSR](https://gatk.broadinstitute.org/hc/en-us/articles/5358826654875-ApplyBQSR) recalibrates the base qualities of the input reads based on the recalibration table produced by the [GATK BaseRecalibrator](#gatk-baserecalibrator) tool.

The resulting recalibrated CRAM files are delivered to the user. Recalibrated CRAM files are usually 2-3 times larger than the duplicate-marked CRAM files.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/Alignment/recalibrated/<sample>/`**

- `<sample>.recal.cram` and `<sample>.recal.cram.crai`
  - CRAM file and index
- if `--save_output_as_bam`:
  - `<sample>.recal.bam` and `<sample>.recal.bam.bai` - BAM file and index
  </details>

## SV_calling

The results regarding structural variant calling are collected in {outdir}/SV_calling/. If some results from a variant caller do not appear here, please check out the `--tools` section to check if the an SV caller was mentioned.

Base Recalibrated CRAM files can used as an input to start the structural variant calling.

### SvABA
SvABA is a method for detecting structural variants in sequencing data using genome-wide local assembly. For reference, check [info](https://github.com/walaj/svaba)

### GRIDSS
GRIDSS is a module software suite containing tools useful for the detection of genomic rearrangements. It can also detect purity and ploidy. For reference, check [info](https://github.com/PapenfussLab/gridss)

## fragCounter
The goal of fragCounter is to correct Whole genome or targeted sequencing data for GC and mappability bias.
The GC bias curve is determined by loess regression of read count by GC and mappability scores. For reference, check [info](https://github.com/mskilab-org/fragCounter)

## Dryclean
Dryclean is a robust principal component analysis (rPCA) based method. It uses a panel of normal (PON) samples to learn the landscape of both biological and technical noise in read depth data. Dryclean then uses this landscape to significantly reduce noise and artifacts in the signal for tumor samples. The input to the algorithm is a GenomicsRanges object containing read depth. 
For reference, check [info](https://github.com/mskilab-org/dryclean)

## CBS
Segmentation is done by circular binary segmentation (CBS) algorithm after getting tumor/normal ratios of corrected read counts. We use a custom module script, check [here](../bin/cbsFH.R)

## HetPileups
Pileup mutational calls are done using a custom module script called HetPileups. We use a custom module script, check [here](../bin/Pileups.R)

## ASCAT
 ASCAT(allele-specific copy number analysis of tumors) is used to accurately dissect the allele-specific copy number of solid tumors, simultaneously estimating and adjusting for both tumor ploidy and nonaberrant cell admixture. We use ASCAT ploidy to supply for JaBbA. For more info regarding ASCAT, check [here](https://github.com/VanLoo-lab/ascat)

## JaBbA
JaBbA builds a genome graph based on junctions and read depth from whole genome sequencing, inferring optimal copy numbers for both vertices (DNA segments) and edges (bonds between segments). It can be used for discovering various patterns of structural variations. For more info regarding JaBbA, check [here](https://github.com/mskilab-org/JaBbA)

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

![MultiQC - FastQC sequence counts plot](images/mqc_fastqc_counts.png)

![MultiQC - FastQC mean quality scores plot](images/mqc_fastqc_quality.png)

![MultiQC - FastQC adapter content plot](images/mqc_fastqc_adapter.png)


### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
