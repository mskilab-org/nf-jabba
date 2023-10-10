/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation
// log.info paramsSummaryLog(workflow)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.ascat_alleles,
    params.ascat_loci,
    params.ascat_loci_gc,
    params.ascat_loci_rt,
    params.bwa,
    params.bwamem2,
    params.cf_chrom_len,
    params.chr_dir,
    //params.cnvkit_reference,
    params.dbnsfp,
    params.dbnsfp_tbi,
    params.dbsnp,
    params.dbsnp_tbi,
    params.dict,
    //params.dragmap,
    params.fasta,
    params.fasta_fai,
    params.germline_resource,
    params.germline_resource_tbi,
    params.input,
    params.intervals,
    params.known_indels,
    params.known_indels_tbi,
    params.known_snps,
    params.known_snps_tbi,
    params.mappability,
    params.multiqc_config,
    params.pon,
    params.pon_tbi,
    params.spliceai_indel,
    params.spliceai_indel_tbi,
    params.spliceai_snv,
    params.spliceai_snv_tbi
]



// Validating the input parameters
WorkflowHeisenbio.initialise(params, log)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Check mandatory parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
println "This is the output of build_only_index: $params.build_only_index"
for (param in checkPathParamList) if (param) file(param, checkIfExists: true)


// Set input, can either be from --input or from automatic retrieval in WorkflowHeisenbio.groovy
if (params.input) {
    ch_from_samplesheet = params.build_only_index ? Channel.empty() : Channel.fromSamplesheet("input")
} else {
    ch_from_samplesheet = params.build_only_index ? Channel.empty() : Channel.fromSamplesheet("input_restart")
}

input_sample = ch_from_samplesheet
        .map{ meta, fastq_1, fastq_2, table, cram, crai, bam, bai, vcf, variantcaller ->
            // generate patient_sample key to group lanes together
            [ meta.patient + meta.sample, [meta, fastq_1, fastq_2, table, cram, crai, bam, bai, vcf, variantcaller] ]
        }

        .tap{ ch_with_patient_sample }                                         // save the channel
        .groupTuple()                                                          //group by patient_sample to get all lanes
        .map { patient_sample, ch_items ->
                                                                               // get number of lanes per sample
            [ patient_sample, ch_items.size() ]
        }
        .combine(ch_with_patient_sample, by: 0)                                // for each entry add numLanes

        .map { patient_sample, num_lanes, ch_items ->

            (meta, fastq_1, fastq_2, table, cram, crai, bam, bai, vcf, variantcaller) = ch_items

            if (meta.lane && fastq_2) {
                meta           = meta + [id: "${meta.sample}-${meta.lane}".toString()]
                def CN         = params.seq_center ? "CN:${params.seq_center}\\t" : ''

                def flowcell   = flowcellLaneFromFastq(fastq_1)
                // Don't use a random element for ID, it breaks resuming
                def read_group = "\"@RG\\tID:${flowcell}.${meta.sample}.${meta.lane}\\t${CN}PU:${meta.lane}\\tSM:${meta.patient}_${meta.sample}\\tLB:${meta.sample}\\tDS:${params.fasta}\\tPL:${params.seq_platform}\""

                meta           = meta - meta.subMap('lane') + [num_lanes: num_lanes.toInteger(), read_group: read_group.toString(), data_type: 'fastq', size: 1]

                if (params.step == 'alignment') return [ meta, [ fastq_1, fastq_2 ] ]
                else {
                    error("Samplesheet contains fastq files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations")
                }

                                                                                 // start from BAM if BAM files are provided
            } else if (meta.lane && bam) {
                if (params.step != 'alignment' && !bai) {
                    error("BAM index (bai) should be provided.")
                }
                meta            = meta + [id: "${meta.sample}-${meta.lane}".toString()]
                def CN          = params.seq_center ? "CN:${params.seq_center}\\t" : ''
                def read_group  = "\"@RG\\tID:${meta.sample}_${meta.lane}\\t${CN}PU:${meta.lane}\\tSM:${meta.patient}_${meta.sample}\\tLB:${meta.sample}\\tDS:${params.fasta}\\tPL:${params.seq_platform}\""

                meta            = meta - meta.subMap('lane') + [num_lanes: num_lanes.toInteger(), read_group: read_group.toString(), data_type: 'bam', size: 1]

                if (params.step != 'annotate') return [ meta - meta.subMap('lane'), bam, bai ]
                else {
                    error("Samplesheet contains bam files but step is `annotate`. The pipeline is expecting vcf files for the annotation. Please check your samplesheet or adjust the step parameter.")
                }

                                                                                  // start from recalibration step
            } else if (table && cram) {
                meta = meta + [id: meta.sample, data_type: 'cram']

                if (!(params.step == 'alignment' || params.step == 'annotate')) return [ meta - meta.subMap('lane'), cram, crai, table ]
                else {
                    error("Samplesheet contains cram files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.")
                }

                                                                                  // starting from recalibration step  when skipping MarkDuplicates
            } else if (table && bam) {
                meta = meta + [id: meta.sample, data_type: 'bam']

                if (!(params.step == 'alignment' || params.step == 'annotate')) return [ meta - meta.subMap('lane'), bam, bai, table ]
                else {
                    error("Samplesheet contains bam files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.")
                }

                                                                                 // when starting from prepare_recalibration or variant_calling
            } else if (cram) {
                meta = meta + [id: meta.sample, data_type: 'cram']

                if (!(params.step == 'alignment' || params.step == 'annotate')) return [ meta - meta.subMap('lane'), cram, crai ]
                else {
                    error("Samplesheet contains cram files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.")
                } 

                                                                                // when starting from prepare_recalibration when skipping MarkDuplicates or `--step markduplicates`
            } else if (bam) {
                meta = meta + [id: meta.sample, data_type: 'bam']

                if (!(params.step == 'alignment' || params.step == 'annotate')) return [ meta - meta.subMap('lane'), bam, bai ]
                else {
                    error("Samplesheet contains bam files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.")
                }

                                                                               // when starting from annotation step
            } else if (vcf) {
                meta = meta + [id: meta.sample, data_type: 'vcf', variantcaller: variantcaller ?: '']

                if (params.step == 'annotate') return [ meta - meta.subMap('lane'), vcf ]
                else {
                    error("Samplesheet contains vcf files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.")
                }
            } else {
                error("Missing or unknown field in csv file header. Please check your samplesheet")
            }
        }


//skipping some lines from sarek that does checks for some downstream steps, not focussing on them now...


// Fails when missing resources for baserecalibrator
// Warns when missing resources for haplotypecaller
if (!params.dbsnp && !params.known_indels) {
    if (params.step in ['alignment', 'markduplicates', 'prepare_recalibration', 'recalibrate'] && (!params.skip_tools || (params.skip_tools && !params.skip_tools.split(',').contains('baserecalibrator')))) {
        error("Base quality score recalibration requires at least one resource file. Please provide at least one of `--dbsnp` or `--known_indels`\nYou can skip this step in the workflow by adding `--skip_tools baserecalibrator` to the command.")
    }
    if (params.tools && (params.tools.split(',').contains('haplotypecaller') || params.tools.split(',').contains('sentieon_haplotyper'))) {
        log.warn "If GATK's Haplotypecaller or Sentieon's Haplotyper is specified, without `--dbsnp` or `--known_indels no filtering will be done. For filtering, please provide at least one of `--dbsnp` or `--known_indels`.\nFor more information see FilterVariantTranches (single-sample, default): https://gatk.broadinstitute.org/hc/en-us/articles/5358928898971-FilterVariantTranches\nFor more information see VariantRecalibration (--joint_germline): https://gatk.broadinstitute.org/hc/en-us/articles/5358906115227-VariantRecalibrator\nFor more information on GATK Best practice germline variant calling: https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-"
    }
}

if (params.joint_germline && (!params.tools || !(params.tools.split(',').contains('haplotypecaller') || params.tools.split(',').contains('sentieon_haplotyper')))) {
    error("The GATK's Haplotypecaller or Sentieon's Haplotyper should be specified as one of the tools when doing joint germline variant calling.) ")
}

if (params.joint_germline && (!params.dbsnp || !params.known_indels || !params.known_snps || params.no_intervals)) {
    log.warn "If GATK's Haplotypecaller or Sentieon's Haplotyper is specified, without `--dbsnp`, `--known_snps`, `--known_indels` or the associated resource labels (ie `known_snps_vqsr`), no variant recalibration will be done. For recalibration you must provide all of these resources.\nFor more information see VariantRecalibration: https://gatk.broadinstitute.org/hc/en-us/articles/5358906115227-VariantRecalibrator \nJoint germline variant calling also requires intervals in order to genotype the samples. As a result, if `--no_intervals` is set to `true` the joint germline variant calling will not be performed."
}


//skipping some regarding mutect2 for now and ascat and variant calling and annotate

if ((params.download_cache) && (params.snpeff_cache || params.vep_cache)) {
    error("Please specify either `--download_cache` or `--snpeff_cache`, `--vep_cache`")
}














/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
//include { INPUT_CHECK } from '../subworkflows/local/input_check'

// Initialize file channels based on params, defined in the params.genomes[params.genome] scope (commenting out some for now for downstream steps)

ascat_alleles      = params.ascat_alleles      ? Channel.fromPath(params.ascat_alleles).collect()     : Channel.empty()
ascat_loci         = params.ascat_loci         ? Channel.fromPath(params.ascat_loci).collect()        : Channel.empty()
ascat_loci_gc      = params.ascat_loci_gc      ? Channel.fromPath(params.ascat_loci_gc).collect()     : Channel.value([])
ascat_loci_rt      = params.ascat_loci_rt      ? Channel.fromPath(params.ascat_loci_rt).collect()     : Channel.value([])
cf_chrom_len       = params.cf_chrom_len       ? Channel.fromPath(params.cf_chrom_len).collect()      : []
chr_dir            = params.chr_dir            ? Channel.fromPath(params.chr_dir).collect()           : Channel.value([])
dbsnp              = params.dbsnp              ? Channel.fromPath(params.dbsnp).collect()             : Channel.value([])
fasta              = params.fasta              ? Channel.fromPath(params.fasta).first()               : Channel.empty()
fasta_fai          = params.fasta_fai          ? Channel.fromPath(params.fasta_fai).collect()         : Channel.empty()
known_indels       = params.known_indels       ? Channel.fromPath(params.known_indels).collect()      : Channel.value([])
known_snps         = params.known_snps         ? Channel.fromPath(params.known_snps).collect()        : Channel.value([])

germline_resource  = params.germline_resource  ? Channel.fromPath(params.germline_resource).collect() : Channel.value([]) // Mutect2 does not require a germline resource, so set to optional input
mappability        = params.mappability        ? Channel.fromPath(params.mappability).collect()       : Channel.value([])
pon                = params.pon                ? Channel.fromPath(params.pon).collect()               : Channel.value([]) // PON is optional for Mutect2 (but highly recommended)
ascat_genome       = params.ascat_genome       ?: Channel.empty()
dbsnp_vqsr         = params.dbsnp_vqsr         ? Channel.value(params.dbsnp_vqsr) : Channel.empty()
known_indels_vqsr  = params.known_indels_vqsr  ? Channel.value(params.known_indels_vqsr) : Channel.empty()
known_snps_vqsr    = params.known_snps_vqsr    ? Channel.value(params.known_snps_vqsr) : Channel.empty()
snpeff_db          = params.snpeff_db          ?: Channel.empty()
vep_cache_version  = params.vep_cache_version  ?: Channel.empty()
vep_genome         = params.vep_genome         ?: Channel.empty()
vep_species        = params.vep_species        ?: Channel.empty()





// Initialize files channels based on params, not defined within the params.genomes[params.genome] scope
if (params.snpeff_cache && params.tools && params.tools.contains("snpeff")) {
    def snpeff_annotation_cache_key = params.use_annotation_cache_keys ? "${params.snpeff_genome}.${params.snpeff_db}/" : ""
    def snpeff_cache_dir =  "${snpeff_annotation_cache_key}${params.snpeff_genome}.${params.snpeff_db}"
    def snpeff_cache_path_full = file("$params.snpeff_cache/$snpeff_cache_dir", type: 'dir')
    if ( !snpeff_cache_path_full.exists() || !snpeff_cache_path_full.isDirectory() ) {
        error("Files within --snpeff_cache invalid. Make sure there is a directory named ${snpeff_cache_dir} in ${params.snpeff_cache}.\n")
    }
    snpeff_cache = Channel.fromPath(file("${params.snpeff_cache}/${snpeff_annotation_cache_key}"), checkIfExists: true).collect()
        .map{ cache -> [ [ id:"${params.snpeff_genome}.${params.snpeff_db}" ], cache ] }
} else snpeff_cache = []

if (params.vep_cache && params.tools && params.tools.contains("vep")) {
    def vep_annotation_cache_key = params.use_annotation_cache_keys ? "${params.vep_cache_version}_${params.vep_genome}/" : ""
    def vep_cache_dir = "${vep_annotation_cache_key}${params.vep_species}/${params.vep_cache_version}_${params.vep_genome}"
    def vep_cache_path_full = file("$params.vep_cache/$vep_cache_dir", type: 'dir')
    if ( !vep_cache_path_full.exists() || !vep_cache_path_full.isDirectory() ) {
        error("Files within --vep_cache invalid. Make sure there is a directory named ${vep_cache_dir} in ${params.vep_cache}.\n")
    }
    vep_cache = Channel.fromPath(file("${params.vep_cache}/${vep_annotation_cache_key}"), checkIfExists: true).collect()
} else vep_cache = []

vep_extra_files = []

//if (params.dbnsfp && params.dbnsfp_tbi) {
//    vep_extra_files.add(file(params.dbnsfp, checkIfExists: true))
//    vep_extra_files.add(file(params.dbnsfp_tbi, checkIfExists: true))
//}

//if (params.spliceai_snv && params.spliceai_snv_tbi && params.spliceai_indel && params.spliceai_indel_tbi) {
//    vep_extra_files.add(file(params.spliceai_indel, checkIfExists: true))
//    vep_extra_files.add(file(params.spliceai_indel_tbi, checkIfExists: true))
//    vep_extra_files.add(file(params.spliceai_snv, checkIfExists: true))
//    vep_extra_files.add(file(params.spliceai_snv_tbi, checkIfExists: true))
//}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL/NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Create samplesheets to restart from different steps
include { CHANNEL_ALIGN_CREATE_CSV                    } from '../subworkflows/local/channel_align_create_csv/main'
include { CHANNEL_MARKDUPLICATES_CREATE_CSV           } from '../subworkflows/local/channel_markduplicates_create_csv/main'
include { CHANNEL_BASERECALIBRATOR_CREATE_CSV         } from '../subworkflows/local/channel_baserecalibrator_create_csv/main'
include { CHANNEL_APPLYBQSR_CREATE_CSV                } from '../subworkflows/local/channel_applybqsr_create_csv/main'
//include { CHANNEL_VARIANT_CALLING_CREATE_CSV          } from '../subworkflows/local/channel_variant_calling_create_csv/main'


// Download annotation cache if needed
include { PREPARE_CACHE                               } from '../subworkflows/local/prepare_cache/main'


// Build indices if needed
include { PREPARE_GENOME                              } from '../subworkflows/local/prepare_genome/main'


// Build intervals if needed
include { PREPARE_INTERVALS                           } from '../subworkflows/local/prepare_intervals/main'


// Convert BAM files to FASTQ files
include { BAM_CONVERT_SAMTOOLS as CONVERT_FASTQ_INPUT } from '../subworkflows/local/bam_convert_samtools/main'
//include { BAM_CONVERT_SAMTOOLS as CONVERT_FASTQ_UMI   } from '../subworkflows/local/bam_convert_samtools/main'


// Map input reads to reference genome
include { FASTQ_ALIGN_BWAMEM_MEM2    } from '../subworkflows/local/fastq_align_bwamem_mem2/main'


// Merge and index BAM files (optional)
include { BAM_MERGE_INDEX_SAMTOOLS                    } from '../subworkflows/local/bam_merge_index_samtools/main'


// Convert BAM files
include { SAMTOOLS_CONVERT as BAM_TO_CRAM             } from '../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_CONVERT as BAM_TO_CRAM_MAPPING     } from '../modules/nf-core/samtools/convert/main'

// Convert CRAM files (optional)
include { SAMTOOLS_CONVERT as CRAM_TO_BAM             } from '../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_CONVERT as CRAM_TO_BAM_RECAL       } from '../modules/nf-core/samtools/convert/main'


// Mark Duplicates (+QC)
include { BAM_MARKDUPLICATES                          } from '../subworkflows/local/bam_markduplicates/main'
//include { BAM_SENTIEON_DEDUP                          } from '../subworkflows/local/bam_sentieon_dedup/main'


// QC on CRAM
include { CRAM_QC_MOSDEPTH_SAMTOOLS as CRAM_QC_NO_MD  } from '../subworkflows/local/cram_qc_mosdepth_samtools/main'
include { CRAM_QC_MOSDEPTH_SAMTOOLS as CRAM_QC_RECAL  } from '../subworkflows/local/cram_qc_mosdepth_samtools/main'


// Create recalibration tables
include { BAM_BASERECALIBRATOR                        } from '../subworkflows/local/bam_baserecalibrator/main'


// Create recalibrated cram files to use for variant calling (+QC)
include { BAM_APPLYBQSR                               } from '../subworkflows/local/bam_applybqsr/main'


//
// MODULE: Installed directly from nf-core/modules
//

// Loading the FASTQC module
include { FASTQC                      } from '../modules/nf-core/fastqc/main'


// TRIM/SPLIT FASTQ Files
include { FASTP                                       } from '../modules/nf-core/fastp/main'


// Loading the MULTIQC module
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'


// Loading the module that dumps the versions of software being used
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'


// SKIPPING SOME FROM SAREK FOR NOW.....


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HEISENBIO() {

	ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
	ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
	ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
	ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

	// To gather all QC reports for MultiQC
	reports  = Channel.empty()
	// To gather used softwares versions for MultiQC
	versions = Channel.empty()

	if (params.download_cache) {
	    PREPARE_CACHE(ensemblvep_info, snpeff_info)
	    snpeff_cache = PREPARE_CACHE.out.snpeff_cache
	    vep_cache    = PREPARE_CACHE.out.ensemblvep_cache.map{ meta, cache -> [ cache ] }

	    versions = versions.mix(PREPARE_CACHE.out.versions)
    }

	// Build indices if needed
	PREPARE_GENOME(
		ascat_alleles,
		ascat_loci,
		ascat_loci_gc,
		ascat_loci_rt,
		chr_dir,
		dbsnp,
		fasta,
		fasta_fai,
		germline_resource,
		known_indels,
		known_snps,
		pon)

    // Gather built indices or get them from the params. If it is supplied it will take that or build it using PREPARE_GENOME() function. Built from the fasta file:
	dict       = params.dict        ? Channel.fromPath(params.dict).map{ it -> [ [id:'dict'], it ] }.collect()
				    : PREPARE_GENOME.out.dict
	fasta_fai  = params.fasta_fai   ? Channel.fromPath(params.fasta_fai).collect()
				    : PREPARE_GENOME.out.fasta_fai
	bwa        = params.bwa         ? Channel.fromPath(params.bwa).collect()
				    : PREPARE_GENOME.out.bwa
	bwamem2    = params.bwamem2     ? Channel.fromPath(params.bwamem2).collect()
				    : PREPARE_GENOME.out.bwamem2
    // Gather index for mapping given the chosen aligner
	index_alignement = (params.aligner == "bwa-mem" || params.aligner == "sentieon-bwamem") ? bwa :
	params.aligner == "bwa-mem2" ? bwamem2
}
    //Skipping some VCFs, ASCAT, and msisensorpro variables loading for now....
    //versions = versions.mix(PREPARE_GENOME.out.versions)
	// ALIGNMENT STARTING

	// MODULE: MultiQC
    //if (!(params.skip_tools && params.skip_tools.split(',').contains('multiqc'))) {
    //    workflow_summary    = WorkflowHeisenbio.paramsSummaryMultiqc(workflow, summary_params)
    //    ch_workflow_summary = Channel.value(workflow_summary)
    //    
    //    methods_description    = WorkflowHeisenbio.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    //    ch_methods_description = Channel.value(methods_description)
    //    
    //    ch_multiqc_files = Channel.empty()
    //    multiqc_files = multiqc_files.mix(version_yaml)
    //    multiqc_files = multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    //    multiqc_files = multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    //    multiqc_files = multiqc_files.mix(reports.collect().ifEmpty([]))
    //    
    //    MULTIQC(multiqc_files.collect(), ch_multiqc_config.collect().ifEmpty([]), ch_multiqc_custom_config.collect().ifEmpty([]), ch_multiqc_logo.collect().ifEmpty([]))
    //    multiqc_report = MULTIQC.out.report.toList()
    //    versions = versions.mix(MULTIQC.out.versions)
    //}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Parse first line of a FASTQ file, return the flowcell id and lane number.
def flowcellLaneFromFastq(path) {
    // expected format:
    // xx:yy:FLOWCELLID:LANE:... (seven fields)
    // or
    // FLOWCELLID:LANE:xx:... (five fields)
    def line
    path.withInputStream {
        InputStream gzipStream = new java.util.zip.GZIPInputStream(it)
        Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
        BufferedReader buffered = new BufferedReader(decoder)
        line = buffered.readLine()
    }
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(':')
    String fcid

    if (fields.size() >= 7) {
        // CASAVA 1.8+ format, from  https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
        // "@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>"
        fcid = fields[2]
    } else if (fields.size() == 5) {
        fcid = fields[0]
    }
    return fcid
}



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
