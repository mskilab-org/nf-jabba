/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'
def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
    params.cnvkit_reference,
    params.dbnsfp,
    params.dbnsfp_tbi,
    params.dbsnp,
    params.dbsnp_tbi,
    params.dict,
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
    params.gcmapdir_frag,
    params.pon_dryclean,
    params.blacklist_coverage_jabba
]

def toolParamMap = [
    "snpeff"     : [params.snpeff_cache],
    "vep"        : [params.vep_cache],
    "svaba"      : [params.indel_mask, params.germ_sv_db, params.simple_seq_db],
    "gridss"     : [params.blacklist_gridss, params.pon_gridss],
    "hetpileups" : [params.hapmap_sites],
    "fusions"    : [params.gencode_fusions],
    "allelic_cn" : [params.mask_non_integer_balance, params.mask_lp_phased_balance]
]

// Check if running tools and add their params to the checkPathParamList
if (params.tools) {
    toolParamMap.each { tool, toolParams ->
        if (params.tools.contains(tool)) {
            checkPathParamList.addAll(toolParams)
        }
    }
}

// Initialise the workflow
WorkflowNfjabba.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Check mandatory parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check if the path parameters exist
checkPathParamList.each { param ->
    if (param) file(param, checkIfExists: true)
}

// Set input, can either be from --input or from automatic retrieval in WorkflowSarek.groovy
def inputType = params.input ? "input" : "input_restart"
def ch_from_samplesheet = params.build_only_index ? Channel.empty() : Channel.fromSamplesheet(inputType)

// Function to handle errors based on step and data type
def handleError(step, dataType) {
    error("Your samplesheet contains ${dataType} files but step is `${step}`. Please check your samplesheet or adjust the step parameter.")
}

input_sample = ch_from_samplesheet
        .map{ meta, fastq_1, fastq_2, table, cram, crai, bam, bai, cov, hets, ploidy, vcf, vcf2, seg, nseg, ggraph, variantcaller ->
            // generate patient_sample key to group lanes together
            [ meta.patient + meta.sample, [meta, fastq_1, fastq_2, table, cram, crai, bam, bai, cov, hets, ploidy, vcf, vcf2, seg, nseg, ggraph, variantcaller] ]
        }
        .tap{ ch_with_patient_sample } // save the channel
        .groupTuple() //group by patient_sample to get all lanes
        .map { patient_sample, ch_items ->
            // get number of lanes per sample
            [ patient_sample, ch_items.size() ]
        }
        .combine(ch_with_patient_sample, by: 0) // for each entry add numLanes
        .map { patient_sample, num_lanes, ch_items ->

            (meta, fastq_1, fastq_2, table, cram, crai, bam, bai, cov, hets, ploidy, vcf, vcf2, seg, nseg, ggraph, variantcaller) = ch_items
            if (meta.lane && fastq_2) {
                meta           = meta + [id: "${meta.sample}-${meta.lane}".toString()]
                def CN         = params.seq_center ? "CN:${params.seq_center}\\t" : ''

                def flowcell   = flowcellLaneFromFastq(fastq_1)
                // Don't use a random element for ID, it breaks resuming
                def read_group = "\"@RG\\tID:${flowcell}.${meta.sample}.${meta.lane}\\t${CN}PU:${meta.lane}\\tSM:${meta.patient}_${meta.sample}\\tLB:${meta.sample}\\tDS:${params.fasta}\\tPL:${params.seq_platform}\""

                meta           = meta - meta.subMap('lane') + [num_lanes: num_lanes.toInteger(), read_group: read_group.toString(), data_type: 'fastq', size: 1]

                if (params.step == 'alignment') return [ meta, [ fastq_1, fastq_2 ] ]
                else {
                    handleError(params.step, 'fastq')
                }

            // start from BAM
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
                    handleError(params.step, 'bam')
                }

            // recalibration
            } else if (table && cram) {
                meta = meta + [id: meta.sample, data_type: 'cram']

                if (!(params.step == 'alignment' || params.step == 'annotate')) return [ meta - meta.subMap('lane'), cram, crai, table ]
                else {
                    handleError(params.step, 'cram')
                }

            // recalibration when skipping MarkDuplicates
            } else if (table && bam) {
                meta = meta + [id: meta.sample, data_type: 'bam']

                if (!(params.step == 'alignment' || params.step == 'annotate')) return [ meta - meta.subMap('lane'), bam, bai, table ]
                else {
                    handleError(params.step, 'bam')
                }

            // prepare_recalibration or variant_calling
            } else if (cram) {
                meta = meta + [id: meta.sample, data_type: 'cram']

                if (!(params.step == 'alignment' || params.step == 'annotate')) return [ meta - meta.subMap('lane'), cram, crai ]
                else {
                    handleError(params.step, 'cram')
                }

            // prepare_recalibration when skipping MarkDuplicates or `--step markduplicates`
            } else if (bam) {
                meta = meta + [id: meta.sample, data_type: 'bam']

                if (!(params.step == 'alignment' || params.step == 'annotate')) return [ meta - meta.subMap('lane'), bam, bai ]
                else {
                    handleError(params.step, 'bam')
                }
            // allelic_cn
            } else if (cov && hets && ggraph) {
                meta = meta + [id: meta.sample, data_type: ['cov', 'hets', 'ggraph']]

                if (params.step == 'allelic_cn') return [ meta - meta.subMap('lane'), cov, hets, ggraph]
                else {
                    handleError(params.step, 'coverage .rds, hets .txt and JaBbA ggraph .rds')
                }
                //jabba
            } else if (cov && hets && ploidy && vcf && vcf2 && seg && nseg) {
                meta = meta + [id: meta.sample, data_type: ['cov', 'vcf', 'hets', 'seg']]

                if (params.step == 'jabba') return [ meta - meta.subMap('lane'), cov, hets, vcf, vcf2, seg, nseg ]
                else {
                    handleError(params.step, 'coverage .rds, hets .txt, SV .vcf and segs .rds')
                }
                // ascat
            } else if (cov && hets) {
                meta = meta + [id: meta.sample, data_type: ['cov', 'hets']]

                if (params.step == 'ascat') return [ meta - meta.subMap('lane'), cov, hets ]
                else {
                    handleError(params.step, 'coverage .rds and hetpileups .txt')
                }
                // dryclean
            } else if (cov) {
                meta = meta + [id: meta.sample, data_type: 'cov']

                if (params.step == 'dryclean') return [ meta - meta.subMap('lane'), cov ]
                else {
                    handleError(params.step, 'coverage .rds')
                }
                // events and fusions
            } else if (ggraph) {
                meta = meta + [id: meta.sample, data_type: 'ggraph']

                if (params.step == 'events' || params.step == 'fusions') return [ meta - meta.subMap('lane'), ggraph ]
                else {
                    handleError(params.step, 'ggraph .rds')
                }
        } else {
                error("Missing or unknown field in csv file header. Please check your samplesheet")
            }
        }

if (!params.dbsnp && !params.known_indels) {
    if (params.step in ['alignment', 'markduplicates', 'prepare_recalibration', 'recalibrate'] && (!params.skip_tools || (params.skip_tools && !params.skip_tools.split(',').contains('baserecalibrator')))) {
        log.warn "Base quality score recalibration requires at least one resource file. Please provide at least one of `--dbsnp` or `--known_indels`\nYou can skip this step in the workflow by adding `--skip_tools baserecalibrator` to the command."
    }
    if (params.tools && (params.tools.split(',').contains('haplotypecaller') || params.tools.split(',').contains('sentieon_haplotyper'))) {
        log.warn "If GATK's Haplotypecaller or Sentieon's Haplotyper is specified, without `--dbsnp` or `--known_indels no filtering will be done. For filtering, please provide at least one of `--dbsnp` or `--known_indels`.\nFor more information see FilterVariantTranches (single-sample, default): https://gatk.broadinstitute.org/hc/en-us/articles/5358928898971-FilterVariantTranches\nFor more information see VariantRecalibration (--joint_germline): https://gatk.broadinstitute.org/hc/en-us/articles/5358906115227-VariantRecalibrator\nFor more information on GATK Best practice germline variant calling: https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-"
    }
}

// Fails when missing sex information for CNV tools
if (params.tools && (params.tools.split(',').contains('ascat'))) {
    input_sample.map{
        if (it[0].sex == 'NA' ) {
            log.warn "Please specify sex information for each sample in your samplesheet when using '--tools' with 'ascat' if known for comparison"
        }
    }
}

if ((params.download_cache) && (params.snpeff_cache || params.vep_cache)) {
    error("Please specify either `--download_cache` or `--snpeff_cache`, `--vep_cache`.")
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Initialize file channels based on params, defined in the params.genomes[params.genome] scope
ascat_alleles      = params.ascat_alleles      ? Channel.fromPath(params.ascat_alleles).collect()     : Channel.empty()
ascat_loci         = params.ascat_loci         ? Channel.fromPath(params.ascat_loci).collect()        : Channel.empty()
ascat_loci_gc      = params.ascat_loci_gc      ? Channel.fromPath(params.ascat_loci_gc).collect()     : Channel.value([])
ascat_loci_rt      = params.ascat_loci_rt      ? Channel.fromPath(params.ascat_loci_rt).collect()     : Channel.value([])
cf_chrom_len       = params.cf_chrom_len       ? Channel.fromPath(params.cf_chrom_len).collect()      : []
chr_dir            = params.chr_dir            ? Channel.fromPath(params.chr_dir).collect()           : Channel.value([])
dbsnp              = params.dbsnp              ? Channel.fromPath(params.dbsnp).collect()             : Channel.value([])
fasta              = params.fasta              ? Channel.fromPath(params.fasta).first()               : Channel.empty()
fasta_fai          = params.fasta_fai          ? Channel.fromPath(params.fasta_fai).collect()         : Channel.empty()
germline_resource  = params.germline_resource  ? Channel.fromPath(params.germline_resource).collect() : Channel.value([]) // Mutect2 does not require a germline resource, so set to optional input
known_indels       = params.known_indels       ? Channel.fromPath(params.known_indels).collect()      : Channel.value([])
known_snps         = params.known_snps         ? Channel.fromPath(params.known_snps).collect()        : Channel.value([])
mappability        = params.mappability        ? Channel.fromPath(params.mappability).collect()       : Channel.value([])
pon                = params.pon                ? Channel.fromPath(params.pon).collect()               : Channel.value([]) // PON is optional for Mutect2 (but highly recommended)

// SVABA
indel_mask         = params.indel_mask         ? Channel.fromPath(params.indel_mask).collect()        : Channel.empty()   // This is the indel mask for SVABA
germ_sv_db         = params.germ_sv_db         ? Channel.fromPath(params.germ_sv_db).collect()        : Channel.empty()   // This is the germline SV mask for Svaba
simple_seq_db      = params.simple_seq_db      ? Channel.fromPath(params.simple_seq_db).collect()     : Channel.empty()   // This is the file containing sites of simple DNA that can confuse the contig re-alignment for SVABA

// GRIDSS
blacklist_gridss   = params.blacklist_gridss   ? Channel.fromPath(params.blacklist_gridss).collect()  : Channel.empty()   // This is the mask for gridss SV calls
pon_gridss         = params.pon_gridss         ? Channel.fromPath(params.pon_gridss).collect()        : Channel.empty()   //This is the pon directory for GRIDSS SOMATIC. (MUST CONTAIN .bed and .bedpe files)

// FragCounter
gcmapdir_frag      = params.gcmapdir_frag      ? Channel.fromPath(params.gcmapdir_frag).collect()     : Channel.empty()   // This is the GC/Mappability directory for fragCounter. (Must contain gc* & map* .rds files)

// HetPileups
hapmap_sites       = params.hapmap_sites       ? Channel.fromPath(params.hapmap_sites).collect()      : Channel.empty()

// Dryclean
pon_dryclean      = params.pon_dryclean      ? Channel.fromPath(params.pon_dryclean).collect()     : Channel.empty()   // This is the path to the PON for Dryclean.
blacklist_path_dryclean      = params.blacklist_path_dryclean      ? Channel.fromPath(params.blacklist_path_dryclean).collect()     : Channel.empty()   // This is the path to the blacklist for Dryclean (optional).
germline_file_dryclean      = params.germline_file_dryclean      ? Channel.fromPath(params.germline_file_dryclean).collect()     : Channel.empty()   // This is the path to the germline mask for dryclean (optional).

// JaBbA
blacklist_coverage_jabba		= params.blacklist_coverage_jabba		  ? Channel.fromPath(params.blacklist_coverage_jabba).collect() : Channel.empty()

// Fusions
gencode_fusions             = params.gencode_fusions        ? Channel.fromPath(params.gencode_fusions).collect() : Channel.empty()

// Allelic CN
mask_non_integer_balance    = params.mask_non_integer_balance   ? Channel.fromPath(params.mask_non_integer_balance).collect() : Channel.empty()
mask_lp_phased_balance    = params.mask_lp_phased_balance   ? Channel.fromPath(params.mask_lp_phased_balance).collect() : Channel.empty()

// Initialize value channels based on params, defined in the params.genomes[params.genome] scope
ascat_genome       = params.ascat_genome       ?: Channel.empty()
dbsnp_vqsr         = params.dbsnp_vqsr         ? Channel.value(params.dbsnp_vqsr) : Channel.empty()
known_indels_vqsr  = params.known_indels_vqsr  ? Channel.value(params.known_indels_vqsr) : Channel.empty()
known_snps_vqsr    = params.known_snps_vqsr    ? Channel.value(params.known_snps_vqsr) : Channel.empty()
snpeff_db          = params.snpeff_db          ?: Channel.empty()
vep_cache_version  = params.vep_cache_version  ?: Channel.empty()
vep_genome         = params.vep_genome         ?: Channel.empty()
vep_species        = params.vep_species        ?: Channel.empty()
error_rate         = params.error_rate         ?: Channel.empty()                                                         // For SVABA

// Hetpileups
filter_hets         = params.filter_hets       ?: Channel.empty()
max_depth           = params.max_depth         ?: Channel.empty()

// FragCounter
windowsize_frag    = params.windowsize_frag    ?: Channel.empty()                                                         // For fragCounter
minmapq_frag       = params.minmapq_frag       ?: Channel.empty()                                                         // For fragCounter
midpoint_frag      = params.midpoint_frag      ?: Channel.empty()                                                         // For fragCounter
paired_frag        = params.paired_frag        ?: Channel.empty()                                                         // For fragCounter
exome_frag         = params.exome_frag         ?: Channel.empty()                                                         // For fragCounter

// Dryclean
centered_dryclean           = params.centered_dryclean          ?: Channel.empty()
cbs_dryclean                = params.cbs_dryclean               ?: Channel.empty()
cnsignif_dryclean           = params.cnsignif_dryclean          ?: Channel.empty()
wholeGenome_dryclean        = params.wholeGenome_dryclean       ?: Channel.empty()
blacklist_dryclean          = params.blacklist_dryclean         ?: Channel.empty()
germline_filter_dryclean    = params.germline_filter_dryclean   ?: Channel.empty()
human_dryclean              = params.human_dryclean             ?: Channel.empty()
field_dryclean              = params.field_dryclean             ?: Channel.empty()
build_dryclean              = params.build_dryclean             ?: Channel.empty()

// ASCAT_seg
field_ascat                 = params.field_ascat                ?: Channel.empty()
hets_thresh_ascat           = params.hets_thresh_ascat          ?: Channel.empty()
penalty_ascat               = params.penalty_ascat              ?: Channel.empty()
gc_correct_ascat            = params.gc_correct_ascat           ?: Channel.empty()
rebin_width_ascat           = params.rebin_width_ascat          ?: Channel.empty()
from_maf_ascat              = params.from_maf_ascat             ?: Channel.empty()

// CBS
cnsignif_cbs                    = params.cnsignif_cbs               ?: Channel.empty()
field_cbs                       = params.field_cbs                  ?: Channel.empty()
name_cbs                        = params.name_cbs                   ?: Channel.empty()

// JaBbA
blacklist_junctions_jabba       = params.blacklist_junctions_jabba      ?: Channel.empty()
geno_jabba					    = params.geno_jabba			            ?: Channel.empty()
indel_jabba					    = params.indel_jabba			        ?: Channel.empty()
tfield_jabba					= params.tfield_jabba			        ?: Channel.empty()
iter_jabba					    = params.iter_jabba			            ?: Channel.empty()
rescue_window_jabba				= params.rescue_window_jabba			?: Channel.empty()
rescue_all_jabba				= params.rescue_all_jabba			    ?: Channel.empty()
nudgebalanced_jabba				= params.nudgebalanced_jabba			?: Channel.empty()
edgenudge_jabba					= params.edgenudge_jabba			    ?: Channel.empty()
strict_jabba					= params.strict_jabba			        ?: Channel.empty()
allin_jabba					    = params.allin_jabba			        ?: Channel.empty()
field_jabba					    = params.field_jabba			        ?: Channel.empty()
maxna_jabba					    = params.maxna_jabba			        ?: Channel.empty()
purity_jabba					= params.purity_jabba                   ?: Channel.empty()
ploidy_jab     					= params.ploidy_jabba                   ?: Channel.empty()
pp_method_jabba					= params.pp_method_jabba                ?: Channel.empty()
cnsignif_jabba					= params.cnsignif_jabba                 ?: Channel.empty()
slack_jabba					    = params.slack_jabba                    ?: Channel.empty()
linear_jabba					= params.linear_jabba                   ?: Channel.empty()
tilim_jabba					    = params.tilim_jabba                    ?: Channel.empty()
epgap_jabba					    = params.epgap_jabba                    ?: Channel.empty()
fix_thres_jabba					= params.fix_thres_jabba			    ?: Channel.empty()
lp_jabba					    = params.lp_jabba			            ?: Channel.empty()
ism_jabba					    = params.ism_jabba			            ?: Channel.empty()
filter_loose_jabba				= params.filter_loose_jabba			    ?: Channel.empty()
gurobi_jabba					= params.gurobi_jabba			        ?: Channel.empty()
nonintegral_jabba				= params.nonintegral_jabba			    ?: Channel.empty()
verbose_jabba					= params.verbose_jabba			        ?: Channel.empty()
help_jabba					    = params.help_jabba			            ?: Channel.empty()

//Alleic CN (Non-integer balance)
field_non_integer_balance = params.field_non_integer_balance ?: Channel.empty()
hets_thresh_non_integer_balance = params.hets_thresh_non_integer_balance ?: Channel.empty()
overwrite_non_integer_balance = params.overwrite_non_integer_balance ?: Channel.empty()
lambda_non_integer_balance = params.lambda_non_integer_balance ?: Channel.empty()
allin_non_integer_balance = params.allin_non_integer_balance ?: Channel.empty()
fix_thresh_non_integer_balance = params.fix_thresh_non_integer_balance ?: Channel.empty()
nodebounds_non_integer_balance = params.nodebounds_non_integer_balance ?: Channel.empty()
ism_non_integer_balance = params.ism_non_integer_balance ?: Channel.empty()
build_non_integer_balance = params.build_non_integer_balance ?: Channel.empty()
epgap_non_integer_balance = params.epgap_non_integer_balance ?: Channel.empty()
tilim_non_integer_balance = params.tilim_non_integer_balance ?: Channel.empty()
gurobi_non_integer_balance = params.gurobi_non_integer_balance ?: Channel.empty()
pad_non_integer_balance = params.pad_non_integer_balance ?: Channel.empty()

// ...(LP Phased balance)
lambda_lp_phased_balance = params.lambda_lp_phased_balance ?: Channel.empty()
cnloh_lp_phased_balance = params.cnloh_lp_phased_balance ?: Channel.empty()
major_lp_phased_balance = params.major_lp_phased_balance ?: Channel.empty()
allin_lp_phased_balance = params.allin_lp_phased_balance ?: Channel.empty()
marginal_lp_phased_balance = params.marginal_lp_phased_balance ?: Channel.empty()
from_maf_lp_phased_balance = params.from_maf_lp_phased_balance ?: Channel.empty()
ism_lp_phased_balance = params.ism_lp_phased_balance ?: Channel.empty()
epgap_lp_phased_balance = params.epgap_lp_phased_balance ?: Channel.empty()
hets_thresh_lp_phased_balance = params.hets_thresh_lp_phased_balance ?: Channel.empty()
min_bins_lp_phased_balance = params.min_bins_lp_phased_balance ?: Channel.empty()
min_width_lp_phased_balance = params.min_width_lp_phased_balance ?: Channel.empty()
trelim_lp_phased_balance = params.trelim_lp_phased_balance ?: Channel.empty()
reward_lp_phased_balance = params.reward_lp_phased_balance ?: Channel.empty()
nodefileind_lp_phased_balance = params.nodefileind_lp_phased_balance ?: Channel.empty()
tilim_lp_phased_balance = params.tilim_lp_phased_balance ?: Channel.empty()

// Initialize files channels based on params, not defined within the params.genomes[params.genome] scope
if (params.snpeff_cache && params.tools && params.tools.contains("snpeff")) {
    def snpeff_annotation_cache_key = params.use_annotation_cache_keys ? "${params.snpeff_genome}.${params.snpeff_db}/" : ""
    def snpeff_cache_dir =  "${snpeff_annotation_cache_key}${params.snpeff_genome}.${params.snpeff_db}"
    def snpeff_cache_path_full = file("$params.snpeff_cache/$snpeff_cache_dir", type: 'dir')
    if ( !snpeff_cache_path_full.exists() || !snpeff_cache_path_full.isDirectory() ) {
        error("Files within --snpeff_cache invalid. Make sure there is a directory named ${snpeff_cache_dir} in ${params.snpeff_cache}.\nhttps://nf-co.re/sarek/dev/usage#how-to-customise-snpeff-and-vep-annotation")
    }
    snpeff_cache = Channel.fromPath(file("${params.snpeff_cache}/${snpeff_annotation_cache_key}"), checkIfExists: true).collect()
        .map{ cache -> [ [ id:"${params.snpeff_genome}.${params.snpeff_db}" ], cache ] }
} else snpeff_cache = []

if (params.vep_cache && params.tools && params.tools.contains("vep")) {
    def vep_annotation_cache_key = params.use_annotation_cache_keys ? "${params.vep_cache_version}_${params.vep_genome}/" : ""
    def vep_cache_dir = "${vep_annotation_cache_key}${params.vep_species}/${params.vep_cache_version}_${params.vep_genome}"
    def vep_cache_path_full = file("$params.vep_cache/$vep_cache_dir", type: 'dir')
    if ( !vep_cache_path_full.exists() || !vep_cache_path_full.isDirectory() ) {
        error("Files within --vep_cache invalid. Make sure there is a directory named ${vep_cache_dir} in ${params.vep_cache}.\nhttps://nf-co.re/sarek/dev/usage#how-to-customise-snpeff-and-vep-annotation")
    }
    vep_cache = Channel.fromPath(file("${params.vep_cache}/${vep_annotation_cache_key}"), checkIfExists: true).collect()
} else vep_cache = []

vep_extra_files = []


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
include { CHANNEL_SVCALLING_CREATE_CSV                } from '../subworkflows/local/channel_svcalling_create_csv/main'
include { CHANNEL_HETPILEUPS_CREATE_CSV               } from '../subworkflows/local/channel_hetpileups_create_csv/main'

// Download annotation cache if needed
include { PREPARE_CACHE                               } from '../subworkflows/local/prepare_cache/main'

// Build indices if needed
include { PREPARE_GENOME                              } from '../subworkflows/local/prepare_genome/main'

// Build intervals if needed
include { PREPARE_INTERVALS                           } from '../subworkflows/local/prepare_intervals/main'

// Convert BAM files to FASTQ files
include { BAM_CONVERT_SAMTOOLS as CONVERT_FASTQ_INPUT } from '../subworkflows/local/bam_convert_samtools/main'

// Run FASTQC
include { FASTQC                                      } from '../modules/nf-core/fastqc/main'

// TRIM/SPLIT FASTQ Files
include { FASTP                                       } from '../modules/nf-core/fastp/main'

// Loading the MULTIQC module
include { MULTIQC                                     } from '../modules/nf-core/multiqc/main'

// Loading the module that dumps the versions of software being used
include { CUSTOM_DUMPSOFTWAREVERSIONS                 } from '../modules/nf-core/custom/dumpsoftwareversions/main'

// Map input reads to reference genome
include { FASTQ_ALIGN_BWAMEM_MEM2                     } from '../subworkflows/local/fastq_align_bwamem_mem2/main'

// Merge and index BAM files (optional)
include { BAM_MERGE_INDEX_SAMTOOLS                    } from '../subworkflows/local/bam_merge_index_samtools/main'

// Convert BAM files
include { SAMTOOLS_CONVERT as BAM_TO_CRAM             } from '../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_CONVERT as BAM_TO_CRAM_MAPPING     } from '../modules/nf-core/samtools/convert/main'

// Convert CRAM files (optional)
include { SAMTOOLS_CONVERT as CRAM_TO_BAM             } from '../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_CONVERT as CRAM_TO_BAM_RECAL       } from '../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_CONVERT as CRAM_TO_BAM_FINAL       } from '../modules/nf-core/samtools/convert/main'

// Mark Duplicates (+QC)
include { BAM_MARKDUPLICATES                          } from '../subworkflows/local/bam_markduplicates/main'

// QC on CRAM
include { CRAM_QC_MOSDEPTH_SAMTOOLS as CRAM_QC_NO_MD  } from '../subworkflows/local/cram_qc_mosdepth_samtools/main'
include { CRAM_QC_MOSDEPTH_SAMTOOLS as CRAM_QC_RECAL  } from '../subworkflows/local/cram_qc_mosdepth_samtools/main'

// Create recalibration tables
include { BAM_BASERECALIBRATOR                        } from '../subworkflows/local/bam_baserecalibrator/main'

// Create recalibrated cram files to use for variant calling (+QC)
include { BAM_APPLYBQSR                               } from '../subworkflows/local/bam_applybqsr/main'

// Svaba
include { BAM_SVCALLING_SVABA                         } from '../subworkflows/local/bam_svcalling_svaba/main'

//GRIDSS
include { BAM_SVCALLING_GRIDSS                        } from '../subworkflows/local/bam_svcalling_gridss/main'
include { BAM_SVCALLING_GRIDSS_SOMATIC                } from '../subworkflows/local/bam_svcalling_gridss/main'

// HETPILEUPS
include { BAM_HETPILEUPS                              } from '../subworkflows/local/bam_hetpileups/main'

// fragCounter
include { BAM_FRAGCOUNTER as TUMOR_FRAGCOUNTER         } from '../subworkflows/local/bam_fragCounter/main'
include { BAM_FRAGCOUNTER as NORMAL_FRAGCOUNTER        } from '../subworkflows/local/bam_fragCounter/main'

// dryclean
include { COV_DRYCLEAN as TUMOR_DRYCLEAN               } from '../subworkflows/local/cov_dryclean/main'
include { COV_DRYCLEAN as NORMAL_DRYCLEAN              } from '../subworkflows/local/cov_dryclean/main'

//ASCAT
include { COV_ASCAT                                   } from '../subworkflows/local/cov_ascat/main'
include { EXTRACT_PURITYPLOIDY                        } from '../modules/local/ascat/main'

// CBS
include { COV_CBS as CBS                              } from '../subworkflows/local/cov_cbs/main'

// JaBbA
include { COV_JUNC_JABBA as JABBA                     } from '../subworkflows/local/jabba/main'
include { COV_JUNC_JABBA as JABBA_WITH_SVABA          } from '../subworkflows/local/jabba/main'
include { COV_JUNC_JABBA as JABBA_WITH_GRIDSS         } from '../subworkflows/local/jabba/main'

// Events
include { GGRAPH_EVENTS as EVENTS                            } from '../subworkflows/local/events/main'
include { GGRAPH_EVENTS as EVENTS_WITH_GRIDSS                } from '../subworkflows/local/events/main'
include { GGRAPH_EVENTS as EVENTS_WITH_SVABA                 } from '../subworkflows/local/events/main'

// Fusions
include { GGRAPH_FUSIONS as FUSIONS                            } from '../subworkflows/local/fusions/main'
include { GGRAPH_FUSIONS as FUSIONS_WITH_GRIDSS                } from '../subworkflows/local/fusions/main'
include { GGRAPH_FUSIONS as FUSIONS_WITH_SVABA                 } from '../subworkflows/local/fusions/main'

// Alleic CN
include { COV_GGRAPH_NON_INTEGER_BALANCE as NON_INTEGER_BALANCE                            } from '../subworkflows/local/allelic_cn/main'
include { COV_GGRAPH_NON_INTEGER_BALANCE as NON_INTEGER_BALANCE_WITH_GRIDSS                } from '../subworkflows/local/allelic_cn/main'
include { COV_GGRAPH_NON_INTEGER_BALANCE as NON_INTEGER_BALANCE_WITH_SVABA                 } from '../subworkflows/local/allelic_cn/main'

include { COV_GGRAPH_LP_PHASED_BALANCE as LP_PHASED_BALANCE                            } from '../subworkflows/local/allelic_cn/main'
include { COV_GGRAPH_LP_PHASED_BALANCE as LP_PHASED_BALANCE_WITH_GRIDSS                } from '../subworkflows/local/allelic_cn/main'
include { COV_GGRAPH_LP_PHASED_BALANCE as LP_PHASED_BALANCE_WITH_SVABA                 } from '../subworkflows/local/allelic_cn/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NFJABBA {

    // MULTIQC
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    // To gather all QC reports for MultiQC
    reports  = Channel.empty()
    // To gather used softwares versions for MultiQC
    versions = Channel.empty()

        // Download cache if needed
    // Assuming that if the cache is provided, the user has already downloaded it
    ensemblvep_info = params.vep_cache    ? [] : Channel.of([ [ id:"${params.vep_cache_version}_${params.vep_genome}" ], params.vep_genome, params.vep_species, params.vep_cache_version ])
    snpeff_info     = params.snpeff_cache ? [] : Channel.of([ [ id:"${params.snpeff_genome}.${params.snpeff_db}" ], params.snpeff_genome, params.snpeff_db ])

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

    // Gather built indices or get them from the params
    // Built from the fasta file:
    dict       = params.dict        ? Channel.fromPath(params.dict).map{ it -> [ [id:'dict'], it ] }.collect()
                                    : PREPARE_GENOME.out.dict
    fasta_fai  = params.fasta_fai   ? Channel.fromPath(params.fasta_fai).collect()
                                    : PREPARE_GENOME.out.fasta_fai
    bwa        = params.bwa         ? Channel.fromPath(params.bwa).collect()
                                    : PREPARE_GENOME.out.bwa
    bwamem2    = params.bwamem2     ? Channel.fromPath(params.bwamem2).collect()
                                    : PREPARE_GENOME.out.bwamem2

    // Gather index for mapping given the chosen aligner
    index_alignement = (params.aligner == "bwa-mem") ? bwa :
        params.aligner == "bwa-mem2" ? bwamem2 : null

    // TODO: add a params for msisensorpro_scan
    msisensorpro_scan      = PREPARE_GENOME.out.msisensorpro_scan

    // For ASCAT, extracted from zip or tar.gz files:
    allele_files           = PREPARE_GENOME.out.allele_files
    chr_files              = PREPARE_GENOME.out.chr_files
    gc_file                = PREPARE_GENOME.out.gc_file
    loci_files             = PREPARE_GENOME.out.loci_files
    rt_file                = PREPARE_GENOME.out.rt_file

    // Tabix indexed vcf files:
    dbsnp_tbi              = params.dbsnp                   ? params.dbsnp_tbi             ? Channel.fromPath(params.dbsnp_tbi).collect()             : PREPARE_GENOME.out.dbsnp_tbi             : Channel.value([])
    germline_resource_tbi  = params.germline_resource       ? params.germline_resource_tbi ? Channel.fromPath(params.germline_resource_tbi).collect() : PREPARE_GENOME.out.germline_resource_tbi : [] //do not change to Channel.value([]), the check for its existence then fails for Getpileupsumamries
    known_indels_tbi       = params.known_indels            ? params.known_indels_tbi      ? Channel.fromPath(params.known_indels_tbi).collect()      : PREPARE_GENOME.out.known_indels_tbi      : Channel.value([])
    known_snps_tbi         = params.known_snps              ? params.known_snps_tbi        ? Channel.fromPath(params.known_snps_tbi).collect()        : PREPARE_GENOME.out.known_snps_tbi        : Channel.value([])
    pon_tbi                = params.pon                     ? params.pon_tbi               ? Channel.fromPath(params.pon_tbi).collect()               : PREPARE_GENOME.out.pon_tbi               : Channel.value([])

    // known_sites is made by grouping both the dbsnp and the known snps/indels resources
    // Which can either or both be optional
    known_sites_indels     = dbsnp.concat(known_indels).collect()
    known_sites_indels_tbi = dbsnp_tbi.concat(known_indels_tbi).collect()

    known_sites_snps       = dbsnp.concat(known_snps).collect()
    known_sites_snps_tbi   = dbsnp_tbi.concat(known_snps_tbi).collect()

    // Build intervals if needed
    PREPARE_INTERVALS(fasta_fai, params.intervals, params.no_intervals)

    // Intervals for speed up preprocessing/variant calling by spread/gather
    // [interval.bed] all intervals in one file
    intervals_bed_combined         = params.no_intervals ? Channel.value([])      : PREPARE_INTERVALS.out.intervals_bed_combined
    intervals_bed_gz_tbi_combined  = params.no_intervals ? Channel.value([])      : PREPARE_INTERVALS.out.intervals_bed_gz_tbi_combined

    // For QC during preprocessing, we don't need any intervals (MOSDEPTH doesn't take them for WGS)
    intervals_for_preprocessing = params.wes ?
        intervals_bed_combined.map{it -> [ [ id:it.baseName ], it ]}.collect() :
        Channel.value([ [ id:'null' ], [] ])

    intervals            = PREPARE_INTERVALS.out.intervals_bed        // [ interval, num_intervals ] multiple interval.bed files, divided by useful intervals for scatter/gather
    intervals_bed_gz_tbi = PREPARE_INTERVALS.out.intervals_bed_gz_tbi // [ interval_bed, tbi, num_intervals ] multiple interval.bed.gz/.tbi files, divided by useful intervals for scatter/gather

    intervals_and_num_intervals = intervals.map{ interval, num_intervals ->
        if ( num_intervals < 1 ) [ [], num_intervals ]
        else [ interval, num_intervals ]
    }

    intervals_bed_gz_tbi_and_num_intervals = intervals_bed_gz_tbi.map{ intervals, num_intervals ->
        if ( num_intervals < 1 ) [ [], [], num_intervals ]
        else [ intervals[0], intervals[1], num_intervals ]
    }

    // Gather used softwares versions
    versions = versions.mix(PREPARE_GENOME.out.versions)
    versions = versions.mix(PREPARE_INTERVALS.out.versions)

    // Initialize all flags to false
    boolean runAlignment = false
    boolean runMarkDuplicates = false
    boolean runPrepareRecalibration = false
    boolean runRecalibrate = false
    boolean runSvCalling = false
    boolean runFragCounter = false
    boolean runHetPileups = false
    boolean runDryClean = false
    boolean runAscat = false
    boolean runJabba = false
    boolean runEvents = false
    boolean runFusions = false
    boolean runAlleicCN = false

    // Set flags based on params.step using a cascading approach
    // Fall through to the next case if the previous case is true
    switch (params.step) {
        case 'alignment':
            runAlignment = true
        case 'markduplicates':
            runMarkDuplicates = true
        case 'prepare_recalibration':
            runPrepareRecalibration = true
        case 'recalibrate':
            runRecalibrate = true
        case 'sv_calling':
            runSvCalling = true
        case 'fragcounter':
            runFragCounter = true
        case 'hetpileups':
            runHetPileups = true
        case 'dryclean':
            runDryClean = true
        case 'ascat':
            runAscat = true
        case 'jabba':
            runJabba = true
        case 'events':
            runEvents = true
        case 'fusions':
            runFusions = true
        case 'allelic_cn':
            runAlleicCN = true
            break
        default:
            error "Invalid step: ${params.step}"
    }

    tools_used = params.tools ? params.tools.split(',') : error("No tools provided")

    if (runAlignment) {

        // Figure out if input is bam or fastq
        input_sample_type = input_sample.branch{
            bam:   it[0].data_type == "bam"
            fastq: it[0].data_type == "fastq"
        }

        // Convert any bam input to fastq
        // fasta are not needed when converting bam to fastq -> [ id:"fasta" ], []
        // No need for fasta.fai -> []
        interleave_input = false // Currently don't allow interleaved input
        CONVERT_FASTQ_INPUT(
            input_sample_type.bam,
            [ [ id:"fasta" ], [] ], // fasta
            [ [ id:'null' ], [] ],  // fasta_fai
            interleave_input)

        // Gather fastq (inputed or converted)
        // Theorically this could work on mixed input (fastq for one sample and bam for another)
        // But not sure how to handle that with the samplesheet
        // Or if we really want users to be able to do that
        input_fastq = input_sample_type.fastq.mix(CONVERT_FASTQ_INPUT.out.reads)
        //input_fastq.view()
        // STEP 0: QC & TRIMMING
        // `--skip_tools fastqc` to skip fastqc
        // Trim only with `--trim_fastq`
        // Additional options to be set up

        // QC
        if (!(params.skip_tools && params.skip_tools.split(',').contains('fastqc'))) {
            FASTQC(input_fastq)

            reports = reports.mix(FASTQC.out.zip.collect{ meta, logs -> logs })
            versions = versions.mix(FASTQC.out.versions.first())
        }

        //skipping the UMI Conscensus calling step for now
        reads_for_fastp = input_fastq

        // Trimming and/or splitting
        if (params.trim_fastq && params.split_fastq > 0) {
            log.warn "You have mentioned trim_fastq to `$params.trim_fastq`, will do trimming"
            save_trimmed_fail = false
            save_merged = false
            FASTP(
                reads_for_fastp,
                [], // we are not using any adapter fastas at the moment
                save_trimmed_fail,
                save_merged
            )

            reports = reports.mix(FASTP.out.json.collect{ meta, json -> json })
            reports = reports.mix(FASTP.out.html.collect{ meta, html -> html })

            if (params.split_fastq) {
                reads_for_alignment = FASTP.out.reads.map{ meta, reads ->
                    read_files = reads.sort(false) { a,b -> a.getName().tokenize('.')[0] <=> b.getName().tokenize('.')[0] }.collate(2)
                    [ meta + [ size:read_files.size() ], read_files ]
                }.transpose()
            } else reads_for_alignment = FASTP.out.reads

            versions = versions.mix(FASTP.out.versions)

        } else {
            println "Skipping trimming since trim_fastq is false"
            reads_for_alignment = reads_for_fastp
        }

        // STEP 1: MAPPING READS TO REFERENCE GENOME
        // reads will be sorted
        reads_for_alignment = reads_for_alignment.map{ meta, reads ->
            // Update meta.id to meta.sample no multiple lanes or splitted fastqs
            if (meta.size * meta.num_lanes == 1) [ meta + [ id:meta.sample ], reads ]
            else [ meta, reads ]
        }

        //reads_for_alignment.view()

        println "Starting Alignment Process..."
        sort_bam = true
        //fasta.view()
        //fasta_fai.view()
        //index_alignement.view()
        FASTQ_ALIGN_BWAMEM_MEM2(reads_for_alignment, index_alignement, sort_bam, fasta, fasta_fai)
        // Grouping the bams from the same samples not to stall the workflow
        bam_mapped = FASTQ_ALIGN_BWAMEM_MEM2.out.bam.map{ meta, bam ->

            // Update meta.id to be meta.sample, ditching sample-lane that is not needed anymore
            // Update meta.data_type
            // Remove no longer necessary fields:
            //   read_group: Now in the BAM header
            //    num_lanes: only needed for mapping
            //         size: only needed for mapping

            // Use groupKey to make sure that the correct group can advance as soon as it is complete
            // and not stall the workflow until all reads from all channels are mapped
            [ groupKey( meta - meta.subMap('num_lanes', 'read_group', 'size') + [ data_type:'bam', id:meta.sample ], (meta.num_lanes ?: 1) * (meta.size ?: 1)), bam ]
        }.groupTuple()

        if (
            params.save_mapped ||
            (
                (params.skip_tools && params.skip_tools.split(',').contains('markduplicates')) &&
                !(tools_used.contains('sentieon_dedup'))
            )
        ) {
            // bams are merged (when multiple lanes from the same sample), indexed and then converted to cram
            BAM_MERGE_INDEX_SAMTOOLS(bam_mapped)

            BAM_TO_CRAM_MAPPING(BAM_MERGE_INDEX_SAMTOOLS.out.bam_bai, fasta, fasta_fai)
            // Create CSV to restart from this step
            params.save_output_as_bam ? CHANNEL_ALIGN_CREATE_CSV(BAM_MERGE_INDEX_SAMTOOLS.out.bam_bai) : CHANNEL_ALIGN_CREATE_CSV(BAM_TO_CRAM_MAPPING.out.alignment_index)

            // Gather used softwares versions
            versions = versions.mix(BAM_MERGE_INDEX_SAMTOOLS.out.versions)
            versions = versions.mix(BAM_TO_CRAM_MAPPING.out.versions)
        }

        // Gather used softwares versions
        versions = versions.mix(CONVERT_FASTQ_INPUT.out.versions)
        versions = versions.mix(FASTQ_ALIGN_BWAMEM_MEM2.out.versions)
    }

    if (runMarkDuplicates) {

        // ch_cram_no_markduplicates_restart = Channel.empty()
        cram_markduplicates_no_spark = Channel.empty()

         // STEP 2: markduplicates (+QC) + convert to CRAM
        // ch_bam_for_markduplicates will contain bam mapped with FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON when step is mapping
        // Or bams that are specified in the samplesheet.csv when step is prepare_recalibration
        cram_for_markduplicates = params.step == 'alignment' ? bam_mapped : input_sample.map{ meta, input, index -> [ meta, input ] }
        // if no MD is done, then run QC on mapped & converted CRAM files
        // or the input BAM (+converted) or CRAM files
        cram_skip_markduplicates = Channel.empty()

        // Should it be possible to restart from converted crams?
        // For now, conversion from bam to cram is only done when skipping markduplicates

        if (
            params.skip_tools &&
            params.skip_tools.split(',').contains('markduplicates')
        ) {
            if (params.step == 'alignment') {
                cram_skip_markduplicates = BAM_TO_CRAM_MAPPING.out.alignment_index
            } else {
                input_markduplicates_convert = input_sample.branch{
                    bam:  it[0].data_type == "bam"
                    cram: it[0].data_type == "cram"
                }
                // Convert any input BAMs to CRAM
                BAM_TO_CRAM(input_markduplicates_convert.bam, fasta, fasta_fai)
                versions = versions.mix(BAM_TO_CRAM.out.versions)

                cram_skip_markduplicates = Channel.empty().mix(input_markduplicates_convert.cram, BAM_TO_CRAM.out.alignment_index)
            }
            CRAM_QC_NO_MD(cram_skip_markduplicates, fasta, intervals_for_preprocessing)

            // Gather QC reports
            reports = reports.mix(CRAM_QC_NO_MD.out.reports.collect{ meta, report -> report })

            // Gather used softwares versions
            versions = versions.mix(CRAM_QC_NO_MD.out.versions)
        } else {
            BAM_MARKDUPLICATES(
                cram_for_markduplicates,
                fasta,
                fasta_fai,
                intervals_for_preprocessing)

            cram_markduplicates_no_spark = BAM_MARKDUPLICATES.out.cram

            // Gather QC reports
            reports = reports.mix(BAM_MARKDUPLICATES.out.reports.collect{ meta, report -> report })

            // Gather used softwares versions
            versions = versions.mix(BAM_MARKDUPLICATES.out.versions)
        }
        // ch_md_cram_for_restart contains either:
        // - crams from markduplicates
        // - crams from sentieon_dedup
        // - crams from markduplicates_spark
        // - crams from input step markduplicates --> from the converted ones only?
        ch_md_cram_for_restart = Channel.empty().mix(cram_markduplicates_no_spark)
            // Make sure correct data types are carried through
            .map{ meta, cram, crai -> [ meta + [data_type: "cram"], cram, crai ] }
        //ch_md_cram_for_restart.view()
        // If params.save_output_as_bam, then convert CRAM files to BAM
        CRAM_TO_BAM(ch_md_cram_for_restart, fasta, fasta_fai)
        versions = versions.mix(CRAM_TO_BAM.out.versions)

        // CSV should be written for the file actually out, either CRAM or BAM
        // Create CSV to restart from this step
        csv_subfolder = (tools_used.contains('sentieon_dedup')) ? 'sentieon_dedup' : 'markduplicates'

        params.save_output_as_bam ? CHANNEL_MARKDUPLICATES_CREATE_CSV(CRAM_TO_BAM.out.alignment_index, csv_subfolder, params.outdir, params.save_output_as_bam) : CHANNEL_MARKDUPLICATES_CREATE_CSV(ch_md_cram_for_restart, csv_subfolder, params.outdir, params.save_output_as_bam)
    }

    if (runPrepareRecalibration) {
        // Run if starting from step "prepare_recalibration"
        if (params.step == 'prepare_recalibration') {
                        // Support if starting from BAM or CRAM files
            input_prepare_recal_convert = input_sample.branch{
                bam: it[0].data_type == "bam"
                cram: it[0].data_type == "cram"
            }
            // BAM files first must be converted to CRAM files since from this step on we base everything on CRAM format
            BAM_TO_CRAM(input_prepare_recal_convert.bam, fasta, fasta_fai)
            versions = versions.mix(BAM_TO_CRAM.out.versions)

            ch_cram_from_bam = BAM_TO_CRAM.out.alignment_index
                // Make sure correct data types are carried through
                .map{ meta, cram, crai -> [ meta + [data_type: "cram"], cram, crai ] }

            ch_cram_for_bam_baserecalibrator = Channel.empty().mix(ch_cram_from_bam, input_prepare_recal_convert.cram)
            ch_md_cram_for_restart = ch_cram_from_bam

        } else {
        // ch_cram_for_bam_baserecalibrator contains either:
            // - crams from markduplicates
            // - crams from markduplicates_spark
            // - crams converted from bam mapped when skipping markduplicates
            // - input cram files, when start from step markduplicates
            ch_cram_for_bam_baserecalibrator = Channel.empty().mix(ch_md_cram_for_restart, cram_skip_markduplicates )
                // Make sure correct data types are carried through
                .map{ meta, cram, crai -> [ meta + [data_type: "cram"], cram, crai ] }

        }
        if (!(params.skip_tools && params.skip_tools.split(',').contains('baserecalibrator'))) {

            ch_table_bqsr_tab    = Channel.empty()

            BAM_BASERECALIBRATOR(
                ch_cram_for_bam_baserecalibrator,
                dict,
                fasta,
                fasta_fai,
                intervals_and_num_intervals,
                known_sites_indels,
                known_sites_indels_tbi)

                ch_table_bqsr_tab = BAM_BASERECALIBRATOR.out.table_bqsr

                versions = versions.mix(BAM_BASERECALIBRATOR.out.versions)
        }
        // ch_table_bqsr contains table from baserecalibrator
        ch_table_bqsr = Channel.empty().mix(ch_table_bqsr_tab)

        reports = reports.mix(ch_table_bqsr.collect{ meta, table -> table })
        cram_applybqsr = ch_cram_for_bam_baserecalibrator.join(ch_table_bqsr, failOnDuplicate: true, failOnMismatch: true)

        // Create CSV to restart from this step
        CHANNEL_BASERECALIBRATOR_CREATE_CSV(ch_md_cram_for_restart.join(ch_table_bqsr, failOnDuplicate: true), params.tools, params.skip_tools, params.save_output_as_bam, params.outdir)
    }

    // STEP 4: RECALIBRATING
    if (runRecalibrate) {
        // Run if starting from step "prepare_recalibration"
        if (params.step == 'recalibrate') {
            // Support if starting from BAM or CRAM files
            input_recal_convert = input_sample.branch{
                bam:  it[0].data_type == "bam"
                cram: it[0].data_type == "cram"
            }

            // If BAM file, split up table and mapped file to convert BAM to CRAM
            input_only_table = input_recal_convert.bam.map{ meta, bam, bai, table -> [ meta, table ] }
            input_only_bam   = input_recal_convert.bam.map{ meta, bam, bai, table -> [ meta, bam, bai ] }

            // BAM files first must be converted to CRAM files since from this step on we base everything on CRAM format
            BAM_TO_CRAM(input_only_bam, fasta, fasta_fai)
            versions = versions.mix(BAM_TO_CRAM.out.versions)

            cram_applybqsr = Channel.empty().mix(
                BAM_TO_CRAM.out.alignment_index.join(input_only_table, failOnDuplicate: true, failOnMismatch: true),
                input_recal_convert.cram)
                // Join together converted cram with input tables
                .map{ meta, cram, crai, table -> [ meta + [data_type: "cram"], cram, crai, table ]}

        }
        if (!(params.skip_tools && params.skip_tools.split(',').contains('baserecalibrator'))) {
            cram_variant_calling_tab    = Channel.empty()

            BAM_APPLYBQSR(
                cram_applybqsr,
                dict,
                fasta,
                fasta_fai,
                intervals_and_num_intervals
            )

            cram_variant_calling_tab = BAM_APPLYBQSR.out.cram

            // Gather used softwares versions
            versions = versions.mix(BAM_APPLYBQSR.out.versions)

            cram_variant_calling = Channel.empty().mix(cram_variant_calling_tab)

            CRAM_QC_RECAL(
                cram_variant_calling,
                fasta,
                intervals_for_preprocessing
                )

            // Gather QC
            reports = reports.mix(CRAM_QC_RECAL.out.reports.collect{ meta, report -> report })

            // Gather software versions
            versions = versions.mix(CRAM_QC_RECAL.out.versions)

            // If params.save_output_as_bam, then convert CRAM files to BAM
            CRAM_TO_BAM_RECAL(cram_variant_calling, fasta, fasta_fai)
            versions = versions.mix(CRAM_TO_BAM_RECAL.out.versions)

            // CSV should be written for the file actually out out, either CRAM or BAM
            csv_recalibration = Channel.empty()
            csv_recalibration = params.save_output_as_bam ?  CRAM_TO_BAM_RECAL.out.alignment_index : cram_variant_calling

            // Create CSV to restart from this step
            CHANNEL_APPLYBQSR_CREATE_CSV(csv_recalibration)

        } else if (params.step == 'recalibrate') {
            // cram_variant_calling contains either:
            // - input bams converted to crams, if started from step recal + skip BQSR
            // - input crams if started from step recal + skip BQSR
            cram_variant_calling = Channel.empty().mix(
                  BAM_TO_CRAM.out.alignment_index,
                  input_recal_convert.cram.map{ meta, cram, crai, table -> [ meta, cram, crai ] })
        } else {
            // cram_variant_calling contains either:
            // - crams from markduplicates = ch_cram_for_bam_baserecalibrator if skip BQSR but not started from step recalibration
            cram_variant_calling = Channel.empty().mix(ch_cram_for_bam_baserecalibrator)
        }
        cram_sv_calling        = cram_variant_calling

        // Converting to BAM files to work downstream (SvABA has very low success rate with CRAMs)
        CRAM_TO_BAM_FINAL(cram_sv_calling, fasta, fasta_fai)
        versions = versions.mix(CRAM_TO_BAM_FINAL.out.versions)

        // Gets the BAM files in a channel (format: [meta, bam, bai]); confirms data type is correct
        bam_sv_calling = Channel.empty().mix(CRAM_TO_BAM_FINAL.out.alignment_index)
                            .map{ meta, bam, bai -> [ meta + [data_type: "bam"], bam, bai ] }

        //cram_fragcounter_calling  = cram_variant_calling
    }



    if (runSvCalling) {
        //when starting from sv_calling
        if (params.step == 'sv_calling') {
            input_sv_calling_convert = input_sample.branch{
                bam:  it[0].data_type == "bam"
                cram: it[0].data_type == "cram"
            }
            // BAM files first must be converted to CRAM files since from this step on we base everything on CRAM format
            CRAM_TO_BAM(input_sv_calling_convert.cram, fasta, fasta_fai)
            versions = versions.mix(CRAM_TO_BAM.out.versions)

            bam_sv_calling = Channel.empty().mix(CRAM_TO_BAM.out.alignment_index, input_sv_calling_convert.bam)
                                .map{ meta, bam, bai -> [ meta + [data_type: "bam"], bam, bai ] }           //making sure that the input data_type is correct

        }

        // getting the tumor and normal cram files separated
        bam_sv_calling_status = bam_sv_calling.branch{
            normal: it[0].status == 0
            tumor:  it[0].status == 1
        }


        // All normal samples
        bam_sv_calling_normal_for_crossing = bam_sv_calling_status.normal.map{ meta, bam, bai -> [ meta.patient, meta, bam, bai ] }

        // All tumor samples
        bam_sv_calling_tumor_for_crossing = bam_sv_calling_status.tumor.map{ meta, bam, bai -> [ meta.patient, meta, bam, bai ] }

        // Crossing the normal and tumor samples to create tumor and normal pairs
        bam_sv_calling_pair = bam_sv_calling_normal_for_crossing.cross(bam_sv_calling_tumor_for_crossing)
            .map { normal, tumor ->
                def meta = [:]

                meta.id         = "${tumor[1].sample}_vs_${normal[1].sample}".toString()
                meta.normal_id  = normal[1].sample
                meta.patient    = normal[0]
                meta.sex        = normal[1].sex
                meta.tumor_id   = tumor[1].sample

                [ meta, normal[2], normal[3], tumor[2], tumor[3] ]
        }

        //bam_sv_calling_pair.view()
        if (tools_used.contains('svaba')) {
            BAM_SVCALLING_SVABA(bam_sv_calling_pair, fasta, fasta_fai, bwa, dbsnp, dbsnp_tbi, indel_mask, germ_sv_db, simple_seq_db, error_rate)

            versions = versions.mix(BAM_SVCALLING_SVABA.out.versions)

            all_sv_vcfs = Channel.empty()
            all_sv_vcfs = all_sv_vcfs.mix(BAM_SVCALLING_SVABA.out.all_output)                     //This one contains multiple files of vcf, to get individual files, call individual output

            vcf_from_sv_calling = Channel.empty().mix(BAM_SVCALLING_SVABA.out.som_sv)

            j_supp = Channel.empty()
            j_supp = j_supp.mix(BAM_SVCALLING_SVABA.out.unfiltered_som_sv)
        }

        if (tools_used.contains('gridss')) {
            BAM_SVCALLING_GRIDSS(bam_sv_calling_pair, fasta, fasta_fai, bwa, blacklist_gridss)     // running GRIDSS

            versions = versions.mix(BAM_SVCALLING_GRIDSS.out.versions)

            vcf_from_gridss_gridss = Channel.empty()
            vcf_from_gridss_gridss = vcf_from_gridss_gridss.mix(BAM_SVCALLING_GRIDSS.out.vcf)             // This one contain only one vcf


            BAM_SVCALLING_GRIDSS_SOMATIC(vcf_from_gridss_gridss, pon_gridss)               //running the somatic filter for GRIDSS
            versions = versions.mix(BAM_SVCALLING_GRIDSS_SOMATIC.out.versions)

            vcf_from_sv_calling_gridss          = Channel.empty().mix(BAM_SVCALLING_GRIDSS_SOMATIC.out.somatic_high_confidence)
            unfiltered_som_sv_gridss          = Channel.empty().mix(BAM_SVCALLING_GRIDSS_SOMATIC.out.somatic_all)
        }

        // TODO: CHANNEL_SVCALLING_CREATE_CSV(vcf_from_sv_calling, params.tools, params.outdir) // Need to fix this!!!!!
        bam_fragcounter_calling = bam_sv_calling
    }

    if (runFragCounter) {

        if (params.step == 'fragcounter') {
            input_fragcounter_convert = input_sample.branch{
                bam:  it[0].data_type == "bam"
                cram: it[0].data_type == "cram"
            }
            // CRAM files first must be converted to BAM files since from this step on we base everything on BAM format
            CRAM_TO_BAM(input_fragcounter_convert.cram, fasta, fasta_fai)
            versions = versions.mix(CRAM_TO_BAM.out.versions)

            bam_fragcounter_calling = Channel.empty().mix(CRAM_TO_BAM.out.alignment_index, input_fragcounter_convert.bam)
                                        .map{ meta, bam, bai -> [ meta + [data_type: "bam"], bam, bai ] }           //making sure that the input data_type is correct

        }

        // getting the tumor and normal bam files separated
        bam_fragcounter_status = bam_fragcounter_calling.branch{
            normal: it[0].status == 0
            tumor:  it[0].status == 1
        }

        if (tools_used.contains('fragcounter')) {
            NORMAL_FRAGCOUNTER(bam_fragcounter_status.normal, midpoint_frag, windowsize_frag, gcmapdir_frag, minmapq_frag, paired_frag, exome_frag)
            normal_frag_cov = Channel.empty().mix(NORMAL_FRAGCOUNTER.out.rebinned_raw_cov)

            TUMOR_FRAGCOUNTER(bam_fragcounter_status.tumor, midpoint_frag, windowsize_frag, gcmapdir_frag, minmapq_frag, paired_frag, exome_frag)
            tumor_frag_cov = Channel.empty().mix(TUMOR_FRAGCOUNTER.out.rebinned_raw_cov)

            // Only need one versions because its just one program (fragcounter)
            versions = versions.mix(NORMAL_FRAGCOUNTER.out.versions)
        }

        bam_hetpileups_calling = bam_fragcounter_calling
            // TODO: Add a subworkflow to write the output file paths into a csv
    }

    if (runHetPileups) {

        if (params.step == 'hetpileups') {
            input_hetpileups_calling_convert = input_sample.branch{
                bam:  it[0].data_type == "bam"
                cram: it[0].data_type == "cram"
            }
            CRAM_TO_BAM(input_hetpileups_calling_convert.cram, fasta, fasta_fai)
            versions = versions.mix(CRAM_TO_BAM.out.versions)

            bam_hetpileups_calling = Channel.empty().mix(CRAM_TO_BAM.out.alignment_index, input_hetpileups_calling_convert.bam)
                                        .map{ meta, bam, bai -> [ meta + [data_type: "bam"], bam, bai ] }           //making sure that the input data_type is correct
        }

        // getting the tumor and normal cram files separated
        bam_hetpileups_status = bam_hetpileups_calling.branch{
            normal: it[0].status == 0
            tumor:  it[0].status == 1
        }

        // All normal samples
        bam_hetpileups_normal_for_crossing = bam_hetpileups_status.normal.map{ meta, bam, bai -> [ meta.patient, meta + [id: meta.sample], bam, bai ] }

        // All tumor samples
        bam_hetpileups_tumor_for_crossing = bam_hetpileups_status.tumor.map{ meta, bam, bai -> [ meta.patient, meta + [id: meta.sample], bam, bai ] }

        // Crossing the normal and tumor samples to create tumor and normal pairs
        bam_hetpileups_pair = bam_hetpileups_normal_for_crossing.cross(bam_hetpileups_tumor_for_crossing)
            .map { normal, tumor ->
                def meta = [:]
                meta.id         = "${tumor[1].sample}_vs_${normal[1].sample}".toString()
                meta.normal_id  = normal[1].sample
                meta.patient    = normal[0]
                meta.sex        = normal[1].sex
                meta.tumor_id   = tumor[1].sample

                [ meta, normal[2], normal[3], tumor[2], tumor[3] ]
        }

        if (tools_used.contains('hetpileups')) {

            BAM_HETPILEUPS(bam_hetpileups_pair, filter_hets, max_depth, hapmap_sites)
            versions = versions.mix(BAM_HETPILEUPS.out.versions)

            sites_from_het_pileups_wgs = Channel.empty().mix(BAM_HETPILEUPS.out.het_pileups_wgs)

        }

    }

    if (runDryClean) {

        if (params.step == 'dryclean') {
            input_dryclean_status = input_sample.branch{
                normal: it[0].status == 0
                tumor:  it[0].status == 1
            }

            tumor_frag_cov  = Channel.empty().mix(input_dryclean_status.tumor)
            normal_frag_cov = Channel.empty().mix(input_dryclean_status.normal)

        }
        if (tools_used.contains('dryclean')) {
            // Dryclean for both tumor & normal
            TUMOR_DRYCLEAN(
                    tumor_frag_cov,
                    pon_dryclean,
                    centered_dryclean,
                    cbs_dryclean,
                    cnsignif_dryclean,
                    wholeGenome_dryclean,
                    blacklist_dryclean,
                    blacklist_path_dryclean,
                    germline_filter_dryclean,
                    germline_file_dryclean,
                    human_dryclean,
                    field_dryclean,
                    build_dryclean
                    )

            tumor_dryclean_cov = Channel.empty().mix(TUMOR_DRYCLEAN.out.dryclean_cov)

            NORMAL_DRYCLEAN(
                    normal_frag_cov,
                    pon_dryclean,
                    centered_dryclean,
                    cbs_dryclean,
                    cnsignif_dryclean,
                    wholeGenome_dryclean,
                    blacklist_dryclean,
                    blacklist_path_dryclean,
                    germline_filter_dryclean,
                    germline_file_dryclean,
                    human_dryclean,
                    field_dryclean,
                    build_dryclean
                    )

            normal_dryclean_cov = Channel.empty().mix(NORMAL_DRYCLEAN.out.dryclean_cov)

            // Only need one versions because it's one program (dryclean)
            versions = versions.mix(TUMOR_DRYCLEAN.out.versions)

        }

        // TODO: Add a subworkflow to write the output file paths into a csv

        if (tools_used.contains('cbs')) {

            tumor_dryclean_cov_for_crossing = tumor_dryclean_cov.map{ meta, cov -> [ meta.patient, meta, cov ] }
            normal_dryclean_cov_for_crossing = normal_dryclean_cov.map{ meta, cov -> [ meta.patient, meta, cov ] }
            cov_cbs = tumor_dryclean_cov_for_crossing.cross(normal_dryclean_cov_for_crossing)
                .map { tumor, normal ->
                    def meta = [:]
                        meta.id             = "${tumor[1].sample}_vs_${normal[1].sample}".toString()
                        meta.sample         = "${tumor[1].sample}".toString()
                        meta.normal_id      = normal[1].sample
                        meta.patient        = normal[0]
                        meta.sex            = normal[1].sex
                        meta.tumor_id       = tumor[1].sample

                        [ meta, tumor[2], normal[2] ]
                }

            CBS(cov_cbs, cnsignif_cbs, field_cbs, name_cbs)

            versions       = versions.mix(CBS.out.versions)
            cbs_cov_rds    = Channel.empty().mix(CBS.out.cbs_cov_rds)
            cbs_seg_rds    = Channel.empty().mix(CBS.out.cbs_seg_rds)
            cbs_nseg_rds   = Channel.empty().mix(CBS.out.cbs_nseg_rds)

        }
    }

    if (runAscat) {
        if (tools_used.contains('ascat')) {

            if (params.step == 'ascat') {
                ascat_inputs = input_sample
            } else {

                het_pileups_for_crossing = sites_from_het_pileups_wgs.map{ meta, hets -> [ meta.patient, meta, hets ] }
                ascat_inputs = tumor_dryclean_cov_for_crossing
                    .cross(het_pileups_for_crossing)
                    .map { cov, hets ->
                        def meta = [:]
                            meta.id             = "${cov[1].sample}".toString()
                            meta.patient        = cov[0]
                            meta.sex            = cov[1].sex

                            [ meta, cov[2], hets[2] ]
                    }
            }

            COV_ASCAT(
                ascat_inputs,
                field_ascat,
                hets_thresh_ascat,
                penalty_ascat,
                gc_correct_ascat,
                rebin_width_ascat,
                from_maf_ascat
            )

            versions = versions.mix(COV_ASCAT.out.versions)
            purityploidy = Channel.empty().mix(COV_ASCAT.out.pp)
            purity = Channel.empty().mix(COV_ASCAT.out.purity)
            ploidy = Channel.empty().mix(COV_ASCAT.out.ploidy)
        }
    }

    if (runJabba) {
        if (tools_used.contains('jabba')) {

            if (params.step == 'jabba') {
                // put all the inputs into a map for easier retrieval
                jabba_input_map = input_sample
                    .map{ meta, cov, hets, ploidy, junction, j_supp, seg, nseg -> [
                            "meta": meta,
                            "junction": junction,
                            "cov": cov,
                            "j_supp": j_supp,
                            "hets": hets,
                            "ploidy": ploidy,
                            "seg": seg,
                            "nseg": nseg
                        ]
                    }

                // Returns a channel with the extracted ploidy when provided with a
                // file, and otherwise just passes the value through.
                jabba_input_ploidy = jabba_input_map.map{ sample ->
                    if (sample.ploidy instanceof File) {
                        EXTRACT_PURITYPLOIDY(sample.meta, sample.ploidy)
                        return [ EXTRACT_PURITYPLOIDY.out.ploidy_val ]
                    } else {
                        // If ploidy is a value, use it directly
                        return [ sample.meta, sample.ploidy ]
                    }
                }

                jabba_inputs = input_sample.join(jabba_input_ploidy)
                    .map{ meta, cov, hets, ploidy, junction, j_supp, seg, nseg, extracted_ploidy ->
                        [
                            meta,
                            junction,
                            cov,
                            j_supp,
                            hets,
                            extracted_ploidy,
                            seg,
                            nseg
                        ]
                    }

                JABBA(
                    jabba_inputs
                    blacklist_junctions_jabba,
                    geno_jabba,
                    indel_jabba,
                    tfield_jabba,
                    iter_jabba,
                    rescue_window_jabba,
                    rescue_all_jabba,
                    nudgebalanced_jabba,
                    edgenudge_jabba,
                    strict_jabba,
                    allin_jabba,
                    field_jabba,
                    maxna_jabba,
                    blacklist_coverage_jabba,
                    purity_jabba,
                    pp_method_jabba,
                    cnsignif_jabba,
                    slack_jabba,
                    linear_jabba,
                    tilim_jabba,
                    epgap_jabba,
                    fix_thres_jabba,
                    lp_jabba,
                    ism_jabba,
                    filter_loose_jabba,
                    gurobi_jabba,
                    verbose_jabba
                )

                jabba_rds           = Channel.empty().mix(JABBA.out.jabba_rds)
                jabba_gg            = Channel.empty().mix(JABBA.out.jabba_gg)
                jabba_vcf           = Channel.empty().mix(JABBA.out.jabba_vcf)
                jabba_raw_rds       = Channel.empty().mix(JABBA.out.jabba_raw_rds)
                opti                = Channel.empty().mix(JABBA.out.opti)
                jabba_seg           = Channel.empty().mix(JABBA.out.jabba_seg)
                karyograph          = Channel.empty().mix(JABBA.out.karyograph)
                versions            = versions.mix(JABBA.out.versions)

            } else {

                tumor_dryclean_cov_for_joining = tumor_dryclean_cov.map { meta, tumor_cov -> [meta.patient, meta, tumor_cov] }

                het_pileups_for_joining = sites_from_het_pileups_wgs.map { meta, hets -> [meta.patient, hets] }

                ploidy_for_joining = ploidy.map{ meta, ploidy -> [ meta.patient, ploidy ] }

                cbs_seg_rds_for_joining = cbs_seg_rds.map{ meta, seg -> [ meta.patient, seg ] }
                cbs_nseg_rds_for_joining = cbs_nseg_rds.map{ meta, nseg -> [ meta.patient, nseg ] }

                // join all previous outputs to be used as input for jabba
                // excluding svs since they can come from either svaba or gridss
                jabba_inputs = tumor_dryclean_cov_for_joining
                    .join(het_pileups_for_joining)
                    .join(ploidy_for_joining)
                    .join(cbs_seg_rds_for_joining)
                    .join(cbs_nseg_rds_for_joining)

                if (tools_used.contains('svaba')) {

                    vcf_from_sv_calling_for_joining = vcf_from_sv_calling.map{ meta, juction -> [ meta.patient, junction ] }
                    unfiltered_som_sv_for_joining = j_supp.map{ meta, j_supp -> [ meta.patient, j_supp ] }

                    jabba_w_svaba_inputs = jabba_inputs
                        .join(vcf_from_sv_calling_for_joining)
                        .join(unfiltered_som_sv_for_joining)
                        .map{ patient, meta, cov, hets, ploidy, seg, nseg, junction, j_supp ->
                            [
                                meta,
                                junction,
                                cov,
                                j_supp,
                                hets,
                                ploidy,
                                seg,
                                nseg
                            ]
                        }

                    JABBA_WITH_SVABA(
                        jabba_w_svaba_inputs,
                        blacklist_junctions_jabba,
                        geno_jabba,
                        indel_jabba,
                        tfield_jabba,
                        iter_jabba,
                        rescue_window_jabba,
                        rescue_all_jabba,
                        nudgebalanced_jabba,
                        edgenudge_jabba,
                        strict_jabba,
                        allin_jabba,
                        field_jabba,
                        maxna_jabba,
                        blacklist_coverage_jabba,
                        purity_jabba,
                        pp_method_jabba,
                        cnsignif_jabba,
                        slack_jabba,
                        linear_jabba,
                        tilim_jabba,
                        epgap_jabba,
                        fix_thres_jabba,
                        lp_jabba,
                        ism_jabba,
                        filter_loose_jabba,
                        gurobi_jabba,
                        verbose_jabba
                    )

                    jabba_rds_with_svaba           = Channel.empty().mix(JABBA_WITH_SVABA.out.jabba_rds)
                    jabba_gg_with_svaba            = Channel.empty().mix(JABBA_WITH_SVABA.out.jabba_gg)
                    jabba_vcf_with_svaba           = Channel.empty().mix(JABBA_WITH_SVABA.out.jabba_vcf)
                    jabba_raw_rds_with_svaba       = Channel.empty().mix(JABBA_WITH_SVABA.out.jabba_raw_rds)
                    opti_with_svaba                = Channel.empty().mix(JABBA_WITH_SVABA.out.opti)
                    jabba_seg_with_svaba           = Channel.empty().mix(JABBA_WITH_SVABA.out.jabba_seg)
                    karyograph_with_svaba          = Channel.empty().mix(JABBA_WITH_SVABA.out.karyograph)
                    versions_with_svaba            = versions.mix(JABBA_WITH_SVABA.out.versions)

                }

                if (tools_used.contains('gridss')) {

                    vcf_from_sv_calling_gridss_for_joining = vcf_from_sv_calling_gridss.map{ meta, junction -> [ meta.patient, junction ] }
                    unfiltered_som_sv_gridss_for_joining = unfiltered_som_sv_gridss.map{ meta, j_supp -> [ meta.patient, j_supp ] }

                    jabba_w_gridss_inputs = jabba_inputs
                        .join(vcf_from_sv_calling_gridss_for_joining)
                        .join(unfiltered_som_sv_gridss_for_joining)
                        .map{ patient, meta, cov, hets, ploidy, seg, nseg, junction, j_supp ->
                            [
                                meta,
                                junction,
                                cov,
                                j_supp,
                                hets,
                                ploidy,
                                seg,
                                nseg
                            ]
                        }

                    JABBA_WITH_GRIDSS(
                        jabba_w_gridss_inputs,
                        blacklist_junctions_jabba,
                        geno_jabba,
                        indel_jabba,
                        tfield_jabba,
                        iter_jabba,
                        rescue_window_jabba,
                        rescue_all_jabba,
                        nudgebalanced_jabba,
                        edgenudge_jabba,
                        strict_jabba,
                        allin_jabba,
                        field_jabba,
                        maxna_jabba,
                        blacklist_coverage_jabba,
                        purity_jabba,
                        pp_method_jabba,
                        cnsignif_jabba,
                        slack_jabba,
                        linear_jabba,
                        tilim_jabba,
                        epgap_jabba,
                        fix_thres_jabba,
                        lp_jabba,
                        ism_jabba,
                        filter_loose_jabba,
                        gurobi_jabba,
                        verbose_jabba
                    )

                    jabba_rds_with_gridss           = Channel.empty().mix(JABBA_WITH_GRIDSS.out.jabba_rds)
                    jabba_gg_with_gridss            = Channel.empty().mix(JABBA_WITH_GRIDSS.out.jabba_gg)
                    jabba_seg_with_gridss           = Channel.empty().mix(JABBA_WITH_GRIDSS.out.jabba_seg)
                    karyograph_with_gridss          = Channel.empty().mix(JABBA_WITH_GRIDSS.out.karyograph)
                    versions_with_gridss            = versions.mix(JABBA_WITH_GRIDSS.out.versions)

                }
            }
        }
    }

    if (runEvents) {
        if (tools_used.contains('events')) {
            if (params.step == 'events') {
                EVENTS(input_sample, fasta)
                versions = Channel.empty().mix(EVENTS.out.versions)

                events_output = Channel.empty().mix(EVENTS.out.events_output)
            } else {
                if (tools_used.contains('gridss')) {
                    EVENTS_WITH_GRIDSS(jabba_rds_with_gridss, fasta)
                    events_w_gridss_versions = Channel.empty().mix(EVENTS_WITH_GRIDSS.out.versions)

                    events_w_gridss_output = Channel.empty().mix(EVENTS_WITH_GRIDSS.out.events_output)
                }
                if (tools_used.contains('svaba')) {
                    EVENTS_WITH_SVABA(jabba_rds_with_svaba, fasta)
                    events_w_svaba_versions = Channel.empty().mix(EVENTS_WITH_SVABA.out.versions)

                    events_w_svaba_output = Channel.empty().mix(EVENTS_WITH_SVABA.out.events_output)
                }
            }
        }
    }

    if (runFusions) {
        if (tools_used.contains('fusions')) {
            if (params.step == 'fusions') {
                FUSIONS(input_sample, gencode_fusions)
                fusions_output = Channel.empty().mix(FUSIONS.out.fusions_output)
                versions = Channel.empty().mix(FUSIONS.out.versions)

            } else {
                if (tools_used.contains('gridss')) {
                    FUSIONS_WITH_GRIDSS(jabba_rds_with_gridss, gencode_fusions)
                    fusions_w_gridss_output = Channel.empty().mix(FUSIONS_WITH_GRIDSS.out.fusions_output)
                    fusions_w_gridss_versions = Channel.empty().mix(FUSIONS_WITH_GRIDSS.out.versions)

                }
                if (tools_used.contains('svaba')) {
                    FUSIONS_WITH_SVABA(jabba_rds_with_svaba, gencode_fusions)
                    fusions_w_svaba_output = Channel.empty().mix(FUSIONS_WITH_SVABA.out.fusions_output)
                    fusions_w_svaba_versions = Channel.empty().mix(FUSIONS_WITH_SVABA.out.versions)

                }
            }
        }
    }

    if (runAlleicCN) {
        if (tools_used.contains('allelic_cn')) {
            if (params.step == 'allelic_cn') {
                non_integer_balance_inputs = input_sample.map{ meta, cov, hets, ggraph -> [ meta, ggraph, cov, hets ] }

                NON_INTEGER_BALANCE(
                    non_integer_balance_inputs,
					field_non_integer_balance,
					hets_thresh_non_integer_balance,
                    mask_non_integer_balance,
					overwrite_non_integer_balance,
					lambda_non_integer_balance,
					allin_non_integer_balance,
					fix_thresh_non_integer_balance,
					nodebounds_non_integer_balance,
					ism_non_integer_balance,
					build_non_integer_balance,
					epgap_non_integer_balance,
					tilim_non_integer_balance,
					gurobi_non_integer_balance,
                    fasta,
					pad_non_integer_balance
                )
                versions = Channel.empty().mix(NON_INTEGER_BALANCE.out.versions)

                non_integer_balance_balanced_gg = Channel.empty().mix(NON_INTEGER_BALANCE.out.non_integer_balance_balanced_gg)
                non_integer_balance_hets_gg = Channel.empty().mix(NON_INTEGER_BALANCE.out.non_integer_balance_hets_gg)

                non_integer_balance_balanced_gg_for_joining = non_integer_balance_balanced_gg.map{ meta, balanced -> [ meta.patient, meta, balanced ] }
                hets_for_joining = input_sample.map{ meta, cov, hets, ggraph -> [ meta.patient, hets ] }
                lp_phased_balance_inputs = non_integer_balance_hets_gg_for_joining.join(hets_for_joining)
                    .map{ patient, meta, balanced, hets -> [ meta, balanced, hets ]
                }

                LP_PHASED_BALANCE(
                    lp_phased_balance_inputs,
                    lambda_lp_phased_balance,
                    cnloh_lp_phased_balance,
                    major_lp_phased_balance,
                    allin_lp_phased_balance,
                    marginal_lp_phased_balance,
                    from_maf_lp_phased_balance,
                    mask_lp_phased_balance,
                    ism_lp_phased_balance,
                    epgap_lp_phased_balance,
                    hets_thresh_lp_phased_balance,
                    min_bins_lp_phased_balance,
                    min_width_lp_phased_balance,
                    trelim_lp_phased_balance,
                    reward_lp_phased_balance,
                    nodefileind_lp_phased_balance,
                    tilim_lp_phased_balance
                )

                lp_phased_balance_balanced_gg = Channel.empty().mix(LP_PHASED_BALANCE.out.lp_phased_balance_balanced_gg)
                lp_phased_balance_binstats_gg = Channel.empty().mix(LP_PHASED_BALANCE.out.lp_phased_balance_binstats_gg)
                lp_phased_balance_unphased_allelic_gg = Channel.empty().mix(LP_PHASED_BALANCE.out.lp_phased_balance_unphased_allelic_gg)
            } else {
                if (tools_used.contains('gridss')) {
                    jabba_rds_with_gridss_for_joining = jabba_rds_with_gridss.map{ meta, rds -> [ meta.patient, meta, rds ] }
                    non_integer_balance_w_gridss_inputs = jabba_rds_with_gridss_for_joining
                        .join(tumor_dryclean_cov_for_joining)
                        .join(het_pileups_for_joining)
                        .map{ patient, meta, rds, cov, hets -> [ meta, rds, cov, hets ] }

                    NON_INTEGER_BALANCE_WITH_GRIDSS(
                        non_integer_balance_w_gridss_inputs,
                        field_non_integer_balance,
                        hets_thresh_non_integer_balance,
                        mask_non_integer_balance,
                        overwrite_non_integer_balance,
                        lambda_non_integer_balance,
                        allin_non_integer_balance,
                        fix_thresh_non_integer_balance,
                        nodebounds_non_integer_balance,
                        ism_non_integer_balance,
                        build_non_integer_balance,
                        epgap_non_integer_balance,
                        tilim_non_integer_balance,
                        gurobi_non_integer_balance,
                        fasta,
                        pad_non_integer_balance
                    )
                    versions_w_gridss = Channel.empty().mix(NON_INTEGER_BALANCE_WITH_GRIDSS.out.versions)

                    non_integer_balance_w_gridss_balanced_gg = Channel.empty().mix(NON_INTEGER_BALANCE_WITH_GRIDSS.out.non_integer_balance_balanced_gg)
                    non_integer_balance_w_gridss_hets_gg = Channel.empty().mix(NON_INTEGER_BALANCE_WITH_GRIDSS.out.non_integer_balance_hets_gg)

                    non_integer_balance_w_gridss_balanced_gg_for_joining = non_integer_balance_w_gridss_balanced_gg.map{ meta, balanced -> [ meta.patient, meta, balanced ] }
                    lp_phased_balance_w_gridss_inputs = non_integer_balance_w_gridss_balanced_gg_for_joining.join(het_pileups_for_joining)
                        .map{ patient, meta, balanced, hets -> [ meta, balanced, hets ] }

                    LP_PHASED_BALANCE_WITH_GRIDSS(
                        lp_phased_balance_w_gridss_inputs,
                        lambda_lp_phased_balance,
                        cnloh_lp_phased_balance,
                        major_lp_phased_balance,
                        allin_lp_phased_balance,
                        marginal_lp_phased_balance,
                        from_maf_lp_phased_balance,
                        mask_lp_phased_balance,
                        ism_lp_phased_balance,
                        epgap_lp_phased_balance,
                        hets_thresh_lp_phased_balance,
                        min_bins_lp_phased_balance,
                        min_width_lp_phased_balance,
                        trelim_lp_phased_balance,
                        reward_lp_phased_balance,
                        nodefileind_lp_phased_balance,
                        tilim_lp_phased_balance
                    )

                    lp_phased_balance_w_gridss_balanced_gg = Channel.empty().mix(LP_PHASED_BALANCE_WITH_GRIDSS.out.lp_phased_balance_balanced_gg)
                    lp_phased_balance_w_gridss_binstats_gg = Channel.empty().mix(LP_PHASED_BALANCE_WITH_GRIDSS.out.lp_phased_balance_binstats_gg)
                    lp_phased_balance_w_gridss_unphased_allelic_gg = Channel.empty().mix(LP_PHASED_BALANCE_WITH_GRIDSS.out.lp_phased_balance_unphased_allelic_gg)
                }
                if (tools_used.contains('svaba')) {
                    jabba_rds_with_svaba_for_joining = jabba_rds_with_svaba.map{ meta, rds -> [ meta.patient, meta, rds ] }
                    non_integer_balance_w_svaba_inputs = jabba_rds_with_svaba_for_joining
                        .join(tumor_dryclean_cov_for_joining)
                        .join(het_pileups_for_joining)
                        .map{ patient, meta, rds, cov, hets -> [ meta, rds, cov, hets ] }
                    NON_INTEGER_BALANCE_WITH_SVABA(
                        non_integer_balance_w_svaba_inputs,
                        field_non_integer_balance,
                        hets_thresh_non_integer_balance,
                        mask_non_integer_balance,
                        overwrite_non_integer_balance,
                        lambda_non_integer_balance,
                        allin_non_integer_balance,
                        fix_thresh_non_integer_balance,
                        nodebounds_non_integer_balance,
                        ism_non_integer_balance,
                        build_non_integer_balance,
                        epgap_non_integer_balance,
                        tilim_non_integer_balance,
                        gurobi_non_integer_balance,
                        fasta,
                        pad_non_integer_balance
                    )
                    versions_w_svaba = Channel.empty().mix(NON_INTEGER_BALANCE_WITH_SVABA.out.versions)

                    non_integer_balance_w_svaba_balanced_gg = Channel.empty().mix(NON_INTEGER_BALANCE_WITH_SVABA.out.non_integer_balance_balanced_gg)
                    non_integer_balance_w_svaba_hets_gg = Channel.empty().mix(NON_INTEGER_BALANCE_WITH_SVABA.out.non_integer_balance_hets_gg)

                    non_integer_balance_w_svaba_balanced_gg_for_joining = non_integer_balance_w_svaba_balanced_gg.map{ meta, balanced -> [ meta.patient, meta, balanced ] }
                    lp_phased_balance_w_svaba_inputs = non_integer_balance_w_svaba_balanced_gg_for_joining.join(het_pileups_for_joining)
                        .map{ patient, meta, balanced, hets -> [ meta, balanced, hets ] }

                    LP_PHASED_BALANCE_WITH_SVABA(
                        lp_phased_balance_w_svaba_inputs,
                        lambda_lp_phased_balance,
                        cnloh_lp_phased_balance,
                        major_lp_phased_balance,
                        allin_lp_phased_balance,
                        marginal_lp_phased_balance,
                        from_maf_lp_phased_balance,
                        mask_lp_phased_balance,
                        ism_lp_phased_balance,
                        epgap_lp_phased_balance,
                        hets_thresh_lp_phased_balance,
                        min_bins_lp_phased_balance,
                        min_width_lp_phased_balance,
                        trelim_lp_phased_balance,
                        reward_lp_phased_balance,
                        nodefileind_lp_phased_balance,
                        tilim_lp_phased_balance
                    )

                    lp_phased_balance_balanced_gg = Channel.empty().mix(LP_PHASED_BALANCE_WITH_SVABA.out.lp_phased_balance_balanced_gg)
                    lp_phased_balance_binstats_gg = Channel.empty().mix(LP_PHASED_BALANCE_WITH_SVABA.out.lp_phased_balance_binstats_gg)
                    lp_phased_balance_unphased_allelic_gg = Channel.empty().mix(LP_PHASED_BALANCE_WITH_SVABA.out.lp_phased_balance_unphased_allelic_gg)
                }
            }
        }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
