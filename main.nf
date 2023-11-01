#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/heisenbio
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://tanubrata/nf-core/heisenbio
    Website: https://nf-co.re/heisenbio
    Slack  : https://nfcore.slack.com/channels/heisenbio
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.ascat_alleles         = WorkflowMain.getGenomeAttribute(params, 'ascat_alleles')
params.ascat_genome          = WorkflowMain.getGenomeAttribute(params, 'ascat_genome')
params.ascat_loci            = WorkflowMain.getGenomeAttribute(params, 'ascat_loci')
params.ascat_loci_gc         = WorkflowMain.getGenomeAttribute(params, 'ascat_loci_gc')
params.ascat_loci_rt         = WorkflowMain.getGenomeAttribute(params, 'ascat_loci_rt')
params.bwa                   = WorkflowMain.getGenomeAttribute(params, 'bwa')
params.bwamem2               = WorkflowMain.getGenomeAttribute(params, 'bwamem2')
params.cf_chrom_len          = WorkflowMain.getGenomeAttribute(params, 'cf_chrom_len')
params.chr_dir               = WorkflowMain.getGenomeAttribute(params, 'chr_dir')
params.dbsnp                 = WorkflowMain.getGenomeAttribute(params, 'dbsnp')
params.dbsnp_tbi             = WorkflowMain.getGenomeAttribute(params, 'dbsnp_tbi')
params.dbsnp_vqsr            = WorkflowMain.getGenomeAttribute(params, 'dbsnp_vqsr')
params.dict                  = WorkflowMain.getGenomeAttribute(params, 'dict')
//params.dragmap               = WorkflowMain.getGenomeAttribute(params, 'dragmap')
params.fasta                 = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.fasta_fai             = WorkflowMain.getGenomeAttribute(params, 'fasta_fai')
//params.germline_resource     = WorkflowMain.getGenomeAttribute(params, 'germline_resource')
//params.germline_resource_tbi = WorkflowMain.getGenomeAttribute(params, 'germline_resource_tbi')
params.intervals             = WorkflowMain.getGenomeAttribute(params, 'intervals')
params.known_snps            = WorkflowMain.getGenomeAttribute(params, 'known_snps')
params.known_snps_tbi        = WorkflowMain.getGenomeAttribute(params, 'known_snps_tbi')
params.known_snps_vqsr       = WorkflowMain.getGenomeAttribute(params, 'known_snps_vqsr')
params.known_indels          = WorkflowMain.getGenomeAttribute(params, 'known_indels')
params.known_indels_tbi      = WorkflowMain.getGenomeAttribute(params, 'known_indels_tbi')
params.known_indels_vqsr     = WorkflowMain.getGenomeAttribute(params, 'known_indels_vqsr')
//params.mappability           = WorkflowMain.getGenomeAttribute(params, 'mappability')
//params.pon                   = WorkflowMain.getGenomeAttribute(params, 'pon')
//params.pon_tbi               = WorkflowMain.getGenomeAttribute(params, 'pon_tbi')
//params.snpeff_db             = WorkflowMain.getGenomeAttribute(params, 'snpeff_db')
//params.snpeff_genome         = WorkflowMain.getGenomeAttribute(params, 'snpeff_genome')
//params.vep_cache_version     = WorkflowMain.getGenomeAttribute(params, 'vep_cache_version')
//params.vep_genome            = WorkflowMain.getGenomeAttribute(params, 'vep_genome')
//params.vep_species           = WorkflowMain.getGenomeAttribute(params, 'vep_species')
params.indel_mask            = WorkflowMain.getGenomeAttribute(params, 'indel_mask')
params.germ_sv_db            = WorkflowMain.getGenomeAttribute(params, 'germ_sv_db')
params.simple_seq_db         = WorkflowMain.getGenomeAttribute(params, 'simple_seq_db')
params.blacklist_gridss      = WorkflowMain.getGenomeAttribute(params, 'blacklist_gridss')
params.pon_gridss            = WorkflowMain.getGenomeAttribute(params, 'pon_gridss')
params.gcmapdir_frag         = WorkflowMain.getGenomeAttribute(params, 'gcmapdir_frag')
params.hapmap_sites          = WorkflowMain.getGenomeAttribute(params, 'hapmap_sites')
params.pon_dryclean          = WorkflowMain.getGenomeAttribute(params, 'pon_dryclean')
params.blacklist_coverage_jabba     = WorkflowMain.getGenomeAttribute(params, 'blacklist_coverage_jabba')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ALTERNATIVE INPUT FILE ON RESTART
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.input_restart = WorkflowHeisenbio.retrieveInput(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsHelp } from 'plugin/nf-validation'

// Print help message if needed
if (params.help) {
    def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
    def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
    def String command = "nextflow run ${workflow.manifest.name} --input samplesheet.csv --genome GRCh37 -profile docker"
    log.info logo + paramsHelp(command) + citation + NfcoreTemplate.dashedLine(params.monochrome_logs)
    System.exit(0)
}

// Validate input parameters
if (params.validate_params) {
    validateParameters()
}

WorkflowMain.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HEISENBIO } from './workflows/heisenbio'

//
// WORKFLOW: Run main nf-core/heisenbio analysis pipeline
//
workflow NFCORE_HEISENBIO {
    HEISENBIO ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_HEISENBIO ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
