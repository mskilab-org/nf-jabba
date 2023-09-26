//
// MAPPING
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { BWAMEM2_MEM            } from '../../../modules/nf-core/bwamem2/mem/main'
include { BWA_MEM as BWAMEM1_MEM } from '../../../modules/nf-core/bwa/mem/main'

workflow FASTQ_ALIGN_BWAMEM_MEM2 {
    take:
    reads // channel: [mandatory] meta, reads
    index // channel: [mandatory] index
    sort  // boolean: [mandatory] true -> sort, false -> don't sort
    fasta
    fasta_fai

    main:

    versions = Channel.empty()
    reports = Channel.empty()

    // Only one of the following should be run
    BWAMEM1_MEM(reads, index.map{ it -> [ [ id:'index' ], it ] }, sort) // If aligner is bwa-mem
    BWAMEM2_MEM(reads, index.map{ it -> [ [ id:'index' ], it ] }, sort) // If aligner is bwa-mem2
    
    // Get the bam files from the aligner
    // Only one aligner is run
    bam = Channel.empty()
    bam = bam.mix(BWAMEM1_MEM.out.bam)
    bam = bam.mix(BWAMEM2_MEM.out.bam)


    // Gather reports of all tools used
    reports = reports.mix(BWAMEM2_MEM.out.log)
    reports = reports.mix(BWAMEM1_MEM.out.log)
   

    // Gather versions of all tools used
    versions = versions.mix(BWAMEM1_MEM.out.versions)
    versions = versions.mix(BWAMEM2_MEM.out.versions)

    emit:
    bam      // channel: [ [meta], bam ]
    bai      // channel: [ [meta], bai ]
    reports
    versions // channel: [ versions.yml ]
}
