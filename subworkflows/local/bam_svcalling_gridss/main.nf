//
// GRIDSS SV CALLING
//
//
//

include { GRIDSS_GRIDSS   } from '../../../modules/local/gridss/gridss/main.nf'
include { GRIDSS_SOMATIC  } from '../../../modules/local/gridss/somaticFilter/main.nf'

workflow BAM_SVCALLING_GRIDSS {
    take:
    cram                                              // channel: [mandatory] [ meta, normalcram, normalcrai, tumorcram, tumorcrai ]
    fasta                                             // channel: [mandatory] reference fasta
    fasta_fai                                         // channel: [mandatory] reference fasta index
    bwa_index                                         // channel: [mandatory] bwa index path
    blacklist_gridss                                  // optional: blacklist bed file for gridss


    main:
    versions               = Channel.empty()
    vcf                    = Channel.empty()
    vcf_index              = Channel.empty()
    assembly_bam           = Channel.empty()

    GRIDSS_GRIDSS(cram, fasta, fasta_fai, bwa_index, blacklist_gridss)

    vcf                    = GRIDSS_GRIDSS.out.vcf
    vcf_index              = GRIDSS_GRIDSS.out.vcf_index
    assembly_bam           = GRIDSS_GRIDSS.out.assembly


    versions = versions.mix(GRIDSS_GRIDSS.out.versions)

    emit:
    vcf
    vcf_index
    assembly_bam


    versions

}


workflow BAM_SVCALLING_GRIDSS_SOMATIC {
    take:
    vcf
    pondir_gridss

    main:
    versions                = Channel.empty()
    somatic_all             = Channel.empty()
    somatic_high_confidence = Channel.empty()

    GRIDSS_SOMATIC(vcf, pondir_gridss)

    somatic_all             = GRIDSS_SOMATIC.out.somatic_all_vcf
    somatic_high_confidence = GRIDSS_SOMATIC.out.somatic_high_vcf

    versions                = GRIDSS_SOMATIC.out.versions

    all_vcf = Channel.empty().mix(somatic_all, somatic_high_confidence)

    emit:
    somatic_all
    somatic_high_confidence
    all_vcf

    versions

}









