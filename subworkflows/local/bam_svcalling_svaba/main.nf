//
// SVABA SV CALLING
//
//
//

include { SVABA } from '../../../modules/local/svaba/main.nf'

workflow BAM_SVCALLING_SVABA {
    take:
    cram                                                  // channel: [mandatory] [ meta, tumorcram, tumorcrai, normalcram, normalcrai ] 
    fasta                                                 // channel: [mandatory] [ fasta ]
    fasta_fai                                             // channel: [mandatory] [ fasta index ]
    bwa_index
    dbsnp
    dbsnp_tbi
    indel_mask
    germ_sv_db
    simple_seq_db
    error_rate

    main:
    versions               = Channel.empty()
    som_sv                 = Channel.empty()
    som_indel              = Channel.empty()
    germ_sv                = Channel.empty()
    germ_indel             = Channel.empty()
    raw_calls              = Channel.empty()
    discordants            = Channel.empty()
    ascii_alignments       = Channel.empty()
    unfiltered_germ_sv     = Channel.empty()
    unfiltered_germ_indel  = Channel.empty()
    unfiltered_som_sv      = Channel.empty()
    unfiltered_som_indel   = Channel.empty()

    //Remapping the input based on Svaba module
    cram_svaba = cram.map { meta, normal_cram, normal_crai, tumor_cram, tumor_crai ->
        [meta, tumor_cram, tumor_crai, normal_cram, normal_crai]
    }
    // Calling Svaba module for run
    SVABA(cram_svaba, fasta, fasta_fai, bwa_index, dbsnp, dbsnp_tbi, indel_mask, germ_sv_db, simple_seq_db, error_rate)

    // getting the outputs in a channel from SVABA
    som_sv                 = SVABA.out.som_sv
    som_indel              = SVABA.out.som_indel
    germ_sv                = SVABA.out.germ_sv
    germ_indel             = SVABA.out.germ_indel
    raw_calls              = SVABA.out.raw_calls
    discordants            = SVABA.out.discordants
    ascii_alignments       = SVABA.out.ascii_alignments
    unfiltered_germ_sv     = SVABA.out.unfiltered_germ_sv
    unfiltered_germ_indel  = SVABA.out.unfiltered_germ_indel
    unfiltered_som_sv      = SVABA.out.unfiltered_som_sv 
    unfiltered_som_indel   = SVABA.out.unfiltered_som_indel

    versions = versions.mix(SVABA.out.versions)
    all_output = Channel.empty().mix(
       som_sv,
       som_indel,
       germ_sv,
       germ_indel,
       unfiltered_som_sv,
       unfiltered_som_indel,unfiltered_germ_sv,
       unfiltered_germ_indel 
    )

    emit:
    som_sv
    som_indel
    germ_sv
    germ_indel
    raw_calls
    discordants
    ascii_alignments
    unfiltered_germ_sv
    unfiltered_germ_indel
    unfiltered_som_sv
    unfiltered_som_indel
    all_output

    versions
}