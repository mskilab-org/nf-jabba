//
// ASCAT
//

include { ASCAT } from '../../../modules/nf-core/ascat/main'


workflow BAM_ASCAT {

    take:
    cram_pair                                         // channel: [mandatory] [ meta, normalcram, normalcrai, tumorcram, tumorcrai ]
    allele_files                                      // channel: [mandatory] zip
    loci_files                                        // channel: [mandatory] zip 
    intervals_bed                                     // channel: [optional]  bed for WES
    fasta                                             // channel: [optional]  fasta needed for cram
    gc_file                                           // channel: [optional]  txt for LogRCorrection
    rt_file                                           // channel: [optional]  txt for LogRCorrection

    main:
    versions     = Channel.empty()
    pp           = Channel.empty()
    baf          = Channel.empty()
    cnv          = Channel.empty()
    logR         = Channel.empty()
    metric       = Channel.empty()
    segment      = Channel.empty()
    alleleFreq   = Channel.empty()

    if (!params.wes) intervals_bed = []                // No intervals needed if not WES
    ASCAT(cram_pair, allele_files, loci_files, intervals_bed, fasta, gc_file, rt_file)


    versions     = versions.mix(ASCAT.out.versions)
    pp           = pp.mix(ASCAT.out.purityploidy)
    baf          = baf.mix(ASCAT.out.bafs)
    cnv          = cnv.mix(ASCAT.out.cnvs)
    logR         = logR.mix(ASCAT.out.logrs)
    metric       = metric.mix(ASCAT.out.metrics)
    segment      = segment.mix(ASCAT.out.segments)
    alleleFreq   = segment.mix(ASCAT.out.allelefreqs)


    emit:
    pp                                                 // only need to emit the purityploidy for JaBbA

    versions
}