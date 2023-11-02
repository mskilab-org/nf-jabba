//
// ASCAT
//

include { ASCAT_SEG } from '../../../modules/local/ascat/main'


workflow COV_ASCAT {

    take:
    hetpileups                                        // channel: [mandatory] [ meta, hets ]
    cbs_cov                                           // channel: [mandatory] [ meta, cbs_cov ]
    field                                             // channel: [mandatory] "foreground" for dryclean/ "ratio"
    hets_threshold                                    // channel: cutoff for hetpileups; default=0.2
    penalty                                           // channel: penalty for ASCAT; default=70
    gc_correct                                        // channel: perform GC correction? Default=TRUE
    rebin_width                                       // channel: width for rebinning, default=5e4
    from_maf                                          // channel: whether to start from MAF, default=FALSE


    main:
    versions     = Channel.empty()
    pp           = Channel.empty()
    seg          = Channel.empty()
    //purity       = Channel.empty()
    //ploidy       = Channel.empty()

    ASCAT_SEG(hetpileups, cbs_cov, field, hets_threshold, penalty, gc_correct, rebin_width, from_maf)


    versions     = versions.mix(ASCAT_SEG.out.versions)
    pp           = pp.mix(ASCAT_SEG.out.purityploidy)
    seg          = seg.mix(ASCAT_SEG.out.segments)
    //purity       = purity.mix(ASCAT_SEG.out.purity)
    //ploidy       = purity.mix(ASCAT_SEG.out.ploidy)

    emit:
    pp                                                 // only need to emit the purityploidy for JaBbA
    //purity
    //ploidy

    versions
}
