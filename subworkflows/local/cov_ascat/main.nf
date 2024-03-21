//
// ASCAT
//

include { ASCAT_SEG } from '../../../modules/local/ascat/main'
include { EXTRACT_PURITYPLOIDY } from '../../../modules/local/ascat/main'


workflow COV_ASCAT {

    take:
    inputs              // [ meta, hets, cbs_cov ]
    field
    hets_threshold
    penalty
    gc_correct
    rebin_width
    from_maf


    main:
    versions     = Channel.empty()
    pp           = Channel.empty()
    seg          = Channel.empty()
    purity       = Channel.empty()
    ploidy       = Channel.empty()

    ASCAT_SEG(inputs, field, hets_threshold, penalty, gc_correct, rebin_width, from_maf)

    ascat_pp = Channel.empty().mix(ASCAT_SEG.out.purityploidy)
    EXTRACT_PURITYPLOIDY(ascat_pp)

    versions     = versions.mix(ASCAT_SEG.out.versions)
    pp           = pp.mix(ASCAT_SEG.out.purityploidy)
    seg          = seg.mix(ASCAT_SEG.out.segments)
    purity       = purity.mix(EXTRACT_PURITYPLOIDY.out.purity_val)
    ploidy       = ploidy.mix(EXTRACT_PURITYPLOIDY.out.ploidy_val)

    emit:
    pp                                                 // only need to emit the purityploidy for JaBbA
    purity
    ploidy

    versions
}
