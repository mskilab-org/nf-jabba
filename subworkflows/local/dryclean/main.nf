//
// DRYCLEAN
//

include { DRYCLEAN } from '../../../modules/local/dryclean/main.nf'

workflow COV_DRYCLEAN {

    take:
    input_dryclean   // channel: [mandatory] [ meta, input ]
    pon_dryclean
    centered_dryclean
    cbs_dryclean
    cnsignif_dryclean
    wholeGenome_dryclean
    blacklist_dryclean
    blacklist_path_dryclean
    germline_filter_dryclean
    germline_file_dryclean
    human_dryclean
    field_dryclean
    build_dryclean

    main:
    versions          = Channel.empty()
    dryclean_cov      = Channel.empty()
    dryclean_obj      = Channel.empty()

    DRYCLEAN(input_dryclean, pon_dryclean, centered_dryclean, cbs_dryclean,
    cnsignif_dryclean, wholeGenome_dryclean, blacklist_dryclean,
    blacklist_path_dryclean, germline_filter_dryclean, germline_file_dryclean,
    human_dryclean, field_dryclean, build_dryclean)

    dryclean_cov      = DRYCLEAN.out.decomposed_cov
    dryclean_obj      = DRYCLEAN.out.dryclean_object

    versions          = DRYCLEAN.out.versions

    emit:
    dryclean_cov    // only need to emit the coverage for JaBbA

    versions
}

