//
// BAM DRYCLEAN
//

include { DRYCLEAN } from '../../../modules/local/dryclean/main.nf'

workflow BAM_DRYCLEAN {

    take:
    tuple val(meta), path(cov)
    path(pon_dryclean)
    val(germline_filter_drycln)
    path(blacklist_drycln)
    path(germline_file_drycln)
    val(collapse_drycln)
    val(field_drycln)
    tuple val(meta), path(bam), path(bai)

    main:
    versions          = Channel.empty()
    dryclean_cov      = Channel.empty()

    DRYCLEAN(cov, pon_dryclean, germline_filter_drycln, blacklist_drycln, germline_file_drycln, collapse_drycln, field_drycln,)





}

