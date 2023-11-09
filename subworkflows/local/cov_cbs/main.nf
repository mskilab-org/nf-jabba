//
// CBS
//
//

include { CBS } from '../../../modules/local/cbs/main.nf'

// Define the main workflow process
workflow COV_CBS {
    // Define the input parameters for the main workflow
    take:
    cov_cbs         // channel: [mandatory] [ meta, tumor_cov, normal_cov ]
    cnsignif_cbs
    field_cbs
    name_cbs

    main:
    cbs_cov_rds             = Channel.empty()
    cbs_seg_rds             = Channel.empty()
    cbs_nseg_rds            = Channel.empty()
    versions                = Channel.empty()

    CBS(cov_cbs, cnsignif_cbs, field_cbs, name_cbs)

    cbs_cov_rds              = CBS.out.cbs_cov_rds
    cbs_seg_rds              = CBS.out.cbs_seg_rds
    cbs_nseg_rds             = CBS.out.cbs_nseg_rds

    versions = versions.mix(CBS.out.versions)

    emit:
    cbs_cov_rds
    cbs_seg_rds
    cbs_nseg_rds

    versions
}
