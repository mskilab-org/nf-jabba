//
// JaBbA
//

include { JABBA } from '../../../modules/local/jabba/main.nf'

workflow COV_JUNC_JABBA {

    take:
    cov_rds_jabba           // [ meta, cov]
    junction_jabba          // [ meta, junction]
    ploidy_jabba            // [ meta, ploidy ]
    het_pileups_wgs_jabba   // [ meta, hets ]
    cbs_seg_rds_jabba       // [ meta, seg_cbs ]
    cbs_nseg_rds_jabba      // [ meta, nseg_cbs ]
    j_supp_jabba            // [ meta, unfiltered_som_sv ]
    blacklist_junctions_jabba
    geno_jabba
    indel_jabba
    tfield_jabba
    iter_jabba
    rescue_window_jabba
    rescue_all_jabba
    nudgebalanced_jabba
    edgenudge_jabba
    strict_jabba
    allin_jabba
    field_jabba
    maxna_jabba
    blacklist_coverage_jabba
    purity_jabba
    pp_method_jabba
    cnsignif_jabba
    slack_jabba
    linear_jabba
    tilim_jabba
    epgap_jabba
    name_jabba
    fix_thres_jabba
    lp_jabba
    ism_jabba
    filter_loose_jabba
    gurobi_jabba
    nonintegral_jabba
    verbose_jabba
    help_jabba

    main:
    versions            = Channel.empty()
    jabba_rds           = Channel.empty()
    jabba_gg            = Channel.empty()
    jabba_vcf           = Channel.empty()
    jabba_raw_rds       = Channel.empty()
    opti                = Channel.empty()
    jabba_seg           = Channel.empty()
    karyograph          = Channel.empty()

    JABBA(cov_rds_jabba, junction_jabba, ploidy_jabba, het_pileups_wgs_jabba,
    cbs_seg_rds_jabba, cbs_nseg_rds_jabba, j_supp_jabba,
    blacklist_junctions_jabba, geno_jabba, indel_jabba, tfield_jabba,
    iter_jabba, rescue_window_jabba, rescue_all_jabba, nudgebalanced_jabba,
    edgenudge_jabba, strict_jabba, allin_jabba, field_jabba, maxna_jabba,
    blacklist_coverage_jabba, purity_jabba, pp_method_jabba, cnsignif_jabba,
    slack_jabba, linear_jabba, tilim_jabba, epgap_jabba, name_jabba,
    fix_thres_jabba, lp_jabba, ism_jabba, filter_loose_jabba, gurobi_jabba,
    nonintegral_jabba, verbose_jabba, help_jabba)

    jabba_rds           = JABBA.out.jabba_rds
    jabba_gg            = JABBA.out.jabba_gg
    jabba_vcf           = JABBA.out.jabba_vcf
    jabba_raw_rds       = JABBA.out.jabba_raw_rds
    opti                = JABBA.out.opti
    jabba_seg           = JABBA.out.jabba_seg
    karyograph          = JABBA.out.karyograph

    versions          = JABBA.out.versions

    emit:
    jabba_rds
    jabba_gg
    jabba_vcf
    jabba_raw_rds
    opti
    jabba_seg
    karyograph

    versions
}

