//
// JaBbA
//

include { JABBA } from '../../../modules/local/jabba/main.nf'
include { COERCE_SEQNAMES as COERCE_SEQNAMES_COV } from '../../../modules/local/jabba/main.nf'
include { COERCE_SEQNAMES as COERCE_SEQNAMES_SOM_SV } from '../../../modules/local/jabba/main.nf'
include { COERCE_SEQNAMES as COERCE_SEQNAMES_UNFIL_SOM_SV } from '../../../modules/local/jabba/main.nf'
include { COERCE_SEQNAMES as COERCE_SEQNAMES_HETS } from '../../../modules/local/jabba/main.nf'

workflow COV_JUNC_JABBA {

    take:
    cov_rds_jabba           // [ meta, cov]
    junction_jabba          // [ meta, junction]
    ploidy_jabba            // [ meta, ploidy ]
    het_pileups_wgs_jabba   // [ meta, hets ]
    cbs_seg_rds_jabba       // [ meta, seg_cbs ]
    cbs_nseg_rds_jabba      // [ meta, nseg_cbs ]
    j_supp_jabba            // [ meta, unfiltered_som_sv ]
    blacklist_junctions_jabba // this is a val, not a path, but is treated like a path to allow for a "NULL" default value
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
    fix_thres_jabba
    lp_jabba
    ism_jabba
    filter_loose_jabba
    gurobi_jabba
    verbose_jabba

    main:
    versions            = Channel.empty()
    jabba_rds           = Channel.empty()
    jabba_gg            = Channel.empty()
    jabba_vcf           = Channel.empty()
    jabba_raw_rds       = Channel.empty()
    opti                = Channel.empty()
    jabba_seg           = Channel.empty()
    karyograph          = Channel.empty()

    // Add channels for the outputs of COERCE_SEQNAMES
    chr_coerced_cov_rds_jabba          = Channel.empty()
    chr_coerced_junction_jabba         = Channel.empty()
    chr_coerced_j_supp_jabba           = Channel.empty()
    chr_coerced_het_pileups_wgs_jabba  = Channel.empty()


    // Run COERCE_SEQNAMES to force inputs to be in common
    COERCE_SEQNAMES_COV(cov_rds_jabba)
    chr_coerced_cov_rds_jabba = COERCE_SEQNAMES_COV.out.file
    chr_coerced_cov_rds_jabba_to_cross = chr_coerced_cov_rds_jabba.map { tuple ->
                                                            def (meta, cov) = tuple
                                                            [meta.patient, meta, cov] }
    COERCE_SEQNAMES_SOM_SV(junction_jabba)
    chr_coerced_junction_jabba = COERCE_SEQNAMES_SOM_SV.out.file
    chr_coerced_junction_jabba_to_cross = chr_coerced_junction_jabba.map { tuple ->
                                                            def (meta, vcf) = tuple
                                                            [meta.patient, meta, vcf] }

    COERCE_SEQNAMES_UNFIL_SOM_SV(j_supp_jabba)
    chr_coerced_j_supp_jabba = COERCE_SEQNAMES_UNFIL_SOM_SV.out.file
    chr_coerced_j_supp_jabba_to_cross = chr_coerced_j_supp_jabba.map { tuple ->
                                                            def (meta, vcf2) = tuple
                                                            [meta.patient, meta, vcf2] }

    COERCE_SEQNAMES_HETS(het_pileups_wgs_jabba)
    chr_coerced_het_pileups_wgs_jabba = COERCE_SEQNAMES_HETS.out.file
    chr_coerced_het_pileups_wgs_jabba_to_cross = chr_coerced_het_pileups_wgs_jabba.map { tuple ->
                                                                       def (meta, hets) = tuple
                                                                       [meta.patient, meta, hets] }

    input_jab = chr_coerced_cov_rds_jabba_to_cross.join(chr_coerced_het_pileups_wgs_jabba_to_cross)
                                                  .join(chr_coerced_junction_jabba_to_cross)
                                                  .join(chr_coerced_j_supp_jabba_to_cross)
                                                  .map{ tuples ->
                                                            [tuples[1]] + [tuples[2]] + [tuples[4]] + [tuples[6]] + [tuples[8]]
                                                        }
    input_coerced_cov = input_jab.map{ meta, cov, hets, vcf, vcf2 -> [ meta, cov ] }     //chr stripped cov
    input_coerced_hets = input_jab.map{ meta, cov, hets, vcf, vcf2 -> [ meta, hets ] }   //chr stripped hetpileups
    input_coerced_vcf = input_jab.map{ meta, cov, hets, vcf, vcf2 -> [ meta, vcf ] }     //chr stripped somatic sv
    input_coerced_vcf2 = input_jab.map{ meta, cov, hets, vcf, vcf2 -> [ meta, vcf2 ] }   //chr stripped unfiltered somatic sv

    JABBA(input_coerced_cov, input_coerced_vcf, ploidy_jabba, input_coerced_hets,
    cbs_seg_rds_jabba, cbs_nseg_rds_jabba, input_coerced_vcf2, blacklist_junctions_jabba,
    geno_jabba, indel_jabba, tfield_jabba,
    iter_jabba, rescue_window_jabba, rescue_all_jabba, nudgebalanced_jabba,
    edgenudge_jabba, strict_jabba, allin_jabba, field_jabba, maxna_jabba,
    blacklist_coverage_jabba, purity_jabba, pp_method_jabba, cnsignif_jabba,
    slack_jabba, linear_jabba, tilim_jabba, epgap_jabba, fix_thres_jabba, lp_jabba,
    ism_jabba, filter_loose_jabba, gurobi_jabba, verbose_jabba)

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
