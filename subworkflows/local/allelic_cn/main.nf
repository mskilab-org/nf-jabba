include { NON_INTEGER_BALANCE } from '../../../modules/local/allelic_cn/main.nf'
include { LP_PHASED_BALANCE } from '../../../modules/local/allelic_cn/main.nf'

workflow COV_GGRAPH_NON_INTEGER_BALANCE {
    take:
    inputs_non_integer_balance          // [ meta, ggraph, cov, hets ]
	field_non_integer_balance
	hets_thresh_non_integer_balance
	mask_non_integer_balance
	overwrite_non_integer_balance
	lambda_non_integer_balance
	allin_non_integer_balance
	fix_thresh_non_integer_balance
	nodebounds_non_integer_balance
	ism_non_integer_balance
	build_non_integer_balance
	epgap_non_integer_balance
	tilim_non_integer_balance
	gurobi_non_integer_balance
	fasta_non_integer_balance           // path to decoy fasta
	fasta_fai_non_integer_balance       // path to decoy fasta.fai
    bwa_non_integer_balance
	pad_non_integer_balance

    main:
    versions            = Channel.empty()
    non_integer_balance_balanced_gg = Channel.empty()
    non_integer_balance_hets_gg = Channel.empty()

    NON_INTEGER_BALANCE(
        inputs_non_integer_balance,
        field_non_integer_balance,
        hets_thresh_non_integer_balance,
        mask_non_integer_balance,
        overwrite_non_integer_balance,
        lambda_non_integer_balance,
        allin_non_integer_balance,
        fix_thresh_non_integer_balance,
        nodebounds_non_integer_balance,
        ism_non_integer_balance,
        build_non_integer_balance,
        epgap_non_integer_balance,
        tilim_non_integer_balance,
        gurobi_non_integer_balance,
        fasta_non_integer_balance,
        fasta_fai_non_integer_balance,
        bwa_non_integer_balance,
        pad_non_integer_balance
    )

    non_integer_balance_balanced_gg = NON_INTEGER_BALANCE.out.non_integer_balance_balanced_gg
    non_integer_balance_hets_gg = NON_INTEGER_BALANCE.out.non_integer_balance_hets_gg
    versions = NON_INTEGER_BALANCE.out.versions

    emit:
    non_integer_balance_balanced_gg
    non_integer_balance_hets_gg

    versions
}

workflow COV_GGRAPH_LP_PHASED_BALANCE {
    take:
    inputs_lp_phased_balance        // [ meta, ggraph, hets ]
    lambda_lp_phased_balance
    cnloh_lp_phased_balance
    major_lp_phased_balance
    allin_lp_phased_balance
    marginal_lp_phased_balance
    from_maf_lp_phased_balance
    mask_lp_phased_balance
    ism_lp_phased_balance
    epgap_lp_phased_balance
    hets_thresh_lp_phased_balance
    min_bins_lp_phased_balance
    min_width_lp_phased_balance
    trelim_lp_phased_balance
    reward_lp_phased_balance
    nodefileind_lp_phased_balance
    tilim_lp_phased_balance

    main:
    versions            = Channel.empty()
    lp_phased_balance_balanced_gg = Channel.empty()
    lp_phased_balance_binstats_gg = Channel.empty()
    lp_phased_balance_unphased_allelic_gg = Channel.empty()

    LP_PHASED_BALANCE(
        inputs_lp_phased_balance,
        lambda_lp_phased_balance,
        cnloh_lp_phased_balance,
        major_lp_phased_balance,
        allin_lp_phased_balance,
        marginal_lp_phased_balance,
        from_maf_lp_phased_balance,
        mask_lp_phased_balance,
        ism_lp_phased_balance,
        epgap_lp_phased_balance,
        hets_thresh_lp_phased_balance,
        min_bins_lp_phased_balance,
        min_width_lp_phased_balance,
        trelim_lp_phased_balance,
        reward_lp_phased_balance,
        nodefileind_lp_phased_balance,
        tilim_lp_phased_balance
    )

    lp_phased_balance_balanced_gg = LP_PHASED_BALANCE.out.lp_phased_balance_balanced_gg
    lp_phased_balance_binstats_gg = LP_PHASED_BALANCE.out.lp_phased_balance_binstats_gg
    lp_phased_balance_unphased_allelic_gg = LP_PHASED_BALANCE.out.lp_phased_balance_unphased_allelic_gg
    versions = LP_PHASED_BALANCE.out.versions

    emit:
    lp_phased_balance_balanced_gg
    lp_phased_balance_binstats_gg
    lp_phased_balance_unphased_allelic_gg

    versions
}
