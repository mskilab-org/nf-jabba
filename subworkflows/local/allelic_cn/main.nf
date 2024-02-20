include { NON_INTEGER_BALANCE } from '../../../modules/local/allelic_cn/main.nf'

workflow NON_INTEGER_BALANCE {
    take:
	jabba_rds_non_integer_balance           // [meta, ggraph]
	decomposed_cov_non_integer_balance      // [meta, cov]
	het_pileups_wgs_non_integer_balance     // [meta, hets]
	id_non_integer_balance
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
	fasta_non_integer_balance    // path to decoy fasta
	pad_non_integer_balance

    main:
    versions            = Channel.empty()
    non_integer_balance_balanced_gg = Channel.empty()
    non_integer_balance_hets_gg = Channel.empty()

    NON_INTEGER_BALANCE(
        jabba_rds_non_integer_balance,
        decomposed_cov_non_integer_balance,
        het_pileups_wgs_non_integer_balance,
        id_non_integer_balance,
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
        pad_non_integer_balance,
    )

    non_integer_balance_balanced_gg = NON_INTEGER_BALANCE.out.non_integer_balance_balanced_gg
    non_integer_balance_hets_gg = NON_INTEGER_BALANCE.out.non_integer_balance_hets_gg
    versions = NON_INTEGER_BALANCE.out.versions

    emit:
    non_integer_balance_balanced_gg
    non_integer_balance_hets_gg

    versions
}
