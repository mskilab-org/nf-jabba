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
    inputs  // [ meta, junction, cov, j_supp, hets, ploidy, seg_cbs, nseg_cbs ]
    blacklist_junctions
    geno
    indel
    tfield
    iter
    rescue_window
    rescue_all
    nudgebalanced
    edgenudge
    strict
    allin
    field
    maxna
    blacklist_coverage
    purity
    pp_method
    cnsignif
    slack
    linear
    tilim
    epgap
    fix_thres
    lp
    ism
    filter_loose
    gurobi
    verbose

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
    chr_coerced_cov_rds          = Channel.empty()
    chr_coerced_junction         = Channel.empty()
    chr_coerced_j_supp           = Channel.empty()
    chr_coerced_het_pileups_wgs  = Channel.empty()

    // put inputs into a map for easier retrieval
    inputs_map = inputs.map { samples ->
        [
            "meta": samples[0],
            "junction": samples[1],
            "cov": samples[2],
            "j_supp": samples[3],
            "hets": samples[4],
            "ploidy": samples[5],
            "seg_cbs": samples[6],
            "nseg_cbs": samples[7]
        ]
    }

    // Run COERCE_SEQNAMES to force inputs to have the same sequence name format
    junction = inputs_map.map { sample ->
        [sample.meta, sample.junction]
    }
    COERCE_SEQNAMES_SOM_SV(junction)
    chr_coerced_junction = COERCE_SEQNAMES_SOM_SV.out.file
    chr_coerced_junction = chr_coerced_junction.map { meta, som_sv ->
        [meta.patient, som_sv]
    }

    cov_rds = inputs_map.map { sample ->
        [sample.meta, sample.cov]
    }
    COERCE_SEQNAMES_COV(cov_rds)
    chr_coerced_cov_rds = COERCE_SEQNAMES_COV.out.file
    chr_coerced_cov_rds = chr_coerced_cov_rds.map { meta, cov ->
        [meta.patient, cov]
    }

    j_supp = inputs_map.map { sample ->
        [sample.meta, sample.j_supp]
    }
    COERCE_SEQNAMES_UNFIL_SOM_SV(j_supp)
    chr_coerced_j_supp = COERCE_SEQNAMES_UNFIL_SOM_SV.out.file
    chr_coerced_j_supp = chr_coerced_j_supp.map { meta, j_supp ->
        [meta.patient, j_supp]
    }

    het_pileups_wgs = inputs_map.map { sample ->
        [sample.meta, sample.hets]
    }
    COERCE_SEQNAMES_HETS(het_pileups_wgs)
    chr_coerced_het_pileups_wgs = COERCE_SEQNAMES_HETS.out.file
    chr_coerced_het_pileups_wgs = chr_coerced_het_pileups_wgs.map { meta, hets ->
        [meta.patient, hets]
    }


    // join coerced inputs into a single channel on key: meta.patient
    inputs = inputs.map { sample ->
        [sample[0].patient] + sample
    }
    final_inputs = inputs
                    .join(chr_coerced_junction)
                    .join(chr_coerced_cov_rds)
                    .join(chr_coerced_j_supp)
                    .join(chr_coerced_het_pileups_wgs)
                    .map { sample ->
                        [
                            sample[1],  // meta (because sample[0] is meta.patient)
                            sample[9],  // chr_coerced_junction
                            sample[10], // chr_coerced_cov_rds
                            sample[11],  // chr_coerced_j_supp
                            sample[12], // chr_coerced_het_pileups_wgs
                            sample[6],  // ploidy
                            sample[7],  // seg_cbs
                            sample[8]  // nseg_cbs
                        ]
                    }

    JABBA(
        final_inputs,
        blacklist_junctions,
        geno,
        indel,
        tfield,
        iter,
        rescue_window,
        rescue_all,
        nudgebalanced,
        edgenudge,
        strict,
        allin,
        field,
        maxna,
        blacklist_coverage,
        purity,
        pp_method,
        cnsignif,
        slack,
        linear,
        tilim,
        epgap,
        fix_thres,
        lp,
        ism,
        filter_loose,
        gurobi,
        verbose
    )

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
