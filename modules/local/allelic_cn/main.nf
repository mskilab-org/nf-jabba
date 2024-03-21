process NON_INTEGER_BALANCE {

    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/allelic_cn:latest':
        'mskilab/allelic_cn:latest' }"

    input:
    tuple val(meta), path(jabba_rds), path(decomposed_cov), path(het_pileups_wgs)
    val(field)
    val(hets_thresh)
    path(mask)
    val(overwrite)
    val(lambda)
    val(allin)
    val(fix_thresh)
    val(nodebounds)
    val(ism)
    val(build)
    val(epgap)
    val(tilim)
    val(gurobi)
    path(fasta)     // path to decoy fasta
    val(pad)

    output:
    tuple val(meta), path("balanced.gg.rds")                , emit: non_integer_balance_balanced_gg, optional: true
    tuple val(meta), path("hets.gg.rds")                    , emit: non_integer_balance_hets_gg, optional: true
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def id          = "${meta.sample}"
    def VERSION    = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """

    export RSCRIPT_PATH=\$(echo "${baseDir}/bin/non_integer_balance.R")

    Rscript \$RSCRIPT_PATH \\
        --id $id \\
        --jab $jabba_rds \\
        --cov $decomposed_cov \\
        --field $field \\
        --hets $het_pileups_wgs \\
        --hets-thresh $hets_thresh \\
        --mask $mask \\
        --overwrite $overwrite \\
        --lambda $lambda \\
        --allin $allin \\
        --fix_thresh $fix_thres \\
        --nodebounds $nodebounds \\
        --ism $ism \\
        --build $build \\
        --epgap $epgap \\
        --tilim $tilim \\
        --gurobi $gurobi \\
        --fasta $fasta \\
        --pad $pad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        non_integer_balance: ${VERSION}
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch balanced.gg.rds hets.gg.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        non_integer_balance: ${VERSION}
    END_VERSIONS
    """
}

process LP_PHASED_BALANCE {

    tag "$id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/allelic_cn:latest':
        'mskilab/allelic_cn:latest' }"

    input:
    tuple val(meta), path(hets_gg), path(hets) // output from non_integer_balance, sites.txt from hetpileups
    val(lambda)
    val(cnloh)
    val(major)
    val(allin)
    val(marginal)
    val(from_maf)
    path(mask)
    val(ism)
    val(epgap)
    val(hets_thresh)
    val(min_bins)
    val(min_width)
    val(trelim)
    val(reward)
    val(nodefileind)
    val(tilim)

    output:
    tuple val(meta), path("balanced.gg.rds")                , emit: lp_phased_balance_balanced_gg, optional: true
    tuple val(meta), path("binstats.gg.rds")                , emit: lp_phased_balance_binstats_gg, optional: true
    tuple val(meta), path("unphased.gg.rds")                , emit: lp_phased_balance_unphased_allelic_gg, optional: true
    path "versions.yml"                                     , emit: versions

    script:
    def args = task.ext.args ?: ''
    def id   = "${meta.sample}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    Rscript \${baseDir}/bin/lp_phased_balance.R \\
        --id $id \\
        --jab $hets_gg \\
        --hets $hets \\
        --lambda $lambda \\
        --cnloh $cnloh \\
        --major $major \\
        --allin $allin \\
        --marginal $marginal \\
        --from_maf $from_maf \\
        --mask $mask \\
        --ism $ism \\
        --epgap $epgap \\
        --hets_thresh $hets_thresh \\
        --min_bins $min_bins \\
        --min_width $min_width \\
        --trelim $trelim \\
        --reward $reward \\
        --nodefileind $nodefileind \\
        --tilim $tilim

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        non_integer_balance: ${VERSION}
    END_VERSIONS
    """
}
