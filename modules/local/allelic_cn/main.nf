process NON_INTEGER_BALANCE {

    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/non_integer_balance:latest':
        'mskilab/non_integer_balance:latest' }"

    input:
    tuple val(meta), path(jabba_rds)
    tuple val(meta), path(decomposed_cov)
    tuple val(meta), path(het_pileups_wgs)
    val(id)
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

process LP_PHASE {

    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/allelic_cn:latest':
        'mskilab/allelic_cn:latest' }"

    input:
    tuple val(meta), path(tumor_bam_wgs, stageAs: "tumor.bam"), path(tumor_bai, stageAs: "tumor.bam.bai"), path(normal_bam_wgs, stageAs: "normal.bam"), path(normal_bai, stageAs: "normal.bam.bai")
    val(filter)
    val(max_depth)
    path(hapmap_sites)

    output:
    tuple val(meta), path("*sites.txt")                     , emit: het_pileups_wgs, optional: true
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def VERSION    = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """

    ## export RSCRIPT_PATH=\$(if [[ ${workflow.containerEngine} == "singularity" && !task.ext.singularity_pull_docker_container ]]; then echo "/Pileup.R"; else echo "\${baseDir}/bin/Pileup.R"; fi)

    ##export RSCRIPT_PATH=\$(if [[ ${workflow.containerEngine} == "singularity" && !task.ext.singularity_pull_docker_container ]]; then
    ##                                echo "/Pileup.R"
    ##                            else
    ##                                echo "${baseDir}/bin/Pileup.R"
    ##                        fi)

    export RSCRIPT_PATH=\$(echo "${baseDir}/bin/Pileups.R")

    Rscript \$RSCRIPT_PATH \\
	--tbam $tumor_bam_wgs \\
	--filter $filter \\
	--nbam $normal_bam_wgs \\
	--ncores $task.cpus \\
	--maxdepth $max_depth \\
	--sites $hapmap_sites

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        allelic_cn: ${VERSION}
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch sites.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        allelic_cn: ${VERSION}
    END_VERSIONS
    """
}
