process ASCAT_SEG {

    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/ascat_seg:latest':
        'mskilab/ascat_seg:latest' }"

    input:
    tuple val(meta), path(hets), path(cbs_cov)                       // channel: [mandatory] [ meta, hets ]
    val(field)                                                       // channel: [mandatory] "foreground" for dryclean/ "ratio"
    val(hets_thresh)                                                 // channel: cutoff for hetpileups; default=0.2
    val(penalty)                                                     // channel: penalty for ASCAT; default=70
    val(gc)                                                          // channel: perform GC correction? Default=TRUE
    val(rebin)                                                       // channel: width for rebinning, default=5e4
    val(from_maf)                                                    // channel: whether to start from MAF, default=FALSE

    output:
    tuple val(meta), path("*ascat_pp.rds")        , emit: purityploidy
    tuple val(meta), path("*ascat_seg.rds")       , emit: segments, optional:true
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.5.2'                              // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    ## RSCRIPT_PATH=\$(if [[ ${workflow.containerEngine} == "singularity" && !task.ext.singularity_pull_docker_container ]]; then echo "/usr/bin/ascat_seg.R"; else echo "\${baseDir}/bin/ascat_seg.R"; fi)
    export RSCRIPT_PATH=\$(echo "${baseDir}/bin/ascat_seg.R")

    Rscript \$RSCRIPT_PATH       \\
    --id          ${meta.id}     \\
    --variants    ${hets}        \\
    --coverage    ${cbs_cov}     \\
    --field       ${field}       \\
    --hets_thresh ${hets_thresh} \\
    --penalty     ${penalty}     \\
    --gc          ${gc}          \\
    --rebin_width ${rebin}       \\
    --from_maf    ${from_maf};

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ASCAT: ${VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.5.2'                              // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ascat_pp.rds
    touch ascat_seg.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ASCAT: ${VERSION}
    END_VERSIONS

    """
}

process EXTRACT_PURITYPLOIDY {
    input:
    tuple val(meta), path(ascat_rds)

    output:
    env purity_val, emit: purity_val
    tuple val(meta), env(ploidy_val), emit: ploidy_val

    script:
    """
	export purity_val=\$(Rscript -e "ascat <- readRDS('${ascat_rds}'); cat(ascat[['purity']])")
	export ploidy_val=\$(Rscript -e "ascat <- readRDS('${ascat_rds}'); cat(ascat[['ploidy']])")
    """
}
