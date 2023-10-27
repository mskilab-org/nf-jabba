process CBS {

    tag "$meta.id"
    label 'process_medium'

    // TODO add fragcounter container
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/cbs:latest':
        'mskilab/cbs:latest' }"

    input:
    tuple val(meta), path(tumor_dryclean_cov), path(normal_dryclean_cov)
    val cnsignif
    val field
    val name

    output:
    path "*.cov.rds", emit: cbs_cov_rds
    path "*seg.rds", emit: cbs_seg_rds
    path "*nseg.rds", emit: cbs_nseg_rds

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def VERSION    = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    RSCRIPT_PATH=\$(if [[ ${workflow.containerEngine} == "singularity" && !task.ext.singularity_pull_docker_container ]]; then echo "/cbsFH.R"; else echo "\${baseDir}/bin/cbsFH.R"; fi)
    Rscript $RSCRIPT_PATH \\
        -t ${tumor_dryclean_cov} \
        -n ${normal_dryclean_cov} \
        --cnsignif ${cnsignif} \
        -m ${name} \
        -f ${field}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cbs: ${VERSION}
    END_VERSIONS

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch cov.rds
    touch seg.rds
    touch nseg.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cbs: ${VERSION}
    END_VERSIONS
    """

}
