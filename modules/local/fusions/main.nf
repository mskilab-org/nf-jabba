process FUSIONS {

    tag "$meta.id"
    label 'process_medium'

    // using events container since the dependencies are the same
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/events:latest':
        'mskilab/events:latest' }"

    input:
    tuple val(meta), path(gGraph)
    path(gencode)

    output:
    tuple val(meta), path("*fusions.rds") , emit: fusions_output
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def id          = "${meta.sample}"
    def VERSION    = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """

    export RSCRIPT_PATH=\$(echo "${baseDir}/bin/fusions.R")

    Rscript \$RSCRIPT_PATH \\
	--id $id \\
	--gGraph $gGraph \\
	--gencode $gencode \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusions: ${VERSION}
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch dummy_output.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusions: ${VERSION}
    END_VERSIONS
    """
}
