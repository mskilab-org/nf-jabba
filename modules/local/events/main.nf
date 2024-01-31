process EVENTS {

    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/events:latest':
        'mskilab/events:latest' }"

    input:
    tuple val(meta), path(gGraph), path(ref)
    val(libdir)
    val(id)
    val(outdir)

    output:
    tuple val(meta), path("${outdir}/*") , emit: events_output
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def VERSION    = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """

    export RSCRIPT_PATH=\$(echo "${baseDir}/bin/Events.R")

    Rscript \$RSCRIPT_PATH \\
	--libdir $libdir \\
	--id $id \\
	--gGraph $gGraph \\
	--ref $ref \\
	--outdir $outdir

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Events: ${VERSION}
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    mkdir -p ${outdir}
    touch ${outdir}/dummy_output.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        events: ${VERSION}
    END_VERSIONS
    """
}
