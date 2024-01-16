process HETPILEUPS {

    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/hetpileups:latest':
        'mskilab/hetpileups:latest' }"

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
        HetPileups: ${VERSION}
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch sites.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hetpileups: ${VERSION}
    END_VERSIONS
    """
}
