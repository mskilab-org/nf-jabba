process FRAGCOUNTER {

    tag "$meta.id"
    label 'process_medium'

    // TODO add fragcounter container
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/fragcounter:latest':
        'mskilab/fragcounter:latest' }"

    input:
    tuple val(meta), path(bam), path(bai)                    // Mandatory: Format should be [meta, bam, bai] : can also provide cram & crai
    val(midpoint)                                            // If TRUE only count midpoint if FALSE then count bin footprint of every fragment interval : Default is TRUE
    val(windowsize)                                          // Bin size: Default is 200
    path(gcmapdir)                                           // Directory containing GC & Mappability bias .rds files
    val(minmapq)                                             // Minimal map quality: Default is 1
    path(fasta)                                              // Mandatory if supplying cram file
    path(fasta_fai)                                          // Fasta index .fai
    val(paired)                                              // If TRUE, will consider the data as paired: Default is TRUE
    val(exome)                                               // Boolean to mention if it is exome; Default is FALSE


    output:
    tuple val(meta), path("*cov.rds")                       , emit: fragcounter_cov, optional: true
    tuple val(meta), path("*cov.corrected.bw")              , emit: corrected_bw, optional: true
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    //def cov        = cov ? "-c ${cov}" : ""
    def midpoint    = midpoint ? "-m ${midpoint}" : ""
    def windowsize  = windowsize ? "-w ${windowsize}" : ""
    def gcmapdir    = gcmapdir ? "-d ${gcmapdir}" : ""
    def minmapq     = minmapq ? "-q ${minmapq}" : ""
    def fasta       = fasta ? "-r ${fasta}" : ""
    def paired      = paired ? "-p ${paired}" : ""
    def exome       = exome ? "-e ${exome}" : ""
    def VERSION    = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    export R_DATATABLE_NUM_THREADS=1

    frag -b $bam \\
    $fasta \\
    $windowsize \\
    $gcmapdir \\
    $minmapq \\
    $paired \\
    $exome

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fragcounter: ${VERSION}
    END_VERSIONS

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch cov.rds
    touch cov.corrected.bw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fragcounter: ${VERSION}
    END_VERSIONS
    """    

}