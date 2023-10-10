
process FRAGCOUNTER {
    tag "$meta.id"
    label 'process_medium'

    // TODO add fragcounter container
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(bam), path(bai)                    // Mandatory: Format should be [meta, bam, bai] : can also provide cram & crai
    path(cov)                                                // Path to existing coverage rds or bedgraph : Default is NULL
    val(midpoint)                                            // If TRUE only count midpoint if FALSE then count bin footprint of every fragment interval : Default is TRUE
    val(windowsize)                                          // Bin size: Default is 200
    path(gcmapdir)                                           // Directory containing GC & Mappability bias .rds files
    val(minmapq)                                             // Minimal map quality: Default is 1
    path(fasta)                                              // Mandatory if supplying cram file
    val(paired)                                              // If TRUE, will consider the data as paired: Default is TRUE
    val(exome)                                               // Boolean to mention if it is exome; Default is FALSE


    output:
    tuple val(meta), path("*.cov.rds")                      , emit: fragcounter_cov, optional: true
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def cov        = cov ? "-c ${cov}" : ""
    def midpoint   = midpoint ? "-m ${midpoint}" : ""
    def windowsize = windowsize ? "-w ${windowsize}" : ""
    def gcmapdir   = gcmapdir ? "-d ${gcmapdir}" : ""
    def minmapq    = minmapq ? "-q ${minmapq}" : ""
    def fasta      = fasta ? "-r ${fasta}" : ""
    def paired     = paired ? "-p ${paired}" : ""
    def exome      = exome ? "-e ${exome}" : ""


    """
    . ~/.bash_profile

    set -x
    R_LIB_PATH="/gpfs/commons/groups/imielinski_lab/lib/R-4.0.2"
    if [ -d "$R_LIB_PATH" ]; then
        export R_LIBS=${R_CUSTOM_LIBS}:$R_LIB_PATH
    fi

    export R_DATATABLE_NUM_THREADS=1

    ## find R installation
    echo "USING LIBRARIES: $(Rscript -e 'print(.libPaths())')"
    export PATH=${PATH}:$(Rscript -e 'cat(paste0(installed.packages()["fragCounter", "LibPath"], "/fragCounter/extdata/"))')
    set +x

    frag -b $bam \\
    $cov \\
    $windowsize \\
    $fasta \\
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
    touch ${meta.id}.cov.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fragcounter: ${VERSION}
    END_VERSIONS
    """    

}