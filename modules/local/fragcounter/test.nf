process FRAGCOUNTER {

    tag "$meta.id"
    label 'process_medium'

    //conda "bioconda::gridss=2.13.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'tanubrata2/fragcounter:0.1':
        'tanubrata2/fragcounter:0.1' }"

    input:
    tuple val(metaid), val(meta), path(bam), path(bai)
    val(windowsize)
    path(gcmapdir)
    val(minmapq)
    path(fasta)
    path(fasta_fai)

    output:
    tuple val(meta), path("*cov.rds")                       , emit: fragcounter_cov, optional: true
    tuple val(meta), path("*cov.raw.rds")                   , emit: fragcounter_raw_cov, optional: true
    tuple val(meta), path("*cov.original.bw")               , emit: raw_bw, optional: true
    tuple val(meta), path("*cov.corrected.bw")              , emit: corrected_bw, optional: true
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def windowsize = windowsize ? "-w ${windowsize}" : ""
    def gcmapdir   = gcmapdir ? "-d ${gcmapdir}" : ""
    def minmapq    = minmapq ? "-q ${minmapq}" : ""
    def fasta      = fasta ? "-r ${fasta}" : ""
    def fasta_fai  = fasta_fai ? "cp -s ${fasta_fai} ."
    def VERSION    = '0.1'

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

    ${fasta_fai}

    frag -b $bam \\
    $windowsize \\
    $fasta \\
    $gcmapdir \\
    $minmapq \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fragcounter: ${VERSION}
    END_VERSIONS


    """

