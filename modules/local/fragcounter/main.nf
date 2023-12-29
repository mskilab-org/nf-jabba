process FRAGCOUNTER {

    tag "$meta.id"
    label 'process_medium'

    // TODO add fragcounter container
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/fragcounter:0.1':
        'mskilab/fragcounter:0.1' }"

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

process REBIN_RAW_FRAGCOUNTER {

    tag "$meta.id"
    label 'process_low'

    // TODO add fragcounter container
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/fragcounter:0.1':
        'mskilab/fragcounter:0.1' }"

    input:
    tuple val(meta), path(cov_raw)
    val(field)
    val(windowsize)

    output:
    tuple val(meta), path("1kb_*"), emit: raw_fragcounter_cov_1kb, optional:true

    script:

    """
    #!/usr/bin/env Rscript

    library(skitools)

    filename = "${cov_raw}"
    outputfn = "1kb_${cov_raw.name}"

    raw_cov = readRDS(filename)
    collapse.cov <- function(cov.gr, bin.size = 1e3, field = "reads.corrected") {
        BINSIZE.ROUGH = bin.size
        cov.gr = cov.gr[, field]
        cov.gr = gr2dt(cov.gr)
        setnames(cov.gr, field, "signal")
        cov.gr = cov.gr[!is.infinite(signal), .(signal = median(signal, na.rm = TRUE)),
                        by = .(seqnames, start = floor(start/BINSIZE.ROUGH)*BINSIZE.ROUGH+1)]
        cov.gr[, end := (start + BINSIZE.ROUGH) - 1]
        setnames(cov.gr, "signal", field)
        cov.gr = dt2gr(cov.gr)
        return(cov.gr)
    }
    rebinned_cov = collapse.cov(raw_cov, bin.size = ${windowsize}, field = "${field}")
    rebinned_cov = rebinned_cov %Q% (!seqnames=="Y")
    saveRDS(rebinned_cov, outputfn)

    """

}
