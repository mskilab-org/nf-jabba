process ASCAT_SEG {

    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/ascat_seg:latest':
        'mskilab/ascat_seg:latest' }"

    input:
    tuple val(meta), path(hets)                                      // channel: [mandatory] [ meta, hets ]
    tuple val(meta), path(cbs_cov)                                   // channel: [mandatory] [ meta, cbs_cov ]
    val(field)                                                       // channel: [mandatory] "foreground" for dryclean/ "ratio" 
    val(hets_thresh)                                                 // channel: cutoff for hetpileups; default=0.2
    val(penalty)                                                     // channel: penalty for ASCAT; default=70
    val(gc)                                                          // channel: perform GC correction? Default=TRUE
    val(rebin)                                                       // channel: width for rebinning, default=5e4
    val(from_maf)                                                    // channel: whether to start from MAF, default=FALSE

    output:
    tuple val(meta), path("*ascat_pp.rds")        , emit: purityploidy, optional:true
    tuple val(meta), path("*ascat_seg.rds")       , emit: segments, optional:true
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.5.2'                              // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    #!/bin/bash
    set -o allexport
    # Check if the environment has the module program installed
    if command -v module &> /dev/null
    then
        # Check if the required modules are available
        if module avail R/4.0.2 &> /dev/null
        then
            ## load correct R and gcc versions
            module unload R
            module load R/4.0.2
        fi
    fi

    set -x
    R_LIB_PATH="~/lab/lib/R-4.0.2"
    if [ -d "\$R_LIB_PATH" ]; then
        export R_LIBS=\$R_LIB_PATH
    fi

    echo "Using R version:"
    cmd="R --version"
    eval \$cmd
    set +x

    ## find R installation and dryclean exec
    echo "USING LIBRARIES: \$(Rscript -e 'print(.libPaths())')"

    RSCRIPT_PATH=\$(if [[ ${workflow.containerEngine} == "singularity" && !task.ext.singularity_pull_docker_container ]]; then echo "./ascat_seg.R"; else echo "\${baseDir}/bin/ascat_seg.R"; fi)
    
    Rscript \$RSCRIPT_PATH       \\
    --id          ${meta.id}     \\
    --variants    ${hets}        \\
    --coverage    ${cbs_cov}     \\
    --field       ${field}       \\
    --hets_thresh ${hets_thresh} \\
    --penalty     ${penalty}     \\
    --gc          ${gc}          \\
    --rebin_width ${rebin}       \\
    --from_maf    ${from_maf}

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
    touch sites.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ASCAT: ${VERSION}
    END_VERSIONS

    """

}