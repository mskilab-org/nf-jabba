process DRYCLEAN {
    tag "$meta.id"
    label 'process_medium'
    publishDir "results", mode: 'copy'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/dryclean:latest':
        'mskilab/dryclean:latest' }"

    input:
    tuple val(meta), path(pon) path(input)
    val centered
    val cbs
    val cnsignif
    val cores
    val wholeGenome
    val blacklist
    path blacklist_path
    val germline_filter
    path germline_file
    val human
    val field
    val build

    output:
    tuple val(meta), path("*.drycleaned.cov.rds")     , emit: decomposed_cov, optional: true
    tuple val(meta), path("*.dryclean.object.rds")    , emit: dryclean_object, optional: true
    //path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
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
    if [ -d "$R_LIB_PATH" ]; then
        export R_LIBS=$R_LIB_PATH
    fi

    echo "Using R version:"
    cmd="R --version"
    eval $cmd
    set +x

    ## find R installation and dryclean exec
    echo "USING LIBRARIES: \$(Rscript -e 'print(.libPaths())')"
    export drycleanPath=\$(Rscript -e 'cat(suppressWarnings(find.package("dryclean")))')
    export drycln=$drycleanPath/extdata/drcln
    echo $drycln
    set +x
    CMD="Rscript $drycln \$@"

    if [ ! -s ./drycleaned.cov.rds ]; then
	if ! { echo "Running:" && echo "${CMD}" && eval ${CMD}; }; then
	    echo "Dryclean broke!"; exit 1;
	fi
    else
	echo "If you wish to rerun Dryclean - please purge directory first"
    fi
    exit 0
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch drycleaned.cov.rds
    """
}
