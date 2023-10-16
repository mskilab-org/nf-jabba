process DRYCLEAN {
    tag "$meta.id"
    label 'process_medium'
    publishDir "results", mode: 'copy'

    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //    'biocontainers/YOUR-TOOL-HERE' }"

    input:
    path tumor_fragcounter_cov from InputCov
    path dryclean_pon from PON
    val pair from name
    val wholeGenome from WholeGenome
    val chromosome from Chromosome
    val cores from cores
    val germline_filter from germlineFilter
    val blacklist from blacklist
    val germline_file from germlineFile
    val collapse from Collapse
    val field from cov_field
    tuple val(meta), path(bam)

    output:
    path "*drycleaned.cov.rds", into: tumor_dryclean_cov
    //path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -o allexport
    . ~/.bash_profile
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
    echo "USING LIBRARIES: $(Rscript -e 'print(.libPaths())')"
    export drycleanPath=$(Rscript -e 'cat(suppressWarnings(find.package("dryclean")))')
    export drycln=$drycleanPath/extdata/drcln
    echo $drycln
    set +x
    CMD="Rscript $drycln $@"

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
