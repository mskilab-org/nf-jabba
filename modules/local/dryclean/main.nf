process DRYCLEAN {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/dryclean:0.0.2':
        'mskilab/dryclean:0.0.2' }"


    input:
    tuple val(meta), path(input)
    path(pon)
    val(centered)
    val(cbs)
    val(cnsignif)
    val(wholeGenome)
    val(blacklist)
    val(blacklist_path)
    val(germline_filter)
    val(germline_file)
    val(human)
    val(field)
    val(build)

    output:
    tuple val(meta), path("*cov.rds")                 , emit: decomposed_cov
    //tuple val(meta), path("*.dryclean.object.rds")    , emit: dryclean_object, optional: true
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.0.2'
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
    export drycleanPath=\$(Rscript -e 'cat(suppressWarnings(find.package("dryclean")))')
    export drycln=\$drycleanPath/extdata/drcln
    echo \$drycln
    set +x

    CMD="Rscript \$drycln \\
        --input             ${input} \\
        --pon               ${pon} \\
        --centered          ${centered} \\
        --cbs               ${cbs} \\
        --cnsignif          ${cnsignif} \\
        --cores             ${task.cpus} \\
        --wholeGenome       ${wholeGenome} \\
        --blacklist         ${blacklist} \\
        --blacklist_path    ${blacklist_path} \\
        --germline.filter   ${germline_filter} \\
        --germline.file     ${germline_file} \\
        --human             ${human} \\
        --field             ${field} \\
        --build             ${build} \\
    "

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dryclean: ${VERSION}
    END_VERSIONS

    if [ ! -s ./drycleaned.cov.rds ]; then
	    if ! { echo "Running:" && echo "\${CMD}" && eval \${CMD}; }; then
	        echo "Dryclean broke!"; exit 1;
	    fi
    else
	    echo "If you wish to rerun Dryclean - please purge directory first"
    fi

    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch drycleaned.cov.rds
    """
}
