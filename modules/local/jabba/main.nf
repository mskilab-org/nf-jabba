process JABBA {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/jabba:latest':
        'mskilab/jabba:latest' }"

    input:
    tuple val(meta), path(cov_rds), path(junctionFilePath)
    val field
    path junctionUnfiltered
    val tfield
    path cbs_nseg_rds
    path cbs_seg_rds
    val slack
    path het_pileups_wgs
    val purity
    val ploidy
    val tilim
    val epgap
    val pp_method
    val maxna
    val flags
    path blacklist_coverage
    path blacklist_junctions
    val iter
    val pair
    val indel
    val cnsignif
    val lp
    val ism
    val treemem
    val fix_thres
    val gurobi
    val nonintegral

    output:
    tuple val(meta), path("*.jabba.simple.rds")      , emit: jabba_rds, optional: true
    tuple val(meta), path("*.jabba.simple.gg.rds")   , emit: jabba_gg, optional: true
    tuple val(meta), path("*.jabba.simple.vcf")      , emit: jabba_vcf, optional: true
    tuple val(meta), path("*.jabba.raw.rds")         , emit: jabba_raw_rds, optional: true
    tuple val(meta), path("*.opt.report.rds")        , emit: opti, optional: true
    tuple val(meta), path("*.jabba.seg")             , emit: jabba_seg, optional: true
    tuple val(meta), path("*.karyograph.rds")        , emit: karyograph, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    """
    set -o allexport

    # Check if the environment has the module program installed
    if command -v module &> /dev/null
    then
        # Check if the required modules are available
        if module avail R/4.0.2 &> /dev/null && module avail gcc/9.2.0 &> /dev/null
        then
            ## load correct R and gcc versions
            module unload R
            module load R/4.0.2
            module unload gcc
            module load gcc/9.2.0
        fi
    fi

    set -x
    R_LIB_PATH="/gpfs/commons/groups/imielinski_lab/lib/R-4.0.2"
    if [ -d "$R_LIB_PATH" ]; then
        export R_LIBS=${R_CUSTOM_LIBS}:$R_LIB_PATH
    fi
    export R_DATATABLE_NUM_THREADS=1
    unset R_HOME
    R_PROFILE_USER="/dev/null"

    ## find R installation
    echo "USING LIBRARIES: $(Rscript -e 'print(.libPaths())')"

    export jabPath=$(Rscript -e 'cat(suppressWarnings(find.package("JaBbA")))')
    export jba=$jabPath/extdata/jba
    echo $jba
    set +x

    export cmd="Rscript $jba $@"

    { echo "Running:" && echo "$(echo $cmd)" && echo && eval $cmd; }
    cmdsig=$?
    if [ "$cmdsig" = 0 ]; then
        echo "Finish!"
    else
        echo "Broke!"
        exit $cmdsig
    fi

    exit 0
    """

    stub:
    """
    touch jabba.simple.rds
    touch jabba.simple.gg.rds
    touch jabba.simple.vcf
    touch jabba.raw.rds
    touch opt.report.rds
    touch jabba.seg
    touch karyograph.rds
    """
}
