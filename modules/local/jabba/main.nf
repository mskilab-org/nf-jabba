process JABBA {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskilab/jabba_cplex:latest':
        'mskilab/jabba_cplex:latest' }"

    input:
    tuple val(meta), path(cov_rds)
    tuple val(meta), path(junction)
    tuple val(meta), val(ploidy)
    tuple val(meta), path(het_pileups_wgs)
    tuple val(meta), path(cbs_seg_rds, stageAs: "cbs_seg.rds")
    tuple val(meta), path(cbs_nseg_rds, stageAs: "cbs_nseg.rds")
    tuple val(meta), path(j_supp)
    val(blacklist_junctions)
    val(geno)
    val(indel)
    val(tfield)
    val(iter)
    val(rescue_window)
    val(rescue_all)
    val(nudgebalanced)
    val(edgenudge)
    val(strict)
    val(allin)
    val(field)
    val(maxna)
    path(blacklist_coverage)
    val(purity)
    val(pp_method)
    val(cnsignif)
    val(slack)
    val(linear)
    val(tilim)
    val(epgap)
    val(name)
    val(fix_thres)
    val(lp)
    val(ism)
    val(filter_loose)
    val(gurobi)
    val(nonintegral)
    val(verbose)
    val(help)

    output:
    tuple val(meta), path("*.jabba.simple.rds")      , emit: jabba_rds, optional: true
    tuple val(meta), path("*.jabba.simple.gg.rds")   , emit: jabba_gg, optional: true
    tuple val(meta), path("*.jabba.simple.vcf")      , emit: jabba_vcf, optional: true
    tuple val(meta), path("*.jabba.raw.rds")         , emit: jabba_raw_rds, optional: true
    tuple val(meta), path("*.opt.report.rds")        , emit: opti, optional: true
    tuple val(meta), path("*.jabba.seg")             , emit: jabba_seg, optional: true
    tuple val(meta), path("*.karyograph.rds")        , emit: karyograph, optional: true
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def VERSION    = '1.1'
    """
    #!/bin/bash

    set -o allexport



    ## find R installation
    echo "USING LIBRARIES: \$(Rscript -e 'print(.libPaths())')"

    export jabPath=\$(Rscript -e 'cat(suppressWarnings(find.package("JaBbA")))')
    export jba=\${jabPath}/extdata/jba
    echo \$jba
    set +x

    export cmd="Rscript \$jba ${junction} ${cov_rds} \\
    --j.supp                ${j_supp} \\
    --indel					${indel} \\
    --tfield				${tfield} \\
    --field					${field} \\
    --seg			        ${cbs_seg_rds} \\
    --blacklist.coverage	${blacklist_coverage} \\
    --nseg			        ${cbs_nseg_rds} \\
    --hets		            ${het_pileups_wgs} \\
    --ploidy				${ploidy} \\
    --purity				${purity} \\
    --ppmethod				${pp_method} \\
    --cores                 ${task.cpus} \\
    "

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        JaBbA: ${VERSION}
    END_VERSIONS

    { echo "Running:" && echo "\$(echo \$cmd)" && echo && eval \$cmd; }
    cmdsig=\$?
    if [ "\$cmdsig" = 0 ]; then
        echo "Finish!"
    else
        echo "Broke!"
        exit \$cmdsig
    fi
    
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    touch jabba.simple.rds
    touch jabba.simple.gg.rds
    touch jabba.simple.vcf
    touch jabba.raw.rds
    touch opt.report.rds
    touch jabba.seg
    touch karyograph.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        JaBbA: ${VERSION}
    END_VERSIONS
    """
}
