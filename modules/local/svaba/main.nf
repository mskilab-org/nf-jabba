
process SVABA {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::svaba=1.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/svaba:1.1.0--h7d7f7ad_2':
        'biocontainers/svaba:1.1.0--h7d7f7ad_2' }"


    input:
    tuple val(meta), path(tumorbam), path(tumorbai), path(normalbam), path(normalbai)
    path fasta
    path fasta_fai
    path bwa_index
    path dbsnp
    path dbsnp_tbi
    path indel_mask
    path germ_sv_db
    path simple_seq_db
    val error_rate

    output:
    tuple val(meta), path("*.svaba.sv.vcf.gz")                        , emit: sv, optional: true
    tuple val(meta), path("*.svaba.indel.vcf.gz")                     , emit: indel, optional: true
    tuple val(meta), path("*.svaba.germline.indel.vcf.gz")            , emit: germ_indel, optional: true
    tuple val(meta), path("*.svaba.germline.sv.vcf.gz")               , emit: germ_sv, optional: true
    tuple val(meta), path("*.svaba.somatic.indel.vcf.gz")             , emit: som_indel,  optional: true
    tuple val(meta), path("*.svaba.somatic.sv.vcf.gz")                , emit: som_sv, optional: true
    tuple val(meta), path("*.svaba.unfiltered.sv.vcf.gz")             , emit: unfiltered_sv, optional: true
    tuple val(meta), path("*.svaba.unfiltered.indel.vcf.gz")          , emit: unfiltered_indel, optional: true
    tuple val(meta), path("*.svaba.unfiltered.germline.indel.vcf.gz") , emit: unfiltered_germ_indel, optional: true
    tuple val(meta), path("*.svaba.unfiltered.germline.sv.vcf.gz")    , emit: unfiltered_germ_sv, optional: true
    tuple val(meta), path("*.svaba.unfiltered.somatic.indel.vcf.gz")  , emit: unfiltered_som_indel,  optional: true
    tuple val(meta), path("*.svaba.unfiltered.somatic.sv.vcf.gz")     , emit: unfiltered_som_sv, optional: true
    tuple val(meta), path("*.bps.txt.gz")                             , emit: raw_calls
    tuple val(meta), path("*.alignments.txt.gz")                      , emit: ascii_alignments, optional:true
    tuple val(meta), path("*.discordants.txt.gz")                     , emit: discordants, optional: true
    tuple val(meta), path("*.log")                                    , emit: log
    path "versions.yml"                                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def bamlist    = normalbam ? "-t ${tumorbam} -n ${normalbam}" : "-t ${tumorbam}"
    def dbsnp      = dbsnp ? "--dbsnp-vcf ${dbsnp}" : ""
    def bwa        = bwa_index ? "cp -s ${bwa_index}/* ." : ""
    def indel_mask = indel_mask ? "--blacklist ${indel_mask}" : ""
    def flags      = germ_sv_db ? "--germline-sv-database ${germ_sv_db} --simple-seq-database ${simple_seq_db}" : ""
    def error_rate = error_rate ? "--error-rate ${error_rate}" : ""

    """
    ${bwa}

    svaba \\
        run \\
        $bamlist \\
        --threads $task.cpus \\
        $dbsnp \\
        $indel_mask \\
        $error_rate \\
        $flags \\
        --id-string $meta.id \\
        --reference-genome $fasta \\
        --g-zip \\
        $args \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svaba: \$(echo \$(svaba --version 2>&1) | sed 's/[^0-9.]*\\([0-9.]*\\).*/\\1/' )
    END_VERSIONS
    """
    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bps.txt.gz
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svaba: \$(echo \$(svaba --version 2>&1) | sed 's/[^0-9.]*\\([0-9.]*\\).*/\\1/' )
    END_VERSIONS
    """

}