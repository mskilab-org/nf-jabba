workflow CHANNEL_COVERAGE_CREATE_CSV {
    take:
        n_cov // channel: [mandatory] meta, normal_cov
        t_cov // channel: [mandatory] meta, tumor_cov
        tools
        outdir

    main:
        if (tools && (tools.split(',').contains('fragcounter'))) {

            // Creating csv files to restart from this step
            n_cov.collectFile(keepHeader: true, skip: 1,sort: true, storeDir: "${params.outdir}/csv"){ meta, normal_cov, normal_bw ->
                patient       = meta.patient
                sample        = meta.id
                normal_cov = "${params.outdir}/Coverages/fragCounter_normal/${meta.id}/${normal_cov.baseName}.rds"
                normal_bw = "${params.outdir}/Coverages/fragCounter_normal/${meta.id}/${normal_bw.baseName}.bw" 
                ["fragcounter_N.csv", "patient,sample,sv_caller,vcf\n${patient},${sample},${normal_cov},${normal_bw}\n"]

            }

            t_cov.collectFile(keepHeader: true, skip: 1,sort: true, storeDir: "${params.outdir}/csv"){ meta, tumor_cov, tumor_bw ->
                patient       = meta.patient
                sample        = meta.id
                tumor_cov = "${params.outdir}/Coverages/fragCounter_tumor/${meta.id}/${tumor_cov.baseName}.rds"
                tumor_bw = "${params.outdir}/Coverages/fragCounter_tumor/${meta.id}/${tumor_bw.baseName}.bw" 
                ["fragcounter_T.csv", "patient,sample,sv_caller,vcf\n${patient},${tumor_cov},${tumor_bw}\n"]

            }
        } 
}