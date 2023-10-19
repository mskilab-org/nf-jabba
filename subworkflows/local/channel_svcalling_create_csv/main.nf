//
// CHANNEL_SV_CALLING_CREATE_CSV
//

workflow CHANNEL_SVCALLING_CREATE_CSV {
    take:
        vcf_from_sv_calling // channel: [mandatory] meta, vcf
        tools
        outdir

    main:
        if (tools && (tools.split(',').contains('svaba'))) {

            // Creating csv files to restart from this step
            vcf_from_sv_calling.collectFile(keepHeader: true, skip: 1,sort: true, storeDir: "${params.outdir}/csv"){ meta, vcf ->
                patient       = meta.patient
                sample        = meta.id
                sv_caller     = "Svaba"
                vcf = "${params.outdir}/SV_calling/SVABA/${meta.id}/${vcf.getName()}"
                ["svaba.csv", "patient,sample,sv_caller,vcf\n${patient},${sample},${sv_caller},${vcf}\n"]

            }
        } 
        
        if (tools && (tools.split(',').contains('gridss'))) {

            // Creating csv files to restart from this step
            vcf_from_sv_calling.collectFile(keepHeader: true, skip: 1,sort: true, storeDir: "${params.outdir}/csv"){ meta, vcf ->
                patient       = meta.patient
                sample        = meta.id
                sv_caller     = "GRIDSS"
                vcf = "${params.outdir}/SV_calling/GRIDSS/${meta.id}/${vcf.baseName}.bgz"
                ["gridss.csv", "patient,sample,sv_caller,vcf\n${patient},${sample},${sv_caller},${vcf}\n"]

            }
        }


}