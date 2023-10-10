//
// CHANNEL_BASERECALIBRATOR_CREATE_CSV
//

workflow CHANNEL_BASERECALIBRATOR_CREATE_CSV {
    take:
        cram_table_bqsr // channel: [mandatory] meta, cram, crai, table
        tools
        skip_tools
        save_output_as_bam
        outdir

    main:
        // Creating csv files to restart from this step
        if (!(skip_tools && (skip_tools.split(',').contains('markduplicates')))) {
            cram_table_bqsr.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${outdir}/csv") { meta, cram, crai, table ->

                patient = meta.patient
                sample  = meta.sample
                sex     = meta.sex
                status  = meta.status
                suffix_aligned = save_output_as_bam ? "bam" : "cram"
                suffix_index   = save_output_as_bam ? "bam.bai" : "cram.crai"
                cram = "${outdir}/Alignment/markduplicates/${sample}/${cram.baseName}.${suffix_aligned}"
                crai = "${outdir}/Alignment/markduplicates/${sample}/${crai.baseName.minus(".cram")}.${suffix_index}"
                table = "${outdir}/Alignment/recal_table/${sample}/${sample}.recal.table"

                type = save_output_as_bam ? "bam" : "cram"
                type_index = save_output_as_bam ? "bai" : "crai"

                ["markduplicates.csv", "patient,sex,status,sample,${type},${type_index},table\n${patient},${sex},${status},${sample},${cram},${crai},${table}\n"]
            }
        } else {
            cram_table_bqsr.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${outdir}/csv") { meta, cram, crai, table ->
                patient = meta.patient
                sample  = meta.sample
                sex     = meta.sex
                status  = meta.status
                suffix_aligned = save_output_as_bam ? "bam" : "cram"
                suffix_index   = save_output_as_bam ? "bam.bai" : "cram.crai"
                cram = "${outdir}/Alignment/${sample}/Mapped/${cram.baseName}.${suffix_aligned}"
                crai = "${outdir}/Alignment/${sample}/Mapped/${crai.baseName.minus(".cram")}.${suffix_index}"
                table = "${outdir}/Alignment/${sample}/recal_table/${sample}.recal.table"

                type = save_output_as_bam ? "bam" : "cram"
                type_index = save_output_as_bam ? "bai" : "crai"

                ["sorted.csv", "patient,sex,status,sample,${type},${type_index},table\n${patient},${sex},${status},${sample},${cram},${crai},${table}\n"]
            }
        }
}
