workflow CHANNEL_HETPILEUPS_CREATE_CSV {
    take:
        sites_from_hetpileups
        tools
        outdir

    main:
        if (tools && (tools.split(',').contains('hetpileups'))) {

            // Creating csv files to restart from this step
            sites_from_hetpileups.collectFile(keepHeader: true, skip: 1,sort: true, storeDir: "${params.outdir}/csv"){ meta, sites ->
                patient       = meta.patient
                sample        = meta.id
                sites = "${params.outdir}/Hetpileups/hetpileups_wgs/${meta.id}/${sites.baseName}.bw"
                ["hetpileups.csv", "patient,sample,sites\n${patient},${sample},${sites}\n"]

            }
        }
}
