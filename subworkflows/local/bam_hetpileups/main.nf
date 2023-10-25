//
// BAM HETPILEUPS
//

include { HETPILEUPS } from '../../../modules/local/hetpileups/main.nf'

workflow BAM_HETPILEUPS {
    // defining inputs
    take:
    input                                             // required: Format should be [meta, tumor bam, normal bam] : can also provide cram
	filter
	max_depth
	hapmap_sites

    //Creating empty channels for output
    main:
    versions            = Channel.empty()
    het_pileups_wgs     = Channel.empty()

    HETPILEUPS(input, filter, max_depth, hapmap_sites)

    // initializing outputs from fragcounter
    versions          = HETPILEUPS.out.versions
    het_pileups_wgs   = HETPILEUPS.out.het_pileups_wgs

    //
    emit:
    het_pileups_wgs

    versions
}
