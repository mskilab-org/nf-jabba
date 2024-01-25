//
// BAM FRAGCOUNTER
//

include { FRAGCOUNTER } from '../../../modules/local/fragcounter/main.nf'
include { REBIN_RAW_FRAGCOUNTER } from '../../../modules/local/fragcounter/main.nf'

workflow BAM_FRAGCOUNTER {
    // defining inputs
    take:
    input                                             // required: Format should be [meta, bam, bai] : can also provide cram & crai
    midpoint
    windowsize
    gcmapdir
    minmapq
    paired
    exome

    //Creating empty channels for output
    main:
    versions          = Channel.empty()
    fragcounter_raw_cov   = Channel.empty()
    fragcounter_cov   = Channel.empty()
    corrected_bw      = Channel.empty()
    rebinned_raw_cov  = Channel.empty()

    FRAGCOUNTER(input, midpoint, windowsize, gcmapdir, minmapq, [], [], paired, exome) // We are keeping cov empty because we don't use any input coverage for fragcounter
    //FRAGCOUNTER(input, midpoint, windowsize, gcmapdir, minmapq, paired, exome)

    // initializing outputs from fragcounter
    fragcounter_raw_cov   = FRAGCOUNTER.out.fragcounter_raw_cov
    fragcounter_cov   = FRAGCOUNTER.out.fragcounter_cov
    versions          = FRAGCOUNTER.out.versions
    corrected_bw      = FRAGCOUNTER.out.corrected_bw

    REBIN_RAW_FRAGCOUNTER(fragcounter_cov, "reads", 1000)

    rebinned_raw_cov  = REBIN_RAW_FRAGCOUNTER.out.raw_fragcounter_cov_1kb

    //
    emit:
    fragcounter_raw_cov
    fragcounter_cov
    rebinned_raw_cov
    corrected_bw

    versions    
}