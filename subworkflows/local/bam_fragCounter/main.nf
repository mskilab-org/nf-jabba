//
// BAM FRAGCOUNTER
//

include { FRAGCOUNTER } from '../../../modules/local/fragcounter/main.nf'

workflow BAM_FRAGCOUNTER {
    // defining inputs
    take:
    input                                             // required: Format should be [meta, bam, bai] : can also provide cram & crai
    midpoint
    windowsize
    gcmapdir
    minmapq
    fasta                                             // Required: if using cram files instead of bam. In our case we are using cram files.
    fasta_fai
    paired
    exome

    //Creating empty channels for output
    main:
    versions          = Channel.empty()
    fragcounter_cov   = Channel.empty()
    corrected_bw      = Channel.empty()

    FRAGCOUNTER(input, midpoint, windowsize, gcmapdir, minmapq, fasta, fasta_fai, paired, exome) // We are keeping cov empty because we don't use any input coverage for fragcounter
    //FRAGCOUNTER(input, midpoint, windowsize, gcmapdir, minmapq, paired, exome)

    // initializing outputs from fragcounter
    fragcounter_cov   = FRAGCOUNTER.out.fragcounter_cov
    versions          = FRAGCOUNTER.out.versions
    corrected_bw      = FRAGCOUNTER.out.corrected_bw

    //
    emit:
    fragcounter_cov
    corrected_bw

    versions    
}