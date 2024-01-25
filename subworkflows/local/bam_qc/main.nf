

include { QUALIMAP_BAMQC } from '../../../modules/nf-core/qualimap/bamqc/main.nf'
include { MOSDEPTH } from '../../../modules/nf-core/mosdepth/main.nf'
include { SAMTOOLS_STATS } from '../../../modules/nf-core/samtools/stats/main.nf'


workflow BAM_QC_QUALIMAP_MOSDEPTH_SAMTOOLS {
    take:
    bam                          // channel: [mandatory] [ meta, cram, crai ]
    fasta                         // channel: [mandatory] [ fasta ]
    intervals

    main:
    versions = Channel.empty()
    reports = Channel.empty()

    qualimap_bam = bam.map { meta, bam, bai ->
        [meta, bam]
    }
    // Reports run on BAM
    QUALIMAP_BAMQC(qualimap_bam, intervals.map{ meta, gff -> gff })

    SAMTOOLS_STATS(bam, fasta.map{ it -> [ [ id:'fasta' ], it ] })

    MOSDEPTH(bam.combine(intervals.map{ meta, bed -> [ bed?:[] ] }), fasta.map{ it -> [ [ id:'fasta' ], it ] })

    // Gather all reports generated
    reports = reports.mix(SAMTOOLS_STATS.out.stats)
    reports = reports.mix(MOSDEPTH.out.global_txt)
    reports = reports.mix(MOSDEPTH.out.regions_txt)
    reports = reports.mix(MOSDEPTH.out.summary_txt)
    //reports = reports.mix()

    // Gather versions of all tools used
    versions = versions.mix(MOSDEPTH.out.versions)
    versions = versions.mix(SAMTOOLS_STATS.out.versions.first())

    emit:
    reports

    versions // channel: [ versions.yml ]
}
