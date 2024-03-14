include { FUSIONS } from '../../../modules/local/fusions/main.nf'

workflow GGRAPH_FUSIONS {
    take:
    input           // required: format [val(meta), path(gGraph)]
    gencode

    main:
    versions            = Channel.empty()
    fusions_output       = Channel.empty()

    FUSIONS(input, gencode)

    fusions_output = FUSIONS.out.fusions_output
    versions = FUSIONS.out.versions

    emit:
    fusions_output

    versions
}
