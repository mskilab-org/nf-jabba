include { EVENTS } from '../../../modules/local/events/main.nf'

workflow GGRAPH_EVENTS {
    take:
    input           // required: format [val(meta), path(gGraph)]
    ref

    main:
    versions            = Channel.empty()
    events_output       = Channel.empty()

    EVENTS(input, ref)

    events_output = EVENTS.out.events_output
    versions = EVENTS.out.versions

    emit:
    events_output

    versions
}
