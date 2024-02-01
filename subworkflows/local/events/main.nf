include { EVENTS } from '../../../modules/local/events/main.nf'

workflow EVENTS {
    take:
    input           // required: format [val(meta), path(gGraph)]
    ref
    id

    main:
    versions            = Channel.empty()
    events_output       = Channel.empty()

    EVENTS(input, ref, id)

    events_output = EVENTS.out.events_output
    versions = EVENTS.out.versions

    emit:
    events_output

    versions
}
