include { EVENTS } from '../../../modules/local/events/main.nf'

workflow EVENTS {
    take:

    main:

    // define channels here for inputs
    // x = Channel.empty()
    EVENTS()


    // gather outputs
    // x = EVENTS.out.x

    emit:
    // x
}
