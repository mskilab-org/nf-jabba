# Refactor
- Cut down on the number of default parameters in the config (possibly at the
  package level). Packages/processes/workflows shouldn't have so many
  parameters (see: JaBbA as the worst offendor), it indicates either
  overparameterization or a dependency tree in the parameter space.
- Add *_create_csv methods for remaining tools to generate samplesheets to
  start from those tools
- Move parameter specification from nextflow.config to module specific config,
  then import them into nextflow config. This keeps the process and its parameter
  configuration tightly coupled--but loses a central interface for modifying
  all the defaults in one place, which is useful if changing one default would
  affect defaults in a different parameter.
- Swtich from the "step" based conditionals to something more declarative (e.g
  cases). Having the program run a block of code by checking if the starting
  step is included in a list of steps is enormously cumbersome, cases should
  make things more parsimoinous and readable.
- Use functions for repetitive declarations/imports and hold their required
  variables in arrays/maps.
- Clean up output channel attribute extraction. Currently very repetitive,
  could be replaced with functions.


## Code Samples
```groovy
// Define the pipeline steps in order
def pipelineSteps = ['alignment', 'sv_calling', 'coverage', 'segmentation', 'ploidy_calling', 'junction_balance']

// User input for the starting step
def startingStep = 'coverage' // This would be provided by the user

// Find the index of the starting step in the pipeline
int startIndex = pipelineSteps.indexOf(startingStep)

// Validate if the starting step is valid
if (startIndex == -1) {
    error "Invalid starting step: $startingStep"
}

// Execute the pipeline steps starting from the starting step
pipelineSteps.drop(startIndex).each { step ->
    switch (step) {
        case 'alignment':
            // run alignment step
            println "Running alignment"
            break
        case 'sv_calling':
            // run sv_calling step
            println "Running sv_calling"
            break
        case 'coverage':
            // run coverage step
            println "Running coverage"
            break
        case 'segmentation':
            // run segmentation step
            println "Running segmentation"
            break
        case 'ploidy_calling':
            // run ploidy_calling step
            println "Running ploidy_calling"
            break
        case 'junction_balance':
            // run junction_balance step
            println "Running junction_balance"
            break
    }
}
```

```groovy
// Define a closure to include a subworkflow from a local path
def includeLocalSubworkflow(String subworkflowName, String alias = null) {
    alias = alias ?: subworkflowName.toUpperCase().replaceAll(/[^A-Z0-9_]/, '_')
    return "include { $alias } from '../subworkflows/local/$subworkflowName/main'"
}

// Define a list of subworkflows to include
def subworkflows = [
    'channel_align_create_csv',
    'channel_markduplicates_create_csv',
    'channel_baserecalibrator_create_csv',
    //...
]

// Include subworkflows using the closure
subworkflows.each { subworkflow ->
    println includeLocalSubworkflow(subworkflow)
}

// For modules that are included multiple times with different aliases, use a map
def samtoolsConvertAliases = [
    'BAM_TO_CRAM': 'bam_convert_samtools',
    'BAM_TO_CRAM_MAPPING': 'bam_convert_samtools',
    // ...
]

samtoolsConvertAliases.each { alias, subworkflow ->
    println includeLocalSubworkflow(subworkflow, alias)
}

// For modules that are included with the same alias but different parameters, include them once
// Then use with different parameters in the workflow logic
println "include { SAMTOOLS_CONVERT } from '../modules/nf-core/samtools/convert/main'"




