/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowDifferentialabundance.initialise(params, log)

def checkPathParamList = [ params.matrix, params.input ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.matrix) { ch_counts = Channel.fromPath(params.matrix) } else { exit 1, 'Gene counts not specified!' }
if (params.input) { ch_input = Channel.fromPath(params.input) } else { exit 1, 'Samplesheet not specified!' }
if (params.contrasts) { ch_contrasts = Channel.fromPath(params.contrasts) } else { exit 1, 'Contrasts not specified!' }

// Check optinal parameters
if (params.control_features) { ch_control_features = file(params.control_features) } else { ch_control_features = [[],[]] }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { IMPORTRNASEQCOUNTS                               } from '../modules/local/importrnaseqcounts/main'

//
// MODULE: Installed directly from nf-core/modules
//
include { GUNZIP as GUNZIP_GTF                              } from '../modules/nf-core/gunzip/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                       } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { SHINYNGS_STATICEXPLORATORY as PLOT_EXPLORATORY    } from '../modules/nf-core/shinyngs/staticexploratory/main'
include { SHINYNGS_STATICDIFFERENTIAL as PLOT_DIFFERENTIAL  } from '../modules/nf-core/shinyngs/staticdifferential/main'
include { SHINYNGS_VALIDATEFOMCOMPONENTS as VALIDATOR       } from '../modules/nf-core/shinyngs/validatefomcomponents/main'
include { DESEQ2_DIFFERENTIAL                               } from '../modules/nf-core/deseq2/differential/main'
include { ATLASGENEANNOTATIONMANIPULATION_GTF2FEATUREANNOTATION as GTF_TO_TABLE } from '../modules/nf-core/atlasgeneannotationmanipulation/gtf2featureannotation/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow DIFFERENTIALABUNDANCE {

    ch_versions = Channel.empty()
    exp_meta = Channel.from([ "id": params.study_name  ])

    switch (params.counts_type) {
        case 'salmon':
            break
        case 'rnaseq_featurecounts':
            IMPORTRNASEQCOUNTS (
                ch_counts,
                ch_input
            )
            ch_counts = IMPORTRNASEQCOUNTS.out.ch_counts
            ch_input = IMPORTRNASEQCOUNTS.out.ch_input
            ch_versions = ch_versions.mix(IMPORTRNASEQCOUNTS.out.versions)
            break
    }

    if (params.gtf) {

        // Get feature annotations from a GTF file, gunzip if necessary

        file_gtf_in = file(params.gtf)
        file_gtf = [ [ "id": file_gtf_in.simpleName ], file_gtf_in ]

        if ( params.gtf.endsWith('.gz') ){
            GUNZIP_GTF(file_gtf)
            file_gtf = GUNZIP_GTF.out.gunzip
            ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
        }

        // Get a features table from the GTF and combine with the counts and sample
        // annotation (fom = features/ observations/ counts)

        GTF_TO_TABLE( file_gtf, [[ "id":""], []])

        // Combine features with the observations and matrices to create a FOM
        // where these things can travel together
        ch_fom = GTF_TO_TABLE.out.feature_annotation
            .map{
                tuple([ "id": params.study_name  ], file(params.input), it[1])               // exp_meta.combine(ch_input).combine(Channel.from(it[1])).combine(ch_counts)
            }
            .combine(ch_counts)


    } else {
        ch_fom = Channel.of([ [ "id": params.study_name  ], file(params.input), file("FALSE")]).merge(ch_counts)
    }

    // Channel for the contrasts file
    ch_contrasts_file = Channel.from([[exp_meta, file(params.contrasts)]])

    // Check compatibility of FOM elements and contrasts
    VALIDATOR(
        ch_fom,
        ch_contrasts_file
    )

    // Split the contrasts up so we can run differential analyses and
    // downstream plots separately.
    // Replace NA strings that might have snuck into the blocking column

    ch_contrasts = VALIDATOR.out.contrasts
        .map{it[1]}
        .splitCsv ( header:true, sep:'\t' )
        .map{
            it.blocking = it.blocking.replace('NA', '')
            it
        }

    // Run the DESeq differential module, which doesn't take the feature
    // annotations

    ch_samples_and_matrix = VALIDATOR.out.fom.map{
        tuple(it[1], it[3])
    }
    ch_contrasts.combine(ch_samples_and_matrix)
    DESEQ2_DIFFERENTIAL (
        ch_contrasts.combine(ch_samples_and_matrix),
        ch_control_features
    )

    // Let's make the simplifying assumption that the processed matrices from
    // the DESeq runs are the same across contrasts. We run the DESeq process
    // with matrices once for each contrast because DESeqDataSetFromMatrix()
    // takes the model, and the model can vary between contrasts if the
    // blocking factors included differ. But the normalised and
    // variance-stabilised matrices are not (IIUC) impacted by the model.

    ch_matrices = DESEQ2_DIFFERENTIAL.out.normalised_counts
        .join(DESEQ2_DIFFERENTIAL.out.vst_counts)
        .map{ it.tail() }
        .first()

    // The exploratory plots are made by coloring by every unique variable used
    // to define contrasts

    ch_contrast_variables = ch_contrasts
        .map{
            [ "id": it.variable ]
        }
        .unique()

    ch_fom_plot_inputs = VALIDATOR.out.fom
        .combine(ch_matrices)                         // Add processed matrices to what we have in the FOM
        .map{
            tuple(it[0], it[1], it[2], [ it[3], it[4], it[5] ]) // Remove the experiment meta and group the matrices
        }


    if (params.gtf) {
        PLOT_EXPLORATORY(
            ch_contrast_variables
                .combine(ch_fom_plot_inputs.map{ it.tail() })
        )

        // Differential analysis using the results of DESeq2

        PLOT_DIFFERENTIAL(
            DESEQ2_DIFFERENTIAL.out.results,
            ch_fom_plot_inputs.first()
        )
    }

    // Gather software versions

    // TODO: Do these have to be separate for correct order? (I.e. the plots at the very end of out.versions?)
    if (params.gtf) {
        ch_versions = ch_versions
            .mix(GTF_TO_TABLE.out.versions)
            .mix(PLOT_EXPLORATORY.out.versions)
            .mix(PLOT_DIFFERENTIAL.out.versions)
    }
    ch_versions = ch_versions
        .mix(VALIDATOR.out.versions)
        .mix(DESEQ2_DIFFERENTIAL.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.adaptivecard(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
