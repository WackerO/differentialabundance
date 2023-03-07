/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowDifferentialabundance.initialise(params, log)

def checkPathParamList = [ params.input ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
def exp_meta = [ "id": params.study_name  ]
if (params.input) { ch_input = Channel.of([ exp_meta, params.input ]) } else { exit 1, 'Input samplesheet not specified!' }
    
if (params.study_type == 'affy_array'){
    if (params.affy_cel_files_archive) { 
        ch_celfiles = Channel.of([ exp_meta, file(params.affy_cel_files_archive, checkIfExists: true) ]) 
    } else { 
        exit 1, 'CEL files archive not specified!' 
    }
} else{
    
    // If this is not an affy array, assume we're reading from a matrix
    
    if (params.matrix) { 
        matrix_file = file(params.matrix, checkIfExists: true)
        ch_in_raw = Channel.of([ exp_meta, matrix_file])
    } else { 
        exit 1, 'Input matrix not specified!' 
    }
}

// Check optional parameters
if (params.control_features) { ch_control_features = file(params.control_features, checkIfExists: true) } else { ch_control_features = [[],[]] } 
if (params.gsea_run) { gene_sets_file = file(params.gsea_gene_sets, checkIfExists: true) } else { gene_sets_file = [] } 

report_file = file(params.report_file, checkIfExists: true)
logo_file = file(params.logo_file, checkIfExists: true)
css_file = file(params.css_file, checkIfExists: true)
citations_file = file(params.citations_file, checkIfExists: true)

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

include { TABULAR_TO_GSEA_CHIP } from '../modules/local/tabular_to_gsea_chip'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { GUNZIP as GUNZIP_GTF                              } from '../modules/nf-core/gunzip/main'
include { UNTAR                                             } from '../modules/nf-core/untar/main.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS                       } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { SHINYNGS_STATICEXPLORATORY as PLOT_EXPLORATORY    } from '../modules/nf-core/shinyngs/staticexploratory/main'  
include { SHINYNGS_STATICDIFFERENTIAL as PLOT_DIFFERENTIAL  } from '../modules/nf-core/shinyngs/staticdifferential/main'  
include { SHINYNGS_VALIDATEFOMCOMPONENTS as VALIDATOR       } from '../modules/nf-core/shinyngs/validatefomcomponents/main'  
include { DESEQ2_DIFFERENTIAL                               } from '../modules/nf-core/deseq2/differential/main'    
include { LIMMA_DIFFERENTIAL                                } from '../modules/nf-core/limma/differential/main'    
include { CUSTOM_MATRIXFILTER                               } from '../modules/nf-core/custom/matrixfilter/main' 
include { ATLASGENEANNOTATIONMANIPULATION_GTF2FEATUREANNOTATION as GTF_TO_TABLE } from '../modules/nf-core/atlasgeneannotationmanipulation/gtf2featureannotation/main' 
include { GSEA_GSEA                                         } from '../modules/nf-core/gsea/gsea/main' 
include { CUSTOM_TABULARTOGSEAGCT                           } from '../modules/nf-core/custom/tabulartogseagct/main'
include { CUSTOM_TABULARTOGSEACLS                           } from '../modules/nf-core/custom/tabulartogseacls/main' 
include { RMARKDOWNNOTEBOOK                                 } from '../modules/nf-core/rmarkdownnotebook/main' 
include { AFFY_JUSTRMA as AFFY_JUSTRMA_RAW                  } from '../modules/nf-core/affy/justrma/main' 
include { AFFY_JUSTRMA as AFFY_JUSTRMA_NORM                 } from '../modules/nf-core/affy/justrma/main' 
include { PROTEUS                                           } from '../modules/nf-core/proteus/main' 

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow DIFFERENTIALABUNDANCE {

    // Set up some basic variables
    ch_versions = Channel.empty()
   
    // If we have affy array data in the form of CEL files we'll be deriving
    // matrix and annotation from them 
    
    if (params.study_type == 'affy_array'){

        // Uncompress the CEL files archive

        UNTAR ( ch_celfiles )

        ch_affy_input = ch_input
            .join(UNTAR.out.untar)

        // Run affy to derive the matrix. Reset the meta so it can be used to
        // define a prefix for different matrix flavours

        AFFY_JUSTRMA_RAW ( 
            ch_affy_input,
            [[],[]] 
        )
        AFFY_JUSTRMA_NORM ( 
            ch_affy_input,
            [[],[]] 
        )

        // Fetch affy outputs and reset the meta

        ch_in_raw = AFFY_JUSTRMA_RAW.out.expression     
        ch_in_norm = AFFY_JUSTRMA_NORM.out.expression
        
        ch_affy_platform_features = AFFY_JUSTRMA_RAW.out.annotation
    } else if (params.study_type == 'Px') {
        proteus_in = Channel.of([ exp_meta, file(params.input), file(params.matrix) ])
        PROTEUS(proteus_in)
    }

    //// Fetch or derive a feature annotation table

    // If user has provided a feature annotation table, use that

    if (params.features){
        ch_features = Channel.of([ exp_meta, file(params.features, checkIfExists: true)]) 
    } else if (params.study_type == 'affy_array'){
        ch_features = ch_affy_platform_features
    } else if (params.gtf){
        // Get feature annotations from a GTF file, gunzip if necessary
     
        file_gtf_in = file(params.gtf)
        file_gtf = [ [ "id": file_gtf_in.simpleName ], file_gtf_in ] 

        if ( params.gtf.endsWith('.gz') ){
            GUNZIP_GTF(file_gtf)
            file_gtf = GUNZIP_GTF.out.gunzip
            ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
        } 

        // Get a features table from the GTF and combine with the matrix and sample
        // annotation (fom = features/ observations/ matrix)

        GTF_TO_TABLE( file_gtf, [[ "id":""], []])
        ch_features = GTF_TO_TABLE.out.feature_annotation
            .map{
                tuple( exp_meta, it[1])
            }

        // Record the version of the GTF -> table tool

        ch_versions = ch_versions
            .mix(GTF_TO_TABLE.out.versions)
    }
    else{

        // Otherwise we can just use the matrix input 
        matrix_as_anno_filename = "matrix_as_anno.${matrix_file.getExtension()}"
        matrix_file.copyTo(matrix_as_anno_filename)
        ch_features = Channel.of([ exp_meta, file(matrix_as_anno_filename)])
    }

    // Channel for the contrasts file
    
   // ch_contrasts_file = Channel.from([[exp_meta, file(params.contrasts)]])

    // Check compatibility of FOM elements and contrasts

    if (params.study_type == 'affy_array'){
        ch_matrices_for_validation = ch_in_raw
            .join(ch_in_norm)
            .map{tuple(it[0], [it[1], it[2]])}
    }
    else{
        ch_matrices_for_validation = ch_in_raw
    }

 //   VALIDATOR(
 //       ch_input.join(ch_matrices_for_validation),
 //       ch_features,
 //       ch_contrasts_file
 //   )

    // For Affy, we've validated multiple input matrices for raw and norm,
    // we'll separate them out again here

    if (params.study_type == 'affy_array'){
        ch_validated_assays = VALIDATOR.out.assays
            .transpose()
            .branch {
                raw: it[1].name.contains('raw')
                normalised: it[1].name.contains('normalised')
            }
        ch_raw = ch_validated_assays.raw
        ch_norm = ch_validated_assays.normalised
        ch_matrix_for_differential = ch_norm        
    } else{
       //ch_raw = VALIDATOR.out.assays 
       // ch_matrix_for_differential = ch_raw
    }

    // Split the contrasts up so we can run differential analyses and
    // downstream plots separately. 
    // Replace NA strings that might have snuck into the blocking column

    // ch_contrasts = VALIDATOR.out.contrasts
    //     .map{it[1]}
    //     .splitCsv ( header:true, sep:'\t' )
    //     .map{
    //         it.blocking = it.blocking.replace('NA', '')
    //         if (!it.id){
    //             it.id = it.values().join('_')
    //         }
    //         it
    //     }

    // // Firstly Filter the input matrix

    // CUSTOM_MATRIXFILTER(
    //     ch_matrix_for_differential,
    //     VALIDATOR.out.sample_meta
    // )

    // // Prepare inputs for differential processes

    // ch_differential_inputs = ch_contrasts.combine(
    //     VALIDATOR.out.sample_meta
    //         .join(CUSTOM_MATRIXFILTER.out.filtered)     // -> meta, samplesheet, filtered matrix
    //         .map{ it.tail() } 
    // )

    // if (params.study_type == 'affy_array'){

    //     LIMMA_DIFFERENTIAL (
    //         ch_differential_inputs
    //     )
    //     ch_differential = LIMMA_DIFFERENTIAL.out.results 
        
    //     ch_versions = ch_versions
    //         .mix(LIMMA_DIFFERENTIAL.out.versions)
    
    //     ch_processed_matrices = ch_norm
    //         .map{ it.tail() }
    //         .first()
    // }
    // else{

    //     // Run the DESeq differential module, which doesn't take the feature
    //     // annotations 

    //     DESEQ2_DIFFERENTIAL (
    //         ch_differential_inputs,
    //         ch_control_features
    //     )

    //     // Let's make the simplifying assumption that the processed matrices from
    //     // the DESeq runs are the same across contrasts. We run the DESeq process
    //     // with matrices once for each contrast because DESeqDataSetFromMatrix()
    //     // takes the model, and the model can vary between contrasts if the
    //     // blocking factors included differ. But the normalised and
    //     // variance-stabilised matrices are not (IIUC) impacted by the model.

    //     ch_norm = DESEQ2_DIFFERENTIAL.out.normalised_counts.first()
    //     ch_vst = DESEQ2_DIFFERENTIAL.out.vst_counts.first()
    //     ch_differential = DESEQ2_DIFFERENTIAL.out.results 
    
    //     ch_versions = ch_versions
    //         .mix(DESEQ2_DIFFERENTIAL.out.versions)
        
    //     ch_processed_matrices = ch_norm
    //         .join(ch_vst)
    //         .map{ it.tail() }
    // }

    // // Run a gene set analysis where directed

    // // Currently, we're letting GSEA work on the expression data. In future we
    // // will allow use of GSEA preranked instead, which will work with the fold
    // // changes/ p values from DESeq2
    
    // if (params.gsea_run){    
    
    //     ch_gene_sets = Channel.from(gene_sets_file)

    //     // For GSEA, we need to convert normalised counts to a GCT format for
    //     // input, and process the sample sheet to generate class definitions
    //     // (CLS) for the variable used in each contrast        
        
    //     CUSTOM_TABULARTOGSEAGCT ( ch_norm )

    //     ch_contrasts_and_samples = ch_contrasts.combine( VALIDATOR.out.sample_meta.map { it[1] } )
    //     CUSTOM_TABULARTOGSEACLS(ch_contrasts_and_samples) 

    //     TABULAR_TO_GSEA_CHIP(
    //         VALIDATOR.out.feature_meta.map{ it[1] },
    //         [params.features_id_col, params.features_name_col]    
    //     )
    
    //     // The normalised matrix does not always have a contrast meta, so we
    //     // need a combine rather than a join here

    //     ch_gsea_inputs = CUSTOM_TABULARTOGSEAGCT.out.gct
    //         .map{ it.tail() }
    //         .combine(CUSTOM_TABULARTOGSEACLS.out.cls)
    //         .map{ tuple(it[1], it[0], it[2]) }
    //         .combine(ch_gene_sets)

    //     GSEA_GSEA( 
    //         ch_gsea_inputs,
    //         ch_gsea_inputs.map{ tuple(it[0].reference, it[0].target) }, // * 
    //         TABULAR_TO_GSEA_CHIP.out.chip.first()
    //     )
        
    //     // * Note: GSEA module currently uses a value channel for the mandatory
    //     // non-file arguments used to define contrasts, hence the indicated
    //     // usage of map to perform that transformation. An active subject of
    //     // debate

    //     ch_gsea_results = GSEA_GSEA.out.report_tsvs_ref
    //         .join(GSEA_GSEA.out.report_tsvs_target)

    //     // Record GSEA versions
    //     ch_versions = ch_versions
    //         .mix(TABULAR_TO_GSEA_CHIP.out.versions)
    //         .mix(GSEA_GSEA.out.versions)
    // }

    // // The exploratory plots are made by coloring by every unique variable used
    // // to define contrasts

    // ch_contrast_variables = ch_contrasts
    //     .map{
    //         [ "id": it.variable ]
    //     }
    //     .unique()

    // ch_all_matrices = VALIDATOR.out.sample_meta                 // meta, samples
    //     .join(VALIDATOR.out.feature_meta)                       // meta, samples, features
    //     .join(ch_raw)                                           // meta, samples, features, raw matrix
    //     .combine(ch_processed_matrices)                         // meta, samples, features, raw, norm, ...
    //     .map{
    //         tuple(it[0], it[1], it[2], it[3..it.size()-1])
    //     }
    //     .first()
 
    // ch_contrast_variables
    //     .combine(ch_all_matrices.map{ it.tail() })

    //     ch_contrast_variables
    //         .combine(ch_all_matrices.map{ it.tail() })
    
    // PLOT_EXPLORATORY(
    //     ch_contrast_variables
    //         .combine(ch_all_matrices.map{ it.tail() })
    // )

    // // Differential analysis using the results of DESeq2

    // PLOT_DIFFERENTIAL(
    //     ch_differential, 
    //     ch_all_matrices
    // )

    // // Gather software versions

    // ch_versions = ch_versions
    //     .mix(VALIDATOR.out.versions)
    //     .mix(PLOT_EXPLORATORY.out.versions)
    //     .mix(PLOT_DIFFERENTIAL.out.versions)

    // CUSTOM_DUMPSOFTWAREVERSIONS (
    //     ch_versions.unique().collectFile(name: 'collated_versions.yml')
    // )
    
    // // Generate a list of files that will be used by the markdown report

    // ch_report_file = Channel.from(report_file)
    //     .map{ tuple(exp_meta, it) }

    // ch_logo_file = Channel.from(logo_file)
    // ch_css_file = Channel.from(css_file)
    // ch_citations_file = Channel.from(citations_file)

    // ch_report_input_files = ch_all_matrices
    //     .map{ it.tail() }
    //     .map{it.flatten()}
    //     .combine(ch_contrasts_file.map{it.tail()})
    //     .combine(CUSTOM_DUMPSOFTWAREVERSIONS.out.yml)
    //     .combine(ch_logo_file)
    //     .combine(ch_css_file)
    //     .combine(ch_citations_file)
    //     .combine(ch_differential.map{it[1]}.toList())
 
    // if (params.gsea_run){ 
    //     ch_report_input_files = ch_report_input_files
    //         .combine(ch_gsea_results
    //             .map{it.tail()}.flatMap().toList()
    //         )
    // }

    // // Make a params list - starting with the input matrices and the relevant
    // // params to use in reporting

    // def report_file_names = [ 'observations', 'features' ] + 
    //     params.exploratory_assay_names.split(',').collect { "${it}_matrix".toString() } +
    //     [ 'contrasts_file', 'versions_file', 'logo', 'css', 'citations' ]

    // // Condition params reported on study type

    // def params_pattern = ~/^(study|observations|features|filtering|exploratory|differential|deseq2|gsea).*/
    // if (params.study_type == 'affy_array'){
    //     params_pattern = ~/^(study|observations|features|filtering|exploratory|differential|affy|limma|gsea).*/
    // } 

    // ch_report_params = ch_report_input_files
    //     .map{
    //         params.findAll{ k,v -> k.matches(params_pattern) } +
    //         [report_file_names, it.collect{ f -> f.name}].transpose().collectEntries() 
    //     }

    // // Render the final report
    
    // RMARKDOWNNOTEBOOK(
    //     ch_report_file,
    //     ch_report_params,
    //     ch_report_input_files
    // )

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
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
