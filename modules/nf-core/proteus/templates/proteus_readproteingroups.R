#!/usr/bin/env Rscript

################################################
################################################
## Functions                                  ##
################################################
################################################

#' Parse out options from a string without recourse to optparse
#'
#' @param x Long-form argument list like --opt1 val1 --opt2 val2
#'
#' @return named list of options and values similar to optparse

parse_args <- function(x){
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- lapply(args_list, function(x) scan(text=x, what='character', quiet = TRUE))

    # Ensure the option vectors are length 2 (key/ value) to catch empty ones
    args_vals <- lapply(args_vals, function(z){ length(z) <- 2; z})

    parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[! is.na(parsed_args)]
}

#' Flexibly read CSV or TSV files
#'
#' @param file Input file
#' @param header Passed to read.delim()
#' @param row.names Passed to read.delim()
#'
#' @return output Data frame

read_delim_flexible <- function(file, header = TRUE, row.names = NULL, check.names = F){

    ext <- tolower(tail(strsplit(basename(file), split = "\\\\.")[[1]], 1))

    if (ext == "tsv" || ext == "txt") {
        separator <- "\\t"
    } else if (ext == "csv") {
        separator <- ","
    } else {
        stop(paste("Unknown separator for", ext))
    }

    read.delim(
        file,
        sep = separator,
        header = header,
        row.names = row.names,
        check.names = check.names
    )
}

#' Round numeric dataframe columns to fixed decimal places by applying
#' formatting and converting back to numerics
#'
#' @param dataframe A data frame
#' @param columns Which columns to round (assumes all of them by default)
#' @param digits How many decimal places to round to?
#'
#' @return output Data frame
# TODO check if this is necessary
round_dataframe_columns <- function(df, columns = NULL, digits = 8){
    if (is.null(columns)){
        columns <- colnames(df)
    }

    df[,columns] <- format(
        data.frame(df[, columns], check.names = FALSE),
        nsmall = digits
    )

    # Convert columns back to numeric

    for (c in columns) {
        df[[c]][grep("^ *NA\$", df[[c]])] <- NA
        df[[c]] <- as.numeric(df[[c]])
    }
    df
}

################################################
################################################
## PARSE PARAMETERS FROM NEXTFLOW             ##
################################################
################################################

# I've defined these in a single array like this so that we could go back to an
# optparse-driven method in future with module bin/ directories, rather than
# the template

# Set defaults and classes

opt <- list(
    quant_file = '$quants',
    sample_file = '$samplesheet',
    contrast_variable = NULL,
    sample_id_col = 'sample',
    measure_col_prefix = 'Intensity',
    normfuns = 'normalizeMedian',
    plotSampleDistributions_method = 'violin',
    plotMV_loess = T,
    palette_name = 'Set1'
)
opt_types <- lapply(opt, class)

# Apply parameter overrides

args_opt <- parse_args('$task.ext.args')
for ( ao in names(args_opt)){
    if (! ao %in% names(opt)){
        stop(paste("Invalid option:", ao))
    }else{

        # Preserve classes from defaults where possible
        if (! is.null(opt[[ao]])){
            args_opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
        }
        opt[[ao]] <- args_opt[[ao]]
    }
}

# Check if required parameters have been provided

required_opts <- c('quant_file', 'sample_file', 'contrast_variable')
missing <- required_opts[unlist(lapply(opt[required_opts], is.null)) | ! required_opts %in% names(opt)]

if (length(missing) > 0){
    stop(paste("Missing required options:", paste(missing, collapse=', ')))
}

# Check file inputs are valid

for (file_input in c('quant_file', 'sample_file')){
    if (is.null(opt[[file_input]])) {
        stop(paste("Please provide", file_input), call. = FALSE)
    }

    if (! file.exists(opt[[file_input]])){
        stop(paste0('Value of ', file_input, ': ', opt[[file_input]], ' is not a valid file'))
    }
}

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################
mytmp <- tempdir()


# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#     install.packages("BiocManager",repos = "http://cran.us.r-project.org", lib=mytmp)
# }
# library("BiocManager", lib.loc=mytmp)
# if (!require("rmarkdown")){
#     BiocManager::install("rmarkdown", lib=mytmp)
#     library("rmarkdown", lib.loc=mytmp)
# }
# if (!require("miniUI")){
#     BiocManager::install("miniUI", lib=mytmp)
#     library("miniUI", lib.loc=mytmp)
# }
# if (!require("pkgdown")){
#     BiocManager::install("pkgdown", lib=mytmp)
#     library("pkgdown", lib.loc=mytmp)
# }
# if (!require("devtools")){
#     BiocManager::install("devtools", lib=mytmp)
#     library("devtools", lib.loc=mytmp)
# }
# if (!require("limma")){
#   BiocManager::install("limma", lib=mytmp)
#   library("limma", lib.loc=mytmp)
# }
# if (!require("ggplot2")){
#     BiocManager::install("ggplot2", lib=mytmp)
#     library("ggplot2", lib.loc=mytmp)
# }
# if (!require("ggplot2")) {
#     install.packages("ggplot2",repos = "http://cran.us.r-project.org", lib=mytmp)
#     library("ggplot2", lib.loc=mytmp)
# }
# devtools::install_github("tidyverse/ggplot2", lib=mytmp)
# library("ggplot2", lib.loc=mytmp)


# if (!require("plotly")){
#     BiocManager::install("plotly", lib=mytmp)
#     library("plotly", lib.loc=mytmp)
# }
# if (!require("proteus")){
#     devtools::install_github("bartongroup/Proteus", lib=mytmp, 
#                             build_opts= c("--no-resave-data", "--no-manual"), build_vignettes=FALSE)
#     library("proteus", lib.loc=mytmp)
# }


library(limma)
library(plotly)
library(proteus)

################################################
################################################
## READ IN COUNTS FILE AND SAMPLE METADATA    ##
################################################
################################################
#opt\$quant_file, 
quant.table <-
    read_delim_flexible(
        file = opt\$quant_file,
        check.names = FALSE,
        row.names = 1
    )
write.table(quant.table, file="/home-link/iivow01/git/differentialabundance/error/quanttablllllllllle.tsv", sep="\t")
sample.sheet <-
    read_delim_flexible(
        file = opt\$sample_file,
        check.names=FALSE
    )

# Deal with spaces that may be in sample column
#opt\$sample_id_col <- make.names(opt\$sample_id_col)

if (! opt\$sample_id_col %in% colnames(sample.sheet)){
    stop(paste0("Specified sample ID column '", opt\$sample_id_col, "' is not in the sample sheet"))
}

# Add metadata columns that are necessary for proteus

sample.sheet\$sample <- sample.sheet[[opt\$sample_id_col]]



#opt\$contrast_variable <- make.names(opt\$contrast_variable)
sample.sheet\$condition <- sample.sheet[[opt\$contrast_variable]]

# Add prefix for proteinGroups measurement columns to the sample IDs from the sampesheet
measure.cols <- setNames(paste0(opt\$measure_col_prefix, sample.sheet[[opt\$sample_id_col]]), sample.sheet[[opt\$sample_id_col]])

# TODO check if this can happen for proteingroups
# Sample sheet can have duplicate rows for multiple sequencing runs, so uniqify
# before assigning row names

#sample.sheet <- sample.sheet[! duplicated(sample.sheet[[opt\$sample_id_col]]), ]
#rownames(sample.sheet) <- sample.sheet[[opt\$sample_id_col]]

# Check that all samples specified in the input sheet are present in the quants
# table

missing_columns <- paste0(opt\$measure_col_prefix, sample.sheet[[opt\$sample_id_col]])
missing_samples <-
    (sample.sheet[[opt\$sample_id_col]])[!missing_columns %in% colnames(quant.table)]

# TODO: Consider if this auto-filter should be kept or removed (probably removed, otherwise I also have to deal with makenames)
#sample.sheet <- sample.sheet[!(rownames(sample.sheet) %in% missing_samples),]

if (length(missing_samples) > 0) {
    stop(paste(
        length(missing_columns),
        'specified samples do not have a(n)',
        opt\$measure_col_prefix,
        'column in quant table. The following columns are missing:',
        paste(missing_columns, collapse = ', ')
    ))
} else{
    # Save any non-quant data, with gene metadata etc we might need later
    # TODO: Maybe just save the whole quant file? (or not; not sure the rest is ever needed)
    nonquant.table <-
        quant.table[, !colnames(quant.table) %in% paste0(opt\$measure_col_prefix, sample.sheet[[opt\$sample_id_col]]), drop = FALSE]
}

################################################
################################################
## CHECK AND FORMAT NORMFUN AND FILTERFUN     ##
################################################
################################################

valid_normfuns <- c("normalizeMedian", "normalizeQuantiles")
normfuns <- opt\$normfuns

# Check validity of normfun(s)
invalid_normfuns <- normfuns[!(normfuns %in% valid_normfuns)]
if (length(invalid_normfuns)>0) {
    stop(paste0("Invalid normfuns argument(s): ",
        paste(invalid_normfuns, collapse=", "),
        ". Valid normfuns are: ",
        paste(valid_normfuns, collapse=", "),
        "."))
}

################################################
################################################
## Run Proteus processes and generate outputs ##
################################################
################################################

# TODO
output_prefix <- "output_prefix"




proteinGroups <- readProteinGroups(
    file=opt\$quant_file,
    meta=sample.sheet,
    measure.cols=measure.cols,
    data.cols=proteus::proteinColumns
)
capture.output(proteinGroups, file="/home-link/iivow01/git/differentialabundance/error/proteingroups")
capture.output(str(proteinGroups), file="/home-link/iivow01/git/differentialabundance/error/proteingroupsstr")

write("1", file="/home-link/iivow01/git/differentialabundance/error/status")

write.table(proteinGroups\$tab, file="/home-link/iivow01/git/differentialabundance/error/tab", quote=F)
write("2", file="/home-link/iivow01/git/differentialabundance/error/status", append=T)

# Generate non-normalized distribution plot
png('proteus.nonnormalized_distributions.png', width = 5*300, height = 5*300, res = 300, pointsize = 8)
print(plotSampleDistributions(proteinGroups, title="Non-normalized sample distributions", fill="condition", method=opt\$plotSampleDistributions_method) + scale_fill_brewer(palette=opt\$palette_name))
dev.off()

# R object for other processes to use
saveRDS(proteinGroups, file = 'proteus.raw_proteingroups.rds')


# Generate plots for all requested normalizations; also, save
# normalized protein groups for limma
for (normfun in normfuns) {
    proteinGroups.normalized <- normalizeData(proteinGroups, norm.fun = eval(parse(text=normfun))) # Proteus also accepts other norm.funs, e.g. from limma

    png(paste0('proteus.', normfun, '_normalized_distributions.png'), width = 5*300, height = 5*300, res = 300, pointsize = 8) 
    print(plotSampleDistributions(proteinGroups.normalized, title=paste("Sample distributions after applying", normfun), fill="condition", method=opt\$plotSampleDistributions_method) + scale_fill_brewer(palette=opt\$palette_name))
    dev.off()
    
    png(paste0('proteus.', normfun, '_normalized_mean_variance_relationship.png'), width = 5*300, height = 5*300, res = 300, pointsize = 8) 
    print(plotMV(proteinGroups.normalized, with.loess=opt\$plotMV_loess) + scale_fill_distiller(palette=opt\$palette_name)) #, title=paste("Sample mean variance relationship after applying", normfun)
    dev.off()
    
    png(paste0('proteus.', normfun, '_normalized_dendrogram.png'), width = 5*300, height = 5*300, res = 300, pointsize = 8)
    print(plotClustering(proteinGroups.normalized), title=paste("Sample clustering after applying", normfun))
    dev.off()

    summary <- summary(proteinGroups.normalized)
    
    # R object for other processes to use
    saveRDS(proteinGroups.normalized, file = paste0('proteus.', normfun, 'normalized_proteingroups.rds'))
write("3", file="/home-link/iivow01/git/differentialabundance/error/status", append=T)

    # Write normalized count matrix
    write.table(
        data.frame(
            gene_id = rownames(proteinGroups.normalized\$tab),
            proteinGroups.normalized\$tab,
            check.names = FALSE
        ),
        file = paste(output_prefix, 'proteus', normfun, 'normalized_proteingroups_tab', 'tsv', sep = '.'),
        col.names = TRUE,
        row.names = FALSE,
        sep = '\t',
        quote = FALSE
    )
    write.table(
        data.frame(
            gene_id = rownames(proteinGroups.normalized\$tab),
            proteinGroups.normalized\$tab,
            check.names = FALSE
        ),
        file = "/home-link/iivow01/git/differentialabundance/error/waaaaas.normalizeMedian.tsv",
        col.names = TRUE,
        row.names = FALSE,
        sep = '\t',
        quote = FALSE
    )
    write.table(
        data.frame(
            gene_id = rownames(proteinGroups.normalized\$tab),
            proteinGroups.normalized\$tab,
            check.names = FALSE
        ),
        file = paste(output_prefix, 'proteus', normfun, 'normalized_proteingroups_tab2', 'tsv', sep = '.'),
        col.names = TRUE,
        row.names = FALSE,
        sep = '\t',
        quote = FALSE
    )
    
 write("5", file="/home-link/iivow01/git/differentialabundance/error/status", append=T)
   
bla <- paste("/home-link/iivow01/git/differentialabundance/error/waaaaas", normfun, 'tsv', sep = '.')
}




matrix_data <- read.delim(
    bla,
    check.names = FALSE,
    header = TRUE,
    sep = "\t",
    row.names = 1
)

write.table(matrix_data, file="/home-link/iivow01/git/differentialabundance/error/gehtdas.tsv", sep = "\t", quote=F)

# Write non-normalized count matrix
write.table(
    data.frame(
        gene_id = rownames(proteinGroups\$tab),
        proteinGroups\$tab,
        check.names = FALSE
    ),
    file = paste(output_prefix, 'proteus', 'raw_proteingroups_tab', 'tsv', sep = '.'),
    col.names = TRUE,
    row.names = FALSE,
    sep = '\t',
    quote = FALSE
)

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink("R_sessionInfo.log")
print(sessionInfo())
sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################
#TODO
r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
limma.version <- as.character(packageVersion('limma'))
plotly.version <- as.character(packageVersion('plotly'))
proteus.version <- as.character(packageVersion('proteus'))
#TODO: change mparker2
# writeLines(
#     c(
#         '"${task.process}":',
#         paste('    r-base:', r.version),
#         paste('    bioconductor-limma:', limma.version),
#         paste('    r-plotly:', plotly.version),
#         paste('    mparker2-proteus:', proteus.version),
#     ),
# 'versions.yml')

################################################
################################################
################################################
################################################
