#!/usr/bin/env Rscript

#
# LMDseq tximport script - Simplified and robust version
# Created by Pierluigi Di Chiaro
# Licensed under MIT License - see LICENSE file for details
#

# Load required libraries
suppressPackageStartupMessages({
  library("optparse")
  library("tidyverse")
  library("tximport")
  library("rtracklayer")
})

# Command line options
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Space-separated list of kallisto output directories", metavar="character"),
  make_option(c("-g", "--gtf"), type="character", default=NULL,
              help="GTF annotation file", metavar="character"),
  make_option(c("-r", "--reference"), type="character", default=NULL,
              help="Reference annotation file (optional)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output file name", metavar="character")
)

optparser <- OptionParser(option_list=option_list)
opt <- parse_args(optparser)

# Validate required arguments
if (is.null(opt$input)) {
  print_help(optparser)
  stop("Input directories must be provided", call.=FALSE)
}
if (is.null(opt$gtf)) {
  print_help(optparser)
  stop("GTF file must be provided", call.=FALSE)
}
if (is.null(opt$output)) {
  print_help(optparser)
  stop("Output file must be provided", call.=FALSE)
}

# Parse input parameters
input_dirs <- trimws(strsplit(opt$input, " ")[[1]])
gtf_file <- opt$gtf
ref_file <- opt$reference
output_file <- opt$output

cat("Input directories:", paste(input_dirs, collapse=", "), "\n")
cat("GTF file:", gtf_file, "\n")
cat("Reference file:", ifelse(is.null(ref_file), "None", ref_file), "\n")
cat("Output file:", output_file, "\n")

# Validate inputs
for(dir in input_dirs) {
    if(!dir.exists(dir)) {
        stop("Directory not found: ", dir, call.=FALSE)
    }
    abundance_file <- file.path(dir, "abundance.h5")
    if(!file.exists(abundance_file)) {
        stop("abundance.h5 not found in: ", dir, call.=FALSE)
    }
}

if(!file.exists(gtf_file)) {
    stop("GTF file not found: ", gtf_file, call.=FALSE)
}

# Read GTF and create transcript-to-gene mapping
cat("Reading GTF file...\n")
gtf_data <- import(gtf_file)
gtf_df <- as.data.frame(gtf_data)

transcripts <- gtf_df[gtf_df$type == "transcript", ]
if(nrow(transcripts) == 0) {
    stop("No transcripts found in GTF file", call.=FALSE)
}

tx2gene <- transcripts[, c("transcript_id", "gene_id")]
tx2gene <- tx2gene[complete.cases(tx2gene), ]

cat("Found", nrow(tx2gene), "transcript-to-gene mappings\n")

# Prepare file paths for tximport
abundance_files <- file.path(input_dirs, "abundance.h5")
names(abundance_files) <- basename(input_dirs)

# Run tximport
cat("Running tximport...\n")
txi <- tximport(abundance_files, type = "kallisto", tx2gene = tx2gene)

# Extract count matrix
counts_matrix <- txi$counts
sample_names <- colnames(counts_matrix)

cat("Imported expression data for", nrow(counts_matrix), "genes and", ncol(counts_matrix), "samples\n")

# Create output matrix
if(!is.null(ref_file) && file.exists(ref_file)) {
    cat("Adding gene annotations from reference file...\n")
    
    # Try to read reference file
    tryCatch({
        ref_data <- read.delim(ref_file, stringsAsFactors = FALSE)
        
        # Create gene annotation dataframe
        gene_info <- data.frame(
            gene_id = rownames(counts_matrix),
            stringsAsFactors = FALSE
        )
        
        # Add basic gene info if available in GTF
        genes_gtf <- gtf_df[gtf_df$type == "gene", ]
        if(nrow(genes_gtf) > 0) {
            gene_names <- setNames(genes_gtf$gene_name, genes_gtf$gene_id)
            gene_info$gene_name <- gene_names[gene_info$gene_id]
        }
        
        # Merge with reference if possible
        if("gene_id" %in% names(ref_data) || "ENSEMBL" %in% names(ref_data)) {
            ref_col <- ifelse("gene_id" %in% names(ref_data), "gene_id", "ENSEMBL")
            gene_info <- merge(gene_info, ref_data, by.x = "gene_id", by.y = ref_col, all.x = TRUE, sort = FALSE)
        }
        
        # Combine annotation with expression data
        final_matrix <- cbind(gene_info, counts_matrix[gene_info$gene_id, , drop = FALSE])
        
    }, error = function(e) {
        cat("Warning: Could not process reference file:", e$message, "\n")
        cat("Creating basic matrix without extended annotations\n")
        final_matrix <- data.frame(gene_id = rownames(counts_matrix), counts_matrix)
    })
    
} else {
    cat("Creating basic gene expression matrix\n")
    final_matrix <- data.frame(gene_id = rownames(counts_matrix), counts_matrix)
}

# Write output
cat("Writing output matrix...\n")
write.table(final_matrix, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Successfully created expression matrix:\n")
cat("  Genes:", nrow(final_matrix), "\n")
cat("  Samples:", ncol(final_matrix) - 1, "\n")  # Subtract 1 for gene_id column
cat("  Output file:", output_file, "\n")

cat("=== tximport completed successfully ===\n")

