#!/usr/bin/env Rscript

# MAF Objects Merger using read.maf from maftools
# This script reads MAF files individually using read.maf, updates sample names, and merges them

library(maftools)
library(data.table)

# Function to merge MAF objects with filename as Tumor_Sample_Barcode
merge_maf_objects <- function(maf_directory = ".", output_file = NULL) {
  
  # Get all MAF files in the directory
  maf_files <- list.files(path = maf_directory, 
                         pattern = "\\.(maf|MAF)$", 
                         full.names = TRUE)
  
  if (length(maf_files) == 0) {
    stop("No MAF files found in the specified directory")
  }
  
  cat("Found", length(maf_files), "MAF files:\n")
  cat(paste(basename(maf_files), collapse = "\n"), "\n\n")
  
  # Initialize list to store MAF objects
  maf_objects <- list()
  
  # Process each MAF file
  for (i in seq_along(maf_files)) {
    file_path <- maf_files[i]
    
    # Extract filename without extension for sample barcode
    sample_name <- tools::file_path_sans_ext(basename(file_path))
    
    cat("Processing:", sample_name, "\n")
    
    # Read MAF file using read.maf
    tryCatch({
      maf_obj <- read.maf(maf = file_path, verbose = FALSE)
      
      # Get the data table from MAF object
      maf_data <- maf_obj@data
      
      # Update Tumor_Sample_Barcode with filename
      maf_data$Tumor_Sample_Barcode <- sample_name
      
      # Create new MAF object with updated sample names
      updated_maf <- read.maf(maf = maf_data, verbose = FALSE)
      
      # Store in list
      maf_objects[[sample_name]] <- updated_maf
      
      cat("  - Successfully processed", nrow(maf_data), "mutations\n")
      
    }, error = function(e) {
      cat("  - Error reading", sample_name, ":", e$message, "\n")
    })
  }
  
  if (length(maf_objects) == 0) {
    stop("No valid MAF files could be processed")
  }
  
  cat("\nMerging", length(maf_objects), "MAF objects...\n")
  
  # Extract data tables from all MAF objects
  all_maf_data <- lapply(maf_objects, function(maf_obj) {
    return(maf_obj@data)
  })
  
  # Combine all data tables
  merged_data <- rbindlist(all_maf_data, use.names = TRUE, fill = TRUE)
  
  # Create final merged MAF object
  cat("Creating merged MAF object...\n")
  merged_maf <- read.maf(maf = merged_data, verbose = FALSE)
  
  # Summary statistics
  cat("\nMerge Summary:\n")
  cat("Total mutations:", nrow(merged_maf@data), "\n")
  cat("Total samples:", length(unique(merged_maf@data$Tumor_Sample_Barcode)), "\n")
  cat("Total genes:", length(unique(merged_maf@data$Hugo_Symbol)), "\n")
  
  # Sample distribution
  sample_counts <- table(merged_maf@data$Tumor_Sample_Barcode)
  cat("\nMutations per sample:\n")
  print(sample_counts)
  
  # Variant classification summary
  cat("\nVariant classification summary:\n")
  print(table(merged_maf@data$Variant_Classification))
  
  # Save merged MAF object if output file specified
  if (!is.null(output_file)) {
    write.mafSummary(maf = merged_maf, basename = tools::file_path_sans_ext(output_file))
    fwrite(merged_maf@data, output_file, sep = "\t", quote = FALSE, na = "")
    cat("\nMerged MAF file written to:", output_file, "\n")
  }
  
  return(merged_maf)
}

merged_maf
plotmafSummary(merged_maf)
oncoplot(merged_maf, top = 20)
# Function to create summary plots for merged MAF
create_summary_plots <- function(merged_maf, output_dir = "maf_plots") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat("Creating summary plots...\n")
  
  # MAF summary plot
  pdf(file.path(output_dir, "maf_summary.pdf"), width = 10, height = 8)
  plotmafSummary(maf = merged_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
  dev.off()
  
  # Oncoplot
  pdf(file.path(output_dir, "oncoplot.pdf"), width = 12, height = 8)
  oncoplot(maf = merged_maf, top = 20, removeNonMutated = TRUE)
  dev.off()
  
  # Transition and transversion plot
  pdf(file.path(output_dir, "titv_plot.pdf"), width = 10, height = 6)
  titv_summary <- titv(maf = merged_maf, plot = TRUE, useSyn = TRUE)
  dev.off()
  
  # Variant classification plot
  pdf(file.path(output_dir, "variant_classification.pdf"), width = 10, height = 6)
  plotVafDistribution(maf = merged_maf, vafCol = 'HGVSp')
  dev.off()
  
  cat("Plots saved to:", output_dir, "\n")
  
  return(invisible(NULL))
}

# Function to compare samples in merged MAF
compare_samples <- function(merged_maf, output_dir = "sample_comparisons") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Get sample names
  samples <- getSampleSummary(merged_maf)$Tumor_Sample_Barcode
  
  if (length(samples) < 2) {
    cat("Need at least 2 samples for comparison\n")
    return(invisible(NULL))
  }
  
  cat("Comparing", length(samples), "samples...\n")
  
  # Sample comparison matrix
  pdf(file.path(output_dir, "sample_comparison.pdf"), width = 12, height = 10)
  somaticInteractions(maf = merged_maf, top = 25, pvalue = c(0.05, 0.1))
  dev.off()
  
  # Mutational signature analysis (if possible)
  tryCatch({
    pdf(file.path(output_dir, "mutational_signatures.pdf"), width = 12, height = 8)
    sig_res <- extractSignatures(maf = merged_maf, nTry = 3, plotBestFit = TRUE)
    dev.off()
  }, error = function(e) {
    cat("Could not perform signature analysis:", e$message, "\n")
  })
  
  cat("Sample comparison plots saved to:", output_dir, "\n")
  
  return(invisible(NULL))
}

# Main execution
if (!interactive()) {
  # Get command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) == 0) {
    maf_dir <- "."
    output_file <- "merged_mutations.maf"
  } else if (length(args) == 1) {
    maf_dir <- args[1]
    output_file <- "merged_mutations.maf"
  } else {
    maf_dir <- args[1]
    output_file <- args[2]
  }
  
  # Merge MAF objects
  merged_maf <- merge_maf_objects(maf_dir, output_file)
  
  # Create summary plots
  create_summary_plots(merged_maf)
  
  # Compare samples
  compare_samples(merged_maf)
  
} else {
  cat("Script loaded. Use merge_maf_objects() to merge MAF files.\n")
  cat("Usage: merged_maf <- merge_maf_objects(maf_directory = '.', output_file = 'merged_mutations.maf')\n")
  cat("Then use: create_summary_plots(merged_maf) for visualization\n")
}

# Example usage after merging:
merged_maf <- merge_maf_objects("/Users/denriquez/Library/CloudStorage/OneDrive-KagoshimaUniversity/maf_htlv")
getSampleSummary(merged_maf)
plotmafSummary(merged_maf)
oncoplot(merged_maf, top = 20)
create_summary_plots(merged_maf)
compare_samples(merged_maf)
library(maftools)
PlotOncogenicPathways(merged_maf)
oncoplot(merged_maf, top = 20, genes = c("TP53", "KRAS", "BRAF", "PIK3CA"))
OncogenicPathways(merged_maf)
