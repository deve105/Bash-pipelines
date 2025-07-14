#!/usr/bin/env Rscript

# MAF Files Merger for maftools
# This script merges multiple MAF files and uses filenames as Tumor_Sample_Barcode

library(data.table)
library(dplyr)

# Function to merge MAF files
merge_maf_files <- function(maf_directory = ".", output_file = "merged_mutations.maf") {
  
  # Get all MAF files in the directory
  maf_files <- list.files(path = maf_directory, 
                         pattern = "\\.(maf|MAF)$", 
                         full.names = TRUE)
  
  if (length(maf_files) == 0) {
    stop("No MAF files found in the specified directory")
  }
  
  cat("Found", length(maf_files), "MAF files:\n")
  cat(paste(basename(maf_files), collapse = "\n"), "\n\n")
  
  # Initialize list to store data
  maf_data_list <- list()
  
  # Process each MAF file
  for (i in seq_along(maf_files)) {
    file_path <- maf_files[i]
    
    # Extract filename without extension for sample barcode
    sample_name <- tools::file_path_sans_ext(basename(file_path))
    
    cat("Processing:", sample_name, "\n")
    
    # Read MAF file
    tryCatch({
      maf_data <- fread(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
      
      # Check if it's a valid MAF file (should have required columns)
      required_cols <- c("Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode")
      missing_cols <- setdiff(required_cols, colnames(maf_data))
      
      if (length(missing_cols) > 0) {
        cat("Warning: Missing required columns in", sample_name, ":", paste(missing_cols, collapse = ", "), "\n")
      }
      
      # Set Tumor_Sample_Barcode to filename
      maf_data$Tumor_Sample_Barcode <- sample_name
      
      # Store in list
      maf_data_list[[sample_name]] <- maf_data
      
    }, error = function(e) {
      cat("Error reading", sample_name, ":", e$message, "\n")
    })
  }
  
  if (length(maf_data_list) == 0) {
    stop("No valid MAF files could be processed")
  }
  
  # Combine all MAF data
  cat("\nMerging", length(maf_data_list), "MAF files...\n")
  
  # Get all unique columns across all files
  all_columns <- unique(unlist(lapply(maf_data_list, colnames)))
  
  # Standardize columns in each data frame
  maf_data_list <- lapply(maf_data_list, function(df) {
    missing_cols <- setdiff(all_columns, colnames(df))
    for (col in missing_cols) {
      df[[col]] <- NA
    }
    return(df[, all_columns, with = FALSE])
  })
  
  # Combine all data
  merged_maf <- rbindlist(maf_data_list, use.names = TRUE, fill = TRUE)
  
  # Ensure essential MAF columns are present with appropriate defaults
  essential_cols <- list(
    "Hugo_Symbol" = "Unknown",
    "Entrez_Gene_Id" = 0,
    "Center" = "Unknown",
    "NCBI_Build" = "GRCh38",
    "Chromosome" = "Unknown",
    "Start_Position" = 0,
    "End_Position" = 0,
    "Strand" = "+",
    "Variant_Classification" = "Unknown",
    "Variant_Type" = "Unknown",
    "Reference_Allele" = "N",
    "Tumor_Seq_Allele1" = "N",
    "Tumor_Seq_Allele2" = "N",
    "Tumor_Sample_Barcode" = "Unknown"
  )
  
  # Add missing essential columns with defaults
  for (col in names(essential_cols)) {
    if (!col %in% colnames(merged_maf)) {
      merged_maf[[col]] <- essential_cols[[col]]
    }
  }
  
  # Reorder columns to put essential ones first
  essential_order <- names(essential_cols)
  other_cols <- setdiff(colnames(merged_maf), essential_order)
  final_order <- c(essential_order, other_cols)
  
  merged_maf <- merged_maf[, final_order, with = FALSE]
  
  # Remove rows with missing Hugo_Symbol or Tumor_Sample_Barcode
  merged_maf <- merged_maf[!is.na(Hugo_Symbol) & Hugo_Symbol != "" & 
                          !is.na(Tumor_Sample_Barcode) & Tumor_Sample_Barcode != ""]
  
  # Summary statistics
  cat("\nMerge Summary:\n")
  cat("Total mutations:", nrow(merged_maf), "\n")
  cat("Total samples:", length(unique(merged_maf$Tumor_Sample_Barcode)), "\n")
  cat("Total genes:", length(unique(merged_maf$Hugo_Symbol)), "\n")
  
  # Sample distribution
  sample_counts <- table(merged_maf$Tumor_Sample_Barcode)
  cat("\nMutations per sample:\n")
  print(sample_counts)
  
  # Write merged MAF file
  fwrite(merged_maf, output_file, sep = "\t", quote = FALSE, na = "")
  cat("\nMerged MAF file written to:", output_file, "\n")
  
  return(merged_maf)
}

# Function to validate MAF for maftools
validate_maf_for_maftools <- function(maf_file) {
  
  cat("Validating MAF file for maftools compatibility...\n")
  
  maf_data <- fread(maf_file, sep = "\t", header = TRUE)
  
  # Check required columns for maftools
  required_for_maftools <- c(
    "Hugo_Symbol",
    "Variant_Classification", 
    "Tumor_Sample_Barcode"
  )
  
  missing_required <- setdiff(required_for_maftools, colnames(maf_data))
  
  if (length(missing_required) > 0) {
    cat("ERROR: Missing required columns for maftools:", paste(missing_required, collapse = ", "), "\n")
    return(FALSE)
  }
  
  # Check for valid variant classifications
  valid_classifications <- c(
    "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
    "Missense_Mutation", "Nonsense_Mutation", "Silent", "Splice_Site",
    "Translation_Start_Site", "Nonstop_Mutation", "3'UTR", "3'Flank",
    "5'UTR", "5'Flank", "IGR", "Intron", "RNA", "Targeted_Region"
  )
  
  invalid_classifications <- setdiff(unique(maf_data$Variant_Classification), valid_classifications)
  
  if (length(invalid_classifications) > 0) {
    cat("WARNING: Non-standard variant classifications found:", paste(invalid_classifications, collapse = ", "), "\n")
  }
  
  cat("MAF file validation complete.\n")
  cat("Ready for maftools processing!\n")
  
  return(TRUE)
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
  
  # Merge MAF files
  merged_data <- merge_maf_files(maf_dir, output_file)
  
  # Validate for maftools
  validate_maf_for_maftools(output_file)
  
} else {
  cat("Script loaded. Use merge_maf_files() to merge MAF files.\n")
  cat("Usage: merge_maf_files(maf_directory = '.', output_file = 'merged_mutations.maf')\n")
}

# Example usage for maftools after merging:
# library(maftools)
# merged_maf <- read.maf("merged_mutations.maf")
# plotmafSummary(merged_maf)
# oncoplot(merged_maf, top = 20)