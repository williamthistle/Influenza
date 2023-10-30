# Load extra RNA functions
source("~/00.setup.R")
home_dir <- "/Genomics/ogtr04/wat2/"

# Helper function to set column to include chr and start/end coordinates in format for UCSC liftover
set_chr_info <- function(differential_analysis_results) {
  chr_loc_info <- strsplit(rownames(differential_analysis_results), "_")
  chr <- sapply(chr_loc_info, function(x) x[1])
  start <- sapply(chr_loc_info, function(x) x[2])
  end <- sapply(chr_loc_info, function(x) x[3])
  differential_analysis_results$chr <- chr
  differential_analysis_results$start <- start
  differential_analysis_results$end <- end
  differential_analysis_results$coordinates <- paste0(chr, ":", start, "-", end)
  new_col_order <- c("chr", "start", "end", "coordinates", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  differential_analysis_results <- differential_analysis_results[new_col_order]
  return(differential_analysis_results)
}


################## SETUP ##################
data_path <- paste0(home_dir, "mintchip/data")
processed_data_path <- paste0(home_dir, "mintchip/processed_data/")
processed_normalized_data_path <- paste0(home_dir, "mintchip/processed_data_normalized/")
output_dir <- paste0(home_dir, "mintchip/analysis/")
analysis_name <- "primary_analysis_11_subject_22_sample"
# data_token is used to choose subset of data that we want to analyze (pre-defined in flu_data_tokens.tsv)
data_token <- "mintchip_paired_sample"
# Grab samples that we want to analyze
data_tokens <- read.table("~/flu_data_tokens.tsv", header = TRUE)
all_sample_id_list <- list.dirs(data_path, recursive = FALSE)
all_sample_id_list <- strsplit(all_sample_id_list, "/")
all_sample_id_list <- unlist(lapply(all_sample_id_list, tail, n = 1L))
samples <- data_tokens[data_tokens$token == data_token,]$samples
samples <-  unlist(strsplit(samples, ","))
sample_id_list <- all_sample_id_list[all_sample_id_list %in% samples]

# Load information about samples
sample_metadata <- read.table("~/all_metadata_sheet.tsv", sep = "\t", header = TRUE)
sample_metadata <- sample_metadata[sample_metadata$aliquot_id %in% sample_id_list,]
sample_metadata_for_SPEEDI_df <- sample_metadata
rownames(sample_metadata_for_SPEEDI_df) <- sample_metadata_for_SPEEDI_df$aliquot_id
sample_metadata_for_SPEEDI_df <- sample_metadata_for_SPEEDI_df[c("subject_id", "time_point", "sex", "viral_load")]
sample_metadata_for_SPEEDI_df$time_point[sample_metadata_for_SPEEDI_df$time_point == 'D-1'] <- 'D_minus_1'

# Break down metadata by category
high_viral_load_samples <- sort(rownames(sample_metadata_for_SPEEDI_df[sample_metadata_for_SPEEDI_df$viral_load == "high",]))
low_viral_load_samples <- sort(rownames(sample_metadata_for_SPEEDI_df[sample_metadata_for_SPEEDI_df$viral_load == "low",]))
all_viral_load_samples <- c(high_viral_load_samples, low_viral_load_samples)
d28_samples <- sort(rownames(sample_metadata_for_SPEEDI_df[sample_metadata_for_SPEEDI_df$time_point == "D28",]))
d_minus_1_samples <- sort(rownames(sample_metadata_for_SPEEDI_df[sample_metadata_for_SPEEDI_df$time_point == "D_minus_1",]))
all_day_samples <- c(d28_samples, d_minus_1_samples)
male_samples <- sort(rownames(sample_metadata_for_SPEEDI_df[sample_metadata_for_SPEEDI_df$sex == "M",]))
female_samples <- sort(rownames(sample_metadata_for_SPEEDI_df[sample_metadata_for_SPEEDI_df$sex == "F",]))
all_sex_samples <- c(male_samples, female_samples)

# Normalize paths (in case user provides relative paths)
data_path <- normalize_dir_path(data_path)
reference_dir <- normalize_dir_path(reference_dir)
output_dir <- normalize_dir_path(output_dir)
# Create output_dir if it doesn't already exist
if (!dir.exists(output_dir)) {dir.create(output_dir)}
# Add "/" to end of output_dir if not already present
last_char_of_output_dir_path <- substr(output_dir, nchar(output_dir), nchar(output_dir))
if(last_char_of_output_dir_path != "/") {
  output_dir <- paste0(output_dir, "/")
}
# Set analysis name
if(is.null(analysis_name)) {
  analysis_name <- paste0(gsub(" ", "_", Sys.time()), "_SPEEDI")
  analysis_name <- gsub(":", "-", analysis_name)
}
# Update our output dir to be the specific analysis directory
output_dir <- paste0(output_dir, analysis_name, "/")
if (!dir.exists(output_dir)) {dir.create(output_dir)}
setwd(output_dir)

# Markers from mintchip
markers <- c("H3K4me1", "H3K4me3", "H3K9me3", "H3K27Ac", "H3K27me3", "H3K36me3")
DAS_matrices <- list()
for(marker in markers) {
  print(marker)
  samples <- read.table(paste0("~/", marker, "_metadata.tsv"), sep = "\t", header = TRUE)
  # Remove subjects that only have one time point (not both)
  samples <- samples[samples$Tissue  %in% names(table(samples$Tissue)[table(samples$Tissue) == 2]),]
  current_peak_info <- DiffBind::dba(sampleSheet=samples)
  current_peak_info <- DiffBind::dba.count(current_peak_info)
  current_peak_info$config$cores <- 32
  current_peaks <- current_peak_info$peaks
  first_peaks <- current_peaks[[1]]
  peak_row_names <- paste0(first_peaks$Chr, "_", first_peaks$Start, "_", first_peaks$End)
  peak_col_names <- current_peak_info$samples$SampleID
  
  # Create a matrix filled with zeroes using nrow and ncol
  numRows <- length(peak_row_names)
  numCols <- length(peak_col_names)
  zeroMatrix <- matrix(0, nrow = numRows, ncol = numCols)
  
  # Create a data frame with the zero matrix and set row and column names
  peak_counts <- data.frame(zeroMatrix)
  rownames(peak_counts) <- peak_row_names
  colnames(peak_counts) <- peak_col_names
  
  for(current_index in 1:length(current_peaks)) {
    peak_counts[,current_index] <- current_peaks[[current_index]]$Reads
  }
  
  rownames(samples) <- samples$SampleID
  
  # Note that Tissue is actually subject ID
  differential_analysis <- DESeq2::DESeqDataSetFromMatrix(countData = peak_counts, colData = samples, design = stats::formula("~ Tissue + Condition"))
  differential_analysis <- DESeq2::DESeq(differential_analysis)
  save(differential_analysis, file = paste0(output_dir, marker, "_", "D28_D1_diff_obj.rds"))
  differential_analysis_results_contrast <- utils::tail(DESeq2::resultsNames(differential_analysis), n=1)
  print(differential_analysis_results_contrast)
  
  # LFC = 0
  differential_analysis_results_0 <- DESeq2::results(differential_analysis, name=differential_analysis_results_contrast)
  differential_analysis_results_0 <- set_chr_info(differential_analysis_results_0)
  
  write.table(differential_analysis_results_0, paste0(output_dir, marker, "_", "D28_D1_diff_unfiltered.tsv"), quote = FALSE, sep = "\t")
  
  differential_analysis_results_0 <- differential_analysis_results_0[rowSums(is.na(differential_analysis_results_0)) == 0, ] # Remove NAs
  differential_analysis_results_0_filtered <- differential_analysis_results_0[differential_analysis_results_0$pvalue < 0.05,]
  
  write.table(differential_analysis_results_0_filtered, paste0(output_dir, marker, "_", "D28_D1_diff_filtered.tsv"), quote = FALSE, sep = "\t")
  
  # LFC = 0.1
  differential_analysis_results_0.1 <- DESeq2::results(differential_analysis, name=differential_analysis_results_contrast, lfcThreshold = 0.1)
  differential_analysis_results_0.1 <- set_chr_info(differential_analysis_results_0.1)
  differential_analysis_results_0.1 <- differential_analysis_results_0.1[rowSums(is.na(differential_analysis_results_0.1)) == 0, ] # Remove NAs
  differential_analysis_results_0.1_filtered <- differential_analysis_results_0.1[differential_analysis_results_0.1$pvalue < 0.05,]
  
  write.table(differential_analysis_results_0.1_filtered, paste0(output_dir, marker, "_", "D28_D1_diff_filtered_0.1.tsv"), quote = FALSE, sep = "\t")
  
  # LFC = 0.2
  differential_analysis_results_0.2 <- DESeq2::results(differential_analysis, name=differential_analysis_results_contrast, lfcThreshold = 0.2)
  differential_analysis_results_0.2 <- set_chr_info(differential_analysis_results_0.2)
  differential_analysis_results_0.2 <- differential_analysis_results_0.2[rowSums(is.na(differential_analysis_results_0.2)) == 0, ] # Remove NAs
  differential_analysis_results_0.2_filtered <- differential_analysis_results_0.2[differential_analysis_results_0.2$pvalue < 0.05,]
  
  write.table(differential_analysis_results_0.2_filtered, paste0(output_dir, marker, "_", "D28_D1_diff_filtered_0.2.tsv"), quote = FALSE, sep = "\t")
  
  # LFC = 0.3
  differential_analysis_results_0.3 <- DESeq2::results(differential_analysis, name=differential_analysis_results_contrast, lfcThreshold = 0.3)
  differential_analysis_results_0.3 <- set_chr_info(differential_analysis_results_0.3)
  differential_analysis_results_0.3 <- differential_analysis_results_0.3[rowSums(is.na(differential_analysis_results_0.3)) == 0, ] # Remove NAs
  differential_analysis_results_0.3_filtered <- differential_analysis_results_0.3[differential_analysis_results_0.3$pvalue < 0.05,]
  
  write.table(differential_analysis_results_0.3_filtered, paste0(output_dir, marker, "_", "D28_D1_diff_filtered_0.3.tsv"), quote = FALSE, sep = "\t")
  
  # LFC = 0.585
  differential_analysis_results_0.585 <- DESeq2::results(differential_analysis, name=differential_analysis_results_contrast, lfcThreshold = 0.585)
  differential_analysis_results_0.585 <- set_chr_info(differential_analysis_results_0.585)
  differential_analysis_results_0.585 <- differential_analysis_results_0.585[rowSums(is.na(differential_analysis_results_0.585)) == 0, ] # Remove NAs
  differential_analysis_results_0.585_filtered <- differential_analysis_results_0.585[differential_analysis_results_0.585$pvalue < 0.05,]
  
  write.table(differential_analysis_results_0.585_filtered, paste0(output_dir, marker, "_", "D28_D1_diff_filtered_0.585.tsv"), quote = FALSE, sep = "\t")
  
  # LFC = 1
  differential_analysis_results_1 <- DESeq2::results(differential_analysis, name=differential_analysis_results_contrast, lfcThreshold = 1)
  differential_analysis_results_1 <- set_chr_info(differential_analysis_results_1)
  differential_analysis_results_1 <- differential_analysis_results_1[rowSums(is.na(differential_analysis_results_1)) == 0, ] # Remove NAs
  differential_analysis_results_1_filtered <- differential_analysis_results_1[differential_analysis_results_1$pvalue < 0.05,]
  
  write.table(differential_analysis_results_1_filtered, paste0(output_dir, marker, "_", "D28_D1_diff_filtered_1.tsv"), quote = FALSE, sep = "\t")
  
  # LFC = 2
  differential_analysis_results_2 <- DESeq2::results(differential_analysis, name=differential_analysis_results_contrast, lfcThreshold = 2)
  differential_analysis_results_2 <- set_chr_info(differential_analysis_results_2)
  differential_analysis_results_2 <- differential_analysis_results_2[rowSums(is.na(differential_analysis_results_2)) == 0, ] # Remove NAs
  differential_analysis_results_2_filtered <- differential_analysis_results_2[differential_analysis_results_2$pvalue < 0.05,]
  
  write.table(differential_analysis_results_2_filtered, paste0(output_dir, marker, "_", "D28_D1_diff_filtered_2.tsv"), quote = FALSE, sep = "\t")
  DAS_matrices[[marker]] <- list(differential_analysis_results_0_filtered, differential_analysis_results_0.1_filtered, 
                                 differential_analysis_results_0.585_filtered, differential_analysis_results_1_filtered, 
                                 differential_analysis_results_2_filtered)
}
