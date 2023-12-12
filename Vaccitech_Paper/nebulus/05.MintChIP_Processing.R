# Load extra RNA functions
source("~/00.setup.R")
home_dir <- "/Genomics/ogtr04/wat2/"

saveDiffBindPeaks <- function(dbObj, fold, output_dir, marker) {
  results_fc_edger <- dba.report(dbObj, method = DBA_EDGER, contrast = 1, fold = fold, bUsePval = TRUE, bNormalized = FALSE)
  results_fc_deseq2 <- dba.report(dbObj, method = DBA_DESEQ2, contrast = 1, fold = fold, bUsePval = TRUE, bNormalized = FALSE)
  consensus_peak_set <- as.data.frame(subset(results_fc_deseq2, names(results_fc_deseq2) %in% intersect(names(results_fc_deseq2), names(results_fc_edger))))
  
  results_fc_edger <- as.data.frame(results_fc_edger)
  results_fc_deseq2 <- as.data.frame(results_fc_deseq2)
  
  if(nrow(results_fc_edger) > 0) {
    results_fc_edger$coordinates <- paste0(results_fc_edger$seqnames, ":", results_fc_edger$start, "-", results_fc_edger$end)
  }
  if(nrow(results_fc_deseq2) > 0) {
    results_fc_deseq2$coordinates <- paste0(results_fc_deseq2$seqnames, ":", results_fc_deseq2$start, "-", results_fc_deseq2$end)
  }
  if(nrow(consensus_peak_set) > 0) {
    consensus_peak_set$coordinates <- paste0(consensus_peak_set$seqnames, ":", consensus_peak_set$start, "-", consensus_peak_set$end)
  }

  write.table(results_fc_edger, file = paste0(output_dir, marker, "_edgeR_FC_", fold, ".tsv"), sep = "\t", quote = FALSE)
  write.table(results_fc_deseq2, file = paste0(output_dir, marker, "_DESeq2_FC_", fold, ".tsv"), sep = "\t", quote = FALSE)
  write.table(consensus_peak_set, file = paste0(output_dir, marker, "_consensus_peak_set_FC_", fold, ".tsv"), sep = "\t", quote = FALSE)
}

saveAllPeaks <- function(dbObj, output_dir, marker) {
  peaks <- dbObj$peaks[[1]]
  peak_set <- data.frame(chr = peaks$Chr, start = peaks$Start, end = peaks$End, strand = "+")
  write.table(peak_set, file = paste0(output_dir, marker, "_all_peaks.tsv"), sep = "\t", quote = FALSE, col.names = FALSE)
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
for(marker in markers) {
  print(marker)
  marker_dir <- paste0(output_dir, marker, "/")
  if (!dir.exists(marker_dir)) {dir.create(marker_dir)}
  samples <- read.table(paste0("~/", marker, "_metadata.tsv"), sep = "\t", header = TRUE)
  # Remove subjects that only have one time point (not both)
  samples <- samples[samples$Tissue  %in% names(table(samples$Tissue)[table(samples$Tissue) == 2]),]
  # Load samples into DiffBind and count peaks
  dbObj <- DiffBind::dba(sampleSheet=samples)
  dbObj$config$cores <- 32
  dbObj <- DiffBind::dba.count(dbObj)
  
  saveAllPeaks(dbObj, output_dir = marker_dir, marker = marker)
  
  # Normalize for both DESeq2 and edgeR
  #dbObj.norm <- dba.normalize(dbObj,normalize=DBA_NORM_NATIVE,
  #                            method=DBA_ALL_METHODS,
  #                            background=TRUE)
  # Find DASs
  #dbObj.norm <- dba.contrast(dbObj.norm, design="~Tissue+Condition")
  #dbObj.norm <- dba.analyze(dbObj.norm, method = DBA_ALL_METHODS, bBlacklist = FALSE, bGreylist = FALSE)
  
  # Save DAS object to file
  #saveRDS(dbObj.norm, file = paste0(marker_dir, marker, "_dbObj.rds"))

  # Save peaks with 0 FC (DESeq2, edgeR, and overlapping)
  #saveDiffBindPeaks(dbObj.norm, fold = 0, output_dir = marker_dir, marker = marker)
  #saveDiffBindPeaks(dbObj.norm, fold = 0.1, output_dir = marker_dir, marker = marker)
  #saveDiffBindPeaks(dbObj.norm, fold = 0.2, output_dir = marker_dir, marker = marker)
  #saveDiffBindPeaks(dbObj.norm, fold = 0.3, output_dir = marker_dir, marker = marker)
  #saveDiffBindPeaks(dbObj.norm, fold = 0.585, output_dir = marker_dir, marker = marker)
  #saveDiffBindPeaks(dbObj.norm, fold = 1, output_dir = marker_dir, marker = marker)
  #saveDiffBindPeaks(dbObj.norm, fold = 2, output_dir = marker_dir, marker = marker)
}
