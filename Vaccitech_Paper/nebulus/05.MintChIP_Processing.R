# Load extra RNA functions
source("~/00.setup.R")
home_dir <- "/Genomics/ogtr04/wat2/"

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


mintchip_data_matrix_file_paths <- list.files(processed_data_path, full.names = TRUE)
DAS_matrices <- list()
for(current_data_file_path in mintchip_data_matrix_file_paths) {
  load(current_data_file_path)
  current_marker <- unique(dbObj.count$samples$Factor)
  # Rename sample IDs to aliquot IDs everywhere
  dbObj.count$class[DBA_ID,] <- mintchip_metadata$Aliquot
  for (i in 1:8) {
    names(dbObj.count$masks[[i]]) <- mintchip_metadata$Aliquot
  }
  colnames(dbObj.count$class) <- mintchip_metadata$Aliquot
  dbObj.count$samples$SampleID <- mintchip_metadata$Aliquot
  dbObj.count$samples$ControlID <- paste0(dbObj.count$samples$SampleID, "_C")
  dbObj.count$class[7,] <- paste0(dbObj.count$samples$SampleID, "_C")
  colnames(dbObj.count$called) <- mintchip_metadata$Aliquot
  colnames(dbObj.count$binding) <- c("CHR", "START", "END", mintchip_metadata$Aliquot)
  # Use tissue field for our subject ID as proxy (maybe should use replicate field?)
  dbObj.count$class[2,] <- mintchip_metadata$Subject
  dbObj.count$samples$Tissue <- mintchip_metadata$Subject
  # Update BAM file names
  for(i in 1:length(dbObj.count$samples$bamReads)) {
    current_bam_token <- strsplit(dbObj.count$samples$bamReads[i], "_")[[1]][2]
    current_bam_token <- strsplit(current_bam_token, "\\.")[[1]][1]
    bam_record <- mintchip_metadata[mintchip_metadata$oldBamToken == current_bam_token,]
    new_bam_token <- bam_record$newBamToken
    dbObj.count$samples$bamReads <- sub(current_bam_token, new_bam_token, dbObj.count$samples$bamReads)
    dbObj.count$samples$Peaks <- sub(current_bam_token, new_bam_token, dbObj.count$samples$Peaks)
  }
  dbObj.count$samples$bamReads <- paste0(mintchip_metadata$bamDir, dbObj.count$samples$bamReads)
  dbObj.count$samples$Peaks <- paste0(mintchip_metadata$bamDir, dbObj.count$samples$Peaks)
  dbObj.count$class[10,] <- dbObj.count$samples$bamReads
  # Grab subjects that have matching pre- and post-exposure and subset DiffBind obj to these subjects
  subject_subset <- table(mintchip_metadata$Subject) == 2
  subject_subset <- names(subject_subset[subject_subset])
  aliquot_subset <- mintchip_metadata[mintchip_metadata$Subject %in% subject_subset,]$Aliquot
  all_aliquots <- dbObj.count$samples$SampleID
  full_subject_aliquot_flag <- all_aliquots %in% aliquot_subset
  dbObj.count$masks$full_subject <- full_subject_aliquot_flag
  names(dbObj.count$masks$full_subject) <- all_aliquots
  dbObj.count_full_subj <- dba(dbObj.count, mask = dbObj.count$masks$full_subject)
  # Finish subsetting
  dbObj.count_full_subj$samples <- dbObj.count_full_subj$samples[dbObj.count_full_subj$samples$SampleID %in% aliquot_subset,]
  rownames(dbObj.count_full_subj$samples) <- NULL
  # Normalize data (using "safest" settings from DiffBind manual)
  dbObj.norm <- dba.normalize(dbObj.count_full_subj,normalize=DBA_NORM_NATIVE,
                              method=DBA_DESEQ2,
                              background=TRUE)
  # Find DASs
  dbObj.norm <- dba.contrast(dbObj.norm,design="~Tissue + Condition")
  dbObj.norm <- dba.analyze(dbObj.norm,bBlacklist = FALSE, bGreylist = FALSE)
  results_fc_0 <- dba.report(dbObj.norm, contrast = 1, fold = 0, bUsePval = TRUE)
  results_fc_0.1 <- dba.report(dbObj.norm, contrast = 1, fold = 0.1)
  results_fc_0.585 <- dba.report(dbObj.norm, contrast = 1, fold = 0.585)
  results_fc_1 <- dba.report(dbObj.norm, contrast = 1, fold = 1)
  results_fc_2 <- dba.report(dbObj.norm, contrast = 1, fold = 2)
  DAS_matrices[[current_marker]] <- list(results_fc_0, results_fc_0.1, results_fc_0.585, 
                                         results_fc_1, results_fc_2)
}
