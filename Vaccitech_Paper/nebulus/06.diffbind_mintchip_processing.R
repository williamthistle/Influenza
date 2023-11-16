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


# TEST FOR ONE MARKER BELOW

marker <- "H3K4me1"
samples <- read.table("~/H3K4me1_metadata.tsv", sep = "\t", header = TRUE)
# Remove subjects that only have one time point (not both)
samples <- samples[samples$Tissue  %in% names(table(samples$Tissue)[table(samples$Tissue) == 2]),]

dbObj <- dba(sampleSheet=samples)
dbObj$config$cores <- 32
dbObj <- dba.count(dbObj)


dbObj.norm <- dba.normalize(dbObj,normalize=DBA_NORM_NATIVE,
                            method=DBA_ALL_METHODS,
                            background=TRUE)
# Find DASs
dbObj.norm <- dba.contrast(dbObj.norm, design="~Tissue+Condition")
dbObj.norm <- dba.analyze(dbObj.norm, method = DBA_ALL_METHODS, bBlacklist = FALSE, bGreylist = FALSE)

results_fc_0_edger <- dba.report(dbObj.norm, method = DBA_EDGER, contrast = 1, fold = 0, bUsePval = TRUE, bNormalized = FALSE)
results_fc_0_deseq2 <- dba.report(dbObj.norm, method = DBA_DESEQ2, contrast = 1, fold = 0, bUsePval = TRUE, bNormalized = FALSE)
consensus_peak_set_0 <- subset(results_fc_0_deseq2, names(results_fc_0_deseq2) %in% intersect(names(results_fc_0_edger), names(results_fc_0_deseq2)))
write.table(consensus_peak_set_0, file = paste0(output_dir, marker, "_consensus_peak_set_0.tsv"), sep = "\t", quote = FALSE)


results_fc_0.1_edger <- dba.report(dbObj.norm, method = DBA_EDGER, contrast = 1, fold = 0.1, bUsePval = TRUE, bNormalized = FALSE)
results_fc_0.1_deseq2 <- dba.report(dbObj.norm, method = DBA_DESEQ2, contrast = 1, fold = 0.1, bUsePval = TRUE, bNormalized = FALSE)
consensus_peak_set_0.1 <- subset(results_fc_0.1_deseq2, names(results_fc_0.1_deseq2) %in% intersect(names(results_fc_0.1_edger), names(results_fc_0.1_deseq2)))

results_fc_0.2_edger <- dba.report(dbObj.norm, method = DBA_EDGER, contrast = 1, fold = 0.2, bUsePval = TRUE, bNormalized = FALSE)
results_fc_0.2_deseq2 <- dba.report(dbObj.norm, method = DBA_DESEQ2, contrast = 1, fold = 0.2, bUsePval = TRUE, bNormalized = FALSE)
consensus_peak_set_0.2 <- subset(results_fc_0.2_deseq2, names(results_fc_0.2_deseq2) %in% intersect(names(results_fc_0.2_edger), names(results_fc_0.2_deseq2)))

results_fc_0.3_edger <- dba.report(dbObj.norm, method = DBA_EDGER, contrast = 1, fold = 0.3, bUsePval = TRUE, bNormalized = FALSE)
results_fc_0.3_deseq2 <- dba.report(dbObj.norm, method = DBA_DESEQ2, contrast = 1, fold = 0.3, bUsePval = TRUE, bNormalized = FALSE)
consensus_peak_set_0.3 <- subset(results_fc_0.3_deseq2, names(results_fc_0.3_deseq2) %in% intersect(names(results_fc_0.3_edger), names(results_fc_0.3_deseq2)))

