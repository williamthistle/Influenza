# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# Add antibody titer info
add_antibody_titer_info <- function(metadata, antibody_titer_data) {
  associated_antibody_titer_data <- antibody_titer_data[antibody_titer_data$subject_id %in% metadata$subject_id,]
  MNT_Baseline_vec <- c()
  MNT_Day_Minus_1_Pre_Challenge_vec <- c()
  MNT_Day_28_Post_Challenge_vec <- c()
  HAI_Day_Minus_1_Pre_Challenge_vec <- c()
  HAI_Day_28_Post_Challenge_vec <- c()
  for(current_row_index in 1:nrow(metadata)) {
    current_row <- metadata[current_row_index,]
    current_subject_id <- current_row$subject_id
    current_antibody_titer_row <- associated_antibody_titer_data[associated_antibody_titer_data$subject_id == current_subject_id,]
    if(nrow(current_antibody_titer_row) == 0) {
      print(current_subject_id)
    }
    MNT_Baseline_vec <- c(MNT_Baseline_vec, current_antibody_titer_row$MNT_Baseline)
    MNT_Day_Minus_1_Pre_Challenge_vec <- c(MNT_Day_Minus_1_Pre_Challenge_vec, current_antibody_titer_row$MNT_Day_Minus_1_Pre_Challenge)
    MNT_Day_28_Post_Challenge_vec <- c(MNT_Day_28_Post_Challenge_vec, current_antibody_titer_row$MNT_Day_28_Post_Challenge)
    HAI_Day_Minus_1_Pre_Challenge_vec <- c(HAI_Day_Minus_1_Pre_Challenge_vec, current_antibody_titer_row$HAI_Day_Minus_1_Pre_Challenge)
    HAI_Day_28_Post_Challenge_vec <- c(HAI_Day_28_Post_Challenge_vec, current_antibody_titer_row$HAI_Day_28_Post_Challenge)
  }
  metadata$MNT_Baseline <- MNT_Baseline_vec
  metadata$MNT_Day_Minus_1_Pre_Challenge <- MNT_Day_Minus_1_Pre_Challenge_vec
  metadata$MNT_Day_28_Post_Challenge <- MNT_Day_28_Post_Challenge_vec
  metadata$HAI_Day_Minus_1_Pre_Challenge <- HAI_Day_Minus_1_Pre_Challenge_vec
  metadata$HAI_Day_28_Post_Challenge <- HAI_Day_28_Post_Challenge_vec
  return(metadata)
}

# High placebo (full time series)
hvl_full_time_series_placebo_second_period_metadata <- hvl_full_time_series_placebo_metadata[hvl_full_time_series_placebo_metadata$period == "2",]
hvl_full_time_series_placebo_second_period_metadata <- hvl_full_time_series_placebo_second_period_metadata[hvl_full_time_series_placebo_second_period_metadata$time_point != "2_D_minus_2",]
hvl_full_time_series_placebo_second_period_metadata <- add_antibody_titer_info(hvl_full_time_series_placebo_second_period_metadata, antibody_titer_data)
hvl_full_time_series_placebo_second_period_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, hvl_full_time_series_placebo_second_period_metadata, contrast = c("2_D_minus_1", "2_D2", "2_D5", "2_D8", "2_D28"))

# High and low placebo (full time series)
# If we subset to just D-1, D8, D28, we see the same trend for D8 prediction but even stronger!
# However, if we do D-1, D8, and D28 for full set of placebo subjects (not just time series), we get insignificant prediction
# This is most likely because we add in the moderate viral load subjects
# But wait! If we look at the full set of HVL and LVL placebo for Day -1, Day 8, Day 28, significance also disappears
# So I think it's a combination of sharp viral load split AND too few data points
both_full_time_series_placebo_second_period_metadata <- both_full_time_series_placebo_metadata[both_full_time_series_placebo_metadata$period == "2",]
both_full_time_series_placebo_second_period_metadata <- both_full_time_series_placebo_second_period_metadata[both_full_time_series_placebo_second_period_metadata$time_point != "2_D_minus_2",]
both_full_time_series_placebo_second_period_metadata <- add_antibody_titer_info(both_full_time_series_placebo_second_period_metadata, antibody_titer_data)
both_placebo_second_period_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, both_full_time_series_placebo_second_period_metadata, contrast = c("2_D_minus_1", "2_D2", "2_D5", "2_D8", "2_D28"))

# Binary D2 vs D-1 for placebo - doesn't work because there are too few observations
matching_all_placebo_2_D2_binary_metadata <- both_full_time_series_placebo_second_period_metadata[both_full_time_series_placebo_second_period_metadata$time_point != "2_D5",]
matching_all_placebo_2_D2_binary_metadata <- matching_all_placebo_2_D2_binary_metadata[matching_all_placebo_2_D2_binary_metadata$time_point != "2_D8",]
matching_all_placebo_2_D2_binary_metadata <- matching_all_placebo_2_D2_binary_metadata[matching_all_placebo_2_D2_binary_metadata$time_point != "2_D28",]
matching_all_placebo_2_D2_binary_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, matching_all_placebo_2_D2_binary_metadata, contrast = c("2_D_minus_1", "2_D2"))

# All subjects (full time series)
all_second_period_metadata <- all_full_time_series_metadata[all_full_time_series_metadata$period == "2",]
all_second_period_metadata <- all_second_period_metadata[all_second_period_metadata$time_point != "2_D_minus_2",]
all_second_period_metadata <- add_antibody_titer_info(all_second_period_metadata, antibody_titer_data)
all_second_period_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, all_second_period_metadata, contrast = c("2_D_minus_1", "2_D2", "2_D5", "2_D8", "2_D28"))

# Binary D2 vs D-1 for all subjects
matching_all_2_D2_binary_metadata <- all_second_period_metadata[all_second_period_metadata$time_point != "2_D5",]
matching_all_2_D2_binary_metadata <- matching_all_2_D2_binary_metadata[matching_all_2_D2_binary_metadata$time_point != "2_D8",]
matching_all_2_D2_binary_metadata <- matching_all_2_D2_binary_metadata[matching_all_2_D2_binary_metadata$time_point != "2_D28",]
matching_all_placebo_2_D2_binary_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, matching_all_2_D2_binary_metadata, contrast = c("2_D_minus_1", "2_D2"))


# High placebo
matching_high_placebo_second_period_metadata <- hvl_placebo_metadata[hvl_placebo_metadata$period == "2",]
matching_high_placebo_second_period_metadata <- matching_high_placebo_second_period_metadata[matching_high_placebo_second_period_metadata$time_point != "2_D_minus_2",]
matching_high_placebo_second_period_metadata <- matching_high_placebo_second_period_metadata[matching_high_placebo_second_period_metadata$time_point != "2_D2",]
matching_high_placebo_second_period_metadata <- matching_high_placebo_second_period_metadata[matching_high_placebo_second_period_metadata$time_point != "2_D5",]
matching_high_placebo_second_period_metadata <- add_antibody_titer_info(matching_high_placebo_second_period_metadata, antibody_titer_data)
matching_high_placebo_second_period_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, matching_high_placebo_second_period_metadata, contrast = c("2_D_minus_1", "2_D8", "2_D28"))

# Low placebo
matching_low_placebo_second_period_metadata <- lvl_placebo_metadata[lvl_placebo_metadata$period == "2",]
matching_low_placebo_second_period_metadata <- matching_low_placebo_second_period_metadata[matching_low_placebo_second_period_metadata$time_point != "2_D_minus_2",]
matching_low_placebo_second_period_metadata <- matching_low_placebo_second_period_metadata[matching_low_placebo_second_period_metadata$time_point != "2_D2",]
matching_low_placebo_second_period_metadata <- matching_low_placebo_second_period_metadata[matching_low_placebo_second_period_metadata$time_point != "2_D5",]
matching_low_placebo_second_period_metadata <- add_antibody_titer_info(matching_low_placebo_second_period_metadata, antibody_titer_data)
matching_low_placebo_second_period_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, matching_low_placebo_second_period_metadata, contrast = c("2_D_minus_1", "2_D8", "2_D28"))

# All placebo (period 2)
matching_all_placebo_second_period_metadata <- placebo_metadata[placebo_metadata$period == "2",]
matching_all_placebo_second_period_metadata <- matching_all_placebo_second_period_metadata[matching_all_placebo_second_period_metadata$time_point != "2_D_minus_2",]
matching_all_placebo_second_period_metadata <- matching_all_placebo_second_period_metadata[matching_all_placebo_second_period_metadata$time_point != "2_D2",]
matching_all_placebo_second_period_metadata <- matching_all_placebo_second_period_metadata[matching_all_placebo_second_period_metadata$time_point != "2_D5",]
matching_all_placebo_second_period_metadata <- add_antibody_titer_info(matching_all_placebo_second_period_metadata, antibody_titer_data)
matching_all_placebo_second_period_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, matching_all_placebo_second_period_metadata, contrast = c("2_D_minus_1", "2_D8", "2_D28"))

# Test HVL and LVL Day -1, Day 8, Day 28 (not full time series)
matching_hvl_and_lvl_placebo_second_period_metadata <- matching_all_placebo_second_period_metadata[matching_all_placebo_second_period_metadata$AVAL <= lvl_threshold | matching_all_placebo_second_period_metadata$AVAL >= hvl_threshold,]
matching_hvl_and_lvl_placebo_second_period_metadata <- add_antibody_titer_info(matching_hvl_and_lvl_placebo_second_period_metadata, antibody_titer_data)
matching_hvl_and_lvl_placebo_second_period_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, matching_hvl_and_lvl_placebo_second_period_metadata, contrast = c("2_D_minus_1", "2_D8", "2_D28"))

# All vaccine (period 2)
matching_all_vaccinated_second_period_metadata <- vaccinated_metadata[vaccinated_metadata$period == "2",]
matching_all_vaccinated_second_period_metadata <- matching_all_vaccinated_second_period_metadata[matching_all_vaccinated_second_period_metadata$time_point != "2_D_minus_2",]
matching_all_vaccinated_second_period_metadata <- matching_all_vaccinated_second_period_metadata[matching_all_vaccinated_second_period_metadata$time_point != "2_D2",]
matching_all_vaccinated_second_period_metadata <- matching_all_vaccinated_second_period_metadata[matching_all_vaccinated_second_period_metadata$time_point != "2_D5",]
matching_all_vaccinated_second_period_metadata <- add_antibody_titer_info(matching_all_vaccinated_second_period_metadata, antibody_titer_data)
matching_all_vaccinated_second_period_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, matching_all_vaccinated_second_period_metadata, contrast = c("2_D_minus_1", "2_D8", "2_D28"))

# All matched samples (period 2)
matching_all_vaccinated_second_period_metadata <- matching_all_vaccinated_second_period_metadata[,-c(51,52,53,54),]
matching_all_samples_second_period_metadata <- rbind(matching_all_placebo_second_period_metadata, matching_all_vaccinated_second_period_metadata)
matching_all_samples_second_period_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, matching_all_samples_second_period_metadata, contrast = c("2_D_minus_1", "2_D8", "2_D28"))

# All matched vaccine (period 1)
matching_all_vaccinated_first_period_metadata <- vaccinated_metadata[vaccinated_metadata$period == "1",]
matching_all_vaccinated_first_period_metadata <- matching_all_vaccinated_first_period_metadata[matching_all_vaccinated_first_period_metadata$time_point != "1_D2",]
matching_all_vaccinated_first_period_metadata <- matching_all_vaccinated_first_period_metadata[matching_all_vaccinated_first_period_metadata$time_point != "1_D28",]
matching_all_vaccinated_first_period_metadata <- add_antibody_titer_info(matching_all_vaccinated_first_period_metadata, antibody_titer_data)
matching_all_vaccinated_first_period_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, matching_all_vaccinated_first_period_metadata, contrast = c("1_D_minus_1", "1_D8"))

# All matched placebo (binary)
matching_all_placebo_2_D28_binary_metadata <- placebo_metadata[placebo_metadata$period == "2",]
matching_all_placebo_2_D28_binary_metadata <- matching_all_placebo_2_D28_binary_metadata[matching_all_placebo_2_D28_binary_metadata$time_point != "2_D_minus_2",]
matching_all_placebo_2_D28_binary_metadata <- matching_all_placebo_2_D28_binary_metadata[matching_all_placebo_2_D28_binary_metadata$time_point != "2_D2",]
matching_all_placebo_2_D28_binary_metadata <- matching_all_placebo_2_D28_binary_metadata[matching_all_placebo_2_D28_binary_metadata$time_point != "2_D5",]
matching_all_placebo_2_D28_binary_metadata <- matching_all_placebo_2_D28_binary_metadata[matching_all_placebo_2_D28_binary_metadata$time_point != "2_D8",]
matching_all_placebo_2_D28_binary_metadata <- add_antibody_titer_info(matching_all_placebo_2_D28_binary_metadata, antibody_titer_data)
matching_all_placebo_2_D28_binary_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, matching_all_placebo_2_D28_binary_metadata, contrast = c("2_D_minus_1", "2_D28"))

matching_all_placebo_2_D8_binary_metadata <- placebo_metadata[placebo_metadata$period == "2",]
matching_all_placebo_2_D8_binary_metadata <- matching_all_placebo_2_D8_binary_metadata[matching_all_placebo_2_D8_binary_metadata$time_point != "2_D_minus_2",]
matching_all_placebo_2_D8_binary_metadata <- matching_all_placebo_2_D8_binary_metadata[matching_all_placebo_2_D8_binary_metadata$time_point != "2_D2",]
matching_all_placebo_2_D8_binary_metadata <- matching_all_placebo_2_D8_binary_metadata[matching_all_placebo_2_D8_binary_metadata$time_point != "2_D5",]
matching_all_placebo_2_D8_binary_metadata <- matching_all_placebo_2_D8_binary_metadata[matching_all_placebo_2_D8_binary_metadata$time_point != "2_D28",]
matching_all_placebo_2_D8_binary_metadata <- add_antibody_titer_info(matching_all_placebo_2_D8_binary_metadata, antibody_titer_data)
matching_all_placebo_2_D8_binary_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, matching_all_placebo_2_D8_binary_metadata, contrast = c("2_D_minus_1", "2_D8"))

matching_all_placebo_2_D_minus_2_binary_metadata <- placebo_metadata[placebo_metadata$period == "2",]
matching_all_placebo_2_D_minus_2_binary_metadata <- matching_all_placebo_2_D_minus_2_binary_metadata[matching_all_placebo_2_D_minus_2_binary_metadata$time_point != "2_D8",]
matching_all_placebo_2_D_minus_2_binary_metadata <- matching_all_placebo_2_D_minus_2_binary_metadata[matching_all_placebo_2_D_minus_2_binary_metadata$time_point != "2_D2",]
matching_all_placebo_2_D_minus_2_binary_metadata <- matching_all_placebo_2_D_minus_2_binary_metadata[matching_all_placebo_2_D_minus_2_binary_metadata$time_point != "2_D5",]
matching_all_placebo_2_D_minus_2_binary_metadata <- matching_all_placebo_2_D_minus_2_binary_metadata[matching_all_placebo_2_D_minus_2_binary_metadata$time_point != "2_D28",]
matching_all_placebo_2_D_minus_2_binary_metadata <- add_antibody_titer_info(matching_all_placebo_2_D_minus_2_binary_metadata, antibody_titer_data)
matching_all_placebo_2_D_minus_2_binary_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, matching_all_placebo_2_D_minus_2_binary_metadata, contrast = c("2_D_minus_1", "2_D_minus_2"))

# All matched vaccinated (binary)
matching_all_vaccinated_2_D28_binary_metadata <- vaccinated_metadata[vaccinated_metadata$period == "2",]
matching_all_vaccinated_2_D28_binary_metadata <- matching_all_vaccinated_2_D28_binary_metadata[matching_all_vaccinated_2_D28_binary_metadata$time_point != "2_D_minus_2",]
matching_all_vaccinated_2_D28_binary_metadata <- matching_all_vaccinated_2_D28_binary_metadata[matching_all_vaccinated_2_D28_binary_metadata$time_point != "2_D2",]
matching_all_vaccinated_2_D28_binary_metadata <- matching_all_vaccinated_2_D28_binary_metadata[matching_all_vaccinated_2_D28_binary_metadata$time_point != "2_D5",]
matching_all_vaccinated_2_D28_binary_metadata <- matching_all_vaccinated_2_D28_binary_metadata[matching_all_vaccinated_2_D28_binary_metadata$time_point != "2_D8",]
matching_all_vaccinated_2_D28_binary_metadata <- add_antibody_titer_info(matching_all_vaccinated_2_D28_binary_metadata, antibody_titer_data)
matching_all_vaccinated_2_D28_binary_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, matching_all_vaccinated_2_D28_binary_metadata, contrast = c("2_D_minus_1", "2_D28"))

matching_all_vaccinated_2_D8_binary_metadata <- vaccinated_metadata[vaccinated_metadata$period == "2",]
matching_all_vaccinated_2_D8_binary_metadata <- matching_all_vaccinated_2_D8_binary_metadata[matching_all_vaccinated_2_D8_binary_metadata$time_point != "2_D_minus_2",]
matching_all_vaccinated_2_D8_binary_metadata <- matching_all_vaccinated_2_D8_binary_metadata[matching_all_vaccinated_2_D8_binary_metadata$time_point != "2_D2",]
matching_all_vaccinated_2_D8_binary_metadata <- matching_all_vaccinated_2_D8_binary_metadata[matching_all_vaccinated_2_D8_binary_metadata$time_point != "2_D5",]
matching_all_vaccinated_2_D8_binary_metadata <- matching_all_vaccinated_2_D8_binary_metadata[matching_all_vaccinated_2_D8_binary_metadata$time_point != "2_D28",]
matching_all_vaccinated_2_D8_binary_metadata <- add_antibody_titer_info(matching_all_vaccinated_2_D8_binary_metadata, antibody_titer_data)
# This is the only one that captures anything reasonable
matching_all_vaccinated_2_D8_binary_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, matching_all_vaccinated_2_D8_binary_metadata, contrast = c("2_D_minus_1", "2_D8"))

# Do subset for testing against placebo - this shows that the signal is pretty random. 
matching_all_vaccinated_2_D8_binary_subjects <- unique(matching_all_vaccinated_2_D8_binary_metadata$subject_id)
matching_all_vaccinated_2_D8_binary_subjects <- sample(matching_all_vaccinated_2_D8_binary_subjects, 45)
matching_all_vaccinated_2_D8_binary_metadata_subset <- matching_all_vaccinated_2_D8_binary_metadata[matching_all_vaccinated_2_D8_binary_metadata$subject_id %in% matching_all_vaccinated_2_D8_binary_subjects,]
matching_all_vaccinated_2_D8_binary_subset_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, matching_all_vaccinated_2_D8_binary_metadata_subset, contrast = c("2_D_minus_1", "2_D8"))


# All samples (binary)
matching_all_samples_2_D28_binary_metadata <- matching_all_samples_second_period_metadata
matching_all_samples_2_D28_binary_metadata <- matching_all_samples_2_D28_binary_metadata[matching_all_samples_2_D28_binary_metadata$time_point != "2_D_minus_2",]
matching_all_samples_2_D28_binary_metadata <- matching_all_samples_2_D28_binary_metadata[matching_all_samples_2_D28_binary_metadata$time_point != "2_D2",]
matching_all_samples_2_D28_binary_metadata <- matching_all_samples_2_D28_binary_metadata[matching_all_samples_2_D28_binary_metadata$time_point != "2_D5",]
matching_all_samples_2_D28_binary_metadata <- matching_all_samples_2_D28_binary_metadata[matching_all_samples_2_D28_binary_metadata$time_point != "2_D8",]
matching_all_samples_2_D28_binary_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, matching_all_samples_2_D28_binary_metadata, contrast = c("2_D_minus_1", "2_D28"))

matching_all_samples_2_D8_binary_metadata <- matching_all_samples_second_period_metadata
matching_all_samples_2_D8_binary_metadata <- matching_all_samples_2_D8_binary_metadata[matching_all_samples_2_D8_binary_metadata$time_point != "2_D_minus_2",]
matching_all_samples_2_D8_binary_metadata <- matching_all_samples_2_D8_binary_metadata[matching_all_samples_2_D8_binary_metadata$time_point != "2_D2",]
matching_all_samples_2_D8_binary_metadata <- matching_all_samples_2_D8_binary_metadata[matching_all_samples_2_D8_binary_metadata$time_point != "2_D5",]
matching_all_samples_2_D8_binary_metadata <- matching_all_samples_2_D8_binary_metadata[matching_all_samples_2_D8_binary_metadata$time_point != "2_D28",]
matching_all_samples_2_D8_binary_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, matching_all_samples_2_D8_binary_metadata, contrast = c("2_D_minus_1", "2_D8"))

matching_all_samples_2_D_minus_2_binary_metadata <- placebo_metadata[placebo_metadata$period == "2",]
matching_all_samples_2_D_minus_2_binary_metadata <- rbind(matching_all_samples_2_D_minus_2_binary_metadata, vaccinated_metadata[vaccinated_metadata$period == "2",])
matching_all_samples_2_D_minus_2_binary_metadata <- matching_all_samples_2_D_minus_2_binary_metadata[matching_all_samples_2_D_minus_2_binary_metadata$time_point != "2_D8",]
matching_all_samples_2_D_minus_2_binary_metadata <- matching_all_samples_2_D_minus_2_binary_metadata[matching_all_samples_2_D_minus_2_binary_metadata$time_point != "2_D2",]
matching_all_samples_2_D_minus_2_binary_metadata <- matching_all_samples_2_D_minus_2_binary_metadata[matching_all_samples_2_D_minus_2_binary_metadata$time_point != "2_D5",]
matching_all_samples_2_D_minus_2_binary_metadata <- matching_all_samples_2_D_minus_2_binary_metadata[matching_all_samples_2_D_minus_2_binary_metadata$time_point != "2_D28",]
matching_all_samples_2_D_minus_2_binary_metadata <- matching_all_samples_2_D_minus_2_binary_metadata[matching_all_samples_2_D_minus_2_binary_metadata$subject_id  %in% names(table(matching_all_samples_2_D_minus_2_binary_metadata$subject_id)[table(matching_all_samples_2_D_minus_2_binary_metadata$subject_id) == 2]),]
matching_all_samples_2_D_minus_2_binary_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, matching_all_samples_2_D_minus_2_binary_metadata, contrast = c("2_D_minus_1", "2_D_minus_2"))






# Other approaches


# HIGH PLACEBO
hvl_placebo_second_period_metadata <- hvl_full_time_series_placebo_metadata[hvl_full_time_series_placebo_metadata$period == "2",]
hvl_placebo_second_period_metadata <- hvl_placebo_second_period_metadata[hvl_placebo_second_period_metadata$time_point != "2_D_minus_2",]
hvl_placebo_second_period_training <- create_multiclassPairs_classifier_for_bulk_data(gene_counts_normalized, hvl_placebo_second_period_metadata)
hvl_placebo_second_period_training_results <- hvl_placebo_second_period_training[[2]]

hvl_placebo_aliquot_d_minus_1_subset <- rownames(hvl_placebo_second_period_metadata[hvl_placebo_second_period_metadata$time_point == "2_D_minus_1",])
hvl_placebo_second_period_training_results_2_D_minus_1 <- hvl_placebo_second_period_training_results[rownames(hvl_placebo_second_period_training_results) %in% hvl_placebo_aliquot_d_minus_1_subset,]
hvl_placebo_2_D_minus_1_metadata <- hvl_placebo_second_period_metadata[hvl_placebo_second_period_metadata$time_point == "2_D_minus_1",]
hvl_placebo_second_period_training_results_2_D_minus_1$AVAL <- hvl_placebo_2_D_minus_1_metadata$AVAL

colnames(hvl_placebo_second_period_training_results_2_D_minus_1) <- c("time_2_D2", "time_2_D8", "time_2_D5", "time_2_D_minus_1", "time_2_D28", "max_score", "tie_flag", "AVAL")

ggplot(data = hvl_placebo_second_period_training_results_2_D_minus_1, mapping = aes(x = time_2_D28, y = AVAL)) +
  geom_point(size = 2) +
  sm_statCorr(corr_method = "spearman") + xlab("Misclassification Prob") + ylab("Viral Load") + labs(title = "Day")

# LOW PLACEBO
lvl_placebo_second_period_metadata <- lvl_placebo_metadata[lvl_placebo_metadata$period == "2",]
# lvl_placebo_second_period_metadata <- lvl_placebo_second_period_metadata[lvl_placebo_second_period_metadata$time_point != "2_D_minus_2",]
lvl_placebo_second_period_training <- create_multiclassPairs_classifier_for_bulk_data(gene_counts_normalized, lvl_placebo_second_period_metadata)
lvl_placebo_second_period_training_results <- lvl_placebo_second_period_training[[2]]

lvl_placebo_second_period_training_results_2_D_minus_1 <- lvl_placebo_second_period_training_results[lvl_placebo_second_period_training_results$max_score == "2_D_minus_1",]
lvl_placebo_2_D_minus_1_metadata <- lvl_placebo_second_period_metadata[lvl_placebo_second_period_metadata$time_point == "2_D_minus_1",]
lvl_placebo_second_period_training_results_2_D_minus_1$AVAL <- lvl_placebo_2_D_minus_1_metadata$AVAL

colnames(lvl_placebo_second_period_training_results_2_D_minus_1) <- c("time_2_D28", "time_2_D8", "time_2_D_minus_1", "time_2_D_minus_2", "time_2_D5", "time_2_D2", "max_score", "tie_flag", "AVAL")

ggplot(data = lvl_placebo_second_period_training_results_2_D_minus_1, mapping = aes(x = time_2_D28, y = AVAL)) +
  geom_point(size = 2) +
  sm_statCorr() + xlab("Misclassification Prob") + ylab("Viral Load") + labs(title = "Day")

# BOTH (HIGH AND LOW) PLACEBO
both_placebo_second_period_metadata <- both_placebo_metadata[both_placebo_metadata$period == "2",]
# both_placebo_second_period_metadata <- both_placebo_second_period_metadata[both_placebo_second_period_metadata$time_point != "2_D_minus_2",]
both_placebo_second_period_training <- create_multiclassPairs_classifier_for_bulk_data(gene_counts_normalized, both_placebo_second_period_metadata)
both_placebo_second_period_training_results <- both_placebo_second_period_training[[2]]

both_placebo_aliquot_d_minus_1_subset <- rownames(both_placebo_second_period_metadata[both_placebo_second_period_metadata$time_point == "2_D_minus_1",])
both_placebo_second_period_training_results_2_D_minus_1 <- both_placebo_second_period_training_results[rownames(both_placebo_second_period_training_results) %in% both_placebo_aliquot_d_minus_1_subset,]
both_placebo_2_D_minus_1_metadata <- both_placebo_second_period_metadata[both_placebo_second_period_metadata$time_point == "2_D_minus_1",]
both_placebo_second_period_training_results_2_D_minus_1$AVAL <- both_placebo_2_D_minus_1_metadata$AVAL

colnames(both_placebo_second_period_training_results_2_D_minus_1) <- c("time_2_D28", "time_2_D8", "time_2_D2", "time_2_D_minus_2", "time_2_D_minus_1", "time_2_D5", "max_score", "tie_flag", "AVAL")

ggplot(data = both_placebo_second_period_training_results_2_D_minus_1, mapping = aes(x = time_2_D5, y = AVAL)) +
  geom_point(size = 2) +
  sm_statCorr() + xlab("Misclassification Score") + ylab("Viral Load") + labs(title = "2_D2")

# ALL PLACEBO
all_placebo_second_period_metadata <- placebo_metadata[placebo_metadata$period == "2",]
all_placebo_second_period_metadata <- all_placebo_second_period_metadata[all_placebo_second_period_metadata$time_point != "2_D_minus_2",]
all_placebo_second_period_training <- create_multiclassPairs_classifier_for_bulk_data(gene_counts_normalized_without_scale, all_placebo_second_period_metadata)
all_placebo_second_period_training_results <- all_placebo_second_period_training[[2]]

all_placebo_aliquot_d_minus_1_subset <- rownames(all_placebo_second_period_metadata[all_placebo_second_period_metadata$time_point == "2_D_minus_1",])
all_placebo_second_period_training_results_2_D_minus_1 <- all_placebo_second_period_training_results[rownames(all_placebo_second_period_training_results) %in% all_placebo_aliquot_d_minus_1_subset,]
all_placebo_2_D_minus_1_metadata <- all_placebo_second_period_metadata[all_placebo_second_period_metadata$time_point == "2_D_minus_1",]
all_placebo_second_period_training_results_2_D_minus_1$AVAL <- all_placebo_2_D_minus_1_metadata$AVAL

colnames(all_placebo_second_period_training_results_2_D_minus_1) <- c("time_2_D_minus_1", "time_2_D28", "time_2_D8", "time_2_D2", "time_2_D5", "max_score", "tie_flag", "AVAL")

ggplot(data = all_placebo_second_period_training_results_2_D_minus_1, mapping = aes(x = time_2_D28, y = AVAL)) +
  geom_point(size = 2) +
  sm_statCorr() + xlab("Misclassification Score") + ylab("Viral Load") + labs(title = "Day")

# ALL PLACEBO BUT MATCHING TIME POINTS
matching_all_placebo_second_period_metadata <- placebo_metadata[placebo_metadata$period == "2",]
matching_all_placebo_second_period_metadata <- matching_all_placebo_second_period_metadata[matching_all_placebo_second_period_metadata$time_point != "2_D_minus_2",]
matching_all_placebo_second_period_metadata <- matching_all_placebo_second_period_metadata[matching_all_placebo_second_period_metadata$time_point != "2_D2",]
matching_all_placebo_second_period_metadata <- matching_all_placebo_second_period_metadata[matching_all_placebo_second_period_metadata$time_point != "2_D5",]
matching_all_placebo_second_period_training <- create_multiclassPairs_classifier_for_bulk_data(gene_counts_normalized, matching_all_placebo_second_period_metadata)
matching_all_placebo_second_period_training_results <- matching_all_placebo_second_period_training[[2]]

matching_all_placebo_aliquot_d_minus_1_subset <- rownames(matching_all_placebo_second_period_metadata[matching_all_placebo_second_period_metadata$time_point == "2_D_minus_1",])
matching_all_placebo_second_period_training_results_2_D_minus_1 <- matching_all_placebo_second_period_training_results[rownames(matching_all_placebo_second_period_training_results) %in% matching_all_placebo_aliquot_d_minus_1_subset,]
matching_all_placebo_2_D_minus_1_metadata <- matching_all_placebo_second_period_metadata[matching_all_placebo_second_period_metadata$time_point == "2_D_minus_1",]
matching_all_placebo_second_period_training_results_2_D_minus_1$AVAL <- matching_all_placebo_2_D_minus_1_metadata$AVAL

colnames(matching_all_placebo_second_period_training_results_2_D_minus_1) <- c("time_2_D_minus_1", "time_2_D28", "time_2_D8", "max_score", "tie_flag", "AVAL")

ggplot(data = matching_all_placebo_second_period_training_results_2_D_minus_1, mapping = aes(x = time_2_D8, y = AVAL)) +
  geom_point(size = 2) +
  sm_statCorr(corr_method = "spearman") + xlab("Misclassification Score") + ylab("Viral Load") + labs(title = "Day")

ggplot(data = matching_all_placebo_second_period_training_results_2_D_minus_1[matching_all_placebo_second_period_training_results_2_D_minus_1$AVAL >= hvl_threshold,], 
       mapping = aes(x = time_2_D28, y = AVAL)) +
  geom_point(size = 2) +
  sm_statCorr(corr_method = "spearman") + xlab("Misclassification Score") + ylab("Viral Load") + labs(title = "Day")

# RANDOM FOREST TEST

# Set up classifier object
matching_all_placebo_second_period_metadata <- placebo_metadata[placebo_metadata$period == "2",]
matching_all_placebo_second_period_metadata <- matching_all_placebo_second_period_metadata[matching_all_placebo_second_period_metadata$time_point != "2_D_minus_2",]
matching_all_placebo_second_period_metadata <- matching_all_placebo_second_period_metadata[matching_all_placebo_second_period_metadata$time_point != "2_D2",]
matching_all_placebo_second_period_metadata <- matching_all_placebo_second_period_metadata[matching_all_placebo_second_period_metadata$time_point != "2_D5",]
matching_all_placebo_second_period_counts <- gene_counts_normalized[,rownames(matching_all_placebo_second_period_metadata)]
rownames(matching_all_placebo_second_period_counts) <- gsub("-", "_", rownames(matching_all_placebo_second_period_counts))
time_points <- matching_all_placebo_second_period_metadata$time_point
classes <- as.vector(unique(time_points))

classifier_object <- ReadData(Data = matching_all_placebo_second_period_counts,
                              Labels = time_points,
                              verbose = TRUE)

# Sort genes
genes_RF <- sort_genes_RF(data_object = classifier_object,
                          # featureNo_altogether, it is better not to specify a number here
                          # featureNo_one_vs_rest, it is better not to specify a number here
                          rank_data = TRUE,
                          platform_wise = FALSE,
                          num.trees = 1000, # more features, more tress are recommended
                          seed=get_speedi_seed(), # for reproducibility
                          verbose = TRUE)
genes_RF # sorted genes object

summary_genes <- summary_genes_RF(sorted_genes_RF = genes_RF,
                                  genes_altogether = c(10,20,50,100,150,200),
                                  genes_one_vs_rest = c(10,20,50,100,150,200))

rules_RF <- sort_rules_RF(data_object = object, 
                          sorted_genes_RF = genes_RF,
                          genes_altogether = 50,
                          genes_one_vs_rest = 50, 
                          num.trees = 1000,# more rules, more tress are recommended 
                          seed=get_speedi_seed(),
                          verbose = TRUE)

parameters <- data.frame(
  gene_repetition=c(3,2,1),
  rules_one_vs_rest=c(2,3,10),
  rules_altogether=c(2,3,10),
  run_boruta=c(FALSE,"make_error",TRUE), # I want to produce error in the 2nd trial
  plot_boruta = FALSE,
  num.trees=c(100,200,300),
  stringsAsFactors = FALSE)

# parameters
# for overall and byclass possible options, check the help files
para_opt <- optimize_RF(data_object = classifier_object,
                        sorted_rules_RF = rules_RF,
                        parameters = parameters,
                        test_object = NULL,
                        overall = c("Accuracy","Kappa"), # wanted overall measurements 
                        byclass = c("F1"), # wanted measurements per class
                        verbose = TRUE,
                        seed = get_speedi_seed())

RF_classifier <- train_RF(data_object = classifier_object,
                          sorted_rules_RF = rules_RF,
                          gene_repetition = 1,
                          rules_altogether = 10,
                          rules_one_vs_rest = 10,
                          run_boruta = TRUE, 
                          plot_boruta = FALSE,
                          probability = TRUE,
                          num.trees = 300,
                          boruta_args = list(),
                          verbose = TRUE)

proximity_matrix_RF(object = classifier_object,
                    classifier = RF_classifier, 
                    plot = TRUE,
                    return_matrix = FALSE, # if we need to get the matrix itself
                    title = "Test",
                    cluster_cols = TRUE)

training_pred <- RF_classifier$RF_scheme$RF_classifier$predictions
if (is.factor(training_pred)) {
  x <- as.character(training_pred)
}

# if the classifier trained using probability   = TRUE
if (is.matrix(training_pred)) {
  x <- colnames(training_pred)[max.col(training_pred)]
}

# training accuracy
caret::confusionMatrix(data =factor(x),
                       reference = factor(classifier_object$data$Labels),
                       mode = "everything")

training_pred <- as.data.frame(training_pred)

training_pred$correct_time_point <- time_points
training_pred$AVAL <- matching_all_placebo_second_period_metadata$AVAL


training_pred_d_minus_1 <- training_pred[training_pred$correct_time_point == "2_D_minus_1",]
colnames(training_pred_d_minus_1) <- c("time_2_D_minus_1", "time_2_D28", "time_2_D8", "correct_time_point", "AVAL")

ggplot(data = training_pred_d_minus_1, mapping = aes(x = time_2_D28, y = AVAL)) +
  geom_point(size = 2) +
  sm_statCorr() + xlab("Misclassification Prob") + ylab("Viral Load") + labs(title = "Day")