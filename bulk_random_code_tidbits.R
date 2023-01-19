high_viral_load_single_cell <- c()

# Count number of assays for each subject
for(subject in high_viral_load_subjects ) {
  subset_metadata <- all_metadata[all_metadata$subject_id == subject,]
  num_scRNA <- sum(subset_metadata$scRNA_seq)
  num_scATAC <- sum(subset_metadata$scATAC_seq)
  num_multiome <- sum(subset_metadata$multiome)
  if(num_scRNA > 0 | num_scATAC > 0) {
    high_viral_load_single_cell <- c(high_viral_load_single_cell, paste0(num_scRNA, " scRNA and ", num_scATAC, " scATAC"))
  } else if(num_multiome > 0) {
    high_viral_load_single_cell <- c(high_viral_load_single_cell, paste0(num_multiome, " multiome"))
  } else {
    high_viral_load_single_cell <- c(high_viral_load_single_cell, "just bulk")
  }
}
    
names(high_viral_load_single_cell) <- high_viral_load_subjects

one_removed_placebo_runs <- list()
for(current_subject in unique(high_placebo_metadata$subject_id)) {
  subset_metadata <- high_placebo_metadata[high_placebo_metadata$subject_id != current_subject,]
  subset_aliquots <- rownames(subset_metadata)
  subset_counts <- high_placebo_counts[,subset_aliquots]
  subset_results <- run_deseq_bulk_analysis("placebo", subset_counts, subset_metadata,
                          "2_D_minus_1", "2_D_minus_2", data_dir, "test")
  one_removed_placebo_runs <- list(one_removed_placebo_runs, subset_results)
}