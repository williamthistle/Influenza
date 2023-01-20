high_viral_load_single_cell <- all_metadata
high_viral_load_single_cell <- high_viral_load_single_cell[0,]
# Count number of assays for each subject
for(subject in high_viral_load_subjects ) {
  subset_metadata <- all_metadata[all_metadata$subject_id == subject,]
  subset_metadata <- subset_metadata[subset_metadata$scRNA_seq == TRUE | subset_metadata$scATAC_seq == TRUE | subset_metadata$multiome == TRUE,]
  high_viral_load_single_cell <- rbind(high_viral_load_single_cell, subset_metadata)
}

high_total_D1_passing_scRNA <- c()
high_total_D1_failing <- c()
high_total_D1_passing_true_multiome <- c()
high_total_D28_passing_scRNA <- c()
high_total_D28_failing <- c()
high_total_D28_passing_true_multiome <- c()
for(subject in unique(high_viral_load_single_cell$subject_id)) {
  print(paste0("Current subject: ", subject))
  subset_metadata <- high_viral_load_single_cell[high_viral_load_single_cell$subject_id == subject,]
  for(aliquot_id in subset_metadata$aliquot_id) {
    aliquot_subset_metadata <- subset_metadata[subset_metadata$aliquot_id == aliquot_id,]
    print(paste0("Current aliquot: ", aliquot_id))
    if(aliquot_subset_metadata$time_point == "2_D_minus_1") {
      if(aliquot_subset_metadata$scRNA_seq == TRUE) {
        if(subset_metadata[subset_metadata$aliquot_id == aliquot_id,]$passed_qc_scRNA_seq  == "Pass" | subset_metadata[subset_metadata$aliquot_id == aliquot_id,]$passed_qc_scRNA_seq  == "Good") {
          high_total_D1_passing_scRNA <- c(high_total_D1_passing_scRNA, aliquot_id)
        } else {
          high_total_D1_failing <- c(high_total_D1_failing, aliquot_id)
        }
      } else if(aliquot_subset_metadata$multiome == TRUE) {
        high_total_D1_passing_true_multiome <- c(high_total_D1_passing_true_multiome, aliquot_id)
      }
    } else if(aliquot_subset_metadata$time_point == "2_D28") {
      if(aliquot_subset_metadata$scRNA_seq == TRUE) {
        if(subset_metadata[subset_metadata$aliquot_id == aliquot_id,]$passed_qc_scRNA_seq  == "Pass" | subset_metadata[subset_metadata$aliquot_id == aliquot_id,]$passed_qc_scRNA_seq  == "Good") {
          high_total_D28_passing_scRNA <- c(high_total_D28_passing_scRNA, aliquot_id)
        } else {
          high_total_D28_failing <- c(high_total_D28_failing, aliquot_id)
        }
      } else if(aliquot_subset_metadata$multiome == TRUE) {
        high_total_D28_passing_true_multiome <- c(high_total_D28_passing_true_multiome, aliquot_id)
      }
    }
  }
}

print(paste0("D1 passing scRNA: ", length(high_total_D1_passing_scRNA)))
print(paste0("D1 passing true multiome: ", length(high_total_D1_passing_true_multiome)))
D1_high_total_passing <- length(high_total_D1_passing_scRNA) + length(high_total_D1_passing_true_multiome)
print(paste0("total D1 passing: ", D1_high_total_passing))
print(paste0("total D1 failing scRNA: ", length(high_total_D1_failing)))

print(paste0("D28 passing scRNA: ", length(high_total_D28_passing_scRNA)))
print(paste0("D28 passing true multiome: ", length(high_total_D28_passing_true_multiome)))
D28_high_total_passing <- length(high_total_D28_passing_scRNA) + length(high_total_D28_passing_true_multiome)
print(paste0("total D28 passing: ", D28_high_total_passing))
print(paste0("total D28 failing scRNA: ", length(high_total_D28_failing)))


low_viral_load_single_cell <- all_metadata
low_viral_load_single_cell <- low_viral_load_single_cell[0,]
# Count number of assays for each subject
for(subject in low_viral_load_subjects ) {
  subset_metadata <- all_metadata[all_metadata$subject_id == subject,]
  subset_metadata <- subset_metadata[subset_metadata$scRNA_seq == TRUE | subset_metadata$scATAC_seq == TRUE | subset_metadata$multiome == TRUE,]
  low_viral_load_single_cell <- rbind(low_viral_load_single_cell, subset_metadata)
}

low_total_D1_passing_scRNA <- c()
low_total_D1_failing <- c()
low_total_D1_passing_true_multiome <- c()
low_total_D28_passing_scRNA <- c()
low_total_D28_failing <- c()
low_total_D28_passing_true_multiome <- c()
for(subject in unique(low_viral_load_single_cell$subject_id)) {
  print(paste0("Current subject: ", subject))
  subset_metadata <- low_viral_load_single_cell[low_viral_load_single_cell$subject_id == subject,]
  for(aliquot_id in subset_metadata$aliquot_id) {
    aliquot_subset_metadata <- subset_metadata[subset_metadata$aliquot_id == aliquot_id,]
    print(paste0("Current aliquot: ", aliquot_id))
    if(aliquot_subset_metadata$time_point == "2_D_minus_1") {
      if(aliquot_subset_metadata$scRNA_seq == TRUE) {
        if(subset_metadata[subset_metadata$aliquot_id == aliquot_id,]$passed_qc_scRNA_seq  == "Pass" | subset_metadata[subset_metadata$aliquot_id == aliquot_id,]$passed_qc_scRNA_seq  == "Good") {
          low_total_D1_passing_scRNA <- c(low_total_D1_passing_scRNA, aliquot_id)
        } else {
          low_total_D1_failing <- c(low_total_D1_failing, aliquot_id)
        }
      } else if(aliquot_subset_metadata$multiome == TRUE) {
        low_total_D1_passing_true_multiome <- c(low_total_D1_passing_true_multiome, aliquot_id)
      }
    } else if(aliquot_subset_metadata$time_point == "2_D28") {
      if(aliquot_subset_metadata$scRNA_seq == TRUE) {
        if(subset_metadata[subset_metadata$aliquot_id == aliquot_id,]$passed_qc_scRNA_seq  == "Pass" | subset_metadata[subset_metadata$aliquot_id == aliquot_id,]$passed_qc_scRNA_seq  == "Good") {
          low_total_D28_passing_scRNA <- c(low_total_D28_passing_scRNA, aliquot_id)
        } else {
          low_total_D28_failing <- c(low_total_D28_failing, aliquot_id)
        }
      } else if(aliquot_subset_metadata$multiome == TRUE) {
        low_total_D28_passing_true_multiome <- c(low_total_D28_passing_true_multiome, aliquot_id)
      }
    }
  }
}

print(paste0("D1 passing scRNA: ", length(low_total_D1_passing_scRNA)))
print(paste0("D1 passing true multiome: ", length(low_total_D1_passing_true_multiome)))
D1_low_total_passing <- length(low_total_D1_passing_scRNA) + length(low_total_D1_passing_true_multiome)
print(paste0("total D1 passing: ", D1_low_total_passing))
print(paste0("total D1 failing scRNA: ", length(low_total_D1_failing)))

print(paste0("D28 passing scRNA: ", length(low_total_D28_passing_scRNA)))
print(paste0("D28 passing true multiome: ", length(low_total_D28_passing_true_multiome)))
D28_low_total_passing <- length(low_total_D28_passing_scRNA) + length(low_total_D28_passing_true_multiome)
print(paste0("total D28 passing: ", D28_low_total_passing))
print(paste0("total D28 failing scRNA: ", length(low_total_D28_failing)))




one_removed_placebo_runs <- list()
for(current_subject in unique(high_placebo_metadata$subject_id)) {
  subset_metadata <- high_placebo_metadata[high_placebo_metadata$subject_id != current_subject,]
  subset_aliquots <- rownames(subset_metadata)
  subset_counts <- high_placebo_counts[,subset_aliquots]
  subset_results <- run_deseq_bulk_analysis("placebo", subset_counts, subset_metadata,
                          "2_D_minus_1", "2_D_minus_2", data_dir, "test")
  one_removed_placebo_runs <- list(one_removed_placebo_runs, subset_results)
}