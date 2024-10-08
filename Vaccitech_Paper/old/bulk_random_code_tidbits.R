find_distribution_of_sc_metadata_for_subjects = function(list_of_subjects, metadata, output_name_prefix=NA) {
  # D1 aliquots (passing scRNA, passing multiome, and failing)
  total_D1_passing_scRNA <- c()
  total_D1_passing_scRNA_sex <- c()
  total_D1_passing_multiome <- c()
  total_D1_passing_multiome_sex <- c()
  total_D1_failing <- c()
  total_D1_failing_sex <- c()
  # D28 aliquots (passing scRNA, passing multiome, and failing)
  total_D28_passing_scRNA <- c()
  total_D28_passing_scRNA_sex <- c()
  total_D28_passing_multiome <- c()
  total_D28_passing_multiome_sex <- c()
  total_D28_failing <- c()
  total_D28_failing_sex <- c()
  # Full pass (D-1 and D28) subjects (passing scRNA and passing multiome)
  total_full_pass_subjects_scRNA_aliquots <- c()
  total_full_pass_subjects_scRNA <- c()
  total_full_pass_subjects_scRNA_sex <- c()
  total_full_pass_subjects_multiome_aliquots <- c()
  total_full_pass_subjects_multiome <- c()
  total_full_pass_subjects_multiome_sex <- c()
  for(subject in list_of_subjects) {
    subject_metadata <- metadata[metadata$subject_id == subject,]
    current_sex <- as.character(subject_metadata$sex[1])
    # Gather information on a subject level
    # Check for full pass for current subject (first for scRNA-seq and then for multiome)
    full_pass = FALSE
    if(sum(subject_metadata$scRNA_seq) == 2 & (sum(subject_metadata$passed_qc_scRNA_seq == "Pass") + sum(subject_metadata$passed_qc_scRNA_seq == "Good")) == 2) {
      total_full_pass_subjects_scRNA <- c(total_full_pass_subjects_scRNA, subject)
      total_full_pass_subjects_scRNA_sex <- c(total_full_pass_subjects_scRNA_sex, current_sex)
      full_pass = TRUE
    } else if(sum(subject_metadata$multiome) == 2) {
      total_full_pass_subjects_multiome <- c(total_full_pass_subjects_multiome, subject)
      total_full_pass_subjects_multiome_sex <- c(total_full_pass_subjects_multiome_sex, current_sex)
      full_pass = TRUE
    }
    # Gather information on an aliquot level
    for(aliquot_id in subject_metadata$aliquot_id) {
      aliquot_subject_metadata <- subject_metadata[subject_metadata$aliquot_id == aliquot_id,]
      # Gather info about D-1
      if(aliquot_subject_metadata$time_point == "2_D_minus_1") {
        if(aliquot_subject_metadata$scRNA_seq == TRUE) {
          if(subject_metadata[subject_metadata$aliquot_id == aliquot_id,]$passed_qc_scRNA_seq  == "Pass" | subject_metadata[subject_metadata$aliquot_id == aliquot_id,]$passed_qc_scRNA_seq  == "Good") {
            total_D1_passing_scRNA <- c(total_D1_passing_scRNA, aliquot_id)
            total_D1_passing_scRNA_sex <- c(total_D1_passing_scRNA_sex, current_sex)
          } else {
            total_D1_failing <- c(total_D1_failing, aliquot_id)
            total_D1_failing_sex <- c(total_D1_failing_sex, current_sex)
          }
          if(full_pass) {
            total_full_pass_subjects_scRNA_aliquots <- c(total_full_pass_subjects_scRNA_aliquots, aliquot_id)
          }
        } else if(aliquot_subject_metadata$multiome == TRUE) {
          total_D1_passing_multiome <- c(total_D1_passing_multiome, aliquot_id)
          total_D1_passing_multiome_sex <- c(total_D1_passing_multiome_sex, current_sex)
          if(full_pass) {
            total_full_pass_subjects_multiome_aliquots <- c(total_full_pass_subjects_multiome_aliquots, aliquot_id)
          }
        }
      # Gather info about D28
      } else if(aliquot_subject_metadata$time_point == "2_D28") {
        if(aliquot_subject_metadata$scRNA_seq == TRUE) {
          if(subject_metadata[subject_metadata$aliquot_id == aliquot_id,]$passed_qc_scRNA_seq  == "Pass" | subject_metadata[subject_metadata$aliquot_id == aliquot_id,]$passed_qc_scRNA_seq  == "Good") {
            total_D28_passing_scRNA <- c(total_D28_passing_scRNA, aliquot_id)
            total_D28_passing_scRNA_sex <- c(total_D28_passing_scRNA_sex, current_sex)
          } else {
            total_D28_failing <- c(total_D28_failing, aliquot_id)
            total_D28_failing_sex <- c(total_D28_failing_sex, current_sex)
          }
          if(full_pass) {
            total_full_pass_subjects_scRNA_aliquots <- c(total_full_pass_subjects_scRNA_aliquots, aliquot_id)
          }
        } else if(aliquot_subject_metadata$multiome == TRUE) {
          total_D28_passing_multiome <- c(total_D28_passing_multiome, aliquot_id)
          total_D28_passing_multiome_sex <- c(total_D28_passing_multiome_sex, current_sex)
          if(full_pass) {
            total_full_pass_subjects_multiome_aliquots <- c(total_full_pass_subjects_multiome_aliquots, aliquot_id)
          }
        }
      }
    }
  }
  # Print report to output
  print("### BREAKDOWN OF SUBJECTS ###")
  print(paste0("Total full pass (D-1 and D28) scRNA subjects: ", length(total_full_pass_subjects_scRNA), " (", sum(total_full_pass_subjects_scRNA_sex == "M"), " Male, ", sum(total_full_pass_subjects_scRNA_sex == "F"), " Female)"))
  print(paste0("Total full pass (D-1 and D28) multiome subjects: ", length(total_full_pass_subjects_multiome), " (", sum(total_full_pass_subjects_multiome_sex == "M"), " Male, ", sum(total_full_pass_subjects_multiome_sex == "F"), " Female)"))
  total_full_pass_subjects <- length(total_full_pass_subjects_scRNA) + length(total_full_pass_subjects_multiome)
  total_full_pass_subjects_male <- sum(total_full_pass_subjects_scRNA_sex == "M") + sum(total_full_pass_subjects_multiome_sex == "M")
  total_full_pass_subjects_female <- sum(total_full_pass_subjects_scRNA_sex == "F") + sum(total_full_pass_subjects_multiome_sex == "F")
  print(paste0("total full pass (D-1 and D28) subjects: ", total_full_pass_subjects, " (", total_full_pass_subjects_male, " Male, ", total_full_pass_subjects_female, " Female)"))
  print("### BREAKDOWN OF ALIQUOTS ###")
  print(paste0("D1 passing scRNA: ", length(total_D1_passing_scRNA), " (", sum(total_D1_passing_scRNA_sex == "M"), " Male, ", sum(total_D1_passing_scRNA_sex == "F"), " Female)"))
  print(paste0("D1 passing true multiome: ", length(total_D1_passing_multiome), " (", sum(total_D1_passing_multiome_sex == "M"), " Male, ", sum(total_D1_passing_multiome_sex == "F"), " Female)"))
  D1_total_passing <- length(total_D1_passing_scRNA) + length(total_D1_passing_multiome)
  D1_total_passing_male <- sum(total_D1_passing_scRNA_sex == "M") + sum(total_D1_passing_multiome_sex == "M")
  D1_total_passing_female <- sum(total_D1_passing_scRNA_sex == "F") + sum(total_D1_passing_multiome_sex == "F")
  print(paste0("total D1 passing: ", D1_total_passing, " (", D1_total_passing_male, " Male, ", D1_total_passing_female, " Female)"))
  print(paste0("total D1 failing scRNA: ", length(total_D1_failing), " (", sum(total_D1_failing_sex == "M"), " Male, ", sum(total_D1_failing_sex == "F"), " Female)"))
  print(paste0("D28 passing scRNA: ", length(total_D28_passing_scRNA), " (", sum(total_D28_passing_scRNA_sex == "M"), " Male, ", sum(total_D28_passing_scRNA_sex == "F"), " Female)"))
  print(paste0("D28 passing true multiome: ", length(total_D28_passing_multiome), " (", sum(total_D28_passing_multiome_sex == "M"), " Male, ", sum(total_D28_passing_multiome_sex == "F"), " Female)"))
  D28_total_passing <- length(total_D28_passing_scRNA) + length(total_D28_passing_multiome)
  D28_total_passing_male <- sum(total_D28_passing_scRNA_sex == "M") + sum(total_D28_passing_multiome_sex == "M")
  D28_total_passing_female <- sum(total_D28_passing_scRNA_sex == "F") + sum(total_D28_passing_multiome_sex == "F")
  print(paste0("total D28 passing: ", D28_total_passing, " (", D28_total_passing_male, " Male, ", D28_total_passing_female, " Female)"))
  print(paste0("total D28 failing scRNA: ", length(total_D28_failing), " (", sum(total_D28_failing_sex == "M"), " Male, ", sum(total_D28_failing_sex == "F"), " Female)"))
  # Print aliquots to file
  aliquot_output <- ""
  aliquot_output <- paste0(aliquot_output, "########\n")
  aliquot_output <- paste0(aliquot_output, "FULL PASS (D-1 AND D28) scRNA SUBJECTS:\n")
  aliquot_index <- 1
  for(index in 1:length(total_full_pass_subjects_scRNA)) {
    first_aliquot <- total_full_pass_subjects_scRNA_aliquots[aliquot_index]
    second_aliquot <- total_full_pass_subjects_scRNA_aliquots[aliquot_index + 1]
    first_aliquot_time <- metadata[metadata$aliquot_id == first_aliquot,]$time_point
    second_aliquot_time <- metadata[metadata$aliquot_id == second_aliquot,]$time_point
    if(first_aliquot_time == "2_D_minus_1") {
      first_aliquot_time <- "D-1"
      second_aliquot_time <- "D28"
    } else if (first_aliquot_time == "2_D28") {
      # Swap aliquots because we want D-1 listed first and D28 listed second
      temp <- first_aliquot
      first_aliquot <- second_aliquot
      second_aliquot <- temp
      first_aliquot_time <- "D-1"
      second_aliquot_time <- "D28"
    }
    aliquot_output <- paste0(aliquot_output, total_full_pass_subjects_scRNA[index], " (" , total_full_pass_subjects_scRNA_sex[index], "): ", first_aliquot, " (", first_aliquot_time, ") and ", second_aliquot, " (", second_aliquot_time, ")\n")
    aliquot_index <- aliquot_index + 2
  }
  aliquot_output <- paste0(aliquot_output, "########\n")
  aliquot_output <- paste0(aliquot_output, "FULL PASS (D-1 AND D28) MULTIOME SUBJECTS:\n")
  aliquot_index <- 1
  for(index in 1:length(total_full_pass_subjects_multiome)) {
    first_aliquot <- total_full_pass_subjects_multiome_aliquots[aliquot_index]
    second_aliquot <- total_full_pass_subjects_multiome_aliquots[aliquot_index + 1]
    first_aliquot_time <- metadata[metadata$aliquot_id == first_aliquot,]$time_point
    second_aliquot_time <- metadata[metadata$aliquot_id == second_aliquot,]$time_point
    if(first_aliquot_time == "2_D_minus_1") {
      first_aliquot_time <- "D-1"
      second_aliquot_time <- "D28"
    } else if (first_aliquot_time == "2_D28") {
      # Swap aliquots because we want D-1 listed first and D28 listed second
      temp <- first_aliquot
      first_aliquot <- second_aliquot
      second_aliquot <- temp
      first_aliquot_time <- "D-1"
      second_aliquot_time <- "D28"
    }
    aliquot_output <- paste0(aliquot_output, total_full_pass_subjects_multiome[index], " (" , total_full_pass_subjects_multiome_sex[index], "): ", first_aliquot, " (", first_aliquot_time, ") and ", second_aliquot, " (", second_aliquot_time, ")\n")
    aliquot_index <- aliquot_index + 2
  }
  aliquot_output <- paste0(aliquot_output, "########\n")
  aliquot_output <- paste0(aliquot_output, "D-1 PASSING scRNA ALIQUOTS:\n")
  for(index in 1:length(total_D1_passing_scRNA)) {
    aliquot_output <- paste0(aliquot_output, total_D1_passing_scRNA[index], " (" , total_D1_passing_scRNA_sex[index], ")\n")
  }
  aliquot_output <- paste0(aliquot_output, "########\n")
  aliquot_output <- paste0(aliquot_output, "D-1 PASSING MULTIOME ALIQUOTS:\n")
  for(index in 1:length(total_D1_passing_multiome)) {
    aliquot_output <- paste0(aliquot_output, total_D1_passing_multiome[index], " (" , total_D1_passing_multiome_sex[index], ")\n")
  }
  aliquot_output <- paste0(aliquot_output, "########\n")
  aliquot_output <- paste0(aliquot_output, "D28 PASSING scRNA ALIQUOTS:\n")
  for(index in 1:length(total_D28_passing_scRNA)) {
    aliquot_output <- paste0(aliquot_output, total_D28_passing_scRNA[index], " (" , total_D28_passing_scRNA_sex[index], ")\n")
  }
  aliquot_output <- paste0(aliquot_output, "########\n")
  aliquot_output <- paste0(aliquot_output, "D28 PASSING MULTIOME ALIQUOTS:\n")
  for(index in 1:length(total_D28_passing_multiome)) {
    aliquot_output <- paste0(aliquot_output, total_D28_passing_multiome[index], " (" , total_D28_passing_multiome_sex[index], ")\n")
  }
  aliquot_output <- paste0(aliquot_output, "########")
  if(!is.na(output_name_prefix)) {
    writeLines(aliquot_output, paste0(data_dir, output_name_prefix, "_aliquot_distribution.txt"))
  } else {
    writeLines(aliquot_output, paste0(data_dir, "aliquot_distribution.txt"))
  }
}


find_distribution_of_sc_metadata_for_subjects(high_viral_load_subjects, all_metadata, "high")
find_distribution_of_sc_metadata_for_subjects(low_viral_load_subjects, all_metadata, "low")




# Is one placebo subject in the high viral load subjects creating a lot of DEGs between Period 2 D-2 and D-1?
one_removed_placebo_runs <- list()
for(current_subject in unique(high_placebo_metadata$subject_id)) {
  subset_metadata <- high_placebo_metadata[high_placebo_metadata$subject_id != current_subject,]
  subset_aliquots <- rownames(subset_metadata)
  subset_counts <- high_placebo_counts[,subset_aliquots]
  subset_results <- run_deseq_bulk_analysis("placebo", subset_counts, subset_metadata,
                          "2_D_minus_1", "2_D_minus_2", data_dir, "test")
  one_removed_placebo_runs <- list(one_removed_placebo_runs, subset_results)
}