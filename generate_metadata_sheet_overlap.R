base_dir <- "C:/Users/wat2/Documents/GitHub/Influenza/"
single_cell_
scRNA_data_list <- paste0(base_dir, "scRNA/scRNA_data_list.txt")
scATAC_data_list <- paste0(base_dir, "scATAC/scATAC_data_list.txt")
multiome_data_list <- paste0(base_dir, "multiome/multiome_data_list.txt")
scRNA_qc_file <- paste0(base_dir, "ECHO_FLU_Vaccitech_PBMC_scrnaseq_coded_qc_report_WT.csv")
multiome_qc_file <- paste0(base_dir, "Stanford_FLU_combined_qc_metric_coded_09015022_qc_data.csv")
overall_metadata_file <- paste0(base_dir, "20220609_metadataECHO_Vaccitech_Coded.csv")
# Read in tables
scRNA_data <- read.table(scRNA_data_list)$V1
scRNA_data <- scRNA_data[1:length(scRNA_data) - 1]
scATAC_data <- read.table(scATAC_data_list)$V1
scATAC_data <- scATAC_data[1:length(scATAC_data) - 1]
multiome_data <- read.table(multiome_data_list)$V1
multiome_data <- multiome_data[1:length(multiome_data) - 1]
scRNA_qc <- read.csv(scRNA_qc_file)
multiome_qc <- read.csv(multiome_qc_file)
overall_metadata <- read.csv(overall_metadata_file)
single_cell_all_metadata_sheet_df <- data.frame(aliquot_id = character(), subject_id = character(), scRNA_seq = character(),
                                      scATAC_seq = character(), multiome = character(), has_metadata = character(),
                                      treatment = character(), time_point = character(), sex = character(), 
                                      race = character(), passed_qc_scRNA_seq = character(), passed_qc_multiome = character())
for(scRNA_entry in scRNA_data) {
  # Add aliquot ID to current row
  current_row <- c()
  current_row <- append(current_row, scRNA_entry)
  # Add subject ID to current row
  current_sample = overall_metadata[overall_metadata$X_aliquot_id == scRNA_entry,]
  if(nrow(current_sample) > 0) {
    current_row <- append(current_row, current_sample$SUBJECT_ID)
  } else {
    current_row <- append(current_row, "N/A")
  }
  # Add presence of sCRNA_seq data to current row
  current_row <- append(current_row, TRUE)
  # Add FALSE for scATAC_seq and multiome for now
  current_row <- append(current_row, FALSE)
  current_row <- append(current_row, FALSE)
  if(nrow(current_sample) > 0) {
    current_row <- append(current_row, TRUE)
    current_row <- append(current_row, current_sample$TREATMENT)
    current_row <- append(current_row, current_sample$Time_Point)
    current_row <- append(current_row, current_sample$SEX)
    current_row <- append(current_row, current_sample$RACE)
  } else {
    current_row <- append(current_row, FALSE)
    current_row <- append(current_row, "N/A")
    current_row <- append(current_row, "N/A")
    current_row <- append(current_row, "N/A")
    current_row <- append(current_row, "N/A")
  }
  if(scRNA_entry %in% scRNA_qc$Aliquot_ID) {
    current_qc <- scRNA_qc[scRNA_qc$Aliquot_ID == scRNA_entry,]
    current_row <- append(current_row, current_qc$QC_result)
  } else {
    current_row <- append(current_row, "N/A")
  }
  # Append N/A for multiome QC for now
  current_row <- append(current_row, "N/A")
  single_cell_all_metadata_sheet_df[nrow(single_cell_all_metadata_sheet_df) + 1,] = current_row
}

#scATAC-seq
for(scATAC_entry in scATAC_data) {
  # We've already captured info about this aliquot above
  if(scATAC_entry %in% single_cell_all_metadata_sheet_df$aliquot_id) {
    single_cell_all_metadata_sheet_df[single_cell_all_metadata_sheet_df$aliquot_id == scATAC_entry,]$scATAC_seq <- TRUE
  } else {
    # Add aliquot ID to current row
    current_row <- c()
    current_row <- append(current_row, scATAC_entry)
    # Add subject ID to current row
    current_sample = overall_metadata[overall_metadata$X_aliquot_id == scATAC_entry,]
    if(nrow(current_sample) > 0) {
      current_row <- append(current_row, current_sample$SUBJECT_ID)
    } else {
      current_row <- append(current_row, "N/A")
    }
    # Add FALSE for scRNA_seq (would have caught it above)
    current_row <- append(current_row, FALSE)
    # Add TRUE For sCATAC_seq and FALSE for multiome for now
    current_row <- append(current_row, TRUE)
    current_row <- append(current_row, FALSE)
    if(nrow(current_sample) > 0) {
      current_row <- append(current_row, TRUE)
      current_row <- append(current_row, current_sample$TREATMENT)
      current_row <- append(current_row, current_sample$Time_Point)
      current_row <- append(current_row, current_sample$SEX)
      current_row <- append(current_row, current_sample$RACE)
    } else {
      current_row <- append(current_row, FALSE)
      current_row <- append(current_row, "N/A")
      current_row <- append(current_row, "N/A")
      current_row <- append(current_row, "N/A")
      current_row <- append(current_row, "N/A")
    }
    # Add N/A for scRNA_seq QC
    current_row <- append(current_row, "N/A")
    # Append N/A for multiome QC for now
    current_row <- append(current_row, "N/A")
    single_cell_all_metadata_sheet_df[nrow(single_cell_all_metadata_sheet_df) + 1,] = current_row
  }
}

# multiome
for(multiome_entry in multiome_data) {
  # We've already captured info about this aliquot above
  if(multiome_entry %in% single_cell_all_metadata_sheet_df$aliquot_id) {
    single_cell_all_metadata_sheet_df[single_cell_all_metadata_sheet_df$aliquot_id == multiome_entry,]$multiome <- TRUE
    current_qc <- multiome_qc[multiome_qc$Aliquot_ID == multiome_entry,]$Submission_QC_result
    single_cell_all_metadata_sheet_df[single_cell_all_metadata_sheet_df$aliquot_id == multiome_entry,]$passed_qc_multiome <- current_qc
  } else {
    # Add aliquot ID to current row
    current_row <- c()
    current_row <- append(current_row, multiome_entry)
    # Add subject ID to current row
    current_sample = overall_metadata[overall_metadata$X_aliquot_id == multiome_entry,]
    if(nrow(current_sample) > 0) {
      current_row <- append(current_row, current_sample$SUBJECT_ID)
    } else {
      current_row <- append(current_row, "N/A")
    }
    # Add FALSE for scRNA_seq (would have caught it above)
    current_row <- append(current_row, FALSE)
    # Add FALSE For sCATAC_seq and TRUE for multiome 
    current_row <- append(current_row, FALSE)
    current_row <- append(current_row, TRUE)
    if(nrow(current_sample) > 0) {
      current_row <- append(current_row, TRUE)
      current_row <- append(current_row, current_sample$TREATMENT)
      current_row <- append(current_row, current_sample$Time_Point)
      current_row <- append(current_row, current_sample$SEX)
      current_row <- append(current_row, current_sample$RACE)
    } else {
      current_row <- append(current_row, FALSE)
      current_row <- append(current_row, "N/A")
      current_row <- append(current_row, "N/A")
      current_row <- append(current_row, "N/A")
      current_row <- append(current_row, "N/A")
    }
    # Add N/A for scRNA_seq QC
    current_row <- append(current_row, "N/A")
    if(multiome_entry %in% multiome_qc$Aliquot_ID) {
      current_qc <- multiome_qc[multiome_qc$Aliquot_ID == multiome_entry,]
      current_row <- append(current_row, current_qc$Submission_QC_result)
    } else {
      current_row <- append(current_row, "N/A")
    }
    single_cell_all_metadata_sheet_df[nrow(single_cell_all_metadata_sheet_df) + 1,] = current_row
  }
}

single_cell_all_metadata_sheet_df <- single_cell_all_metadata_sheet_df[order(single_cell_all_metadata_sheet_df$subject_id),]
write.table(single_cell_all_metadata_sheet_df, paste0(base_dir, "single_cell_all_metadata_sheet.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)

single_cell_placebo_metadata_sheet_df <- single_cell_all_metadata_sheet_df[single_cell_all_metadata_sheet_df$treatment == "PLACEBO",]
write.table(single_cell_placebo_metadata_sheet_df, paste0(base_dir, "single_cell_placebo_metadata_sheet.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)

sexes <- c()
visited_subjects <- c()
for (current_row in 1:nrow(single_cell_placebo_metadata_sheet_df)) {
  if(!single_cell_placebo_metadata_sheet_df[current_row,]$subject_id %in% visited_subjects ) {
    sexes <- append(sexes, single_cell_placebo_metadata_sheet_df[current_row,]$sex)
    visited_subjects <- append(visited_subjects, single_cell_placebo_metadata_sheet_df[current_row,]$subject_id)
  }
}

races <- c()
visited_subjects <- c()
for (current_row in 1:nrow(single_cell_placebo_metadata_sheet_df)) {
  if(!single_cell_placebo_metadata_sheet_df[current_row,]$subject_id %in% visited_subjects ) {
    races <- append(races, single_cell_placebo_metadata_sheet_df[current_row,]$race)
    visited_subjects <- append(visited_subjects, single_cell_placebo_metadata_sheet_df[current_row,]$subject_id)
  }
}


print(paste0("The total number of aliquots is: ", nrow(single_cell_all_metadata_sheet_df)))
print(paste0("The total number of placebo aliquots is: ", nrow(single_cell_placebo_metadata_sheet_df)))
print(paste0("The total number of D-1 aliquots is: ", nrow(single_cell_placebo_metadata_sheet_df[single_cell_placebo_metadata_sheet_df$time_point == "D-1",])))
print(paste0("The total number of D28 aliquots is: ", nrow(single_cell_placebo_metadata_sheet_df[single_cell_placebo_metadata_sheet_df$time_point == "D28",])))
print(paste0("The total number of unique placebo subjects is: ", length(unique(single_cell_placebo_metadata_sheet_df$subject_id))))
print(paste0("The total number of male subjects is: ", sum(sexes == "M")))
print(paste0("The total number of female subjects is: ", sum(sexes == "F")))
print(paste0("The total number of white subjects is: ", sum(races == "WHITE")))
print(paste0("The total number of other race subjects is: ", length(races) - sum(races == "WHITE")))

paired_single_cell_placebo_metadata_sheet_df <- single_cell_placebo_metadata_sheet_df[(single_cell_placebo_metadata_sheet_df$scRNA_seq == TRUE & 
                                                                 single_cell_placebo_metadata_sheet_df$scATAC_seq == TRUE) | 
                                                                single_cell_placebo_metadata_sheet_df$multiome == TRUE,]

print(paste0("The total number of subjects with paired scRNA-seq / scATAC-seq data or multiome data (for MAGICAL) is: ", length(unique(paired_single_cell_placebo_metadata_sheet_df$subject_id))))
print(paste0("The total number of subjects with paired scRNA-seq / scATAC-seq data or multiome data (for MAGICAL) and both D-1 and D28 timepoints is: ", length(unique(paired_single_cell_placebo_metadata_sheet_df$subject_id)) - sum(table(paired_single_cell_placebo_metadata_sheet_df$subject_id) < 2)))




