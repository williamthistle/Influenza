# Function to find total number of each sex in dataframe
find_sex_count = function(curr_df) {
  sexes <- c()
  visited_subjects <- c()
  for (current_row in 1:nrow(curr_df)) {
    if(!curr_df[current_row,]$subject_id %in% visited_subjects ) {
      sexes <- append(sexes, curr_df[current_row,]$sex)
      visited_subjects <- append(visited_subjects, curr_df[current_row,]$subject_id)
    }
  }
  return(sexes)
}

# Function to find total number of each race in dataframe
find_race_count = function(curr_df) {
  races <- c()
  visited_subjects <- c()
  for (current_row in 1:nrow(curr_df)) {
    if(!curr_df[current_row,]$subject_id %in% visited_subjects ) {
      races <- append(races, curr_df[current_row,]$race)
      visited_subjects <- append(visited_subjects, curr_df[current_row,]$subject_id)
    }
  }
  return(races)
}

# Function to find total number of each age in dataframe
find_age_count = function(curr_df) {
  ages <- c()
  visited_subjects <- c()
  for (current_row in 1:nrow(curr_df)) {
    if(!curr_df[current_row,]$subject_id %in% visited_subjects ) {
      ages <- append(ages, curr_df[current_row,]$age)
      visited_subjects <- append(visited_subjects, curr_df[current_row,]$subject_id)
    }
  }
  return(ages)
}


base_dir <- "C:/Users/willi/Documents/GitHub/Influenza/"
scRNA_data_list <- paste0(base_dir, "scRNA/scRNA_data_list.txt")
scATAC_data_list <- paste0(base_dir, "scATAC/scATAC_data_list.txt")
multiome_data_list <- paste0(base_dir, "multiome/multiome_data_list.txt")
bulkRNA_data_list <- paste0(base_dir, "bulkRNA/bulkRNA_data_list.txt")
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
bulkRNA_data <- read.table(bulkRNA_data_list)$V1
scRNA_qc <- read.csv(scRNA_qc_file)
multiome_qc <- read.csv(multiome_qc_file)
overall_metadata <- read.csv(overall_metadata_file)
all_metadata_sheet_df <- data.frame(aliquot_id = character(), subject_id = character(), scRNA_seq = character(),
                                      scATAC_seq = character(), multiome = character(), bulkRNA_seq = character(),
                                      has_metadata = character(), specimen_prep = character(), treatment = character(), 
                                      period = character(), time_point = character(), sex = character(), age = character(), 
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
  # Add FALSE for scATAC_seq, multiome, and bulkRNA_seq for now
  current_row <- append(current_row, FALSE)
  current_row <- append(current_row, FALSE)
  current_row <- append(current_row, FALSE)
  if(nrow(current_sample) > 0) {
    current_row <- append(current_row, TRUE)
    current_row <- append(current_row, current_sample$specimen_prep)
    current_row <- append(current_row, current_sample$TREATMENT)
    current_row <- append(current_row, current_sample$Period)
    current_row <- append(current_row, current_sample$Time_Point)
    current_row <- append(current_row, current_sample$SEX)
    current_row <- append(current_row, current_sample$AGE)
    current_row <- append(current_row, current_sample$RACE)
  } else {
    current_row <- append(current_row, FALSE)
    current_row <- append(current_row, "N/A")   
    current_row <- append(current_row, "N/A")
    current_row <- append(current_row, "N/A")
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
  all_metadata_sheet_df[nrow(all_metadata_sheet_df) + 1,] = current_row
}

#scATAC-seq
for(scATAC_entry in scATAC_data) {
  # We've already captured info about this aliquot above
  if(scATAC_entry %in% all_metadata_sheet_df$aliquot_id) {
    all_metadata_sheet_df[all_metadata_sheet_df$aliquot_id == scATAC_entry,]$scATAC_seq <- TRUE
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
    # Add TRUE For sCATAC_seq and FALSE for multiome and bulkRNA_seq for now
    current_row <- append(current_row, TRUE)
    current_row <- append(current_row, FALSE)
    current_row <- append(current_row, FALSE)
    if(nrow(current_sample) > 0) {
      current_row <- append(current_row, TRUE)
      current_row <- append(current_row, current_sample$specimen_prep)
      current_row <- append(current_row, current_sample$TREATMENT)
      current_row <- append(current_row, current_sample$Period)
      current_row <- append(current_row, current_sample$Time_Point)
      current_row <- append(current_row, current_sample$SEX)
      current_row <- append(current_row, current_sample$AGE)
      current_row <- append(current_row, current_sample$RACE)
    } else {
      current_row <- append(current_row, FALSE)
      current_row <- append(current_row, "N/A")
      current_row <- append(current_row, "N/A")
      current_row <- append(current_row, "N/A")
      current_row <- append(current_row, "N/A")
      current_row <- append(current_row, "N/A")
      current_row <- append(current_row, "N/A")
      current_row <- append(current_row, "N/A")
    }
    # Add N/A for scRNA_seq QC
    current_row <- append(current_row, "N/A")
    # Append N/A for multiome QC for now
    current_row <- append(current_row, "N/A")
    all_metadata_sheet_df[nrow(all_metadata_sheet_df) + 1,] = current_row
  }
}

# multiome
for(multiome_entry in multiome_data) {
  # We've already captured info about this aliquot above
  if(multiome_entry %in% all_metadata_sheet_df$aliquot_id) {
    all_metadata_sheet_df[all_metadata_sheet_df$aliquot_id == multiome_entry,]$multiome <- TRUE
    current_qc <- multiome_qc[multiome_qc$Aliquot_ID == multiome_entry,]$Submission_QC_result
    all_metadata_sheet_df[all_metadata_sheet_df$aliquot_id == multiome_entry,]$passed_qc_multiome <- current_qc
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
    # Add FALSE For sCATAC_seq and TRUE for multiome and FALSE for bulkRNA_seq for now
    current_row <- append(current_row, FALSE)
    current_row <- append(current_row, TRUE)
    current_row <- append(current_row, FALSE)
    if(nrow(current_sample) > 0) {
      current_row <- append(current_row, TRUE)
      current_row <- append(current_row, current_sample$specimen_prep)
      current_row <- append(current_row, current_sample$TREATMENT)
      current_row <- append(current_row, current_sample$Period)
      current_row <- append(current_row, current_sample$Time_Point)
      current_row <- append(current_row, current_sample$SEX)
      current_row <- append(current_row, current_sample$AGE)
      current_row <- append(current_row, current_sample$RACE)
    } else {
      current_row <- append(current_row, FALSE)
      current_row <- append(current_row, "N/A")
      current_row <- append(current_row, "N/A")
      current_row <- append(current_row, "N/A")
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
    all_metadata_sheet_df[nrow(all_metadata_sheet_df) + 1,] = current_row
  }
}

# bulkRNA-seq
for(bulkRNA_entry in bulkRNA_data) {
  # We've already captured info about this aliquot above
  if(bulkRNA_entry %in% all_metadata_sheet_df$aliquot_id) {
    all_metadata_sheet_df[all_metadata_sheet_df$aliquot_id == bulkRNA_entry,]$bulkRNA_seq <- TRUE
  } else {
    # Add aliquot ID to current row
    current_row <- c()
    current_row <- append(current_row, bulkRNA_entry)
    # Add subject ID to current row
    current_sample = overall_metadata[overall_metadata$X_aliquot_id == bulkRNA_entry,]
    if(nrow(current_sample) > 0) {
      current_row <- append(current_row, current_sample$SUBJECT_ID)
    } else {
      current_row <- append(current_row, "N/A")
    }
    # Add FALSE for scRNA_seq (would have caught it above)
    current_row <- append(current_row, FALSE)
    # Add FALSE For sCATAC_seq and multiome and TRUE for bulkRNA_seq for now
    current_row <- append(current_row, FALSE)
    current_row <- append(current_row, FALSE)
    current_row <- append(current_row, TRUE)
    if(nrow(current_sample) > 0) {
      current_row <- append(current_row, TRUE)
      current_row <- append(current_row, current_sample$specimen_prep)
      current_row <- append(current_row, current_sample$TREATMENT)
      current_row <- append(current_row, current_sample$Period)
      current_row <- append(current_row, current_sample$Time_Point)
      current_row <- append(current_row, current_sample$SEX)
      current_row <- append(current_row, current_sample$AGE)
      current_row <- append(current_row, current_sample$RACE)
    } else {
      current_row <- append(current_row, FALSE)
      current_row <- append(current_row, "N/A")
      current_row <- append(current_row, "N/A")
      current_row <- append(current_row, "N/A")
      current_row <- append(current_row, "N/A")
      current_row <- append(current_row, "N/A")
      current_row <- append(current_row, "N/A")
      current_row <- append(current_row, "N/A")
    }
    # Add N/A for scRNA_seq QC and multiome QC
    current_row <- append(current_row, "N/A")
    current_row <- append(current_row, "N/A")
    all_metadata_sheet_df[nrow(all_metadata_sheet_df) + 1,] = current_row
  }
}

# Convert age to different categories
# 1: 18-27
# 2: 28-37
# 3: 38-47
# 4: 48-57
all_metadata_sheet_df$age <- as.numeric(all_metadata_sheet_df$age)
all_metadata_sheet_df$age[all_metadata_sheet_df$age >= 18 & all_metadata_sheet_df$age <= 27 & !is.na(all_metadata_sheet_df$age)] <- 1
all_metadata_sheet_df$age[all_metadata_sheet_df$age >= 28 & all_metadata_sheet_df$age <= 37 & !is.na(all_metadata_sheet_df$age)] <- 2
all_metadata_sheet_df$age[all_metadata_sheet_df$age >= 38 & all_metadata_sheet_df$age <= 47 & !is.na(all_metadata_sheet_df$age)] <- 3
all_metadata_sheet_df$age[all_metadata_sheet_df$age >= 48 & all_metadata_sheet_df$age <= 57 & !is.na(all_metadata_sheet_df$age)] <- 4





# Write complete metadata DF to file
all_metadata_sheet_df <- all_metadata_sheet_df[order(all_metadata_sheet_df$subject_id),]
write.table(all_metadata_sheet_df, paste0(base_dir, "all_metadata_sheet.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)

# Write placebo metadata DF to file
placebo_metadata_sheet_df <- all_metadata_sheet_df[all_metadata_sheet_df$treatment == "PLACEBO",]
write.table(placebo_metadata_sheet_df, paste0(base_dir, "placebo_metadata_sheet.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)

# Write vaccinated metadata DF to file
vaccinated_metadata_sheet_df <- all_metadata_sheet_df[all_metadata_sheet_df$treatment == "MVA-NP+M1",]
write.table(vaccinated_metadata_sheet_df, paste0(base_dir, "vaccinated_metadata_sheet.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)




# Print stats
print("*********** OVERALL ***********")
print(paste0("The total number of aliquots is: ", nrow(all_metadata_sheet_df)))
print(paste0("The total number of placebo aliquots is: ", nrow(placebo_metadata_sheet_df)))
print(paste0("The total number of unique placebo subjects is: ", length(unique(placebo_metadata_sheet_df$subject_id))))
print(paste0("The total number of male placebo subjects is: ", sum(find_sex_count(placebo_metadata_sheet_df) == "M")))
print(paste0("The total number of female placebo subjects is: ", sum(find_sex_count(placebo_metadata_sheet_df) == "F")))
print(paste0("The total number of white placebo subjects is: ", sum(find_race_count(placebo_metadata_sheet_df) == "WHITE")))
print(paste0("The total number of other race placebo subjects is: ", length(find_race_count(placebo_metadata_sheet_df)) - sum(find_race_count(placebo_metadata_sheet_df) == "WHITE")))
# Since the vast majority of subjects are white, we don't need to keep tracking race
print("*********** SCRNA-SEQ (PLACEBO) ***********")
scRNAseq_placebo_metadata_sheet_df <- placebo_metadata_sheet_df[placebo_metadata_sheet_df$scRNA_seq == TRUE,]
print(paste0("The total number of aliquots is: ", nrow(scRNAseq_placebo_metadata_sheet_df)))
print(paste0("The total number of D-1 aliquots is: ", nrow(scRNAseq_placebo_metadata_sheet_df[scRNAseq_placebo_metadata_sheet_df$time_point == "D-1",])))
print(paste0("The total number of D28 aliquots is: ", nrow(scRNAseq_placebo_metadata_sheet_df[scRNAseq_placebo_metadata_sheet_df$time_point == "D28",])))
print(paste0("The total number of unique subjects is: ", length(unique(scRNAseq_placebo_metadata_sheet_df$subject_id))))
print(paste0("The total number of male subjects is: ", sum(find_sex_count(scRNAseq_placebo_metadata_sheet_df) == "M")))
print(paste0("The total number of female subjects is: ", sum(find_sex_count(scRNAseq_placebo_metadata_sheet_df) == "F")))
# Note terrible sex ratio (3 male vs 9 female)
scRNAseq_paired_timepoints_placebo_metadata_sheet_df <- scRNAseq_placebo_metadata_sheet_df[scRNAseq_placebo_metadata_sheet_df$subject_id %in% names(table(scRNAseq_placebo_metadata_sheet_df$subject_id)[table(scRNAseq_placebo_metadata_sheet_df$subject_id) == 2]),]
print(paste0("The total number of unique subjects with both D-1 and D28 timepoints is: ", length(unique(scRNAseq_paired_timepoints_placebo_metadata_sheet_df$subject_id))))
print(paste0("The total number of male subjects with both D-1 and D28 timepoints is: ", sum(find_sex_count(scRNAseq_paired_timepoints_placebo_metadata_sheet_df) == "M")))
print(paste0("The total number of female subjects with both D-1 and D28 timepoints is: ", sum(find_sex_count(scRNAseq_paired_timepoints_placebo_metadata_sheet_df) == "F")))
# Note terrible sex ratio (2 male vs 7 female)
print("*********** SCRNA-SEQ (PLACEBO) THAT PASSED QC ***********")
passed_qc_scRNAseq_placebo_metadata_sheet_df <- scRNAseq_placebo_metadata_sheet_df[scRNAseq_placebo_metadata_sheet_df$passed_qc_scRNA_seq  == "Good" | scRNAseq_placebo_metadata_sheet_df$passed_qc_scRNA_seq  == "Pass",]
print(paste0("The total number of aliquots is: ", nrow(passed_qc_scRNAseq_placebo_metadata_sheet_df)))
print(paste0("The total number of D-1 aliquots is: ", nrow(passed_qc_scRNAseq_placebo_metadata_sheet_df[passed_qc_scRNAseq_placebo_metadata_sheet_df$time_point == "D-1",])))
print(paste0("The total number of D28 aliquots is: ", nrow(passed_qc_scRNAseq_placebo_metadata_sheet_df[passed_qc_scRNAseq_placebo_metadata_sheet_df$time_point == "D28",])))
print(paste0("The total number of unique subjects is: ", length(unique(passed_qc_scRNAseq_placebo_metadata_sheet_df$subject_id))))
print(paste0("The total number of male subjects is: ", sum(find_sex_count(passed_qc_scRNAseq_placebo_metadata_sheet_df) == "M")))
print(paste0("The total number of female subjects is: ", sum(find_sex_count(passed_qc_scRNAseq_placebo_metadata_sheet_df) == "F")))
# Note terrible sex ratio (2 male vs 6 female)
scRNAseq_passed_qc_paired_timepoints_placebo_metadata_sheet_df <- passed_qc_scRNAseq_placebo_metadata_sheet_df[passed_qc_scRNAseq_placebo_metadata_sheet_df$subject_id %in% names(table(passed_qc_scRNAseq_placebo_metadata_sheet_df$subject_id)[table(passed_qc_scRNAseq_placebo_metadata_sheet_df$subject_id) == 2]),]
print(paste0("The total number of unique subjects with both D-1 and D28 timepoints is: ", length(unique(scRNAseq_passed_qc_paired_timepoints_placebo_metadata_sheet_df$subject_id))))
print(paste0("The total number of male subjects with both D-1 and D28 timepoints is: ", sum(find_sex_count(scRNAseq_passed_qc_paired_timepoints_placebo_metadata_sheet_df) == "M")))
print(paste0("The total number of female subjects with both D-1 and D28 timepoints is: ", sum(find_sex_count(scRNAseq_passed_qc_paired_timepoints_placebo_metadata_sheet_df) == "F")))
# Note terrible sex ratio (1 male vs 3 female)
print("*********** SCATAC-SEQ (PLACEBO) ***********")
scATACseq_placebo_metadata_sheet_df <- placebo_metadata_sheet_df[placebo_metadata_sheet_df$scATAC_seq == TRUE,]
print(paste0("The total number of aliquots is: ", nrow(scATACseq_placebo_metadata_sheet_df)))
print(paste0("The total number of D-1 aliquots is: ", nrow(scATACseq_placebo_metadata_sheet_df[scATACseq_placebo_metadata_sheet_df$time_point == "D-1",])))
print(paste0("The total number of D28 aliquots is: ", nrow(scATACseq_placebo_metadata_sheet_df[scATACseq_placebo_metadata_sheet_df$time_point == "D28",])))
print(paste0("The total number of unique subjects is: ", length(unique(scATACseq_placebo_metadata_sheet_df$subject_id))))
print(paste0("The total number of male subjects is: ", sum(find_sex_count(scATACseq_placebo_metadata_sheet_df) == "M")))
print(paste0("The total number of female subjects is: ", sum(find_sex_count(scATACseq_placebo_metadata_sheet_df) == "F")))
# Note terrible sex ratio (4 male vs 9 female)
scATACseq_paired_timepoints_placebo_metadata_sheet_df <- scATACseq_placebo_metadata_sheet_df[scATACseq_placebo_metadata_sheet_df$subject_id %in% names(table(scATACseq_placebo_metadata_sheet_df$subject_id)[table(scATACseq_placebo_metadata_sheet_df$subject_id) == 2]),]
print(paste0("The total number of unique subjects with both D-1 and D28 timepoints is: ", length(unique(scATACseq_paired_timepoints_placebo_metadata_sheet_df$subject_id))))
print(paste0("The total number of male subjects with both D-1 and D28 timepoints is: ", sum(find_sex_count(scATACseq_paired_timepoints_placebo_metadata_sheet_df) == "M")))
print(paste0("The total number of female subjects with both D-1 and D28 timepoints is: ", sum(find_sex_count(scATACseq_paired_timepoints_placebo_metadata_sheet_df) == "F")))
# Note terrible sex ratio (3 male vs 6 female)
print("*********** PAIRED SCRNA-SEQ AND SCATAC-SEQ (PLACEBO) FOR MAGICAL ***********")
paired_placebo_metadata_sheet_df <- placebo_metadata_sheet_df[(placebo_metadata_sheet_df$scRNA_seq == TRUE & 
                                                                 placebo_metadata_sheet_df$scATAC_seq == TRUE),]
print(paste0("The total number of aliquots is: ", nrow(paired_placebo_metadata_sheet_df)))
print(paste0("The total number of D-1 aliquots is: ", nrow(paired_placebo_metadata_sheet_df[paired_placebo_metadata_sheet_df$time_point == "D-1",])))
print(paste0("The total number of D28 aliquots is: ", nrow(paired_placebo_metadata_sheet_df[paired_placebo_metadata_sheet_df$time_point == "D28",])))
print(paste0("The total number of unique subjects is: ", length(unique(paired_placebo_metadata_sheet_df$subject_id))))
print(paste0("The total number of male subjects is: ", sum(find_sex_count(paired_placebo_metadata_sheet_df) == "M")))
print(paste0("The total number of female subjects is: ", sum(find_sex_count(paired_placebo_metadata_sheet_df) == "F")))
# Note terrible sex ratio (3 male vs 9 female)
paired_paired_timepoints_placebo_metadata_sheet_df <- paired_placebo_metadata_sheet_df[paired_placebo_metadata_sheet_df$subject_id %in% names(table(paired_placebo_metadata_sheet_df$subject_id)[table(paired_placebo_metadata_sheet_df$subject_id) == 2]),]
print(paste0("The total number of unique subjects with both D-1 and D28 timepoints is: ", length(unique(paired_paired_timepoints_placebo_metadata_sheet_df$subject_id))))
print(paste0("The total number of male subjects with both D-1 and D28 timepoints is: ", sum(find_sex_count(paired_paired_timepoints_placebo_metadata_sheet_df) == "M")))
print(paste0("The total number of female subjects with both D-1 and D28 timepoints is: ", sum(find_sex_count(paired_paired_timepoints_placebo_metadata_sheet_df) == "F")))
# Note terrible sex ratio (2 male vs 4 female)
print("*********** PAIRED SCRNA-SEQ THAT PASSED QC AND SCATAC-SEQ (PLACEBO) FOR MAGICAL ***********")
passed_qc_paired_placebo_metadata_sheet_df <- paired_placebo_metadata_sheet_df[paired_placebo_metadata_sheet_df$passed_qc_scRNA_seq  == "Good" | paired_placebo_metadata_sheet_df$passed_qc_scRNA_seq  == "Pass",]
print(paste0("The total number of aliquots is: ", nrow(passed_qc_paired_placebo_metadata_sheet_df)))
print(paste0("The total number of D-1 aliquots is: ", nrow(passed_qc_paired_placebo_metadata_sheet_df[passed_qc_paired_placebo_metadata_sheet_df$time_point == "D-1",])))
print(paste0("The total number of D28 aliquots is: ", nrow(passed_qc_paired_placebo_metadata_sheet_df[passed_qc_paired_placebo_metadata_sheet_df$time_point == "D28",])))
print(paste0("The total number of unique subjects is: ", length(unique(passed_qc_paired_placebo_metadata_sheet_df$subject_id))))
print(paste0("The total number of male subjects is: ", sum(find_sex_count(passed_qc_paired_placebo_metadata_sheet_df) == "M")))
print(paste0("The total number of female subjects is: ", sum(find_sex_count(passed_qc_paired_placebo_metadata_sheet_df) == "F")))
# Note terrible sex ratio (2 male vs 5 female)
passed_qc_paired_paired_timepoints_placebo_metadata_sheet_df <- passed_qc_paired_placebo_metadata_sheet_df[passed_qc_paired_placebo_metadata_sheet_df$subject_id %in% names(table(passed_qc_paired_placebo_metadata_sheet_df$subject_id)[table(passed_qc_paired_placebo_metadata_sheet_df$subject_id) == 2]),]
print(paste0("The total number of unique subjects with both D-1 and D28 timepoints is: ", length(unique(passed_qc_paired_paired_timepoints_placebo_metadata_sheet_df$subject_id))))
print(paste0("The total number of male subjects with both D-1 and D28 timepoints is: ", sum(find_sex_count(passed_qc_paired_paired_timepoints_placebo_metadata_sheet_df) == "M")))
print(paste0("The total number of female subjects with both D-1 and D28 timepoints is: ", sum(find_sex_count(passed_qc_paired_paired_timepoints_placebo_metadata_sheet_df) == "F")))
# Note sex ratio (1 male vs 2 female)
print("*********** MULTIOME (PLACEBO) - NOTE THAT ALL SAMPLES PASSED QC ***********")
multiome_placebo_metadata_sheet_df <- placebo_metadata_sheet_df[placebo_metadata_sheet_df$multiome == TRUE,]
print(paste0("The total number of aliquots is: ", nrow(multiome_placebo_metadata_sheet_df)))
print(paste0("The total number of D-1 aliquots is: ", nrow(multiome_placebo_metadata_sheet_df[multiome_placebo_metadata_sheet_df$time_point == "D-1",])))
print(paste0("The total number of D28 aliquots is: ", nrow(multiome_placebo_metadata_sheet_df[multiome_placebo_metadata_sheet_df$time_point == "D28",])))
print(paste0("The total number of unique subjects is: ", length(unique(multiome_placebo_metadata_sheet_df$subject_id))))
print(paste0("The total number of male subjects is: ", sum(find_sex_count(multiome_placebo_metadata_sheet_df) == "M")))
print(paste0("The total number of female subjects is: ", sum(find_sex_count(multiome_placebo_metadata_sheet_df) == "F")))
# Note mediocre sex ratio (4 male vs 7 female)
multiome_paired_timepoints_placebo_metadata_sheet_df <- multiome_placebo_metadata_sheet_df[multiome_placebo_metadata_sheet_df$subject_id %in% names(table(multiome_placebo_metadata_sheet_df$subject_id)[table(multiome_placebo_metadata_sheet_df$subject_id) == 2]),]
print(paste0("The total number of unique subjects with both D-1 and D28 timepoints is: ", length(unique(multiome_paired_timepoints_placebo_metadata_sheet_df$subject_id))))
print(paste0("The total number of male subjects with both D-1 and D28 timepoints is: ", sum(find_sex_count(multiome_paired_timepoints_placebo_metadata_sheet_df) == "M")))
print(paste0("The total number of female subjects with both D-1 and D28 timepoints is: ", sum(find_sex_count(multiome_paired_timepoints_placebo_metadata_sheet_df) == "F")))
# Note mediocre sex ratio (2 male vs 6 female)
print("*********** BULK RNA-SEQ (PLACEBO) ***********")
bulkRNAseq_placebo_metadata_sheet_df <- placebo_metadata_sheet_df[placebo_metadata_sheet_df$bulkRNA_seq == TRUE,]
print(paste0("The total number of aliquots is: ", nrow(bulkRNAseq_placebo_metadata_sheet_df)))
print(paste0("The total number of D-1 (Period 1) aliquots is: ", nrow(bulkRNAseq_placebo_metadata_sheet_df[bulkRNAseq_placebo_metadata_sheet_df$time_point == "D1 predose" & bulkRNAseq_placebo_metadata_sheet_df$period == "1",])))
print(paste0("The total number of D2 (Period 1) aliquots is: ", nrow(bulkRNAseq_placebo_metadata_sheet_df[bulkRNAseq_placebo_metadata_sheet_df$time_point == "D2" & bulkRNAseq_placebo_metadata_sheet_df$period == "1",])))
print(paste0("The total number of D8 (Period 1) aliquots is: ", nrow(bulkRNAseq_placebo_metadata_sheet_df[bulkRNAseq_placebo_metadata_sheet_df$time_point == "D8" & bulkRNAseq_placebo_metadata_sheet_df$period == "1",])))
print(paste0("The total number of D28 (Period 1) aliquots is: ", nrow(bulkRNAseq_placebo_metadata_sheet_df[bulkRNAseq_placebo_metadata_sheet_df$time_point == "D28" & bulkRNAseq_placebo_metadata_sheet_df$period == "1",])))
print(paste0("The total number of D-2 (Period 2) aliquots is: ", nrow(bulkRNAseq_placebo_metadata_sheet_df[bulkRNAseq_placebo_metadata_sheet_df$time_point == "D-2" & bulkRNAseq_placebo_metadata_sheet_df$period == "2",])))
print(paste0("The total number of D-1 (Period 2) aliquots is: ", nrow(bulkRNAseq_placebo_metadata_sheet_df[bulkRNAseq_placebo_metadata_sheet_df$time_point == "D-1" & bulkRNAseq_placebo_metadata_sheet_df$period == "2",])))
print(paste0("The total number of D2 (Period 2) aliquots is: ", nrow(bulkRNAseq_placebo_metadata_sheet_df[bulkRNAseq_placebo_metadata_sheet_df$time_point == "D2" & bulkRNAseq_placebo_metadata_sheet_df$period == "2",])))
print(paste0("The total number of D5 (Period 2) aliquots is: ", nrow(bulkRNAseq_placebo_metadata_sheet_df[bulkRNAseq_placebo_metadata_sheet_df$time_point == "D5" & bulkRNAseq_placebo_metadata_sheet_df$period == "2",])))
print(paste0("The total number of D8 (Period 2) aliquots is: ", nrow(bulkRNAseq_placebo_metadata_sheet_df[bulkRNAseq_placebo_metadata_sheet_df$time_point == "D8" & bulkRNAseq_placebo_metadata_sheet_df$period == "2",])))
print(paste0("The total number of D28 (Period 2) aliquots is: ", nrow(bulkRNAseq_placebo_metadata_sheet_df[bulkRNAseq_placebo_metadata_sheet_df$time_point == "D28" & bulkRNAseq_placebo_metadata_sheet_df$period == "2",])))
print(paste0("The total number of unique subjects is: ", length(unique(bulkRNAseq_placebo_metadata_sheet_df$subject_id))))
print(paste0("The total number of male subjects is: ", sum(find_sex_count(bulkRNAseq_placebo_metadata_sheet_df) == "M")))
print(paste0("The total number of female subjects is: ", sum(find_sex_count(bulkRNAseq_placebo_metadata_sheet_df) == "F")))
# Good sex ratio (22 male vs 24 female)
print(paste0("The total number of age 1 subjects is: ", sum(find_age_count(bulkRNAseq_placebo_metadata_sheet_df) == 1)))
print(paste0("The total number of age 2 subjects is: ", sum(find_age_count(bulkRNAseq_placebo_metadata_sheet_df) == 2)))
print(paste0("The total number of age 3 subjects is: ", sum(find_age_count(bulkRNAseq_placebo_metadata_sheet_df) == 3)))
print(paste0("The total number of age 4 subjects is: ", sum(find_age_count(bulkRNAseq_placebo_metadata_sheet_df) == 4)))
bulkRNAseq_all_timepoints_placebo_metadata_sheet_df <- bulkRNAseq_placebo_metadata_sheet_df[bulkRNAseq_placebo_metadata_sheet_df$subject_id %in% names(table(bulkRNAseq_placebo_metadata_sheet_df$subject_id)[table(bulkRNAseq_placebo_metadata_sheet_df$subject_id) == 10]),]
print(paste0("The total number of unique subjects with both D-1 and D28 timepoints is: ", length(unique(bulkRNAseq_all_timepoints_placebo_metadata_sheet_df$subject_id))))
print(paste0("The total number of male subjects with both D-1 and D28 timepoints is: ", sum(find_sex_count(bulkRNAseq_all_timepoints_placebo_metadata_sheet_df) == "M")))
print(paste0("The total number of female subjects with both D-1 and D28 timepoints is: ", sum(find_sex_count(bulkRNAseq_all_timepoints_placebo_metadata_sheet_df) == "F")))
# Note bad sex ratio (7 male vs 16 female)
print("*********** BULK RNA-SEQ (VACCINATED) ***********")
bulkRNAseq_vaccinated_metadata_sheet_df <- vaccinated_metadata_sheet_df[vaccinated_metadata_sheet_df$bulkRNA_seq == TRUE,]
print(paste0("The total number of aliquots is: ", nrow(bulkRNAseq_vaccinated_metadata_sheet_df)))
print(paste0("The total number of D-1 (Period 1) aliquots is: ", nrow(bulkRNAseq_vaccinated_metadata_sheet_df[bulkRNAseq_vaccinated_metadata_sheet_df$time_point == "D1 predose" & bulkRNAseq_vaccinated_metadata_sheet_df$period == "1",])))
print(paste0("The total number of D2 (Period 1) aliquots is: ", nrow(bulkRNAseq_vaccinated_metadata_sheet_df[bulkRNAseq_vaccinated_metadata_sheet_df$time_point == "D2" & bulkRNAseq_vaccinated_metadata_sheet_df$period == "1",])))
print(paste0("The total number of D8 (Period 1) aliquots is: ", nrow(bulkRNAseq_vaccinated_metadata_sheet_df[bulkRNAseq_vaccinated_metadata_sheet_df$time_point == "D8" & bulkRNAseq_vaccinated_metadata_sheet_df$period == "1",])))
print(paste0("The total number of D28 (Period 1) aliquots is: ", nrow(bulkRNAseq_vaccinated_metadata_sheet_df[bulkRNAseq_vaccinated_metadata_sheet_df$time_point == "D28" & bulkRNAseq_vaccinated_metadata_sheet_df$period == "1",])))
print(paste0("The total number of D-2 (Period 2) aliquots is: ", nrow(bulkRNAseq_vaccinated_metadata_sheet_df[bulkRNAseq_vaccinated_metadata_sheet_df$time_point == "D-2" & bulkRNAseq_vaccinated_metadata_sheet_df$period == "2",])))
print(paste0("The total number of D-1 (Period 2) aliquots is: ", nrow(bulkRNAseq_vaccinated_metadata_sheet_df[bulkRNAseq_vaccinated_metadata_sheet_df$time_point == "D-1" & bulkRNAseq_vaccinated_metadata_sheet_df$period == "2",])))
print(paste0("The total number of D2 (Period 2) aliquots is: ", nrow(bulkRNAseq_vaccinated_metadata_sheet_df[bulkRNAseq_vaccinated_metadata_sheet_df$time_point == "D2" & bulkRNAseq_vaccinated_metadata_sheet_df$period == "2",])))
print(paste0("The total number of D5 (Period 2) aliquots is: ", nrow(bulkRNAseq_vaccinated_metadata_sheet_df[bulkRNAseq_vaccinated_metadata_sheet_df$time_point == "D5" & bulkRNAseq_vaccinated_metadata_sheet_df$period == "2",])))
print(paste0("The total number of D8 (Period 2) aliquots is: ", nrow(bulkRNAseq_vaccinated_metadata_sheet_df[bulkRNAseq_vaccinated_metadata_sheet_df$time_point == "D8" & bulkRNAseq_vaccinated_metadata_sheet_df$period == "2",])))
print(paste0("The total number of D28 (Period 2) aliquots is: ", nrow(bulkRNAseq_vaccinated_metadata_sheet_df[bulkRNAseq_vaccinated_metadata_sheet_df$time_point == "D28" & bulkRNAseq_vaccinated_metadata_sheet_df$period == "2",])))
print(paste0("The total number of unique subjects is: ", length(unique(bulkRNAseq_vaccinated_metadata_sheet_df$subject_id))))
print(paste0("The total number of male subjects is: ", sum(find_sex_count(bulkRNAseq_vaccinated_metadata_sheet_df) == "M")))
print(paste0("The total number of female subjects is: ", sum(find_sex_count(bulkRNAseq_vaccinated_metadata_sheet_df) == "F")))
print(paste0("The total number of white subjects is: ", sum(find_race_count(bulkRNAseq_vaccinated_metadata_sheet_df) == "WHITE")))
print(paste0("The total number of Asian subjects is: ", sum(find_race_count(bulkRNAseq_vaccinated_metadata_sheet_df) == "ASIAN")))
print(paste0("The total number of black subjects is: ", sum(find_race_count(bulkRNAseq_vaccinated_metadata_sheet_df) == "BLACK OR AFRICAN AMERICAN")))
# Good sex ratio (22 male vs 24 female)
print(paste0("The total number of age 1 subjects is: ", sum(find_age_count(bulkRNAseq_vaccinated_metadata_sheet_df) == 1)))
print(paste0("The total number of age 2 subjects is: ", sum(find_age_count(bulkRNAseq_vaccinated_metadata_sheet_df) == 2)))
print(paste0("The total number of age 3 subjects is: ", sum(find_age_count(bulkRNAseq_vaccinated_metadata_sheet_df) == 3)))
print(paste0("The total number of age 4 subjects is: ", sum(find_age_count(bulkRNAseq_vaccinated_metadata_sheet_df) == 4)))
# Overall conclusions:
# In general, really bad sex ratios throughout (way more females than males)
# But Vincy didn't correct for it in her analysis, so maybe it's OK? 
