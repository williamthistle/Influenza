base_dir <- "C:/Users/wat2/Documents/GitHub/Influenza/"
bulkRNA_data_list <- paste0(base_dir, "bulkRNA/bulkRNA_data_list.txt")
overall_metadata_file <- paste0(base_dir, "20220609_metadataECHO_Vaccitech_Coded.csv")
# Read in tables
bulkRNA_data <- read.table(bulkRNA_data_list)$V1
overall_metadata <- read.csv(overall_metadata_file)
bulkRNA_all_metadata_sheet_df <- data.frame(aliquot_id = character(), subject_id = character(), has_metadata = character(),
                                      treatment = character(), period = character(), time_point = character(), 
                                    sex = character())
for(bulkRNA_entry in bulkRNA_data) {
  # Add aliquot ID to current row
  current_row <- c()
  current_row <- append(current_row, bulkRNA_entry)
  # Add subject ID to current row
  current_sample = overall_metadata[overall_metadata$X_aliquot_id == bulkRNA_entry,]
  if(nrow(current_sample) > 0) {
    current_row <- append(current_row, current_sample$SUBJECT_ID)
    current_row <- append(current_row, TRUE)
    current_row <- append(current_row, current_sample$TREATMENT)
    current_row <- append(current_row, current_sample$Period)
    current_row <- append(current_row, current_sample$Time_Point)
    current_row <- append(current_row, current_sample$SEX)
  } else {
    current_row <- append(current_row, "N/A")
    current_row <- append(current_row, FALSE)
    current_row <- append(current_row, "N/A")
    current_row <- append(current_row, "N/A")
    current_row <- append(current_row, "N/A")
    current_row <- append(current_row, "N/A")
  }
  bulkRNA_all_metadata_sheet_df[nrow(bulkRNA_all_metadata_sheet_df) + 1,] = current_row
}

bulkRNA_all_metadata_sheet_df <- bulkRNA_all_metadata_sheet_df[order(bulkRNA_all_metadata_sheet_df$subject_id),]
write.table(bulkRNA_all_metadata_sheet_df, paste0(base_dir, "bulkRNA_all_metadata_sheet.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)

bulkRNA_placebo_metadata_sheet_df <- bulkRNA_all_metadata_sheet_df[bulkRNA_all_metadata_sheet_df$treatment == "PLACEBO",]
write.table(bulkRNA_placebo_metadata_sheet_df, paste0(base_dir, "bulkRNA_placebo_metadata_sheet.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)

sexes <- c()
visited_subjects <- c()
for (current_row in 1:nrow(bulkRNA_placebo_metadata_sheet_df)) {
  if(!bulkRNA_placebo_metadata_sheet_df[current_row,]$subject_id %in% visited_subjects ) {
    sexes <- append(sexes, bulkRNA_placebo_metadata_sheet_df[current_row,]$sex)
    visited_subjects <- append(visited_subjects, bulkRNA_placebo_metadata_sheet_df[current_row,]$subject_id)
  }
}

print(paste0("The total number of aliquots is: ", nrow(bulkRNA_all_metadata_sheet_df)))
print(paste0("The total number of placebo aliquots is: ", nrow(bulkRNA_placebo_metadata_sheet_df)))
print(paste0("The total number of D-1 aliquots is: ", nrow(bulkRNA_placebo_metadata_sheet_df[bulkRNA_placebo_metadata_sheet_df$time_point == "D-1",])))
print(paste0("The total number of D28 aliquots is: ", nrow(bulkRNA_placebo_metadata_sheet_df[bulkRNA_placebo_metadata_sheet_df$time_point == "D28",])))
print(paste0("The total number of unique placebo subjects is: ", length(unique(bulkRNA_placebo_metadata_sheet_df$subject_id))))
print(paste0("The total number of male subjects is: ", sum(sexes == "M")))
print(paste0("The total number of female subjects is: ", sum(sexes == "F")))
print(paste0("The total number of subjects with all timepoints is: ", length(unique(bulkRNA_placebo_metadata_sheet_df$subject_id)) - sum(table(bulkRNA_placebo_metadata_sheet_df$subject_id) < 10)))




