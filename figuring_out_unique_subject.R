single_cell_stuff <- all_metadata[all_metadata$scRNA_seq == TRUE | all_metadata$scATAC_seq == TRUE | all_metadata$multiome == TRUE,]
single_cell_stuff <- single_cell_stuff[single_cell_stuff$treatment == "PLACEBO",]
single_cell_subjects <- unique(single_cell_stuff$subject_id)
for(subject in single_cell_subjects) {
  current_metadata <- all_metadata[all_metadata$subject_id == subject,]
  print(nrow(current_metadata))
  if(nrow(current_metadata) == 5) {
    print(current_metadata)
  }
}