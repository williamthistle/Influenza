my_datasets <- c()
for(dataset in bacteria_list) {
  my_datasets <- c(my_datasets, dataset$formattedName)
}

bacteria_metadata <- study_metadata[study_metadata$Study %in% my_datasets,]


final_sample_sizes <- c()
for(bacteria_dataset_id in bacteria_metadata$Study) {
  for(bacteria_dataset in bacteria_list) {
    if(bacteria_dataset$formattedName == bacteria_dataset_id) {
      final_sample_sizes <- c(final_sample_sizes, nrow(bacteria_dataset$pheno))
    }
  }
}
bacteria_metadata$final_sample_size <- final_sample_sizes

# FILTERING STRATEGY 1
flagged_bacteria_datasets <- c(5, 6, 11, 12, 14, 17, 20, 22, 27, 36, 37, 38, 39, 46, 47, 48, 52)
filtered_bacteria_metadata <- bacteria_metadata[-flagged_bacteria_datasets,]
filtered_bacteria_metadata <- filtered_bacteria_metadata[filtered_bacteria_metadata$Standardized.Exposure != "Unknown",]
filtered_bacteria_metadata <- filtered_bacteria_metadata[!grepl("Influenza virus", filtered_bacteria_metadata$Standardized.Exposure),]
filtered_bacteria_metadata <- filtered_bacteria_metadata[filtered_bacteria_metadata$final_sample_size >= 10,]
filtered_bacteria_list <- bacteria_list[my_datasets %in% filtered_bacteria_metadata$Study]
 
# Fix time series
filtered_bacteria_list[[27]]$pheno <- filtered_bacteria_list[[27]]$pheno[filtered_bacteria_list[[27]]$pheno$`Time.Point` == "HC" | filtered_bacteria_list[[27]]$pheno$`Time.Point` == "-T1",]
filtered_bacteria_list[[27]]$pheno <- filtered_bacteria_list[[27]]$pheno[filtered_bacteria_list[[27]]$pheno$`SubjectID` != "HC3" & filtered_bacteria_list[[27]]$pheno$`SubjectID` != "HC4"
                                                                         & filtered_bacteria_list[[27]]$pheno$`SubjectID` != "HC5" & filtered_bacteria_list[[27]]$pheno$`SubjectID` != "HC6",]
kept_accessions_27 <- rownames(filtered_bacteria_list[[27]]$pheno)
filtered_bacteria_list[[27]]$expr <- filtered_bacteria_list[[27]]$expr[,colnames(filtered_bacteria_list[[27]]$expr) %in% kept_accessions_27]
filtered_bacteria_list[[27]]$class <- filtered_bacteria_list[[27]]$class[names(filtered_bacteria_list[[27]]$class) %in% kept_accessions_27]

# Capture updated sample sizes
final_sample_sizes <- c()
for(bacteria_dataset_id in filtered_bacteria_metadata$Study) {
  for(bacteria_dataset in filtered_bacteria_list) {
    if(bacteria_dataset$formattedName == bacteria_dataset_id) {
      final_sample_sizes <- c(final_sample_sizes, nrow(bacteria_dataset$pheno))
    }
  }
}
filtered_bacteria_metadata$final_sample_size <- final_sample_sizes

total_sample_sizes <- -1
platform_check <- FALSE
pbmc_check <- FALSE
adult_children_ratio <- -1
while((total_sample_sizes < 1300 | total_sample_sizes > 1500) | platform_check == FALSE | pbmc_check == FALSE) {
  potential_discovery_metadata <- filtered_bacteria_metadata[sample(nrow(filtered_bacteria_metadata), 18), ]
  total_sample_sizes <- sum(potential_discovery_metadata$final_sample_size)
  if(sum(potential_discovery_metadata$Manufacturer == "affy", na.rm = TRUE) < 3) {
    platform_check <- TRUE 
  } else {
    platform_check <- FALSE
  }
  if(sum(potential_discovery_metadata$Tissue == "PBMCs", na.rm = TRUE) == 1) {
    pbmc_check <- TRUE 
  } else {
    pbmc_check <- FALSE
  }
}

print(sum(potential_discovery_metadata$final_sample_size))
potential_validation_metadata <- filtered_bacteria_metadata[!(filtered_bacteria_metadata$Study %in% potential_discovery_metadata$Study),]
print(sum(potential_validation_metadata$final_sample_size))
bacteria_discovery_metadata <- potential_discovery_metadata
bacteria_validation_metadata <- potential_validation_metadata
save(bacteria_discovery_metadata, file = "bacteria_discovery_metadata.rda")
save(bacteria_validation_metadata, file = "bacteria_validation_metadata.rda")


