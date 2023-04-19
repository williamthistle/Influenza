my_datasets <- c()
for(dataset in noninfectious_list) {
  my_datasets <- c(my_datasets, dataset$formattedName)
}

noninfectious_metadata <- lapply(noninfectious_list, function(x){
  data.frame('Study' = x$formattedName, 'Condition' = unique(x$pheno$Standardized.Exposure))
}) %>%
  bind_rows()

platform <- c()
sample_size <- c()
for(current_row_index in 1:nrow(noninfectious_metadata)) {
  current_row <- noninfectious_metadata[current_row_index,]
  current_study <- current_row$Study
  current_platform <- strsplit(current_study, "_")[[1]][2]
  print(current_platform)
  if(current_platform == "GPL10558") {
    platform <- c(platform, "ilmn")
  } else if(current_platform == "GPL7020") {
    platform <- c(platform, "affy")
  } else if(current_platform == "GPL15988") {
    platform <- c(platform, "ilmn")
  } else if(current_platform == "GPL570") {
    platform <- c(platform, "affy")
  } else if(current_platform == "GPL6947") {
    platform <- c(platform, "ilmn")
  } else if(current_platform == "GPL11532") {
    platform <- c(platform, "affy")
  } else if(current_platform == "GPL13667") {
    platform <- c(platform, "affy")
  } else if(current_platform == "GPL6244") {
    platform <- c(platform, "affy")
  } else if(current_platform == "GPL6480") {
    platform <- c(platform, "agilent")
  } else if(current_platform == "GPL4133") {
    platform <- c(platform, "agilent")
  } else if(current_platform == "GPL10775") {
    platform <- c(platform, "Other")
  } else if(current_platform == "GPL6947") {
    platform <- c(platform, "ilmn")
  } else if(current_platform == "GPL13607") {
    platform <- c(platform, "agilent")
  } else if(current_platform == "GPL10904") {
    platform <- c(platform, "ilmn")
  } else if(current_platform == "GPL10558") {
    platform <- c(platform, "ilmn")   
  } else if(current_platform == "GPL13158") {
    platform <- c(platform, "affy")   
  } else if(current_platform == "GPL80") {
    platform <- c(platform, "affy")   
  } else if(current_platform == "GPL20844") {
    platform <- c(platform, "agilent")   
  } else {
    print(current_row_index)
  }
  sample_size <- c(sample_size, nrow(noninfectious_list[[current_row_index]]$pheno))
}
noninfectious_metadata$platform <- platform
noninfectious_metadata$sample_size <- sample_size

total_sample_sizes <- -1
while(total_sample_sizes < 1800 | total_sample_sizes > 2000) {
  potential_discovery_metadata <- noninfectious_metadata[sample(nrow(noninfectious_metadata), 23), ]
  total_sample_sizes <- sum(potential_discovery_metadata$sample_size)
}

print(sum(potential_discovery_metadata$sample_size))
potential_validation_metadata <- noninfectious_metadata[!(noninfectious_metadata$Study %in% potential_discovery_metadata$Study),]
print(sum(potential_validation_metadata$sample_size))
noninfectious_discovery_metadata <- potential_discovery_metadata
noninfectious_validation_metadata <- potential_validation_metadata
noninfectious_discovery_metadata
noninfectious_validation_metadata

save(noninfectious_discovery_metadata, file = "noninfectious_discovery_metadata.rda")
save(noninfectious_validation_metadata, file = "noninfectious_validation_metadata.rda")
