# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# NOTE: We currently remove 0 qPCRAUC samples and any samples that don't have viral loads from consideration
# If we're uploading these samples, we may want to include them in these tables
summarize_stats <- function(current_table) {
  ### Subject based stats
  current_table_subject_subset <- current_table[!duplicated(current_table$subject_id), ]
  female_current_table_subject_subset <- current_table_subject_subset[current_table_subject_subset$sex == "F",]
  male_current_table_subject_subset <- current_table_subject_subset[current_table_subject_subset$sex == "M",]
  
  ## Total
  print("TOTAL SUBJECTS")
  print(nrow(current_table_subject_subset))
  
  ## Sex
  print("SEX DISTRIBUTION")
  print(table(current_table_subject_subset$sex))
  
  ## Age
  print("AGE DISTRIBUTION")
  print("ALL")
  bulk_all_mean_age <- mean(current_table_subject_subset$raw_age)
  bulk_all_sd_age <- sd(current_table_subject_subset$raw_age)
  print(paste0(round(bulk_all_mean_age, 2), " \u00B1 ", round(bulk_all_sd_age, 2)))
  
  print("FEMALE")
  bulk_female_mean_age <- mean(female_current_table_subject_subset$raw_age)
  bulk_female_sd_age <- sd(female_current_table_subject_subset$raw_age)
  print(paste0(round(bulk_female_mean_age, 2), " \u00B1 ", round(bulk_female_sd_age, 2)))
  
  print("MALE")
  bulk_male_mean_age <- mean(male_current_table_subject_subset$raw_age)
  bulk_male_sd_age <- sd(male_current_table_subject_subset$raw_age)
  print(paste0(round(bulk_male_mean_age, 2), " \u00B1 ", round(bulk_male_sd_age, 2)))
  
  ## Race
  print("RACE DISTRIBUTION")
  print(table(current_table_subject_subset$race))
  
  ## Viral Load
  print("VIRAL LOAD DISTRIBUTION")
  print(table(current_table_subject_subset$viral_load_category))
  
  ### Sample based stats
  
  ## Total
  print("TOTAL SAMPLES")
  print(current_table)
  print(nrow(current_table))
  
  ## Sex
  print("SEX DISTRIBUTION")
  print(table(current_table$sex))
  
  
  ## Race
  print("RACE DISTRIBUTION")
  print(table(current_table$race))
  
  ## Viral Load
  print("VIRAL LOAD DISTRIBUTION")
  print(table(current_table$viral_load_category))
}

### PLACEBO REPORT
placebo_report <- all_metadata[all_metadata$treatment == "PLACEBO",]
# Remove any samples that don't have viral load (or 0 qPCRAUC)
placebo_report <- placebo_report[placebo_report$viral_load > 0,]

# Bulk RNA-Seq
bulk_placebo_report <- placebo_report[placebo_report$bulkRNA_seq == TRUE,]
summarize_stats(bulk_placebo_report)

# scRNA-seq
scRNA_placebo_report <- placebo_report[placebo_report$scRNA_seq == TRUE,]
summarize_stats(scRNA_placebo_report)

# scATAC-seq
scATAC_placebo_report <- placebo_report[placebo_report$scATAC_seq == TRUE,]
summarize_stats(scATAC_placebo_report)

# snMultiome
snMultiome_placebo_report <- placebo_report[placebo_report$multiome == TRUE,]
summarize_stats(snMultiome_placebo_report)

# snMethylation
snMethylation_placebo_report <- placebo_report[placebo_report$snME == TRUE,]
summarize_stats(snMethylation_placebo_report)

# Methylation
methylation_placebo_report <- placebo_report[placebo_report$bulk_methylation == TRUE,]
summarize_stats(methylation_placebo_report)

# miRNA
miRNA_placebo_report <- placebo_report[placebo_report$miRNA == TRUE,]
summarize_stats(miRNA_placebo_report)

# MintChIP
mintchip_placebo_report <- placebo_report[placebo_report$mintchip == TRUE,]
summarize_stats(mintchip_placebo_report)

### VACCINATED REPORT
vaccinated_report <- all_metadata[all_metadata$treatment == "MVA-NP+M1",]
# Remove any samples that don't have viral load info (marked as -1 viral load) or have no viral load (0 qPCRAUC)
vaccinated_report <- vaccinated_report[vaccinated_report$viral_load > 0,]

# Bulk RNA-Seq
bulk_vaccinated_report <- vaccinated_report[vaccinated_report$bulkRNA_seq == TRUE,]
summarize_stats(bulk_vaccinated_report)

# scRNA-seq
scRNA_vaccinated_report <- vaccinated_report[vaccinated_report$scRNA_seq == TRUE,]
summarize_stats(scRNA_vaccinated_report)

# scATAC-seq
scATAC_vaccinated_report <- vaccinated_report[vaccinated_report$scATAC_seq == TRUE,]
summarize_stats(scATAC_vaccinated_report)

# snMultiome
snMultiome_vaccinated_report <- vaccinated_report[vaccinated_report$multiome == TRUE,]
summarize_stats(snMultiome_vaccinated_report)

# Methylation
methylation_vaccinated_report <- vaccinated_report[vaccinated_report$bulk_methylation == TRUE,]
summarize_stats(methylation_vaccinated_report)
