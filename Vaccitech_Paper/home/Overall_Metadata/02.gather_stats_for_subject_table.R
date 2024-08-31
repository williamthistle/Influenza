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
# Subset to the 11 subjects that we actually analyze
mintchip_paired_samples <- c("17f82eccf700879f","ac2a58addd62cd54","551b6363face39ec","561e4c6a66f8b1e7",
                             "54079742b6a3cc8c","db7f26dc0f695274","1e5232f7d77ffbcf","83b86144423425e6",
                             "97053ff8c6b8623c","ce08d2c374a48d3c","00e4bf2b268c136e","18962e531bca3732",
                             "546421cbf94adc40","d51dc1a1431f87d6","18a86a5137ad6bf4","6892e26fc1a12f37",
                             "e0e22d32d1fbf6fa","f8e9a81fc752cc09","bbb598e6ac174c1d","cf68fb9f4b61f719",
                             "5a2aeba907c7baf4","e767550290c513c2")
mintchip_placebo_paired_samples_report <- mintchip_placebo_report[mintchip_placebo_report$aliquot_id %in% mintchip_paired_samples,]
summarize_stats(mintchip_placebo_paired_samples_report)

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
