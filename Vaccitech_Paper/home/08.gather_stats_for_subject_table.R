# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

summarize_stats <- function(current_table) {
  current_table_subject_subset <- current_table[!duplicated(current_table$subject_id), ]
  female_current_table_subject_subset <- current_table_subject_subset[current_table_subject_subset$sex == "F",]
  male_current_table_subject_subset <- current_table_subject_subset[current_table_subject_subset$sex == "M",]
  
  ## Sex
  print("SEX DISTRIBUTION")
  print(table(current_table_subject_subset$sex))
  
  ## Age
  print("AGE DISTRIBUTION")
  print("ALL")
  bulk_all_mean_age <- mean(current_table_subject_subset$raw_age)
  print(bulk_all_mean_age)
  bulk_all_sd_age <- sd(current_table_subject_subset$raw_age)
  print(bulk_all_sd_age)
  
  print("FEMALE")
  bulk_female_mean_age <- mean(female_current_table_subject_subset$raw_age)
  print(bulk_female_mean_age)
  bulk_female_sd_age <- sd(female_current_table_subject_subset$raw_age)
  print(bulk_female_sd_age)
  
  print("MALE")
  bulk_male_mean_age <- mean(male_current_table_subject_subset$raw_age)
  print(bulk_male_mean_age)
  bulk_male_sd_age <- sd(male_current_table_subject_subset$raw_age)
  print(bulk_male_sd_age)
  
  ##  Race
  print("RACE DISTRIBUTION")
  print(table(current_table_subject_subset$race))
}

placebo_metadata_for_table <- sample_metadata[sample_metadata$treatment == "PLACEBO",]

placebo_metadata_for_table$time_point <- paste0(placebo_metadata_for_table$period, "_", placebo_metadata_for_table$time_point)
placebo_metadata_for_table$time_point[placebo_metadata_for_table$time_point == '1_D1 predose'] <- '1_D_minus_1'
placebo_metadata_for_table$time_point[placebo_metadata_for_table$time_point == '2_D-2'] <- '2_D_minus_2'
placebo_metadata_for_table$time_point[placebo_metadata_for_table$time_point == '2_D-1'] <- '2_D_minus_1'

### Bulk RNA-Seq
bulk_placebo_metadata_for_table <- placebo_metadata_for_table[placebo_metadata_for_table$bulkRNA_seq == TRUE,]
summarize_stats(bulk_placebo_metadata_for_table)

### scRNA-seq / snATAC-seq
scRNA_aliquots <- flu_tokens[flu_tokens$token == "single_cell_paired_sample",]$samples
scRNA_aliquots <- unlist(strsplit(scRNA_aliquots, ","))
scRNA_placebo_metadata_for_table <- placebo_metadata_for_table[placebo_metadata_for_table$aliquot_id %in% scRNA_snATAC_aliquots,]
summarize_stats(scRNA_placebo_metadata_for_table)

### sn Multiome
multiome_aliquots <- flu_tokens[flu_tokens$token == "all_multiome_paired_minus_0_sample",]$samples
multiome_aliquots <- unlist(strsplit(multiome_aliquots, ","))
multiome_placebo_metadata_for_table <- placebo_metadata_for_table[placebo_metadata_for_table$aliquot_id %in% multiome_aliquots,]
summarize_stats(multiome_placebo_metadata_for_table)

# snME
snME_placebo_metadata_for_table <- placebo_metadata_for_table[placebo_metadata_for_table$snME == TRUE,]
summarize_stats(snME_placebo_metadata_for_table)

# Bulk methylation
bulk_methyl_placebo_metadata_for_table <- placebo_metadata_for_table[placebo_metadata_for_table$bulk_methylation == TRUE,]
summarize_stats(bulk_methyl_placebo_metadata_for_table)

# miRNA-seq is same as snME

# Mint-ChIP
mintchip_metadata_for_table <- mintchip_metadata[mintchip_metadata$Subject  %in% names(table(mintchip_metadata$Subject)[table(mintchip_metadata$Subject) == 2]),]
mintchip_placebo_metadata_for_table <- placebo_metadata_for_table[placebo_metadata_for_table$aliquot_id %in% mintchip_metadata_for_table$Aliquot,]
summarize_stats(mintchip_placebo_metadata_for_table)

