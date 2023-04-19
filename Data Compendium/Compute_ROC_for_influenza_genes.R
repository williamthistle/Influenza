library(MetaIntegrator)
library(dplyr)
library(openxlsx)
library(ggplot2)

# Set working directory to wherever this script is
# (because Compendium_Functions.R will be in the same place)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("Compendium_Functions.R")
# Now set working directory to where compendium data files are
setwd("C:/Users/willi/Desktop/Influenza Work (April 2023)/Data Compendium")

# Grab study metadata (will be used to confirm that discovery / validation splits are balanced)
study_metadata <- read.xlsx('metadata_tables/darpa_compendium_metadata_v1.xlsx', sheet = 1)
# Grab Daniel's used influenza datasets - this will provide the IDs for the influenza and non-influenza virus datasets
# (We still need to grab IDs for bacteria and non-infectious)
daniel_paper_datasets <- read.xlsx('metadata_tables/Daniel_Paper_Datasets.xlsx', sheet = 1)

# All infectious datasets are found in data_list
data_list <- readRDS('Processed, Manuscript Files/data_list_v1.RDS')

# FIXING CERTAIN DATASETS
current_dataset <- data_list[[127]]
current_dataset$pheno[grepl("acute", current_dataset$pheno$title),]$Study <- "GSE97741_GPL10558"
current_dataset$pheno[grepl("acute", current_dataset$pheno$title),]$Class <- "Virus"
current_dataset$pheno[grepl("acute", current_dataset$pheno$title),]$Pathogen <- "Respiratory syncytial virus;Rhinovirus"
current_dataset$pheno[grepl("discharge", current_dataset$pheno$title),]$Study <- "GSE97741_GPL10558"
current_dataset$pheno[grepl("discharge", current_dataset$pheno$title),]$Class <- "Convalescent"
current_dataset$pheno[grepl("discharge", current_dataset$pheno$title),]$Pathogen <- "Healthy"
data_list[[127]] <- current_dataset

# We have to load noninfectious datasets separately
noninfectious_list <- readRDS('Processed, Manuscript Files/noninfectious_data_list.RDS')

# Single-cell RNA-seq MAGICAL genes (paired)
MAGICAL_single_cell_genes <- read.table("../Single Cell RNA-Seq/MAGICAL_flu_genes.csv")$V1

# True multiome RNA-seq MAGICAL genes (paired)
MAGICAL_multiome_14_genes <- read.table("../True Multiome/MAGICAL_flu_genes_14.csv")$V1

# True multiome RNA-seq MAGICAL genes (all)
MAGICAL_multiome_19_genes <- read.table("../True Multiome/MAGICAL_flu_genes_19.csv")$V1

# Find all influenza datasets
flu_samples <- sapply(data_list, findPathogen, pathogen = 'Influenza virus')
flu_list <- data_list[flu_samples]
flu_list <- lapply(flu_list, formatPathogenContrast, pathogen = 'Influenza virus')
filtered_samples <- sapply(flu_list, filterSampleSizes, N = 4)
flu_list <- flu_list[filtered_samples]

# Find all non-influenza virus datasets
non_flu_virus_list <- lapply(data_list, removePathogen, pathogen = 'Influenza virus', target_class = "Virus")
filtered_samples <- sapply(non_flu_virus_list, filterSampleSizes, N = 4)
non_flu_virus_list <- non_flu_virus_list[filtered_samples]

# Find all bacteria datasets
bacteria_list <- lapply(data_list, removePathogen, pathogen = 'Influenza virus', target_class = "Bacteria")
filtered_samples <- sapply(bacteria_list, filterSampleSizes, N = 4)
bacteria_list <- bacteria_list[filtered_samples]

# Find all non-infectious datasets
noninfectious_list <- lapply(noninfectious_list, formatContrast, target_class = "Other.NonInfectious")
filtered_samples <- sapply(noninfectious_list, filterSampleSizes, N = 4)
noninfectious_list <- noninfectious_list[filtered_samples]

# Find influenza discovery datasets
flu_discovery_list <- list()
flu_discovery_dataset_names <- c()
current_index <- 1
daniel_paper_discovery_datasets <- daniel_paper_datasets[daniel_paper_datasets$Usage == "Discovery",]
for(current_flu_dataset in flu_list) {
  current_accession <- current_flu_dataset$formattedName
  for(daniel_accession_index in 1:nrow(daniel_paper_discovery_datasets)) {
    if(grepl(daniel_paper_discovery_datasets[daniel_accession_index,]$Accession, current_accession)) {
      flu_discovery_dataset_names <- c(flu_discovery_dataset_names, current_accession)
      flu_discovery_list[[current_index]] <- current_flu_dataset
      current_index <- current_index + 1
    }
  }
}
flu_discovery_metadata <- study_metadata[study_metadata$Study %in% flu_discovery_dataset_names,]

# Find influenza validation datasets
flu_validation_list <- list()
flu_validation_dataset_names <- c()
current_index <- 1
daniel_paper_validation_datasets <- daniel_paper_datasets[daniel_paper_datasets$Usage == "Validation (robustness)",]
for(current_flu_dataset in flu_list) {
  current_accession <- current_flu_dataset$formattedName
  for(daniel_accession_index in 1:nrow(daniel_paper_validation_datasets)) {
    if(grepl(daniel_paper_validation_datasets[daniel_accession_index,]$Accession, current_accession)) {
      flu_validation_dataset_names <- c(flu_validation_dataset_names, current_accession)
      flu_validation_list[[current_index]] <- current_flu_dataset
      current_index <- current_index + 1
    }
  }
}
flu_validation_metadata <- study_metadata[study_metadata$Study %in% flu_validation_dataset_names,]

# FIXING TIME SERIES FOR INFLUENZA DISCOVERY
# Need to filter time series data so it only contains most acute time point and maybe one healthy time point (convalescent / healthy)
# Dataset 5
kept_indices_5 <- c(1, 5, 10, 11, 12, 13, 16, 19, 22, 26, 27, 28, 29, 30, 31, 32, 33, 35, 36, 37, 38, 39, 40)
kept_accessions_5 <- flu_discovery_list[[5]]$pheno$geo_accession[kept_indices_5]
flu_discovery_list[[5]]$pheno <- flu_discovery_list[[5]]$pheno[flu_discovery_list[[5]]$pheno$geo_accession %in% kept_accessions_5,]
flu_discovery_list[[5]]$expr <- flu_discovery_list[[5]]$expr[,colnames(flu_discovery_list[[5]]$expr) %in% kept_accessions_5]
flu_discovery_list[[5]]$class <- flu_discovery_list[[5]]$class[names(flu_discovery_list[[5]]$class) %in% kept_accessions_5]
# Dataset 6
kept_indices_6 <- c(1, 5, 10, 15, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37)
kept_accessions_6 <- flu_discovery_list[[6]]$pheno$geo_accession[kept_indices_6]
flu_discovery_list[[6]]$pheno <- flu_discovery_list[[6]]$pheno[flu_discovery_list[[6]]$pheno$geo_accession %in% kept_accessions_6,]
kept_accessions_6_alternate <- rownames(flu_discovery_list[[6]]$pheno)
flu_discovery_list[[6]]$expr <- flu_discovery_list[[6]]$expr[,colnames(flu_discovery_list[[6]]$expr) %in% kept_accessions_6_alternate]
flu_discovery_list[[6]]$class <- flu_discovery_list[[6]]$class[names(flu_discovery_list[[6]]$class) %in% kept_accessions_6_alternate]
# Dataset 7
flu_discovery_list[[7]]$pheno <- flu_discovery_list[[7]]$pheno[flu_discovery_list[[7]]$pheno$`Time.Point` != "late period",]
kept_accessions_7 <- rownames(flu_discovery_list[[7]]$pheno)
flu_discovery_list[[7]]$expr <- flu_discovery_list[[7]]$expr[,colnames(flu_discovery_list[[7]]$expr) %in% kept_accessions_7]
flu_discovery_list[[7]]$class <- flu_discovery_list[[7]]$class[names(flu_discovery_list[[7]]$class) %in% kept_accessions_7]
# Dataset 8
flu_discovery_list[[8]]$pheno <- flu_discovery_list[[8]]$pheno[flu_discovery_list[[8]]$pheno$`Time.Point` == "Less than 72 hour after symptoms",]
flu_discovery_list[[8]]$pheno <- flu_discovery_list[[8]]$pheno[flu_discovery_list[[8]]$pheno$Pathogen == "Influenza virus" | flu_discovery_list[[8]]$pheno$Pathogen == "Healthy",]
kept_accessions_8 <- rownames(flu_discovery_list[[8]]$pheno)
flu_discovery_list[[8]]$expr <- flu_discovery_list[[8]]$expr[,colnames(flu_discovery_list[[8]]$expr) %in% kept_accessions_8]
flu_discovery_list[[8]]$class <- flu_discovery_list[[8]]$class[names(flu_discovery_list[[8]]$class) %in% kept_accessions_8]
# Dataset 9
flu_discovery_list[[9]]$pheno <- flu_discovery_list[[9]]$pheno[flu_discovery_list[[9]]$pheno$`Time.Point` == "Baseline" | flu_discovery_list[[9]]$pheno$`Time.Point` == "Day0",]
kept_accessions_9 <- rownames(flu_discovery_list[[9]]$pheno)
flu_discovery_list[[9]]$expr <- flu_discovery_list[[9]]$expr[,colnames(flu_discovery_list[[9]]$expr) %in% kept_accessions_9]
flu_discovery_list[[9]]$class <- flu_discovery_list[[9]]$class[names(flu_discovery_list[[9]]$class) %in% kept_accessions_9]
# Dataset 10
flu_discovery_list[[10]]$pheno <- flu_discovery_list[[10]]$pheno[flu_discovery_list[[10]]$pheno$`Time.Point` == "Pre-Exposure" | flu_discovery_list[[10]]$pheno$`Time.Point` == "48 hours",]
kept_accessions_10 <- rownames(flu_discovery_list[[10]]$pheno)
flu_discovery_list[[10]]$expr <- flu_discovery_list[[10]]$expr[,colnames(flu_discovery_list[[10]]$expr) %in% kept_accessions_10]
flu_discovery_list[[10]]$class <- flu_discovery_list[[10]]$class[names(flu_discovery_list[[10]]$class) %in% kept_accessions_10]

flu_discovery_final_sample_sizes <- c()
for(flu_dataset_id in flu_discovery_metadata$Study) {
  for(flu_dataset in flu_discovery_list) {
    if(flu_dataset$formattedName == flu_dataset_id) {
      flu_discovery_final_sample_sizes <- c(flu_discovery_final_sample_sizes, nrow(flu_dataset$pheno))
    }
  }
}
flu_discovery_metadata$final_sample_size <- flu_discovery_final_sample_sizes

# FIXING TIME SERIES FOR INFLUENZA VALIDATION
# Need to filter time series data so it only contains most acute time point and maybe one healthy time point (convalescent / healthy)
# Dataset 3
flu_validation_list[[3]]$pheno <- flu_validation_list[[3]]$pheno[flu_validation_list[[3]]$pheno$`Time.Point` != "T2" & flu_validation_list[[3]]$pheno$`Time.Point` != "T3",]
kept_accessions_3 <- rownames(flu_validation_list[[3]]$pheno)
flu_validation_list[[3]]$expr <- flu_validation_list[[3]]$expr[,colnames(flu_validation_list[[3]]$expr) %in% kept_accessions_3]
flu_validation_list[[3]]$class <- flu_validation_list[[3]]$class[names(flu_validation_list[[3]]$class) %in% kept_accessions_3]
# Dataset 4
flu_validation_list[[4]]$pheno <- flu_validation_list[[4]]$pheno[flu_validation_list[[4]]$pheno$`Time.Point` == 1,]
kept_accessions_4 <- rownames(flu_validation_list[[4]]$pheno)
flu_validation_list[[4]]$expr <- flu_validation_list[[4]]$expr[,colnames(flu_validation_list[[4]]$expr) %in% kept_accessions_4]
flu_validation_list[[4]]$class <- flu_validation_list[[4]]$class[names(flu_validation_list[[4]]$class) %in% kept_accessions_4]

flu_validation_final_sample_sizes <- c()
for(flu_dataset_id in flu_validation_metadata$Study) {
  for(flu_dataset in flu_validation_list) {
    if(flu_dataset$formattedName == flu_dataset_id) {
      flu_validation_final_sample_sizes <- c(flu_validation_final_sample_sizes, nrow(flu_dataset$pheno))
    }
  }
}
flu_validation_metadata$final_sample_size <- flu_validation_final_sample_sizes

# Find non-influenza virus Daniel datasets
non_flu_virus_Daniel_list <- list()
non_flu_virus_Daniel_dataset_names <- c()
current_index <- 1
daniel_paper_cross_reactivity_datasets <- daniel_paper_datasets[daniel_paper_datasets$Usage == "Validation (cross-reactivity)",]
for(current_non_flu_virus_dataset in non_flu_virus_list) {
  current_accession <- current_non_flu_virus_dataset$formattedName
  for(daniel_accession_index in 1:nrow(daniel_paper_cross_reactivity_datasets)) {
    if(grepl(daniel_paper_cross_reactivity_datasets[daniel_accession_index,]$Accession, current_accession) & !(current_accession %in% non_flu_virus_Daniel_dataset_names)) {
      non_flu_virus_Daniel_dataset_names <- c(non_flu_virus_Daniel_dataset_names, current_accession)
      non_flu_virus_Daniel_list[[current_index]] <- current_non_flu_virus_dataset
      current_index <- current_index + 1
    }
  }
}
non_flu_virus_Daniel_metadata <- study_metadata[study_metadata$Study %in% non_flu_virus_Daniel_dataset_names,]

# # FIXING TIME SERIES FOR NON INFLUENZA VIRUS
# Need to filter time series data so it only contains most acute time point and maybe one healthy time point (convalescent / healthy)
# Dataset 10 (does it matter that this is also in flu_discovery? Maybe make sure it's in discovery set for non-influenza virus?)
non_flu_virus_Daniel_list[[10]]$pheno <- non_flu_virus_Daniel_list[[10]]$pheno[non_flu_virus_Daniel_list[[10]]$pheno$`Time.Point` == "Less than 72 hour after symptoms",]
kept_accessions_10 <- rownames(non_flu_virus_Daniel_list[[10]]$pheno)
non_flu_virus_Daniel_list[[10]]$expr <- non_flu_virus_Daniel_list[[10]]$expr[,colnames(non_flu_virus_Daniel_list[[10]]$expr) %in% kept_accessions_10]
non_flu_virus_Daniel_list[[10]]$class <- non_flu_virus_Daniel_list[[10]]$class[names(non_flu_virus_Daniel_list[[10]]$class) %in% kept_accessions_10]
# Dataset 11 (does it matter that this is also in flu_discovery? Maybe make sure it's in discovery set for non-influenza virus?)
non_flu_virus_Daniel_list[[11]]$pheno <- non_flu_virus_Daniel_list[[11]]$pheno[non_flu_virus_Daniel_list[[11]]$pheno$`Time.Point` == "Baseline" | non_flu_virus_Daniel_list[[11]]$pheno$`Time.Point` == "Day0",]
non_flu_virus_Daniel_list[[11]]$pheno <- non_flu_virus_Daniel_list[[11]]$pheno[!(grepl("Influenza virus", non_flu_virus_Daniel_list[[11]]$pheno$`Pathogen`)),]
kept_accessions_11 <- rownames(non_flu_virus_Daniel_list[[11]]$pheno)
non_flu_virus_Daniel_list[[11]]$expr <- non_flu_virus_Daniel_list[[11]]$expr[,colnames(non_flu_virus_Daniel_list[[11]]$expr) %in% kept_accessions_11]
non_flu_virus_Daniel_list[[11]]$class <- non_flu_virus_Daniel_list[[11]]$class[names(non_flu_virus_Daniel_list[[11]]$class) %in% kept_accessions_11]

# Add sample sizes and figure out which datasets overlap with influenza discovery / validation
non_flu_virus_Daniel_final_sample_sizes <- c()
non_flu_virus_Daniel_final_overlapping_with_flu <- c()
for(dataset_id in non_flu_virus_Daniel_metadata$Study) {
  for(non_flu_virus_Daniel_dataset in non_flu_virus_Daniel_list) {
    if(non_flu_virus_Daniel_dataset$formattedName == dataset_id) {
      non_flu_virus_Daniel_final_sample_sizes <- c(non_flu_virus_Daniel_final_sample_sizes, nrow(non_flu_virus_Daniel_dataset$pheno))
    }
  } 
  if(dataset_id %in% flu_discovery_metadata$Study) {
    non_flu_virus_Daniel_final_overlapping_with_flu <- c(non_flu_virus_Daniel_final_overlapping_with_flu, TRUE)
  } else {
    non_flu_virus_Daniel_final_overlapping_with_flu <- c(non_flu_virus_Daniel_final_overlapping_with_flu, FALSE)
  }
}
non_flu_virus_Daniel_metadata$final_sample_size <- non_flu_virus_Daniel_final_sample_sizes
non_flu_virus_Daniel_metadata$overlap_with_flu_discovery <- non_flu_virus_Daniel_final_overlapping_with_flu

# create discovery and validation datasets for non-influenza virus
non_flu_discovery_indices <- c(3, 4, 5, 6, 10, 11, 12)
non_flu_validation_indices <- c(1, 2, 7, 8, 9)
non_flu_virus_discovery_list <- non_flu_virus_Daniel_list[non_flu_discovery_indices]
non_flu_virus_validation_list <- non_flu_virus_Daniel_list[non_flu_validation_indices]

# BACTERIA
load("bacteria_discovery_metadata.rda")
load("bacteria_validation_metadata.rda")

bacteria_discovery_list <- bacteria_list[sort(match(bacteria_discovery_metadata$Study, bacteria_metadata$Study))]
bacteria_validation_list <- bacteria_list[sort(match(bacteria_validation_metadata$Study, bacteria_metadata$Study))]

# NON INFECTIOUS
load("noninfectious_discovery_metadata.rda")
load("noninfectious_validation_metadata.rda")

noninfectious_discovery_list <- noninfectious_list[sort(match(noninfectious_discovery_metadata$Study, noninfectious_metadata$Study))]
noninfectious_validation_list <- noninfectious_list[sort(match(noninfectious_validation_metadata$Study, noninfectious_metadata$Study))]

# Calculate flu (discovery) AUCs for individual genes
sc_discovery_flu_aucs <- test_individual_genes_on_datasets(MAGICAL_single_cell_genes, flu_discovery_list, "Single_Cell_Paired", "flu_discovery")
multiome_discovery_paired_flu_aucs <- test_individual_genes_on_datasets(MAGICAL_multiome_14_genes, flu_discovery_list, "Multiome_Paired", "flu_discovery")
multiome_discovery_all_flu_aucs <- test_individual_genes_on_datasets(MAGICAL_multiome_19_genes, flu_discovery_list, "Multiome_All", "flu_discovery")

# Calculate non-influenza virus (discovery) AUCs for individual genes
sc_discovery_non_flu_virus_aucs <- test_individual_genes_on_datasets(MAGICAL_single_cell_genes, non_flu_virus_discovery_list, "Single_Cell_Paired", "non_flu_virus_discovery")
multiome_discovery_paired_non_flu_virus_aucs <- test_individual_genes_on_datasets(MAGICAL_multiome_14_genes, non_flu_virus_discovery_list, "Multiome_Paired", "non_flu_virus_discovery")
multiome_discovery_all_non_flu_virus_aucs <- test_individual_genes_on_datasets(MAGICAL_multiome_19_genes, non_flu_virus_discovery_list, "Multiome_All", "non_flu_virus_discovery")

# Calculate bacteria (discovery) AUCs for individual genes
sc_discovery_bacteria_aucs <- test_individual_genes_on_datasets(MAGICAL_single_cell_genes, bacteria_discovery_list, "Single_Cell_Paired", "bacteria_discovery")
multiome_discovery_paired_bacteria_aucs <- test_individual_genes_on_datasets(MAGICAL_multiome_14_genes, bacteria_discovery_list, "Multiome_Paired", "bacteria_discovery")
multiome_discovery_all_bacteria_aucs <- test_individual_genes_on_datasets(MAGICAL_multiome_19_genes, bacteria_discovery_list, "Multiome_All", "bacteria_discovery")

# Calculate non-infectious (discovery) AUCs for individual genes
sc_discovery_noninfectious_aucs <- test_individual_genes_on_datasets(MAGICAL_single_cell_genes, noninfectious_discovery_list, "Single_Cell_Paired", "noninfectious_discovery")
multiome_discovery_paired_noninfectious_aucs <- test_individual_genes_on_datasets(MAGICAL_multiome_14_genes, noninfectious_discovery_list, "Multiome_Paired", "noninfectious_discovery")
multiome_discovery_all_noninfectious_aucs <- test_individual_genes_on_datasets(MAGICAL_multiome_19_genes, noninfectious_discovery_list, "Multiome_All", "noninfectious_discovery")

# Combine results from the different gene lists into comprehensive discovery dataframes
all_discovery_flu_aucs <- rbind(sc_discovery_flu_aucs, multiome_discovery_paired_flu_aucs, multiome_discovery_all_flu_aucs)
all_discovery_non_flu_virus_aucs <- rbind(sc_discovery_non_flu_virus_aucs, multiome_discovery_paired_non_flu_virus_aucs, multiome_discovery_all_non_flu_virus_aucs)
all_discovery_bacteria_aucs <- rbind(sc_discovery_bacteria_aucs, multiome_discovery_paired_bacteria_aucs, multiome_discovery_all_bacteria_aucs)
all_discovery_noninfectious_aucs <- rbind(sc_discovery_noninfectious_aucs, multiome_discovery_paired_noninfectious_aucs, multiome_discovery_all_noninfectious_aucs)

# Create final, overall discovery dataframe
all_discovery_aucs <- all_discovery_flu_aucs
all_discovery_aucs$non_flu_virus_discovery_gene_auc <- all_discovery_non_flu_virus_aucs$non_flu_virus_discovery_gene_auc
all_discovery_aucs$bacteria_discovery_gene_auc <- all_discovery_bacteria_aucs$bacteria_discovery_gene_auc
all_discovery_aucs$noninfectious_discovery_gene_auc <- all_discovery_noninfectious_aucs$noninfectious_discovery_gene_auc
all_discovery_aucs <- all_discovery_aucs[,c(1,2,4,5,6,3)]

# Plot non-influenza virus AUC vs influenza AUC
# NOTE - check for overlap between different data types (e.g., noninfectious vs flu)
# NOTE - one row is removed because we got an NA value - this is probably fine
# No labels
ggplot(all_discovery_aucs, aes(x=flu_discovery_gene_auc, y=non_flu_virus_discovery_gene_auc, group = source)) + geom_point(aes(col = source)) + 
  coord_fixed(xlim = c(0,1), ylim = c(0,1)) + geom_vline(xintercept=0.3, linetype=2) + geom_vline(xintercept=0.7, linetype=2) + 
  geom_abline(slope = 1, intercept = 0.075, linetype=2) + geom_abline(slope = 1, intercept = -0.075, linetype=2)

# AUC > 0.7 genes
ggplot(all_discovery_aucs, aes(x=flu_discovery_gene_auc, y=non_flu_virus_discovery_gene_auc, group = source, label = gene_name)) + geom_point(aes(col = source)) + 
  coord_fixed(xlim = c(0,1), ylim = c(0,1)) + geom_vline(xintercept=0.3, linetype=2) + geom_vline(xintercept=0.7, linetype=2) + 
  geom_abline(slope = 1, intercept = 0.075, linetype=2) + geom_abline(slope = 1, intercept = -0.075, linetype=2) + geom_text(aes(label=ifelse(flu_discovery_gene_auc>0.7,as.character(gene_name),'')))
current_df <- all_discovery_aucs[all_discovery_aucs$flu_discovery_gene_auc > 0.7,]
for(current_row_index in 1:nrow(current_df)) {
  current_row <- current_df[current_row_index,]
  line_coord <- current_row$flu_discovery_gene_auc - 0.075
  if(current_row$non_flu_virus_discovery_gene_auc < line_coord) {
    print(current_row$gene_name)
  }
}

# AUC < 0.3 genes
ggplot(all_discovery_aucs, aes(x=flu_discovery_gene_auc, y=non_flu_virus_discovery_gene_auc, group = source, label = gene_name)) + geom_point(aes(col = source)) + 
  coord_fixed(xlim = c(0,1), ylim = c(0,1)) + geom_vline(xintercept=0.3, linetype=2) + geom_vline(xintercept=0.7, linetype=2) + 
  geom_abline(slope = 1, intercept = 0.075, linetype=2) + geom_abline(slope = 1, intercept = -0.075, linetype=2) + geom_text(aes(label=ifelse(flu_discovery_gene_auc<0.3,as.character(gene_name),'')))
all_discovery_aucs[all_discovery_aucs$flu_discovery_gene_auc < 0.3,]
current_df <- all_discovery_aucs[all_discovery_aucs$flu_discovery_gene_auc < 0.3,]
for(current_row_index in 1:nrow(current_df)) {
  current_row <- current_df[current_row_index,]
  line_coord <- current_row$flu_discovery_gene_auc + 0.075
  if(current_row$non_flu_virus_discovery_gene_auc > line_coord) {
    print(current_row$gene_name)
  }
}

# Plot non-influenza virus AUC vs influenza AUC for viral specificity
ggplot(all_discovery_aucs, aes(x=flu_discovery_gene_auc, y=non_flu_virus_discovery_gene_auc, group = source)) + geom_point(aes(col = source)) + 
  coord_fixed(xlim = c(0,1), ylim = c(0,1)) + geom_hline(yintercept=0.3, linetype=2) + geom_hline(yintercept=0.7, linetype=2) + 
  geom_abline(slope = 1, intercept = 0.075, linetype=2) + geom_abline(slope = 1, intercept = -0.075, linetype=2)

# AUC > 0.7 genes
ggplot(all_discovery_aucs, aes(x=flu_discovery_gene_auc, y=non_flu_virus_discovery_gene_auc, group = source, label = gene_name)) + geom_point(aes(col = source)) + 
  coord_fixed(xlim = c(0,1), ylim = c(0,1)) + geom_hline(yintercept=0.3, linetype=2) + geom_hline(yintercept=0.7, linetype=2) + 
  geom_abline(slope = 1, intercept = 0.075, linetype=2) + geom_abline(slope = 1, intercept = -0.075, linetype=2) + geom_text(aes(label=ifelse(non_flu_virus_discovery_gene_auc>0.7,as.character(gene_name),'')))
all_discovery_aucs[all_discovery_aucs$non_flu_virus_discovery_gene_auc > 0.7,]
current_df <- na.omit(all_discovery_aucs[all_discovery_aucs$non_flu_virus_discovery_gene_auc > 0.7,])
for(current_row_index in 1:nrow(current_df)) {
  current_row <- current_df[current_row_index,]
  line_coord <- current_row$flu_discovery_gene_auc + 0.075
  if(current_row$non_flu_virus_discovery_gene_auc > line_coord) {
    print(current_row$gene_name)
  }
}



# AUC < 0.3 genes
ggplot(all_discovery_aucs, aes(x=flu_discovery_gene_auc, y=non_flu_virus_discovery_gene_auc, group = source, label = gene_name)) + geom_point(aes(col = source)) + 
  coord_fixed(xlim = c(0,1), ylim = c(0,1)) + geom_hline(yintercept=0.3, linetype=2) + geom_hline(yintercept=0.7, linetype=2) + 
  geom_abline(slope = 1, intercept = 0.075, linetype=2) + geom_abline(slope = 1, intercept = -0.075, linetype=2) + geom_text(aes(label=ifelse(non_flu_virus_discovery_gene_auc<0.3,as.character(gene_name),'')))
all_discovery_aucs[all_discovery_aucs$non_flu_virus_discovery_gene_auc < 0.3,]
current_df <- na.omit(all_discovery_aucs[all_discovery_aucs$non_flu_virus_discovery_gene_auc < 0.3,])
for(current_row_index in 1:nrow(current_df)) {
  current_row <- current_df[current_row_index,]
  line_coord <- current_row$flu_discovery_gene_auc - 0.075
  if(current_row$non_flu_virus_discovery_gene_auc < line_coord) {
    print(current_row$gene_name)
  }
}

# Plot bacteria AUC vs influenza AUC for bacterial specificity
ggplot(all_discovery_aucs, aes(x=flu_discovery_gene_auc, y=bacteria_discovery_gene_auc , group = source)) + geom_point(aes(col = source)) + 
  coord_fixed(xlim = c(0,1), ylim = c(0,1)) + geom_hline(yintercept=0.3, linetype=2) + geom_hline(yintercept=0.7, linetype=2) + 
  geom_abline(slope = 1, intercept = 0.125, linetype=2) + geom_abline(slope = 1, intercept = -0.125, linetype=2)

# AUC > 0.7 genes
ggplot(all_discovery_aucs, aes(x=flu_discovery_gene_auc, y=bacteria_discovery_gene_auc, group = source, label = gene_name)) + geom_point(aes(col = source)) + 
  coord_fixed(xlim = c(0,1), ylim = c(0,1)) + geom_hline(yintercept=0.3, linetype=2) + geom_hline(yintercept=0.7, linetype=2) + 
  geom_abline(slope = 1, intercept = 0.125, linetype=2) + geom_abline(slope = 1, intercept = -0.125, linetype=2) + geom_text(aes(label=ifelse(bacteria_discovery_gene_auc>0.7,as.character(gene_name),'')))
all_discovery_aucs[all_discovery_aucs$bacteria_discovery_gene_auc > 0.7,]
current_df <- na.omit(all_discovery_aucs[all_discovery_aucs$bacteria_discovery_gene_auc > 0.7,])
for(current_row_index in 1:nrow(current_df)) {
  current_row <- current_df[current_row_index,]
  line_coord <- current_row$flu_discovery_gene_auc + 0.125
  if(current_row$bacteria_discovery_gene_auc > line_coord) {
    print(paste0(current_row$gene_name, " - ", current_row$source))
  }
}

# AUC < 0.3 genes
ggplot(all_discovery_aucs, aes(x=flu_discovery_gene_auc, y=bacteria_discovery_gene_auc, group = source, label = gene_name)) + geom_point(aes(col = source)) + 
  coord_fixed(xlim = c(0,1), ylim = c(0,1)) + geom_hline(yintercept=0.3, linetype=2) + geom_hline(yintercept=0.7, linetype=2) + 
  geom_abline(slope = 1, intercept = 0.125, linetype=2) + geom_abline(slope = 1, intercept = -0.125, linetype=2) + geom_text(aes(label=ifelse(bacteria_discovery_gene_auc<0.3,as.character(gene_name),'')))
all_discovery_aucs[all_discovery_aucs$bacteria_discovery_gene_auc < 0.3,]
current_df <- na.omit(all_discovery_aucs[all_discovery_aucs$bacteria_discovery_gene_auc < 0.3,])
for(current_row_index in 1:nrow(current_df)) {
  current_row <- current_df[current_row_index,]
  line_coord <- current_row$flu_discovery_gene_auc - 0.125
  if(current_row$bacteria_discovery_gene_auc < line_coord) {
    print(paste0(current_row$gene_name, " - ", current_row$source))
  }
}

# Plot noninfectious AUC vs influenza AUC for non-infectious specificity
ggplot(all_discovery_aucs, aes(x=flu_discovery_gene_auc, y=noninfectious_discovery_gene_auc, group = source)) + geom_point(aes(col = source)) + 
  coord_fixed(xlim = c(0,1), ylim = c(0,1)) + geom_hline(yintercept=0.3, linetype=2) + geom_hline(yintercept=0.7, linetype=2) + 
  geom_abline(slope = 1, intercept = 0.075, linetype=2) + geom_abline(slope = 1, intercept = -0.075, linetype=2)