library(MetaIntegrator)
library(dplyr)
library(openxlsx)

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

# NON INFLUENZA VIRUS
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

non_flu_virus_Daniel_final_sample_sizes <- c()
for(dataset_id in non_flu_virus_Daniel_metadata$Study) {
  for(non_flu_virus_Daniel_dataset in non_flu_virus_Daniel_list) {
    if(non_flu_virus_Daniel_dataset$formattedName == dataset_id) {
      non_flu_virus_Daniel_final_sample_sizes <- c(non_flu_virus_Daniel_final_sample_sizes, nrow(non_flu_virus_Daniel_dataset$pheno))
    }
  }
}
non_flu_virus_Daniel_metadata$final_sample_size <- non_flu_virus_Daniel_final_sample_sizes

# INFLUENZA VIRUS
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

# Need to create non_influenza_discovery_list, bacteria_discovery_list, and non_infectious_discovery_list




sc_flu_aucs <- test_individual_genes_on_datasets(MAGICAL_single_cell_genes, flu_discovery_list, "Single_Cell_Paired")
multiome_paired_flu_aucs <- test_individual_genes_on_datasets(MAGICAL_multiome_14_genes, flu_discovery_list, "Multiome_Paired")
all_flu_aucs <- rbind(sc_flu_aucs, multiome_paired_flu_aucs)

# Create gene set signature (using MetaIntegrator format)
sig <- list()
sig$posGeneNames <- '7124'
sig$negGeneNames <- ''
sig$filterDescription <- 'test'
sig$FDRThresh <- 0
sig$effectSizeThresh <- 0
sig$numberStudiesThresh <- 1
sig$isLeaveOneOut <- F
sig$heterogeneityPvalThresh <- 0
sig$timestamp <- Sys.time()

# Test gene set signature on influenza samples and print median AUROC
flu_aucs <- sapply(X = flu_discovery_list, FUN = calculateAUROC, signature = sig)
print(median(flu_aucs, na.rm = TRUE))

non_flu_virus_aucs <- sapply(X = non_flu_virus_list, FUN = calculateAUROC, signature = sig)
print(median(non_flu_virus_aucs, na.rm = TRUE))

bacteria_aucs <- sapply(X = bacteria_list, FUN = calculateAUROC, signature = sig)
print(median(bacteria_aucs, na.rm = TRUE))

noninfectious_aucs <- sapply(X = noninfectious_list, FUN = calculateAUROC, signature = sig)
print(median(noninfectious_aucs, na.rm = TRUE))







sapply(flu_list, countSampleSizes)

chosen_flu_dataset_sizes <- c(159, 101, 197, 39, 42, 19, 97)
chosen_flu_indices <- c(1,3,4,9,11,13,15)

discovery_influenza_list <- flu_list[chosen_flu_indices]