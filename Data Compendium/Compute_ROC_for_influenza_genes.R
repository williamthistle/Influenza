library(MetaIntegrator)
library(dplyr)
library(openxlsx)
library(ggplot2)

# Set working directory to wherever this script is
# (because Compendium_Functions.R will be in the same place)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("Compendium_Functions.R")
# Now set working directory to where compendium data files are
setwd("~/")
setwd("../OneDrive - Princeton University/Influenza Analysis/Data Compendium")

# Grab study metadata (will be used to confirm that discovery / validation splits are balanced)
study_metadata <- read.xlsx('metadata_tables/darpa_compendium_metadata_v1.xlsx', sheet = 1)
# Grab the accessions from Daniel's publication for influenza and non-influenza virus datasets
daniel_paper_datasets <- read.xlsx('metadata_tables/Daniel_Paper_Datasets.xlsx', sheet = 1)

# All infectious datasets are loaded into data_list
data_list <- readRDS('Processed, Manuscript Files/data_list_v1.RDS')
data_list <- fix_data_list(data_list)

# We have to load noninfectious datasets separately
noninfectious_list <- readRDS('Processed, Manuscript Files/noninfectious_data_list.RDS')

# Read in MAGICAL gene lists
# TODO: Add single cell RNA-seq gene lists (that pass pseudobulk)

# Single-cell RNA-seq MAGICAL genes (paired)
MAGICAL_single_cell_genes <- read.table("../Single Cell RNA-Seq/MAGICAL_flu_genes.csv")$V1

# True multiome RNA-seq MAGICAL genes (paired)
MAGICAL_multiome_14_genes <- read.table("../True Multiome/MAGICAL_flu_genes_14.csv")$V1

# True multiome RNA-seq MAGICAL genes (all)
MAGICAL_multiome_19_genes <- read.table("../True Multiome/MAGICAL_flu_genes_19.csv")$V1

# MintCHiP genes (all) - probably not analyzed here
mintchip_genes <- read.table("../MintChIP/mintchip_gene_list.csv")$V1

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
flu_discovery_info <- split_flu_data(flu_list, study_metadata, daniel_paper_datasets, "Discovery")
flu_discovery_list <- flu_discovery_info[[1]]
flu_discovery_metadata <- flu_discovery_info[[2]]

# Find influenza validation datasets
flu_validation_info <- split_flu_data(flu_list, study_metadata, daniel_paper_datasets, "Validation (robustness)")
flu_validation_list <- flu_validation_info[[1]]
flu_validation_metadata <- flu_validation_info[[2]]

# Fix time series data for influenza discovery datasets
flu_discovery_list <- fix_time_series_for_flu_discovery(flu_discovery_list)

# Append final sample size to influenza discovery metadata
flu_discovery_metadata <- append_final_sample_size(flu_discovery_metadata, flu_discovery_list)

# # Fix time series data for influenza validation datasets
flu_validation_list <- fix_time_series_for_flu_validation(flu_validation_list)

# Append final sample size to influenza validation metadata
flu_validation_metadata <- append_final_sample_size(flu_validation_metadata, flu_validation_list)

# Find non-influenza virus Daniel datasets
non_flu_virus_Daniel_info <- split_flu_data(non_flu_virus_list, study_metadata, daniel_paper_datasets, "Validation (cross-reactivity)")
non_flu_virus_Daniel_list <- non_flu_virus_Daniel_info[[1]]
non_flu_virus_Daniel_metadata <- non_flu_virus_Daniel_info[[2]]

# Fix time series data for non influenza virus datasets
non_flu_virus_Daniel_list <- fix_time_series_for_non_flu_virus(non_flu_virus_Daniel_list)

# Add sample sizes and figure out which datasets overlap with influenza discovery / validation
non_flu_virus_Daniel_metadata <- append_final_sample_size(non_flu_virus_Daniel_metadata, non_flu_virus_Daniel_list)
non_flu_virus_Daniel_metadata <- append_overlap_with_other_metadata(non_flu_virus_Daniel_metadata, flu_discovery_metadata)

# create discovery and validation datasets for non-influenza virus
non_flu_discovery_indices <- c(3, 4, 5, 6, 10, 11, 12)
non_flu_validation_indices <- c(1, 2, 7, 8, 9)
non_flu_virus_discovery_list <- non_flu_virus_Daniel_list[non_flu_discovery_indices]
non_flu_virus_validation_list <- non_flu_virus_Daniel_list[non_flu_validation_indices]

# BACTERIA
load("bacteria_discovery_metadata.rda")
load("bacteria_validation_metadata.rda")

bacteria_dataset_names <- c()
for(current_bacteria_dataset in bacteria_list) {
  current_accession <- current_bacteria_dataset$formattedName
  bacteria_dataset_names <- c(bacteria_dataset_names, current_accession)
}

bacteria_metadata <- study_metadata[study_metadata$Study %in% bacteria_dataset_names,]
bacteria_discovery_list <- bacteria_list[sort(match(bacteria_discovery_metadata$Study, bacteria_metadata$Study))]
bacteria_validation_list <- bacteria_list[sort(match(bacteria_validation_metadata$Study, bacteria_metadata$Study))]

# NON INFECTIOUS
load("noninfectious_discovery_metadata.rda")
load("noninfectious_validation_metadata.rda")

noninfectious_dataset_names <- c()
for(current_noninfectious_dataset in noninfectious_list) {
  current_accession <- current_noninfectious_dataset$formattedName
  noninfectious_dataset_names <- c(noninfectious_dataset_names, current_accession)
}

# TODO: Check overlap between noninfectious and other stuff (want to make sure discovery stuff is together and validation is together)
noninfectious_discovery_list <- noninfectious_list[sort(match(noninfectious_discovery_metadata$Study, noninfectious_dataset_names))]
noninfectious_validation_list <- noninfectious_list[sort(match(noninfectious_validation_metadata$Study, noninfectious_dataset_names))]

# Get discovery AUCs
all_discovery_aucs <- get_discovery_aucs(MAGICAL_single_cell_genes, MAGICAL_multiome_14_genes, MAGICAL_multiome_19_genes, flu_discovery_list, non_flu_virus_discovery_list, 
                                         bacteria_discovery_list, noninfectious_discovery_list)

# Gene lists to store positive and negative genes for gene signature
flu_pos_genes <- c()
flu_neg_genes <- c()
vir_pos_genes <- c()
vir_neg_genes <- c()
bac_pos_genes <- c()
bac_neg_genes <- c()
nif_pos_genes <- c()
nif_neg_genes <- c()

# Gather flu positive and negative genes
flu_pos_and_neg_info <- capture_flu_pos_and_neg(all_discovery_aucs)
flu_pos_df <- flu_pos_and_neg_info[[1]]
flu_neg_df <- flu_pos_and_neg_info[[2]]
flu_pos_genes <- flu_pos_and_neg_info[[3]]
flu_neg_genes <- flu_pos_and_neg_info[[4]]

# Gather non-flu virus positive and negative genes
vir_pos_and_neg_info <- capture_vir_pos_and_neg(all_discovery_aucs)
vir_pos_df <- vir_pos_and_neg_info[[1]]
vir_neg_df <- vir_pos_and_neg_info[[2]]
vir_pos_genes <- vir_pos_and_neg_info[[3]]
vir_neg_genes <- vir_pos_and_neg_info[[4]]

# Plot bacteria AUC vs influenza AUC for bacterial specificity
bac_pos_and_neg_info <- capture_bac_pos_and_neg(all_discovery_aucs)
bac_pos_df <- bac_pos_and_neg_info[[1]]
bac_neg_df <- bac_pos_and_neg_info[[2]]
bac_pos_genes <- bac_pos_and_neg_info[[3]]
bac_neg_genes <- bac_pos_and_neg_info[[4]]

# Plot noninfectious AUC vs influenza AUC for non-infectious specificity
# No genes pass criteria
nif_pos_and_neg_info <- capture_nif_pos_and_neg(all_discovery_aucs)
nif_pos_df <- nif_pos_and_neg_info[[1]]
nif_neg_df <- nif_pos_and_neg_info[[2]]
nif_pos_genes <- nif_pos_and_neg_info[[3]]
nif_neg_genes <- nif_pos_and_neg_info[[4]]




# Create final gene signature

# Using difference in AUC to create a smaller sig
flu_pos_genes <- c("ELF1", "CAPN2", "RAB8B", "ARIH1", "MEF2A")
flu_neg_genes <- c("AUTS2", "USP36", "LPCAT1", "SECISBP2", "UQCR11", "ETS1", "HLA-DRA", "PRKCA", "HNRNPDL", "RASSF1", "CD247", "ZNF831", "BRD1")
vir_pos_genes <- c("IRAK3")
vir_neg_genes <- c("SLC38A1")
bac_pos_genes <- c("ADGRE5", "TUBA1A", "LRRK2", "SIRPA", "HCAR3")
bac_neg_genes <- c("CD69", "EZR", "CEBPZ", "SERBP1", "SLC38A1", "PPP3CC", "PRMT1", "NCL")

# My signature
sig_pos_genes <- c(flu_pos_genes, vir_neg_genes, bac_neg_genes, nif_neg_genes)
sig_neg_genes <- c(flu_neg_genes, vir_pos_genes, bac_pos_genes, nif_pos_genes)

# Daniel's signature
old_flu_pos_genes <- c("CAPN2", "CITED2", "DUSP1", "ELF1", "IER2", "KLF6", "MAP3K8", "MOB1A", "RAB8B", "SAMHD1", "USP3", "USP8")
old_flu_neg_genes <- c("HNRNPDL", "KRT10", "OAZ2", "OST4", "PRMT1", "TOMM6", "UQCR11")
old_vir_pos_genes <- c("S100A4", "UBE2J1")
old_vir_neg_genes <- c("JUN")
old_bac_pos_genes <- c("APLP2", "HCAR3", "IL1RAP", "IQSEC1", "LRRK2", "VCAN")
old_bac_neg_genes <- c("CD69", "CEBPZ", "SERBP1")
old_nif_pos_genes <- c()
old_nif_neg_genes <- c("CHURC1")

old_sig_pos_genes <- c(old_flu_pos_genes, old_vir_neg_genes, old_bac_neg_genes, old_nif_neg_genes)
old_sig_neg_genes <- c(old_flu_neg_genes, old_vir_pos_genes, old_bac_pos_genes, old_nif_pos_genes)

# Signature 1 (28619954)
other_sig_1_pos_genes <- c("IFI27")
other_sig_1_neg_genes <- c()

# Signature 2 (26682989)
other_sig_2_pos_genes <- c("CD38", "HERC5", "HERC6", "IFI6", "IFIH1", "LGALS3BP", "LY6E", "MX1", "PARP12", "RTP4", "ZBP1")
other_sig_2_neg_genes <- c()


plot_pooled_auc(sig_pos_genes, sig_neg_genes, flu_validation_list, "Influenza Pooled AUC (New)")
plot_pooled_auc(old_sig_pos_genes, old_sig_neg_genes, flu_validation_list, "Influenza Pooled AUC (Old)")
plot_pooled_auc(other_sig_1_pos_genes, other_sig_1_neg_genes, flu_validation_list, "Influenza Pooled AUC (Other Sig 1)")
plot_pooled_auc(other_sig_2_pos_genes, other_sig_2_neg_genes, flu_validation_list, "Influenza Pooled AUC (Other Sig 2)")

plot_pooled_auc(sig_pos_genes, sig_neg_genes, non_flu_virus_validation_list, "Non-Influenza Virus Pooled AUC (New)")
plot_pooled_auc(old_sig_pos_genes, old_sig_neg_genes, non_flu_virus_validation_list, "Non-Influenza Virus Pooled AUC (Old)")
plot_pooled_auc(other_sig_1_pos_genes, other_sig_1_neg_genes, non_flu_virus_validation_list, "Non-Influenza Virus Pooled AUC (Other Sig 1)")
plot_pooled_auc(other_sig_2_pos_genes, other_sig_2_neg_genes, non_flu_virus_validation_list, "Non-Influenza Virus Pooled AUC (Other Sig 2)")

plot_pooled_auc(sig_pos_genes, sig_neg_genes, bacteria_validation_list, "Bacteria Pooled AUC (New)")
plot_pooled_auc(old_sig_pos_genes, old_sig_neg_genes, bacteria_validation_list, "Bacteria Pooled AUC (Old)")

plot_pooled_auc(sig_pos_genes, sig_neg_genes, noninfectious_validation_list, "Noninfectious Pooled AUC (New)")
plot_pooled_auc(old_sig_pos_genes, old_sig_neg_genes, noninfectious_validation_list, "Noninfectious Pooled AUC (Old)")








