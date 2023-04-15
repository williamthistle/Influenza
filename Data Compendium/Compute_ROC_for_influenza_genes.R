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

# All infectious datasets are found in data_list
data_list <- readRDS('Processed, Manuscript Files/data_list_v1.RDS')
# We have to load noninfectious datasets separately
noninfectious_list <- readRDS('Processed, Manuscript Files/noninfectious_data_list.RDS')

# Find influenza datasets 
flu_samples <- sapply(data_list, findPathogen, pathogen = 'Influenza virus')
flu_list <- data_list[flu_samples]
flu_list <- lapply(flu_list, formatPathogenContrast, pathogen = 'Influenza virus')
filtered_samples <- sapply(flu_list, filterSampleSizes, N = 4)
flu_list <- flu_list[filtered_samples]

# Find non-influenza virus datasets
non_flu_virus_list <- lapply(data_list, removePathogen, pathogen = 'Influenza virus', target_class = "Virus")
filtered_samples <- sapply(non_flu_virus_list, filterSampleSizes, N = 4)
non_flu_virus_list <- non_flu_virus_list[filtered_samples]

# Find bacteria datasets
bacteria_list <- lapply(data_list, removePathogen, pathogen = 'Influenza virus', target_class = "Bacteria")
filtered_samples <- sapply(bacteria_list, filterSampleSizes, N = 4)
bacteria_list <- bacteria_list[filtered_samples]

# Find non-infectious datasets
noninfectious_list <- lapply(noninfectious_list, formatContrast, target_class = "Other.NonInfectious")
filtered_samples <- sapply(noninfectious_list, filterSampleSizes, N = 4)
noninfectious_list <- noninfectious_list[filtered_samples]

# Create gene set signature (using MetaIntegrator format)
sig <- list()
sig$posGeneNames <- '6280'
sig$negGeneNames <- ''
sig$filterDescription <- 'test'
sig$FDRThresh <- 0
sig$effectSizeThresh <- 0
sig$numberStudiesThresh <- 1
sig$isLeaveOneOut <- F
sig$heterogeneityPvalThresh <- 0
sig$timestamp <- Sys.time()

# Test gene set signature on influenza samples and print median AUROC
flu_aucs <- sapply(X = flu_list, FUN = calculateAUROC, signature = sig)
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