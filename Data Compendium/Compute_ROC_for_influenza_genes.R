library(MetaIntegrator)
library(dplyr)
library(openxlsx)

setwd("C:/Users/willi/Desktop/Influenza Work (April 2023)/Data Compendium")
source("Compendium_Functions.R")

study_metadata <- read.xlsx('metadata_tables/darpa_compendium_metadata_v1.xlsx', sheet = 1)
printable_df <- study_metadata %>%
  dplyr::select(Accession, Platform, Standardized.Exposure, Classes, Control.Type, Age)

data_list <- readRDS('Processed, Manuscript Files/data_list_v1.RDS')


flu_samples <- sapply(data_list, findPathogen, pathogen = 'Influenza virus')
flu_list <- data_list[flu_samples]

flu_list <- lapply(flu_list, formatPathogenContrast, pathogen = 'Influenza virus')
# filter out studies that do not have at least N cases and N controls
flu_N = 4
filtered_samples <- sapply(flu_list, filterSampleSizes, N = flu_N)
flu_list <- flu_list[filtered_samples]

sapply(flu_list, countSampleSizes)

chosen_flu_dataset_sizes <- c(159, 101, 197, 39, 42, 19, 97)
chosen_flu_indices <- c(1,3,4,9,11,13,15)

discovery_influenza_list <- flu_list[chosen_flu_indices]


sig <- list()
sig$posGeneNames <- '8840'
sig$negGeneNames <- ''
sig$filterDescription <- 'test'
sig$FDRThresh <- 0
sig$effectSizeThresh <- 0
sig$numberStudiesThresh <- 1
sig$isLeaveOneOut <- F
sig$heterogeneityPvalThresh <- 0
sig$timestamp <- Sys.time()

flu_aucs <- sapply(X = flu_list, FUN = calculateAUROC, signature = sig)
print(median(flu_aucs, na.rm = TRUE))
