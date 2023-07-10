library(data.table)
library(DESeq2)
library(MetaIntegrator)

base_dir <- "~/GitHub/Influenza/"
source(paste0(base_dir, "bulk_RNA_analysis_helper.R"))
source(paste0(base_dir, "pseudobulk_analysis_helper.R"))
source(paste0(base_dir, "Data Compendium/Compendium_Functions.R"))
setup_bulk_analysis()
sample_metadata <- read.table(paste0(base_dir, "all_metadata_sheet.tsv"), sep = "\t", header = TRUE)
cell_types <- c("CD4_Naive", "CD8_Naive", "CD4_Memory", "CD8_Memory", "cDC", "HSPC", "pDC", "Platelet", "Plasmablast", "Proliferating", "NK", "T_Naive", "CD14_Mono", "CD16_Mono", "MAIT")
single_cell_magical_dir <- "C:/Users/willi/OneDrive - Princeton University/Influenza Analysis/Single Cell RNA-Seq/MAGICAL Analyses/Placebo 4 Subject HVL (SPEEDI) - SCT/"
single_cell_pseudobulk_dir <- paste0(single_cell_magical_dir, "scRNA_pseudobulk/")
multiome_magical_dir <- "C:/Users/willi/OneDrive - Princeton University/Influenza Analysis/True Multiome/MAGICAL Analyses/14 Placebo Sample (Final)/"
multiome_pseudobulk_dir <- paste0(multiome_magical_dir, "scRNA_pseudobulk/")
set.seed(2000)