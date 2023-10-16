library(SPEEDI)
library(Seurat)
library(parallel)
library(doMC)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(openxlsx)
library(future)
# Load extra RNA functions
home_dir <- "~/"
source(paste0(home_dir, "extra_functions/rna/preprocessing_and_qc.R"))
source(paste0(home_dir, "extra_functions/rna/get_stats.R"))
source(paste0(home_dir, "extra_functions/rna/manipulate_data.R"))
source(paste0(home_dir, "extra_functions/rna/differential_expression.R"))
source(paste0(home_dir, "extra_functions/rna/visualization.R"))
source(paste0(home_dir, "extra_functions/atac/preprocessing_and_qc.R"))
source(paste0(home_dir, "extra_functions/atac/processing.R"))

date <- Sys.Date()
hostname <- system("hostname", intern = TRUE)
if(hostname == "lumos.Princeton.EDU") {
  addArchRThreads(16)
}

addArchRGenome("hg38")