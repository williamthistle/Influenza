library(SPEEDI)
library(Seurat)
library(parallel)
library(doMC)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(openxlsx)
library(data.table)
library(DESeq2)
# library(MetaIntegrator)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(DiffBind)
library(pheatmap)
library(ArchRtoSignac)
library(EnsDb.Hsapiens.v86)
library(chromVARmotifs)
library(Signac)
# Load extra RNA functions
home_dir <- "~/"
source(paste0(home_dir, "extra_functions/rna/preprocessing_and_qc.R"))
source(paste0(home_dir, "extra_functions/rna/get_stats.R"))
source(paste0(home_dir, "extra_functions/rna/manipulate_data.R"))
source(paste0(home_dir, "extra_functions/rna/differential_expression.R"))
source(paste0(home_dir, "extra_functions/rna/visualization.R"))
source(paste0(home_dir, "extra_functions/atac/preprocessing_and_qc.R"))
source(paste0(home_dir, "extra_functions/atac/processing.R"))

# Load motif data (CIS-BP V2)
data("human_pwms_v2")

# Load mintchip metadata
mintchip_metadata <- read.table(paste0(home_dir, "mintchip_metadata.tsv"), sep = "\t", header = TRUE)

date <- Sys.Date()
addArchRThreads(16)

addArchRGenome("hg38")