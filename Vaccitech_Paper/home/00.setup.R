library(data.table)
library(DESeq2)
library(MetaIntegrator)
library(openxlsx)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPseeker)
library(DiffBind)
library(SPEEDI)
library(pheatmap)
library(dplyr)
library(GenomicRanges)
library(SummarizedExperiment)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
library(monaLisa)
library(ComplexHeatmap)
library(circlize)
library(chromVARmotifs)
library(sva)
library(corrr)
library(enrichR)
library(smplot2)

# Load bulk RNA-seq analysis results
onedrive_dir <- "~"
setwd(onedrive_dir)
setwd("../")
laptop_use <- grep("OneDrive", getwd())
if(length(laptop_use == 1)) {
  setwd("../")
}
onedrive_dir <- getwd()
onedrive_dir <- paste0(onedrive_dir, "/OneDrive - Princeton University/")
load(paste0(onedrive_dir, "Vaccitech_Paper/Analyses/Bulk_RNA-Seq/Results/bulk_RNA_analysis.RData"))

# Reload onedrive dir
onedrive_dir <- "~"
setwd(onedrive_dir)
setwd("../")
laptop_use <- grep("OneDrive", getwd())
if(length(laptop_use == 1)) {
  setwd("../")
}
onedrive_dir <- getwd()
onedrive_dir <- paste0(onedrive_dir, "/OneDrive - Princeton University/")

# Misc settings
set.seed(get_speedi_seed())
options(max.print=10000)

# Gene databases for hg19 (MintChIP) and hg38 (everything else)
txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene
txdb_hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Set base script dir 
script_base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"

# Source various scripts
source(paste0(script_base_dir, "extra_functions/bulk_RNA_analysis_helper.R"))
source(paste0(script_base_dir, "extra_functions/Compendium_Functions.R"))
source(paste0(script_base_dir, "extra_functions/humanbase_functions.R"))
source(paste0(script_base_dir, "extra_functions/pseudobulk_analysis_helper.R"))
source(paste0(script_base_dir, "extra_functions/MAGICAL_functions.R"))

# Set overall metadata dir
overall_metadata_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/Overall_Metadata/")

# Set bulk RNA-seq dirs
bulk_rna_data_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/Bulk_RNA-Seq/Data/")
bulk_rna_results_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/Bulk_RNA-Seq/Results/")

# Read in bulk RNA-seq cell type proportions
cibersort_cell_type_proportions <- read.table(paste0(bulk_rna_data_dir, "cibertsortX_rsem_genes_count.processed.CPM.tsv"), sep = "\t", header = TRUE)
# Rename Mixture to aliquot_id for merging
cibersort_cell_type_proportions$aliquot_id <- cibersort_cell_type_proportions$Mixture
cibersort_cell_type_proportions <- cibersort_cell_type_proportions[,-c(1)]

setup_bulk_rna_analysis(metadata_dir = overall_metadata_dir, data_dir = bulk_rna_data_dir)

bulk_rna_cell_types <- c("B.cells.naive", "T.cells.CD8", "T.cells.CD4.naive", "T.cells.CD4.memory.resting", "T.cells.regulatory..Tregs.", "NK.cells.resting",
                         "Monocytes", "Mast.cells.resting", "Neutrophils")
innate_cell_types <- c("CD16 Mono","CD14 Mono","cDC","pDC","NK","NK_CD56bright")
adaptive_cell_types <- c("CD4 Naive", "CD8 Naive", "CD4 Memory", "CD8 Memory", "B memory", "MAIT", "B naive", "B")
scRNA_cell_types <- c("CD16 Mono","CD14 Mono","cDC","pDC","NK","NK_CD56bright", "CD4 Naive", "CD8 Naive", "CD4 Memory", "CD8 Memory", "B memory", "MAIT", "B naive")
snME_cell_types <- c("B-Mem", "B-Naive", "Monocyte", "NK-cell2", "Tc-Mem", "Tc-Naive", "Th-Mem", "Th-Naive")
scATAC_cell_types <- c("B", "CD4 Memory", "CD8 Memory", "CD14 Mono", "CD16 Mono", "NK", "CD4 Naive", "CD8 Naive", "pDC", "cDC", "MAIT", "Proliferating")
mintchip_markers <- c("H3K4me1", "H3K4me3", "H3K9me3", "H3K27Ac", "H3K27me3", "H3K36me3")

# scRNA-seq dirs
# HVL
scRNA_hvl_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/scRNA-Seq/Results/HVL/")
scRNA_hvl_placebo_deg_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/scRNA-Seq/Results/HVL/Placebo/DEGs/Normal/")
scRNA_hvl_placebo_MAGICAL_cell_types_deg_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/scRNA-Seq/Results/HVL/Placebo/DEGs/MAGICAL/")
scRNA_hvl_placebo_downstream_analysis_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/scRNA-Seq/Results/HVL/Placebo/Downstream_Analysis/")
scRNA_hvl_vaccinated_deg_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/scRNA-Seq/Results/HVL/Vaccinated/DEGs/Normal/")
scRNA_hvl_vaccinated_MAGICAL_cell_types_deg_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/scRNA-Seq/Results/HVL/Vaccinated/DEGs/MAGICAL/")
scRNA_hvl_vaccinated_downstream_analysis_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/scRNA-Seq/Results/HVL/Vaccinated/Downstream_Analysis/")

# LVL
scRNA_lvl_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/scRNA-Seq/Results/LVL/")
scRNA_lvl_placebo_deg_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/scRNA-Seq/Results/LVL/Placebo/DEGs/Normal/")
scRNA_lvl_placebo_MAGICAL_cell_types_deg_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/scRNA-Seq/Results/LVL/Placebo/DEGs/MAGICAL/")
scRNA_lvl_placebo_downstream_analysis_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/scRNA-Seq/Results/LVL/Placebo/Downstream_Analysis/")

# scATAC-seq dirs
# HVL
scATAC_hvl_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/scATAC-Seq/Results/HVL/")
scATAC_hvl_placebo_das_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/scATAC-Seq/Results/HVL/Placebo/DASs/")
scATAC_hvl_placebo_das_annotated_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/scATAC-Seq/Results/HVL/Placebo/DASs/Annotated/")
scATAC_hvl_vaccinated_das_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/scATAC-Seq/Results/HVL/Vaccinated/DASs/")
scATAC_hvl_vaccinated_das_annotated_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/scATAC-Seq/Results/HVL/Vaccinated/DASs/Annotated/")

# LVL
scATAC_lvl_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/scATAC-Seq/Results/LVL/")
scATAC_lvl_placebo_das_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/scATAC-Seq/Results/LVL/Placebo/DASs/")
scATAC_lvl_placebo_das_annotated_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/scATAC-Seq/Results/LVL/Placebo/DASs/Annotated/")

# MAGICAL analysis dirs
# HVL
MAGICAL_hvl_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/MAGICAL/Results/HVL/")
MAGICAL_hvl_placebo_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/MAGICAL/Results/HVL/Placebo/")
MAGICAL_hvl_placebo_output_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/MAGICAL/Results/HVL/Placebo/Output/")
MAGICAL_hvl_vaccinated_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/MAGICAL/Results/HVL/Vaccinated/")
MAGICAL_hvl_vaccinated_output_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/MAGICAL/Results/HVL/Vaccinated/Output/")

# LVL
MAGICAL_lvl_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/MAGICAL/Results/LVL/")
MAGICAL_lvl_placebo_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/MAGICAL/Results/LVL/Placebo/")
MAGICAL_lvl_placebo_output_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/MAGICAL/Results/LVL/Placebo/Output/")

# MintChIP analysis dirs
mintchip_metadata_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/MintChIP/Metadata/")
mintchip_das_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/MintChIP/Results/differentially_accessible_sites/")
homer_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/MintChIP/Results/HOMER/")

# snME analysis dirs
snME_data_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/snMethylation/Data/")
snME_results_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/snMethylation/Results/")

# miRNA analysis dirs
miRNA_data_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/miRNA/Data/")
miRNA_metadata_dir <-paste0(onedrive_dir, "Vaccitech_Paper/Analyses/miRNA/Metadata/")

# total RNA analysis dirs
totalRNA_data_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/Total_RNA/Data/")
totalRNA_metadata_dir <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/Total_RNA/Metadata/")

# Read in data lists for creating comprehensive metadata sheet
scRNA_data_list <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/scRNA-Seq/Metadata/scRNA_data_list.txt")
scATAC_data_list <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/scATAC-Seq/Metadata/snATAC_data_list.txt")
multiome_data_list <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/snMultiome/Metadata/multiome_data_list.txt")
bulkRNA_data_list <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/Bulk_RNA-Seq/Metadata/bulkRNA_data_list.txt")
mintchip_data_list <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/MintChIP/Metadata/mintchip_data_list.txt")
miRNA_data_list <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/miRNA/Metadata/miRNA_data_list.txt")
totalRNA_data_list <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/Total_RNA/Metadata/totalRNA_data_list.txt")
snme_data_list <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/snMethylation/Metadata/snME_data_list.txt")
bulk_methylation_metadata_file <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/Bulk_Methylation/Metadata/Bulk_Methylation_Sample_Metadata.csv")
scRNA_qc_file <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/scRNA-Seq/Metadata/ECHO_FLU_Vaccitech_PBMC_scrnaseq_coded_qc_report_WT.csv")
multiome_qc_file <- paste0(onedrive_dir, "Vaccitech_Paper/Analyses/snMultiome/Metadata/Stanford_FLU_combined_qc_metric_coded_09015022_qc_data.csv")
original_unfiltered_metadata_file <- paste0(metadata_dir, "Vaccitech_Original_Metadata_Sheet.csv")
scRNA_data <- read.table(scRNA_data_list)$V1
scRNA_data <- scRNA_data[1:length(scRNA_data) - 1]
scATAC_data <- read.table(scATAC_data_list)$V1
scATAC_data <- scATAC_data[1:length(scATAC_data) - 1]
multiome_data <- read.table(multiome_data_list)$V1
multiome_data <- multiome_data[1:length(multiome_data) - 1]
bulkRNA_data <- read.table(bulkRNA_data_list)$V1
mintchip_data <- read.table(mintchip_data_list)$V1
miRNA_data <- read.table(miRNA_data_list)$V1
totalRNA_data <- read.table(totalRNA_data_list)$V1
snme_data <- read.table(snme_data_list)$V1
bulk_methylation_data <- read.table(bulk_methylation_metadata_file, sep = ",", header = TRUE)$Sample.ID
scRNA_qc <- read.csv(scRNA_qc_file)
multiome_qc <- read.csv(multiome_qc_file)
original_unfiltered_metadata <- read.csv(original_unfiltered_metadata_file)

# Read in metadata and/or data for different assays

# Overall metadata
all_metadata <- read.table(paste0(overall_metadata_dir, "all_metadata_sheet.tsv"), sep = "\t", header = TRUE)
flu_ID_tokens <- read.table(paste0(overall_metadata_dir, "flu_data_tokens.tsv"), sep = "\t", header = TRUE)

# MintChIP
mintchip_metadata <- read.table(paste0(mintchip_metadata_dir, "mintchip_metadata.tsv"), sep = "\t", header = TRUE)

# MAGICAL
MAGICAL_hvl_placebo_results <- read.table(file = paste0(MAGICAL_hvl_placebo_output_dir, "MAGICAL_overall_output.tsv"), sep = "\t", header = TRUE)
# MAGICAL_hvl_vaccinated_results <- read.table(file = paste0(MAGICAL_hvl_vaccinated_output_dir, "MAGICAL_overall_output.tsv"), sep = "\t", header = TRUE)
# MAGICAL_lvl_placebo_results <- read.table(file = paste0(MAGICAL_lvl_output_dir, "MAGICAL_overall_output.tsv"), sep = "\t", header = TRUE)

# scRNA-seq - TODO: verify that SC and pseudobulk FC are in same direction
# HVL
scRNA_hvl_placebo_degs <- read.table(paste0(scRNA_hvl_placebo_deg_dir, "D28-vs-D_minus_1-degs-time_point.final.list.tsv"), sep = "\t", header = TRUE)
scRNA_hvl_placebo_degs <- scRNA_hvl_placebo_degs[scRNA_hvl_placebo_degs$Cell_Type != "Platelet",]
scRNA_hvl_placebo_MAGICAL_cell_types_degs <-  read.table(paste0(scRNA_hvl_placebo_MAGICAL_cell_types_deg_dir, "D28-vs-D_minus_1-degs-time_point.final.list.tsv"), sep = "\t", header = TRUE)
scRNA_hvl_placebo_MAGICAL_cell_types_degs <- scRNA_hvl_placebo_MAGICAL_cell_types_degs[scRNA_hvl_placebo_MAGICAL_cell_types_degs$Cell_Type != "Platelet",]
# TODO: Fix this!
scRNA_hvl_placebo_cell_metadata <- read.table(paste0(scRNA_hvl_dir, "HVL_RNA_cell_metadata.tsv"), sep = "\t", comment.char = "", header = TRUE)

scRNA_hvl_vaccinated_degs <- read.table(paste0(scRNA_hvl_vaccinated_deg_dir, "D28-vs-D_minus_1-degs-time_point.final.list.tsv"), sep = "\t", header = TRUE)
scRNA_hvl_vaccinated_degs <- scRNA_hvl_vaccinated_degs[scRNA_hvl_vaccinated_degs$Cell_Type != "Platelet",]
scRNA_hvl_vaccinated_MAGICAL_cell_types_degs <-  read.table(paste0(scRNA_hvl_vaccinated_MAGICAL_cell_types_deg_dir, "D28-vs-D_minus_1-degs-time_point.final.list.tsv"), sep = "\t", header = TRUE)
scRNA_hvl_vaccinated_MAGICAL_cell_types_degs <- scRNA_hvl_vaccinated_MAGICAL_cell_types_degs[scRNA_hvl_vaccinated_MAGICAL_cell_types_degs$Cell_Type != "Platelet",]

# LVL
scRNA_lvl_placebo_degs <- read.table(paste0(scRNA_lvl_placebo_deg_dir, "D28-vs-D_minus_1-degs-time_point.final.list.tsv"), sep = "\t", header = TRUE)
scRNA_lvl_placebo_degs <- scRNA_lvl_placebo_degs[scRNA_lvl_placebo_degs$Cell_Type != "Platelet",]
scRNA_lvl_placebo_MAGICAL_cell_types_degs <-  read.table(paste0(scRNA_lvl_placebo_MAGICAL_cell_types_deg_dir, "D28-vs-D_minus_1-degs-time_point.final.list.tsv"), sep = "\t", header = TRUE)
scRNA_lvl_placebo_MAGICAL_cell_types_degs <- scRNA_lvl_placebo_MAGICAL_cell_types_degs[scRNA_lvl_placebo_MAGICAL_cell_types_degs$Cell_Type != "Platelet",]





# Tables containing DEG results for single cell RNA-seq processing (LVL)
# Includes genes that passed pseudobulk filtering (remove platelets)
sc_pseudobulk_deg_lvl_table <- read.table(paste0(sc_deg_lvl_dir, "D28-vs-D_minus_1-degs-time_point.final.list.tsv"), sep = "\t", header = TRUE)
sc_pseudobulk_deg_lvl_table <- sc_pseudobulk_deg_lvl_table[sc_pseudobulk_deg_lvl_table$Cell_Type != "Platelet",]
sc_pseudobulk_deg_lvl_table <- subset(sc_pseudobulk_deg_lvl_table, sign(sc_log2FC) == sign(pseudo_bulk_log2FC))
sc_pseudobulk_deg_lvl_combined_cell_types_table <- read.table(paste0(sc_deg_lvl_combined_cell_types_dir, "D28-vs-D_minus_1-degs-time_point.final.list.tsv"), sep = "\t", header = TRUE)
sc_pseudobulk_deg_lvl_combined_cell_types_table <- sc_pseudobulk_deg_lvl_combined_cell_types_table[sc_pseudobulk_deg_lvl_combined_cell_types_table$Cell_Type != "Platelet",]
sc_pseudobulk_deg_lvl_combined_cell_types_table <- subset(sc_pseudobulk_deg_lvl_combined_cell_types_table, sign(sc_log2FC) == sign(pseudo_bulk_log2FC))
innate_sc_pseudobulk_deg_lvl_table <- sc_pseudobulk_deg_lvl_table[sc_pseudobulk_deg_lvl_table$Cell_Type %in% innate_cell_types,]
adaptive_sc_pseudobulk_deg_lvl_table <- sc_pseudobulk_deg_lvl_table[sc_pseudobulk_deg_lvl_table$Cell_Type %in% adaptive_cell_types,]

# Peak related tables
sc_peaks <- read.table(paste0(sc_das_dir, "HVL_peaks_info.txt"), sep = "\t", header = TRUE)
sc_motifs <- read.table(paste0(sc_das_dir, "peak_motif_matches.txt"), sep = "\t", header = TRUE)
atac_cell_metadata <- read.table(paste0(sc_das_dir, "HVL_ATAC_cell_metadata.tsv"), sep = "\t", comment.char = "")

# Create more lenient versions of sc_peaks (adds 250 or 500 to start / end of coordinates)
colnames(sc_peaks)[1] <- "seqnames"
sc_peaks_lenient <- sc_peaks
sc_peaks_lenient$start <- sc_peaks_lenient$start - 250 
sc_peaks_lenient$end <- sc_peaks_lenient$end + 250
sc_peaks_more_lenient <- sc_peaks
sc_peaks_more_lenient$start <- sc_peaks_lenient$start - 250 
sc_peaks_more_lenient$end <- sc_peaks_lenient$end + 250

sc_peaks_granges <- makeGRangesFromDataFrame(df = sc_peaks, keep.extra.columns = TRUE)
sc_peaks_lenient_granges <- makeGRangesFromDataFrame(df = sc_peaks_lenient, keep.extra.columns = TRUE)
sc_peaks_more_lenient_granges <- makeGRangesFromDataFrame(df = sc_peaks_more_lenient, keep.extra.columns = TRUE)

# snME related tables
snME_dms <- read.table(paste0(snME_data_dir, "snME_dms_processed.tsv"), sep = "\t", header = TRUE)
# Possible snME cell types
snME_dms_cell_types <- unique(snME_dms$celltype)

# miRNA tables
miRNA_raw_counts <- read.table(paste0(miRNA_data_dir, "miRNA_raw_counts.csv"), sep =",", header = TRUE)
miRNA_metadata_table <- read.table(paste0(miRNA_metadata_dir, "miRNA_metadata.tsv"), sep = "\t", header = TRUE)

# totalRNA tables
totalRNA_raw_counts <- read.table(paste0(totalRNA_data_dir, "STAREstCountFileT_Vaccitech_SampleNum_Corrected.csv"), sep =",", header = TRUE)
totalRNA_metadata_table <- read.table(paste0(totalRNA_metadata_dir, "totalRNA_metadata.tsv"), sep = "\t", header = TRUE)

# Grab gene lists from result tables and report number of genes
sc_pseudobulk_degs <- unique(sc_pseudobulk_deg_table$Gene_Name)
print(paste0("Number of DEGs that pass pseudobulk (scRNA): ", length(sc_pseudobulk_degs)))

pos_sc_pseudobulk_degs <- unique(sc_pseudobulk_deg_table[sc_pseudobulk_deg_table$sc_log2FC > 0,]$Gene_Name)
print(paste0("Number of positive DEGs that pass pseudobulk (scRNA): ", length(pos_sc_pseudobulk_degs)))

neg_sc_pseudobulk_degs <- unique(sc_pseudobulk_deg_table[sc_pseudobulk_deg_table$sc_log2FC < 0,]$Gene_Name)
print(paste0("Number of negative DEGs that pass pseudobulk (scRNA): ", length(neg_sc_pseudobulk_degs)))

innate_sc_pseudobulk_degs <- unique(innate_sc_pseudobulk_deg_table$Gene_Name)
print(paste0("Number of DEGs that pass pseudobulk in innate cell types (scRNA): ", length(innate_sc_pseudobulk_degs)))

pos_innate_sc_pseudobulk_degs <- unique(innate_sc_pseudobulk_deg_table[innate_sc_pseudobulk_deg_table$sc_log2FC > 0,]$Gene_Name)
print(paste0("Number of upregulated DEGs that pass pseudobulk in innate cell types (scRNA): ", length(pos_innate_sc_pseudobulk_degs)))

neg_innate_sc_pseudobulk_degs <- unique(innate_sc_pseudobulk_deg_table[innate_sc_pseudobulk_deg_table$sc_log2FC < 0,]$Gene_Name)
print(paste0("Number of downregulated DEGs that pass pseudobulk in innate cell types (scRNA): ", length(neg_innate_sc_pseudobulk_degs)))

adaptive_sc_pseudobulk_degs <- unique(adaptive_sc_pseudobulk_deg_table$Gene_Name)
print(paste0("Number of DEGs that pass pseudobulk in adaptive cell types (scRNA): ", length(adaptive_sc_pseudobulk_degs)))

pos_adaptive_sc_pseudobulk_degs <- unique(adaptive_sc_pseudobulk_deg_table[adaptive_sc_pseudobulk_deg_table$sc_log2FC > 0,]$Gene_Name)
print(paste0("Number of upregulated DEGs that pass pseudobulk in adaptive cell types (scRNA): ", length(pos_adaptive_sc_pseudobulk_degs)))

neg_adaptive_sc_pseudobulk_degs <- unique(adaptive_sc_pseudobulk_deg_table[adaptive_sc_pseudobulk_deg_table$sc_log2FC < 0,]$Gene_Name)
print(paste0("Number of downregulated DEGs that pass pseudobulk in adaptive cell types (scRNA): ", length(neg_adaptive_sc_pseudobulk_degs)))

# Next, we can select bulk RNA-seq associated with the specific subjects that we used for single-cell / multiome
sc_aliquots <- c("91910a04711aa3dd","3731a6247ae23831","2232300b0e3a3d06","76ea83ff9293871a","5fdfdbaeb3c8eee8","981520e7e138d460","bb3d7b309cb5fc58","8338411dc3e181e9","da4fe21a89c8f7f4","41d248a6ec3b87e2","e3e01c75894ef461","4534496c580cb408") # 12 samples - 6 paired
sc_subjects <- as.character(unique(all_metadata[all_metadata$aliquot_id %in% sc_aliquots,]$subject_id))
multiome_aliquots <- c("717579a2ae2fb6c2", "dde63f8ca98af665", "a464019298ae6682", "d554be0e36e4d789", "e43db0f72b9c2e31", "b82bb7c75d47dac1", "6f609a68dca1261f", "9c6ec1b704700c7d", "7b54cfac7e67b0fa", "575d74707585856a", "c1eb160d7bd1f29f", "8832fff8247b18b9", "abf6d19ee03be1e8", "216bb226181591dd") # 14 samples - 7 paired
multiome_subjects <- as.character(unique(all_metadata[all_metadata$aliquot_id %in% multiome_aliquots,]$subject_id))

sc_placebo_metadata <- placebo_metadata[(placebo_metadata$subject_id %in% sc_subjects),]
sc_placebo_counts <- placebo_counts[,(colnames(placebo_counts) %in% rownames(sc_placebo_metadata))]

multiome_placebo_metadata <- placebo_metadata[(placebo_metadata$subject_id %in% multiome_subjects),]
multiome_placebo_counts <- placebo_counts[,(colnames(placebo_counts) %in% rownames(multiome_placebo_metadata))]

high_sc_placebo_metadata <- high_placebo_metadata[(high_placebo_metadata$subject_id %in% sc_subjects),]
high_sc_placebo_counts <- high_placebo_counts[,(colnames(high_placebo_counts) %in% rownames(high_sc_placebo_metadata))]

high_non_sc_placebo_metadata <- high_placebo_metadata[!(high_placebo_metadata$subject_id %in% sc_subjects),]
high_non_sc_placebo_counts <- high_placebo_counts[,(colnames(high_placebo_counts) %in% rownames(high_non_sc_placebo_metadata))]

high_multiome_placebo_metadata <- high_placebo_metadata[(high_placebo_metadata$subject_id %in% multiome_subjects),]
high_multiome_placebo_counts <- high_placebo_counts[,(colnames(high_placebo_counts) %in% rownames(high_multiome_placebo_metadata))]

high_non_multiome_placebo_metadata <- high_placebo_metadata[!(high_placebo_metadata$subject_id %in% multiome_subjects),]
high_non_multiome_placebo_counts <- high_placebo_counts[,(colnames(high_placebo_counts) %in% rownames(high_non_multiome_placebo_metadata))]

low_sc_placebo_metadata <- low_placebo_metadata[(low_placebo_metadata$subject_id %in% sc_subjects),]
low_sc_placebo_counts <- low_placebo_counts[,(colnames(low_placebo_counts) %in% rownames(low_sc_placebo_metadata))]

low_non_sc_placebo_metadata <- low_placebo_metadata[!(low_placebo_metadata$subject_id %in% sc_subjects),]
low_non_sc_placebo_counts <- low_placebo_counts[,(colnames(low_placebo_counts) %in% rownames(low_non_sc_placebo_metadata))]

low_multiome_placebo_metadata <- low_placebo_metadata[(low_placebo_metadata$subject_id %in% multiome_subjects),]
low_multiome_placebo_counts <- low_placebo_counts[,(colnames(low_placebo_counts) %in% rownames(low_multiome_placebo_metadata))]

low_non_multiome_placebo_metadata <- low_placebo_metadata[!(low_placebo_metadata$subject_id %in% multiome_subjects),]
low_non_multiome_placebo_counts <- low_placebo_counts[,(colnames(low_placebo_counts) %in% rownames(low_non_multiome_placebo_metadata))]