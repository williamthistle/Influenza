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

# Load bulk RNA-seq analysis results
onedrive_dir <- "~"
setwd(onedrive_dir)
setwd("../")
onedrive_dir <- getwd()
onedrive_dir <- paste0(onedrive_dir, "/OneDrive - Princeton University/")
load(paste0(onedrive_dir, "Influenza Analysis/bulk_RNA_analysis.RData"))

set.seed(get_speedi_seed())
options(max.print=10000)

base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
data_dir <- "~/local_data_files/"
source(paste0(base_dir, "extra_functions/bulk_RNA_analysis_helper.R"))
source(paste0(base_dir, "extra_functions/Compendium_Functions.R"))
source(paste0(base_dir, "extra_functions/humanbase_functions.R"))
source(paste0(base_dir, "extra_functions/pseudobulk_analysis_helper.R"))
source(paste0(base_dir, "extra_functions/MAGICAL_functions.R"))
setup_bulk_analysis(base_dir = base_dir, data_dir = data_dir)
sample_metadata <- read.table(paste0(base_dir, "metadata/all_metadata_sheet.tsv"), sep = "\t", header = TRUE)
flu_tokens <- read.table(paste0(base_dir, "metadata/flu_data_tokens.tsv"), sep = "\t", header = TRUE)
possible_cell_types <- c("CD4_Naive", "CD8_Naive", "CD4_Memory", "CD8_Memory", "cDC", "HSPC", "pDC", "Platelet", "Plasmablast", "Proliferating", "NK", "NK_CD56bright", "T_Naive", "CD14_Mono", "CD16_Mono", "MAIT", "B", "B_naive", "B_memory")
innate_cell_types <- c("CD16_Mono","CD14_Mono","cDC","pDC","NK","NK_CD56bright")
# Reload one drive dir
onedrive_dir <- "~"
setwd(onedrive_dir)
setwd("../")
onedrive_dir <- getwd()
onedrive_dir <- paste0(onedrive_dir, "/OneDrive - Princeton University/")
# Set other dirs
sc_rna_dir <- paste0(onedrive_dir, "Influenza Analysis/Vaccitech Paper Analysis/Current Analyses/Single Cell/RNA/HVL/")
sc_deg_dir <- paste0(onedrive_dir, "Influenza Analysis/Vaccitech Paper Analysis/Current Analyses/Single Cell/RNA/HVL/DEGs/")
sc_deg_combined_cell_types_dir <- paste0(onedrive_dir, "Influenza Analysis/Vaccitech Paper Analysis/Current Analyses/Single Cell/RNA/HVL/DEGs_with_ATAC_cell_types/")
sc_magical_dir <- paste0(onedrive_dir, "Influenza Analysis/Vaccitech Paper Analysis/Current Analyses/Single Cell/MAGICAL/HVL/")
sc_humanbase_dir <- paste0(onedrive_dir, "Influenza Analysis/Vaccitech Paper Analysis/Current Analyses/Single Cell/RNA/HVL/HumanBase/")
sc_das_dir <- paste0(onedrive_dir, "Influenza Analysis/Vaccitech Paper Analysis/Current Analyses/Single Cell/ATAC/HVL/")
mintchip_metadata_dir <- paste0(onedrive_dir, "Influenza Analysis/MintChIP/")
mintchip_das_dir <- paste0(onedrive_dir, "Influenza Analysis/Vaccitech Paper Analysis/Current Analyses/MintChIP/differentially_accessible_sites/")

# Read in Mint-ChIP metadata
mintchip_metadata <- read.table(paste0(mintchip_metadata_dir, "mintchip_metadata.tsv"), sep = "\t", header = TRUE)
# Tables containing DEG results for single cell RNA-seq processing
# Includes genes that passed pseudobulk filtering (remove platelets)
# TODO: Remove genes that have different sign for sc FC and pseudobulk FC? Not relevant for my current SC dataset at least
sc_pseudobulk_deg_table <- read.table(paste0(sc_deg_dir, "D28-vs-D_minus_1-degs-time_point.final.list.tsv"), sep = "\t", header = TRUE)
sc_pseudobulk_deg_table <- sc_pseudobulk_deg_table[sc_pseudobulk_deg_table$Cell_Type != "Platelet",]
sc_pseudobulk_deg_combined_cell_types_table <- read.table(paste0(sc_deg_combined_cell_types_dir, "D28-vs-D_minus_1-degs-time_point.final.list.tsv"), sep = "\t", header = TRUE)
sc_pseudobulk_deg_combined_cell_types_table <- sc_pseudobulk_deg_combined_cell_types_table[sc_pseudobulk_deg_combined_cell_types_table$Cell_Type != "Platelet",]
sc_pseudobulk_deg_combined_cell_types_table$Cell_Type[sc_pseudobulk_deg_combined_cell_types_table$Cell_Type == "NK_MAGICAL"] <- "NK"
innate_sc_pseudobulk_deg_table <- sc_pseudobulk_deg_table[sc_pseudobulk_deg_table$Cell_Type %in% innate_cell_types,]
rna_cell_metadata <- read.table(paste0(sc_rna_dir, "HVL_RNA_cell_metadata.tsv"), sep = "\t", comment.char = "", header = TRUE)

# Peak related tables
sc_peaks <- read.table(paste0(sc_das_dir, "HVL_peaks_info.txt"), sep = "\t", header = TRUE)
sc_motifs <- read.table(paste0(sc_das_dir, "peak_motif_matches.txt"), sep = "\t", header = TRUE)
sc_motifs$chr <- paste0("chr", sc_motifs$chr)
sc_das_lenient <- read.table(paste0(sc_das_dir, "diff_peaks/D28_D1_diff_lenient_final.tsv"), sep = "\t", header = TRUE)
sc_motif_lenient_subset <- sc_motifs[sc_motifs$chr %in% sc_das_lenient$chr & sc_motifs$point1 %in% sc_das_lenient$start & sc_motifs$point2 %in% sc_das_lenient$end,]
sc_motif_lenient_subset$chr <- as.numeric(substr(sc_motif_lenient_subset$chr, 4, 6))
sc_motifs$chr <- as.numeric(substr(sc_motifs$chr, 4, 6))
sc_das_stricter <- read.table(paste0(sc_das_dir, "diff_peaks/D28_D1_diff_stricter_final.tsv"), sep = "\t", header = TRUE)
sc_das_strictest <- read.table(paste0(sc_das_dir, "diff_peaks/D28_D1_diff_strictest_final.tsv"), sep = "\t", header = TRUE)
atac_cell_metadata <- read.table(paste0(sc_das_dir, "HVL_cell_metadata.tsv"), sep = "\t", comment.char = "", header = TRUE)

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