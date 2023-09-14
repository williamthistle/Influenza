library(data.table)
library(DESeq2)
library(MetaIntegrator)
library(openxlsx)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(DiffBind)
library(SPEEDI)

set.seed(get_speedi_seed())

base_dir <- "~/GitHub/Influenza/"
source(paste0(base_dir, "bulk_RNA_analysis_helper.R"))
source(paste0(base_dir, "pseudobulk_analysis_helper.R"))
source(paste0(base_dir, "Data Compendium/Compendium_Functions.R"))
setup_bulk_analysis()
sample_metadata <- read.table(paste0(base_dir, "all_metadata_sheet.tsv"), sep = "\t", header = TRUE)
possible_cell_types <- c("CD4_Naive", "CD8_Naive", "CD4_Memory", "CD8_Memory", "cDC", "HSPC", "pDC", "Platelet", "Plasmablast", "Proliferating", "NK", "T_Naive", "CD14_Mono", "CD16_Mono", "MAIT")
onedrive_dir <- "~"
setwd(onedrive_dir)
setwd("../")
onedrive_dir <- getwd()
onedrive_dir <- paste0(onedrive_dir, "/OneDrive - Princeton University/")
sc_pseudobulk_dir <- paste0(onedrive_dir, "Influenza Analysis/Vaccitech Paper Analysis/Current Analyses/Single Cell/RNA/HVL/DEGs/")
sc_humanbase_dir <- paste0(onedrive_dir, "Influenza Analysis/Vaccitech Paper Analysis/Current Analyses/Single Cell/RNA/HVL/HumanBase/")
multiome_pseudobulk_dir <- paste0(onedrive_dir, "Influenza Analysis/Vaccitech Paper Analysis/Current Analyses/Multiome/RNA/HVL/DEGs/")
mintchip_dir <- paste0(onedrive_dir, "Influenza Analysis/MintChIP/")
sc_peak_dir <- paste0(onedrive_dir, "Influenza Analysis/Vaccitech Paper Analysis/Current Analyses/Single Cell/ATAC/HVL/")

# Read in Mint-ChIP metadata
mintchip_metadata <- read.table(paste0(mintchip_dir, "mintchip_metadata.tsv"), sep = "\t", header = TRUE)
# Tables containing results for single cell and multiome RNA-seq processing
# Includes genes that passed pseudobulk filtering and genes that passed MAGICAL filtering (remove platelets)
sc_pseudobulk_gene_table <- read.table(paste0(sc_pseudobulk_dir, "D28-vs-D_minus_1-degs-time_point.final.list.tsv"), sep = "\t", header = TRUE)
sc_pseudobulk_gene_table <- sc_pseudobulk_gene_table[sc_pseudobulk_gene_table$Cell_Type != "Platelet",]
multiome_pseudobulk_gene_table <- read.table(paste0(multiome_pseudobulk_dir, "D28-vs-D_minus_1-degs-time_point.final.list.tsv"), sep = "\t", header = TRUE)
multiome_pseudobulk_gene_table <- multiome_pseudobulk_gene_table[multiome_pseudobulk_gene_table$Cell_Type != "Platelet",]

# Peak related tables
mintchip_table <- read.xlsx(xlsxFile = paste0(mintchip_dir, "MultiMarkerBestSignatureAnnotations_All3000.xlsx"), sheet = 1)
sc_peaks <- read.table(paste0(sc_peak_dir, "HVL_peaks_info.txt"), sep = "\t", header = TRUE)
sc_motifs <- read.table(paste0(sc_peak_dir, "peak_motif_matches.txt"), sep = "\t", header = TRUE)
sc_motifs$chr <- paste0("chr", sc_motifs$chr)
sc_das_lenient <- read.table(paste0(sc_peak_dir, "diff_peaks/D28_D1_diff_lenient_final.tsv"), sep = "\t", header = TRUE)
sc_motif_lenient_subset <- sc_motifs[sc_motifs$chr %in% sc_das_lenient$chr & sc_motifs$point1 %in% sc_das_lenient$start & sc_motifs$point2 %in% sc_das_lenient$end,]
sc_motif_lenient_subset$chr <- as.numeric(substr(sc_motif_lenient_subset$chr, 4, 6))
sc_motifs$chr <- as.numeric(substr(sc_motifs$chr, 4, 6))
sc_das_stricter <- read.table(paste0(sc_peak_dir, "diff_peaks/D28_D1_diff_stricter_final.tsv"), sep = "\t", header = TRUE)
sc_das_strictest <- read.table(paste0(sc_peak_dir, "diff_peaks/D28_D1_diff_strictest_final.tsv"), sep = "\t", header = TRUE)

# Grab gene lists from result tables and report number of genes
sc_pseudobulk_genes <- unique(sc_pseudobulk_gene_table$Gene_Name)
print(paste0("Number of genes that pass pseudobulk (scRNA): ", length(sc_pseudobulk_genes)))
multiome_pseudobulk_genes <- unique(multiome_pseudobulk_gene_table$Gene_Name)
print(paste0("Number of genes that pass pseudobulk (multiome): ", length(multiome_pseudobulk_genes)))

# Next, let's test our gene lists on the actual bulk RNA-seq data! (Should I add extra samples for subjects that have them? How to label as HVL or LVL?)
# First, let's remove the 0 PCR sample from low because it's questionable
removed_low_viral_aliquots <- rownames(placebo_metadata[placebo_metadata$subject_id == "f18c54d93cef4a4e",])
placebo_metadata <- placebo_metadata[!(placebo_metadata$subject_id %in% "f18c54d93cef4a4e"),]
placebo_counts <- placebo_counts[,!(colnames(placebo_counts) %in% removed_low_viral_aliquots)]
low_placebo_metadata <- low_placebo_metadata[!(low_placebo_metadata$subject_id %in% "f18c54d93cef4a4e"),]
low_placebo_counts <- low_placebo_counts[,!(colnames(low_placebo_counts) %in% removed_low_viral_aliquots)]
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