library(data.table)
library(openxlsx)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)

# Set up 
base_dir <- "~/GitHub/Influenza/"
single_cell_magical_dir <- "C:/Users/willi/OneDrive - Princeton University/Influenza Analysis/Single Cell RNA-Seq/MAGICAL Analyses/Placebo 4 Subject HVL (SPEEDI) - SCT/"
multiome_magical_dir <- "C:/Users/willi/OneDrive - Princeton University/Influenza Analysis/True Multiome/MAGICAL Analyses/14 Placebo Sample (Final)/"
old_wd <- getwd()
setwd("~/")
setwd("../OneDrive - Princeton University/Influenza Analysis/Data Compendium")
mintchip_gene_table <- read.xlsx("../MintChIP/MultiMarkerBestSignatureAnnotations_All3000.xlsx", sheet = 1)
setwd(old_wd)
set.seed(2000)

# Tables containing results for single cell and multiome RNA-seq processing
# Includes genes that passed pseudobulk filtering and genes that passed MAGICAL filtering
# Note that LR includes latent variable subject in differential expression analysis
single_cell_pseudobulk_pass_peak_table <- read.table(paste0(single_cell_magical_dir, "D28_D1_MAGICAL_peaks_passing_pseudobulk_unadjusted.txt"), sep = "\t", header = TRUE)
#single_cell_pseudobulk_pass_peak_FDR_table <- read.table(paste0(single_cell_magical_dir, "D28_D1_MAGICAL_peaks_passing_pseudobulk_FDR.txt"), sep = "\t", header = TRUE)
cell_types <- unique(single_cell_pseudobulk_pass_peak_table$cell_type)
#cell_types_FDR <- unique(single_cell_pseudobulk_pass_peak_FDR_table$cell_type)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

cell_type_peak_annotations <- list()
cell_type_peak_annotations_FDR <- list()
index <- 1

for(cell_type in cell_types) {
  cell_type_pseudobulk_pass_peak_table <- single_cell_pseudobulk_pass_peak_table[single_cell_pseudobulk_pass_peak_table$cell_type == cell_type,]
  chipseeker_input_peaks <- cell_type_pseudobulk_pass_peak_table[,c(2,3,4)]
  
  chipseeker_input_peaks$chr[chipseeker_input_peaks$chr == "23"] <- "X"
  chipseeker_input_peaks$chr <- paste0('chr', chipseeker_input_peaks$chr)
  chipseeker_input_peaks$length <- chipseeker_input_peaks$end - chipseeker_input_peaks$start
  
  chipseeker_output <- annotatePeak(makeGRangesFromDataFrame(chipseeker_input_peaks), TxDb = txdb, annoDb = "org.Hs.eg.db")
  chipseeker_output_df <- as.data.frame(chipseeker_output)
  cell_type_peak_annotations[[index]] <- chipseeker_output_df
  index <- index + 1
}

index <- 1

for(cell_type in cell_types_FDR) {
  cell_type_pseudobulk_pass_peak_table_FDR <- single_cell_pseudobulk_pass_peak_FDR_table[single_cell_pseudobulk_pass_peak_FDR_table$cell_type == cell_type,]
  chipseeker_input_peaks_FDR <- cell_type_pseudobulk_pass_peak_table_FDR[,c(2,3,4)]
  
  chipseeker_input_peaks_FDR$chr[chipseeker_input_peaks_FDR$chr == "23"] <- "X"
  chipseeker_input_peaks_FDR$chr <- paste0('chr', chipseeker_input_peaks_FDR$chr)
  chipseeker_input_peaks_FDR$length <- chipseeker_input_peaks_FDR$end - chipseeker_input_peaks_FDR$start
  
  chipseeker_output_FDR <- annotatePeak(makeGRangesFromDataFrame(chipseeker_input_peaks_FDR), TxDb = txdb, annoDb = "org.Hs.eg.db")
  chipseeker_output_df_FDR <- as.data.frame(chipseeker_output_FDR)
  cell_type_peak_annotations_FDR[[index]] <- chipseeker_output_df_FDR
  index <- index + 1
}

# Figure out why NA is there
mint_cell_type_peak_annotations <- list()
index <- 1
total_curated_peaks <- 0
total_genes <- c()
for(cell_type in cell_types) {
  mint_cell_type_peak_annotations[[index]] <- cell_type_peak_annotations[[index]][cell_type_peak_annotations[[index]]$SYMBOL %in% mintchip_gene_table$SYMBOL,]
  print(cell_type)
  print(nrow(mint_cell_type_peak_annotations[[index]]))
  total_curated_peaks <- total_curated_peaks + nrow(mint_cell_type_peak_annotations[[index]])
  total_genes <- c(total_genes, mint_cell_type_peak_annotations[[index]]$SYMBOL)
  index <- index + 1
}

index <- 1
for(mint_cell_type_peak_annotation in mint_cell_type_peak_annotations) {
  current_cell_type <- cell_types[index]
  for(row_index in 1:nrow(mint_cell_type_peak_annotation)) {
    seqname <- mint_cell_type_peak_annotation[row_index,]$seqnames
    start <- mint_cell_type_peak_annotation[row_index,]$start
    end <- mint_cell_type_peak_annotation[row_index,]$end
    current_peak_table <- read.xlsx(paste0(single_cell_magical_dir, "scATAC_DAPs/", current_cell_type, "_D28_D1_diff.xlsx"), sheet = 1)
  }
  index <- index + 1
}












print(length(cell_type_peak_annotations[[1]]$SYMBOL))
print(length(intersect(cell_type_peak_annotations[[1]]$SYMBOL, mintchip_gene_table$SYMBOL)))

mint_cell_type_peak_annotations <- cell_type_peak_annotations[[1]][cell_type_peak_annotations[[1]]$SYMBOL %in% mintchip_gene_table$SYMBOL,]
mint_cell_type_genes <- unique(mint_cell_type_peak_annotations$SYMBOL)
cell_type_peak_annotations_subset <- cell_type_peak_annotations[[1]]
