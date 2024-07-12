# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

search_terms <- c("histone", "interferon", "interleukin", "AP-1", "methyltransferase", "acetyltransferase", "demethylase")
gene_terms <- c("CCL3", "CX3CR1", "CCL3L1", "CXCL16", "IL32", "CASP1", "NFIL3", "IRAK3", "IL1RAP", "RIPK1", "PTGES", "CEBPB", 
                "IRF2", "IRF7", "IFNG", "OAS1", "NMRAL1", "MNDA", "GBP1", "PSMB9", "IFNGR1", "DNAJC3", "GBP5", "CEMIP2", "USP38",  
                "JUN", "JUNB", "FOSB", "FOSL2", "JDP2", 
                "MAP3K11", "CSK", "DUSP1", "DUSP2", "DUSP6", "TRAF3IP3", "ABHD17A", "CSF1R", "MAPK7", "MAPK8", "MAP2K1", "MAP3K8", "MAP3K20", "MAP4K3", "MAPKAPK2",
                "PTK2B", "RELL1", "MINK1", "BRAF",
                "NFKB1", "JAK1", "STAT3", "STAT4",
                "IRF", "STAT", "JAK", "MAP", "DUSP")
gene_terms <- c(gene_terms, c("HIST1H1C", "HIST1H1D", "HIST1H1E", "H2AFZ", "KAT6A", "BRD1", "MSL2", "HDAC5", "HDAC7", "HDAC9", "METTL23", "ASH1L", "PRDM2", "SETD2","KDM2B", "KDM3A",
                              "KAT", "HDAC", "KDM"))

txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene
txdb_hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene
promoter_terms <- c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)")

# MintChIP - "H3K4me1", "H3K4me3", "H3K9me3", "H3K27Ac", "H3K27me3", "H3K36me3"
# H2K27ac
mintchip_H3K27ac_peaks <- read.table(paste0(mintchip_das_dir, "H3K27Ac/H3K27Ac_DESeq2_FC_0.1.tsv"), sep = "\t", header = TRUE)
mintchip_H3K27ac_peaks_annotated <- annotatePeak(makeGRangesFromDataFrame(mintchip_H3K27ac_peaks), TxDb = txdb_hg19, annoDb = "org.Hs.eg.db")
mintchip_H3K27ac_peaks_annotated <- as.data.frame(mintchip_H3K27ac_peaks_annotated)
mintchip_H3K27ac_peaks_annotated$Fold <- mintchip_H3K27ac_peaks$Fold
mintchip_H3K27ac_peaks_annotated$p.value <- mintchip_H3K27ac_peaks$p.value

mintchip_H3K27ac_peaks_annotated_final <- mintchip_H3K27ac_peaks_annotated %>%
  filter(
    str_detect(SYMBOL, paste(gene_terms, collapse = "|")) |
      str_detect(GENENAME, paste(search_terms, collapse = "|"))
  )

mintchip_H3K27ac_peaks_annotated_final_pos <- mintchip_H3K27ac_peaks_annotated_final[mintchip_H3K27ac_peaks_annotated_final$Fold > 0,]
mintchip_H3K27ac_peaks_annotated_final_neg <- mintchip_H3K27ac_peaks_annotated_final[mintchip_H3K27ac_peaks_annotated_final$Fold < 0,]

# H3K4me1
mintchip_H3K4me1_peaks <- read.table(paste0(mintchip_das_dir, "H3K4me1/H3K4me1_DESeq2_FC_0.3.tsv"), sep = "\t", header = TRUE)
mintchip_H3K4me1_peaks_annotated <- annotatePeak(makeGRangesFromDataFrame(mintchip_H3K4me1_peaks), TxDb = txdb_hg19, annoDb = "org.Hs.eg.db")
mintchip_H3K4me1_peaks_annotated <- as.data.frame(mintchip_H3K4me1_peaks_annotated)
mintchip_H3K4me1_peaks_annotated$Fold <- mintchip_H3K4me1_peaks$Fold
mintchip_H3K4me1_peaks_annotated$p.value <- mintchip_H3K4me1_peaks$p.value

mintchip_H3K4me1_peaks_annotated_final <- mintchip_H3K4me1_peaks_annotated %>%
  filter(
    str_detect(SYMBOL, paste(gene_terms, collapse = "|")) |
      str_detect(GENENAME, paste(search_terms, collapse = "|"))
  )

mintchip_H3K4me1_peaks_annotated_final_pos <- mintchip_H3K4me1_peaks_annotated_final[mintchip_H3K4me1_peaks_annotated_final$Fold > 0,]
mintchip_H3K4me1_peaks_annotated_final_neg <- mintchip_H3K4me1_peaks_annotated_final[mintchip_H3K4me1_peaks_annotated_final$Fold < 0,]

# H3K4me3
mintchip_H3K4me3_peaks <- read.table(paste0(mintchip_das_dir, "H3K4me3/H3K4me3_DESeq2_FC_0.1.tsv"), sep = "\t", header = TRUE)
mintchip_H3K4me3_peaks_annotated <- annotatePeak(makeGRangesFromDataFrame(mintchip_H3K4me3_peaks), TxDb = txdb_hg19, annoDb = "org.Hs.eg.db")
mintchip_H3K4me3_peaks_annotated <- as.data.frame(mintchip_H3K4me3_peaks_annotated)
mintchip_H3K4me3_peaks_annotated$Fold <- mintchip_H3K4me3_peaks$Fold
mintchip_H3K4me3_peaks_annotated$p.value <- mintchip_H3K4me3_peaks$p.value

mintchip_H3K4me3_peaks_annotated_final <- mintchip_H3K4me3_peaks_annotated %>%
  filter(
    str_detect(SYMBOL, paste(gene_terms, collapse = "|")) |
      str_detect(GENENAME, paste(search_terms, collapse = "|"))
  )

mintchip_H3K4me3_peaks_annotated_final_pos <- mintchip_H3K4me3_peaks_annotated_final[mintchip_H3K4me3_peaks_annotated_final$Fold > 0,]
mintchip_H3K4me3_peaks_annotated_final_neg <- mintchip_H3K4me3_peaks_annotated_final[mintchip_H3K4me3_peaks_annotated_final$Fold < 0,]

# H3K9me3
mintchip_H3K9me3_peaks <- read.table(paste0(mintchip_das_dir, "H3K9me3/H3K9me3_DESeq2_FC_0.3.tsv"), sep = "\t", header = TRUE)
mintchip_H3K9me3_peaks_annotated <- annotatePeak(makeGRangesFromDataFrame(mintchip_H3K9me3_peaks), TxDb = txdb_hg19, annoDb = "org.Hs.eg.db")
mintchip_H3K9me3_peaks_annotated <- as.data.frame(mintchip_H3K9me3_peaks_annotated)
mintchip_H3K9me3_peaks_annotated$Fold <- mintchip_H3K9me3_peaks$Fold
mintchip_H3K9me3_peaks_annotated$p.value <- mintchip_H3K9me3_peaks$p.value

mintchip_H3K9me3_peaks_annotated_final <- mintchip_H3K9me3_peaks_annotated %>%
  filter(
    str_detect(SYMBOL, paste(gene_terms, collapse = "|")) |
      str_detect(GENENAME, paste(search_terms, collapse = "|"))
  )

mintchip_H3K9me3_peaks_annotated_final_pos <- mintchip_H3K9me3_peaks_annotated_final[mintchip_H3K9me3_peaks_annotated_final$Fold > 0,]
mintchip_H3K9me3_peaks_annotated_final_neg <- mintchip_H3K9me3_peaks_annotated_final[mintchip_H3K9me3_peaks_annotated_final$Fold < 0,]

# H3K27me3
mintchip_H3K27me3_peaks <- read.table(paste0(mintchip_das_dir, "H3K27me3/H3K27me3_DESeq2_FC_0.1.tsv"), sep = "\t", header = TRUE)
mintchip_H3K27me3_peaks_annotated <- annotatePeak(makeGRangesFromDataFrame(mintchip_H3K27me3_peaks), TxDb = txdb_hg19, annoDb = "org.Hs.eg.db")
mintchip_H3K27me3_peaks_annotated <- as.data.frame(mintchip_H3K27me3_peaks_annotated)
mintchip_H3K27me3_peaks_annotated$Fold <- mintchip_H3K27me3_peaks$Fold
mintchip_H3K27me3_peaks_annotated$p.value <- mintchip_H3K27me3_peaks$p.value

mintchip_H3K27me3_peaks_annotated_final <- mintchip_H3K27me3_peaks_annotated %>%
  filter(
    str_detect(SYMBOL, paste(gene_terms, collapse = "|")) |
      str_detect(GENENAME, paste(search_terms, collapse = "|"))
  )

mintchip_H3K27me3_peaks_annotated_final_pos <- mintchip_H3K27me3_peaks_annotated_final[mintchip_H3K27me3_peaks_annotated_final$Fold > 0,]
mintchip_H3K27me3_peaks_annotated_final_neg <- mintchip_H3K27me3_peaks_annotated_final[mintchip_H3K27me3_peaks_annotated_final$Fold < 0,]

# H3K36me3
mintchip_H3K36me3_peaks <- read.table(paste0(mintchip_das_dir, "H3K36me3/H3K36me3_DESeq2_FC_0.3.tsv"), sep = "\t", header = TRUE)
mintchip_H3K36me3_peaks_annotated <- annotatePeak(makeGRangesFromDataFrame(mintchip_H3K36me3_peaks), TxDb = txdb_hg19, annoDb = "org.Hs.eg.db")
mintchip_H3K36me3_peaks_annotated <- as.data.frame(mintchip_H3K36me3_peaks_annotated)
mintchip_H3K36me3_peaks_annotated$Fold <- mintchip_H3K36me3_peaks$Fold
mintchip_H3K36me3_peaks_annotated$p.value <- mintchip_H3K36me3_peaks$p.value

mintchip_H3K36me3_peaks_annotated_final <- mintchip_H3K36me3_peaks_annotated %>%
  filter(
    str_detect(SYMBOL, paste(gene_terms, collapse = "|")) |
      str_detect(GENENAME, paste(search_terms, collapse = "|"))
  )

mintchip_H3K36me3_peaks_annotated_final_pos <- mintchip_H3K36me3_peaks_annotated_final[mintchip_H3K36me3_peaks_annotated_final$Fold > 0,]
mintchip_H3K36me3_peaks_annotated_final_neg <- mintchip_H3K36me3_peaks_annotated_final[mintchip_H3K36me3_peaks_annotated_final$Fold < 0,]

# snME

# Monocytes
snME_monocytes_peaks <- snME_dms[snME_dms$celltype == "Monocyte",]
snME_monocytes_peaks_annotated <- annotatePeak(makeGRangesFromDataFrame(snME_monocytes_peaks), TxDb = txdb_hg38, annoDb = "org.Hs.eg.db")
snME_monocytes_peaks_annotated <- as.data.frame(snME_monocytes_peaks_annotated)
snME_monocytes_peaks <- snME_monocytes_peaks[snME_monocytes_peaks$chr %in% snME_monocytes_peaks_annotated$seqnames &
                                               snME_monocytes_peaks$start %in% snME_monocytes_peaks_annotated$start & 
                                               snME_monocytes_peaks$end %in% snME_monocytes_peaks_annotated$end,]
snME_monocytes_peaks_annotated$methylation <- snME_monocytes_peaks$methylation
snME_monocytes_peaks_annotated$D1 <- snME_monocytes_peaks$D1
snME_monocytes_peaks_annotated$D28 <- snME_monocytes_peaks$D28

snME_monocytes_peaks_annotated_final <- snME_monocytes_peaks_annotated %>%
  filter(
    str_detect(SYMBOL, paste(gene_terms, collapse = "|")) |
      str_detect(GENENAME, paste(search_terms, collapse = "|"))
  )

snME_monocytes_peaks_annotated_final_pos <- snME_monocytes_peaks_annotated_final[snME_monocytes_peaks_annotated_final$methylation == "upregulated",]
snME_monocytes_peaks_annotated_final_pos_promoters <- snME_monocytes_peaks_annotated_final_pos[snME_monocytes_peaks_annotated_final_pos$annotation %in% promoter_terms,]
snME_monocytes_peaks_annotated_final_neg <- snME_monocytes_peaks_annotated_final[snME_monocytes_peaks_annotated_final$methylation == "downregulated",]
snME_monocytes_peaks_annotated_final_neg_promoters <- snME_monocytes_peaks_annotated_final_neg[snME_monocytes_peaks_annotated_final_neg$annotation %in% promoter_terms,]

# NK
snME_NK_peaks <- snME_dms[snME_dms$celltype == "NK-cell2",]
snME_NK_peaks_annotated <- annotatePeak(makeGRangesFromDataFrame(snME_NK_peaks), TxDb = txdb_hg38, annoDb = "org.Hs.eg.db")
snME_NK_peaks_annotated <- as.data.frame(snME_NK_peaks_annotated)
snME_NK_peaks <- snME_NK_peaks[snME_NK_peaks$chr %in% snME_NK_peaks_annotated$seqnames &
                                               snME_NK_peaks$start %in% snME_NK_peaks_annotated$start & 
                                               snME_NK_peaks$end %in% snME_NK_peaks_annotated$end,]
snME_NK_peaks_annotated$methylation <- snME_NK_peaks$methylation
snME_NK_peaks_annotated$D1 <- snME_NK_peaks$D1
snME_NK_peaks_annotated$D28 <- snME_NK_peaks$D28

snME_NK_peaks_annotated_final <- snME_NK_peaks_annotated %>%
  filter(
    str_detect(SYMBOL, paste(gene_terms, collapse = "|")) |
      str_detect(GENENAME, paste(search_terms, collapse = "|"))
  )

snME_NK_peaks_annotated_final_pos <- snME_NK_peaks_annotated_final[snME_NK_peaks_annotated_final$methylation == "upregulated",]
snME_NK_peaks_annotated_final_pos_promoters <- snME_NK_peaks_annotated_final_pos[snME_NK_peaks_annotated_final_pos$annotation %in% promoter_terms,]
snME_NK_peaks_annotated_final_neg <- snME_NK_peaks_annotated_final[snME_NK_peaks_annotated_final$methylation == "downregulated",]
snME_NK_peaks_annotated_final_neg_promoters <- snME_NK_peaks_annotated_final_neg[snME_NK_peaks_annotated_final_neg$annotation %in% promoter_terms,]

# CD14 mono peaks
scATAC_CD14_mono_peaks <- read.table(paste0(scATAC_hvl_placebo_das_dir, "D28-vs-D_minus_1-degs-CD14_Mono-time_point-controlling_for_subject_id_final_pct_0.01.tsv"), sep = "\t", header = TRUE)

chromosomes <- sapply(strsplit(scATAC_CD14_mono_peaks$Peak_Name, "-"), `[`, 1)
start_coords <- as.numeric(sapply(strsplit(scATAC_CD14_mono_peaks$Peak_Name, "-"), `[`, 2))
end_coords <- as.numeric(sapply(strsplit(scATAC_CD14_mono_peaks$Peak_Name, "-"), `[`, 3))
scATAC_CD14_mono_peaks$seqnames <- chromosomes
scATAC_CD14_mono_peaks$start <- start_coords
scATAC_CD14_mono_peaks$end <- end_coords

scATAC_CD14_mono_peaks_annotated <- annotatePeak(makeGRangesFromDataFrame(scATAC_CD14_mono_peaks), TxDb = txdb_hg38, annoDb = "org.Hs.eg.db")
scATAC_CD14_mono_peaks_annotated <- as.data.frame(scATAC_CD14_mono_peaks_annotated)
scATAC_CD14_mono_peaks_annotated$sc_log2FC <- scATAC_CD14_mono_peaks$sc_log2FC
scATAC_CD14_mono_peaks_annotated$sc_pval <- scATAC_CD14_mono_peaks$sc_pval

scATAC_CD14_mono_peaks_annotated_final <- scATAC_CD14_mono_peaks_annotated %>%
  filter(
    str_detect(SYMBOL, paste(gene_terms, collapse = "|")) |
      str_detect(GENENAME, paste(search_terms, collapse = "|"))
  )

scATAC_CD14_mono_peaks_annotated_final_pos <- scATAC_CD14_mono_peaks_annotated_final[scATAC_CD14_mono_peaks_annotated_final$sc_log2FC > 0,]
scATAC_CD14_mono_peaks_annotated_final_pos_promoters <- scATAC_CD14_mono_peaks_annotated_final_pos[scATAC_CD14_mono_peaks_annotated_final_pos$annotation %in% promoter_terms,]
scATAC_CD14_mono_peaks_annotated_final_neg <- scATAC_CD14_mono_peaks_annotated_final[scATAC_CD14_mono_peaks_annotated_final$sc_log2FC < 0,]
scATAC_CD14_mono_peaks_annotated_final_neg_promoters <- scATAC_CD14_mono_peaks_annotated_final_neg[scATAC_CD14_mono_peaks_annotated_final_neg$annotation %in% promoter_terms,]

# CD16 mono peaks
scATAC_CD16_mono_peaks <- read.table(paste0(scATAC_hvl_placebo_das_dir, "D28-vs-D_minus_1-degs-CD16_Mono-time_point-controlling_for_subject_id_final_pct_0.01.tsv"), sep = "\t", header = TRUE)

chromosomes <- sapply(strsplit(scATAC_CD16_mono_peaks$Peak_Name, "-"), `[`, 1)
start_coords <- as.numeric(sapply(strsplit(scATAC_CD16_mono_peaks$Peak_Name, "-"), `[`, 2))
end_coords <- as.numeric(sapply(strsplit(scATAC_CD16_mono_peaks$Peak_Name, "-"), `[`, 3))
scATAC_CD16_mono_peaks$seqnames <- chromosomes
scATAC_CD16_mono_peaks$start <- start_coords
scATAC_CD16_mono_peaks$end <- end_coords

scATAC_CD16_mono_peaks_annotated <- annotatePeak(makeGRangesFromDataFrame(scATAC_CD16_mono_peaks), TxDb = txdb_hg38, annoDb = "org.Hs.eg.db")
scATAC_CD16_mono_peaks_annotated <- as.data.frame(scATAC_CD16_mono_peaks_annotated)
scATAC_CD16_mono_peaks_annotated$sc_log2FC <- scATAC_CD16_mono_peaks$sc_log2FC
scATAC_CD16_mono_peaks_annotated$sc_pval <- scATAC_CD16_mono_peaks$sc_pval

scATAC_CD16_mono_peaks_annotated_final <- scATAC_CD16_mono_peaks_annotated %>%
  filter(
    str_detect(SYMBOL, paste(gene_terms, collapse = "|")) |
      str_detect(GENENAME, paste(search_terms, collapse = "|"))
  )

scATAC_CD16_mono_peaks_annotated_final_pos <- scATAC_CD16_mono_peaks_annotated_final[scATAC_CD16_mono_peaks_annotated_final$sc_log2FC > 0,]
scATAC_CD16_mono_peaks_annotated_final_pos_promoters <- scATAC_CD16_mono_peaks_annotated_final_pos[scATAC_CD16_mono_peaks_annotated_final_pos$annotation %in% promoter_terms,]
scATAC_CD16_mono_peaks_annotated_final_neg <- scATAC_CD16_mono_peaks_annotated_final[scATAC_CD16_mono_peaks_annotated_final$sc_log2FC < 0,]
scATAC_CD16_mono_peaks_annotated_final_neg_promoters <- scATAC_CD16_mono_peaks_annotated_final_neg[scATAC_CD16_mono_peaks_annotated_final_neg$annotation %in% promoter_terms,]

# cDC peaks
scATAC_cDC_peaks <- read.table(paste0(scATAC_hvl_placebo_das_dir, "D28-vs-D_minus_1-degs-cDC-time_point-controlling_for_subject_id_final_pct_0.01.tsv"), sep = "\t", header = TRUE)

chromosomes <- sapply(strsplit(scATAC_cDC_peaks$Peak_Name, "-"), `[`, 1)
start_coords <- as.numeric(sapply(strsplit(scATAC_cDC_peaks$Peak_Name, "-"), `[`, 2))
end_coords <- as.numeric(sapply(strsplit(scATAC_cDC_peaks$Peak_Name, "-"), `[`, 3))
scATAC_cDC_peaks$seqnames <- chromosomes
scATAC_cDC_peaks$start <- start_coords
scATAC_cDC_peaks$end <- end_coords

scATAC_cDC_peaks_annotated <- annotatePeak(makeGRangesFromDataFrame(scATAC_cDC_peaks), TxDb = txdb_hg38, annoDb = "org.Hs.eg.db")
scATAC_cDC_peaks_annotated <- as.data.frame(scATAC_cDC_peaks_annotated)
scATAC_cDC_peaks_annotated$sc_log2FC <- scATAC_cDC_peaks$sc_log2FC
scATAC_cDC_peaks_annotated$sc_pval <- scATAC_cDC_peaks$sc_pval

scATAC_cDC_peaks_annotated_final <- scATAC_cDC_peaks_annotated %>%
  filter(
    str_detect(SYMBOL, paste(gene_terms, collapse = "|")) |
      str_detect(GENENAME, paste(search_terms, collapse = "|"))
  )

scATAC_cDC_peaks_annotated_final_pos <- scATAC_cDC_peaks_annotated_final[scATAC_cDC_peaks_annotated_final$sc_log2FC > 0,]
scATAC_cDC_peaks_annotated_final_pos_promoters <- scATAC_cDC_peaks_annotated_final_pos[scATAC_cDC_peaks_annotated_final_pos$annotation %in% promoter_terms,]
scATAC_cDC_peaks_annotated_final_neg <- scATAC_cDC_peaks_annotated_final[scATAC_cDC_peaks_annotated_final$sc_log2FC < 0,]
scATAC_cDC_peaks_annotated_final_neg_promoters <- scATAC_cDC_peaks_annotated_final_neg[scATAC_cDC_peaks_annotated_final_neg$annotation %in% promoter_terms,]

# NK peaks
scATAC_NK_peaks <- read.table(paste0(scATAC_hvl_placebo_das_dir, "D28-vs-D_minus_1-degs-NK-time_point-controlling_for_subject_id_final_pct_0.01.tsv"), sep = "\t", header = TRUE)

chromosomes <- sapply(strsplit(scATAC_NK_peaks$Peak_Name, "-"), `[`, 1)
start_coords <- as.numeric(sapply(strsplit(scATAC_NK_peaks$Peak_Name, "-"), `[`, 2))
end_coords <- as.numeric(sapply(strsplit(scATAC_NK_peaks$Peak_Name, "-"), `[`, 3))
scATAC_NK_peaks$seqnames <- chromosomes
scATAC_NK_peaks$start <- start_coords
scATAC_NK_peaks$end <- end_coords

scATAC_NK_peaks_annotated <- annotatePeak(makeGRangesFromDataFrame(scATAC_NK_peaks), TxDb = txdb_hg38, annoDb = "org.Hs.eg.db")
scATAC_NK_peaks_annotated <- as.data.frame(scATAC_NK_peaks_annotated)
scATAC_NK_peaks_annotated$sc_log2FC <- scATAC_NK_peaks$sc_log2FC
scATAC_NK_peaks_annotated$sc_pval <- scATAC_NK_peaks$sc_pval

scATAC_NK_peaks_annotated_final <- scATAC_NK_peaks_annotated %>%
  filter(
    str_detect(SYMBOL, paste(gene_terms, collapse = "|")) |
      str_detect(GENENAME, paste(search_terms, collapse = "|"))
  )

scATAC_NK_peaks_annotated_final_pos <- scATAC_NK_peaks_annotated_final[scATAC_NK_peaks_annotated_final$sc_log2FC > 0,]
scATAC_NK_peaks_annotated_final_pos_promoters <- scATAC_NK_peaks_annotated_final_pos[scATAC_NK_peaks_annotated_final_pos$annotation %in% promoter_terms,]
scATAC_NK_peaks_annotated_final_neg <- scATAC_NK_peaks_annotated_final[scATAC_NK_peaks_annotated_final$sc_log2FC < 0,]
scATAC_NK_peaks_annotated_final_neg_promoters <- scATAC_NK_peaks_annotated_final_neg[scATAC_NK_peaks_annotated_final_neg$annotation %in% promoter_terms,]



