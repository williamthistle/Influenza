# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

find_atac_and_rna_correlation <- function(atac_df, rna_df) {
  rna_df <- rna_df[rownames(rna_df) %in% atac_df$SYMBOL,]
  atac_df <- atac_df[atac_df$SYMBOL %in% rownames(rna_df),]
  
  # Grab FC associated with each DAS symbol
  corresponding_rna_fc <- c()
  
  for(current_row_index in 1:nrow(atac_df)) {
    current_symbol <- atac_df[current_row_index,]$SYMBOL
    current_rna_fc <- rna_df[rownames(rna_df) == current_symbol,]$avg_log2FC
    corresponding_rna_fc <- c(corresponding_rna_fc, current_rna_fc)
  }
  
  comparing_first_vs_second_df <- data.frame(gene_name = atac_df$SYMBOL, first_fc = atac_df$sc_log2FC,
                                             second_fc = corresponding_rna_fc)
  return(comparing_first_vs_second_df)
}

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

# Get AP-1 loci
peak_motif_matches <- read.table(paste0(scATAC_hvl_dir, "Placebo/peak_motif_matches.txt"), sep = "\t", header = TRUE)

peak_motif_matches_AP1 <- peak_motif_matches[peak_motif_matches$FOS > 0 | peak_motif_matches$FOSL1 > 0 | peak_motif_matches$FOSL2 > 0
                                             | peak_motif_matches$JUN > 0 | peak_motif_matches$JUNB > 0 | peak_motif_matches$JUND > 0,]
peak_motif_matches_AP1$chr <- paste0("chr", peak_motif_matches_AP1$chr)
peak_motif_matches_AP1$chr <- gsub("chr23", "chrX", peak_motif_matches_AP1$chr)

# Get IRF1 loci
peak_motif_matches <- read.table(paste0(scATAC_hvl_dir, "Placebo/peak_motif_matches.txt"), sep = "\t", header = TRUE)

peak_motif_matches_IRF1 <- peak_motif_matches[peak_motif_matches$IRF1 > 0,]
peak_motif_matches_IRF1$chr <- paste0("chr", peak_motif_matches_IRF1$chr)
peak_motif_matches_IRF1$chr <- gsub("chr23", "chrX", peak_motif_matches_IRF1$chr)

correlation_cell_types <- c("CD14_Mono", "CD16_Mono", "cDC", "pDC", "NK")

sc_annotated_genes <- list()
sc_correlations <- list()
sc_correlation_plots <- list()

for(cell_type in correlation_cell_types) {
  print(cell_type)
  sc_correlations[[cell_type]] <- list()
  sc_correlation_plots[[cell_type]] <- list()
  cell_type_no_underscore <- gsub("_", " ", cell_type)
  # Unfiltered scRNA
  unfiltered_cell_type_sc_deg <- read.table(paste0(scRNA_hvl_placebo_deg_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_sc_unfiltered.tsv"),
                                            sep = "\t", header = TRUE)
  
  # Get AP-1 loci for current cell type
  scATAC_peaks <- read.table(paste0(scATAC_hvl_placebo_das_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_final_pct_0.01.tsv"), sep = "\t", header = TRUE)
  
  chromosomes <- sapply(strsplit(scATAC_peaks$Peak_Name, "-"), `[`, 1)
  start_coords <- as.numeric(sapply(strsplit(scATAC_peaks$Peak_Name, "-"), `[`, 2))
  end_coords <- as.numeric(sapply(strsplit(scATAC_peaks$Peak_Name, "-"), `[`, 3))
  
  scATAC_peaks$seqnames <- chromosomes
  scATAC_peaks$start <- start_coords
  scATAC_peaks$end <- end_coords
  
  scATAC_peaks_AP1 <- scATAC_peaks[scATAC_peaks$seqnames %in% peak_motif_matches_AP1$chr &
                                     scATAC_peaks$start %in% peak_motif_matches_AP1$point1 &
                                     scATAC_peaks$end %in% peak_motif_matches_AP1$point2,]
  
  # Annotate with genes
  scATAC_peaks_AP1_annotated <- annotatePeak(makeGRangesFromDataFrame(scATAC_peaks_AP1), TxDb = txdb_hg38, annoDb = "org.Hs.eg.db")
  scATAC_peaks_AP1_annotated <- as.data.frame(scATAC_peaks_AP1_annotated)
  scATAC_peaks_AP1_annotated$sc_log2FC <- scATAC_peaks_AP1$sc_log2FC
  scATAC_peaks_AP1_annotated$sc_pval <- scATAC_peaks_AP1$sc_pval
  
  sc_annotated_genes[[cell_type]] <- scATAC_peaks_AP1_annotated
  
  # Check correlation for primary significant SC genes in validation set
  comparing_first_vs_second_df <- find_atac_and_rna_correlation(scATAC_peaks_AP1_annotated, unfiltered_cell_type_sc_deg)
  
  # Calculate correlation
  correlation_val <- cor.test(comparing_first_vs_second_df$first_fc, comparing_first_vs_second_df$second_fc,
                              method = "spearman")
  print(correlation_val$estimate)
  print(correlation_val$p.value)
  sc_correlations[[cell_type]][["AP1_atac_vs_rna"]] <- correlation_val
  
  # Plot correlation
  sc_correlation_plots[[cell_type]][["AP1_atac_vs_rna"]] <- ggplot(data = comparing_first_vs_second_df, mapping = aes(x = first_fc, y = second_fc)) +
    geom_point(size = 2) +
    sm_statCorr(corr_method = "spearman") + xlab("ATAC FC") + ylab("RNA FC") + labs(title = cell_type_no_underscore) + xlim(-6, 6) + ylim(-6, 6)
}

pseudobulk_corrected_plots <- lapply(sc_correlation_plots, function(x) x[[1]])
ggsave("C:/Users/wat2/Desktop/AP1_atac_vs_rna_correlations.png", plot = patchwork::wrap_plots(pseudobulk_corrected_plots, ncol = 2, nrow = 3), height = 10, width = 10)

# CD14 Mono
scATAC_CD14_mono_peaks_AP1_annotated <- sc_annotated_genes[["CD14_Mono"]]
scATAC_CD14_mono_peaks_AP1_annotated_pos <- scATAC_CD14_mono_peaks_AP1_annotated[scATAC_CD14_mono_peaks_AP1_annotated$sc_log2FC > 0,]
scATAC_CD14_mono_peaks_AP1_annotated_neg <- scATAC_CD14_mono_peaks_AP1_annotated[scATAC_CD14_mono_peaks_AP1_annotated$sc_log2FC < 0,]

scATAC_CD14_mono_peaks_AP1_annotated_final <- scATAC_CD14_mono_peaks_AP1_annotated %>%
  filter(
    str_detect(SYMBOL, paste(gene_terms, collapse = "|")) |
      str_detect(GENENAME, paste(search_terms, collapse = "|"))
  )

scATAC_CD14_mono_peaks_AP1_annotated_final_pos <- scATAC_CD14_mono_peaks_AP1_annotated_final[scATAC_CD14_mono_peaks_AP1_annotated_final$sc_log2FC > 0,]
scATAC_CD14_mono_peaks_AP1_annotated_final_neg <- scATAC_CD14_mono_peaks_AP1_annotated_final[scATAC_CD14_mono_peaks_AP1_annotated_final$sc_log2FC < 0,]

# IRF1
for(cell_type in correlation_cell_types) {
  print(cell_type)
  sc_correlations[[cell_type]] <- list()
  sc_correlation_plots[[cell_type]] <- list()
  cell_type_no_underscore <- gsub("_", " ", cell_type)
  # Unfiltered scRNA
  unfiltered_cell_type_sc_deg <- read.table(paste0(scRNA_hvl_placebo_deg_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_sc_unfiltered.tsv"),
                                            sep = "\t", header = TRUE)
  
  # Get AP-1 loci for current cell type
  scATAC_peaks <- read.table(paste0(scATAC_hvl_placebo_das_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_final_pct_0.01.tsv"), sep = "\t", header = TRUE)
  
  chromosomes <- sapply(strsplit(scATAC_peaks$Peak_Name, "-"), `[`, 1)
  start_coords <- as.numeric(sapply(strsplit(scATAC_peaks$Peak_Name, "-"), `[`, 2))
  end_coords <- as.numeric(sapply(strsplit(scATAC_peaks$Peak_Name, "-"), `[`, 3))
  
  scATAC_peaks$seqnames <- chromosomes
  scATAC_peaks$start <- start_coords
  scATAC_peaks$end <- end_coords
  
  scATAC_peaks_IRF1 <- scATAC_peaks[scATAC_peaks$seqnames %in% peak_motif_matches_IRF1$chr &
                                     scATAC_peaks$start %in% peak_motif_matches_IRF1$point1 &
                                     scATAC_peaks$end %in% peak_motif_matches_IRF1$point2,]
  
  # Check number of upregulated and downregulated loci
  nrow(scATAC_peaks_IRF1[scATAC_peaks_IRF1$sc_log2FC > 0,])
  nrow(scATAC_peaks_IRF1[scATAC_peaks_IRF1$sc_log2FC < 0,])
  
  # Now, let's do SC (no pseudobulk correction)
  #scATAC_peaks <- read.table(paste0(scATAC_hvl_placebo_das_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_sc_pct_0.01.tsv"), sep = "\t", header = TRUE)
  
  #chromosomes <- sapply(strsplit(rownames(scATAC_peaks), "-"), `[`, 1)
  #start_coords <- as.numeric(sapply(strsplit(rownames(scATAC_peaks), "-"), `[`, 2))
  #end_coords <- as.numeric(sapply(strsplit(rownames(scATAC_peaks), "-"), `[`, 3))
  
  #scATAC_peaks$seqnames <- chromosomes
  #scATAC_peaks$start <- start_coords
  #scATAC_peaks$end <- end_coords
  
  #scATAC_peaks_IRF1 <- scATAC_peaks[scATAC_peaks$seqnames %in% peak_motif_matches_IRF1$chr &
  #                                    scATAC_peaks$start %in% peak_motif_matches_IRF1$point1 &
  #                                    scATAC_peaks$end %in% peak_motif_matches_IRF1$point2,]
  
  # Check number of upregulated and downregulated loci
  #nrow(scATAC_peaks_IRF1[scATAC_peaks_IRF1$avg_log2FC > 0,])
  #nrow(scATAC_peaks_IRF1[scATAC_peaks_IRF1$avg_log2FC < 0,])
  
  # Annotate with genes
  scATAC_peaks_IRF1_annotated <- annotatePeak(makeGRangesFromDataFrame(scATAC_peaks_IRF1), TxDb = txdb_hg38, annoDb = "org.Hs.eg.db")
  scATAC_peaks_IRF1_annotated <- as.data.frame(scATAC_peaks_IRF1_annotated)
  scATAC_peaks_IRF1_annotated$sc_log2FC <- scATAC_peaks_IRF1$sc_log2FC
  scATAC_peaks_IRF1_annotated$sc_pval <- scATAC_peaks_IRF1$sc_pval
  
  sc_annotated_genes[[cell_type]] <- scATAC_peaks_AP1_annotated
  
  # Check correlation for primary significant SC genes in validation set
  comparing_first_vs_second_df <- find_atac_and_rna_correlation(scATAC_peaks_AP1_annotated, unfiltered_cell_type_sc_deg)
  
  # Calculate correlation
  correlation_val <- cor.test(comparing_first_vs_second_df$first_fc, comparing_first_vs_second_df$second_fc,
                              method = "spearman")
  print(correlation_val$estimate)
  print(correlation_val$p.value)
  sc_correlations[[cell_type]][["AP1_atac_vs_rna"]] <- correlation_val
  
  # Plot correlation
  sc_correlation_plots[[cell_type]][["AP1_atac_vs_rna"]] <- ggplot(data = comparing_first_vs_second_df, mapping = aes(x = first_fc, y = second_fc)) +
    geom_point(size = 2) +
    sm_statCorr(corr_method = "spearman") + xlab("ATAC FC") + ylab("RNA FC") + labs(title = cell_type_no_underscore) + xlim(-6, 6) + ylim(-6, 6)
}
