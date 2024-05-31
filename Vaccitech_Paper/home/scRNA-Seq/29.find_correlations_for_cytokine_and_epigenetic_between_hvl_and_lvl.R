# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# I picked my 9 favorite cell types so I can have a 3x3 grid of correlation plots
correlation_cell_types <- c("CD16_Mono", "CD14_Mono", "NK", "NK_CD56bright", "cDC", "pDC")

sc_correlations <- list()
sc_correlation_plots <- list()

cytokine_genes <- c("CCL3", "CX3CR1", "CXCL16", "NFIL3", "IL32", "IRAK3", "IL1RAP", "IRF2", "IRF7", "IFNGR1", "IFNG", "JUN", "JUNB", 
                   "JUND", "FOS", "FOSL1", "FOSL2")

epigenetic_genes <- c("KAT6A", "HDAC5", "HDAC7", "HDAC9", "ASH1L", "PRDM2", "SETD2", "KDM2B", "KDM3A")

for(cell_type in correlation_cell_types) {
  print(cell_type)
  sc_correlations[[cell_type]] <- list()
  sc_correlation_plots[[cell_type]] <- list()
  cell_type_no_underscore <- gsub("_", " ", cell_type)
  # Unfiltered SC
  unfiltered_cell_type_sc_degs <- read.table(paste0(scRNA_hvl_placebo_deg_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_sc_unfiltered.tsv"),
                                             sep = "\t", header = TRUE)
  unfiltered_cell_type_lvl_sc_degs <- read.table(paste0(scRNA_lvl_placebo_deg_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_sc_unfiltered.tsv"),
                                                        sep = "\t", header = TRUE)
  
  # Check correlation for cytokine genes
  cytokine_unfiltered_cell_type_sc_degs <- unfiltered_cell_type_sc_degs[rownames(unfiltered_cell_type_sc_degs) %in% cytokine_genes,]
  cytokine_unfiltered_cell_type_lvl_sc_degs <- unfiltered_cell_type_lvl_sc_degs[rownames(unfiltered_cell_type_lvl_sc_degs) %in% cytokine_genes,]
  
  cytokine_unfiltered_cell_type_sc_degs <- cytokine_unfiltered_cell_type_sc_degs[order(rownames(cytokine_unfiltered_cell_type_sc_degs)),]
  cytokine_unfiltered_cell_type_lvl_sc_degs <- cytokine_unfiltered_cell_type_lvl_sc_degs[order(rownames(cytokine_unfiltered_cell_type_lvl_sc_degs)),]
  
  comparing_first_vs_second_df <- data.frame(gene_name = rownames(cytokine_unfiltered_cell_type_sc_degs), first_fc = cytokine_unfiltered_cell_type_sc_degs$avg_log2FC,
                                             second_fc = cytokine_unfiltered_cell_type_lvl_sc_degs$avg_log2FC)
  # Calculate correlation
  correlation_val <- cor.test(comparing_first_vs_second_df$first_fc, comparing_first_vs_second_df$second_fc)
  print(correlation_val$estimate)
  print(correlation_val$p.value)
  sc_correlations[[cell_type]][["sc_primary_vs_sc_lvl"]] <- correlation_val
  
  # Check correlation for epigenetic genes
  epigenetic_unfiltered_cell_type_sc_degs <- unfiltered_cell_type_sc_degs[rownames(unfiltered_cell_type_sc_degs) %in% epigenetic_genes,]
  epigenetic_unfiltered_cell_type_lvl_sc_degs <- unfiltered_cell_type_lvl_sc_degs[rownames(unfiltered_cell_type_lvl_sc_degs) %in% epigenetic_genes,]
  
  epigenetic_unfiltered_cell_type_sc_degs <- epigenetic_unfiltered_cell_type_sc_degs[order(rownames(epigenetic_unfiltered_cell_type_sc_degs)),]
  epigenetic_unfiltered_cell_type_lvl_sc_degs <- epigenetic_unfiltered_cell_type_lvl_sc_degs[order(rownames(epigenetic_unfiltered_cell_type_lvl_sc_degs)),]
  
  comparing_first_vs_second_df <- data.frame(gene_name = rownames(epigenetic_unfiltered_cell_type_sc_degs), first_fc = epigenetic_unfiltered_cell_type_sc_degs$avg_log2FC,
                                             second_fc = epigenetic_unfiltered_cell_type_lvl_sc_degs$avg_log2FC)
  # Calculate correlation
  correlation_val <- cor.test(comparing_first_vs_second_df$first_fc, comparing_first_vs_second_df$second_fc)
  print(correlation_val$estimate)
  print(correlation_val$p.value)
  
  # Plot correlation
  sc_correlation_plots[[cell_type]][["sc_primary_vs_sc_lvl"]] <- ggplot(data = comparing_first_vs_second_df, mapping = aes(x = first_fc, y = second_fc)) +
    geom_point(size = 2) +
    sm_statCorr() + xlab("High Viral Load FC") + ylab("Low Viral Load FC") + labs(title = cell_type_no_underscore)
  
  # Check correlation for primary significant SC genes (pseudobulk_corrected) in LVL data
  primary_sc_pseudobulk_corrected_degs <- rownames(cell_type_sc_pseudobulk_corrected_degs)
  print(paste0("Number of SC DEGs (pseudobulk corrected) in primary data is: ", length(primary_sc_pseudobulk_corrected_degs)))
  compare_first_df <- cell_type_sc_pseudobulk_corrected_degs
  compare_second_df <- unfiltered_cell_type_lvl_sc_degs[rownames(unfiltered_cell_type_lvl_sc_degs) %in% primary_sc_pseudobulk_corrected_degs,]
  compare_first_df <- compare_first_df[rownames(compare_first_df) %in% rownames(compare_second_df),]
  compare_second_df <- compare_second_df[rownames(compare_second_df) %in% rownames(compare_first_df),]
  compare_first_df <- compare_first_df[order(rownames(compare_first_df)),]
  compare_second_df <- compare_second_df[order(rownames(compare_second_df)),]
  
  comparing_first_vs_second_df <- data.frame(gene_name = rownames(compare_first_df), first_fc = compare_first_df$sc_log2FC,
                                             second_fc = compare_second_df$avg_log2FC)
  # Calculate correlation
  correlation_val <- cor.test(comparing_first_vs_second_df$first_fc, comparing_first_vs_second_df$second_fc)
  print(correlation_val$estimate)
  print(correlation_val$p.value)
  sc_correlations[[cell_type]][["sc_primary_pseudobulk_corrected_vs_sc_lvl"]] <- correlation_val
  
  # Plot correlation
  sc_correlation_plots[[cell_type]][["sc_primary_pseudobulk_corrected_vs_sc_lvl"]] <- ggplot(data = comparing_first_vs_second_df, mapping = aes(x = first_fc, y = second_fc)) +
    geom_point(size = 2) +
    sm_statCorr(corr_method = "spearman") + xlab("HVL Naive FC") + ylab("LVL Naive FC") + labs(title = cell_type_no_underscore) +  xlim(-1.5, 1.5) + ylim(-1.5, 1.5)
}

hvl_vs_lvl_pseudobulk_corrected_plots <- lapply(sc_correlation_plots, function(x) x[[2]])
ggsave("C:/Users/willi/Desktop/hvl_vs_lvl_fc_correlation.png", plot = patchwork::wrap_plots(hvl_vs_lvl_pseudobulk_corrected_plots, ncol = 3, nrow = 3), height = 10, width = 10)

# Other plotting attempts
#n <- length(pseudobulk_corrected_plots)
#nCol <- floor(sqrt(n))
# do.call("grid.arrange", c(pseudobulk_corrected_plots, ncol=nCol))
# test <- list(pseudobulk_corrected_plots[[1]], pseudobulk_corrected_plots[[2]], pseudobulk_corrected_plots[[3]])
# do.call("grid.arrange", c(test, ncol=nCol))