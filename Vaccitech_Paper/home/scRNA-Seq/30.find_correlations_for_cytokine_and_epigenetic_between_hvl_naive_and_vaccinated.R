# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

correlation_cell_types <- c("CD14_Mono", "CD16_Mono", "NK", "NK_CD56bright", "cDC", "pDC")

hvl_vs_vaccinated_cytoepi_correlations <- list()
hvl_vs_vaccinated_cytoepi_correlation_plots <- list()

cytokine_genes <- c("CCL3", "CX3CR1", "CXCL16", "NFIL3", "IL32", "IRAK3", "IL1RAP", "IRF2", "IRF7", "IFNGR1", "IFNG", "JUN", "JUNB", 
                   "JUND", "FOS", "FOSL1", "FOSL2")

epigenetic_genes <- c("KAT6A", "HDAC5", "HDAC7", "HDAC9", "ASH1L", "PRDM2", "SETD2", "KDM2B", "KDM3A")

for(cell_type in correlation_cell_types) {
  print(cell_type)
  hvl_vs_vaccinated_cytoepi_correlations[[cell_type]] <- list()
  hvl_vs_vaccinated_cytoepi_correlation_plots[[cell_type]] <- list()
  cell_type_no_underscore <- gsub("_", " ", cell_type)
  if(cell_type_no_underscore == "NK CD56bright") {
    cell_type_no_underscore <- "NK_CD56bright"
  }
  # Unfiltered SC
  unfiltered_cell_type_sc_degs <- read.table(paste0(scRNA_hvl_placebo_deg_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_sc_unfiltered.tsv"),
                                             sep = "\t", header = TRUE)
  unfiltered_cell_type_vaccinated_sc_degs <- read.table(paste0(scRNA_hvl_vaccinated_deg_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_sc_unfiltered.tsv"),
                                                        sep = "\t", header = TRUE)
  
  # Check correlation for cytokine genes
  cytokine_unfiltered_cell_type_sc_degs <- unfiltered_cell_type_sc_degs[rownames(unfiltered_cell_type_sc_degs) %in% cytokine_genes,]
  cytokine_unfiltered_cell_type_vaccinated_sc_degs <- unfiltered_cell_type_vaccinated_sc_degs[rownames(unfiltered_cell_type_vaccinated_sc_degs) %in% cytokine_genes,]
  
  cytokine_unfiltered_cell_type_sc_degs <- cytokine_unfiltered_cell_type_sc_degs[order(rownames(cytokine_unfiltered_cell_type_sc_degs)),]
  cytokine_unfiltered_cell_type_vaccinated_sc_degs <- cytokine_unfiltered_cell_type_vaccinated_sc_degs[order(rownames(cytokine_unfiltered_cell_type_vaccinated_sc_degs)),]
  
  comparing_first_vs_second_df <- data.frame(gene_name = rownames(cytokine_unfiltered_cell_type_sc_degs), first_fc = cytokine_unfiltered_cell_type_sc_degs$avg_log2FC,
                                             second_fc = cytokine_unfiltered_cell_type_vaccinated_sc_degs$avg_log2FC)
  # Calculate correlation
  correlation_val <- cor.test(comparing_first_vs_second_df$first_fc, comparing_first_vs_second_df$second_fc)
  print(correlation_val$estimate)
  print(correlation_val$p.value)
  hvl_vs_vaccinated_cytoepi_correlations[[cell_type]][["cytokine"]] <- comparing_first_vs_second_df
  
  # Plot correlation
  hvl_vs_vaccinated_cytoepi_correlation_plots[[cell_type]][["cytokine"]] <- ggplot(data = comparing_first_vs_second_df, mapping = aes(x = first_fc, y = second_fc)) +
    geom_point(size = 2) +
    sm_statCorr(corr_method = "spearman") + xlab("High Viral Load FC") + ylab("Low Viral Load FC") + labs(title = cell_type_no_underscore) + xlim(-2, 2) + ylim(-2, 2)
  
  # Check correlation for epigenetic genes
  epigenetic_unfiltered_cell_type_sc_degs <- unfiltered_cell_type_sc_degs[rownames(unfiltered_cell_type_sc_degs) %in% epigenetic_genes,]
  epigenetic_unfiltered_cell_type_vaccinated_sc_degs <- unfiltered_cell_type_vaccinated_sc_degs[rownames(unfiltered_cell_type_vaccinated_sc_degs) %in% epigenetic_genes,]
  
  epigenetic_unfiltered_cell_type_sc_degs <- epigenetic_unfiltered_cell_type_sc_degs[order(rownames(epigenetic_unfiltered_cell_type_sc_degs)),]
  epigenetic_unfiltered_cell_type_vaccinated_sc_degs <- epigenetic_unfiltered_cell_type_vaccinated_sc_degs[order(rownames(epigenetic_unfiltered_cell_type_vaccinated_sc_degs)),]
  
  comparing_first_vs_second_df <- data.frame(gene_name = rownames(epigenetic_unfiltered_cell_type_sc_degs), first_fc = epigenetic_unfiltered_cell_type_sc_degs$avg_log2FC,
                                             second_fc = epigenetic_unfiltered_cell_type_vaccinated_sc_degs$avg_log2FC)
  # Calculate correlation
  correlation_val <- cor.test(comparing_first_vs_second_df$first_fc, comparing_first_vs_second_df$second_fc)
  print(correlation_val$estimate)
  print(correlation_val$p.value)
  hvl_vs_vaccinated_cytoepi_correlations[[cell_type]][["epigenetic"]] <- comparing_first_vs_second_df
  
  # Plot correlation
  hvl_vs_vaccinated_cytoepi_correlation_plots[[cell_type]][["epigenetic"]] <- ggplot(data = comparing_first_vs_second_df, mapping = aes(x = first_fc, y = second_fc)) +
    geom_point(size = 2) +
    sm_statCorr(corr_method = "spearman") + xlab("High Viral Load FC") + ylab("Low Viral Load FC") + labs(title = cell_type_no_underscore) + xlim(-2, 2) + ylim(-2, 2)
  
}

cyto_hvl_vs_vaccinated_cytoepi_correlation_plots <- lapply(hvl_vs_vaccinated_cytoepi_correlation_plots, function(x) x[[1]])
ggsave("C:/Users/willi/Desktop/hvl_vs_vaccinated_fc_cyto_correlation.png", plot = patchwork::wrap_plots(cyto_hvl_vs_vaccinated_cytoepi_correlation_plots, ncol = 2, nrow = 3), height = 10, width = 10)

epigenetic_hvl_vs_vaccinated_cytoepi_correlation_plots <- lapply(hvl_vs_vaccinated_cytoepi_correlation_plots, function(x) x[[2]])
ggsave("C:/Users/willi/Desktop/hvl_vs_vaccinated_fc_epigenetic_correlation.png", plot = patchwork::wrap_plots(epigenetic_hvl_vs_vaccinated_cytoepi_correlation_plots, ncol = 2, nrow = 3), height = 10, width = 10)


# Other plotting attempts
#n <- length(pseudobulk_corrected_plots)
#nCol <- floor(sqrt(n))
# do.call("grid.arrange", c(pseudobulk_corrected_plots, ncol=nCol))
# test <- list(pseudobulk_corrected_plots[[1]], pseudobulk_corrected_plots[[2]], pseudobulk_corrected_plots[[3]])
# do.call("grid.arrange", c(test, ncol=nCol))