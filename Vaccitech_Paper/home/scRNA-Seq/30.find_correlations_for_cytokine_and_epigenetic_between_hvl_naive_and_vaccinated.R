# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

correlation_cell_types <- c("pDC", "NK_CD56bright", "CD16_Mono", "cDC", "NK", "CD14_Mono")

hvl_vs_lvl_cytoepi_correlations <- list()
hvl_vs_lvl_cytoepi_correlation_plots <- list()

cytokine_genes <- c("CCL3", "CX3CR1", "CCL3L1", "CXCL16", 
                    "IL32", "CASP1", "CSF1R", "NFIL3", "IRAK3", "IL1RAP", "PTGES",   
                    "IRF2", "IRF7", "IFNG", "OAS1", "MNDA", "PSMB9", "IFNGR1", "DNAJC3", "GBP5", "USP38",  
                    "JUN", "JUNB", "FOSB", "FOSL2", "JDP2", 
                    "MAP3K11", "CSK", "DUSP1", "DUSP2", "DUSP6", "ABHD17A", "RIPK1", "MAPK7", "MAPK8", "MAP2K1", "MAP3K8", "MAP3K20", "MAP4K3", "MAPKAPK2",
                    "PTK2B", "RELL1", "MINK1", "BRAF",
                    "NMRAL1", "NFKBIA", "NFKBIZ", "NFKB1",
                    "JAK1", "STAT3", "STAT4")

cytokine_genes_hvl_vaccinated <- c("CX3CR1", "CCR2", "CCL3", "CXCL16", "CXCL2", "CXCL3", "CXCL8", "CXCR4",
                                   "CSF1R", "NFIL3", "IL10", "IL10RA", "IL1R1", "IL1R2", "IL1RAP",
                                   "IFNG", "IRF1", "IRF2", "OAS1", "IFI16", "IFIH1", "IFITM1", "IFNAR1", "IRF8", "IFNGR1", "IRF2BP2", 
                                   "FOS", "JUN", "FOSL2", "JDP2",
                                   "DUSP1", "MAPK14", "MAP2K3", "MAP3K8", "MAPK6", "MAPK8", "MAPKAPK2", "DUSP4",
                                   "NFKBIA", "NFKBIZ",
                                   "JAK2", "STAT1")

cytokine_genes <- union(cytokine_genes, cytokine_genes_hvl_vaccinated)

epigenetic_genes <- c("HIST1H1C", "HIST1H1D", "HIST1H1E", "H2AFZ", "KAT6A", "BRD1", "MSL2", "HDAC5", "HDAC7", "HDAC9", "METTL23", "ASH1L", "PRDM2", "SETD2","KDM2B", "KDM3A")

epigenetic_genes_hvl_vaccinated <- c("HDAC2", "HDAC9", "SETD3", "KAT6B")
  
epigenetic_genes <- union(epigenetic_genes, epigenetic_genes_hvl_vaccinated)

for(cell_type in correlation_cell_types) {
  print(cell_type)
  hvl_vs_lvl_cytoepi_correlations[[cell_type]] <- list()
  hvl_vs_lvl_cytoepi_correlation_plots[[cell_type]] <- list()
  cell_type_no_underscore <- gsub("_", " ", cell_type)
  if(cell_type_no_underscore == "NK CD56bright") {
    cell_type_no_underscore <- "NK_CD56bright"
  }
  # Unfiltered SC
  unfiltered_cell_type_sc_degs <- read.table(paste0(scRNA_hvl_placebo_deg_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_sc_unfiltered.tsv"),
                                             sep = "\t", header = TRUE)
  unfiltered_cell_type_lvl_sc_degs <- read.table(paste0(scRNA_hvl_vaccinated_deg_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_sc_unfiltered.tsv"),
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
  hvl_vs_lvl_cytoepi_correlations[[cell_type]][["cytokine"]] <- comparing_first_vs_second_df
  
  # Plot correlation
  hvl_vs_lvl_cytoepi_correlation_plots[[cell_type]][["cytokine"]] <- ggplot(data = comparing_first_vs_second_df, mapping = aes(x = second_fc, y = first_fc)) +
    geom_point(size = 2) +
    sm_statCorr(corr_method = "spearman", text_size = 6) + xlab("Vaccinated HVL FC") + ylab("HVL HVL FC") + labs(title = cell_type_no_underscore) + 
    xlim(-3, 3) + ylim(-3, 3) + theme(aspect.ratio = 1, text=element_text(size=15)) + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")
  
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
  hvl_vs_lvl_cytoepi_correlations[[cell_type]][["epigenetic"]] <- comparing_first_vs_second_df
  
  # Plot correlation
  hvl_vs_lvl_cytoepi_correlation_plots[[cell_type]][["epigenetic"]] <- ggplot(data = comparing_first_vs_second_df, mapping = aes(x = second_fc, y = first_fc)) +
    geom_point(size = 2) +
    sm_statCorr(corr_method = "spearman", text_size = 6) + xlab("Vaccinated HVL FC") + ylab("Naive HVL FC") + labs(title = cell_type_no_underscore) + 
    xlim(-2, 2) + ylim(-2, 2) + theme(aspect.ratio = 1, text=element_text(size=15)) + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")
}

cyto_hvl_vs_lvl_cytoepi_correlation_plots <- lapply(hvl_vs_lvl_cytoepi_correlation_plots, function(x) x[[1]])
ggsave("C:/Users/willi/Desktop/hvl_vs_vaccinated_fc_cyto_correlation.png", plot = patchwork::wrap_plots(cyto_hvl_vs_lvl_cytoepi_correlation_plots, ncol = 3, nrow = 2), height = 10, width = 10)

epigenetic_hvl_vs_lvl_cytoepi_correlation_plots <- lapply(hvl_vs_lvl_cytoepi_correlation_plots, function(x) x[[2]])
ggsave("C:/Users/willi/Desktop/hvl_vs_vaccinated_fc_epigenetic_correlation.png", plot = patchwork::wrap_plots(epigenetic_hvl_vs_lvl_cytoepi_correlation_plots, ncol = 3, nrow = 2), height = 10, width = 10)


# Other plotting attempts
#n <- length(pseudobulk_corrected_plots)
#nCol <- floor(sqrt(n))
# do.call("grid.arrange", c(pseudobulk_corrected_plots, ncol=nCol))
# test <- list(pseudobulk_corrected_plots[[1]], pseudobulk_corrected_plots[[2]], pseudobulk_corrected_plots[[3]])
# do.call("grid.arrange", c(test, ncol=nCol))