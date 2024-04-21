all_sc_das_pseudobulk_corrected_table <- data.frame(Cell_Type = character(), Site = character(), log2FC = numeric())

for(cell_type in scATAC_cell_types) {
  cell_type_for_file_name <- sub(" ", "_", cell_type)
  full_sc_das <- read.table(paste0(scATAC_hvl_placebo_das_dir, "D28-vs-D_minus_1-degs-", cell_type_for_file_name, "-time_point-controlling_for_subject_id_final_pct_0.01.tsv"),
                                    sep = "\t", header = TRUE)
  all_sc_das_table_for_cell_type <- data.frame(Cell_Type = cell_type, Site = full_sc_das$Peak_Name, log2FC = full_sc_das$sc_log2FC)
  all_sc_das_pseudobulk_corrected_table <- rbind(all_sc_das_pseudobulk_corrected_table, all_sc_das_table_for_cell_type)
}

all_sc_das_pseudobulk_corrected_table <- all_sc_das_pseudobulk_corrected_table  %>%
  mutate(direction = ifelse(log2FC > 0, "Upregulated", "Downregulated"))

pos_all_sc_das_pseudobulk_corrected <- all_sc_das_pseudobulk_corrected_table[all_sc_das_pseudobulk_corrected_table$log2FC > 0,]
neg_all_sc_das_pseudobulk_corrected <- all_sc_das_pseudobulk_corrected_table[all_sc_das_pseudobulk_corrected_table$log2FC < 0,]

cell_type_order <- c("pDC", "CD4 Memory", "cDC", "Proliferating", "CD4 Naive", "CD16 Mono", "MAIT", "CD8 Naive", "B", "CD8 Memory", "NK", "CD14 Mono")
all_sc_das_pseudobulk_corrected_table$Cell_Type <- factor(all_sc_das_pseudobulk_corrected_table$Cell_Type, levels = cell_type_order)
all_sc_das_pseudobulk_corrected_table$direction <- factor(all_sc_das_pseudobulk_corrected_table$direction, levels = c("Upregulated", "Downregulated"))

all_sc_das_plot <- ggplot(all_sc_das_pseudobulk_corrected_table, aes(Cell_Type, log2FC, color = direction)) +
  geom_jitter() + geom_hline(yintercept = 0.1, linetype = "dashed", color = "black") + 
  geom_hline(yintercept = -0.1, linetype = "dashed", color = "black") + theme_minimal(base_size = 14) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("Cell Type") +
  ylab("log2FC (Site Expression)") + ggtitle("scATAC DAS by Cell Type (28 Days Post-Exposure vs Pre-Exposure)") +
  theme(plot.title = element_text(hjust = 0.5)) + guides(color=guide_legend(title="Fold Change Direction"))

ggsave(filename = paste0("C:/Users/willi/Desktop/", "scatac_das_jitterplot.tiff"), plot = all_sc_das_plot, device='tiff', dpi=300, width = 15)
