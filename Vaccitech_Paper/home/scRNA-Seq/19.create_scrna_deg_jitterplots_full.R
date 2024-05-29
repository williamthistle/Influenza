# HVL
hvl_scrna_jitterplot_table <- data.frame(Cell_Type = character(), Gene = character(), log2FC = numeric())

innate_cell_types_jitterplot <- c("CD14 Mono", "CD16 Mono", "cDC", "NK", "NK_CD56bright", "pDC")

for(cell_type in innate_cell_types_jitterplot) {
  #cell_type_for_file_name <- sub(" ", "_", cell_type)
  cell_type_scRNA_hvl_placebo_degs <- scRNA_hvl_placebo_degs[scRNA_hvl_placebo_degs$Cell_Type == cell_type,]
  #write.table(cell_type_scRNA_hvl_placebo_degs, 
  #            file = paste0("C:/Users/wat2/Desktop/sc_DEGs/", cell_type_for_file_name, "_pseudobulk_corrected_DEGs.tsv"), 
  #            sep = "\t", quote = FALSE, row.names = FALSE)
  #full_sc_degs <- read.table(paste0(sc_deg_dir, "D28-vs-D_minus_1-degs-", cell_type_for_file_name, "-time_point-controlling_for_subject_id_sc.tsv"),
  #                                  sep = "\t", header = TRUE)
  #write.table(full_sc_degs, 
  #            file = paste0("C:/Users/wat2/Desktop/sc_DEGs/", cell_type_for_file_name, "_full_sc_DEGs.tsv"), 
  #            sep = "\t", quote = FALSE)
  hvl_scrna_jitterplot_table_for_cell_type <- data.frame(Cell_Type = cell_type, Gene = cell_type_scRNA_hvl_placebo_degs$Gene_Name, log2FC = cell_type_scRNA_hvl_placebo_degs$sc_log2FC)
  hvl_scrna_jitterplot_table <- rbind(hvl_scrna_jitterplot_table, hvl_scrna_jitterplot_table_for_cell_type)
}

#hvl_scrna_jitterplot_table <- subset(hvl_scrna_jitterplot_table, !(Cell_Type == "HSPC"))
#hvl_scrna_jitterplot_table <- subset(hvl_scrna_jitterplot_table, !(Cell_Type == "NK_CD56bright"))

#hvl_scrna_jitterplot_table <- hvl_scrna_jitterplot_table %>%
  #mutate(direction = ifelse(log2FC > 0, "Upregulated", "Downregulated")) %>%
  #mutate(Cell_Type = ifelse(Cell_Type == "B naive", "B Naive", Cell_Type)) %>%
  #mutate(Cell_Type = ifelse(Cell_Type == "B memory", "B Memory", Cell_Type))

hvl_scrna_jitterplot_table <- hvl_scrna_jitterplot_table %>%
mutate(direction = ifelse(log2FC > 0, "Upregulated", "Downregulated"))


pos_all_sc_degs <- hvl_scrna_jitterplot_table[hvl_scrna_jitterplot_table$log2FC > 0,]
neg_all_sc_degs <- hvl_scrna_jitterplot_table[hvl_scrna_jitterplot_table$log2FC < 0,]



cell_type_order <- c("pDC", "NK_CD56bright", "CD16 Mono", "cDC", "NK", "CD14 Mono")
hvl_scrna_jitterplot_table$Cell_Type <- factor(hvl_scrna_jitterplot_table$Cell_Type, levels = cell_type_order)
hvl_scrna_jitterplot_table$direction <- factor(hvl_scrna_jitterplot_table$direction, levels = c("Upregulated", "Downregulated"))

all_sc_degs_plot <- ggplot(hvl_scrna_jitterplot_table, aes(Cell_Type, log2FC, color = direction)) +
  geom_jitter() + geom_hline(yintercept = 0.1, linetype = "dashed", color = "black") + 
  geom_hline(yintercept = -0.1, linetype = "dashed", color = "black") + theme_minimal(base_size = 14) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("Cell Type") +
  ylab("log2FC (Gene Expression)") + ggtitle("scRNA DEGs by Cell Type (28 Days Post-Exposure vs Pre-Exposure)") +
  theme(plot.title = element_text(hjust = 0.5)) + guides(color=guide_legend(title="Fold Change Direction"))

ggsave(filename = paste0("C:/Users/willi/Desktop/", "scrna_deg_jitterplot.tiff"), plot = all_sc_degs_plot, device='tiff', dpi=300)

# LVL
lvl_sc_degs_table <- data.frame(Cell_Type = character(), Gene = character(), log2FC = numeric())

for(cell_type in unique(scRNA_hvl_placebo_degs$Cell_Type)) {
  cell_type_for_file_name <- sub(" ", "_", cell_type)
  full_sc_degs <- read.table(paste0(sc_deg_lvl_dir, "D28-vs-D_minus_1-degs-", cell_type_for_file_name, "-time_point-controlling_for_subject_id_sc.tsv"),
                             sep = "\t", header = TRUE)
  lvl_sc_degs_table_for_cell_type <- data.frame(Cell_Type = cell_type, Gene = rownames(full_sc_degs), log2FC = full_sc_degs$avg_log2FC)
  lvl_sc_degs_table <- rbind(lvl_sc_degs_table, lvl_sc_degs_table_for_cell_type)
}

lvl_sc_degs_table <- subset(lvl_sc_degs_table, !(Cell_Type == "HSPC"))
lvl_sc_degs_table <- subset(lvl_sc_degs_table, !(Cell_Type == "NK_CD56bright"))

lvl_sc_degs_table <- lvl_sc_degs_table %>%
  mutate(direction = ifelse(log2FC > 0, "Upregulated", "Downregulated")) %>%
  mutate(Cell_Type = ifelse(Cell_Type == "B naive", "B Naive", Cell_Type)) %>%
  mutate(Cell_Type = ifelse(Cell_Type == "B memory", "B Memory", Cell_Type))

pos_all_sc_degs <- lvl_sc_degs_table[lvl_sc_degs_table$log2FC > 0,]
neg_all_sc_degs <- lvl_sc_degs_table[lvl_sc_degs_table$log2FC < 0,]



cell_type_order <- c("CD4 Naive", "pDC", "CD8 Naive", "cDC", "MAIT", "CD16 Mono", "B Memory", "B Naive", "CD4 Memory", "CD14 Mono", "NK", "CD8 Memory")
lvl_sc_degs_table$Cell_Type <- factor(lvl_sc_degs_table$Cell_Type, levels = cell_type_order)
lvl_sc_degs_table$direction <- factor(lvl_sc_degs_table$direction, levels = c("Upregulated", "Downregulated"))

all_sc_degs_plot <- ggplot(lvl_sc_degs_table, aes(Cell_Type, log2FC, color = direction)) +
  geom_jitter() + geom_hline(yintercept = 0.1, linetype = "dashed", color = "black") + 
  geom_hline(yintercept = -0.1, linetype = "dashed", color = "black") + theme_minimal(base_size = 14) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("Cell Type") +
  ylab("log2FC (Gene Expression)") + ggtitle("scRNA DEGs by Cell Type (28 Days Post-Exposure vs Pre-Exposure)") +
  theme(plot.title = element_text(hjust = 0.5)) + guides(color=guide_legend(title="Fold Change Direction"))

ggsave(filename = paste0("C:/Users/willi/Desktop/", "scrna_deg_jitterplot.tiff"), plot = all_sc_degs_plot, device='tiff', dpi=300)
