all_sc_degs_table <- data.frame(Cell_Type = character(), Gene = character(), log2FC = numeric())

for(cell_type in unique(sc_pseudobulk_deg_table$Cell_Type)) {
  cell_type_for_file_name <- sub(" ", "_", cell_type)
  cell_type_sc_pseudobulk_deg_table <- sc_pseudobulk_deg_table[sc_pseudobulk_deg_table$Cell_Type == cell_type,]
  write.table(cell_type_sc_pseudobulk_deg_table, 
              file = paste0("C:/Users/willi/Desktop/sc_DEGs/", cell_type_for_file_name, "_pseudobulk_corrected_DEGs.tsv"), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  full_sc_degs <- read.table(paste0(sc_deg_dir, "D28-vs-D_minus_1-degs-", cell_type_for_file_name, "-time_point-controlling_for_subject_id_sc.tsv"),
                                    sep = "\t", header = TRUE)
  write.table(full_sc_degs, 
              file = paste0("C:/Users/willi/Desktop/sc_DEGs/", cell_type_for_file_name, "_full_sc_DEGs.tsv"), 
              sep = "\t", quote = FALSE)
  all_sc_degs_table_for_cell_type <- data.frame(Cell_Type = cell_type, Gene = rownames(full_sc_degs), log2FC = full_sc_degs$avg_log2FC)
  all_sc_degs_table <- rbind(all_sc_degs_table, all_sc_degs_table_for_cell_type)
}

all_sc_degs_table <- subset(all_sc_degs_table, !(Cell_Type == "HSPC"))
all_sc_degs_table <- subset(all_sc_degs_table, !(Cell_Type == "NK_CD56bright"))

all_sc_degs_table <- all_sc_degs_table %>%
  mutate(direction = ifelse(log2FC > 0, "Upregulated", "Downregulated")) %>%
  mutate(Cell_Type = ifelse(Cell_Type == "B naive", "B Naive", Cell_Type)) %>%
  mutate(Cell_Type = ifelse(Cell_Type == "B memory", "B Memory", Cell_Type))

pos_all_sc_degs <- all_sc_degs_table[all_sc_degs_table$log2FC > 0,]
neg_all_sc_degs <- all_sc_degs_table[all_sc_degs_table$log2FC < 0,]



cell_type_order <- c("CD4 Naive", "pDC", "CD8 Naive", "cDC", "MAIT", "CD16 Mono", "B Memory", "B Naive", "CD4 Memory", "CD14 Mono", "NK", "CD8 Memory")
all_sc_degs_table$Cell_Type <- factor(all_sc_degs_table$Cell_Type, levels = cell_type_order)
all_sc_degs_table$direction <- factor(all_sc_degs_table$direction, levels = c("Upregulated", "Downregulated"))

all_sc_degs_plot <- ggplot(all_sc_degs_table, aes(Cell_Type, log2FC, color = direction)) +
  geom_jitter() + geom_hline(yintercept = 0.1, linetype = "dashed", color = "black") + 
  geom_hline(yintercept = -0.1, linetype = "dashed", color = "black") + theme_minimal(base_size = 14) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("Cell Type") +
  ylab("log2FC (Gene Expression)") + ggtitle("scRNA DEGs by Cell Type (28 Days Post-Exposure vs Pre-Exposure)") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = paste0("C:/Users/willi/Desktop/", "scrna_deg_jitterplot.tiff"), plot = all_sc_degs_plot, device='tiff', dpi=300)


