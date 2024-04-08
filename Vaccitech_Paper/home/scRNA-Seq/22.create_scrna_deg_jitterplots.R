# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# HVL
hvl_jitterplot_data <- scRNA_hvl_placebo_degs

hvl_jitterplot_data <- subset(hvl_jitterplot_data, !(Cell_Type == "HSPC"))
hvl_jitterplot_data <- subset(hvl_jitterplot_data, !(Cell_Type == "NK_CD56bright"))
hvl_jitterplot_data <- subset(hvl_jitterplot_data, !(Cell_Type == "pDC"))
hvl_jitterplot_data <- subset(hvl_jitterplot_data, !(Cell_Type == "Proliferating"))


hvl_jitterplot_data <- hvl_jitterplot_data %>%
  mutate(direction = ifelse(sc_log2FC > 0, "Upregulated", "Downregulated")) %>%
  mutate(Cell_Type = ifelse(Cell_Type == "B naive", "B Naive", Cell_Type)) %>%
  mutate(Cell_Type = ifelse(Cell_Type == "B memory", "B Memory", Cell_Type))

pos_all_sc_degs <- hvl_jitterplot_data[hvl_jitterplot_data$sc_log2FC > 0,]
neg_all_sc_degs <- hvl_jitterplot_data[hvl_jitterplot_data$sc_log2FC < 0,]


cell_type_order <- c("CD16 Mono", "CD8 Naive", "B Memory", "cDC", "MAIT", "B Naive", "CD4 Naive", "NK", "CD14 Mono", "CD4 Memory", "CD8 Memory")
hvl_jitterplot_data$Cell_Type <- factor(hvl_jitterplot_data$Cell_Type, levels = cell_type_order)
hvl_jitterplot_data$direction <- factor(hvl_jitterplot_data$direction, levels = c("Upregulated", "Downregulated"))

all_sc_degs_plot <- ggplot(hvl_jitterplot_data, aes(Cell_Type, sc_log2FC, color = direction)) +
  geom_jitter() + geom_hline(yintercept = 0.1, linetype = "dashed", color = "black") + 
  geom_hline(yintercept = -0.1, linetype = "dashed", color = "black") + theme_minimal(base_size = 14) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("Cell Type") +
  ylab("log2FC (Gene Expression)") + ggtitle("scRNA DEGs by Cell Type (28 Days Post-Exposure vs Pre-Exposure)") +
  theme(plot.title = element_text(hjust = 0.5)) + guides(color=guide_legend(title="Fold Change Direction"))

ggsave(filename = paste0("C:/Users/willi/Desktop/", "scrna_pseudobulk_deg_jitterplot.tiff"), plot = all_sc_degs_plot, device='tiff', dpi=300)

# LVL
lvl_sc_degs_table <- data.frame(Cell_Type = character(), Gene = character(), log2FC = numeric())

for(cell_type in unique(sc_pseudobulk_deg_table$Cell_Type)) {
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
