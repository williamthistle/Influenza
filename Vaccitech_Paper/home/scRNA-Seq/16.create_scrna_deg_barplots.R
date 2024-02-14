# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

sc_cell_types <- c("CD14 Mono", "CD16 Mono", "cDC", "pDC", "NK", "NK_CD56bright", "CD4 Naive", "CD8 Naive", "CD4 Memory", "CD8 Memory",
                   "MAIT", "B memory", "B naive")
sc_cell_types <- rep(sc_cell_types, each = 4)
fc_dir_vector <- rep(c("Upregulated", "Downregulated"), 26)
pseudobulk_correction_vector <- rep(c("Pseudobulk Corrected", "Pseudobulk Corrected", "Raw", "Raw"), 13)
count_vector <- c()

for(cell_type in unique(sc_cell_types)) {
  cell_type_specific_sc_pseudobulk_deg_table <- sc_pseudobulk_deg_table[sc_pseudobulk_deg_table$Cell_Type == cell_type,]
  pos_pseudobulk_genes <- nrow(cell_type_specific_sc_pseudobulk_deg_table[cell_type_specific_sc_pseudobulk_deg_table$sc_log2FC > 0,])
  neg_pseudobulk_genes <- nrow(cell_type_specific_sc_pseudobulk_deg_table[cell_type_specific_sc_pseudobulk_deg_table$sc_log2FC < 0,])
  cell_type_for_file_name <- sub(" ", "_", cell_type)
  current_deg_table <- read.table(paste0(sc_deg_dir, "D28-vs-D_minus_1-degs-", cell_type_for_file_name, "-time_point-controlling_for_subject_id_sc.tsv"), sep = "\t", header = TRUE)
  pos_genes <- nrow(current_deg_table[current_deg_table$avg_log2FC > 0,])
  neg_genes <- nrow(current_deg_table[current_deg_table$avg_log2FC < 0,])
  current_count_vector <- c(pos_pseudobulk_genes, neg_pseudobulk_genes, pos_genes, neg_genes)
  count_vector <- c(count_vector, current_count_vector)
}

scrna_deg_plot_df <- data.frame(cell_types = sc_cell_types, pseudobulk_correction = pseudobulk_correction_vector, 
                                fc_direction = fc_dir_vector, count = count_vector)

scrna_deg_plot_df$cell_types <- replace(scrna_deg_plot_df$cell_types, scrna_deg_plot_df$cell_types=="B naive", "B Naive")
scrna_deg_plot_df$cell_types <- replace(scrna_deg_plot_df$cell_types, scrna_deg_plot_df$cell_types=="B memory", "B Memory")

scrna_deg_plot_df$cell_types <- factor(scrna_deg_plot_df$cell_types, levels = unique(scrna_deg_plot_df$cell_types))
scrna_deg_plot_df$fc_direction <- factor(scrna_deg_plot_df$fc_direction, levels = c("Upregulated", "Downregulated"))
scrna_deg_plot_df$interaction <- interaction(scrna_deg_plot_df$pseudobulk_correction, scrna_deg_plot_df$fc_direction)
scrna_deg_plot_df$interaction <- factor(scrna_deg_plot_df$interaction, c("Raw.Upregulated", "Pseudobulk Corrected.Upregulated",
                                                                         "Raw.Downregulated", 
                                                                         "Pseudobulk Corrected.Downregulated"))


ggplot(scrna_deg_plot_df, aes(x = cell_types, y = count, fill = interaction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  labs(title = "scRNA-Seq DEGs (28 Days Post-Exposure vs. Pre-Exposure to H3N2)",
       x = "Cell Types",
       y = "Count",
       fill = "Fold Change Direction\n") +
  scale_fill_manual(values = c("Pseudobulk Corrected.Upregulated" = "maroon", 
                               "Raw.Upregulated" = "maroon4",
                               "Pseudobulk Corrected.Downregulated" = "steelblue",
                               "Raw.Downregulated" = "steelblue4"),
                    labels=c("Upregulated (Raw)","Upregulated (Pseudobulk Corrected)","Downregulated (Raw)","Downregulated (Pseudobulk Corrected)")) +
  theme_minimal(base_size = 14) + theme(plot.title = element_text(hjust = 0.5))

raw_scrna_deg_plot_df <- scrna_deg_plot_df
raw_scrna_deg_plot_df <- raw_scrna_deg_plot_df[raw_scrna_deg_plot_df$pseudobulk_correction == "Raw",]

ggplot(raw_scrna_deg_plot_df, aes(x = cell_types, y = count, fill = fc_direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Upregulated" = "maroon", "Downregulated" = "steelblue")) +
  labs(title = "scRNA-Seq DEGs (28 Days Post-Exposure vs. Pre-Exposure to H3N2)",
       x = "Cell Type",
       y = "Count") +
  theme_minimal(base_size = 14) + guides(fill=guide_legend(title="Fold Change Direction")) + 
  theme(plot.title = element_text(hjust = 0.5))

