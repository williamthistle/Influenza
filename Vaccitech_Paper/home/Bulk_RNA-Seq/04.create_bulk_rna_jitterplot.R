# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

bulk_degs_jitter_table <- data.frame(Day = character(), Gene = character(), log2FC = numeric())

all_degs_table_for_d2 <- data.frame(Day = "Day 2", Gene = rownames(raw_high_placebo_period_2_D2_vs_D_minus_1_results), log2FC = raw_high_placebo_period_2_D2_vs_D_minus_1_results$log2FoldChange)
bulk_degs_jitter_table <- rbind(bulk_degs_jitter_table, all_degs_table_for_d2)

all_degs_table_for_d5 <- data.frame(Day = "Day 5", Gene = rownames(raw_high_placebo_period_2_D5_vs_D_minus_1_results), log2FC = raw_high_placebo_period_2_D5_vs_D_minus_1_results$log2FoldChange)
bulk_degs_jitter_table <- rbind(bulk_degs_jitter_table, all_degs_table_for_d5)

all_degs_table_for_d8 <- data.frame(Day = "Day 8", Gene = rownames(raw_high_placebo_period_2_D8_vs_D_minus_1_results), log2FC = raw_high_placebo_period_2_D8_vs_D_minus_1_results$log2FoldChange)
bulk_degs_jitter_table <- rbind(bulk_degs_jitter_table, all_degs_table_for_d8)

all_degs_table_for_d28 <- data.frame(Day = "Day 28", Gene = rownames(raw_high_placebo_period_2_D28_vs_D_minus_1_results), log2FC = raw_high_placebo_period_2_D28_vs_D_minus_1_results$log2FoldChange)
bulk_degs_jitter_table <- rbind(bulk_degs_jitter_table, all_degs_table_for_d28)

bulk_degs_jitter_table <- bulk_degs_jitter_table %>%
  mutate(direction = ifelse(log2FC > 0, "Upregulated", "Downregulated"))

pos_bulk_jitter_degs <- bulk_degs_jitter_table[bulk_degs_jitter_table$log2FC > 0,]
neg_bulk_jitter_degs <- bulk_degs_jitter_table[bulk_degs_jitter_table$log2FC < 0,]

bulk_degs_jitter_table$Day <- factor(bulk_degs_jitter_table$Day, levels = c("Day 2", "Day 5", "Day 8", "Day 28"))
bulk_degs_jitter_table$direction <- factor(bulk_degs_jitter_table$direction, levels = c("Upregulated", "Downregulated"))

bulk_degs_jitter_plot <- ggplot(bulk_degs_jitter_table, aes(Day, log2FC, color = direction)) +
  geom_jitter() + theme_minimal(base_size = 14) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("Day (Post-Exposure)") +
  ylab("log2FC (Gene Expression)") + ggtitle("RNA-Seq DEGs by Day Post-Exposure") +
  theme(plot.title = element_text(hjust = 0.5)) + guides(color=guide_legend(title="Fold Change Direction"))

ggsave(filename = paste0("C:/Users/willi/Desktop/", "bulk_deg_jitterplot.tiff"), plot = bulk_degs_jitter_plot, device='tiff', dpi=300)