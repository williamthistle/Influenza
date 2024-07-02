all_mintchip_das_table <- data.frame(Marker = character(), Site = character(), log2FC = numeric())

for(marker in mintchip_markers) {
  full_mintchip_das <- read.table(paste0(mintchip_das_dir, marker, "/", marker, "_DESeq2_FC_0.1.tsv"),
                                    sep = "\t", header = TRUE)
  full_mintchip_das_for_cell_type <- data.frame(Marker = marker, Site = full_mintchip_das$coordinates, log2FC = full_mintchip_das$Fold)
  all_mintchip_das_table <- rbind(all_mintchip_das_table, full_mintchip_das_for_cell_type)
}

all_mintchip_das_table <- all_mintchip_das_table %>%
  mutate(direction = ifelse(log2FC > 0, "Upregulated", "Downregulated"))
  
pos_all_mintchip_das <- all_mintchip_das_table[all_mintchip_das_table$log2FC > 0,]
neg_all_mintchip_das <- all_mintchip_das_table[all_mintchip_das_table$log2FC < 0,]

marker_order <- c("H3K27me3", "H3K27Ac", "H3K4me3", "H3K36me3", "H3K9me3", "H3K4me1")
all_mintchip_das_table$Marker <- factor(all_mintchip_das_table$Marker, levels = marker_order)
all_mintchip_das_table$direction <- factor(all_mintchip_das_table$direction, levels = c("Upregulated", "Downregulated"))

mintchip_das_plot <- ggplot(all_mintchip_das_table, aes(Marker, log2FC, color = direction)) +
  geom_jitter() + geom_hline(yintercept = 0.1, linetype = "dashed", color = "black") + 
  geom_hline(yintercept = -0.1, linetype = "dashed", color = "black") + theme_minimal(base_size = 18) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("Histone Marker") +
  ylab("log2FC (Marker Enrichment)") + theme(plot.title = element_text(hjust = 0.5)) + 
  guides(color=guide_legend(title="Fold Change Direction")) + theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1, size = 18))

ggsave(filename = paste0("C:/Users/willi/Desktop/", "mintchip_das_jitterplot.png"), plot = mintchip_das_plot, device='png', dpi=300, width = 4, height = 5, units = "in")


