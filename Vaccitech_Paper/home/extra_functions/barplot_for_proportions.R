# RNA
test <- read.table("C:/Users/willi/Desktop/RNA_cell_type_proportion_time_point.csv", sep = ",", header = TRUE)

#manova_results <- manova(cbind(CD16_Mono, T_Naive, CD14_Mono, CD4_Memory, Platelet, B, cDC, NK_MAGICAL, MAIT, CD8_Memory, pDC, Proliferating, HSPC) ~ Condition, data = test)
#one.way <- aov(CD16_Mono ~ Condition, data = test)
#one.way <- aov(CD14_Mono ~ Condition, data = test)
#one.way <- aov(NK_MAGICAL ~ Condition, data = test)

# Only CD14 Mono has any significance
kruskal.test(CD14_Mono ~ Condition, data = test)

# Get ready for stacked barplot
results <- aggregate(cbind(CD16_Mono, T_Naive, CD14_Mono, CD4_Memory, Platelet, B, cDC, NK_MAGICAL, MAIT, CD8_Memory, pDC, Proliferating, HSPC) ~ Condition, data = test, sum)
for(colname in colnames(results)) {
  if(colname != "Condition") {
    results[[colname]] <- results[[colname]] / 4
  }
}

colnames(results) <- c("Condition", "CD16 Monocyte", "T Naive", "CD14 Monocyte", "CD4 Memory", "Platelet", "B", "cDC", "NK", "MAIT", "CD8 Memory", "pDC", "Proliferating", "HSPC")
results$Dendritic <- results$pDC + results$cDC
results <- results[,-c(8,12)]

results_long <- tidyr::gather(results, key = "Cell_Type", value = "Expression", -Condition)

results_long <- results_long %>% 
  mutate(Condition = case_when(
    Condition == 'D_minus_1' ~ 'Day -1',
    Condition == 'D28' ~ 'Day 28',
    # Add more conditions as needed
    TRUE ~ Condition  # Keep other values unchanged
  ))

rna_barplot <- ggplot(results_long, aes(fill=Cell_Type, y=Expression, x=Condition)) + 
  geom_bar(position="stack", stat="identity") + 
  xlab("Group") + ylab("Cell Type Proportion") + scale_fill_brewer(palette = "Paired") +
  labs(fill = "Cell Type") +
  ggtitle("Cell Type Proportions for scRNA-Seq (High Viral Load)") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + coord_fixed(ratio = 1.5) + theme_bw()

ggsave("C:/Users/willi/Desktop/rna_barplot.png", plot = rna_barplot, device = "tiff", dpi = 300)






# ATAC
test <- read.table("C:/Users/willi/Desktop/ATAC_cell_type_proportion_time_point.csv", sep = ",", header = TRUE)

#manova_results <- manova(cbind(T.Naive, CD8.Memory, CD4.Memory, B, NK, CD14.Mono, MAIT, CD16.Mono, Proliferating) ~ Condition, data = test)
#one.way <- aov(CD16_Mono ~ Condition, data = test)
#one.way <- aov(CD14_Mono ~ Condition, data = test)
#one.way <- aov(NK_MAGICAL ~ Condition, data = test)

# Nothing has any significance
kruskal.test(B ~ Condition, data = test)

# Get ready for stacked barplot
results_ATAC <- aggregate(cbind(T.Naive, CD8.Memory, CD4.Memory, B, NK, CD14.Mono, MAIT, CD16.Mono, Proliferating) ~ Condition, data = test, sum)
for(colname in colnames(results_ATAC)) {
  if(colname != "Condition") {
    results_ATAC[[colname]] <- results_ATAC[[colname]] / 4
  }
}

colnames(results_ATAC) <- c("Condition", "T Naive", "CD8 Memory", "CD4 Memory", "B", "NK", "CD14 Mono", "MAIT", "CD16 Mono", "Proliferating")

results_ATAC_long <- tidyr::gather(results_ATAC, key = "Cell_Type", value = "Expression", -Condition)

results_ATAC_long <- results_ATAC_long %>% 
  mutate(Condition = case_when(
    Condition == 'D_minus_1' ~ 'Day -1',
    Condition == 'D28' ~ 'Day 28',
    # Add more conditions as needed
    TRUE ~ Condition  # Keep other values unchanged
  ))

atac_barplot <- ggplot(results_ATAC_long, aes(fill=Cell_Type, y=Expression, x=Condition)) + 
  geom_bar(position="stack", stat="identity") + 
  xlab("Group") + ylab("Cell Type Proportion") + scale_fill_brewer(palette = "Paired") +
  labs(fill = "Cell Type") +
  ggtitle("Cell Type Proportions for scATAC-Seq (High Viral Load)") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + coord_fixed(ratio = 1.5) + theme_bw()

ggsave("C:/Users/willi/Desktop/atac_barplot.png", plot = atac_barplot, device = "tiff", dpi = 300)