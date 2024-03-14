# RNA
scRNA_cell_type_proportions <- read.table("C:/Users/willi/Desktop/RNA_cell_type_proportion_time_point.csv", sep = ",", header = TRUE)
scRNA_cell_type_proportions$Condition <- factor(scRNA_cell_type_proportions$Condition, levels = c("D_minus_1", "D28"))

cell_type_proportion_p_values <- c()
cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(CD16_Mono ~ Condition, data = scRNA_cell_type_proportions)))
cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(T_Naive ~ Condition, data = scRNA_cell_type_proportions)))
cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(CD14_Mono ~ Condition, data = scRNA_cell_type_proportions)))
cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(CD4_Memory ~ Condition, data = scRNA_cell_type_proportions)))
cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(B ~ Condition, data = scRNA_cell_type_proportions)))
cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(cDC ~ Condition, data = scRNA_cell_type_proportions)))
cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(NK_MAGICAL ~ Condition, data = scRNA_cell_type_proportions)))
cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(MAIT ~ Condition, data = scRNA_cell_type_proportions)))
cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(CD8_Memory ~ Condition, data = scRNA_cell_type_proportions)))
cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(pDC ~ Condition, data = scRNA_cell_type_proportions)))
cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(Proliferating ~ Condition, data = scRNA_cell_type_proportions)))
cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(HSPC ~ Condition, data = scRNA_cell_type_proportions)))
# Adjust for multiple hypothesis testing
cell_type_proportion_p_values_adjusted <- p.adjust(cell_type_proportion_p_values, method = "BH")
# names(cell_type_proportion_p_values_adjusted) <- bulk_cell_types







# Get ready for stacked barplot
results <- aggregate(cbind(CD16_Mono, T_Naive, CD14_Mono, CD4_Memory, Platelet, B, cDC, NK_MAGICAL, MAIT, CD8_Memory, pDC, Proliferating, HSPC) ~ Condition, data = test, sum)
for(colname in colnames(results)) {
  if(colname != "Condition") {
    results[[colname]] <- results[[colname]] / 6
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
  labs(fill = "Cell Type") + theme_bw() +
  ggtitle("Cell Type Proportions for scRNA-Seq") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + coord_fixed(ratio = 1.5)

ggsave("C:/Users/wat2/Desktop/rna_barplot.tiff", plot = rna_barplot, device = "tiff", dpi = 300)






# ATAC
test <- read.table("C:/Users/wat2/Desktop/ATAC_cell_type_proportion_time_point.csv", sep = ",", header = TRUE)

#manova_results <- manova(cbind(T.Naive, CD8.Memory, CD4.Memory, B, NK, CD14.Mono, MAIT, CD16.Mono, Proliferating) ~ Condition, data = test)
#one.way <- aov(CD16_Mono ~ Condition, data = test)
#one.way <- aov(CD14_Mono ~ Condition, data = test)
#one.way <- aov(NK_MAGICAL ~ Condition, data = test)

# Nothing has any significance
kruskal.test(CD14.Mono ~ Condition, data = test)

# Get ready for stacked barplot
results_ATAC <- aggregate(cbind(T.Naive, CD8.Memory, CD4.Memory, B, NK, CD14.Mono, MAIT, CD16.Mono, Proliferating) ~ Condition, data = test, sum)
for(colname in colnames(results_ATAC)) {
  if(colname != "Condition") {
    results_ATAC[[colname]] <- results_ATAC[[colname]] / 6
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
  labs(fill = "Cell Type") + theme_bw() +
  ggtitle("Cell Type Proportions for scATAC-Seq") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + coord_fixed(ratio = 1.5) 

ggsave("C:/Users/wat2/Desktop/atac_barplot.tiff", plot = atac_barplot, device = "tiff", dpi = 300)