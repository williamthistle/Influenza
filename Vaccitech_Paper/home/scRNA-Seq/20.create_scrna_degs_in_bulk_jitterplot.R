# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

#hvl_upregulated_sc_genes_in_bulk_unfiltered <- fill_in_sc_deg_info_for_time_series(innate_sc_pseudobulk_deg_table, high_placebo_counts, high_placebo_metadata,
#                                                                             paste0(bulk_results_dir, "hvl_upregulated_sc_genes_found_in_bulk_unfiltered/"), "up", alpha = 0.05,
#                                                                             filter_D28 = FALSE)
#saveRDS(hvl_upregulated_sc_genes_in_bulk_unfiltered, file = paste0(bulk_results_dir, "hvl_upregulated_sc_genes_found_in_bulk/hvl_upregulated_sc_genes_in_bulk_unfiltered.RDS"))

hvl_upregulated_sc_genes_in_bulk_unfiltered <- readRDS(paste0(bulk_results_dir, "hvl_upregulated_sc_genes_found_in_bulk/hvl_upregulated_sc_genes_in_bulk_unfiltered.RDS"))

# Remove cell types
hvl_upregulated_sc_genes_in_bulk_unfiltered <- hvl_upregulated_sc_genes_in_bulk_unfiltered %>%
  mutate(Gene = sapply(Gene, remove_text_in_parentheses))

# Change Positive / Negative to Significant
hvl_upregulated_sc_genes_in_bulk_unfiltered$Fold.Change.Direction <- replace(hvl_upregulated_sc_genes_in_bulk_unfiltered$Fold.Change.Direction,
                                                                       hvl_upregulated_sc_genes_in_bulk_unfiltered$Fold.Change.Direction == "Positive",
                                                                       "Significant")

hvl_upregulated_sc_genes_in_bulk_unfiltered$Fold.Change.Direction <- replace(hvl_upregulated_sc_genes_in_bulk_unfiltered$Fold.Change.Direction,
                                                                       hvl_upregulated_sc_genes_in_bulk_unfiltered$Fold.Change.Direction == "Negative",
                                                                       "Significant")

# Change absolute fold change to include direction of fold change (and rename column to remove Abs)
hvl_upregulated_sc_genes_in_bulk_unfiltered <- hvl_upregulated_sc_genes_in_bulk_unfiltered %>%
  mutate(Fold.Change = ifelse(Fold.Change.Direction.Raw == "Negative", -Fold.Change.Abs, Fold.Change.Abs)) %>%
  select(-Fold.Change.Abs)

hvl_upregulated_sc_genes_in_bulk_unfiltered$Day <- factor(hvl_upregulated_sc_genes_in_bulk_unfiltered$Day, levels = c("Day.2", "Day.5", "Day.8", "Day.28"))
hvl_upregulated_sc_genes_in_bulk_unfiltered$Fold.Change.Direction.Raw <- factor(hvl_upregulated_sc_genes_in_bulk_unfiltered$Fold.Change.Direction.Raw, 
                                                                                levels = c("Positive", "Negative"))
hvl_upregulated_sc_genes_in_bulk_unfiltered$Fold.Change.Direction <- factor(hvl_upregulated_sc_genes_in_bulk_unfiltered$Fold.Change.Direction, 
                                                                                levels = c("Significant", "Not Significant"))

hvl_upregulated_sc_genes_in_bulk_unfiltered_plot <- ggplot(hvl_upregulated_sc_genes_in_bulk_unfiltered, aes(Day, Fold.Change, color = Fold.Change.Direction.Raw, shape = Fold.Change.Direction, label = ifelse(abs(Fold.Change) >= 1 == "Significant",as.character(Gene),''))) +
  geom_jitter(size=3, position = position_jitter(seed = 1)) + ggrepel::geom_text_repel(position = position_jitter(seed = 1)) +
  scale_shape_manual(values = c(8, 16)) + theme_minimal(base_size = 14) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("Day (Post-Exposure)") +
  ylab("log2FC (Gene Expression)") + ggtitle("Fold Change Over Time for Single Cell DEGs Upregulated at Day 28") +
  theme(plot.title = element_text(hjust = 0.5)) + guides(color=guide_legend(title="Fold Change"))

ggsave(filename = paste0("C:/Users/willi/Desktop/", "bulk_deg_jitterplot.tiff"), plot = bulk_degs_jitter_plot, device='tiff', dpi=300)

hvl_upregulated_sc_genes_in_bulk_unfiltered_day_2 <- hvl_upregulated_sc_genes_in_bulk_unfiltered[hvl_upregulated_sc_genes_in_bulk_unfiltered$Day == "Day.2",]
hvl_upregulated_sc_genes_in_bulk_unfiltered_day_2 <- hvl_upregulated_sc_genes_in_bulk_unfiltered_day_2[order(hvl_upregulated_sc_genes_in_bulk_unfiltered_day_2$Fold.Change),]

hvl_upregulated_sc_genes_in_bulk_unfiltered_day_5 <- hvl_upregulated_sc_genes_in_bulk_unfiltered[hvl_upregulated_sc_genes_in_bulk_unfiltered$Day == "Day.5",]
hvl_upregulated_sc_genes_in_bulk_unfiltered_day_5 <- hvl_upregulated_sc_genes_in_bulk_unfiltered_day_5[order(hvl_upregulated_sc_genes_in_bulk_unfiltered_day_5$Fold.Change),]

hvl_upregulated_sc_genes_in_bulk_unfiltered_day_8 <- hvl_upregulated_sc_genes_in_bulk_unfiltered[hvl_upregulated_sc_genes_in_bulk_unfiltered$Day == "Day.8",]
hvl_upregulated_sc_genes_in_bulk_unfiltered_day_8 <- hvl_upregulated_sc_genes_in_bulk_unfiltered_day_8[order(hvl_upregulated_sc_genes_in_bulk_unfiltered_day_8$Fold.Change),]

hvl_upregulated_sc_genes_in_bulk_unfiltered_day_28 <- hvl_upregulated_sc_genes_in_bulk_unfiltered[hvl_upregulated_sc_genes_in_bulk_unfiltered$Day == "Day.28",]

# Highest upregulated day 5 genes
tail(hvl_upregulated_sc_genes_in_bulk_unfiltered_day_2, n = 20)
# Highest downregulated day 5 genes
hvl_upregulated_sc_genes_in_bulk_unfiltered_day_2[1:20,]

# Highest upregulated day 5 genes
tail(hvl_upregulated_sc_genes_in_bulk_unfiltered_day_5, n = 20)
# Highest downregulated day 5 genes
hvl_upregulated_sc_genes_in_bulk_unfiltered_day_5[1:20,]

# Highest upregulated day 8 genes
tail(hvl_upregulated_sc_genes_in_bulk_unfiltered_day_8, n = 20)
# Highest downregulated day 8 genes
hvl_upregulated_sc_genes_in_bulk_unfiltered_day_8[1:20,]

# Highest upregulated day 28 genes
tail(hvl_upregulated_sc_genes_in_bulk_unfiltered_day_28, n = 20)

# Create a new data frame with additional rows
new_rows <- data.frame(
  Gene = rep(unique(hvl_upregulated_sc_genes_in_bulk_unfiltered$Gene), each = 1),  # Repeat each gene once
  Day = "Day.Minus.1",
  Fold.Change.Direction.Raw = "Positive",
  Fold.Change.Direction = "Not Significant",
  Adjusted.P.Value = 1,
  Fold.Change = 0
)

# Concatenate the original data frame with the new rows
hvl_upregulated_sc_genes_in_bulk_unfiltered_with_pre_exposure <- rbind(hvl_upregulated_sc_genes_in_bulk_unfiltered, new_rows)

hvl_upregulated_sc_genes_in_bulk_unfiltered_with_pre_exposure <- hvl_upregulated_sc_genes_in_bulk_unfiltered_with_pre_exposure %>%
  arrange(Gene)

hvl_upregulated_sc_genes_in_bulk_unfiltered_with_pre_exposure <- hvl_upregulated_sc_genes_in_bulk_unfiltered_with_pre_exposure[order(hvl_upregulated_sc_genes_in_bulk_unfiltered_with_pre_exposure$Gene, match(hvl_upregulated_sc_genes_in_bulk_unfiltered_with_pre_exposure$Day, c("Day.Minus.1", "Day.2", "Day.5", "Day.8", "Day.28"))), ]



gene_sequences <- c()

# Iterate over unique genes
for (gene in unique(hvl_upregulated_sc_genes_in_bulk_unfiltered_with_pre_exposure$Gene)) {
  # Subset the dataframe for the current gene
  gene_data <- hvl_upregulated_sc_genes_in_bulk_unfiltered_with_pre_exposure[hvl_upregulated_sc_genes_in_bulk_unfiltered_with_pre_exposure$Gene == gene, ]
  
  # Initialize an empty vector to store the directions
  directions <- c()
  
  # Iterate over rows for the current gene
  for (i in 1:(nrow(gene_data) - 1)) {
    # Check if the fold change is increasing or decreasing
    direction <- ifelse(gene_data$Fold.Change[i + 1] > gene_data$Fold.Change[i], "UP", "DOWN")
    directions <- c(directions, direction)
  }
  
  # Combine the directions and store in the vector
  gene_sequence <- paste(directions, collapse = " ")
  gene_sequences <- c(gene_sequences, gene_sequence)
}

# Create a new dataframe with gene names and sequences
result <- data.frame(Gene = unique(hvl_upregulated_sc_genes_in_bulk_unfiltered_with_pre_exposure$Gene), Fold.Change.Sequence = gene_sequences)

result <- result %>% 
  arrange(Fold.Change.Sequence)

print(result)

