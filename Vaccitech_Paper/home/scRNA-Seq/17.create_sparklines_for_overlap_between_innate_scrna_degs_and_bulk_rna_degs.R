# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# Load RDS files
hvl_upregulated_sc_genes_in_bulk_0.05 <- readRDS(paste0(bulk_results_dir, "hvl_upregulated_sc_genes_found_in_bulk/hvl_upregulated_sc_genes_in_bulk_0.05.RDS"))

# Remove cell types
hvl_upregulated_sc_genes_in_bulk_0.05 <- hvl_upregulated_sc_genes_in_bulk_0.05 %>%
  mutate(Gene = sapply(Gene, remove_text_in_parentheses))

# Change Positive / Negative to Significant
hvl_upregulated_sc_genes_in_bulk_0.05$Fold.Change.Direction <- replace(hvl_upregulated_sc_genes_in_bulk_0.05$Fold.Change.Direction,
                                                                       hvl_upregulated_sc_genes_in_bulk_0.05$Fold.Change.Direction == "Positive",
                                                                       "Significant")

hvl_upregulated_sc_genes_in_bulk_0.05$Fold.Change.Direction <- replace(hvl_upregulated_sc_genes_in_bulk_0.05$Fold.Change.Direction,
                                                                       hvl_upregulated_sc_genes_in_bulk_0.05$Fold.Change.Direction == "Negative",
                                                                       "Significant")

# Change absolute fold change to include direction of fold change (and rename column to remove Abs)
hvl_upregulated_sc_genes_in_bulk_0.05 <- hvl_upregulated_sc_genes_in_bulk_0.05 %>%
  mutate(Fold.Change = ifelse(Fold.Change.Direction.Raw == "Negative", -Fold.Change.Abs, Fold.Change.Abs)) %>%
  select(-Fold.Change.Abs)

# Unique genes
unique_genes <- unique(hvl_upregulated_sc_genes_in_bulk_0.05$Gene)

# Create a new data frame with additional rows
new_rows <- data.frame(
  Gene = rep(unique_genes, each = 1),  # Repeat each gene once
  Day = "Day.Minus.1",
  Fold.Change.Direction.Raw = "Positive",
  Fold.Change.Direction = "Not Significant",
  Adjusted.P.Value = 1,
  Fold.Change = 0
)

# Concatenate the original data frame with the new rows
hvl_upregulated_sc_genes_in_bulk_0.05 <- rbind(hvl_upregulated_sc_genes_in_bulk_0.05, new_rows)

hvl_upregulated_sc_genes_in_bulk_0.05 <- hvl_upregulated_sc_genes_in_bulk_0.05 %>%
  arrange(Gene)

hvl_upregulated_sc_genes_in_bulk_0.05 <- hvl_upregulated_sc_genes_in_bulk_0.05[order(hvl_upregulated_sc_genes_in_bulk_0.05$Gene, match(hvl_upregulated_sc_genes_in_bulk_0.05$Day, c("Day.Minus.1", "Day.2", "Day.5", "Day.8", "Day.28"))), ]

for(current_gene in unique(hvl_upregulated_sc_genes_in_bulk_0.05$Gene)) {
  current_gene_stats <- hvl_upregulated_sc_genes_in_bulk_0.05[hvl_upregulated_sc_genes_in_bulk_0.05$Gene == current_gene,]
  current_gene_fc <- c()
  value_spots_list <- list()
  for(current_row_index in 1:nrow(current_gene_stats)) {
    current_row <- current_gene_stats[current_row_index,]
    current_gene_fc <- c(current_gene_fc, round(current_row$Fold.Change, digits = 7))
    if(current_row$Fold.Change.Direction == "Significant") {
      value_spots_list[[as.character(round(current_row$Fold.Change, digits = 7))]] <- "red"
    }
  }
  # Save as 1100 x 230 in RStudio window
  sparklines::sparkline(current_gene_fc, "line", list(fillColor="white", spotRadius = 16, lineWidth = 8, minSpotColor = "", maxSpotColor = "", spotColor = "", valueSpots = value_spots_list), 
                 width = "1000px", height = "200px")
  #saveWidget(sparkline_plot, file = paste0("C:/Users/willi/Desktop/sparkline_test/", current_gene, ".html"))
  #webshot(paste0("C:/Users/willi/Desktop/sparkline_test/", current_gene, ".html"), 
  #        file = paste0("C:/Users/willi/Desktop/sparkline_test/", current_gene, ".png"), vwidth = 1000, vheight = 200)
}


gene_sequences <- c()

# Iterate over unique genes
for (gene in unique(hvl_upregulated_sc_genes_in_bulk_0.05$Gene)) {
  # Subset the dataframe for the current gene
  gene_data <- hvl_upregulated_sc_genes_in_bulk_0.05[hvl_upregulated_sc_genes_in_bulk_0.05$Gene == gene, ]
  
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
result <- data.frame(Gene = unique(hvl_upregulated_sc_genes_in_bulk_0.05$Gene), Fold.Change.Sequence = gene_sequences)

result <- result %>% 
  arrange(Fold.Change.Sequence)

print(result)

hvl_upregulated_sc_genes_in_bulk_0.05$Day <- factor(hvl_upregulated_sc_genes_in_bulk_0.05$Day, levels = c("Day.Minus.1",
                                                                                                          "Day.2", "Day.5",
                                                                                                          "Day.8", "Day.28"))

# Line plot for all genes
ggplot(data=hvl_upregulated_sc_genes_in_bulk_0.05, aes(x=Day, y=Fold.Change, group = Gene, color = Gene)) +
  geom_line() +
  geom_point()

# Line plot for up up down down genes
up_up_down_down_genes <- result[result$Fold.Change.Sequence == "UP UP DOWN DOWN",]$Gene
hvl_upregulated_sc_genes_in_bulk_0.05_up_up_down_down_subset <- hvl_upregulated_sc_genes_in_bulk_0.05[hvl_upregulated_sc_genes_in_bulk_0.05$Gene %in% up_up_down_down_genes,]

ggplot(data=hvl_upregulated_sc_genes_in_bulk_0.05_up_up_down_down_subset, aes(x=Day, y=Fold.Change, group = Gene, color = Gene)) +
  geom_line() +
  geom_point()

# Line plot for up down up up genes
up_down_up_up_genes <- result[result$Fold.Change.Sequence == "UP DOWN UP UP",]$Gene
hvl_upregulated_sc_genes_in_bulk_0.05_up_down_up_up_subset <- hvl_upregulated_sc_genes_in_bulk_0.05[hvl_upregulated_sc_genes_in_bulk_0.05$Gene %in% up_down_up_up_genes,]

ggplot(data=hvl_upregulated_sc_genes_in_bulk_0.05_up_down_up_up_subset, aes(x=Day, y=Fold.Change, group = Gene, color = Gene)) +
  geom_line() +
  geom_point()

# Line plot for up down down up genes
up_down_down_up_genes <- result[result$Fold.Change.Sequence == "UP DOWN DOWN UP",]$Gene
hvl_upregulated_sc_genes_in_bulk_0.05_up_down_down_up_subset <- hvl_upregulated_sc_genes_in_bulk_0.05[hvl_upregulated_sc_genes_in_bulk_0.05$Gene %in% up_down_down_up_genes,]

ggplot(data=hvl_upregulated_sc_genes_in_bulk_0.05_up_down_down_up_subset, aes(x=Day, y=Fold.Change, group = Gene, color = Gene)) +
  geom_line() +
  geom_point()

# Line plot for up down up down genes
up_down_up_down_genes <- result[result$Fold.Change.Sequence == "UP DOWN UP DOWN",]$Gene
hvl_upregulated_sc_genes_in_bulk_0.05_up_down_up_down_subset <- hvl_upregulated_sc_genes_in_bulk_0.05[hvl_upregulated_sc_genes_in_bulk_0.05$Gene %in% up_down_up_down_genes,]

ggplot(data=hvl_upregulated_sc_genes_in_bulk_0.05_up_down_up_down_subset, aes(x=Day, y=Fold.Change, group = Gene, color = Gene)) +
  geom_line() +
  geom_point()

# Line plot for up up down up genes
up_up_down_up_genes <- result[result$Fold.Change.Sequence == "UP UP DOWN UP",]$Gene
hvl_upregulated_sc_genes_in_bulk_0.05_up_up_down_up_subset <- hvl_upregulated_sc_genes_in_bulk_0.05[hvl_upregulated_sc_genes_in_bulk_0.05$Gene %in% up_up_down_up_genes,]

ggplot(data=hvl_upregulated_sc_genes_in_bulk_0.05_up_up_down_up_subset, aes(x=Day, y=Fold.Change, group = Gene, color = Gene)) +
  geom_line() +
  geom_point()

