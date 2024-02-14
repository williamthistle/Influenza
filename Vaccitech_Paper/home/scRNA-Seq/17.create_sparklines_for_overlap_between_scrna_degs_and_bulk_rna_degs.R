# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# Load RDS files
hvl_upregulated_sc_genes_in_bulk_0.05 <- readRDS(paste0(bulk_results_dir, "hvl_upregulated_sc_genes_found_in_bulk/hvl_upregulated_sc_genes_in_bulk_0.05.RDS"))

# Change Positive / Negative to Significant
hvl_upregulated_sc_genes_in_bulk_0.05$Fold.Change.Direction <- replace(hvl_upregulated_sc_genes_in_bulk_0.05$Fold.Change.Direction,
                                                                       hvl_upregulated_sc_genes_in_bulk_0.05$Fold.Change.Direction == "Positive",
                                                                       "Significant")

hvl_upregulated_sc_genes_in_bulk_0.05$Fold.Change.Direction <- replace(hvl_upregulated_sc_genes_in_bulk_0.05$Fold.Change.Direction,
                                                                       hvl_upregulated_sc_genes_in_bulk_0.05$Fold.Change.Direction == "Negative",
                                                                       "Significant")

# Change absolute fold change to include direction of fold change (and rename column to remover Abs)
hvl_upregulated_sc_genes_in_bulk_0.05 <- hvl_upregulated_sc_genes_in_bulk_0.05 %>%
  mutate(Fold.Change = ifelse(Fold.Change.Direction.Raw == "Negative", -Fold.Change.Abs, Fold.Change.Abs)) %>%
  select(-Fold.Change.Abs)

hvl_upregulated_sc_genes_in_bulk_0.05 <- hvl_upregulated_sc_genes_in_bulk_0.05 %>%
  arrange(Gene)

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

