### THIS FUNCTION:
### RUNS MAGICAL AND CREATES 1) OVERALL MAGICAL DF, AND 2) TF-CENTRIC MAGICAL DF

# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

source(paste0(base_dir, "extra_functions/MAGICAL_functions.R"))
source(paste0(base_dir, "extra_functions/my_MAGICAL_helper_functions.R"))

# Grab starting DFs
overall_magical_df <- read.table(paste0(MAGICAL_hvl_placebo_output_dir, "MAGICAL_overall_output_with_stage_1_info.tsv"), sep = "\t", header = TRUE)
magical_gene_overlap_df <- create_magical_gene_overlap_df(overall_magical_df, hvl_full_time_series_placebo_period_2_D5_vs_D_minus_1_results, hvl_full_time_series_placebo_period_2_D8_vs_D_minus_1_results)
hg38_ref_seq <- read.table(paste0(MAGICAL_hvl_placebo_dir, "hg38_Refseq.txt"), sep = "\t", header = TRUE)

# Parse to important genes
important_magical_genes <- c("IL32", "CASP1", "CSF1R", "IRAK3", "IL1RAP", "NFIL3", "PTGES", "IL34", "CSF1R", "IL1B", "IL15RA", "IL21R", "IL1RN", "IL6R", "IL6R", "IL5RA", "IL4R", "IL6", "IL37", "IL17F", "IL20RB", "IRAK4", "IL10RA", "IL17B", "IL11", "IRF2", "IRF7", "IFNG", "OAS1", "PSMB9", "MNDA", "DNAJC3", "IFNGR1", "GBP5", "USP38", "IRF2", "IRF1", "IFNGR1", "IFI35", "IRF8", "IFITM3", "IRF5", "IFI6", "IRF1", "STING1", "IRF5", "IFI27", "IFNAR2", "LUARIS", "STING1", "ABHD17A", "CSK", "MAP3K11", "DUSP1", "DUSP2", "DUSP6", "MAPK7", "MAPK8", "MAP2K1", "MAP3K8", "MAP3K20", "MAP4K3", "MAPKAPK2", "PTK2B", "RELL1", "MINK1", "BRAF", "RIPK1", "MAPK13", "MAP2K5", "MAP3K11", "DUSP16", "DUSP16", "MAP3K15", "MAPK3", "MAP3K11", "MAP4K4", "MAP3K6", "MAPKBP1", "MAP3K3", "MAPK8IP1", "MINK1", "MAP4K4", "MAP4K3", "RIPK1", "MAP2K3", "MAP4K4", "MAP4K1", "MAP2K6", "MAP2K3", "MAPK15", "MAP2K2", "MAP4K2", "MAP2K1", "MAPKAP1", "MAPK10", "MAP4K1", "FOSL2", "FOSL1", "JUND", "JUN", "JUNB", "JUND", "JDP2", "FOSL2", "FOS", "FOSB", "JAK1", "STAT3", "STAT4", "JAK1", "STAT2", "STAT4", "STAT6", "JAK2", "JAK3", "STAT4", "JAK2", "JAK2", "JAK1", "STAT5A", "CCL3", "CX3CR1", "CCL3L1", "CXCL16", "CXCR2", "CCL22", "CCL26", "CXCR3", "CCL17", "CCR7", "CXCL12", "XCR1", "HIST1H1C", "HIST1H1D", "HIST1H1E", "H2AFZ", "H2BC5", "H2BC18", "H4C13", "H2BC15", "H3C2", "H2AC14", "H2AC4", "H4C8", "H2AC13", "H3-3A", "H2BC10", "H4C5", "H3-3B", "H2AC25", "H1-4", "KAT6A", "BRD1", "MSL2", "KAT2B", "KAT2A", "KAT2B", "BRD1", "MSL2", "KAT5", "HDAC5", "HDAC7", "HDAC9", "HDAC10", "HDAC4", "HDAC4", "HDAC11", "METTL23", "ASH1L", "PRDM2", "SETD2", "ASH1L", "SUV39H1", "ASH2L", "KMT2B", "SETD7", "PRMT2", "KMT2C", "KDM3A", "KDM2B", "KDM2B", "KDM5B", "KDM2A", "KDM7A", "KDM2A", "KDM8", "KDM4A", "KDM4B", "KDM2B", "KDM6A", "KDM4C", "KDM4A", "KDM4B", "KDM4E", "KDM5A")
important_overall_magical_df <- overall_magical_df[overall_magical_df$Gene_symbol %in% important_magical_genes,]
# Remove some singleton genes from categories
important_overall_magical_df <- important_overall_magical_df[!(important_overall_magical_df$Gene_symbol %in% c("CCL3L1", "STAT4")),]
important_magical_gene_overlap_df <- magical_gene_overlap_df[magical_gene_overlap_df$Gene_Name %in% important_magical_genes,]

# Mapping gene names to gene types
heatmap_gene_types <- list()

heatmap_gene_types[["IRAK3"]] <- "Interleukin"
heatmap_gene_types[["IL1RAP"]] <- "Interleukin"
heatmap_gene_types[["CASP1"]] <- "Interleukin"

heatmap_gene_types[["DNAJC3"]] <- "Interferon"
heatmap_gene_types[["GBP5"]] <- "Interferon"
heatmap_gene_types[["IFNG"]] <- "Interferon"
heatmap_gene_types[["PSMB9"]] <- "Interferon"
heatmap_gene_types[["IRF2"]] <- "Interferon"

heatmap_gene_types[["FOSB"]] <- "AP-1"
heatmap_gene_types[["FOSL2"]] <- "AP-1"
heatmap_gene_types[["JDP2"]] <- "AP-1"

heatmap_gene_types[["DUSP1"]] <- "MAP Kinase"
heatmap_gene_types[["DUSP2"]] <- "MAP Kinase"
heatmap_gene_types[["DUSP6"]] <- "MAP Kinase"
heatmap_gene_types[["MAP3K11"]] <- "MAP Kinase"
heatmap_gene_types[["CSF1R"]] <- "MAP Kinase"
heatmap_gene_types[["MAPK7"]] <- "MAP Kinase"
heatmap_gene_types[["MINK1"]] <- "MAP Kinase"
heatmap_gene_types[["MAP4K3"]] <- "MAP Kinase"
heatmap_gene_types[["ABHD17A"]] <- "MAP Kinase"
heatmap_gene_types[["CSK"]] <- "MAP Kinase"
heatmap_gene_types[["RIPK1"]] <- "MAP Kinase"

heatmap_gene_types[["SETD2"]] <- "Lysine Methyltransferase"
heatmap_gene_types[["METTL23"]] <- "Lysine Methyltransferase"
heatmap_gene_types[["ASH1L"]] <- "Lysine Methyltransferase"

heatmap_gene_types[["KDM2B"]] <- "Lysine Demethylase"
heatmap_gene_types[["KDM3A"]] <- "Lysine Demethylase"

heatmap_gene_types[["BRD1"]] <- "Lysine Acetyltransferase"

# Vectors to store info for DF
gene_name <- c()
gene_fc <- c()
cell_type <- c()
site_fc <- c()
binding_tfs <- c()
other_assays <- c()
strand <- c()
dist_to_TSS <- c()

# Add info to DF
for(current_row_index in 1:nrow(important_overall_magical_df)) {
  current_row <- important_overall_magical_df[current_row_index,]
  gene_name <- c(gene_name, current_row$Gene_symbol)
  current_magical_gene_overlap <- important_magical_gene_overlap_df[important_magical_gene_overlap_df$Gene_Name == current_row$Gene_symbol,][1,]
  gene_fc <- c(gene_fc, current_row$rna_sc_fc)
  cell_type <- c(cell_type, current_row$Cell_Type)
  site_fc <- c(site_fc, current_row$atac_sc_fc)
  binding_tfs <- c(binding_tfs, current_row$TFs.binding.prob.)
  gene_strand_info <- unique(hg38_ref_seq[hg38_ref_seq$gene_name == current_row$Gene_symbol,]$strand)
  strand <- c(strand, gene_strand_info)
  current_tss <- current_row$Gene_TSS
  current_peak_midpoint <- floor((current_row$Peak_start + current_row$Peak_end) / 2)
  if(gene_strand_info == "+") {
    dist_to_TSS <- c(dist_to_TSS, current_peak_midpoint - current_tss)
  } else {
    dist_to_TSS <- c(dist_to_TSS, current_tss - current_peak_midpoint)
  }
  current_other_assays <- colnames(current_magical_gene_overlap)[apply(current_magical_gene_overlap, 2, function(col) any(col == TRUE))]
  if(length(current_other_assays) > 0) {
    other_assays <- c(other_assays, paste(current_other_assays, collapse = "/"))
  } else {
    other_assays <- c(other_assays, "N/A")
  }
}

# Create final DF
final_magical_df_for_paper <- data.frame(gene_name = gene_name, gene_fc = gene_fc, strand = strand, cell_type = cell_type,
                                         site_fc = site_fc, binding_tfs = binding_tfs, dist_to_TSS = dist_to_TSS, other_assays = other_assays)

# Add gene types
gene_types <- c()
for(current_row_index in 1:nrow(final_magical_df_for_paper)) {
  current_row <- final_magical_df_for_paper[current_row_index,]
  gene_types <- c(gene_types, heatmap_gene_types[[current_row$gene_name]])
}

final_magical_df_for_paper$gene_type <- gene_types

# Add categories for distance to TSS
dist_to_tss_categories <- c()
for(current_row_index in 1:nrow(final_magical_df_for_paper)) {
  current_row <- final_magical_df_for_paper[current_row_index,]
  current_dist_to_tss <- current_row$dist_to_TSS
  if(abs(current_dist_to_tss) <= 3000) {
    dist_to_tss_categories <- c(dist_to_tss_categories, "Promoter")
  } else if(abs(current_dist_to_tss) <= 10000) {
    dist_to_tss_categories <- c(dist_to_tss_categories, "Distal")
  } else {
    dist_to_tss_categories <- c(dist_to_tss_categories, "Trans")
  }
}

final_magical_df_for_paper$dist_to_tss_category <- dist_to_tss_categories

write.table(final_magical_df_for_paper, file = "C:/Users/wat2/Desktop/final_magical_df_for_paper.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)

magical_heatmap_plots <- list()

for(cell_type in unique(final_magical_df_for_paper$cell_type)) {
  cell_type_final_magical_df_for_paper <- final_magical_df_for_paper[final_magical_df_for_paper$cell_type == cell_type,]
  
  cell_type_final_magical_df_for_paper$gene_type <- factor(cell_type_final_magical_df_for_paper$gene_type, levels = c("Interleukin", "Interferon", "AP-1", "MAP Kinase", "Lysine Methyltransferase", "Lysine Demethylase", "Lysine Acetyltransferase"))
  
  # Organize the rows
  cell_type_final_magical_df_for_paper <- cell_type_final_magical_df_for_paper %>%
    arrange(gene_type, desc(gene_fc), gene_name, desc(site_fc))
  
  cell_type_final_magical_df_for_paper <- cell_type_final_magical_df_for_paper %>%
    group_by(gene_name) %>%
    mutate(entry_id = row_number()) %>%
    ungroup()
  
  cell_type_final_magical_df_for_paper_long <- cell_type_final_magical_df_for_paper %>%
    pivot_longer(cols = c(gene_fc, site_fc), names_to = "type", values_to = "value")
  
  cell_type_final_magical_df_for_paper_long <- cell_type_final_magical_df_for_paper_long %>%
    mutate(gene_entry = paste(gene_name, entry_id, sep = "_"))
  
  cell_type_final_magical_df_for_paper_long$gene_type <- factor(cell_type_final_magical_df_for_paper_long$gene_type, levels = c("Interleukin", "Interferon", "AP-1", "MAP Kinase",
                                                  "Lysine Methyltransferase", "Lysine Demethylase"))
  
  # Sort the data frame by the gene_type column
  cell_type_final_magical_df_for_paper_long <- cell_type_final_magical_df_for_paper_long %>%
    arrange(gene_type)
  
  cell_type_final_magical_df_for_paper_long$gene_entry <- factor(cell_type_final_magical_df_for_paper_long$gene_entry, levels = rev(unique(cell_type_final_magical_df_for_paper_long$gene_entry)))
  

  magical_heatmap_plot <- ggplot() + 
    geom_tile(data = cell_type_final_magical_df_for_paper_long, aes(x = type, y = gene_entry, fill = value), color = "black") +
    scale_fill_gradient2(low="navy", mid="white", high="red") +
    theme_classic(base_size = 14) + labs(title = NULL,
                                         x = NULL,
                                         y = NULL, fill = "Fold Change") + theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank())
  
  magical_heatmap_plots[[cell_type]] <- magical_heatmap_plot
}

ggsave(filename = paste0("C:/Users/wat2/Desktop/CD14_Mono_magical_heatmap.png"), plot = magical_heatmap_plots[["CD14_Mono"]], device='png', dpi=300, width = 3.5, height = 6, units = "in")

# Without legend
magical_heatmap_plots_without_legend <- list()

for(cell_type in unique(final_magical_df_for_paper$cell_type)) {
  cell_type_final_magical_df_for_paper <- final_magical_df_for_paper[final_magical_df_for_paper$cell_type == cell_type,]
  
  cell_type_final_magical_df_for_paper$gene_type <- factor(cell_type_final_magical_df_for_paper$gene_type, levels = c("Interleukin", "Interferon", "AP-1", "MAP Kinase", "Lysine Methyltransferase", "Lysine Demethylase", "Lysine Acetyltransferase"))
  
  # Organize the rows
  cell_type_final_magical_df_for_paper <- cell_type_final_magical_df_for_paper %>%
    arrange(gene_type, desc(gene_fc), gene_name, desc(site_fc))
  
  cell_type_final_magical_df_for_paper <- cell_type_final_magical_df_for_paper %>%
    group_by(gene_name) %>%
    mutate(entry_id = row_number()) %>%
    ungroup()
  
  cell_type_final_magical_df_for_paper_long <- cell_type_final_magical_df_for_paper %>%
    pivot_longer(cols = c(gene_fc, site_fc), names_to = "type", values_to = "value")
  
  cell_type_final_magical_df_for_paper_long <- cell_type_final_magical_df_for_paper_long %>%
    mutate(gene_entry = paste(gene_name, entry_id, sep = "_"))
  
  cell_type_final_magical_df_for_paper_long$gene_type <- factor(cell_type_final_magical_df_for_paper_long$gene_type, levels = c("Interleukin", "Interferon", "AP-1", "MAP Kinase",
                                                                                                                                "Lysine Methyltransferase", "Lysine Demethylase"))
  
  # Sort the data frame by the gene_type column
  cell_type_final_magical_df_for_paper_long <- cell_type_final_magical_df_for_paper_long %>%
    arrange(gene_type)
  
  cell_type_final_magical_df_for_paper_long$gene_entry <- factor(cell_type_final_magical_df_for_paper_long$gene_entry, levels = rev(unique(cell_type_final_magical_df_for_paper_long$gene_entry)))
  
  
  magical_heatmap_plot <- ggplot() + 
    geom_tile(data = cell_type_final_magical_df_for_paper_long, aes(x = type, y = gene_entry, fill = value), color = "black") +
    scale_fill_gradient2(low="navy", mid="white", high="red") +
    theme_classic(base_size = 14) + labs(title = NULL,
                                         x = NULL,
                                         y = NULL, fill = "Fold Change") + theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank()) + guides(fill="none")
  
  magical_heatmap_plots_without_legend[[cell_type]] <- magical_heatmap_plot
}

ggsave(filename = paste0("C:/Users/wat2/Desktop/CD14_Mono_magical_heatmap_without_legend.png"), plot = magical_heatmap_plots_without_legend[["CD14_Mono"]], device='png', dpi=300, width = 1, height = 6, units = "in")

# Other TF-based stuff
magical_tf_vs_gene_df <- create_tf_targets_df(overall_magical_df)
important_magical_tf_vs_gene_df <- magical_tf_vs_gene_df[magical_tf_vs_gene_df$Gene_symbol %in% important_overall_magical_df$Gene_symbol,]

# Add gene types
gene_types <- c()
for(current_row_index in 1:nrow(important_magical_tf_vs_gene_df)) {
  current_row <- important_magical_tf_vs_gene_df[current_row_index,]
  gene_types <- c(gene_types, heatmap_gene_types[[current_row$Gene_symbol]])
}

important_magical_tf_vs_gene_df$gene_type <- gene_types

important_magical_tf_vs_gene_df <- important_magical_tf_vs_gene_df[important_magical_tf_vs_gene_df$Cell_Type == "CD14_Mono",]
important_magical_tf_vs_gene_df <- important_magical_tf_vs_gene_df[important_magical_tf_vs_gene_df$gene_type == "MAP Kinase",]
