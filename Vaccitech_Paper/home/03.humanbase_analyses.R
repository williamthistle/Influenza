# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# Grab upregulated and downregulated cell types for innate (trained immunity) cell types
sc_trained_immunity_upregulated_genes <- unique(innate_sc_pseudobulk_deg_table[innate_sc_pseudobulk_deg_table$sc_log2FC > 0,]$Gene_Name)
sc_trained_immunity_downregulated_genes <- unique(innate_sc_pseudobulk_deg_table[innate_sc_pseudobulk_deg_table$sc_log2FC < 0,]$Gene_Name)

# Assign cell types to upregulated genes from scRNA-seq (blood network)
HB_innate_upregulated_genes_blood <- assign_cell_types_to_humanbase_results(paste0(sc_humanbase_dir, "HB_innate_upregulated_genes_blood.tsv"), sc_pseudobulk_gene_table)
HB_innate_upregulated_genes_blood_modules <- HB_innate_upregulated_genes_blood[[1]]
HB_innate_upregulated_genes_blood_go_terms <- HB_innate_upregulated_genes_blood[[2]]

output_file <- paste0(sc_humanbase_dir, "HB_innate_upregulated_genes_blood_with_cell_types.tsv")
cat(paste0("# GENES ASSOCIATED WITH EACH MODULE \n"), file=output_file)
utils::write.table(HB_innate_upregulated_genes_blood_modules, file = output_file, append=TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
cat(paste0("# GENES ASSOCIATED WITH EACH GO TERM WITHIN EACH MODULE \n"), file=output_file, append = TRUE)
utils::write.table(HB_innate_upregulated_genes_blood_go_terms, file = output_file, append=TRUE, sep = "\t", quote = FALSE, row.names = FALSE)

# Assign cell types to downregulated genes from scRNA-seq (blood network)
HB_innate_downregulated_genes_blood <- assign_cell_types_to_humanbase_results(paste0(sc_humanbase_dir, "HB_innate_downregulated_genes_blood.tsv"), sc_pseudobulk_gene_table)
HB_innate_downregulated_genes_blood_modules <- HB_innate_downregulated_genes_blood[[1]]
HB_innate_downregulated_genes_blood_go_terms <- HB_innate_downregulated_genes_blood[[2]]

output_file <- paste0(sc_humanbase_dir, "HB_innate_downregulated_genes_blood_with_cell_types.tsv")
cat(paste0("# GENES ASSOCIATED WITH EACH MODULE \n"), file=output_file)
utils::write.table(HB_innate_downregulated_genes_blood_modules, file = output_file, append=TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
cat(paste0("# GENES ASSOCIATED WITH EACH GO TERM WITHIN EACH MODULE \n"), file=output_file, append = TRUE)
utils::write.table(HB_innate_downregulated_genes_blood_go_terms, file = output_file, append=TRUE, sep = "\t", quote = FALSE, row.names = FALSE)

# Assign cell types to genes in pathway analyses (from Natalie's script)


up_modules <- c("M1", "M2", "M3", "M4", "M5", "M6")
up_pathway_plot <- data.frame(Term = character(), Count = numeric(), Adjusted.P.value = numeric(), Genes = character(), Cell_Types = character(), Module = character())
for(module in up_modules) {   
  module_pathway_info <- read.table(paste0(sc_humanbase_dir, "HB_innate_upregulated_genes_blood_", module, "_gsea_info_output.tsv"), sep = "\t", header = TRUE)
  module_pathway_info_with_cell_types <- assign_cell_types_to_pathway_results(module_pathway_info, innate_sc_pseudobulk_deg_table)
  write.table(module_pathway_info_with_cell_types, paste0(sc_humanbase_dir, "HB_innate_upregulated_genes_blood_", module, "_gsea_info_output_with_cell_types.tsv"), sep = "\t", quote = FALSE)
  module_pathway_info_with_cell_types_edited <- module_pathway_info_with_cell_types[,c("Term", "Overlap", "Adjusted.P.value", "Genes", "Cell_Types")]
  colnames(module_pathway_info_with_cell_types_edited) <- c("Term", "Count", "Adjusted.P.Value", "Genes", "Cell_Types")
  module_pathway_info_with_cell_types_edited <- module_pathway_info_with_cell_types_edited[1:5,]
  module_pathway_info_with_cell_types_edited$Count <- unlist(strsplit(module_pathway_info_with_cell_types_edited$Count, "/"))[c(1, 3, 5, 7, 9)]
  module_pathway_info_with_cell_types_edited$Module <- module
  up_pathway_plot <- rbind(up_pathway_plot, module_pathway_info_with_cell_types_edited)
}

write.table(up_pathway_plot, paste0(sc_humanbase_dir, "HB_innate_upregulated_genes_blood_gsea_top_upregulated_pathways.tsv"), sep = "\t", quote = FALSE)

ggplot(data = up_pathway_plot, aes(x = Module, y = Term, size = Count, color = Adjusted.P.Value)) +
  geom_point() +
  scale_size_discrete(range = c(2, 16)) +
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "orange")) +
  theme_minimal() +
  labs(
    title = "Scatter Plot",
    x = "Module",
    y = "Term",
    size = "Count",
    color = "Adjusted P-Value"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



down_modules <- c("M1", "M2", "M3", "M4")
for(module in down_modules) {
  module_pathway_info <- read.table(paste0(sc_humanbase_dir, "HB_innate_downregulated_genes_blood_", module, "_gsea_info_output.tsv"), sep = "\t", header = TRUE)
  module_pathway_info_with_cell_types <- assign_cell_types_to_pathway_results(module_pathway_info, innate_sc_pseudobulk_deg_table)
  write.table(module_pathway_info_with_cell_types, paste0(sc_humanbase_dir, "HB_innate_downregulated_genes_blood_", module, "_gsea_info_output_with_cell_types.tsv"), sep = "\t", quote = FALSE)
}

# Create plot for upregulated module pathway analysis (Figure 2)



# Create plot for downregulated module pathway analysis (Figure 2)


