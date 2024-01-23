### THIS FUNCTION:
### 1) ASSIGNS CELL TYPES TO HUMANBASE MODULES (UP AND DOWNREGULATED SC DEGS)
### 2) ASSIGNS CELL TYPES TO REACTOME PATHWAYS (MODULES FROM HB ANALYSIS)
### 3) CREATES REACTOME PATHWAY PLOTS FOR FIGURE (NOTE: DOESN'T INCLUDE CELL TYPES YET)

# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

monocyte_dir <- paste0(sc_humanbase_dir, "Monocytes/")

# Grab upregulated and downregulated cell types for adaptive (trained immunity) cell types
sc_monocytes <- sc_pseudobulk_deg_table[sc_pseudobulk_deg_table$Cell_Type == "CD14 Mono" | sc_pseudobulk_deg_table$Cell_Type == "CD16 Mono",]
sc_monocytes_upregulated <- sc_monocytes[sc_monocytes$sc_log2FC > 0,]
sc_monocytes_downregulated <- sc_monocytes[sc_monocytes$sc_log2FC < 0,]

# Assign cell types to upregulated genes from scRNA-seq (blood network)
HB_monocyte_upregulated_genes_monocyte <- assign_cell_types_to_humanbase_results(paste0(monocyte_dir, "HB_monocytes_upregulated_genes_monocyte.tsv"), sc_monocytes_upregulated)
HB_monocyte_upregulated_genes_monocyte_modules <- HB_monocyte_upregulated_genes_monocyte[[1]]
HB_monocyte_upregulated_genes_monocyte_go_terms <- HB_monocyte_upregulated_genes_monocyte[[2]]

output_file <- paste0(monocyte_dir, "HB_monocyte_upregulated_genes_monocyte_with_cell_types.tsv")
cat(paste0("# GENES ASSOCIATED WITH EACH MODULE \n"), file=output_file)
utils::write.table(HB_monocyte_upregulated_genes_monocyte_modules, file = output_file, append=TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
cat(paste0("# GENES ASSOCIATED WITH EACH GO TERM WITHIN EACH MODULE \n"), file=output_file, append = TRUE)
utils::write.table(HB_monocyte_upregulated_genes_monocyte_go_terms, file = output_file, append=TRUE, sep = "\t", quote = FALSE, row.names = FALSE)

# Assign cell types to downregulated genes from scRNA-seq (monocyte network)
HB_monocyte_downregulated_genes_monocyte <- assign_cell_types_to_humanbase_results(paste0(monocyte_dir, "HB_monocytes_downregulated_genes_monocyte.tsv"), sc_monocytes_downregulated)
HB_monocyte_downregulated_genes_monocyte_modules <- HB_monocyte_downregulated_genes_monocyte[[1]]
HB_monocyte_downregulated_genes_monocyte_go_terms <- HB_monocyte_downregulated_genes_monocyte[[2]]

output_file <- paste0(monocyte_dir, "HB_monocyte_downregulated_genes_monocyte_with_cell_types.tsv")
cat(paste0("# GENES ASSOCIATED WITH EACH MODULE \n"), file=output_file)
utils::write.table(HB_monocyte_downregulated_genes_monocyte_modules, file = output_file, append=TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
cat(paste0("# GENES ASSOCIATED WITH EACH GO TERM WITHIN EACH MODULE \n"), file=output_file, append = TRUE)
utils::write.table(HB_monocyte_downregulated_genes_monocyte_go_terms, file = output_file, append=TRUE, sep = "\t", quote = FALSE, row.names = FALSE)

# Assign cell types to genes in pathway analyses (from Natalie's script)
# Note that M2 didn't have any significant results
up_modules <- c("M1", "M3", "M4", "M5")
up_pathway_plot <- data.frame(Term = character(), Count = numeric(), Adjusted.P.value = numeric(), Genes = character(), Cell_Types = character(), Module = character())
for(module in up_modules) {   
  module_pathway_info <- read.table(paste0(monocyte_dir, "HB_monocyte_upregulated_genes_monocyte_", module, "_gsea_info_output.tsv"), sep = "\t", header = TRUE)
  module_pathway_info_with_cell_types <- assign_cell_types_to_pathway_results(module_pathway_info, monocyte_sc_pseudobulk_deg_table[monocyte_sc_pseudobulk_deg_table$sc_log2FC > 0,])
  write.table(module_pathway_info_with_cell_types, paste0(monocyte_dir, "HB_monocyte_upregulated_genes_monocyte_", module, "_gsea_info_output_with_cell_types.tsv"), sep = "\t", quote = FALSE)
  module_pathway_info_with_cell_types_edited <- module_pathway_info_with_cell_types[,c("Term", "Overlap", "Adjusted.P.value", "Genes", "Cell_Types")]
  colnames(module_pathway_info_with_cell_types_edited) <- c("Term", "Count", "Adjusted.P.Value", "Genes", "Cell_Types")
  module_pathway_info_with_cell_types_edited <- module_pathway_info_with_cell_types_edited[1:5,]
  module_pathway_info_with_cell_types_edited$Count <- unlist(strsplit(module_pathway_info_with_cell_types_edited$Count, "/"))[c(1, 3, 5, 7, 9)]
  module_pathway_info_with_cell_types_edited$Module <- module
  up_pathway_plot <- rbind(up_pathway_plot, module_pathway_info_with_cell_types_edited)
}

write.table(up_pathway_plot, paste0(monocyte_dir, "HB_monocyte_upregulated_genes_monocyte_gsea_top_upregulated_pathways.tsv"), sep = "\t", quote = FALSE)

down_modules <- c("M1", "M2", "M3", "M4")
down_pathway_plot <- data.frame(Term = character(), Count = numeric(), Adjusted.P.value = numeric(), Genes = character(), Cell_Types = character(), Module = character())
for(module in down_modules) {
  module_pathway_info <- read.table(paste0(monocyte_dir, "HB_monocyte_downregulated_genes_monocyte_", module, "_gsea_info_output.tsv"), sep = "\t", header = TRUE)
  module_pathway_info_with_cell_types <- assign_cell_types_to_pathway_results(module_pathway_info, monocyte_sc_pseudobulk_deg_table[monocyte_sc_pseudobulk_deg_table$sc_log2FC < 0,])
  write.table(module_pathway_info_with_cell_types, paste0(monocyte_dir, "HB_monocyte_downregulated_genes_monocyte_", module, "_gsea_info_output_with_cell_types.tsv"), sep = "\t", quote = FALSE)
  module_pathway_info_with_cell_types_edited <- module_pathway_info_with_cell_types[,c("Term", "Overlap", "Adjusted.P.value", "Genes", "Cell_Types")]
  colnames(module_pathway_info_with_cell_types_edited) <- c("Term", "Count", "Adjusted.P.Value", "Genes", "Cell_Types")
  module_pathway_info_with_cell_types_edited <- module_pathway_info_with_cell_types_edited[1:5,]
  module_pathway_info_with_cell_types_edited$Count <- unlist(strsplit(module_pathway_info_with_cell_types_edited$Count, "/"))[c(1, 3, 5, 7, 9)]
  module_pathway_info_with_cell_types_edited$Module <- module
  down_pathway_plot <- rbind(down_pathway_plot, module_pathway_info_with_cell_types_edited)
}

write.table(down_pathway_plot, paste0(monocyte_dir, "HB_monocyte_upregulated_genes_monocyte_gsea_top_downregulated_pathways.tsv"), sep = "\t", quote = FALSE)


# Create plot for upregulated module pathway analysis (Figure 2)
colnames(up_pathway_plot) <- c("Reactome.Pathway", "Gene.Count", "Adjusted.P.Value", "Genes", "Cell.Types", "Module")
up_pathway_plot$Reactome.Pathway <- sub("\\R-HSA-.*", "", up_pathway_plot$Reactome.Pathway)
up_pathway_plot$Gene.Count <- as.numeric(up_pathway_plot$Gene.Count)
up_pathway_plot$Reactome.Pathway <- factor(up_pathway_plot$Reactome.Pathway, levels = unique(up_pathway_plot$Reactome.Pathway))

up_pathway_plot_file <- ggplot(data = up_pathway_plot, aes(x = Module, y = Reactome.Pathway, size = Gene.Count, color = Adjusted.P.Value)) +
  geom_point() +
  theme_minimal() +
  labs(
    title = "Reactome Pathways for Upregulated Genes from monocyte Immune Cells",
    x = "Module",
    y = "Reactome Pathway",
    size = "Count",
    color = "Adjusted P Value"
  ) +
  theme(plot.title = element_text(hjust = 1)) + theme(aspect.ratio = 3/1)

ggsave(filename = paste0(monocyte_dir, "HB_monocyte_upregulated_genes_monocyte_gsea_top_upregulated_pathways.tiff"), plot = up_pathway_plot_file, device='tiff', dpi=300)

# Create plot for downregulated module pathway analysis (Figure 2)
colnames(down_pathway_plot) <- c("Reactome.Pathway", "Gene.Count", "Adjusted.P.Value", "Genes", "Cell.Types", "Module")
down_pathway_plot$Reactome.Pathway <- sub("\\R-HSA-.*", "", down_pathway_plot$Reactome.Pathway)
down_pathway_plot$Gene.Count <- as.numeric(down_pathway_plot$Gene.Count)
down_pathway_plot$Reactome.Pathway <- factor(down_pathway_plot$Reactome.Pathway, levels = unique(down_pathway_plot$Reactome.Pathway))

down_pathway_plot_file <- ggplot(data = down_pathway_plot, aes(x = Module, y = Reactome.Pathway, size = Gene.Count, color = Adjusted.P.Value)) +
  geom_point() +
  theme_minimal() +
  labs(
    title = "Reactome Pathways for Downregulated Genes from monocyte Immune Cells",
    x = "Module",
    y = "Reactome Pathway",
    size = "Count",
    color = "Adjusted P Value"
  ) +
  theme(plot.title = element_text(hjust = 1)) + theme(aspect.ratio = 3/1)

ggsave(filename = paste0(monocyte_dir, "HB_monocyte_downregulated_genes_monocyte_gsea_top_downregulated_pathways.tiff"), plot = down_pathway_plot_file, device='tiff', dpi=300)

