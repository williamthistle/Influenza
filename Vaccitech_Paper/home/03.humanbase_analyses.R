# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# Assign cell types to upregulated genes from scRNA-seq (blood network)
HB_all_upregulated_genes_blood <- assign_cell_types_to_humanbase_results(paste0(sc_humanbase_dir, "HB_all_upregulated_genes_blood.tsv"), sc_pseudobulk_gene_table)
HB_all_upregulated_genes_blood_modules <- HB_all_upregulated_genes_blood[[1]]
HB_all_upregulated_genes_blood_go_terms <- HB_all_upregulated_genes_blood[[2]]

output_file <- paste0(sc_humanbase_dir, "HB_all_upregulated_genes_blood_with_cell_types.tsv")
cat(paste0("# GENES ASSOCIATED WITH EACH MODULE \n"), file=output_file)
utils::write.table(HB_all_upregulated_genes_blood_modules, file = output_file, append=TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
cat(paste0("# GENES ASSOCIATED WITH EACH GO TERM WITHIN EACH MODULE \n"), file=output_file, append = TRUE)
utils::write.table(HB_all_upregulated_genes_blood_go_terms, file = output_file, append=TRUE, sep = "\t", quote = FALSE, row.names = FALSE)

# Assign cell types to downregulated genes from scRNA-seq (blood network)
HB_all_downregulated_genes_blood <- assign_cell_types_to_humanbase_results(paste0(sc_humanbase_dir, "HB_all_downregulated_genes_blood.tsv"), sc_pseudobulk_gene_table)
HB_all_downregulated_genes_blood_modules <- HB_all_downregulated_genes_blood[[1]]
HB_all_downregulated_genes_blood_go_terms <- HB_all_downregulated_genes_blood[[2]]

output_file <- paste0(sc_humanbase_dir, "HB_all_downregulated_genes_blood_with_cell_types.tsv")
cat(paste0("# GENES ASSOCIATED WITH EACH MODULE \n"), file=output_file)
utils::write.table(HB_all_downregulated_genes_blood_modules, file = output_file, append=TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
cat(paste0("# GENES ASSOCIATED WITH EACH GO TERM WITHIN EACH MODULE \n"), file=output_file, append = TRUE)
utils::write.table(HB_all_downregulated_genes_blood_go_terms, file = output_file, append=TRUE, sep = "\t", quote = FALSE, row.names = FALSE)

# Assign cell types to genes in pathway analyses (from Natalie's script)
up_modules <- c("M1", "M2", "M3", "M4", "M5")
down_modules <- c("M1", "M2", "M3", "M4", "M5", "M6")
for(module in up_modules) {
  module_pathway_info <- read.table(paste0(sc_humanbase_dir, "HB_all_upregulated_genes_blood_", module, "_gsea_info_output.tsv"), sep = "\t", header = TRUE)
  module_pathway_info_with_cell_types <- assign_cell_types_to_pathway_results(module_pathway_info, sc_pseudobulk_gene_table)
  write.table(module_pathway_info_with_cell_types, paste0(sc_humanbase_dir, "HB_all_upregulated_genes_blood_", module, "_gsea_info_output_with_cell_types.tsv"), sep = "\t", quote = FALSE)
}

for(module in down_modules) {
  module_pathway_info <- read.table(paste0(sc_humanbase_dir, "HB_all_downregulated_genes_blood_", module, "_gsea_info_output.tsv"), sep = "\t", header = TRUE)
  module_pathway_info_with_cell_types <- assign_cell_types_to_pathway_results(module_pathway_info, sc_pseudobulk_gene_table)
  write.table(module_pathway_info_with_cell_types, paste0(sc_humanbase_dir, "HB_all_downregulated_genes_blood_", module, "_gsea_info_output_with_cell_types.tsv"), sep = "\t", quote = FALSE)
}