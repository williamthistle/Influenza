# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# Assign cell types to upregulated genes from scRNA-seq (blood network)
HB_all_upregulated_genes_blood <- assign_cell_types_to_humanbase_results(paste0(sc_humanbase_dir, "HB_all_upregulated_genes_blood.tsv"), sc_pseudobulk_gene_table)
HB_all_upregulated_genes_blood_modules <- HB_all_upregulated_genes_blood[[1]]
HB_all_upregulated_genes_blood_go_terms <- HB_all_upregulated_genes_blood[[2]]

# Assign cell types to downregulated genes from scRNA-seq (blood network)
HB_all_downregulated_genes_blood <- assign_cell_types_to_humanbase_results(paste0(sc_humanbase_dir, "HB_all_downregulated_genes_blood.tsv"), sc_pseudobulk_gene_table)
HB_all_downregulated_genes_blood_modules <- HB_all_downregulated_genes_blood[[1]]
HB_all_downregulated_genes_blood_go_terms <- HB_all_downregulated_genes_blood[[2]]

