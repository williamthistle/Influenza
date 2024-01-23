### THIS FUNCTION:
### 1) ASSIGNS CELL TYPES TO HUMANBASE MODULES (UP AND DOWNREGULATED SC DEGS)
### 2) ASSIGNS CELL TYPES TO REACTOME PATHWAYS (MODULES FROM HB ANALYSIS)
### 3) CREATES REACTOME PATHWAY PLOTS FOR FIGURE (NOTE: DOESN'T INCLUDE CELL TYPES YET)

# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

nk_dir <- paste0(sc_humanbase_dir, "NK/")

# Grab upregulated and downregulated cell types for adaptive (trained immunity) cell types
sc_nks <- sc_pseudobulk_deg_table[sc_pseudobulk_deg_table$Cell_Type == "NK" | sc_pseudobulk_deg_table$Cell_Type == "NK_CD56bright",]
sc_nks_upregulated <- sc_nks[sc_nks$sc_log2FC > 0,]
sc_nks_downregulated <- sc_nks[sc_nks$sc_log2FC < 0,]

# Assign cell types to upregulated genes from scRNA-seq (blood network)
HB_nk_upregulated_genes_nk <- assign_cell_types_to_humanbase_results(paste0(nk_dir, "HB_NKs_upregulated_genes_natural_killer_cell.tsv"), sc_nks_upregulated)
HB_nk_upregulated_genes_nk_modules <- HB_nk_upregulated_genes_nk[[1]]
HB_nk_upregulated_genes_nk_go_terms <- HB_nk_upregulated_genes_nk[[2]]

output_file <- paste0(nk_dir, "HB_nk_upregulated_genes_nk_with_cell_types.tsv")
cat(paste0("# GENES ASSOCIATED WITH EACH MODULE \n"), file=output_file)
utils::write.table(HB_nk_upregulated_genes_nk_modules, file = output_file, append=TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
cat(paste0("# GENES ASSOCIATED WITH EACH GO TERM WITHIN EACH MODULE \n"), file=output_file, append = TRUE)
utils::write.table(HB_nk_upregulated_genes_nk_go_terms, file = output_file, append=TRUE, sep = "\t", quote = FALSE, row.names = FALSE)

# Assign cell types to downregulated genes from scRNA-seq (nk network)
HB_nk_downregulated_genes_nk <- assign_cell_types_to_humanbase_results(paste0(nk_dir, "HB_NKs_downregulated_genes_natural_killer_cell.tsv"), sc_nks_downregulated)
HB_nk_downregulated_genes_nk_modules <- HB_nk_downregulated_genes_nk[[1]]
HB_nk_downregulated_genes_nk_go_terms <- HB_nk_downregulated_genes_nk[[2]]

output_file <- paste0(nk_dir, "HB_nk_downregulated_genes_nk_with_cell_types.tsv")
cat(paste0("# GENES ASSOCIATED WITH EACH MODULE \n"), file=output_file)
utils::write.table(HB_nk_downregulated_genes_nk_modules, file = output_file, append=TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
cat(paste0("# GENES ASSOCIATED WITH EACH GO TERM WITHIN EACH MODULE \n"), file=output_file, append = TRUE)
utils::write.table(HB_nk_downregulated_genes_nk_go_terms, file = output_file, append=TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
