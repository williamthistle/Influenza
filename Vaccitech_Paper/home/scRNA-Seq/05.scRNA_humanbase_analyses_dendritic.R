### THIS FUNCTION:
### 1) ASSIGNS CELL TYPES TO HUMANBASE MODULES (UP AND DOWNREGULATED SC DEGS)
### 2) ASSIGNS CELL TYPES TO REACTOME PATHWAYS (MODULES FROM HB ANALYSIS)
### 3) CREATES REACTOME PATHWAY PLOTS FOR FIGURE (NOTE: DOESN'T INCLUDE CELL TYPES YET)

# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

dendritic_dir <- paste0(sc_humanbase_dir, "Dendritic/")

# Grab upregulated and downregulated cell types for adaptive (trained immunity) cell types
sc_Dendritic <- sc_pseudobulk_deg_table[sc_pseudobulk_deg_table$Cell_Type == "cDC" | sc_pseudobulk_deg_table$Cell_Type == "pDC",]
sc_Dendritic_upregulated <- sc_Dendritic[sc_Dendritic$sc_log2FC > 0,]
sc_Dendritic_downregulated <- sc_Dendritic[sc_Dendritic$sc_log2FC < 0,]

# Assign cell types to upregulated genes from scRNA-seq (blood network)
HB_dendritic_upregulated_genes_dendritic <- assign_cell_types_to_humanbase_results(paste0(dendritic_dir, "HB_dendritic_upregulated_genes_dendritic.tsv"), sc_Dendritic_upregulated)
HB_dendritic_upregulated_genes_dendritic_modules <- HB_dendritic_upregulated_genes_dendritic[[1]]
HB_dendritic_upregulated_genes_dendritic_go_terms <- HB_dendritic_upregulated_genes_dendritic[[2]]

output_file <- paste0(dendritic_dir, "HB_dendritic_upregulated_genes_dendritic_with_cell_types.tsv")
cat(paste0("# GENES ASSOCIATED WITH EACH MODULE \n"), file=output_file)
utils::write.table(HB_dendritic_upregulated_genes_dendritic_modules, file = output_file, append=TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
cat(paste0("# GENES ASSOCIATED WITH EACH GO TERM WITHIN EACH MODULE \n"), file=output_file, append = TRUE)
utils::write.table(HB_dendritic_upregulated_genes_dendritic_go_terms, file = output_file, append=TRUE, sep = "\t", quote = FALSE, row.names = FALSE)

# Assign cell types to downregulated genes from scRNA-seq (dendritic network)
HB_dendritic_downregulated_genes_dendritic <- assign_cell_types_to_humanbase_results(paste0(dendritic_dir, "HB_Dendritic_downregulated_genes_dendritic.tsv"), sc_Dendritic_downregulated)
HB_dendritic_downregulated_genes_dendritic_modules <- HB_dendritic_downregulated_genes_dendritic[[1]]
HB_dendritic_downregulated_genes_dendritic_go_terms <- HB_dendritic_downregulated_genes_dendritic[[2]]

output_file <- paste0(dendritic_dir, "HB_dendritic_downregulated_genes_dendritic_with_cell_types.tsv")
cat(paste0("# GENES ASSOCIATED WITH EACH MODULE \n"), file=output_file)
utils::write.table(HB_dendritic_downregulated_genes_dendritic_modules, file = output_file, append=TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
cat(paste0("# GENES ASSOCIATED WITH EACH GO TERM WITHIN EACH MODULE \n"), file=output_file, append = TRUE)
utils::write.table(HB_dendritic_downregulated_genes_dendritic_go_terms, file = output_file, append=TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
