# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# CD14 Mono
cd14_mono_upregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "CD14_Mono_Upregulated.tsv"))
cd14_mono_downregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "CD14_Mono_Downregulated.tsv"))
# CD16 Mono
cd16_mono_upregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "CD16_Mono_Upregulated.tsv"))
cd16_mono_downregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "CD16_Mono_Downregulated.tsv"))
# NK
nk_upregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "NK_Upregulated.tsv"))
nk_downregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "NK_Downregulated.tsv"))
# cDC
cDC_upregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "cDC_Upregulated.tsv"))
cDC_downregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "cDC_Downregulated.tsv"))
# CD4 Memory
cd4_memory_upregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "CD4_Memory_Upregulated.tsv"))
cd4_memory_downregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "CD4_Memory_Downregulated.tsv"))
# CD8 Memory
cd8_memory_upregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "CD8_Memory_Upregulated.tsv"))
cd8_memory_downregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "CD8_Memory_Downregulated.tsv"))
# CD8 Naive
cd8_naive_upregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "CD8_Naive_Upregulated.tsv"))
cd8_naive_downregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "CD8_Naive_Downregulated.tsv"))
# MAIT
mait_upregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "MAIT_Upregulated.tsv"))
mait_downregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "MAIT_Downregulated.tsv"))
# B Memory
b_memory_upregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "B_Memory_Upregulated.tsv"))
b_memory_downregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "B_Memory_Downregulated.tsv"))
# B Naive
b_naive_upregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "B_Naive_Upregulated.tsv"))
b_naive_downregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "B_Naive_Downregulated.tsv"))