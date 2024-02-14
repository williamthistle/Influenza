base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# UPREGULATED SCRNA-SEQ
analysis_dir <- "C:/Users/willi/Desktop/full_sc_test/scRNA_upregulated_comparisons/"
fc_direction <- "upregulated"

get_lists_of_commands(analysis_dir, fc_direction)

# DOWNREGULATED SCRNA-SEQ
analysis_dir <- "C:/Users/willi/Desktop/full_sc_test/scRNA_downregulated_comparisons/"
fc_direction <- "downregulated"

get_lists_of_commands(analysis_dir, fc_direction)

# UPREGULATED BULK RNA-SEQ
analysis_dir <- "C:/Users/willi/OneDrive - Princeton University/Vaccitech_Paper/network_comparisons/bulk_upregulated_comparisons/"
fc_direction <- "upregulated"

get_lists_of_commands(analysis_dir, fc_direction)

# DOWNREGULATED BULK RNA-SEQ
analysis_dir <- "C:/Users/willi/OneDrive - Princeton University/Vaccitech_Paper/network_comparisons/bulk_downregulated_comparisons/"
fc_direction <- "downregulated"

get_lists_of_commands(analysis_dir, fc_direction)











