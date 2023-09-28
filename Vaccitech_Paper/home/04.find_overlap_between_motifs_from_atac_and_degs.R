# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

CD16_Mono_cell_motifs_up <- read.table(paste0(sc_peak_dir, "diff_peaks/CD16_Mono_D28_D1_motif_up_all_positive.tsv"), sep = "\t", header = TRUE)
CD16_Mono_cell_motifs_up <- CD16_Mono_cell_motifs_up[CD16_Mono_cell_motifs_up$p_adj < 0.05,]
CD16_Mono_cell_tfs_up <- CD16_Mono_cell_motifs_up$TF
CD16_Mono_cell_tfs_up <- sub("_.*", "", CD16_Mono_cell_tfs_up)

CD14_Mono_cell_motifs_down <- read.table(paste0(sc_peak_dir, "diff_peaks/CD14_Mono_D28_D1_motif_down_all_negative.tsv"), sep = "\t", header = TRUE)
CD14_Mono_cell_motifs_down <- CD14_Mono_cell_motifs_down[CD14_Mono_cell_motifs_down$p_adj < 0.05,]
CD14_Mono_cell_tfs_down <- CD14_Mono_cell_motifs_down$TF
CD14_Mono_cell_tfs_down <- sub("_.*", "", CD14_Mono_cell_tfs_down)

