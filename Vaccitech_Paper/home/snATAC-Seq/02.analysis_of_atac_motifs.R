# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

pos_overlapping_peaks <- read.table("C:/Users/wat2/OneDrive - Princeton University/Vaccitech_Paper/Analyses/snATAC-Seq/Results/HVL/motif_enrichment/CD14_Mono/overlapping_peak/0.01/with_bg/D28-vs-D_minus_1-degs-CD14_Mono-overlapping_peak_pct_0.01_total_peaks_1041_pos_motifs_with_bg.tsv", sep = "\t", header = TRUE)
pos_overlapping_peaks <- pos_overlapping_peaks[pos_overlapping_peaks$p.adjust < 0.05,]
neg_overlapping_peaks <- read.table("C:/Users/wat2/OneDrive - Princeton University/Vaccitech_Paper/Analyses/snATAC-Seq/Results/HVL/motif_enrichment/CD14_Mono/overlapping_peak/0.01/with_bg/D28-vs-D_minus_1-degs-CD14_Mono-overlapping_peak_pct_0.01_total_peaks_1442_neg_motifs_with_bg.tsv", sep = "\t", header = TRUE)
neg_overlapping_peaks <- neg_overlapping_peaks[neg_overlapping_peaks$p.adjust < 0.05,]

pos_sc_peaks <- read.table("C:/Users/wat2/OneDrive - Princeton University/Vaccitech_Paper/Analyses/snATAC-Seq/Results/HVL/motif_enrichment/CD14_Mono/sc_filtered/0.01/with_bg/D28-vs-D_minus_1-degs-CD14_Mono-sc_filtered_pct_0.01_FC_1_total_peaks_2166_pos_motifs_with_bg.tsv", sep = "\t", header = TRUE)
pos_sc_peaks <- pos_sc_peaks[pos_sc_peaks$p.adjust < 0.05,]

neg_sc_peaks <- read.table("C:/Users/wat2/OneDrive - Princeton University/Vaccitech_Paper/Analyses/snATAC-Seq/Results/HVL/motif_enrichment/CD14_Mono/sc_filtered/0.01/with_bg/D28-vs-D_minus_1-degs-CD14_Mono-sc_filtered_pct_0.01_FC_-1_total_peaks_2397_neg_motifs_with_bg.tsv", sep = "\t", header = TRUE)
neg_sc_peaks <- neg_sc_peaks[neg_sc_peaks$p.adjust < 0.05,]

