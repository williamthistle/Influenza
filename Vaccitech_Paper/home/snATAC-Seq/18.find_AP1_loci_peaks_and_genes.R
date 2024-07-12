# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

txdb_hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene
promoter_terms <- c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)")

peak_motif_matches <- read.table(paste0(scATAC_hvl_dir, "Placebo/peak_motif_matches.txt"), sep = "\t", header = TRUE)

peak_motif_matches_AP1 <- peak_motif_matches[peak_motif_matches$FOS > 0 | peak_motif_matches$FOSL1 > 0 | peak_motif_matches$FOSL2 > 0
                                             | peak_motif_matches$JUN > 0 | peak_motif_matches$JUNB > 0 | peak_motif_matches$JUND > 0,]
peak_motif_matches_AP1$chr <- paste0("chr", peak_motif_matches_AP1$chr)
peak_motif_matches_AP1$chr <- gsub("chr23", "chrX", peak_motif_matches_AP1$chr)
# CD14 Mono
scATAC_CD14_mono_peaks <- read.table(paste0(scATAC_hvl_placebo_das_dir, "D28-vs-D_minus_1-degs-CD14_Mono-time_point-controlling_for_subject_id_final_pct_0.01.tsv"), sep = "\t", header = TRUE)

chromosomes <- sapply(strsplit(scATAC_CD14_mono_peaks$Peak_Name, "-"), `[`, 1)
start_coords <- as.numeric(sapply(strsplit(scATAC_CD14_mono_peaks$Peak_Name, "-"), `[`, 2))
end_coords <- as.numeric(sapply(strsplit(scATAC_CD14_mono_peaks$Peak_Name, "-"), `[`, 3))
scATAC_CD14_mono_peaks$seqnames <- chromosomes
scATAC_CD14_mono_peaks$start <- start_coords
scATAC_CD14_mono_peaks$end <- end_coords

# CD14 Mono AP-1 subset
scATAC_CD14_mono_peaks_AP1 <- scATAC_CD14_mono_peaks[scATAC_CD14_mono_peaks$seqnames %in% peak_motif_matches_AP1$chr &
                                                       scATAC_CD14_mono_peaks$start %in% peak_motif_matches_AP1$point1 &
                                                       scATAC_CD14_mono_peaks$end %in% peak_motif_matches_AP1$point2,]

scATAC_CD14_mono_peaks_AP1_annotated <- annotatePeak(makeGRangesFromDataFrame(scATAC_CD14_mono_peaks_AP1), TxDb = txdb_hg38, annoDb = "org.Hs.eg.db")
scATAC_CD14_mono_peaks_AP1_annotated <- as.data.frame(scATAC_CD14_mono_peaks_AP1_annotated)
scATAC_CD14_mono_peaks_AP1_annotated$sc_log2FC <- scATAC_CD14_mono_peaks_AP1$sc_log2FC
scATAC_CD14_mono_peaks_AP1_annotated$sc_pval <- scATAC_CD14_mono_peaks_AP1$sc_pval