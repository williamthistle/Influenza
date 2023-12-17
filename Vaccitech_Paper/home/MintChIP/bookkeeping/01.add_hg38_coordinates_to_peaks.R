# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

mintchip_markers <- c("H3K4me1", "H3K4me3", "H3K9me3", "H3K27Ac", "H3K27me3", "H3K36me3")
for(marker in mintchip_markers) {
  current_hg38_sites <- read.table(paste0(mintchip_das_dir, marker, "/", marker, "_all_peaks_hg38.tsv"), sep = "\t", header = TRUE)
  # colnames(current_hg38_sites) <- c("seqnames", "start", "end", "hg38_coordinates", "overlaps")
  write.table(current_hg38_sites, file = paste0(mintchip_das_dir, marker, "/", marker, "_all_peaks_with_hg38_coordinates.tsv"), sep = "\t", quote = FALSE,
              row.names = FALSE)
}