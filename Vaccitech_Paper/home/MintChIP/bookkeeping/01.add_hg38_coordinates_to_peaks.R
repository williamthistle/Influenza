# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

mintchip_markers <- c("H3K4me1", "H3K4me3", "H3K9me3", "H3K27Ac", "H3K27me3", "H3K36me3")
for(marker in mintchip_markers) {
  print(marker)
  # Read in hg19 (original) mintchip peaks
  current_hg19_peaks <- read.table(paste0(mintchip_das_dir, marker, "/", marker, "_all_peaks.tsv"), sep = "\t")
  colnames(current_hg19_peaks) <- c("seqnames", "start", "end")
  # Read in hg38 (new) associated peaks
  current_hg38_peaks <- read.table(paste0(mintchip_das_dir, marker, "/", marker, "_all_peaks_hg38.tsv"), sep = "\t")
  colnames(current_hg38_peaks) <- c("seqnames", "start", "end")
  # Read in errors that occurred while mapping peaks
  peak_transfer_errors <- read.table(paste0(mintchip_das_dir, marker, "/", marker, "_all_peaks_error.tsv"), sep = "\t")
  colnames(peak_transfer_errors) <- c("seqnames", "start", "end")
  # Remove any peaks that had issues from current_hg19_peaks
  current_hg19_peaks <- anti_join(current_hg19_peaks, peak_transfer_errors, by = c("seqnames", "start", "end"))
  # Now, all hg19 / hg38 peaks should be aligned row-by-row
  current_hg19_peaks$hg38_start <- current_hg38_peaks$start
  current_hg19_peaks$hg38_end <- current_hg38_peaks$end
  write.table(current_hg19_peaks, file = paste0(mintchip_das_dir, marker, "/", marker, "_all_peaks_with_hg38_coordinates.tsv"), sep = "\t", quote = FALSE,
              row.names = FALSE)
}