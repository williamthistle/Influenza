# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# We will use consensus peak FC 0 threshold for liftover (hg19 -> hg38)

mintchip_markers <- c("H3K4me1", "H3K4me3", "H3K9me3", "H3K27Ac", "H3K27me3", "H3K36me3")

for(marker in mintchip_markers) {
  current_hg19_sites <- read.table(paste0(mintchip_das_dir, marker, "/", marker, "_consensus_peak_set_FC_0.tsv"), sep = "\t", header = TRUE)
  current_hg19_sites <- current_hg19_sites$coordinates
  write.table(current_hg19_sites, file = paste0(mintchip_das_dir, marker, "/", marker, "_consensus_peak_set_FC_0_coordinates.tsv"),
                                                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Go to UCSC site to perform liftover

# Add hg38 coordinates to each file
for(marker in mintchip_markers) {
  current_marker_sites <- read.table(paste0(mintchip_das_dir, marker, "/", marker, "_consensus_peak_set_FC_0.tsv"), sep = "\t", header = TRUE)
  associated_hg38_coordinates <- read.table(paste0(mintchip_das_dir, marker, "/", marker, "_consensus_peak_set_FC_0_coordinates_hg38.tsv"), 
                                            sep = "\t")$V1
  if(file.size(paste0(mintchip_das_dir, marker, "/", marker, "_consensus_peak_set_FC_0_coordinates_hg38_error.tsv")) > 2) {
    associated_hg38_errors <- read.table(paste0(mintchip_das_dir, marker, "/", marker, "_consensus_peak_set_FC_0_coordinates_hg38_error.tsv"), 
                                         sep = "\t")$V1
    # Remove any site that didn't have associated hg38 site in liftover
    current_marker_sites <- current_marker_sites[current_marker_sites$coordinates != associated_hg38_errors, ]
  }
  # Add hg38 info
  current_marker_sites$hg38_coordinates <- associated_hg38_coordinates
  write.table(current_marker_sites, file = paste0(mintchip_das_dir, marker, "/", marker, "_consensus_peak_set_FC_0_with_hg38_coordinates.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# Next, we want to find overlap for each 
