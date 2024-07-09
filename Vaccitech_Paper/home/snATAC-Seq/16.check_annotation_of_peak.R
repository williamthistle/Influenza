# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# Peak to check
checking_peak <- data.frame(seqnames = "chr11", start = "66069371", end = "66069871")
checking_peak_annotated <- annotatePeak(makeGRangesFromDataFrame(checking_peak), TxDb = txdb, annoDb = "org.Hs.eg.db")
checking_peak_annotated <- as.data.frame(checking_peak_annotated)
checking_peak_annotated
