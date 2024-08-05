# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# Peak to check
# KDM2B
# First peak: chr12:121546911-121547411 (Intron 6 for standard transcript ENST00000377071, promoter for ENST00000538503)
# Second peak: chr12:121580095-121580595 (Promoter for standard transcript ENST00000377071, promoter for ENST00000377069)
# H3K4me1 negative: chr12:121469731-121470131 (Intron 12 for standard transcript ENST00000377071, promoter for ENST00000538243) - 110.9kb
# H3K4me3 negative: chr12:121491198-121491598 (Intron 12 for standard transcript ENST00000377071) - 89.4 kb
# H3K9me3 negative: chr12:121436478-121436878 (Intron 22 for standard transcript ENST00000377071) - 144.1 kb
# H3K36me3 positive: chr12:121464402-121464802 (Intron 12 for standard transcript ENST00000377071) - 116.2 kb
# snME positive: chr12:121461800-121461800 (Intron 12 for standard transcript ENST00000377071) - 119.2 kb
# Gene TSS: 121581023

checking_peak <- data.frame(seqnames = "chr12", start = "121436478", end = "121436878")
checking_peak_annotated <- annotatePeak(makeGRangesFromDataFrame(checking_peak), TxDb = txdb, annoDb = "org.Hs.eg.db")
checking_peak_annotated <- as.data.frame(checking_peak_annotated)
checking_peak_annotated
checking_peak_annotated <- annotatePeak(makeGRangesFromDataFrame(checking_peak), TxDb = txdb, annoDb = "org.Hs.eg.db",
                                        genomicAnnotationPriority = "Intron")
checking_peak_annotated <- as.data.frame(checking_peak_annotated)
checking_peak_annotated