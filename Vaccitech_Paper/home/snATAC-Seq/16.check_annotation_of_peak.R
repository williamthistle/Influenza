# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# Peak to check
# chr12:121469731-121470131: H3K4me1
# chr12:121491198-121491598: H3K4me3
# chr12:121464402-121464802: H3K36me3
# chr12:121461800-121461800: hypermethylation
# chr12:121436478-121436878: H3K9me3

checking_peak <- data.frame(seqnames = "chr12", start = "121436478", end = "121436878")
checking_peak_annotated <- annotatePeak(makeGRangesFromDataFrame(checking_peak), TxDb = txdb, annoDb = "org.Hs.eg.db")
checking_peak_annotated <- as.data.frame(checking_peak_annotated)
checking_peak_annotated
checking_peak_annotated <- annotatePeak(makeGRangesFromDataFrame(checking_peak), TxDb = txdb, annoDb = "org.Hs.eg.db",
                                        genomicAnnotationPriority = "Intron")
checking_peak_annotated <- as.data.frame(checking_peak_annotated)
checking_peak_annotated