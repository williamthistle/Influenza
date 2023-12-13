library(GenomicRanges)
library(SummarizedExperiment)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg19)
library(monaLisa)
library(ComplexHeatmap)
library(circlize)
library(chromVARmotifs)

# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# Load motifs
data("human_pwms_v2")

set.seed(SPEEDI::get_speedi_seed())
testing_marker <- read.table(paste0(mintchip_das_dir, "H3K4me1/H3K4me1_consensus_peak_set_FC_0.1.tsv"), sep = "\t", header = TRUE)
testing_marker <- makeGRangesFromDataFrame(df = testing_marker, keep.extra.columns = TRUE)

hist(testing_marker$Fold, 100, col = "gray", main = "",
     xlab = "Change of accessibility (post - pre)", ylab = "Number of peaks")


bins <- bin(x = testing_marker$Fold, binmode = "equalN", nElement = 500, minAbsX = 0.5)
plotBinDensity(testing_marker$Fold, bins, legend = "topleft")

testing_marker_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, testing_marker)