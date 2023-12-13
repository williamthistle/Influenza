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
testing_marker <- read.table(paste0(sc_das_dir, "diff_peaks/D28-vs-D_minus_1-degs-CD16_Mono-time_point-controlling_for_subject_id_sc_filtered_pct_0.01.tsv"), sep = "\t", header = TRUE)
testing_marker <- testing_marker[abs(testing_marker$avg_log2FC) >= 1,]
chromosomes <- sapply(strsplit(rownames(testing_marker), "-"), `[`, 1)
start_coords <- as.numeric(sapply(strsplit(rownames(testing_marker), "-"), `[`, 2))
end_coords <- as.numeric(sapply(strsplit(rownames(testing_marker), "-"), `[`, 3))
testing_marker$seqnames <- chromosomes
testing_marker$start <- start_coords
testing_marker$end <- end_coords
testing_marker$width <- testing_marker$end - testing_marker$start 

# Remove peaks that aren't length 400 (why do these exist, anyway?)
#testing_marker <- testing_marker[testing_marker$width == 401,]

testing_marker <- makeGRangesFromDataFrame(df = testing_marker, keep.extra.columns = TRUE)

hist(testing_marker$avg_log2FC, 50, col = "gray", main = "",
     xlab = "Change of accessibility (post - pre)", ylab = "Number of peaks")

# Perform binned motif analysis
bins <- bin(x = testing_marker$avg_log2FC, binmode = "equalN", nElement = 500)
plotBinDensity(testing_marker$avg_log2FC, bins, legend = "topleft")

testing_marker_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, testing_marker)

plotBinDiagnostics(seqs = testing_marker_seqs, bins = bins, aspect = "GCfrac")
plotBinDiagnostics(seqs = testing_marker_seqs, bins = bins, aspect = "dinucfreq")

se <- calcBinnedMotifEnrR(seqs = testing_marker_seqs, bins = bins, pwmL = human_pwms_v2 )
#colData(se)
#assay(se, "negLog10Padj")
# max(assay(se, "negLog10Padj"), na.rm = TRUE)

sel <- apply(assay(se, "negLog10Padj"), 1, 
             function(x) max(abs(x), 0, na.rm = TRUE)) > 4.0

seSel <- se[sel, ]

plotMotifHeatmaps(x = seSel, which.plots = c("log2enr", "negLog10Padj"), 
                  width = 2.0, cluster = TRUE, maxEnr = 2, maxSig = 10, 
                  show_motif_GC = TRUE)

### Perform two bin motif analysis (up and down-regulated) ###
testing_marker_negative <- testing_marker[testing_marker$avg_log2FC < 1]
testing_marker_positive <- testing_marker[testing_marker$avg_log2FC > 1]

testing_markers_combined <- c(testing_marker_negative, testing_marker_positive)
testing_markers_combined_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, testing_markers_combined)

bins2 <- rep(c("down", "up"), c(length(testing_marker_negative), length(testing_marker_positive)))
bins2 <- factor(bins2)

se2 <- calcBinnedMotifEnrR(seqs = testing_markers_combined_seqs, bins = bins2,
                           pwmL = human_pwms_v2)
sel2 <- apply(assay(se2, "negLog10Padj"), 1, 
             function(x) max(abs(x), 0, na.rm = TRUE)) > -log10(0.05)

seSel2 <- se2[sel2, ]

plotMotifHeatmaps(x = se2[sel2,], which.plots = c("log2enr", "negLog10Padj"), 
                  width = 1.8, cluster = TRUE, maxEnr = 2, maxSig = 10,
                  show_seqlogo = TRUE)

# Test against genome background
testing_marker_negative_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, testing_marker_negative)
testing_marker_positive_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, testing_marker_positive)

# Negative
se3 <- calcBinnedMotifEnrR(seqs = testing_marker_negative_seqs,
                           pwmL = human_pwms_v2,
                           background = "genome",
                           genome = BSgenome.Hsapiens.UCSC.hg19,
                           genome.regions = NULL, # sample from full genome
                           genome.oversample = 2, 
                           BPPARAM = BiocParallel::SerialParam(RNGseed = 42),
                           verbose = TRUE)

assay(se3, "negLog10P")

# Positive
se4 <- calcBinnedMotifEnrR(seqs = testing_marker_positive_seqs,
                           pwmL = human_pwms_v2,
                           background = "genome",
                           genome = BSgenome.Hsapiens.UCSC.hg19,
                           genome.regions = NULL, # sample from full genome
                           genome.oversample = 2, 
                           BPPARAM = BiocParallel::SerialParam(RNGseed = 42),
                           verbose = TRUE)

assay(se4, "negLog10P")

test <- assay(se4, "negLog10Padj")[, 1]
test <- test > 4
sum(test, na.rm = TRUE)