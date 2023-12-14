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
background_marker <- read.table(paste0(mintchip_das_dir, "H3K4me1/H3K4me1_all_peaks.tsv"), sep = "\t")
colnames(background_marker) <- c("idx", "seqnames", "start", "end", "strand")
background_marker$width <- background_marker$end - background_marker$start + 1

# Remove peaks that aren't length 400 (why do these exist, anyway?)
testing_marker_df <- testing_marker[testing_marker$width == 401,]
background_marker_df <- background_marker[background_marker$width == 401,]

testing_marker_negative_df <- testing_marker_df[testing_marker_df$Fold < 0,]
testing_marker_positive_df <- testing_marker_df[testing_marker_df$Fold > 0,]

testing_marker <- makeGRangesFromDataFrame(df = testing_marker_df, keep.extra.columns = TRUE)
background_marker <- makeGRangesFromDataFrame(df = background_marker_df, keep.extra.columns = TRUE)

hist(testing_marker$Fold, 100, col = "gray", main = "",
     xlab = "Change of accessibility (post - pre)", ylab = "Number of peaks")

# Perform binned motif analysis
bins <- bin(x = testing_marker$Fold, binmode = "equalN", nElement = 500, minAbsX = 0.8)
plotBinDensity(testing_marker$Fold, bins, legend = "topleft")

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
testing_marker_negative <- testing_marker[testing_marker$Fold <= 0]
testing_marker_positive <- testing_marker[testing_marker$Fold >= 0]

testing_markers_combined <- c(testing_marker_negative, testing_marker_positive)
testing_markers_combined_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, testing_markers_combined)

bins2 <- rep(c("down", "up"), c(length(testing_marker_negative), length(testing_marker_positive)))
bins2 <- factor(bins2)

se2 <- calcBinnedMotifEnrR(seqs = testing_markers_combined_seqs, bins = bins2,
                           pwmL = human_pwms_v2)
sel2 <- apply(assay(se2, "negLog10Padj"), 1, 
             function(x) max(abs(x), 0, na.rm = TRUE)) > -log(0.05)

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













# New strategy!
testing_marker <- read.table(paste0(mintchip_das_dir, "H3K4me1/H3K4me1_consensus_peak_set_FC_0.3.tsv"), sep = "\t", header = TRUE)
background_marker <- read.table(paste0(mintchip_das_dir, "H3K4me1/H3K4me1_all_peaks.tsv"), sep = "\t")
colnames(background_marker) <- c("idx", "seqnames", "start", "end", "strand")
background_marker$width <- background_marker$end - background_marker$start + 1
background_marker$strand <- "*"

background_marker <- anti_join(background_marker, testing_marker, by = c("seqnames", "start", "end"))

rownames(testing_marker) <- c(1:nrow(testing_marker))
marker_offset_index <- nrow(testing_marker)
rownames(background_marker) <- c((marker_offset_index + 1):(marker_offset_index + nrow(background_marker)))

# Remove peaks that aren't length 400 (why do these exist, anyway?)
testing_marker_df <- testing_marker[testing_marker$width == 401,]
background_marker_df <- background_marker[background_marker$width == 401,]

testing_marker_negative_df <- testing_marker_df[testing_marker_df$Fold < 0,]
testing_marker_positive_df <- testing_marker_df[testing_marker_df$Fold > 0,]

testing_marker <- makeGRangesFromDataFrame(df = testing_marker_df, keep.extra.columns = TRUE)
background_marker <- makeGRangesFromDataFrame(df = background_marker_df, keep.extra.columns = TRUE)

testing_marker_negative <- testing_marker[testing_marker$Fold < 0]
testing_marker_negative_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, testing_marker_negative)
testing_marker_positive <- testing_marker[testing_marker$Fold > 0]
testing_marker_positive_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, testing_marker_positive)
background_marker_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, background_marker)

gc_content_negative <- letterFrequency(testing_marker_negative_seqs, letters = "GC", as.prob = TRUE)
gc_content_positive <- letterFrequency(testing_marker_positive_seqs, letters = "GC", as.prob = TRUE)
gc_content_background <- letterFrequency(background_marker_seqs, letters = "GC", as.prob = TRUE)

testing_marker_negative_df$GC.percent <- gc_content_negative[,1]
testing_marker_positive_df$GC.percent <- gc_content_positive[,1]
background_marker_df$GC.percent <- gc_content_background[,1]

testing_marker_negative <- makeGRangesFromDataFrame(df = testing_marker_negative_df, keep.extra.columns = TRUE)
testing_marker_positive <- makeGRangesFromDataFrame(df = testing_marker_positive_df, keep.extra.columns = TRUE)

neg.peaks.matched <- MatchRegionStats(
  meta.feature = background_marker_df,
  query.feature = testing_marker_negative_df,
  n = 40000
)

neg_background_df <- background_marker_df[neg.peaks.matched,]
neg_background <- makeGRangesFromDataFrame(df = neg_background_df, keep.extra.columns = TRUE)
neg_background_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, neg_background)

# combine the two sets or genomic regions
neg_markers_combined <- c(neg_background, testing_marker_negative)

# extract sequences from the genome
neg_markers_combined_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, neg_markers_combined)

bins2 <- rep(c("background", "down"), c(length(neg_background_seqs), length(testing_marker_negative_seqs)))
bins2 <- factor(bins2)
table(bins2)

se2 <- calcBinnedMotifEnrR(seqs = neg_markers_combined_seq, bins = bins2,
                           pwmL = human_pwms_v2)

assay(se2, "negLog10Padj")

sel2 <- apply(assay(se2, "negLog10Padj"), 1, 
              function(x) max(abs(x), 0, na.rm = TRUE)) > -log(0.05)
