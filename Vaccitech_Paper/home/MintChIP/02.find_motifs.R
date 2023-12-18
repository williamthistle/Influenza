# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# Load motifs
data("human_pwms_v2")

binned_motif_analyses <- list()
two_bin_motif_analyses <- list()
two_bin_background_motif_analyses <- list()

for(marker in mintchip_markers) {
  print(marker)
  # Read in background peaks
  background_marker_df <- read.table(paste0(mintchip_das_dir, marker, "/", marker, "_all_peaks.tsv"), sep = "\t")
  colnames(background_marker_df) <- c("seqnames", "start", "end")
  background_marker_df$width <- background_marker_df$end - background_marker_df$start + 1
  # Remove background peaks that aren't length 400
  background_marker_df <- background_marker_df[background_marker_df$width == 401,]
  # Set up GRanges and DNAStringSet objects for background peaks
  background_marker_granges <- makeGRangesFromDataFrame(df = background_marker_df, keep.extra.columns = TRUE)
  background_marker_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, background_marker_granges)
  # Traverse different FC thresholds for marker peaks
  for(fc_threshold in c(0, 0.1, 0.3)) {
    print(fc_threshold)
    list_name <- paste0(marker, "_", fc_threshold)
    list_name_pos <- paste0(list_name, "_", "pos")
    list_name_neg <- paste0(list_name, "_", "neg")
    # Read in current marker peak DF
    testing_marker_df <- read.table(paste0(mintchip_das_dir, marker, "/", marker, "_DESeq2_FC_", fc_threshold, ".tsv"), sep = "\t", header = TRUE)
    # Remove peaks that aren't length 400
    testing_marker_df <- testing_marker_df[testing_marker_df$width == 401,]
    # Row names are set here so they don't conflict with background row names later on
    rownames(testing_marker_df) <- c(1:nrow(testing_marker_df))
    
    # Set up GRanges and DNAStringSet objects for marker peaks
    testing_marker_granges <- makeGRangesFromDataFrame(df = testing_marker_df, keep.extra.columns = TRUE)
    testing_marker_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, testing_marker_granges)
   
    # Set up DF, GRanges, and DNAStringSet objects for negative marker peaks
    testing_marker_negative_df <- testing_marker_df[testing_marker_df$Fold < 0,]
    testing_marker_negative_granges <- makeGRangesFromDataFrame(df = testing_marker_negative_df, keep.extra.columns = TRUE)
    testing_marker_negative_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, testing_marker_negative_granges)
    gc_content_negative <- letterFrequency(testing_marker_negative_seqs, letters = "GC", as.prob = TRUE)
    testing_marker_negative_df$GC.percent <- gc_content_negative[,1]
    testing_marker_negative_granges <- makeGRangesFromDataFrame(df = testing_marker_negative_df, keep.extra.columns = TRUE)
    testing_marker_negative_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, testing_marker_negative_granges)
    
    # Set up DF, GRanges, and DNAStringSet objects for positive marker peaks
    testing_marker_positive_df <- testing_marker_df[testing_marker_df$Fold > 0,]
    testing_marker_positive_granges <- makeGRangesFromDataFrame(df = testing_marker_positive_df, keep.extra.columns = TRUE)
    testing_marker_positive_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, testing_marker_positive_granges)
    gc_content_positive <- letterFrequency(testing_marker_positive_seqs, letters = "GC", as.prob = TRUE)
    testing_marker_positive_df$GC.percent <- gc_content_positive[,1]
    testing_marker_positive_granges <- makeGRangesFromDataFrame(df = testing_marker_positive_df, keep.extra.columns = TRUE)
    testing_marker_positive_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, testing_marker_positive_granges)
    
    # Set up DF, GRanges, and DNAStringSet objects for combined (negative then positive) marker peaks
    # Reorganizing in this way makes it easier to do two bin analysis
    testing_markers_combined_df <- rbind(testing_marker_negative_df, testing_marker_positive_df)
    testing_markers_combined_granges <- c(testing_marker_negative_granges, testing_marker_positive_granges)
    testing_markers_combined_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, testing_markers_combined_granges)
    
    # Set up DF, GRanges, and DNAStringSet objects for background minus marker peaks
    # We use this for Strategy 3 (actually using the background with our marker peaks)
    background_marker_without_candidate_peaks_df <- anti_join(background_marker_df, testing_marker_df, by = c("seqnames", "start", "end"))
    # Candidate and background can't have same row names, so make sure they aren't the same
    marker_offset_index <- nrow(testing_marker_df)
    rownames(background_marker_without_candidate_peaks_df) <- c((marker_offset_index + 1):(marker_offset_index + nrow(background_marker_without_candidate_peaks_df)))
    background_marker_without_candidate_peaks_granges <- makeGRangesFromDataFrame(df = background_marker_without_candidate_peaks_df, keep.extra.columns = TRUE)
    background_marker_without_candidate_peaks_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, background_marker_without_candidate_peaks_granges)
    gc_content_background <- letterFrequency(background_marker_without_candidate_peaks_seqs, letters = "GC", as.prob = TRUE)
    background_marker_without_candidate_peaks_df$GC.percent <- gc_content_background[,1]
    background_marker_without_candidate_peaks_granges <- makeGRangesFromDataFrame(df = background_marker_without_candidate_peaks_df, keep.extra.columns = TRUE)
    background_marker_without_candidate_peaks_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, background_marker_without_candidate_peaks_granges)

    # Grab subset of GC matched negative peaks from background
    neg.peaks.matched <- MatchRegionStats(
      meta.feature = background_marker_without_candidate_peaks_df,
      query.feature = testing_marker_negative_df,
      n = 40000
    )
    
    # Set up negative peak background DF, GRanges, and DNAStringSet objects
    neg_background_marker_without_candidate_peaks_df <- background_marker_without_candidate_peaks_df[neg.peaks.matched,]
    neg_background_marker_without_candidate_peaks_granges <- makeGRangesFromDataFrame(df = neg_background_marker_without_candidate_peaks_df, keep.extra.columns = TRUE)
    neg_background_marker_without_candidate_peaks_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, neg_background_marker_without_candidate_peaks_granges)
    
    # Set up DNAStringSet objects for combined (background then negative) marker peaks
    # Reorganizing in this way makes it easier to do two bin analysis
    neg_markers_combined_seqs <- c(neg_background_marker_without_candidate_peaks_seqs, testing_marker_negative_seqs)
    
    # Grab subset of GC matched positive peaks from background
    pos.peaks.matched <- MatchRegionStats(
      meta.feature = background_marker_without_candidate_peaks_df,
      query.feature = testing_marker_positive_df,
      n = 40000
    )
    
    # Set up positive peak background DF, GRanges, and DNAStringSet objects
    pos_background_marker_without_candidate_peaks_df <- background_marker_without_candidate_peaks_df[pos.peaks.matched,]
    pos_background_marker_without_candidate_peaks_granges <- makeGRangesFromDataFrame(df = pos_background_marker_without_candidate_peaks_df, keep.extra.columns = TRUE)
    pos_background_marker_without_candidate_peaks_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, pos_background_marker_without_candidate_peaks_granges)
    
    # Set up DNAStringSet objects for combined (background then positive) marker peaks
    # Reorganizing in this way makes it easier to do two bin analysis
    pos_markers_combined_seqs <- c(pos_background_marker_without_candidate_peaks_seqs, testing_marker_positive_seqs)
    
    # NOTE CD16 Mono pattern - what does that mean?
    #hist(testing_marker_df$Fold, 100, col = "gray", main = "",
    #     xlab = "Change of accessibility (post - pre)", ylab = "Number of peaks")
    
    # Strategy 1: perform binned motif analysis
    bins <- bin(x = testing_marker_granges$Fold, binmode = "equalN")
    # plotBinDensity(testing_marker$Fold, bins, legend = "topleft")
    # plotBinDiagnostics(seqs = testing_marker_seqs, bins = bins, aspect = "GCfrac")
    # plotBinDiagnostics(seqs = testing_marker_seqs, bins = bins, aspect = "dinucfreq")
    
    se <- calcBinnedMotifEnrR(seqs = testing_marker_seqs, bins = bins, pwmL = human_pwms_v2 )
    binned_motif_analyses[[list_name]] <- se
    
    #colData(se)
    #assay(se, "negLog10Padj")
    # max(assay(se, "negLog10Padj"), na.rm = TRUE)
    
    #sel <- apply(assay(se, "negLog10Padj"), 1, 
    #             function(x) max(abs(x), 0, na.rm = TRUE)) > 4.0
    
    #seSel <- se[sel, ]
    
    #plotMotifHeatmaps(x = seSel, which.plots = c("log2enr", "negLog10Padj"), 
    #                  width = 2.0, cluster = TRUE, maxEnr = 2, maxSig = 10, 
    #                  show_motif_GC = TRUE)
    
    ### Strategy 2: perform two bin motif analysis (up and down-regulated) ###
    bins2 <- rep(c("down", "up"), c(length(testing_marker_negative_granges), length(testing_marker_positive_granges)))
    bins2 <- factor(bins2)
    
    se <- calcBinnedMotifEnrR(seqs = testing_markers_combined_seqs, bins = bins2,
                               pwmL = human_pwms_v2)
    two_bin_motif_analyses[[list_name]] <- se
    
    #sel2 <- apply(assay(se2, "negLog10Padj"), 1, 
    #              function(x) max(abs(x), 0, na.rm = TRUE)) > -log(0.05)
    
    #seSel2 <- se2[sel2, ]
    
    #plotMotifHeatmaps(x = se2[sel2,], which.plots = c("log2enr", "negLog10Padj"), 
    #                  width = 1.8, cluster = TRUE, maxEnr = 2, maxSig = 10,
    #                  show_seqlogo = TRUE)
    
    # Strategy 3 - use GC matched background from histone marker
    
    
    bins3 <- rep(c("background", "down"), c(length(neg_background_marker_without_candidate_peaks_granges), length(testing_marker_negative_granges)))
    bins3 <- factor(bins3)
    
    se <- calcBinnedMotifEnrR(seqs = neg_markers_combined_seqs, bins = bins3,
                               pwmL = human_pwms_v2)
    two_bin_background_motif_analyses[[list_name_neg]] <- se
    
    bins4 <- rep(c("background", "down"), c(length(pos_background_marker_without_candidate_peaks_granges), length(testing_marker_positive_granges)))
    bins4 <- factor(bins4)
    
    se <- calcBinnedMotifEnrR(seqs = pos_markers_combined_seqs, bins = bins4,
                              pwmL = human_pwms_v2)
    two_bin_background_motif_analyses[[list_name_pos]] <- se
    
    #assay(se2, "negLog10Padj")
    
    #sel2 <- apply(assay(se2, "negLog10Padj"), 1, 
    #              function(x) max(abs(x), 0, na.rm = TRUE)) > -log(0.05)
  }
}


























# EXTRA JUNK


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
