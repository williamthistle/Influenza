# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# Load motifs
data("human_pwms_v2")

binned_motif_analyses <- list()
two_bin_motif_analyses <- list()
two_bin_background_motif_analyses <- list()

for(cell_type in snME_dms_cell_types) {
  print(cell_type)
  if(cell_type == "B-Naive" | cell_type == "B-Mem") {
    background_cell_type <- "B"
  } else if(cell_type == "NK-cell2") {
    background_cell_type <- "NK-cell"
  } else {
    background_cell_type <- cell_type
  }
  print(paste0("Background cell type: ", background_cell_type))
  # Read in background peaks
  background_marker_df <- read.table("C:/Users/willi/Desktop/Monocyte_merged_subset_filtered.tsv", sep = "\t")
  colnames(background_marker_df) <- c("seqnames", "start", "strand", "context", "reads", "coverage", "test")
  background_marker_df$end <- background_marker_df$start
  background_marker_df$start <- background_marker_df$start - 100
  background_marker_df$end <- background_marker_df$start + 200
  background_marker_df$width <- background_marker_df$end - background_marker_df$start + 1
  # Set up DF for query peaks
  snME_dms_query_df <- snME_dms[snME_dms$celltype == cell_type,]
  snME_dms_query_df <- snME_dms_query_df[,1:6]
  colnames(snME_dms_query_df) <- c("celltype", "seqnames", "start", "end", "methylation", "count")
  snME_dms_query_df$start <- snME_dms_query_df$start - 100
  snME_dms_query_df$end <- snME_dms_query_df$end + 100
  snME_dms_query_df$width <- snME_dms_query_df$end - snME_dms_query_df$start + 1
  # Row names are set here so they don't conflict with background row names later on
  rownames(snME_dms_query_df) <- c(1:nrow(snME_dms_query_df))
  # Set up GRanges and DNAStringSet objects for query peaks
  snME_dms_query_granges <- makeGRangesFromDataFrame(df = snME_dms_query_df, keep.extra.columns = TRUE)
  snME_dms_query_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, snME_dms_query_granges)
  # Set up DF, GRanges, and DNAStringSet objects for negative marker peaks
  snME_dms_query_negative_df <- snME_dms_query_df[snME_dms_query_df$methylation == "downregulated",]
  snME_dms_query_negative_granges <- makeGRangesFromDataFrame(df = snME_dms_query_negative_df, keep.extra.columns = TRUE)
  snME_dms_query_negative_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, snME_dms_query_negative_granges)
  gc_content_negative <- letterFrequency(snME_dms_query_negative_seqs, letters = "GC", as.prob = TRUE)
  snME_dms_query_negative_df$GC.percent <- gc_content_negative[,1]
  snME_dms_query_negative_granges <- makeGRangesFromDataFrame(df = snME_dms_query_negative_df, keep.extra.columns = TRUE)
  snME_dms_query_negative_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, snME_dms_query_negative_granges)
    
  # Set up DF, GRanges, and DNAStringSet objects for positive marker peaks
  snME_dms_query_positive_df <- snME_dms_query_df[snME_dms_query_df$methylation == "upregulated",]
  snME_dms_query_positive_granges <- makeGRangesFromDataFrame(df = snME_dms_query_positive_df, keep.extra.columns = TRUE)
  snME_dms_query_positive_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, snME_dms_query_positive_granges)
  gc_content_positive <- letterFrequency(snME_dms_query_positive_seqs, letters = "GC", as.prob = TRUE)
  snME_dms_query_positive_df$GC.percent <- gc_content_positive[,1]
  snME_dms_query_positive_granges <- makeGRangesFromDataFrame(df = snME_dms_query_positive_df, keep.extra.columns = TRUE)
  snME_dms_query_positive_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, snME_dms_query_positive_granges)
    
  # Set up DF, GRanges, and DNAStringSet objects for combined (negative then positive) marker peaks
  # Reorganizing in this way makes it easier to do two bin analysis
  snME_dms_querys_combined_df <- rbind(snME_dms_query_negative_df, snME_dms_query_positive_df)
  snME_dms_querys_combined_granges <- c(snME_dms_query_negative_granges, snME_dms_query_positive_granges)
  snME_dms_querys_combined_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, snME_dms_querys_combined_granges)
    
  # Set up DF, GRanges, and DNAStringSet objects for background minus marker peaks
  # We use this for Strategy 3 (actually using the background with our marker peaks)
  background_marker_without_candidate_peaks_df <- anti_join(background_marker_df, snME_dms_query_df, by = c("seqnames", "start", "end"))
  # Candidate and background can't have same row names, so make sure they aren't the same
  marker_offset_index <- nrow(snME_dms_query_df)
  rownames(background_marker_without_candidate_peaks_df) <- c((marker_offset_index + 1):(marker_offset_index + nrow(background_marker_without_candidate_peaks_df)))
  # Subset 500000 peaks (don't need to sample from more than that)
  background_marker_without_candidate_peaks_df <- background_marker_without_candidate_peaks_df[sample(nrow(background_marker_without_candidate_peaks_df), 500000), ]
  background_marker_without_candidate_peaks_granges <- makeGRangesFromDataFrame(df = background_marker_without_candidate_peaks_df, keep.extra.columns = TRUE)
  background_marker_without_candidate_peaks_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, background_marker_without_candidate_peaks_granges)
  gc_content_background <- letterFrequency(background_marker_without_candidate_peaks_seqs, letters = "GC", as.prob = TRUE)
  background_marker_without_candidate_peaks_df$GC.percent <- gc_content_background[,1]
  
  # Grab subset of GC matched negative peaks from background
  neg.peaks.matched <- Signac::MatchRegionStats(
    meta.feature = background_marker_without_candidate_peaks_df,
    query.feature = snME_dms_query_negative_df,
    n = 40000
  )
    
  # Set up negative peak background DF, GRanges, and DNAStringSet objects
  neg_background_marker_without_candidate_peaks_df <- background_marker_without_candidate_peaks_df[neg.peaks.matched,]
  neg_background_marker_without_candidate_peaks_granges <- makeGRangesFromDataFrame(df = neg_background_marker_without_candidate_peaks_df, keep.extra.columns = TRUE)
  neg_background_marker_without_candidate_peaks_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, neg_background_marker_without_candidate_peaks_granges)
    
  # Set up DNAStringSet objects for combined (background then negative) marker peaks
  # Reorganizing in this way makes it easier to do two bin analysis
  neg_markers_combined_seqs <- c(neg_background_marker_without_candidate_peaks_seqs, snME_dms_query_negative_seqs)
    
  # Grab subset of GC matched positive peaks from background
  pos.peaks.matched <- Signac::MatchRegionStats(
    meta.feature = background_marker_without_candidate_peaks_df,
    query.feature = snME_dms_query_positive_df,
    n = 40000
  )
    
  # Set up positive peak background DF, GRanges, and DNAStringSet objects
  pos_background_marker_without_candidate_peaks_df <- background_marker_without_candidate_peaks_df[pos.peaks.matched,]
  pos_background_marker_without_candidate_peaks_granges <- makeGRangesFromDataFrame(df = pos_background_marker_without_candidate_peaks_df, keep.extra.columns = TRUE)
  pos_background_marker_without_candidate_peaks_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, pos_background_marker_without_candidate_peaks_granges)
    
  # Set up DNAStringSet objects for combined (background then positive) marker peaks
  # Reorganizing in this way makes it easier to do two bin analysis
  pos_markers_combined_seqs <- c(pos_background_marker_without_candidate_peaks_seqs, snME_dms_query_positive_seqs)
    
  ### Strategy 1: perform two bin motif analysis (up and down-regulated) ###
  bins2 <- rep(c("down", "up"), c(length(snME_dms_query_negative_granges), length(snME_dms_query_positive_granges)))
  bins2 <- factor(bins2)
    
  se2 <- calcBinnedMotifEnrR(seqs = snME_dms_querys_combined_seqs, bins = bins2,
                              pwmL = human_pwms_v2)
  two_bin_motif_analyses[[list_name]] <- se
    
    sel2 <- apply(assay(se2, "negLog10Padj"), 1, 
                  function(x) max(abs(x), 0, na.rm = TRUE)) > -log(0.05)
    
    #seSel2 <- se2[sel2, ]
    
    #plotMotifHeatmaps(x = se2[sel2,], which.plots = c("log2enr", "negLog10Padj"), 
    #                  width = 1.8, cluster = TRUE, maxEnr = 2, maxSig = 10,
    #                  show_seqlogo = TRUE)
    
    # Strategy 3 - use GC matched background from histone marker
    
    
    bins3 <- rep(c("background", "down"), c(length(neg_background_marker_without_candidate_peaks_granges), length(snME_dms_query_negative_granges)))
    bins3 <- factor(bins3)
    
    se <- calcBinnedMotifEnrR(seqs = neg_markers_combined_seqs, bins = bins3,
                               pwmL = human_pwms_v2)
    two_bin_background_motif_analyses[[list_name_neg]] <- se
    
    bins4 <- rep(c("background", "up"), c(length(pos_background_marker_without_candidate_peaks_granges), length(snME_dms_query_positive_granges)))
    bins4 <- factor(bins4)
    
    se <- calcBinnedMotifEnrR(seqs = pos_markers_combined_seqs, bins = bins4,
                              pwmL = human_pwms_v2)
    two_bin_background_motif_analyses[[list_name_pos]] <- se
    
    se3 <- calcBinnedMotifEnrR(seqs = snME_dms_query_negative_seqs,
                               pwmL = human_pwms_v2,
                               background = "genome",
                               genome = BSgenome.Hsapiens.UCSC.hg38,
                               genome.regions = NULL, # sample from full genome
                               genome.oversample = 2, 
                               verbose = TRUE)
    sel3 <- assay(se3, "negLog10P")[, 1] > 2
    
    
    #assay(se2, "negLog10Padj")
    
    #sel2 <- apply(assay(se2, "negLog10Padj"), 1, 
    #              function(x) max(abs(x), 0, na.rm = TRUE)) > -log(0.05)
}

#saveRDS(binned_motif_analyses, paste0(mintchip_das_dir, "binned_motif_analyses.RDS" ))
#saveRDS(two_bin_motif_analyses, paste0(mintchip_das_dir, "two_bin_motif_analyses.RDS" ))
#saveRDS(two_bin_background_motif_analyses, paste0(mintchip_das_dir, "two_bin_background_motif_analyses.RDS" ))

binned_motif_analyses <- readRDS(paste0(mintchip_das_dir, "binned_motif_analyses.RDS"))
two_bin_motif_analyses <- readRDS(paste0(mintchip_das_dir, "two_bin_motif_analyses.RDS"))
two_bin_background_motif_analyses <- readRDS(paste0(mintchip_das_dir, "two_bin_background_motif_analyses.RDS"))

for(binned_motif_analysis in binned_motif_analyses) {
  sel <- apply(assay(binned_motif_analysis, "negLog10Padj"), 1, 
                            function(x) max(abs(x), 0, na.rm = TRUE)) > -log(0.05)
               
  seSel <- binned_motif_analysis[sel, ]
  if(nrow(seSel) > 0) {
    print(nrow(seSel))
  }
}


for(binned_motif_analysis in two_bin_motif_analyses) {
  sel <- apply(assay(binned_motif_analysis, "negLog10Padj"), 1, 
               function(x) max(abs(x), 0, na.rm = TRUE)) > -log(0.05)
  
  seSel <- binned_motif_analysis[sel, ]
  if(nrow(seSel) > 0) {
    print(nrow(seSel))
  }
}

for(index in 1:length(two_bin_background_motif_analyses)) {
  binned_motif_analysis <- two_bin_background_motif_analyses[[index]]
  sel <- apply(assay(binned_motif_analysis, "negLog10Padj"), 1, 
               function(x) max(abs(x), 0, na.rm = TRUE)) > -log(0.05)
  
  seSel <- binned_motif_analysis[sel, ]
  if(nrow(seSel) > 0) {
    print(index)
    print(nrow(seSel))
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
