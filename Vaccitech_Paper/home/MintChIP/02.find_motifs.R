# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# Load motifs
data("human_pwms_v2")

for(marker in mintchip_markers) {
  print(marker)
  # Read in background peaks
  background_marker_df <- read.table(paste0(mintchip_das_dir, marker, "/", marker, "_all_peaks.tsv"), sep = "\t")
  colnames(background_marker_df) <- c("seqnames", "start", "end")
  background_marker_df$width <- background_marker_df$end - background_marker_df$start + 1
  # Remove background peaks that aren't length 400
  background_marker_df <- background_marker_df[background_marker_df$width == 401,]
  rownames(background_marker_df) <- paste0("background_peak_", c(1:nrow(background_marker_df)))
  # Set up GRanges and DNAStringSet objects for background peaks
  background_marker_granges <- makeGRangesFromDataFrame(df = background_marker_df, keep.extra.columns = TRUE)
  background_marker_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, background_marker_granges)
  # Traverse different FC thresholds for marker peaks
  for(analysis_type in c("consensus_peak_set", "DESeq2", "edgeR")) {
    for(fc_threshold in c(0, 0.1, 0.2, 0.3, 0.585)) {
      print(fc_threshold)
      list_name <- paste0(marker, "_", fc_threshold)
      list_name_pos <- paste0(list_name, "_", "pos")
      list_name_neg <- paste0(list_name, "_", "neg")
      # Read in current marker peak DF
      marker_file_path <- paste0(mintchip_das_dir, marker, "/", marker, "_", analysis_type, "_FC_", fc_threshold, ".tsv")
      if(file.size(marker_file_path) > 1) {
        testing_marker_df <- read.table(marker_file_path, sep = "\t", header = TRUE)
        # Remove peaks that aren't length 400
        testing_marker_df <- testing_marker_df[testing_marker_df$width == 401,]
        # Row names are set here so they don't conflict with background row names later on
        rownames(testing_marker_df) <- paste0("peak_", c(1:nrow(testing_marker_df)))
        
        # Set up GRanges and DNAStringSet objects for marker peaks
        testing_marker_negative_df <- testing_marker_df[testing_marker_df$Fold < 0,]
        testing_marker_positive_df <- testing_marker_df[testing_marker_df$Fold > 0,]
        
        # Set up DF, GRanges, and DNAStringSet objects for background minus marker peaks
        background_marker_without_candidate_peaks_df <- anti_join(background_marker_df, testing_marker_df, by = c("seqnames", "start", "end"))
        background_marker_without_candidate_peaks_granges <- makeGRangesFromDataFrame(df = background_marker_without_candidate_peaks_df, keep.extra.columns = TRUE)
        background_marker_without_candidate_peaks_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, background_marker_without_candidate_peaks_granges)
        gc_content_background <- letterFrequency(background_marker_without_candidate_peaks_seqs, letters = "GC", as.prob = TRUE)
        background_marker_without_candidate_peaks_df$GC.percent <- gc_content_background[,1] * 100
        background_marker_without_candidate_peaks_granges <- makeGRangesFromDataFrame(df = background_marker_without_candidate_peaks_df, keep.extra.columns = TRUE)
        background_marker_without_candidate_peaks_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, background_marker_without_candidate_peaks_granges)
        
        testing_marker_granges <- makeGRangesFromDataFrame(df = testing_marker_df, keep.extra.columns = TRUE)
        testing_marker_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, testing_marker_granges)
        gc_content <- letterFrequency(testing_marker_seqs, letters = "GC", as.prob = TRUE)
        testing_marker_df$GC.percent <- gc_content[,1] * 100
        testing_marker_granges <- makeGRangesFromDataFrame(df = testing_marker_df, keep.extra.columns = TRUE)
        testing_marker_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, testing_marker_granges)
        
        # Grab subset of GC matched negative peaks from background
        peaks.matched <- Signac::MatchRegionStats(
          meta.feature = background_marker_without_candidate_peaks_df,
          query.feature = testing_marker_df,
          features.match = c("GC.percent"),
          n = 100000
        )
        
        background_marker_without_candidate_peaks_df <- background_marker_without_candidate_peaks_df[peaks.matched, , drop = FALSE]
        background_marker_without_candidate_peaks_granges <- makeGRangesFromDataFrame(df = background_marker_without_candidate_peaks_df, keep.extra.columns = TRUE)
        background_marker_without_candidate_peaks_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, background_marker_without_candidate_peaks_granges)
        
        if(nrow(testing_marker_negative_df) > 0) {
          # Set up DF, GRanges, and DNAStringSet objects for negative marker peaks
          testing_marker_negative_granges <- makeGRangesFromDataFrame(df = testing_marker_negative_df, keep.extra.columns = TRUE)
          testing_marker_negative_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, testing_marker_negative_granges)
          gc_content_negative <- letterFrequency(testing_marker_negative_seqs, letters = "GC", as.prob = TRUE)
          testing_marker_negative_df$GC.percent <- gc_content_negative[,1] * 100
          testing_marker_negative_granges <- makeGRangesFromDataFrame(df = testing_marker_negative_df, keep.extra.columns = TRUE)
          testing_marker_negative_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, testing_marker_negative_granges)
          
          # Negative peaks
          query_motifs <- motifmatchr::matchMotifs(human_pwms_v2, testing_marker_negative_seqs, genome = "hg19") 
          query_motifs_matches <- motifmatchr::motifMatches(query_motifs)
          
          background_motifs <- motifmatchr::matchMotifs(human_pwms_v2, background_marker_without_candidate_peaks_seqs, genome = "hg19") 
          background_motifs_matches <- motifmatchr::motifMatches(background_motifs)
          
          query.counts <- colSums(x = query_motifs_matches)
          background.counts <- colSums(x = background_motifs_matches)
          percent.observed <- query.counts/length(x = testing_marker_negative_seqs) * 100
          percent.background <- background.counts/length(x = background_marker_without_candidate_peaks_seqs) * 
            100
          fold.enrichment <- percent.observed/percent.background
          p.list <- vector(mode = "numeric")
          for (i in seq_along(along.with = query.counts)) {
            p.list[[i]] <- phyper(q = query.counts[[i]] - 1, m = background.counts[[i]], 
                                  n = length(x = background_marker_without_candidate_peaks_seqs) - background.counts[[i]], 
                                  k = length(x = testing_marker_negative_seqs), lower.tail = FALSE)
          }
          
          results_neg <- data.frame(motif = names(x = query.counts), observed = query.counts, 
                                    background = background.counts, percent.observed = percent.observed, 
                                    percent.background = percent.background, fold.enrichment = fold.enrichment, 
                                    pvalue = p.list, motif.name = name(human_pwms_v2), 
                                    p.adjust = p.adjust(p = p.list, method = "BH"), 
                                    stringsAsFactors = FALSE)
          results_neg <- results_neg[order(results_neg$p.adjust),]
          write.table(results_neg, file = paste0(mintchip_das_dir, marker, "/", marker, "_motifs_", analysis_type, "_", fc_threshold, "_neg.tsv"), sep = "\t", quote = FALSE)
        }
        
        if(nrow(testing_marker_positive_df) > 0) {
          # Set up DF, GRanges, and DNAStringSet objects for positive marker peaks
          testing_marker_positive_granges <- makeGRangesFromDataFrame(df = testing_marker_positive_df, keep.extra.columns = TRUE)
          testing_marker_positive_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, testing_marker_positive_granges)
          gc_content_positive <- letterFrequency(testing_marker_positive_seqs, letters = "GC", as.prob = TRUE)
          testing_marker_positive_df$GC.percent <- gc_content_positive[,1] * 100
          testing_marker_positive_granges <- makeGRangesFromDataFrame(df = testing_marker_positive_df, keep.extra.columns = TRUE)
          testing_marker_positive_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, testing_marker_positive_granges)
          
          # POS
          query_motifs <- motifmatchr::matchMotifs(human_pwms_v2, testing_marker_positive_seqs, genome = "hg19") 
          query_motifs_matches <- motifmatchr::motifMatches(query_motifs)
          
          background_motifs <- motifmatchr::matchMotifs(human_pwms_v2, background_marker_without_candidate_peaks_seqs, genome = "hg19") 
          background_motifs_matches <- motifmatchr::motifMatches(background_motifs)
          
          query.counts <- colSums(x = query_motifs_matches)
          background.counts <- colSums(x = background_motifs_matches)
          percent.observed <- query.counts/length(x = testing_marker_positive_seqs) * 100
          percent.background <- background.counts/length(x = background_marker_without_candidate_peaks_seqs) * 
            100
          fold.enrichment <- percent.observed/percent.background
          p.list <- vector(mode = "numeric")
          for (i in seq_along(along.with = query.counts)) {
            p.list[[i]] <- phyper(q = query.counts[[i]] - 1, m = background.counts[[i]], 
                                  n = length(x = background_marker_without_candidate_peaks_seqs) - background.counts[[i]], 
                                  k = length(x = testing_marker_positive_seqs), lower.tail = FALSE)
          }
          
          results_pos <- data.frame(motif = names(x = query.counts), observed = query.counts, 
                                    background = background.counts, percent.observed = percent.observed, 
                                    percent.background = percent.background, fold.enrichment = fold.enrichment, 
                                    pvalue = p.list, motif.name = name(human_pwms_v2), 
                                    p.adjust = p.adjust(p = p.list, method = "BH"), 
                                    stringsAsFactors = FALSE)
          results_pos <- results_pos[order(results_pos$p.adjust),]
          write.table(results_pos, file = paste0(mintchip_das_dir, marker, "/", marker, "_motifs_", analysis_type, "_", fc_threshold, "_pos.tsv"), sep = "\t", quote = FALSE)
        }
      }
    }
  }
}



# nullranges approach - pretty cool and definitely works, but only matches up to # of query peaks
background_marker_without_candidate_peaks_df
testing_marker_negative_df <- testing_marker_negative_df[,c(1,2,3,4,13)]
background_marker_without_candidate_peaks_df$query <- FALSE
testing_marker_negative_df$query <- TRUE

finding_negative_matching_background_df <- rbind(testing_marker_negative_df, background_marker_without_candidate_peaks_df)
rownames(finding_negative_matching_background_df) <- NULL

testing_marker_negative_granges <- makeGRangesFromDataFrame(df = finding_negative_matching_background_df, keep.extra.columns = TRUE)

mgr <- matchRanges(focal = testing_marker_negative_granges[testing_marker_negative_granges$query],
                   pool = testing_marker_negative_granges[!testing_marker_negative_granges$query],
                   covar = ~GC.percent)
mgr
