### MAP3K11 in CD14 Monocytes ###

# CD14_Mono     MAP3K11    chr11 65614221    chr11   65615225 65615725     1.000000
# CD14_Mono     MAP3K11    chr11 65614221    chr11   65706207 65706707     0.961039
# CD14_Mono     MAP3K11    chr11 65614221    chr11   65827775 65828275     0.976024
# CD14_Mono     MAP3K11    chr11 65614221    chr11   66069371 66069871     0.958042
# Downregulation of H3K27ac: 	chr11:65612955-65613355
# Hypermethylation: chr11:65613039-65613039
idxPass <- which(hvl_placebo_sc_obj$predicted_celltype_majority_vote %in% "CD14 Mono")
cellsPass <- names(hvl_placebo_sc_obj$orig.ident[idxPass])
hvl_placebo_CD14_mono_sc_obj <- subset(x = hvl_placebo_sc_obj, subset = cell_name %in% cellsPass)

hvl_placebo_CD14_mono_sc_obj$time_point <- gsub("D_minus_1", "Pre-Challenge", hvl_placebo_CD14_mono_sc_obj$time_point)
hvl_placebo_CD14_mono_sc_obj$time_point <- gsub("D28", "28 Days Post-Challenge", hvl_placebo_CD14_mono_sc_obj$time_point)

Idents(hvl_placebo_CD14_mono_sc_obj) <- "time_point"
levels(hvl_placebo_CD14_mono_sc_obj) <- c("Pre-Challenge", "28 Days Post-Challenge")

# Plot everything
#cov_plot <- CoveragePlot(
#  object = hvl_placebo_CD14_mono_sc_obj,
#  region = c("chr11-65611955-65615725", "chr11-65706207-65706707", "chr11-65827775-65828275", "chr11-66069371-66069871"),
#  annotation = TRUE,
#  peaks = FALSE
#)

# Alternative plot of everything
cov_plot <- CoveragePlot(
  object = hvl_placebo_CD14_mono_sc_obj,
  region = c("chr11-65611955-66069871"),
  annotation = TRUE,
  peaks = FALSE
)


# cov_plot_build <- ggplot_build(cov_plot[[1]][[1]])

# cov_plot <- cov_plot & scale_fill_manual(values = c("grey", "grey"))

ggsave("/home/wat2/everything_plot.png", plot = cov_plot, width = 8)

# Highlight regions of interest so I can annotate them

# Acetylation (and hypermethylation)
# First site
# Second site 
# Third site
# Fourth site
ranges.show <- StringToGRanges(c("chr11-65612955-65613355", "chr11-65615225-65615725", "chr11-65706207-65706707", "chr11-65827775-65828275", "chr11-66069371-66069871"))
ranges.show$color <- "orange"

cov_plot <- CoveragePlot(
  object = hvl_placebo_CD14_mono_sc_obj,
  region = c("chr11-65611955-66069871"),
  annotation = TRUE,
  peaks = FALSE,
  region.highlight = ranges.show
)

# cov_plot <- cov_plot & scale_fill_manual(values = c("grey", "grey"))

ggsave("/home/wat2/highlighted_plot.png", plot = cov_plot, width = 8)

### FOSB in CD16 Monocytes ###
# CD16_Mono        FOSB    chr19 45467995    chr19   45584223 45584723     0.999001
# Hypermethylation: chr19-45464252-45464252
# Transcript Location: 45467996-45475179
idxPass <- which(hvl_placebo_sc_obj$predicted_celltype_majority_vote %in% "CD16 Mono")
cellsPass <- names(hvl_placebo_sc_obj$orig.ident[idxPass])
hvl_placebo_CD16_mono_sc_obj <- subset(x = hvl_placebo_sc_obj, subset = cell_name %in% cellsPass)

hvl_placebo_CD16_mono_sc_obj$time_point <- gsub("D_minus_1", "Pre-Challenge", hvl_placebo_CD16_mono_sc_obj$time_point)
hvl_placebo_CD16_mono_sc_obj$time_point <- gsub("D28", "28 Days Post-Challenge", hvl_placebo_CD16_mono_sc_obj$time_point)

Idents(hvl_placebo_CD16_mono_sc_obj) <- "time_point"
levels(hvl_placebo_CD16_mono_sc_obj) <- c("Pre-Challenge", "28 Days Post-Challenge")

# Plot everything
cov_plot <- CoveragePlot(
  object = hvl_placebo_CD16_mono_sc_obj,
  region = c("chr19-45464252-45475179", "chr19-45584223-45584723"),
  annotation = TRUE,
  peaks = FALSE
)

cov_plot <- cov_plot & scale_fill_manual(values = c("#00bfc4", "grey"))

ggsave("/home/wat2/everything_plot_FOSB.png", plot = cov_plot, width = 20)

# Highlight regions of interest so I can annotate them
ranges.show <- StringToGRanges(c("chr19-45464152-45464352"))
ranges.show$color <- "orange"

# Plot everything
cov_plot <- CoveragePlot(
  object = hvl_placebo_CD16_mono_sc_obj,
  region = c("chr19-45464252-45475179", "chr19-45584223-45584723"),
  annotation = TRUE,
  peaks = FALSE,
  region.highlight = ranges.show
)

ggsave("/home/wat2/highlighted_plot_FOSB.png", plot = cov_plot, width = 20)

### IL1RAP in CD14 Monocytes ###
# CD14_Mono      IL1RAP     chr3 190514084     chr3  190586645 190587145     0.984016
# Downregulation of H3K36me3: chr3:190634549-190634949
# Transcript Location: chr3-190514085-190651514
idxPass <- which(hvl_placebo_sc_obj$predicted_celltype_majority_vote %in% "CD14 Mono")
cellsPass <- names(hvl_placebo_sc_obj$orig.ident[idxPass])
hvl_placebo_CD14_mono_sc_obj <- subset(x = hvl_placebo_sc_obj, subset = cell_name %in% cellsPass)

hvl_placebo_CD14_mono_sc_obj$time_point <- gsub("D_minus_1", "Pre-Challenge", hvl_placebo_CD14_mono_sc_obj$time_point)
hvl_placebo_CD14_mono_sc_obj$time_point <- gsub("D28", "28 Days Post-Challenge", hvl_placebo_CD14_mono_sc_obj$time_point)

Idents(hvl_placebo_CD14_mono_sc_obj) <- "time_point"
levels(hvl_placebo_CD14_mono_sc_obj) <- c("Pre-Challenge", "28 Days Post-Challenge")

# Plot everything
cov_plot <- CoveragePlot(
  object = hvl_placebo_CD14_mono_sc_obj,
  region = c("chr3-190514085-190651514"),
  annotation = TRUE,
  peaks = FALSE
)

cov_plot <- cov_plot & scale_fill_manual(values = c("#00bfc4", "grey"))

ggsave("/home/wat2/everything_plot_IL1RAP.png", plot = cov_plot, width = 20)

# Highlight regions of interest so I can annotate them
ranges.show <- StringToGRanges(c("chr3-190634549-190634949", "chr3-190586645-190587145"))
ranges.show$color <- "orange"

# Plot everything
cov_plot <- CoveragePlot(
  object = hvl_placebo_CD14_mono_sc_obj,
  region = c("chr3-190514085-190651514"),
  annotation = TRUE,
  peaks = FALSE,
  region.highlight = ranges.show
)

ggsave("/home/wat2/highlighted_plot_IL1RAP.png", plot = cov_plot, width = 20)

### SETD2 in CD14 Monocytes ###
# CD14_Mono       SETD2     chr3 47164113     chr3   46949083 46949583     0.952048
# CD14_Mono       SETD2     chr3 47164113     chr3   47030952 47031452     0.959041
# CD14_Mono       SETD2     chr3 47164113     chr3   47163885 47164385     0.991009
# CD14_Mono       SETD2     chr3 47164113     chr3   47282071 47282571     0.960040
# CD14_Mono       SETD2     chr3 47164113     chr3   47282636 47283136     0.978022
# Transcript Location: chr3-47016428-47164113
idxPass <- which(hvl_placebo_sc_obj$predicted_celltype_majority_vote %in% "CD14 Mono")
cellsPass <- names(hvl_placebo_sc_obj$orig.ident[idxPass])
hvl_placebo_CD14_mono_sc_obj <- subset(x = hvl_placebo_sc_obj, subset = cell_name %in% cellsPass)

hvl_placebo_CD14_mono_sc_obj$time_point <- gsub("D_minus_1", "Pre-Challenge", hvl_placebo_CD14_mono_sc_obj$time_point)
hvl_placebo_CD14_mono_sc_obj$time_point <- gsub("D28", "28 Days Post-Challenge", hvl_placebo_CD14_mono_sc_obj$time_point)

Idents(hvl_placebo_CD14_mono_sc_obj) <- "time_point"
levels(hvl_placebo_CD14_mono_sc_obj) <- c("Pre-Challenge", "28 Days Post-Challenge")

# Plot everything
cov_plot <- CoveragePlot(
  object = hvl_placebo_CD14_mono_sc_obj,
  region = c("chr3-46949083-46949583", "chr3-47030952-47031452", "chr3-47163885-47164385",
             "chr3-47282071-47282571", "chr3-47282636-47283136"),
  annotation = TRUE,
  peaks = FALSE
)

cov_plot <- cov_plot & scale_fill_manual(values = c("#00bfc4", "grey"))

ggsave("/home/wat2/everything_plot_SETD2_CD14_Mono.png", plot = cov_plot, width = 20)

### SETD2 in cDCs ###
# cDC       SETD2     chr3 47164113     chr3   47381267 47381767     0.975025
# Transcript Location: chr3-47016428-47164113
idxPass <- which(hvl_placebo_sc_obj$predicted_celltype_majority_vote %in% "cDC")
cellsPass <- names(hvl_placebo_sc_obj$orig.ident[idxPass])
hvl_placebo_cDC_sc_obj <- subset(x = hvl_placebo_sc_obj, subset = cell_name %in% cellsPass)

hvl_placebo_cDC_sc_obj$time_point <- gsub("D_minus_1", "Pre-Challenge", hvl_placebo_cDC_sc_obj$time_point)
hvl_placebo_cDC_sc_obj$time_point <- gsub("D28", "28 Days Post-Challenge", hvl_placebo_cDC_sc_obj$time_point)

Idents(hvl_placebo_cDC_sc_obj) <- "time_point"
levels(hvl_placebo_cDC_sc_obj) <- c("Pre-Challenge", "28 Days Post-Challenge")

# Plot everything
cov_plot <- CoveragePlot(
  object = hvl_placebo_cDC_sc_obj,
  region = c("chr3-47381267-47381767"),
  annotation = TRUE,
  peaks = FALSE
)

cov_plot <- cov_plot & scale_fill_manual(values = c("#00bfc4", "grey"))

ggsave("/home/wat2/everything_plot_SETD2_cDC.png", plot = cov_plot, width = 20)

### ASH1L in cDCs ###
#       cDC       ASH1L     chr1 155562803     chr1  155934239 155934739     0.964036
#       cDC       ASH1L     chr1 155562803     chr1  155944375 155944875     0.957043
#       cDC       ASH1L     chr1 155562803     chr1  156504834 156505334     0.972028
# Transcript Location: chr1-155335268-155562803
idxPass <- which(hvl_placebo_sc_obj$predicted_celltype_majority_vote %in% "cDC")
cellsPass <- names(hvl_placebo_sc_obj$orig.ident[idxPass])
hvl_placebo_cDC_sc_obj <- subset(x = hvl_placebo_sc_obj, subset = cell_name %in% cellsPass)

hvl_placebo_cDC_sc_obj$time_point <- gsub("D_minus_1", "Pre-Challenge", hvl_placebo_cDC_sc_obj$time_point)
hvl_placebo_cDC_sc_obj$time_point <- gsub("D28", "28 Days Post-Challenge", hvl_placebo_cDC_sc_obj$time_point)

Idents(hvl_placebo_cDC_sc_obj) <- "time_point"
levels(hvl_placebo_cDC_sc_obj) <- c("Pre-Challenge", "28 Days Post-Challenge")

# Plot everything
cov_plot <- CoveragePlot(
  object = hvl_placebo_cDC_sc_obj,
  region = c("chr1-155934239-155934739", "chr1-155944375-155944875", "chr1-156504834-156505334"),
  annotation = TRUE,
  peaks = FALSE
)

cov_plot <- cov_plot & scale_fill_manual(values = c("#00bfc4", "grey"))

ggsave("/home/wat2/everything_plot_ASH1L_cDC.png", plot = cov_plot, width = 20)



