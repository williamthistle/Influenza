# MAP3K11
# CD14_Mono     MAP3K11    chr11 65614221    chr11   65615225 65615725     1.000000
# CD14_Mono     MAP3K11    chr11 65614221    chr11   65706207 65706707     0.961039
# CD14_Mono     MAP3K11    chr11 65614221    chr11   65827775 65828275     0.976024
# CD14_Mono     MAP3K11    chr11 65614221    chr11   66069371 66069871     0.958042
idxPass <- which(hvl_placebo_sc_obj$predicted_celltype_majority_vote %in% "CD14 Mono")
cellsPass <- names(hvl_placebo_sc_obj$orig.ident[idxPass])
hvl_placebo_CD14_mono_sc_obj <- subset(x = hvl_placebo_sc_obj, subset = cell_name %in% cellsPass)

Idents(hvl_placebo_CD14_mono_sc_obj) <- "time_point"

hvl_placebo_CD14_mono_sc_obj$time_point <- gsub("D_minus_1", "Pre-Challenge", hvl_placebo_CD14_mono_sc_obj$time_point)
hvl_placebo_CD14_mono_sc_obj$time_point <- gsub("D28", "28 Days Post-Challenge", hvl_placebo_CD14_mono_sc_obj$time_point)

# We do each subplot with scale.factor = 1e7, and we see what the highest possible ymax is
# Then, we use that for each plot and then combine them together

# Plot gene - original ymax is 410
cov_plot <- CoveragePlot(
  object = hvl_placebo_CD14_mono_sc_obj,
  idents = c('Pre-Challenge', '28 Days Post-Challenge'),
  region = "chr11-65613721-65632852",
  annotation = TRUE,
  peaks = FALSE,
  scale.factor = 1e7,
  ymax = 420
)

ggsave("/home/wat2/gene_plot_1.png", plot = cov_plot)

# Plot first peak - original ymax is 420
cov_plot <- CoveragePlot(
  object = hvl_placebo_CD14_mono_sc_obj,
  idents = c('D_minus_1', 'D28'),
  region = "chr11-65615225-65615725",
  annotation = FALSE,
  peaks = FALSE,
  scale.factor = 1e7,
  ymax = 420
)

ggsave("/home/wat2/peak_plot_1.png", plot = cov_plot)

# Plot second peak - original ymax is 21
cov_plot <- CoveragePlot(
  object = hvl_placebo_CD14_mono_sc_obj,
  idents = c('D_minus_1', 'D28'),
  region = "chr11-65706207-65706707",
  annotation = FALSE,
  peaks = FALSE,
  scale.factor = 1e7,
  ymax = 420
)

ggsave("/home/wat2/peak_plot_2.png", plot = cov_plot)

# Plot third peak - original ymax is 63
cov_plot <- CoveragePlot(
  object = hvl_placebo_CD14_mono_sc_obj,
  idents = c('D_minus_1', 'D28'),
  region = "chr11-65827775-65828275",
  annotation = FALSE,
  peaks = FALSE,
  scale.factor = 1e7,
  ymax = 420
)

ggsave("/home/wat2/peak_plot_3.png", plot = cov_plot)

# Plot fourth peak - original ymax is 9.1
cov_plot <- CoveragePlot(
  object = hvl_placebo_CD14_mono_sc_obj,
  idents = c('D_minus_1', 'D28'),
  region = "chr11-66069371-66069871",
  annotation = FALSE,
  peaks = FALSE,
  scale.factor = 1e7,
  ymax = 420
)

ggsave("/home/wat2/peak_plot_4.png", plot = cov_plot)

# Overall winner max is: 420
