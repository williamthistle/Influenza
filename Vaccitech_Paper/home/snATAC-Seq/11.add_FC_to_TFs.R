# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

snATAC_cell_types <- c("B", "CD4 Memory", "CD8 Memory", "CD14 Mono", "CD16 Mono", "NK", "CD4 Naive", "CD8 Naive", "cDC", "MAIT", "Proliferating", "pDC")

# motif_dirs <- c(scATAC_hvl_placebo_das_motif_dir, scATAC_hvl_vaccinated_das_motif_dir, scATAC_lvl_placebo_das_motif_dir)
motif_dirs <- c(scATAC_lvl_placebo_das_motif_dir)

# Step 1: Add FC to motif files
for(motif_dir in motif_dirs) {
  for(snATAC_cell_type in snATAC_cell_types) {
    snATAC_cell_type_for_file_name <- sub(" ", "_", snATAC_cell_type)
    for(analysis_type in c("sc", "final")) {
      for(pct in c(0.01, 0.05, 0.1)) {
        current_motif_dir <- paste0(motif_dir, snATAC_cell_type_for_file_name, "/", analysis_type, "/", pct, "/with_bg/")
        if(dir.exists(current_motif_dir)) {
          current_output_dir <- paste0(current_motif_dir, "with_fc_added/")
          if (!dir.exists(current_output_dir)) {dir.create(current_output_dir)}
          for(fc in c(0.1, 0.3, 0.585, 1, 2)) {
            # Grab associated motif files using pattern
            upregulated_motifs_file <- list.files(path = current_motif_dir, pattern = paste0("_FC_", fc, "_"), full.names = TRUE)
            downregulated_motifs_file <- list.files(path = current_motif_dir, pattern = paste0("_FC_-", fc, "_"), full.names = TRUE)
            if(length(upregulated_motifs_file) == 1 && length(downregulated_motifs_file) == 1) {
              upregulated_motifs <- read.table(upregulated_motifs_file, sep = "\t", header = TRUE)
              downregulated_motifs <- read.table(downregulated_motifs_file, sep = "\t", header = TRUE)
              # Add FC values for upregulated motifs
              motif_fc_values <- c()
              for(current_row in 1:nrow(upregulated_motifs)) {
                current_motif <- upregulated_motifs[current_row,]
                current_motif_percent_observed <- current_motif$percent.observed
                alternative_motif <- downregulated_motifs[downregulated_motifs$motif %in% current_motif$motif,]
                alternative_motif_percent_observed <- alternative_motif$percent.observed
                motif_fc_value <- log(current_motif_percent_observed / alternative_motif_percent_observed)
                motif_fc_values <- c(motif_fc_values, motif_fc_value)
              }
              upregulated_motifs$fc_value <- motif_fc_values
              upregulated_motifs_sorted_by_fc <- upregulated_motifs
              upregulated_motifs_sorted_by_fc <- upregulated_motifs_sorted_by_fc[upregulated_motifs_sorted_by_fc$p.adjust < 0.05,]
              upregulated_motifs_sorted_by_fc <- upregulated_motifs_sorted_by_fc[rev(order(upregulated_motifs_sorted_by_fc$fc_value)),]
              # Add FC values for downregulated motifs
              motif_fc_values <- c()
              for(current_row in 1:nrow(downregulated_motifs)) {
                current_motif <- downregulated_motifs[current_row,]
                current_motif_percent_observed <- current_motif$percent.observed
                alternative_motif <- upregulated_motifs[upregulated_motifs$motif %in% current_motif$motif,]
                alternative_motif_percent_observed <- alternative_motif$percent.observed
                motif_fc_value <- -log(current_motif_percent_observed / alternative_motif_percent_observed)
                motif_fc_values <- c(motif_fc_values, motif_fc_value)
              }
              downregulated_motifs$fc_value <- motif_fc_values
              downregulated_motifs_sorted_by_fc <- downregulated_motifs
              downregulated_motifs_sorted_by_fc <- downregulated_motifs_sorted_by_fc[downregulated_motifs_sorted_by_fc$p.adjust < 0.05,]
              downregulated_motifs_sorted_by_fc <- downregulated_motifs_sorted_by_fc[order(downregulated_motifs_sorted_by_fc$fc_value),]
              # Save to file
              write.table(upregulated_motifs, file = paste0(current_output_dir, "D28-vs-D_minus_1-degs-", snATAC_cell_type_for_file_name, "-", analysis_type, "_pct_", pct, "_FC_", fc, "_with_fc_values.tsv"),
                          quote = FALSE, sep = "\t")
              write.table(upregulated_motifs_sorted_by_fc, file = paste0(current_output_dir, "D28-vs-D_minus_1-degs-", snATAC_cell_type_for_file_name, "-", analysis_type, "_pct_", pct, "_FC_", fc, "_with_fc_values_sorted.tsv"),
                          quote = FALSE, sep = "\t")
              write.table(downregulated_motifs, file = paste0(current_output_dir, "D28-vs-D_minus_1-degs-", snATAC_cell_type_for_file_name, "-", analysis_type, "_pct_", pct, "_FC_-", fc, "_with_fc_values.tsv"),
                          quote = FALSE, sep = "\t")
              write.table(downregulated_motifs_sorted_by_fc, file = paste0(current_output_dir, "D28-vs-D_minus_1-degs-", snATAC_cell_type_for_file_name, "-", analysis_type, "_pct_", pct, "_FC_-", fc, "_with_fc_values_sorted.tsv"),
                          quote = FALSE, sep = "\t")
            }
          }
        }
      }
    }
  }
}


# Step 2: Create plots

# CD14 Mono (Naive HVL)
cd14_mono_upregulated_motifs <- read.table(paste0(scATAC_hvl_placebo_das_motif_dir, "CD14_Mono/sc/0.01/with_bg/with_fc_added/D28-vs-D_minus_1-degs-CD14_Mono-sc_pct_0.01_FC_1_with_fc_values.tsv"),
                                           sep = "\t", header = TRUE)
cd14_mono_upregulated_motifs <- cd14_mono_upregulated_motifs[rowSums(is.na(cd14_mono_upregulated_motifs)) == 0, ] # Remove NAs
cd14_mono_upregulated_motifs <- cd14_mono_upregulated_motifs[cd14_mono_upregulated_motifs$fc_value > 0,]
cd14_mono_upregulated_motifs$p.adjust.log <- -log(cd14_mono_upregulated_motifs$p.adjust, base = 10)

cd14_mono_downregulated_motifs <-  read.table(paste0(scATAC_hvl_placebo_das_motif_dir, "CD14_Mono/sc/0.01/with_bg/with_fc_added/D28-vs-D_minus_1-degs-CD14_Mono-sc_pct_0.01_FC_-1_with_fc_values.tsv"),
                                              sep = "\t", header = TRUE)
cd14_mono_downregulated_motifs <- cd14_mono_downregulated_motifs[rowSums(is.na(cd14_mono_downregulated_motifs)) == 0, ] # Remove NAs
cd14_mono_downregulated_motifs <- cd14_mono_downregulated_motifs[cd14_mono_downregulated_motifs$fc_value < 0,]
cd14_mono_downregulated_motifs$p.adjust.log <- -log(cd14_mono_downregulated_motifs$p.adjust, base = 10)

cd14_mono_upregulated_motif_df_for_plotting <- data.frame(tf = cd14_mono_upregulated_motifs$motif.name,
                                                          fc = cd14_mono_upregulated_motifs$fc_value,
                                                          p.adjust.log = cd14_mono_upregulated_motifs$p.adjust.log,
                                                          peak_direction = "upregulated")

cd14_mono_downregulated_motif_df_for_plotting <- data.frame(tf = cd14_mono_downregulated_motifs$motif.name,
                                                          fc = cd14_mono_downregulated_motifs$fc_value,
                                                          p.adjust.log = cd14_mono_downregulated_motifs$p.adjust.log,
                                                          peak_direction = "downregulated")

combined_cd14_mono_motif_df_for_plotting <- rbind(cd14_mono_upregulated_motif_df_for_plotting, cd14_mono_downregulated_motif_df_for_plotting)

combined_cd14_mono_motif_df_for_plotting$color <- with(combined_cd14_mono_motif_df_for_plotting, 
                                                       ifelse(p.adjust.log >= 7.5 & abs(fc) >= 0.5 & peak_direction == "upregulated", "red",
                                                              ifelse(p.adjust.log >= 7.5 & abs(fc) >= 0.5 & peak_direction == "downregulated", "blue", "grey")))

cd14_mono_motif_plot <- ggplot(combined_cd14_mono_motif_df_for_plotting, aes(x = fc, y = p.adjust.log, label = tf)) +
  geom_point(aes(color = color), size = 3) +
  scale_color_identity() +
  # geom_text(aes(label = ifelse(color %in% c("red", "blue"), tf, '')), vjust = 2, hjust = 0.5, size = 3) +
  theme_classic(base_size = 30) +
  labs(title = "CD14 Mono",
       x = "Log2(FC)",
       y = "-Log10(Adjusted P-value)") +
  theme(plot.title = element_text(hjust = 0.5)) + geom_vline(xintercept=c(0), linetype="dotted", color = "grey")

ggsave(filename = "C:/Users/willi/Desktop/CD14_Mono_Motif.png", plot = cd14_mono_motif_plot, width = 1800, height = 1703, units = "px")

# CD14 Mono (Naive LVL)
cd14_mono_upregulated_motifs <- read.table(paste0(scATAC_lvl_placebo_das_motif_dir, "CD14_Mono/sc/0.01/with_bg/with_fc_added/D28-vs-D_minus_1-degs-CD14_Mono-sc_pct_0.01_FC_1_with_fc_values.tsv"),
                                           sep = "\t", header = TRUE)
cd14_mono_upregulated_motifs <- cd14_mono_upregulated_motifs[rowSums(is.na(cd14_mono_upregulated_motifs)) == 0, ] # Remove NAs
cd14_mono_upregulated_motifs <- cd14_mono_upregulated_motifs[cd14_mono_upregulated_motifs$fc_value > 0,]
cd14_mono_upregulated_motifs$p.adjust.log <- -log(cd14_mono_upregulated_motifs$p.adjust, base = 10)

cd14_mono_downregulated_motifs <-  read.table(paste0(scATAC_lvl_placebo_das_motif_dir, "CD14_Mono/sc/0.01/with_bg/with_fc_added/D28-vs-D_minus_1-degs-CD14_Mono-sc_pct_0.01_FC_-1_with_fc_values.tsv"),
                                              sep = "\t", header = TRUE)
cd14_mono_downregulated_motifs <- cd14_mono_downregulated_motifs[rowSums(is.na(cd14_mono_downregulated_motifs)) == 0, ] # Remove NAs
cd14_mono_downregulated_motifs <- cd14_mono_downregulated_motifs[cd14_mono_downregulated_motifs$fc_value < 0,]
cd14_mono_downregulated_motifs$p.adjust.log <- -log(cd14_mono_downregulated_motifs$p.adjust, base = 10)

cd14_mono_upregulated_motif_df_for_plotting <- data.frame(tf = cd14_mono_upregulated_motifs$motif.name,
                                                          fc = cd14_mono_upregulated_motifs$fc_value,
                                                          p.adjust.log = cd14_mono_upregulated_motifs$p.adjust.log,
                                                          peak_direction = "upregulated")

cd14_mono_downregulated_motif_df_for_plotting <- data.frame(tf = cd14_mono_downregulated_motifs$motif.name,
                                                            fc = cd14_mono_downregulated_motifs$fc_value,
                                                            p.adjust.log = cd14_mono_downregulated_motifs$p.adjust.log,
                                                            peak_direction = "downregulated")

combined_cd14_mono_motif_df_for_plotting <- rbind(cd14_mono_upregulated_motif_df_for_plotting, cd14_mono_downregulated_motif_df_for_plotting)

combined_cd14_mono_motif_df_for_plotting$color <- with(combined_cd14_mono_motif_df_for_plotting, 
                                                       ifelse(p.adjust.log >= 7.5 & abs(fc) >= 0.5 & peak_direction == "upregulated", "red",
                                                              ifelse(p.adjust.log >= 7.5 & abs(fc) >= 0.5 & peak_direction == "downregulated", "blue", "grey")))

cd14_mono_motif_plot <- ggplot(combined_cd14_mono_motif_df_for_plotting, aes(x = fc, y = p.adjust.log, label = tf)) +
  geom_point(aes(color = color), size = 3) +
  scale_color_identity() +
  # geom_text(aes(label = ifelse(color %in% c("red", "blue"), tf, '')), vjust = 2, hjust = 0.5, size = 3) +
  theme_classic(base_size = 30) +
  labs(title = "CD14 Mono",
       x = "Log2(FC)",
       y = "-Log10(Adjusted P-value)") +
  theme(plot.title = element_text(hjust = 0.5)) + geom_vline(xintercept=c(0), linetype="dotted", color = "grey")

ggsave(filename = "C:/Users/willi/Desktop/CD14_Mono_Motif.png", plot = cd14_mono_motif_plot, width = 1800, height = 1703, units = "px")


# CD16 Mono
cd16_mono_upregulated_motifs <- read.table(paste0(scATAC_hvl_placebo_das_motif_dir, "CD16_Mono/sc/0.01/with_bg/with_fc_added/D28-vs-D_minus_1-degs-cd16_Mono-sc_pct_0.01_FC_1_with_fc_values.tsv"),
                                           sep = "\t", header = TRUE)
cd16_mono_upregulated_motifs <- cd16_mono_upregulated_motifs[rowSums(is.na(cd16_mono_upregulated_motifs)) == 0, ] # Remove NAs
cd16_mono_upregulated_motifs <- cd16_mono_upregulated_motifs[cd16_mono_upregulated_motifs$fc_value > 0,]
cd16_mono_upregulated_motifs$p.adjust.log <- -log(cd16_mono_upregulated_motifs$p.adjust, base = 10)

cd16_mono_downregulated_motifs <-  read.table(paste0(scATAC_hvl_placebo_das_motif_dir, "CD16_Mono/sc/0.01/with_bg/with_fc_added/D28-vs-D_minus_1-degs-cd16_Mono-sc_pct_0.01_FC_-1_with_fc_values.tsv"),
                                              sep = "\t", header = TRUE)
cd16_mono_downregulated_motifs <- cd16_mono_downregulated_motifs[rowSums(is.na(cd16_mono_downregulated_motifs)) == 0, ] # Remove NAs
cd16_mono_downregulated_motifs <- cd16_mono_downregulated_motifs[cd16_mono_downregulated_motifs$fc_value < 0,]
cd16_mono_downregulated_motifs$p.adjust.log <- -log(cd16_mono_downregulated_motifs$p.adjust, base = 10)

cd16_mono_upregulated_motif_df_for_plotting <- data.frame(tf = cd16_mono_upregulated_motifs$motif.name,
                                                          fc = cd16_mono_upregulated_motifs$fc_value,
                                                          p.adjust.log = cd16_mono_upregulated_motifs$p.adjust.log,
                                                          peak_direction = "upregulated")

cd16_mono_downregulated_motif_df_for_plotting <- data.frame(tf = cd16_mono_downregulated_motifs$motif.name,
                                                            fc = cd16_mono_downregulated_motifs$fc_value,
                                                            p.adjust.log = cd16_mono_downregulated_motifs$p.adjust.log,
                                                            peak_direction = "downregulated")

combined_cd16_mono_motif_df_for_plotting <- rbind(cd16_mono_upregulated_motif_df_for_plotting, cd16_mono_downregulated_motif_df_for_plotting)

combined_cd16_mono_motif_df_for_plotting$color <- with(combined_cd16_mono_motif_df_for_plotting, 
                                                       ifelse(p.adjust.log >= 5 & abs(fc) >= 0.5 & peak_direction == "upregulated", "red",
                                                              ifelse(p.adjust.log >= 5 & abs(fc) >= 0.5 & peak_direction == "downregulated", "blue", "grey")))

cd16_mono_motif_plot <- ggplot(combined_cd16_mono_motif_df_for_plotting, aes(x = fc, y = p.adjust.log, label = tf)) +
  geom_point(aes(color = color), size = 4) +
  scale_color_identity() +
  # geom_text(aes(label = ifelse(color %in% c("red", "blue"), tf, '')), vjust = 2, hjust = 0.5, size = 3) +
  theme_classic(base_size = 30) +
  labs(title = "CD16 Mono",
       x = "Log2(FC)",
       y = "-Log10(Adjusted P-value)") +
  theme(plot.title = element_text(hjust = 0.5)) + geom_vline(xintercept=c(0), linetype="dotted", color = "grey")

ggsave(filename = "C:/Users/willi/Desktop/CD16_Mono_Motif.png", plot = cd16_mono_motif_plot, width = 1800, height = 1703, units = "px")

# cDC
cDC_upregulated_motifs <- read.table(paste0(scATAC_hvl_placebo_das_motif_dir, "cDC/sc/0.01/with_bg/with_fc_added/D28-vs-D_minus_1-degs-cDC-sc_pct_0.01_FC_1_with_fc_values.tsv"),
                                           sep = "\t", header = TRUE)
cDC_upregulated_motifs <- cDC_upregulated_motifs[rowSums(is.na(cDC_upregulated_motifs)) == 0, ] # Remove NAs
cDC_upregulated_motifs <- cDC_upregulated_motifs[cDC_upregulated_motifs$fc_value > 0,]
cDC_upregulated_motifs$p.adjust.log <- -log(cDC_upregulated_motifs$p.adjust, base = 10)

cDC_downregulated_motifs <-  read.table(paste0(scATAC_hvl_placebo_das_motif_dir, "cDC/sc/0.01/with_bg/with_fc_added/D28-vs-D_minus_1-degs-cDC-sc_pct_0.01_FC_-1_with_fc_values.tsv"),
                                              sep = "\t", header = TRUE)
cDC_downregulated_motifs <- cDC_downregulated_motifs[rowSums(is.na(cDC_downregulated_motifs)) == 0, ] # Remove NAs
cDC_downregulated_motifs <- cDC_downregulated_motifs[cDC_downregulated_motifs$fc_value < 0,]
cDC_downregulated_motifs$p.adjust.log <- -log(cDC_downregulated_motifs$p.adjust, base = 10)

cDC_upregulated_motif_df_for_plotting <- data.frame(tf = cDC_upregulated_motifs$motif.name,
                                                          fc = cDC_upregulated_motifs$fc_value,
                                                          p.adjust.log = cDC_upregulated_motifs$p.adjust.log,
                                                          peak_direction = "upregulated")

cDC_downregulated_motif_df_for_plotting <- data.frame(tf = cDC_downregulated_motifs$motif.name,
                                                            fc = cDC_downregulated_motifs$fc_value,
                                                            p.adjust.log = cDC_downregulated_motifs$p.adjust.log,
                                                            peak_direction = "downregulated")

combined_cDC_motif_df_for_plotting <- rbind(cDC_upregulated_motif_df_for_plotting, cDC_downregulated_motif_df_for_plotting)

combined_cDC_motif_df_for_plotting$color <- with(combined_cDC_motif_df_for_plotting, 
                                                       ifelse(p.adjust.log >= 10 & abs(fc) >= 0.5 & peak_direction == "upregulated", "red",
                                                              ifelse(p.adjust.log >= 10 & abs(fc) >= 0.5 & peak_direction == "downregulated", "blue", "grey")))

cDC_motif_plot <- ggplot(combined_cDC_motif_df_for_plotting, aes(x = fc, y = p.adjust.log, label = tf)) +
  geom_point(aes(color = color), size = 4) +
  scale_color_identity() +
  # geom_text(aes(label = ifelse(color %in% c("red", "blue"), tf, '')), vjust = 2, hjust = 0.5, size = 3) +
  theme_classic(base_size = 30) +
  labs(title = "cDCs",
       x = "Log2(FC)",
       y = "-Log10(Adjusted P-value)") +
  theme(plot.title = element_text(hjust = 0.5)) + geom_vline(xintercept=c(0), linetype="dotted", color = "grey")

ggsave(filename = "C:/Users/willi/Desktop/cDC_Motif.png", plot = cDC_motif_plot, width = 1800, height = 1703, units = "px")

# pDC
pDC_upregulated_motifs <- read.table(paste0(scATAC_hvl_placebo_das_motif_dir, "pDC/sc/0.01/with_bg/with_fc_added/D28-vs-D_minus_1-degs-pDC-sc_pct_0.01_FC_1_with_fc_values.tsv"),
                                     sep = "\t", header = TRUE)
pDC_upregulated_motifs <- pDC_upregulated_motifs[rowSums(is.na(pDC_upregulated_motifs)) == 0, ] # Remove NAs
pDC_upregulated_motifs <- pDC_upregulated_motifs[pDC_upregulated_motifs$fc_value > 0,]
pDC_upregulated_motifs$p.adjust.log <- -log(pDC_upregulated_motifs$p.adjust, base = 10)

pDC_downregulated_motifs <-  read.table(paste0(scATAC_hvl_placebo_das_motif_dir, "pDC/sc/0.01/with_bg/with_fc_added/D28-vs-D_minus_1-degs-pDC-sc_pct_0.01_FC_-1_with_fc_values.tsv"),
                                        sep = "\t", header = TRUE)
pDC_downregulated_motifs <- pDC_downregulated_motifs[rowSums(is.na(pDC_downregulated_motifs)) == 0, ] # Remove NAs
pDC_downregulated_motifs <- pDC_downregulated_motifs[pDC_downregulated_motifs$fc_value < 0,]
pDC_downregulated_motifs$p.adjust.log <- -log(pDC_downregulated_motifs$p.adjust, base = 10)

pDC_upregulated_motif_df_for_plotting <- data.frame(tf = pDC_upregulated_motifs$motif.name,
                                                    fc = pDC_upregulated_motifs$fc_value,
                                                    p.adjust.log = pDC_upregulated_motifs$p.adjust.log,
                                                    peak_direction = "upregulated")

pDC_downregulated_motif_df_for_plotting <- data.frame(tf = pDC_downregulated_motifs$motif.name,
                                                      fc = pDC_downregulated_motifs$fc_value,
                                                      p.adjust.log = pDC_downregulated_motifs$p.adjust.log,
                                                      peak_direction = "downregulated")

combined_pDC_motif_df_for_plotting <- rbind(pDC_upregulated_motif_df_for_plotting, pDC_downregulated_motif_df_for_plotting)

combined_pDC_motif_df_for_plotting$color <- with(combined_pDC_motif_df_for_plotting, 
                                                 ifelse(p.adjust.log >= 3.5 & abs(fc) >= 0.5 & peak_direction == "upregulated", "red",
                                                        ifelse(p.adjust.log >= 3.5 & abs(fc) >= 0.5 & peak_direction == "downregulated", "blue", "grey")))

pDC_motif_plot <- ggplot(combined_pDC_motif_df_for_plotting, aes(x = fc, y = p.adjust.log, label = tf)) +
  geom_point(aes(color = color), size = 4) +
  scale_color_identity() +
  # geom_text(aes(label = ifelse(color %in% c("red", "blue"), tf, '')), vjust = 2, hjust = 0.5, size = 3) +
  theme_classic(base_size = 30) +
  labs(title = "pDCs",
       x = "Log2(FC)",
       y = "-Log10(Adjusted P-value)") +
  theme(plot.title = element_text(hjust = 0.5)) + geom_vline(xintercept=c(0), linetype="dotted", color = "grey")

ggsave(filename = "C:/Users/willi/Desktop/pDC_Motif.png", plot = pDC_motif_plot, width = 1800, height = 1703, units = "px")

