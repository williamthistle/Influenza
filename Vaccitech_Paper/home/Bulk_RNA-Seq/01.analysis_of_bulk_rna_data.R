### THIS FUNCTION RUNS DESEQ2 ON ALL BULK RNA-SEQ GROUPINGS
### WE CAN LOOK AT HIGH EFFECT SIZE GENES INDIVIDUALLY (NO NEED TO DO FC THRESHOLD?)
### IT ALSO RUNS FUNCTIONAL MODULE DISCOVERY ON DIFFERENT FOLD CHANGE THRESHOLDS
### FOR EACH BULK RNA-SEQ GROUPING
### THESE HB PLOTS (OR DOWNSTREAM REACTOME PLOTS) MAY GO IN SUPPLEMENTARY FIGURES
### IN MAIN TEXT, DISCUSS EACH DAY AND WHAT WE SEE (HIGH EFFECT SIZE GENES, MODULES, PATHWAYS)

# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

################ PLACEBO ################ 
#### PERIOD 1 HIGH VIRAL LOAD ####
# We expect ~0 DEGs because placebo was used (no actual vaccination)
# 0/0/0/0/0/0/0 DEGs
high_placebo_period_1_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", high_placebo_counts, high_placebo_metadata,
                                                                                "1_D2", "1_D_minus_1", paste0(bulk_results_dir, "hvl_bulk_placebo_period_1_D2_vs_D_minus_1/"), "high")
raw_high_placebo_period_1_D2_vs_D_minus_1_results <- high_placebo_period_1_D2_vs_D_minus_1_results[[1]]
# 0/0/0/0/0/0/0 DEGs
high_placebo_period_1_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", high_placebo_counts, high_placebo_metadata,
                                                                                "1_D8", "1_D_minus_1", paste0(bulk_results_dir, "hvl_bulk_period_1_D8_vs_D_minus_1/"), "high")
raw_high_placebo_period_1_D8_vs_D_minus_1_results <- high_placebo_period_1_D8_vs_D_minus_1_results[[1]]
# 16/1/1/1/1/1/1 DEGs
high_placebo_period_1_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", high_placebo_counts, high_placebo_metadata,
                                                                                "1_D28", "1_D_minus_1", paste0(bulk_results_dir, "hvl_bulk_period_1_D28_vs_D_minus_1/"), "high")
raw_high_placebo_period_1_D28_vs_D_minus_1_results <- high_placebo_period_1_D28_vs_D_minus_1_results[[1]]
#### PERIOD 1 LOW VIRAL LOAD ####
# We expect ~0 DEGs because placebo was used (no actual vaccination)
# 1/1/1/1/1/0/0 DEGs
low_placebo_period_1_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", low_placebo_counts, low_placebo_metadata,
                                                                         "1_D2", "1_D_minus_1", paste0(bulk_results_dir, "lvl_bulk_placebo_period_1_D2_vs_D_minus_1/"), "low")
raw_low_placebo_period_1_D2_vs_D_minus_1_results <- low_placebo_period_1_D2_vs_D_minus_1_results[[1]]
# 0/0/0/0/0/0/0 DEGs
low_placebo_period_1_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", low_placebo_counts, low_placebo_metadata,
                                                                         "1_D8", "1_D_minus_1", paste0(bulk_results_dir, "lvl_bulk_placebo_period_1_D8_vs_D_minus_1/"), "low")
raw_low_placebo_period_1_D8_vs_D_minus_1_results <- low_placebo_period_1_D8_vs_D_minus_1_results[[1]]
# 63/19/0/0/0/0/0 DEGs
low_placebo_period_1_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", low_placebo_counts, low_placebo_metadata,
                                                                          "1_D28", "1_D_minus_1", paste0(bulk_results_dir, "lvl_bulk_placebo_period_1_D28_vs_D_minus_1/"), "low")
raw_low_placebo_period_1_D28_vs_D_minus_1_results <- low_placebo_period_1_D28_vs_D_minus_1_results[[1]]

#### PERIOD 2 HIGH VIRAL LOAD ####
# 2 D minus 1 vs 2 D minus 2 - should be virtually zero unless some weird stuff happened between blood draws
# 2154/609/219/87/1/0/0 DEGs found - seems weird that so many are found
# Could there be something that happened shortly before infection that made subjects more susceptible to getting sick?
# Or is there something about these subjects that make them more predisposed to getting sick?
# Could it be a single patient (or two) messing everything up?
# Is it just inflammation from getting blood drawn?
# If so, why would this signal exist for HVL and not LVL?
high_placebo_period_2_D_minus_1_vs_D_minus_2_results <- run_deseq_bulk_analysis_time_series("placebo", high_placebo_counts, high_placebo_metadata,
                                                                                "2_D_minus_1", "2_D_minus_2", paste0(bulk_results_dir, "hvl_bulk_placebo_period_2_D_minus_1_vs_D_minus_2/"), "high")
raw_high_placebo_period_2_D_minus_1_vs_D_minus_2_results <- high_placebo_period_2_D_minus_1_vs_D_minus_2_results[[1]]

# 2 D2 vs 2 D minus 1 - 1384/299/17/4/0/0/0 DEGs
high_placebo_period_2_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", high_placebo_counts, high_placebo_metadata,
                                                                         "2_D2", "2_D_minus_1", paste0(bulk_results_dir, "hvl_bulk_placebo_period_2_D2_vs_D_minus_1/"), "high")
raw_high_placebo_period_2_D2_vs_D_minus_1_results <- high_placebo_period_2_D2_vs_D_minus_1_results[[1]]

# 2 D5 vs 2 D minus 1 - 5276/3432/2230/1556/600/228/51 DEGs
high_placebo_period_2_D5_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", high_placebo_counts, high_placebo_metadata,
                                                                         "2_D5", "2_D_minus_1", paste0(bulk_results_dir, "hvl_bulk_placebo_period_2_D5_vs_D_minus_1/"), "high")
raw_high_placebo_period_2_D5_vs_D_minus_1_results <- high_placebo_period_2_D5_vs_D_minus_1_results[[1]]

# 2 D8 vs 2 D minus 1 - 3314/1940/1210/813/269/59/8 DEGs
high_placebo_period_2_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", high_placebo_counts, high_placebo_metadata,
                                                                         "2_D8", "2_D_minus_1", paste0(bulk_results_dir, "hvl_bulk_placebo_period_2_D8_vs_D_minus_1/"), "high")
raw_high_placebo_period_2_D8_vs_D_minus_1_results <- high_placebo_period_2_D8_vs_D_minus_1_results[[1]]

# 2 D28 vs 2 D minus 1 - 1394/34/3/2/2/0/0 DEGs
high_placebo_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", high_placebo_counts, high_placebo_metadata,
                                                                          "2_D28", "2_D_minus_1", paste0(bulk_results_dir, "hvl_bulk_placebo_period_2_D28_vs_D_minus_1/"), "high")
raw_high_placebo_period_2_D28_vs_D_minus_1_results <- high_placebo_period_2_D28_vs_D_minus_1_results[[1]]

#### PERIOD 2 LOW VIRAL LOAD ####
# 2 D minus 2 vs 2 D minus 1 - should be virtually zero unless some weird stuff happened between blood draws
# 0/0/0/0/0/0/0 DEGs found
low_placebo_period_2_D_minus_1_vs_D_minus_2_results <- run_deseq_bulk_analysis_time_series("placebo", low_placebo_counts, low_placebo_metadata,
                                                                               "2_D_minus_1", "2_D_minus_2", paste0(bulk_results_dir, "lvl_bulk_placebo_period_2_D_minus_1_vs_D_minus_2/"), "low")
raw_low_placebo_period_2_D_minus_1_vs_D_minus_2_results <- low_placebo_period_2_D_minus_1_vs_D_minus_2_results[[1]]
# 2 D2 vs 2 D minus 1 - 2/1/1/1/1/1/1 DEGs
low_placebo_period_2_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", low_placebo_counts, low_placebo_metadata,
                                                                        "2_D2", "2_D_minus_1", paste0(bulk_results_dir, "lvl_bulk_placebo_period_2_D2_vs_D_minus_1/"), "low")
raw_low_placebo_period_2_D2_vs_D_minus_1_results <- low_placebo_period_2_D2_vs_D_minus_1_results[[1]]
# 2 D5 vs 2 D minus 1 - 1/1/1/0/0/0/0 DEGs
low_placebo_period_2_D5_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", low_placebo_counts, low_placebo_metadata,
                                                                        "2_D5", "2_D_minus_1", paste0(bulk_results_dir, "lvl_bulk_placebo_period_2_D5_vs_D_minus_1/"), "low")
raw_low_placebo_period_2_D5_vs_D_minus_1_results <- low_placebo_period_2_D5_vs_D_minus_1_results[[1]]
# 2 D8 vs 2 D minus 1 - 0/0/0/0/0/0/0 DEGs
low_placebo_period_2_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", low_placebo_counts, low_placebo_metadata,
                                                                        "2_D8", "2_D_minus_1", paste0(bulk_results_dir, "lvl_bulk_placebo_period_2_D8_vs_D_minus_1/"), "low")
raw_low_placebo_period_2_D8_vs_D_minus_1_results <- low_placebo_period_2_D8_vs_D_minus_1_results[[1]]
# 2 D28 vs 2 D minus 1 - 647/74/0/0/0/0/0 DEGs
low_placebo_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", low_placebo_counts, low_placebo_metadata,
                                                                         "2_D28", "2_D_minus_1", paste0(bulk_results_dir, "lvl_bulk_placebo_period_2_D28_vs_D_minus_1/"), "low")
raw_low_placebo_period_2_D28_vs_D_minus_1_results <- low_placebo_period_2_D28_vs_D_minus_1_results[[1]]

################ VACCINATED (HVL AND LVL) ################
# 2 D minus 2 vs 2 D minus 1 - should be virtually zero unless some weird stuff happened between blood draws
# 1272/458/169/67/3/2/0 DEGs
hvl_vaccinated_period_2_D_minus_1_vs_D_minus_2_results <- run_deseq_bulk_analysis_time_series("vaccinated", hvl_vaccinated_counts, hvl_vaccinated_metadata,
                                                                                            "2_D_minus_1", "2_D_minus_2", paste0(bulk_results_dir, "hvl_bulk_vaccinated_period_2_D_minus_1_vs_D_minus_2/"), "high")
raw_hvl_vaccinated_period_2_D_minus_1_vs_D_minus_2_results <- hvl_vaccinated_period_2_D_minus_1_vs_D_minus_2_results[[1]]

# 2 D2 vs 2 D minus 1 - 0/0/0/0/0/0/0 DEGs
hvl_vaccinated_period_2_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", hvl_vaccinated_counts, hvl_vaccinated_metadata,
                                                                                     "2_D2", "2_D_minus_1", paste0(bulk_results_dir, "hvl_bulk_vaccinated_period_2_D2_vs_D_minus_1/"), "high")
raw_hvl_vaccinated_period_2_D2_vs_D_minus_1_results <- hvl_vaccinated_period_2_D2_vs_D_minus_1_results[[1]]

# 2 D5 vs 2 D minus 1 - 72/60/49/45/25/11/0 DEGs
hvl_vaccinated_period_2_D5_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", hvl_vaccinated_counts, hvl_vaccinated_metadata,
                                                                                     "2_D5", "2_D_minus_1", paste0(bulk_results_dir, "hvl_bulk_vaccinated_period_2_D5_vs_D_minus_1/"), "high")
raw_hvl_vaccinated_period_2_D5_vs_D_minus_1_results <- hvl_vaccinated_period_2_D5_vs_D_minus_1_results[[1]]

# 2 D8 vs 2 D minus 1 - 1020/300/157/93/33/19/1 DEGs
hvl_vaccinated_period_2_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", hvl_vaccinated_counts, hvl_vaccinated_metadata,
                                                                                     "2_D8", "2_D_minus_1", paste0(bulk_results_dir, "hvl_bulk_vaccinated_period_2_D8_vs_D_minus_1/"), "high")
raw_hvl_vaccinated_period_2_D8_vs_D_minus_1_results <- hvl_vaccinated_period_2_D8_vs_D_minus_1_results[[1]]

# 2 D28 vs 2 D minus 1 - 5860/2485/580/101/1/0/0 DEGs
hvl_vaccinated_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", hvl_vaccinated_counts, hvl_vaccinated_metadata,
                                                                                      "2_D28", "2_D_minus_1", paste0(bulk_results_dir, "hvl_bulk_vaccinated_period_2_D28_vs_D_minus_1/"), "high")
raw_hvl_vaccinated_period_2_D28_vs_D_minus_1_results <- hvl_vaccinated_period_2_D28_vs_D_minus_1_results[[1]]

#### PERIOD 2 LOW VIRAL LOAD ####
# 2 D minus 2 vs 2 D minus 1 - should be virtually zero unless some weird stuff happened between blood draws
# 194/95/46/4/1/0/0 DEGs found
lvl_vaccinated_period_2_D_minus_1_vs_D_minus_2_results <- run_deseq_bulk_analysis_time_series("vaccinated", lvl_vaccinated_counts, lvl_vaccinated_metadata,
                                                                                           "2_D_minus_1", "2_D_minus_2", paste0(bulk_results_dir, "lvl_bulk_vaccinated_period_2_D_minus_1_vs_D_minus_2/"), "low")
raw_lvl_vaccinated_period_2_D_minus_1_vs_D_minus_2_results <- lvl_vaccinated_period_2_D_minus_1_vs_D_minus_2_results[[1]]
# 2 D2 vs 2 D minus 1 - 16/14/9/7/3/2/1 DEGs
lvl_vaccinated_period_2_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", lvl_vaccinated_counts, lvl_vaccinated_metadata,
                                                                                    "2_D2", "2_D_minus_1", paste0(bulk_results_dir, "lvl_bulk_vaccinated_period_2_D2_vs_D_minus_1/"), "low")
raw_lvl_vaccinated_period_2_D2_vs_D_minus_1_results <- lvl_vaccinated_period_2_D2_vs_D_minus_1_results[[1]]
# 2 D5 vs 2 D minus 1 - 133/30/4/4/2/1/0 DEGs
lvl_vaccinated_period_2_D5_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", lvl_vaccinated_counts, lvl_vaccinated_metadata,
                                                                                    "2_D5", "2_D_minus_1", paste0(bulk_results_dir, "lvl_bulk_vaccinated_period_2_D5_vs_D_minus_1/"), "low")
raw_lvl_vaccinated_period_2_D5_vs_D_minus_1_results <- lvl_vaccinated_period_2_D5_vs_D_minus_1_results[[1]]
# 2 D8 vs 2 D minus 1 - 202/2/0/0/0/0/0 DEGs
lvl_vaccinated_period_2_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", lvl_vaccinated_counts, lvl_vaccinated_metadata,
                                                                                    "2_D8", "2_D_minus_1", paste0(bulk_results_dir, "lvl_bulk_vaccinated_period_2_D8_vs_D_minus_1/"), "low")
raw_lvl_vaccinated_period_2_D8_vs_D_minus_1_results <- lvl_vaccinated_period_2_D8_vs_D_minus_1_results[[1]]
# 2 D28 vs 2 D minus 1 - 2386/768/166/23/5/2/1 DEGs
lvl_vaccinated_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", lvl_vaccinated_counts, lvl_vaccinated_metadata,
                                                                                     "2_D28", "2_D_minus_1", paste0(bulk_results_dir, "lvl_bulk_vaccinated_period_2_D28_vs_D_minus_1/"), "low")
raw_lvl_vaccinated_period_2_D28_vs_D_minus_1_results <- lvl_vaccinated_period_2_D28_vs_D_minus_1_results[[1]]

# Save results to disk
save.image(paste0(bulk_results_dir, "bulk_RNA_analysis.RData"))

# Compare fold changes for overlapping genes between vaccinated and placebo for D28
overlapping_genes <- intersect(rownames(raw_hvl_vaccinated_period_2_D28_vs_D_minus_1_results), rownames(raw_high_placebo_period_2_D28_vs_D_minus_1_results))
compare_placebo_D28_df <- raw_high_placebo_period_2_D28_vs_D_minus_1_results[rownames(raw_high_placebo_period_2_D28_vs_D_minus_1_results) %in% overlapping_genes,]
compare_vaccinated_D28_df <- raw_hvl_vaccinated_period_2_D28_vs_D_minus_1_results[rownames(raw_hvl_vaccinated_period_2_D28_vs_D_minus_1_results) %in% overlapping_genes,]
compare_placebo_D28_df <- compare_placebo_D28_df[order(rownames(compare_placebo_D28_df)),]
compare_vaccinated_D28_df <- compare_vaccinated_D28_df[order(rownames(compare_vaccinated_D28_df)),]

comparing_placebo_vs_vaccinated_D28_df <- data.frame(gene_name = rownames(compare_placebo_D28_df), placebo_fc = compare_placebo_D28_df$log2FoldChange,
                                                     vaccinated_fc = compare_vaccinated_D28_df$log2FoldChange)
# Calculate correlation - very high correlation for those overlapping genes
cor(comparing_placebo_vs_vaccinated_D28_df$placebo_fc, comparing_placebo_vs_vaccinated_D28_df$vaccinated_fc)

# Plot correlation
ggplot(comparing_placebo_vs_vaccinated_D28_df, aes(x=placebo_fc, y=vaccinated_fc)) + 
  geom_point()+
  geom_smooth(method=lm)

# HVL PLACEBO VS HVL VACCINATED
# Many more DEGs for D8 vs D28 instead of D28 vs D28 (expected, but nice to see)

high_placebo_metadata_for_comparison <- high_placebo_metadata[high_placebo_metadata$time_point == "2_D28",]
high_placebo_metadata_for_comparison$status <- "placebo"
hvl_vaccinated_metadata_for_comparison <- hvl_vaccinated_metadata[hvl_vaccinated_metadata$time_point == "2_D8",]
hvl_vaccinated_metadata_for_comparison$status <- "vaccinated"
hvl_placebo_and_vaccinated_metadata <- rbind(high_placebo_metadata_for_comparison, hvl_vaccinated_metadata_for_comparison)
hvl_placebo_and_vaccinated_counts <- cbind(high_placebo_counts, hvl_vaccinated_counts)

# Select the relevant time point from our metadata
hvl_placebo_and_vaccinated_metadata_subset <- hvl_placebo_and_vaccinated_metadata
# Select subset of counts associated with subjects
hvl_placebo_and_vaccinated_counts_subset <- hvl_placebo_and_vaccinated_counts[rownames(hvl_placebo_and_vaccinated_metadata_subset)]
# Run DESeq2
current_analysis <- DESeqDataSetFromMatrix(countData = hvl_placebo_and_vaccinated_counts_subset, 
                                           colData = hvl_placebo_and_vaccinated_metadata_subset, 
                                           design = ~ sex + age + status)
current_analysis <- DESeq(current_analysis)
current_analysis_results <- results(current_analysis, contrast = c("status", "placebo", "vaccinated"), alpha = 0.05, lfcThreshold = 0.1)
current_analysis_results <- current_analysis_results[order(current_analysis_results$padj),]
current_analysis_results <- subset(current_analysis_results, padj < 0.05)

# HIGH VS LOW VIRAL LOAD
# 2 D minus 2
# 144/143/142/138/151
# Note: If I don't control for age and sex, I find nothing.
placebo_period_2_D_minus_2_high_vs_low_results <- run_deseq_bulk_analysis_viral_load("both", both_placebo_counts, both_placebo_metadata,
                                                                              "2_D_minus_2", "HIGH", "LOW", paste0(bulk_results_dir, "placebo_period_2_D_minus_2_high_vs_low/"), "2_D_minus_2")
# 2 D minus 1
# 167/157/155/151/178
# Note: If I don't control for age and sex, I find nothing.
placebo_period_2_D_minus_1_high_vs_low_results <- run_deseq_bulk_analysis_viral_load("both", both_placebo_counts, both_placebo_metadata,
                                                                              "2_D_minus_1", "HIGH", "LOW", paste0(bulk_results_dir, "placebo_period_2_D_minus_1_high_vs_low/"), "2_D_minus_1")
# 2 D2
# 139/138/137/134/139
placebo_period_2_D2_high_vs_low_results <- run_deseq_bulk_analysis_viral_load("both", both_placebo_counts, both_placebo_metadata,
                                                                                     "2_D2", "HIGH", "LOW", paste0(bulk_results_dir, "placebo_period_2_D2_high_vs_low/"), "2_D2")
# 2 D5
# Note: If I don't control for age and sex, I still find a lot of DEGs.
# https://hb.flatironinstitute.org/module/overview/?body_tag=9a43682013acc9404b08417d64f9e60f2436bfba all upregulated
# 883/269/155/88/1653
placebo_period_2_D5_high_vs_low_results <- run_deseq_bulk_analysis_viral_load("both", both_placebo_counts, both_placebo_metadata,
                                                                              "2_D5", "HIGH", "LOW", paste0(bulk_results_dir, "placebo_period_2_D5_high_vs_low/"), "2_D5")
# 2 D8
# https://hb.flatironinstitute.org/module/overview/?body_tag=e7397a8502c94a074b03dedb68c473ef8912a8e9 all upregulated
# 298/207/196/187/365
placebo_period_2_D8_high_vs_low_results <- run_deseq_bulk_analysis_viral_load("both", both_placebo_counts, both_placebo_metadata,
                                                                              "2_D8", "HIGH", "LOW", paste0(bulk_results_dir, "placebo_period_2_D8_high_vs_low/"), "2_D8")
# 2 D28
# 193/180/172/170/209
placebo_period_2_D28_high_vs_low_results <- run_deseq_bulk_analysis_viral_load("both", both_placebo_counts, both_placebo_metadata,
                                                                              "2_D28", "HIGH", "LOW", paste0(bulk_results_dir, "placebo_period_2_D8_high_vs_low/"), "2_D28")

# VACCINATED
#### PERIOD 1 HIGH T CELL RESPONSE ####
# 6234/1686/678/163/7799 DEGs
high_vaccinated_period_1_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", high_vaccinated_counts, high_vaccinated_metadata,
                                                                                     "1_D2", "1_D_minus_1", paste0(bulk_results_dir, "high_vaccinated_period_1_D2_vs_D_minus_1/"), "high")
raw_high_vaccinated_period_1_D2_vs_D_minus_1_results <- high_vaccinated_period_1_D2_vs_D_minus_1_results[[5]]
# 59/15/2/0/96 DEGs
high_vaccinated_period_1_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", high_vaccinated_counts, high_vaccinated_metadata,
                                                                                     "1_D8", "1_D_minus_1", paste0(bulk_results_dir, "high_vaccinated_period_1_D8_vs_D_minus_1/"), "high")
raw_high_vaccinated_period_1_D8_vs_D_minus_1_results <- high_vaccinated_period_1_D8_vs_D_minus_1_results[[5]]
# 0/0/0/0/0 DEGs
high_vaccinated_period_1_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", high_vaccinated_counts, high_vaccinated_metadata,
                                                                                      "1_D28", "1_D_minus_1", paste0(bulk_results_dir, "high_vaccinated_period_1_D28_vs_D_minus_1/"), "high")
raw_high_vaccinated_period_1_D28_vs_D_minus_1_results <- high_vaccinated_period_1_D28_vs_D_minus_1_results[[5]]
#### PERIOD 1 LOW VIRAL LOAD ####
# 2280/312/129/18/3854 DEGs
low_vaccinated_period_1_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", low_vaccinated_counts, low_vaccinated_metadata,
                                                                                    "1_D2", "1_D_minus_1", paste0(bulk_results_dir, "low_vaccinated_period_1_D2_vs_D_minus_1/"), "low")
raw_low_vaccinated_period_1_D2_vs_D_minus_1_results <- low_vaccinated_period_1_D2_vs_D_minus_1_results[[5]]
# 70/0/0/0/521 DEGs
low_vaccinated_period_1_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", low_vaccinated_counts, low_vaccinated_metadata,
                                                                                    "1_D8", "1_D_minus_1", paste0(bulk_results_dir, "low_vaccinated_period_1_D8_vs_D_minus_1/"), "low")
raw_low_vaccinated_period_1_D8_vs_D_minus_1_results <- low_vaccinated_period_1_D8_vs_D_minus_1_results[[5]]
# 2/0/0/0/3 DEGs
low_vaccinated_period_1_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", low_vaccinated_counts, low_vaccinated_metadata,
                                                                                     "1_D28", "1_D_minus_1", paste0(bulk_results_dir, "low_vaccinated_period_1_D28_vs_D_minus_1/"), "low")
raw_low_vaccinated_period_1_D28_vs_D_minus_1_results <- low_vaccinated_period_1_D28_vs_D_minus_1_results[[5]]

#### PERIOD 2 HIGH VIRAL LOAD ####
# 2 D minus 2 vs 2 D minus 1 - should be virtually zero unless some weird stuff happened between blood draws
# 294/0/0/0/1225 DEGs found - seems weird that so many are found
high_vaccinated_period_2_D_minus_1_vs_D_minus_2_results <- run_deseq_bulk_analysis_time_series("vaccinated", high_vaccinated_counts, high_vaccinated_metadata,
                                                                                            "2_D_minus_1", "2_D_minus_2", paste0(bulk_results_dir, "high_vaccinated_period_2_D_minus_1_vs_D_minus_2/"), "high")
raw_high_vaccinated_period_2_D_minus_1_vs_D_minus_2_results <- high_vaccinated_period_2_D_minus_1_vs_D_minus_2_results[[5]]
# 0.585 is good, and -0.585 is good
#collected_fmd_high_vaccinated_period_2_D_minus_1_vs_D_minus_2 <- run_fmd_on_flu_data(raw_high_vaccinated_period_2_D_minus_1_vs_D_minus_2_results)

# 2 D2 vs 2 D minus 1 - 0/0/0/0/2 DEGs
high_vaccinated_period_2_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", high_vaccinated_counts, high_vaccinated_metadata,
                                                                                     "2_D2", "2_D_minus_1", paste0(bulk_results_dir, "high_vaccinated_period_2_D2_vs_D_minus_1/"), "high")
raw_high_vaccinated_period_2_D2_vs_D_minus_1_results <- high_vaccinated_period_2_D2_vs_D_minus_1_results[[5]]
# 0.585 is good, and -0.585 is good
#collected_fmd_high_vaccinated_period_2_D2_vs_D_minus_1 <- run_fmd_on_flu_data(raw_high_vaccinated_period_2_D2_vs_D_minus_1_results)

# 2 D5 vs 2 D minus 1 - 704/20/0/0/1588 DEGs
high_vaccinated_period_2_D5_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", high_vaccinated_counts, high_vaccinated_metadata,
                                                                                     "2_D5", "2_D_minus_1", paste0(bulk_results_dir, "high_vaccinated_period_2_D5_vs_D_minus_1/"), "high")
raw_high_vaccinated_period_2_D5_vs_D_minus_1_results <- high_vaccinated_period_2_D5_vs_D_minus_1_results[[5]]
# 0.585 or 1 is good, and -0.585 or -1 is good 
#collected_fmd_high_vaccinated_period_2_D5_vs_D_minus_1 <- run_fmd_on_flu_data(raw_high_vaccinated_period_2_D5_vs_D_minus_1_results, c(0.585, 1, 2))

# 2 D8 vs 2 D minus 1 - 1359/85/1/0/2441 DEGs
high_vaccinated_period_2_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", high_vaccinated_counts, high_vaccinated_metadata,
                                                                                     "2_D8", "2_D_minus_1", paste0(bulk_results_dir, "high_vaccinated_period_2_D8_vs_D_minus_1/"), "high")
raw_high_vaccinated_period_2_D8_vs_D_minus_1_results <- high_vaccinated_period_2_D8_vs_D_minus_1_results[[5]]
# ???
# collected_fmd_high_vaccinated_period_2_D8_vs_D_minus_1 <- run_fmd_on_flu_data(raw_high_vaccinated_period_2_D8_vs_D_minus_1_results, c(0.585, 1, 2))
# 2 D28 vs 2 D minus 1 - 6/4/4/3/21 DEGs
high_vaccinated_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", high_vaccinated_counts, high_vaccinated_metadata,
                                                                                      "2_D28", "2_D_minus_1", paste0(bulk_results_dir, "high_vaccinated_period_2_D28_vs_D_minus_1/"), "high")
raw_high_vaccinated_period_2_D28_vs_D_minus_1_results <- high_vaccinated_period_2_D28_vs_D_minus_1_results[[5]]
# 0.3 is good, and -0.3 is good
#collected_fmd_high_vaccinated_period_2_D28_vs_D_minus_1 <- run_fmd_on_flu_data(raw_high_vaccinated_period_2_D28_vs_D_minus_1_results)

#### PERIOD 2 LOW VIRAL LOAD ####
# 2 D minus 2 vs 2 D minus 1 - should be virtually zero unless some weird stuff happened between blood draws
# 1/0/0/0/9 DEGs found - seems weird that so many are found
low_vaccinated_period_2_D_minus_1_vs_D_minus_2_results <- run_deseq_bulk_analysis_time_series("vaccinated", low_vaccinated_counts, low_vaccinated_metadata,
                                                                                               "2_D_minus_1", "2_D_minus_2", paste0(bulk_results_dir, "low_vaccinated_period_2_D_minus_1_vs_D_minus_2/"), "low")
raw_low_vaccinated_period_2_D_minus_1_vs_D_minus_2_results <- low_vaccinated_period_2_D_minus_1_vs_D_minus_2_results[[5]]
# 0.585 is good, and -0.585 is good
#collected_fmd_low_vaccinated_period_2_D_minus_1_vs_D_minus_2 <- run_fmd_on_flu_data(raw_low_vaccinated_period_2_D_minus_1_vs_D_minus_2_results)

# 2 D2 vs 2 D minus 1 - 2/0/0/0/2 DEGs
low_vaccinated_period_2_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", low_vaccinated_counts, low_vaccinated_metadata,
                                                                                        "2_D2", "2_D_minus_1", paste0(bulk_results_dir, "low_vaccinated_period_2_D2_vs_D_minus_1/"), "low")
raw_low_vaccinated_period_2_D2_vs_D_minus_1_results <- low_vaccinated_period_2_D2_vs_D_minus_1_results[[5]]
# 0.585 is good, and -0.585 is good
#collected_fmd_low_vaccinated_period_2_D2_vs_D_minus_1 <- run_fmd_on_flu_data(raw_low_vaccinated_period_2_D2_vs_D_minus_1_results)

# 2 D5 vs 2 D minus 1 - 101/20/0/0/366 DEGs
low_vaccinated_period_2_D5_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", low_vaccinated_counts, low_vaccinated_metadata,
                                                                                        "2_D5", "2_D_minus_1", paste0(bulk_results_dir, "low_vaccinated_period_2_D5_vs_D_minus_1/"), "low")
raw_low_vaccinated_period_2_D5_vs_D_minus_1_results <- low_vaccinated_period_2_D5_vs_D_minus_1_results[[5]]
# 0.585 or 1 is good, and -0.585 or -1 is good 
#collected_fmd_low_vaccinated_period_2_D5_vs_D_minus_1 <- run_fmd_on_flu_data(raw_low_vaccinated_period_2_D5_vs_D_minus_1_results, c(0.585, 1, 2))

# 2 D8 vs 2 D minus 1 - 4/0/0/0/15 DEGs
low_vaccinated_period_2_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", low_vaccinated_counts, low_vaccinated_metadata,
                                                                                        "2_D8", "2_D_minus_1", paste0(bulk_results_dir, "low_vaccinated_period_2_D8_vs_D_minus_1/"), "low")
#raw_low_vaccinated_period_2_D8_vs_D_minus_1_results <- low_vaccinated_period_2_D8_vs_D_minus_1_results[[5]]
# ???
# collected_fmd_low_vaccinated_period_2_D8_vs_D_minus_1 <- run_fmd_on_flu_data(raw_low_vaccinated_period_2_D8_vs_D_minus_1_results, c(0.585, 1, 2))
# 2 D28 vs 2 D minus 1 - 106/0/0/0/800 DEGs
low_vaccinated_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", low_vaccinated_counts, low_vaccinated_metadata,
                                                                                         "2_D28", "2_D_minus_1", paste0(bulk_results_dir, "low_vaccinated_period_2_D28_vs_D_minus_1/"), "low")
raw_low_vaccinated_period_2_D28_vs_D_minus_1_results <- low_vaccinated_period_2_D28_vs_D_minus_1_results[[5]]
# 0.3 is good, and -0.3 is good
#collected_fmd_low_vaccinated_period_2_D28_vs_D_minus_1 <- run_fmd_on_flu_data(raw_low_vaccinated_period_2_D28_vs_D_minus_1_results)


















# LRT test for period 2 (HVL)
# 6431 DEGs found (note that LRT doesn't have a fold change threshold)
high_placebo_period_2_LRT_metadata <- high_placebo_metadata[high_placebo_metadata$time_point == "2_D28" | high_placebo_metadata$time_point == "2_D8" | 
                                                              high_placebo_metadata$time_point == "2_D5" | high_placebo_metadata$time_point == "2_D2" |
                                                              high_placebo_metadata$time_point == "2_D_minus_1",]
high_placebo_period_2_LRT_metadata <- high_placebo_period_2_LRT_metadata[high_placebo_period_2_LRT_metadata$subject_id 
                                                                         %in% names(table(high_placebo_period_2_LRT_metadata$subject_id)
                                                                                    [table(high_placebo_period_2_LRT_metadata$subject_id) == 5]),]
high_placebo_period_2_LRT_counts <- high_placebo_counts[rownames(high_placebo_period_2_LRT_metadata)]
high_placebo_period_2_LRT_analysis <- DESeqDataSetFromMatrix(countData = high_placebo_period_2_LRT_counts,
                                                             colData = high_placebo_period_2_LRT_metadata,
                                                             design = ~ subject_id + Monocytes + Neutrophils + time_point)
high_placebo_period_2_LRT_analysis <- DESeq(high_placebo_period_2_LRT_analysis, test = "LRT", reduced = ~ subject_id + Monocytes + Neutrophils)
high_placebo_period_2_LRT_analysis_results <- results(high_placebo_period_2_LRT_analysis, alpha = 0.05)
high_placebo_period_2_LRT_analysis_results <- high_placebo_period_2_LRT_analysis_results[order(high_placebo_period_2_LRT_analysis_results$padj),]
high_placebo_period_2_LRT_analysis_results <- subset(high_placebo_period_2_LRT_analysis_results, padj < 0.05)
# Write top 2000 genes to file
high_placebo_period_2_LRT_analysis_2000_topGenes <- rownames(high_placebo_period_2_LRT_analysis_results)[1:2000]
if (!dir.exists(paste0(bulk_data_dir, "high_placebo_period_2_LRT"))) {dir.create(paste0(bulk_data_dir, "high_placebo_period_2_LRT")) }
write.table(high_placebo_period_2_LRT_analysis_2000_topGenes, file = paste0(bulk_data_dir, "high_placebo_period_2_LRT/LRT_top_200_genes.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
# Plot heatmap with top 20 genes
high_placebo_period_2_LRT_analysis_betas <- coef(high_placebo_period_2_LRT_analysis)
high_placebo_period_2_LRT_analysis_results_for_plotting <- results(high_placebo_period_2_LRT_analysis, alpha = 0.05)
high_placebo_period_2_LRT_analysis_20_topGenes <- head(order(high_placebo_period_2_LRT_analysis_results_for_plotting$padj),20)
high_placebo_period_2_LRT_analysis_20_mat <- high_placebo_period_2_LRT_analysis_betas[high_placebo_period_2_LRT_analysis_20_topGenes, -c(1:13)]
high_placebo_period_2_LRT_analysis_thr <- 6
high_placebo_period_2_LRT_analysis_20_mat[high_placebo_period_2_LRT_analysis_20_mat < -high_placebo_period_2_LRT_analysis_thr] <- -high_placebo_period_2_LRT_analysis_thr
high_placebo_period_2_LRT_analysis_20_mat[high_placebo_period_2_LRT_analysis_20_mat > high_placebo_period_2_LRT_analysis_thr] <- high_placebo_period_2_LRT_analysis_thr
colnames(high_placebo_period_2_LRT_analysis_20_mat) <- c("Day 2 vs Day -1", "Day 5 vs Day -1", "Day 8 vs Day -1", "Day 28 vs Day -1")
pheatmap(high_placebo_period_2_LRT_analysis_20_mat, breaks=seq(from=-high_placebo_period_2_LRT_analysis_thr, to=high_placebo_period_2_LRT_analysis_thr, length=101),
         cluster_col=FALSE, fontsize_col=14, filename = paste0(bulk_data_dir, "high_placebo_period_2_LRT/high_placebo_period_2_LRT_analysis_top_20_genes.png"))
# Plot heatmap with top 200 genes (and clean up heatmap plot by removing row names and clustering)
high_placebo_period_2_LRT_analysis_200_topGenes <- head(order(high_placebo_period_2_LRT_analysis_results_for_plotting$padj),200)
high_placebo_period_2_LRT_analysis_200_mat <- high_placebo_period_2_LRT_analysis_betas[high_placebo_period_2_LRT_analysis_200_topGenes, -c(1:13)]
colnames(high_placebo_period_2_LRT_analysis_200_mat) <- c("Day 2 vs Day -1", "Day 5 vs Day -1", "Day 8 vs Day -1", "Day 28 vs Day -1")
high_placebo_period_2_LRT_analysis_200_mat[high_placebo_period_2_LRT_analysis_200_mat < -high_placebo_period_2_LRT_analysis_thr] <- -high_placebo_period_2_LRT_analysis_thr
high_placebo_period_2_LRT_analysis_200_mat[high_placebo_period_2_LRT_analysis_200_mat > high_placebo_period_2_LRT_analysis_thr] <- high_placebo_period_2_LRT_analysis_thr
pheatmap(high_placebo_period_2_LRT_analysis_200_mat, breaks=seq(from=-high_placebo_period_2_LRT_analysis_thr, to=high_placebo_period_2_LRT_analysis_thr, length=101),
         cluster_col=FALSE, fontsize_col=14, show_rownames = FALSE, cluster_row = FALSE, filename = paste0(bulk_data_dir, "high_placebo_period_2_LRT/high_placebo_period_2_LRT_analysis_top_200_genes.png"))


# LRT test for period 2 (LVL)
# 0 DEGs found (note that LRT doesn't have a fold change threshold)
low_placebo_period_2_LRT_metadata <- low_placebo_metadata[low_placebo_metadata$time_point == "2_D28" | low_placebo_metadata$time_point == "2_D8" | 
                                                            low_placebo_metadata$time_point == "2_D5" | low_placebo_metadata$time_point == "2_D2" |
                                                            low_placebo_metadata$time_point == "2_D_minus_1",]
low_placebo_period_2_LRT_metadata <- low_placebo_period_2_LRT_metadata[low_placebo_period_2_LRT_metadata$subject_id 
                                                                       %in% names(table(low_placebo_period_2_LRT_metadata$subject_id)
                                                                                  [table(low_placebo_period_2_LRT_metadata$subject_id) == 5]),]
low_placebo_period_2_LRT_counts <- low_placebo_counts[rownames(low_placebo_period_2_LRT_metadata)]
low_placebo_period_2_LRT_analysis <- DESeqDataSetFromMatrix(countData = low_placebo_period_2_LRT_counts,
                                                            colData = low_placebo_period_2_LRT_metadata,
                                                            design = ~ subject_id + time_point)
low_placebo_period_2_LRT_analysis <- DESeq(low_placebo_period_2_LRT_analysis, test = "LRT", reduced = ~ subject_id)
low_placebo_period_2_LRT_analysis_results <- results(low_placebo_period_2_LRT_analysis, alpha = 0.05)
low_placebo_period_2_LRT_analysis_results <- low_placebo_period_2_LRT_analysis_results[order(low_placebo_period_2_LRT_analysis_results$padj),]
low_placebo_period_2_LRT_analysis_results <- subset(low_placebo_period_2_LRT_analysis_results, padj < 0.05)