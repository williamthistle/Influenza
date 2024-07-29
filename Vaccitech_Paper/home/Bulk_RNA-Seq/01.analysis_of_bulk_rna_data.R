### THIS FUNCTION RUNS DESEQ2 ON ALL BULK RNA-SEQ GROUPINGS
### WE CAN LOOK AT HIGH EFFECT SIZE GENES INDIVIDUALLY (NO NEED TO DO FC THRESHOLD?)
### IT ALSO RUNS FUNCTIONAL MODULE DISCOVERY ON DIFFERENT FOLD CHANGE THRESHOLDS
### FOR EACH BULK RNA-SEQ GROUPING
### THESE HB PLOTS (OR DOWNSTREAM REACTOME PLOTS) MAY GO IN SUPPLEMENTARY FIGURES
### IN MAIN TEXT, DISCUSS EACH DAY AND WHAT WE SEE (HIGH EFFECT SIZE GENES, MODULES, PATHWAYS)

# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# load(paste0(bulk_rna_results_dir, "bulk_RNA_analysis.RData"))

################ PLACEBO FULL TIME SERIES ################

#### PERIOD 1 HIGH VIRAL LOAD ####
# We expect ~0 DEGs because placebo was used (no actual vaccination)
# 0/0/0/0/0/0/0 DEGs
hvl_full_time_series_placebo_period_1_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", hvl_full_time_series_placebo_counts, hvl_full_time_series_placebo_metadata,
                                                                                "1_D2", "1_D_minus_1", paste0(bulk_rna_results_dir, "hvl_bulk_placebo_period_1_D2_vs_D_minus_1/"), "high")
# 5/0/0/0/0/0/0 DEGs
hvl_full_time_series_placebo_period_1_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", hvl_full_time_series_placebo_counts, hvl_full_time_series_placebo_metadata,
                                                                                "1_D8", "1_D_minus_1", paste0(bulk_rna_results_dir, "hvl_bulk_period_1_D8_vs_D_minus_1/"), "high")
# 103/0/0/0/0/0/0 DEGs
hvl_full_time_series_placebo_period_1_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", hvl_full_time_series_placebo_counts, hvl_full_time_series_placebo_metadata,
                                                                                "1_D28", "1_D_minus_1", paste0(bulk_rna_results_dir, "hvl_bulk_period_1_D28_vs_D_minus_1/"), "high")
#### PERIOD 1 MODERATE VIRAL LOAD ####
# No samples

#### PERIOD 1 LOW VIRAL LOAD ####
# We expect ~0 DEGs because placebo was used (no actual vaccination)
# 4/4/4/4/4/4/4 DEGs
lvl_full_time_series_placebo_period_1_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", lvl_full_time_series_placebo_counts, lvl_full_time_series_placebo_metadata,
                                                                         "1_D2", "1_D_minus_1", paste0(bulk_rna_results_dir, "lvl_bulk_placebo_period_1_D2_vs_D_minus_1/"), "low")
# 24/0/0/0/0/0/0 DEGs
lvl_full_time_series_placebo_period_1_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", lvl_full_time_series_placebo_counts, lvl_full_time_series_placebo_metadata,
                                                                         "1_D8", "1_D_minus_1", paste0(bulk_rna_results_dir, "lvl_bulk_placebo_period_1_D8_vs_D_minus_1/"), "low")
# 205/62/5/0/0/0/0 DEGs
lvl_full_time_series_placebo_period_1_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", lvl_full_time_series_placebo_counts, lvl_full_time_series_placebo_metadata,
                                                                          "1_D28", "1_D_minus_1", paste0(bulk_rna_results_dir, "lvl_bulk_placebo_period_1_D28_vs_D_minus_1/"), "low")
#### PERIOD 2 HIGH VIRAL LOAD ####
# 2 D minus 1 vs 2 D minus 2 - should be virtually zero unless some weird stuff happened between blood draws
# 458/114/37/1/0/0/0 DEGs found - seems weird that so many are found
# Could there be something that happened shortly before infection that made subjects more susceptible to getting sick?
# Or is there something about these subjects that make them more predisposed to getting sick?
# Could it be a single patient (or two) messing everything up?
# Is it just inflammation from getting blood drawn?
# If so, why would this signal exist for HVL and not LVL?
hvl_full_time_series_placebo_period_2_D_minus_1_vs_D_minus_2_results <- run_deseq_bulk_analysis_time_series("placebo", hvl_full_time_series_placebo_counts, hvl_full_time_series_placebo_metadata,
                                                                                "2_D_minus_1", "2_D_minus_2", paste0(bulk_rna_results_dir, "hvl_bulk_placebo_period_2_D_minus_1_vs_D_minus_2/"), "high")
# 2 D2 vs 2 D minus 1 - 2/0/0/0/0/0/0 DEGs
hvl_full_time_series_placebo_period_2_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", hvl_full_time_series_placebo_counts, hvl_full_time_series_placebo_metadata,
                                                                         "2_D2", "2_D_minus_1", paste0(bulk_rna_results_dir, "hvl_bulk_placebo_period_2_D2_vs_D_minus_1/"), "high")
# 2 D5 vs 2 D minus 1 - 3799/2337/1501/1024/403/178/41 DEGs
hvl_full_time_series_placebo_period_2_D5_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", hvl_full_time_series_placebo_counts, hvl_full_time_series_placebo_metadata,
                                                                         "2_D5", "2_D_minus_1", paste0(bulk_rna_results_dir, "hvl_bulk_placebo_period_2_D5_vs_D_minus_1/"), "high")
# 2 D8 vs 2 D minus 1 - 3587/2056/1268/816/267/66/8 DEGs
hvl_full_time_series_placebo_period_2_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", hvl_full_time_series_placebo_counts, hvl_full_time_series_placebo_metadata,
                                                                         "2_D8", "2_D_minus_1", paste0(bulk_rna_results_dir, "hvl_bulk_placebo_period_2_D8_vs_D_minus_1/"), "high")
# 2 D28 vs 2 D minus 1 - 77/0/0/0/0/0/0 DEGs
hvl_full_time_series_placebo_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", hvl_full_time_series_placebo_counts, hvl_full_time_series_placebo_metadata,
                                                                          "2_D28", "2_D_minus_1", paste0(bulk_rna_results_dir, "hvl_bulk_placebo_period_2_D28_vs_D_minus_1/"), "high")

# Run LRT on period 2 (Day -1, Day 2, Day 5, Day 8, Day 28)
hvl_full_time_series_placebo_lrt <- run_deseq2_LRT(hvl_full_time_series_placebo_counts, hvl_full_time_series_placebo_metadata)

#### PERIOD 2 MODERATE VIRAL LOAD ####
# No samples

#### PERIOD 2 LOW VIRAL LOAD ####
# 2 D minus 2 vs 2 D minus 1 - should be virtually zero unless some weird stuff happened between blood draws
# 2/2/1/1/1/1/1 DEGs found
lvl_full_time_series_placebo_period_2_D_minus_1_vs_D_minus_2_results <- run_deseq_bulk_analysis_time_series("placebo", lvl_full_time_series_placebo_counts, lvl_full_time_series_placebo_metadata,
                                                                               "2_D_minus_1", "2_D_minus_2", paste0(bulk_rna_results_dir, "lvl_bulk_placebo_period_2_D_minus_1_vs_D_minus_2/"), "low")
# 2 D2 vs 2 D minus 1 - 2/1/1/1/1/1/1 DEGs
lvl_full_time_series_placebo_period_2_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", lvl_full_time_series_placebo_counts, lvl_full_time_series_placebo_metadata,
                                                                        "2_D2", "2_D_minus_1", paste0(bulk_rna_results_dir, "lvl_bulk_placebo_period_2_D2_vs_D_minus_1/"), "low")
# 2 D5 vs 2 D minus 1 - 4/4/4/4/4/4/4 DEGs
lvl_full_time_series_placebo_period_2_D5_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", lvl_full_time_series_placebo_counts, lvl_full_time_series_placebo_metadata,
                                                                        "2_D5", "2_D_minus_1", paste0(bulk_rna_results_dir, "lvl_bulk_placebo_period_2_D5_vs_D_minus_1/"), "low")
# 2 D8 vs 2 D minus 1 - 2/2/2/2/2/2/2 DEGs
lvl_full_time_series_placebo_period_2_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", lvl_full_time_series_placebo_counts, lvl_full_time_series_placebo_metadata,
                                                                        "2_D8", "2_D_minus_1", paste0(bulk_rna_results_dir, "lvl_bulk_placebo_period_2_D8_vs_D_minus_1/"), "low")
# 2 D28 vs 2 D minus 1 - 192/14/2/2/1/1/1 DEGs
lvl_full_time_series_placebo_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", lvl_full_time_series_placebo_counts, lvl_full_time_series_placebo_metadata,
                                                                         "2_D28", "2_D_minus_1", paste0(bulk_rna_results_dir, "lvl_bulk_placebo_period_2_D28_vs_D_minus_1/"), "low")

################ PLACEBO MAX NUMBER OF SAMPLES ################

#### PERIOD 2 HIGH VIRAL LOAD ####
# 2 D8 vs 2 D minus 1 - 3587/2056/1268/816/267/66/8 DEGs
hvl_placebo_period_2_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", hvl_placebo_counts, hvl_placebo_metadata,
                                                                                    "2_D8", "2_D_minus_1", paste0(bulk_rna_results_dir, "hvl_bulk_placebo_period_2_D8_vs_D_minus_1/"), "high")
# 2 D28 vs 2 D minus 1 - 77/0/0/0/0/0/0 DEGs
hvl_placebo_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", hvl_placebo_counts, hvl_placebo_metadata,
                                                                                     "2_D28", "2_D_minus_1", paste0(bulk_rna_results_dir, "hvl_bulk_placebo_period_2_D28_vs_D_minus_1/"), "high")
#### PERIOD 2 MODERATE VIRAL LOAD ####
# 2 D8 vs 2 D minus 1 - 77/0/0/0/0/0/0 DEGs
mvl_placebo_period_2_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", mvl_placebo_counts, mvl_placebo_metadata,
                                                                                    "2_D8", "2_D_minus_1", paste0(bulk_rna_results_dir, "hvl_bulk_placebo_period_2_D8_vs_D_minus_1/"), "moderate")
# 2 D28 vs 2 D minus 1 - 77/0/0/0/0/0/0 DEGs
mvl_placebo_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", mvl_placebo_counts, mvl_placebo_metadata,
                                                                                     "2_D28", "2_D_minus_1", paste0(bulk_rna_results_dir, "hvl_bulk_placebo_period_2_D28_vs_D_minus_1/"), "moderate")
#### PERIOD 2 LOW VIRAL LOAD ####
# 2 D8 vs 2 D minus 1 - 2/2/2/2/2/2/2 DEGs
lvl_placebo_period_2_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", lvl_placebo_counts, lvl_placebo_metadata,
                                                                                    "2_D8", "2_D_minus_1", paste0(bulk_rna_results_dir, "lvl_bulk_placebo_period_2_D8_vs_D_minus_1/"), "low")
# 2 D28 vs 2 D minus 1 - 192/14/2/2/1/1/1 DEGs
lvl_placebo_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", lvl_placebo_counts, lvl_placebo_metadata,
                                                                                     "2_D28", "2_D_minus_1", paste0(bulk_rna_results_dir, "lvl_bulk_placebo_period_2_D28_vs_D_minus_1/"), "low")

################ VACCINATED ################

#### PERIOD 1 ALL VIRAL LOAD
all_vaccinated_period_1_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", vaccinated_counts, vaccinated_metadata,
                                                                                       "1_D8", "1_D_minus_1", paste0(bulk_rna_results_dir, "hvl_bulk_vaccinated_period_2_D8_vs_D_minus_1/"), "high")
#### PERIOD 1 HIGH VIRAL LOAD ####
# D8 vs D minus 1 - 
hvl_vaccinated_period_1_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", hvl_vaccinated_counts, hvl_vaccinated_metadata,
                                                                                       "1_D8", "1_D_minus_1", paste0(bulk_rna_results_dir, "hvl_bulk_vaccinated_period_2_D8_vs_D_minus_1/"), "high")
#### PERIOD 1 MODERATE VIRAL LOAD ####
# D8 vs D-1
mvl_vaccinated_period_1_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", mvl_vaccinated_counts, mvl_vaccinated_metadata,
                                                                                       "1_D8", "1_D_minus_1", paste0(bulk_rna_results_dir, "mvl_bulk_vaccinated_period_2_D8_vs_D_minus_1/"), "moderate")
#### PERIOD 1 LOW VIRAL LOAD ####
# D8 vs D-1
lvl_vaccinated_period_1_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", lvl_vaccinated_counts, lvl_vaccinated_metadata,
                                                                                       "1_D8", "1_D_minus_1", paste0(bulk_rna_results_dir, "lvl_bulk_vaccinated_period_2_D8_vs_D_minus_1/"), "low")
#### PERIOD 2 HIGH VIRAL LOAD ####
# D minus 1 vs D minus 2
hvl_vaccinated_period_2_D_minus_1_vs_D_minus_2_results <- run_deseq_bulk_analysis_time_series("vaccinated", hvl_vaccinated_counts, hvl_vaccinated_metadata,
                                                                                       "2_D_minus_1", "2_D_minus_2", paste0(bulk_rna_results_dir, "hvl_bulk_vaccinated_period_2_D8_vs_D_minus_1/"), "high")
# D8 vs D-1
hvl_vaccinated_period_2_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", hvl_vaccinated_counts, hvl_vaccinated_metadata,
                                                                                       "2_D8", "2_D_minus_1", paste0(bulk_rna_results_dir, "hvl_bulk_vaccinated_period_2_D8_vs_D_minus_1/"), "high")
# D28 vs D-1
hvl_vaccinated_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", hvl_vaccinated_counts, hvl_vaccinated_metadata,
                                                                                        "2_D28", "2_D_minus_1", paste0(bulk_rna_results_dir, "hvl_bulk_vaccinated_period_2_D28_vs_D_minus_1/"), "high")
#### PERIOD 2 HIGH VIRAL LOAD (FULL TIME SERIES) ####
# 2 D minus 1 vs 2 D minus 2
hvl_full_time_series_vaccinated_period_2_D_minus_1_vs_D_minus_2_results <- run_deseq_bulk_analysis_time_series("vaccinated", hvl_full_time_series_vaccinated_counts, hvl_full_time_series_vaccinated_metadata,
                                                                                                            "2_D_minus_1", "2_D_minus_2", paste0(bulk_rna_results_dir, "hvl_bulk_vaccinated_period_2_D_minus_1_vs_D_minus_2/"), "high")
# 2 D2 vs 2 D minus 1 
hvl_full_time_series_vaccinated_period_2_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", hvl_full_time_series_vaccinated_counts, hvl_full_time_series_vaccinated_metadata,
                                                                                                     "2_D2", "2_D_minus_1", paste0(bulk_rna_results_dir, "hvl_bulk_vaccinated_period_2_D2_vs_D_minus_1/"), "high")
# 2 D5 vs 2 D minus 1 
hvl_full_time_series_vaccinated_period_2_D5_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", hvl_full_time_series_vaccinated_counts, hvl_full_time_series_vaccinated_metadata,
                                                                                                     "2_D5", "2_D_minus_1", paste0(bulk_rna_results_dir, "hvl_bulk_vaccinated_period_2_D5_vs_D_minus_1/"), "high")
# 2 D8 vs 2 D minus 1 
hvl_full_time_series_vaccinated_period_2_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", hvl_full_time_series_vaccinated_counts, hvl_full_time_series_vaccinated_metadata,
                                                                                                     "2_D8", "2_D_minus_1", paste0(bulk_rna_results_dir, "hvl_bulk_vaccinated_period_2_D8_vs_D_minus_1/"), "high")
# 2 D28 vs 2 D minus 1
hvl_full_time_series_vaccinated_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", hvl_full_time_series_vaccinated_counts, hvl_full_time_series_vaccinated_metadata,
                                                                                                      "2_D28", "2_D_minus_1", paste0(bulk_rna_results_dir, "hvl_bulk_vaccinated_period_2_D28_vs_D_minus_1/"), "high")

#### PERIOD 2 MODERATE VIRAL LOAD ####
# D8 vs D-1
mvl_vaccinated_period_2_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", mvl_vaccinated_counts, mvl_vaccinated_metadata,
                                                                                       "2_D8", "2_D_minus_1", paste0(bulk_rna_results_dir, "mvl_bulk_vaccinated_period_2_D8_vs_D_minus_1/"), "moderate")
# D28 vs D-1
mvl_vaccinated_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", mvl_vaccinated_counts, mvl_vaccinated_metadata,
                                                                                        "2_D28", "2_D_minus_1", paste0(bulk_rna_results_dir, "mvl_bulk_vaccinated_period_2_D28_vs_D_minus_1/"), "moderate")
#### PERIOD 2 LOW VIRAL LOAD ####
# D minus 1 vs D minus 2
lvl_vaccinated_period_2_D_minus_1_vs_D_minus_2_results <- run_deseq_bulk_analysis_time_series("vaccinated", lvl_vaccinated_counts, lvl_vaccinated_metadata,
                                                                                              "2_D_minus_1", "2_D_minus_2", paste0(bulk_rna_results_dir, "hvl_bulk_vaccinated_period_2_D8_vs_D_minus_1/"), "high")
# D8 vs D-1
lvl_vaccinated_period_2_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", lvl_vaccinated_counts, lvl_vaccinated_metadata,
                                                                                       "2_D8", "2_D_minus_1", paste0(bulk_rna_results_dir, "lvl_bulk_vaccinated_period_2_D8_vs_D_minus_1/"), "low")
# D28 vs D-1
lvl_vaccinated_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", lvl_vaccinated_counts, lvl_vaccinated_metadata,
                                                                                        "2_D28", "2_D_minus_1", paste0(bulk_rna_results_dir, "lvl_bulk_vaccinated_period_2_D28_vs_D_minus_1/"), "low")
# Save results to disk
save.image(paste0(bulk_rna_results_dir, "bulk_RNA_analysis.RData"))

# Compare fold changes for overlapping genes between vaccinated and placebo for D28
overlapping_genes <- intersect(rownames(raw_hvl_vaccinated_period_2_D28_vs_D_minus_1_results), rownames(raw_hvl_placebo_period_2_D28_vs_D_minus_1_results))

compare_placebo_D28_df <- raw_hvl_placebo_period_2_D28_vs_D_minus_1_results[rownames(raw_hvl_placebo_period_2_D28_vs_D_minus_1_results) %in% overlapping_genes,]
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

# More general framework for correlation
# Weird: using vac-only DEGS for HVL and testing for correlation in LVL results in 0.56 - higher than expected?
# Kind of suggests that these HVL genes are consistent even in LVL
# Note 2: log approach seems to result in higher correlation between placebo and vac (HVL) for Day 28 with cell type proportion
# correction
first <- hvl_placebo_period_2_D8_vs_D_minus_1_results[[1]]
second <- hvl_vaccinated_period_2_D8_vs_D_minus_1_results[[1]]
unfiltered_first <- unfiltered_hvl_placebo_period_2_D8_vs_D_minus_1_results[[1]]
unfiltered_second <- unfiltered_hvl_vaccinated_period_2_D8_vs_D_minus_1_results[[1]]

# 1: Check overlapping DEGs
overlapping_genes <- intersect(rownames(second), rownames(first))
print(length(overlapping_genes))
compare_first_df <- first[rownames(first) %in% overlapping_genes,]
compare_second_df <- second[rownames(second) %in% overlapping_genes,]
compare_second_df <- compare_second_df[rownames(compare_second_df) %in% rownames(compare_first_df),]
compare_first_df <- compare_first_df[order(rownames(compare_first_df)),]
compare_second_df <- compare_second_df[order(rownames(compare_second_df)),]

comparing_first_vs_second_df <- data.frame(gene_name = rownames(compare_first_df), first_fc = compare_first_df$log2FoldChange,
                                           second_fc = compare_second_df$log2FoldChange)
# Calculate correlation
correlation_val <- cor.test(comparing_first_vs_second_df$first_fc, comparing_first_vs_second_df$second_fc)
print(correlation_val$estimate)
print(correlation_val$p.value)

# Plot correlation
ggplot(comparing_first_vs_second_df, aes(x=first_fc, y=second_fc)) + 
  geom_point()+
  geom_smooth(method=lm) + xlab("Placebo FC") + ylab("Vaccinated FC")

# 2. Check placebo DEGs
placebo_degs <- rownames(first)
print(length(placebo_degs))
compare_first_df <- unfiltered_first[rownames(unfiltered_first) %in% placebo_degs,]
compare_second_df <- unfiltered_second[rownames(unfiltered_second) %in% placebo_degs,]
compare_first_df <- compare_first_df[rownames(compare_first_df) %in% rownames(compare_second_df),]
compare_second_df <- compare_second_df[rownames(compare_second_df) %in% rownames(compare_first_df),]
compare_first_df <- compare_first_df[order(rownames(compare_first_df)),]
compare_second_df <- compare_second_df[order(rownames(compare_second_df)),]

comparing_first_vs_second_df <- data.frame(gene_name = rownames(compare_first_df), first_fc = compare_first_df$log2FoldChange,
                                           second_fc = compare_second_df$log2FoldChange)
# Calculate correlation
correlation_val <- cor.test(comparing_first_vs_second_df$first_fc, comparing_first_vs_second_df$second_fc)
print(correlation_val$estimate)
print(correlation_val$p.value)

# Plot correlation
ggplot(comparing_first_vs_second_df, aes(x=first_fc, y=second_fc)) + 
  geom_point()+
  geom_smooth(method=lm) + xlab("Placebo FC") + ylab("Vaccinated FC")

# 3. Check vaccinated DEGs
vaccinated_degs <- rownames(second)
print(length(vaccinated_degs))
compare_first_df <- unfiltered_first[rownames(unfiltered_first) %in% vaccinated_degs,]
compare_second_df <- unfiltered_second[rownames(unfiltered_second) %in% vaccinated_degs,]
compare_first_df <- compare_first_df[rownames(compare_first_df) %in% rownames(compare_second_df),]
compare_second_df <- compare_second_df[rownames(compare_second_df) %in% rownames(compare_first_df),]
compare_first_df <- compare_first_df[order(rownames(compare_first_df)),]
compare_second_df <- compare_second_df[order(rownames(compare_second_df)),]

comparing_first_vs_second_df <- data.frame(gene_name = rownames(compare_first_df), first_fc = compare_first_df$log2FoldChange,
                                           second_fc = compare_second_df$log2FoldChange)
# Calculate correlation
correlation_val <- cor.test(comparing_first_vs_second_df$first_fc, comparing_first_vs_second_df$second_fc)
print(correlation_val$estimate)
print(correlation_val$p.value)

# Plot correlation
ggplot(comparing_first_vs_second_df, aes(x=first_fc, y=second_fc)) + 
  geom_point()+
  geom_smooth(method=lm) + xlab("Placebo FC") + ylab("Vaccinated FC")


# 4. Check SC DEGs
# TODO: Check correlation with 1) log and 2) scaling
sc_first <- unfiltered_first[rownames(unfiltered_first) %in% sc_pseudobulk_deg_table$Gene_Name,]
sc_second <- unfiltered_second[rownames(unfiltered_second) %in% sc_pseudobulk_deg_table$Gene_Name,]
overlapping_genes <- intersect(rownames(sc_second), rownames(sc_first))
print(length(overlapping_genes))
sc_second <- sc_second[rownames(sc_second) %in% rownames(sc_first),]
sc_first <- sc_first[order(rownames(sc_first)),]
sc_second <- sc_second[order(rownames(sc_second)),]

comparing_first_vs_second_df <- data.frame(gene_name = rownames(sc_first), first_fc = sc_first$log2FoldChange,
                                           second_fc = sc_second$log2FoldChange)
# Calculate correlation
correlation_val <- cor.test(comparing_first_vs_second_df$first_fc, comparing_first_vs_second_df$second_fc)
print(correlation_val$estimate)
print(correlation_val$p.value)

# Plot correlation
ggplot(comparing_first_vs_second_df, aes(x=first_fc, y=second_fc)) + 
  geom_point()+
  geom_smooth(method=lm) + xlab("Placebo FC") + ylab("Vaccinated FC")

# 5. Check all genes
overlapping_genes <- intersect(rownames(unfiltered_second), rownames(unfiltered_first))
print(length(overlapping_genes))
compare_first_df <- unfiltered_first[rownames(unfiltered_first) %in% overlapping_genes,]
compare_second_df <- unfiltered_second[rownames(unfiltered_second) %in% overlapping_genes,]
compare_second_df <- compare_second_df[rownames(compare_second_df) %in% rownames(compare_first_df),]
compare_first_df <- compare_first_df[order(rownames(compare_first_df)),]
compare_second_df <- compare_second_df[order(rownames(compare_second_df)),]

comparing_first_vs_second_df <- data.frame(gene_name = rownames(compare_first_df), first_fc = compare_first_df$log2FoldChange,
                                           second_fc = compare_second_df$log2FoldChange)
# Calculate correlation
correlation_val <- cor.test(comparing_first_vs_second_df$first_fc, comparing_first_vs_second_df$second_fc)
print(correlation_val$estimate)
print(correlation_val$p.value)

# Plot correlation
ggplot(comparing_first_vs_second_df, aes(x=first_fc, y=second_fc)) + 
  geom_point()+
  geom_smooth(method=lm) + xlab("Placebo FC") + ylab("Vaccinated FC")


# Run correlation for FC in all genes
unfiltered_hvl_placebo_period_2_D_minus_2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", hvl_placebo_counts, hvl_placebo_metadata,
                                                                                            "2_D_minus_2", "2_D_minus_1", paste0(bulk_rna_results_dir, "hvl_bulk_placebo_period_2_D_minus_1_vs_D_minus_2/"), "high", alpha = 0.999999)
unfiltered_hvl_placebo_period_2_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", hvl_placebo_counts, hvl_placebo_metadata,
                                                                                     "2_D2", "2_D_minus_1", paste0(bulk_rna_results_dir, "hvl_bulk_placebo_period_2_D2_vs_D_minus_1/"), "high", alpha = 0.999999)
unfiltered_hvl_placebo_period_2_D5_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", hvl_placebo_counts, hvl_placebo_metadata,
                                                                                     "2_D5", "2_D_minus_1", paste0(bulk_rna_results_dir, "hvl_bulk_placebo_period_2_D2_vs_D_minus_1/"), "high", alpha = 0.999999)
unfiltered_hvl_placebo_period_2_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", hvl_placebo_counts, hvl_placebo_metadata,
                                                                                     "2_D8", "2_D_minus_1", paste0(bulk_rna_results_dir, "hvl_bulk_placebo_period_2_D2_vs_D_minus_1/"), "high", alpha = 0.999999)
unfiltered_hvl_placebo_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", hvl_placebo_counts, hvl_placebo_metadata,
                                                                                     "2_D28", "2_D_minus_1", paste0(bulk_rna_results_dir, "hvl_bulk_placebo_period_2_D2_vs_D_minus_1/"), "high", alpha = 0.999999)
unfiltered_hvl_vaccinated_period_2_D_minus_2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", hvl_vaccinated_counts, hvl_vaccinated_metadata,
                                                                                                         "2_D_minus_2", "2_D_minus_1", paste0(bulk_rna_results_dir, "hvl_bulk_vaccinated_period_2_D_minus_1_vs_D_minus_2/"), "high", alpha = 0.999999)
unfiltered_hvl_vaccinated_period_2_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", hvl_vaccinated_counts, hvl_vaccinated_metadata,
                                                                                       "2_D2", "2_D_minus_1", paste0(bulk_rna_results_dir, "hvl_bulk_vaccinated_period_2_D2_vs_D_minus_1/"), "high", alpha = 0.999999)
unfiltered_hvl_vaccinated_period_2_D5_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", hvl_vaccinated_counts, hvl_vaccinated_metadata,
                                                                                       "2_D5", "2_D_minus_1", paste0(bulk_rna_results_dir, "hvl_bulk_vaccinated_period_2_D2_vs_D_minus_1/"), "high", alpha = 0.999999)
unfiltered_hvl_vaccinated_period_2_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", hvl_vaccinated_counts, hvl_vaccinated_metadata,
                                                                                       "2_D8", "2_D_minus_1", paste0(bulk_rna_results_dir, "hvl_bulk_vaccinated_period_2_D2_vs_D_minus_1/"), "high", alpha = 0.999999)
unfiltered_hvl_vaccinated_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", hvl_vaccinated_counts, hvl_vaccinated_metadata,
                                                                                       "2_D28", "2_D_minus_1", paste0(bulk_rna_results_dir, "hvl_bulk_vaccinated_period_2_D2_vs_D_minus_1/"), "high", alpha = 0.999999)

unfiltered_lvl_placebo_period_2_D_minus_2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", lvl_placebo_counts, lvl_placebo_metadata,
                                                                                                       "2_D_minus_2", "2_D_minus_1", paste0(bulk_rna_results_dir, "lvl_bulk_placebo_period_2_D_minus_1_vs_D_minus_2/"), "low", alpha = 0.999999)
unfiltered_lvl_placebo_period_2_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", lvl_placebo_counts, lvl_placebo_metadata,
                                                                                                "2_D2", "2_D_minus_1", paste0(bulk_rna_results_dir, "lvl_bulk_placebo_period_2_D2_vs_D_minus_1/"), "low", alpha = 0.999999)
unfiltered_lvl_placebo_period_2_D5_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", lvl_placebo_counts, lvl_placebo_metadata,
                                                                                                "2_D5", "2_D_minus_1", paste0(bulk_rna_results_dir, "lvl_bulk_placebo_period_2_D2_vs_D_minus_1/"), "low", alpha = 0.999999)
unfiltered_lvl_placebo_period_2_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", lvl_placebo_counts, lvl_placebo_metadata,
                                                                                                "2_D8", "2_D_minus_1", paste0(bulk_rna_results_dir, "lvl_bulk_placebo_period_2_D2_vs_D_minus_1/"), "low", alpha = 0.999999)
unfiltered_lvl_placebo_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", lvl_placebo_counts, lvl_placebo_metadata,
                                                                                                 "2_D28", "2_D_minus_1", paste0(bulk_rna_results_dir, "lvl_bulk_placebo_period_2_D2_vs_D_minus_1/"), "low", alpha = 0.999999)
unfiltered_lvl_vaccinated_period_2_D_minus_2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", lvl_vaccinated_counts, lvl_vaccinated_metadata,
                                                                                                         "2_D_minus_2", "2_D_minus_1", paste0(bulk_rna_results_dir, "lvl_bulk_vaccinated_period_2_D_minus_1_vs_D_minus_2/"), "low", alpha = 0.999999)
unfiltered_lvl_vaccinated_period_2_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", lvl_vaccinated_counts, lvl_vaccinated_metadata,
                                                                                                  "2_D2", "2_D_minus_1", paste0(bulk_rna_results_dir, "lvl_bulk_vaccinated_period_2_D2_vs_D_minus_1/"), "low", alpha = 0.999999)
unfiltered_lvl_vaccinated_period_2_D5_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", lvl_vaccinated_counts, lvl_vaccinated_metadata,
                                                                                                  "2_D5", "2_D_minus_1", paste0(bulk_rna_results_dir, "lvl_bulk_vaccinated_period_2_D2_vs_D_minus_1/"), "low", alpha = 0.999999)
unfiltered_lvl_vaccinated_period_2_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", lvl_vaccinated_counts, lvl_vaccinated_metadata,
                                                                                                  "2_D8", "2_D_minus_1", paste0(bulk_rna_results_dir, "lvl_bulk_vaccinated_period_2_D2_vs_D_minus_1/"), "low", alpha = 0.999999)
unfiltered_lvl_vaccinated_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", lvl_vaccinated_counts, lvl_vaccinated_metadata,
                                                                                                   "2_D28", "2_D_minus_1", paste0(bulk_rna_results_dir, "lvl_bulk_vaccinated_period_2_D2_vs_D_minus_1/"), "low", alpha = 0.999999)




# HVL PLACEBO VS HVL VACCINATED
# 150 DEGs for D28 (placebo) vs D28 (vaccinated)
hvl_placebo_metadata_for_comparison <- hvl_placebo_metadata[hvl_placebo_metadata$time_point == "2_D28",]
hvl_placebo_metadata_for_comparison$status <- "placebo"
hvl_vaccinated_metadata_for_comparison <- hvl_vaccinated_metadata[hvl_vaccinated_metadata$time_point == "2_D28",]
hvl_vaccinated_metadata_for_comparison$status <- "vaccinated"
hvl_placebo_and_vaccinated_metadata <- rbind(hvl_placebo_metadata_for_comparison, hvl_vaccinated_metadata_for_comparison)
hvl_placebo_and_vaccinated_counts <- cbind(hvl_placebo_counts, hvl_vaccinated_counts)

# Select the relevant time point from our metadata
hvl_placebo_and_vaccinated_metadata_subset <- hvl_placebo_and_vaccinated_metadata
# Select subset of counts associated with subjects
hvl_placebo_and_vaccinated_counts_subset <- hvl_placebo_and_vaccinated_counts[rownames(hvl_placebo_and_vaccinated_metadata_subset)]
# Run DESeq2
current_analysis <- DESeqDataSetFromMatrix(countData = hvl_placebo_and_vaccinated_counts_subset, 
                                           colData = hvl_placebo_and_vaccinated_metadata_subset, 
                                           design = ~ sex + age + status)
current_analysis <- DESeq(current_analysis)
current_analysis_results <- results(current_analysis, contrast = c("status", "placebo", "vaccinated"), alpha = 0.05)
current_analysis_results <- current_analysis_results[order(current_analysis_results$padj),]
current_analysis_results <- subset(current_analysis_results, padj < 0.05)

# HIGH VS LOW VIRAL LOAD
# OVERALL, everything significant for D-2 / D-1 / D2 / D28 are very low base mean and kind of nonsense
# D5 / D8 are higher base mean and actually have real DEG
# This makes sense
# 2 D minus 2
# 144/143/142/138/151
# Note: If I don't control for age and sex, I find nothing.
placebo_period_2_D_minus_2_high_vs_low_results <- run_deseq_bulk_analysis_viral_load("both", both_full_time_series_placebo_counts, both_full_time_series_placebo_metadata,
                                                                              "2_D_minus_2", "high", "low", paste0(bulk_rna_results_dir, "placebo_period_2_D_minus_2_high_vs_low/"), "2_D_minus_2")
# 2 D minus 1
# 167/157/155/151/178
# Note: If I don't control for age and sex, I find nothing.
placebo_period_2_D_minus_1_high_vs_low_results <- run_deseq_bulk_analysis_viral_load("both", both_full_time_series_placebo_counts, both_full_time_series_placebo_metadata,
                                                                              "2_D_minus_1", "high", "low", paste0(bulk_rna_results_dir, "placebo_period_2_D_minus_1_high_vs_low/"), "2_D_minus_1")
# 2 D2
# 139/138/137/134/139
placebo_period_2_D2_high_vs_low_results <- run_deseq_bulk_analysis_viral_load("both", both_full_time_series_placebo_counts, both_full_time_series_placebo_metadata,
                                                                                     "2_D2", "high", "low", paste0(bulk_rna_results_dir, "placebo_period_2_D2_high_vs_low/"), "2_D2")
# 2 D5
# Note: If I don't control for age and sex, I still find a lot of DEGs.
# https://hb.flatironinstitute.org/module/overview/?body_tag=9a43682013acc9404b08417d64f9e60f2436bfba all upregulated
# 883/269/155/88/1653
placebo_period_2_D5_high_vs_low_results <- run_deseq_bulk_analysis_viral_load("both", both_full_time_series_placebo_counts, both_full_time_series_placebo_metadata,
                                                                              "2_D5", "high", "low", paste0(bulk_rna_results_dir, "placebo_period_2_D5_high_vs_low/"), "2_D5")
# 2 D8
# https://hb.flatironinstitute.org/module/overview/?body_tag=e7397a8502c94a074b03dedb68c473ef8912a8e9 all upregulated
# 298/207/196/187/365
placebo_period_2_D8_high_vs_low_results <- run_deseq_bulk_analysis_viral_load("both", both_full_time_series_placebo_counts, both_full_time_series_placebo_metadata,
                                                                              "2_D8", "high", "low", paste0(bulk_rna_results_dir, "placebo_period_2_D8_high_vs_low/"), "2_D8")
# 2 D28
# 193/180/172/170/209
placebo_period_2_D28_high_vs_low_results <- run_deseq_bulk_analysis_viral_load("both", both_full_time_series_placebo_counts, both_full_time_series_placebo_metadata,
                                                                              "2_D28", "high", "low", paste0(bulk_rna_results_dir, "placebo_period_2_D8_high_vs_low/"), "2_D28")

# VACCINATED
#### PERIOD 1 HIGH T CELL RESPONSE ####
# 6234/1686/678/163/7799 DEGs
high_vaccinated_period_1_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", high_vaccinated_counts, high_vaccinated_metadata,
                                                                                     "1_D2", "1_D_minus_1", paste0(bulk_rna_results_dir, "high_vaccinated_period_1_D2_vs_D_minus_1/"), "high")
raw_high_vaccinated_period_1_D2_vs_D_minus_1_results <- high_vaccinated_period_1_D2_vs_D_minus_1_results[[5]]
# 59/15/2/0/96 DEGs
high_vaccinated_period_1_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", high_vaccinated_counts, high_vaccinated_metadata,
                                                                                     "1_D8", "1_D_minus_1", paste0(bulk_rna_results_dir, "high_vaccinated_period_1_D8_vs_D_minus_1/"), "high")
raw_high_vaccinated_period_1_D8_vs_D_minus_1_results <- high_vaccinated_period_1_D8_vs_D_minus_1_results[[5]]
# 0/0/0/0/0 DEGs
high_vaccinated_period_1_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", high_vaccinated_counts, high_vaccinated_metadata,
                                                                                      "1_D28", "1_D_minus_1", paste0(bulk_rna_results_dir, "high_vaccinated_period_1_D28_vs_D_minus_1/"), "high")
raw_high_vaccinated_period_1_D28_vs_D_minus_1_results <- high_vaccinated_period_1_D28_vs_D_minus_1_results[[5]]
#### PERIOD 1 LOW VIRAL LOAD ####
# 2280/312/129/18/3854 DEGs
low_vaccinated_period_1_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", low_vaccinated_counts, low_vaccinated_metadata,
                                                                                    "1_D2", "1_D_minus_1", paste0(bulk_rna_results_dir, "low_vaccinated_period_1_D2_vs_D_minus_1/"), "low")
raw_low_vaccinated_period_1_D2_vs_D_minus_1_results <- low_vaccinated_period_1_D2_vs_D_minus_1_results[[5]]
# 70/0/0/0/521 DEGs
low_vaccinated_period_1_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", low_vaccinated_counts, low_vaccinated_metadata,
                                                                                    "1_D8", "1_D_minus_1", paste0(bulk_rna_results_dir, "low_vaccinated_period_1_D8_vs_D_minus_1/"), "low")
raw_low_vaccinated_period_1_D8_vs_D_minus_1_results <- low_vaccinated_period_1_D8_vs_D_minus_1_results[[5]]
# 2/0/0/0/3 DEGs
low_vaccinated_period_1_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", low_vaccinated_counts, low_vaccinated_metadata,
                                                                                     "1_D28", "1_D_minus_1", paste0(bulk_rna_results_dir, "low_vaccinated_period_1_D28_vs_D_minus_1/"), "low")
raw_low_vaccinated_period_1_D28_vs_D_minus_1_results <- low_vaccinated_period_1_D28_vs_D_minus_1_results[[5]]

#### PERIOD 2 HIGH VIRAL LOAD ####
# 2 D minus 2 vs 2 D minus 1 - should be virtually zero unless some weird stuff happened between blood draws
# 294/0/0/0/1225 DEGs found - seems weird that so many are found
high_vaccinated_period_2_D_minus_1_vs_D_minus_2_results <- run_deseq_bulk_analysis_time_series("vaccinated", high_vaccinated_counts, high_vaccinated_metadata,
                                                                                            "2_D_minus_1", "2_D_minus_2", paste0(bulk_rna_results_dir, "high_vaccinated_period_2_D_minus_1_vs_D_minus_2/"), "high")
raw_high_vaccinated_period_2_D_minus_1_vs_D_minus_2_results <- high_vaccinated_period_2_D_minus_1_vs_D_minus_2_results[[5]]
# 0.585 is good, and -0.585 is good
#collected_fmd_high_vaccinated_period_2_D_minus_1_vs_D_minus_2 <- run_fmd_on_flu_data(raw_high_vaccinated_period_2_D_minus_1_vs_D_minus_2_results)

# 2 D2 vs 2 D minus 1 - 0/0/0/0/2 DEGs
high_vaccinated_period_2_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", high_vaccinated_counts, high_vaccinated_metadata,
                                                                                     "2_D2", "2_D_minus_1", paste0(bulk_rna_results_dir, "high_vaccinated_period_2_D2_vs_D_minus_1/"), "high")
raw_high_vaccinated_period_2_D2_vs_D_minus_1_results <- high_vaccinated_period_2_D2_vs_D_minus_1_results[[5]]
# 0.585 is good, and -0.585 is good
#collected_fmd_high_vaccinated_period_2_D2_vs_D_minus_1 <- run_fmd_on_flu_data(raw_high_vaccinated_period_2_D2_vs_D_minus_1_results)

# 2 D5 vs 2 D minus 1 - 704/20/0/0/1588 DEGs
high_vaccinated_period_2_D5_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", high_vaccinated_counts, high_vaccinated_metadata,
                                                                                     "2_D5", "2_D_minus_1", paste0(bulk_rna_results_dir, "high_vaccinated_period_2_D5_vs_D_minus_1/"), "high")
raw_high_vaccinated_period_2_D5_vs_D_minus_1_results <- high_vaccinated_period_2_D5_vs_D_minus_1_results[[5]]
# 0.585 or 1 is good, and -0.585 or -1 is good 
#collected_fmd_high_vaccinated_period_2_D5_vs_D_minus_1 <- run_fmd_on_flu_data(raw_high_vaccinated_period_2_D5_vs_D_minus_1_results, c(0.585, 1, 2))

# 2 D8 vs 2 D minus 1 - 1359/85/1/0/2441 DEGs
high_vaccinated_period_2_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", high_vaccinated_counts, high_vaccinated_metadata,
                                                                                     "2_D8", "2_D_minus_1", paste0(bulk_rna_results_dir, "high_vaccinated_period_2_D8_vs_D_minus_1/"), "high")
raw_high_vaccinated_period_2_D8_vs_D_minus_1_results <- high_vaccinated_period_2_D8_vs_D_minus_1_results[[5]]
# ???
# collected_fmd_high_vaccinated_period_2_D8_vs_D_minus_1 <- run_fmd_on_flu_data(raw_high_vaccinated_period_2_D8_vs_D_minus_1_results, c(0.585, 1, 2))
# 2 D28 vs 2 D minus 1 - 6/4/4/3/21 DEGs
high_vaccinated_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", high_vaccinated_counts, high_vaccinated_metadata,
                                                                                      "2_D28", "2_D_minus_1", paste0(bulk_rna_results_dir, "high_vaccinated_period_2_D28_vs_D_minus_1/"), "high")
raw_high_vaccinated_period_2_D28_vs_D_minus_1_results <- high_vaccinated_period_2_D28_vs_D_minus_1_results[[5]]
# 0.3 is good, and -0.3 is good
#collected_fmd_high_vaccinated_period_2_D28_vs_D_minus_1 <- run_fmd_on_flu_data(raw_high_vaccinated_period_2_D28_vs_D_minus_1_results)

#### PERIOD 2 LOW VIRAL LOAD ####
# 2 D minus 2 vs 2 D minus 1 - should be virtually zero unless some weird stuff happened between blood draws
# 1/0/0/0/9 DEGs found - seems weird that so many are found
low_vaccinated_period_2_D_minus_1_vs_D_minus_2_results <- run_deseq_bulk_analysis_time_series("vaccinated", low_vaccinated_counts, low_vaccinated_metadata,
                                                                                               "2_D_minus_1", "2_D_minus_2", paste0(bulk_rna_results_dir, "low_vaccinated_period_2_D_minus_1_vs_D_minus_2/"), "low")
raw_low_vaccinated_period_2_D_minus_1_vs_D_minus_2_results <- low_vaccinated_period_2_D_minus_1_vs_D_minus_2_results[[5]]
# 0.585 is good, and -0.585 is good
#collected_fmd_low_vaccinated_period_2_D_minus_1_vs_D_minus_2 <- run_fmd_on_flu_data(raw_low_vaccinated_period_2_D_minus_1_vs_D_minus_2_results)

# 2 D2 vs 2 D minus 1 - 2/0/0/0/2 DEGs
low_vaccinated_period_2_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", low_vaccinated_counts, low_vaccinated_metadata,
                                                                                        "2_D2", "2_D_minus_1", paste0(bulk_rna_results_dir, "low_vaccinated_period_2_D2_vs_D_minus_1/"), "low")
raw_low_vaccinated_period_2_D2_vs_D_minus_1_results <- low_vaccinated_period_2_D2_vs_D_minus_1_results[[5]]
# 0.585 is good, and -0.585 is good
#collected_fmd_low_vaccinated_period_2_D2_vs_D_minus_1 <- run_fmd_on_flu_data(raw_low_vaccinated_period_2_D2_vs_D_minus_1_results)

# 2 D5 vs 2 D minus 1 - 101/20/0/0/366 DEGs
low_vaccinated_period_2_D5_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", low_vaccinated_counts, low_vaccinated_metadata,
                                                                                        "2_D5", "2_D_minus_1", paste0(bulk_rna_results_dir, "low_vaccinated_period_2_D5_vs_D_minus_1/"), "low")
raw_low_vaccinated_period_2_D5_vs_D_minus_1_results <- low_vaccinated_period_2_D5_vs_D_minus_1_results[[5]]
# 0.585 or 1 is good, and -0.585 or -1 is good 
#collected_fmd_low_vaccinated_period_2_D5_vs_D_minus_1 <- run_fmd_on_flu_data(raw_low_vaccinated_period_2_D5_vs_D_minus_1_results, c(0.585, 1, 2))

# 2 D8 vs 2 D minus 1 - 4/0/0/0/15 DEGs
low_vaccinated_period_2_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", low_vaccinated_counts, low_vaccinated_metadata,
                                                                                        "2_D8", "2_D_minus_1", paste0(bulk_rna_results_dir, "low_vaccinated_period_2_D8_vs_D_minus_1/"), "low")
#raw_low_vaccinated_period_2_D8_vs_D_minus_1_results <- low_vaccinated_period_2_D8_vs_D_minus_1_results[[5]]
# ???
# collected_fmd_low_vaccinated_period_2_D8_vs_D_minus_1 <- run_fmd_on_flu_data(raw_low_vaccinated_period_2_D8_vs_D_minus_1_results, c(0.585, 1, 2))
# 2 D28 vs 2 D minus 1 - 106/0/0/0/800 DEGs
low_vaccinated_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", low_vaccinated_counts, low_vaccinated_metadata,
                                                                                         "2_D28", "2_D_minus_1", paste0(bulk_rna_results_dir, "low_vaccinated_period_2_D28_vs_D_minus_1/"), "low")
raw_low_vaccinated_period_2_D28_vs_D_minus_1_results <- low_vaccinated_period_2_D28_vs_D_minus_1_results[[5]]
# 0.3 is good, and -0.3 is good
#collected_fmd_low_vaccinated_period_2_D28_vs_D_minus_1 <- run_fmd_on_flu_data(raw_low_vaccinated_period_2_D28_vs_D_minus_1_results)


















# LRT test for period 2 (HVL)
# 6431 DEGs found (note that LRT doesn't have a fold change threshold)
hvl_placebo_period_2_LRT_metadata <- hvl_placebo_metadata[hvl_placebo_metadata$time_point == "2_D28" | hvl_placebo_metadata$time_point == "2_D8" | 
                                                              hvl_placebo_metadata$time_point == "2_D5" | hvl_placebo_metadata$time_point == "2_D2" |
                                                              hvl_placebo_metadata$time_point == "2_D_minus_1",]
hvl_placebo_period_2_LRT_metadata <- hvl_placebo_period_2_LRT_metadata[hvl_placebo_period_2_LRT_metadata$subject_id 
                                                                         %in% names(table(hvl_placebo_period_2_LRT_metadata$subject_id)
                                                                                    [table(hvl_placebo_period_2_LRT_metadata$subject_id) == 5]),]
hvl_placebo_period_2_LRT_counts <- hvl_placebo_counts[rownames(hvl_placebo_period_2_LRT_metadata)]
hvl_placebo_period_2_LRT_analysis <- DESeqDataSetFromMatrix(countData = hvl_placebo_period_2_LRT_counts,
                                                             colData = hvl_placebo_period_2_LRT_metadata,
                                                             design = ~ subject_id + Monocytes + Neutrophils + time_point)
hvl_placebo_period_2_LRT_analysis <- DESeq(hvl_placebo_period_2_LRT_analysis, test = "LRT", reduced = ~ subject_id + Monocytes + Neutrophils)
hvl_placebo_period_2_LRT_analysis_results <- results(hvl_placebo_period_2_LRT_analysis, alpha = 0.05)
hvl_placebo_period_2_LRT_analysis_results <- hvl_placebo_period_2_LRT_analysis_results[order(hvl_placebo_period_2_LRT_analysis_results$padj),]
hvl_placebo_period_2_LRT_analysis_results <- subset(hvl_placebo_period_2_LRT_analysis_results, padj < 0.05)
# Write top 2000 genes to file
hvl_placebo_period_2_LRT_analysis_2000_topGenes <- rownames(hvl_placebo_period_2_LRT_analysis_results)[1:2000]
if (!dir.exists(paste0(bulk_data_dir, "hvl_placebo_period_2_LRT"))) {dir.create(paste0(bulk_data_dir, "hvl_placebo_period_2_LRT")) }
write.table(hvl_placebo_period_2_LRT_analysis_2000_topGenes, file = paste0(bulk_data_dir, "hvl_placebo_period_2_LRT/LRT_top_200_genes.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
# Plot heatmap with top 20 genes
hvl_placebo_period_2_LRT_analysis_betas <- coef(hvl_placebo_period_2_LRT_analysis)
hvl_placebo_period_2_LRT_analysis_results_for_plotting <- results(hvl_placebo_period_2_LRT_analysis, alpha = 0.05)
hvl_placebo_period_2_LRT_analysis_20_topGenes <- head(order(hvl_placebo_period_2_LRT_analysis_results_for_plotting$padj),20)
hvl_placebo_period_2_LRT_analysis_20_mat <- hvl_placebo_period_2_LRT_analysis_betas[hvl_placebo_period_2_LRT_analysis_20_topGenes, -c(1:13)]
hvl_placebo_period_2_LRT_analysis_thr <- 6
hvl_placebo_period_2_LRT_analysis_20_mat[hvl_placebo_period_2_LRT_analysis_20_mat < -hvl_placebo_period_2_LRT_analysis_thr] <- -hvl_placebo_period_2_LRT_analysis_thr
hvl_placebo_period_2_LRT_analysis_20_mat[hvl_placebo_period_2_LRT_analysis_20_mat > hvl_placebo_period_2_LRT_analysis_thr] <- hvl_placebo_period_2_LRT_analysis_thr
colnames(hvl_placebo_period_2_LRT_analysis_20_mat) <- c("Day 2 vs Day -1", "Day 5 vs Day -1", "Day 8 vs Day -1", "Day 28 vs Day -1")
pheatmap(hvl_placebo_period_2_LRT_analysis_20_mat, breaks=seq(from=-hvl_placebo_period_2_LRT_analysis_thr, to=hvl_placebo_period_2_LRT_analysis_thr, length=101),
         cluster_col=FALSE, fontsize_col=14, filename = paste0(bulk_data_dir, "hvl_placebo_period_2_LRT/hvl_placebo_period_2_LRT_analysis_top_20_genes.png"))
# Plot heatmap with top 200 genes (and clean up heatmap plot by removing row names and clustering)
hvl_placebo_period_2_LRT_analysis_200_topGenes <- head(order(hvl_placebo_period_2_LRT_analysis_results_for_plotting$padj),200)
hvl_placebo_period_2_LRT_analysis_200_mat <- hvl_placebo_period_2_LRT_analysis_betas[hvl_placebo_period_2_LRT_analysis_200_topGenes, -c(1:13)]
colnames(hvl_placebo_period_2_LRT_analysis_200_mat) <- c("Day 2 vs Day -1", "Day 5 vs Day -1", "Day 8 vs Day -1", "Day 28 vs Day -1")
hvl_placebo_period_2_LRT_analysis_200_mat[hvl_placebo_period_2_LRT_analysis_200_mat < -hvl_placebo_period_2_LRT_analysis_thr] <- -hvl_placebo_period_2_LRT_analysis_thr
hvl_placebo_period_2_LRT_analysis_200_mat[hvl_placebo_period_2_LRT_analysis_200_mat > hvl_placebo_period_2_LRT_analysis_thr] <- hvl_placebo_period_2_LRT_analysis_thr
pheatmap(hvl_placebo_period_2_LRT_analysis_200_mat, breaks=seq(from=-hvl_placebo_period_2_LRT_analysis_thr, to=hvl_placebo_period_2_LRT_analysis_thr, length=101),
         cluster_col=FALSE, fontsize_col=14, show_rownames = FALSE, cluster_row = FALSE, filename = paste0(bulk_data_dir, "hvl_placebo_period_2_LRT/hvl_placebo_period_2_LRT_analysis_top_200_genes.png"))


# LRT test for period 2 (LVL)
# 0 DEGs found (note that LRT doesn't have a fold change threshold)
lvl_placebo_period_2_LRT_metadata <- lvl_placebo_metadata[lvl_placebo_metadata$time_point == "2_D28" | lvl_placebo_metadata$time_point == "2_D8" | 
                                                            lvl_placebo_metadata$time_point == "2_D5" | lvl_placebo_metadata$time_point == "2_D2" |
                                                            lvl_placebo_metadata$time_point == "2_D_minus_1",]
lvl_placebo_period_2_LRT_metadata <- lvl_placebo_period_2_LRT_metadata[lvl_placebo_period_2_LRT_metadata$subject_id 
                                                                       %in% names(table(lvl_placebo_period_2_LRT_metadata$subject_id)
                                                                                  [table(lvl_placebo_period_2_LRT_metadata$subject_id) == 5]),]
lvl_placebo_period_2_LRT_counts <- lvl_placebo_counts[rownames(lvl_placebo_period_2_LRT_metadata)]
lvl_placebo_period_2_LRT_analysis <- DESeqDataSetFromMatrix(countData = lvl_placebo_period_2_LRT_counts,
                                                            colData = lvl_placebo_period_2_LRT_metadata,
                                                            design = ~ subject_id + time_point)
lvl_placebo_period_2_LRT_analysis <- DESeq(lvl_placebo_period_2_LRT_analysis, test = "LRT", reduced = ~ subject_id)
lvl_placebo_period_2_LRT_analysis_results <- results(lvl_placebo_period_2_LRT_analysis, alpha = 0.05)
lvl_placebo_period_2_LRT_analysis_results <- lvl_placebo_period_2_LRT_analysis_results[order(lvl_placebo_period_2_LRT_analysis_results$padj),]
lvl_placebo_period_2_LRT_analysis_results <- subset(lvl_placebo_period_2_LRT_analysis_results, padj < 0.05)




