# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# High placebo
hvl_placebo_second_period_metadata <- hvl_placebo_metadata[hvl_placebo_metadata$period == "2",]
hvl_placebo_second_period_metadata <- hvl_placebo_second_period_metadata[hvl_placebo_second_period_metadata$time_point != "2_D_minus_2",]
hvl_placebo_second_period_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, hvl_placebo_second_period_metadata, contrast = c("2_D_minus_1", "2_D2", "2_D5", "2_D8", "2_D28"))

# All matched placebo (period 2)
matching_all_placebo_second_period_metadata <- placebo_metadata[placebo_metadata$period == "2",]
matching_all_placebo_second_period_metadata <- matching_all_placebo_second_period_metadata[matching_all_placebo_second_period_metadata$time_point != "2_D_minus_2",]
matching_all_placebo_second_period_metadata <- matching_all_placebo_second_period_metadata[matching_all_placebo_second_period_metadata$time_point != "2_D2",]
matching_all_placebo_second_period_metadata <- matching_all_placebo_second_period_metadata[matching_all_placebo_second_period_metadata$time_point != "2_D5",]
matching_all_placebo_second_period_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, matching_all_placebo_second_period_metadata, contrast = c("2_D_minus_1", "2_D8", "2_D28"))

# All matched vaccine (period 2)
matching_all_vaccinated_second_period_metadata <- vaccinated_metadata[vaccinated_metadata$period == "2",]
matching_all_vaccinated_second_period_metadata <- matching_all_vaccinated_second_period_metadata[matching_all_vaccinated_second_period_metadata$time_point != "2_D_minus_2",]
matching_all_vaccinated_second_period_metadata <- matching_all_vaccinated_second_period_metadata[matching_all_vaccinated_second_period_metadata$time_point != "2_D2",]
matching_all_vaccinated_second_period_metadata <- matching_all_vaccinated_second_period_metadata[matching_all_vaccinated_second_period_metadata$time_point != "2_D5",]
matching_all_vaccinated_second_period_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, matching_all_vaccinated_second_period_metadata, contrast = c("2_D_minus_1", "2_D8", "2_D28"))

# All matched samples (period 2)
matching_all_samples_second_period_metadata <- rbind(matching_all_placebo_second_period_metadata, matching_all_vaccinated_second_period_metadata)
matching_all_samples_second_period_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, matching_all_samples_second_period_metadata, contrast = c("2_D_minus_1", "2_D8", "2_D28"))

# All matched vaccine (period 1)
matching_all_vaccinated_first_period_metadata <- vaccinated_metadata[vaccinated_metadata$period == "1",]
matching_all_vaccinated_first_period_metadata <- matching_all_vaccinated_first_period_metadata[matching_all_vaccinated_first_period_metadata$time_point != "1_D2",]
matching_all_vaccinated_first_period_metadata <- matching_all_vaccinated_first_period_metadata[matching_all_vaccinated_first_period_metadata$time_point != "1_D28",]
matching_all_vaccinated_first_period_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, matching_all_vaccinated_first_period_metadata, contrast = c("1_D_minus_1", "1_D8"))

# All matched placebo (binary)
matching_all_placebo_2_D28_binary_metadata <- placebo_metadata[placebo_metadata$period == "2",]
matching_all_placebo_2_D28_binary_metadata <- matching_all_placebo_2_D28_binary_metadata[matching_all_placebo_2_D28_binary_metadata$time_point != "2_D_minus_2",]
matching_all_placebo_2_D28_binary_metadata <- matching_all_placebo_2_D28_binary_metadata[matching_all_placebo_2_D28_binary_metadata$time_point != "2_D2",]
matching_all_placebo_2_D28_binary_metadata <- matching_all_placebo_2_D28_binary_metadata[matching_all_placebo_2_D28_binary_metadata$time_point != "2_D5",]
matching_all_placebo_2_D28_binary_metadata <- matching_all_placebo_2_D28_binary_metadata[matching_all_placebo_2_D28_binary_metadata$time_point != "2_D8",]
matching_all_placebo_2_D28_binary_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, matching_all_placebo_2_D28_binary_metadata, contrast = c("2_D_minus_1", "2_D28"))

matching_all_placebo_2_D8_binary_metadata <- placebo_metadata[placebo_metadata$period == "2",]
matching_all_placebo_2_D8_binary_metadata <- matching_all_placebo_2_D8_binary_metadata[matching_all_placebo_2_D8_binary_metadata$time_point != "2_D_minus_2",]
matching_all_placebo_2_D8_binary_metadata <- matching_all_placebo_2_D8_binary_metadata[matching_all_placebo_2_D8_binary_metadata$time_point != "2_D2",]
matching_all_placebo_2_D8_binary_metadata <- matching_all_placebo_2_D8_binary_metadata[matching_all_placebo_2_D8_binary_metadata$time_point != "2_D5",]
matching_all_placebo_2_D8_binary_metadata <- matching_all_placebo_2_D8_binary_metadata[matching_all_placebo_2_D8_binary_metadata$time_point != "2_D28",]
matching_all_placebo_2_D8_binary_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, matching_all_placebo_2_D8_binary_metadata, contrast = c("2_D_minus_1", "2_D8"))

matching_all_placebo_2_D_minus_2_binary_metadata <- placebo_metadata[placebo_metadata$period == "2",]
matching_all_placebo_2_D_minus_2_binary_metadata <- matching_all_placebo_2_D_minus_2_binary_metadata[matching_all_placebo_2_D_minus_2_binary_metadata$time_point != "2_D8",]
matching_all_placebo_2_D_minus_2_binary_metadata <- matching_all_placebo_2_D_minus_2_binary_metadata[matching_all_placebo_2_D_minus_2_binary_metadata$time_point != "2_D2",]
matching_all_placebo_2_D_minus_2_binary_metadata <- matching_all_placebo_2_D_minus_2_binary_metadata[matching_all_placebo_2_D_minus_2_binary_metadata$time_point != "2_D5",]
matching_all_placebo_2_D_minus_2_binary_metadata <- matching_all_placebo_2_D_minus_2_binary_metadata[matching_all_placebo_2_D_minus_2_binary_metadata$time_point != "2_D28",]
matching_all_placebo_2_D_minus_2_binary_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, matching_all_placebo_2_D_minus_2_binary_metadata, contrast = c("2_D_minus_1", "2_D_minus_2"))

# All matched vaccinated (binary)
matching_all_vaccinated_2_D28_binary_metadata <- vaccinated_metadata[vaccinated_metadata$period == "2",]
matching_all_vaccinated_2_D28_binary_metadata <- matching_all_vaccinated_2_D28_binary_metadata[matching_all_vaccinated_2_D28_binary_metadata$time_point != "2_D_minus_2",]
matching_all_vaccinated_2_D28_binary_metadata <- matching_all_vaccinated_2_D28_binary_metadata[matching_all_vaccinated_2_D28_binary_metadata$time_point != "2_D2",]
matching_all_vaccinated_2_D28_binary_metadata <- matching_all_vaccinated_2_D28_binary_metadata[matching_all_vaccinated_2_D28_binary_metadata$time_point != "2_D5",]
matching_all_vaccinated_2_D28_binary_metadata <- matching_all_vaccinated_2_D28_binary_metadata[matching_all_vaccinated_2_D28_binary_metadata$time_point != "2_D8",]
matching_all_vaccinated_2_D28_binary_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, matching_all_vaccinated_2_D28_binary_metadata, contrast = c("2_D_minus_1", "2_D28"))

matching_all_vaccinated_2_D8_binary_metadata <- vaccinated_metadata[vaccinated_metadata$period == "2",]
matching_all_vaccinated_2_D8_binary_metadata <- matching_all_vaccinated_2_D8_binary_metadata[matching_all_vaccinated_2_D8_binary_metadata$time_point != "2_D_minus_2",]
matching_all_vaccinated_2_D8_binary_metadata <- matching_all_vaccinated_2_D8_binary_metadata[matching_all_vaccinated_2_D8_binary_metadata$time_point != "2_D2",]
matching_all_vaccinated_2_D8_binary_metadata <- matching_all_vaccinated_2_D8_binary_metadata[matching_all_vaccinated_2_D8_binary_metadata$time_point != "2_D5",]
matching_all_vaccinated_2_D8_binary_metadata <- matching_all_vaccinated_2_D8_binary_metadata[matching_all_vaccinated_2_D8_binary_metadata$time_point != "2_D28",]
matching_all_vaccinated_2_D8_binary_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, matching_all_vaccinated_2_D8_binary_metadata, contrast = c("2_D_minus_1", "2_D8"))

# All samples (binary)
matching_all_samples_2_D28_binary_metadata <- matching_all_samples_second_period_metadata
matching_all_samples_2_D28_binary_metadata <- matching_all_samples_2_D28_binary_metadata[matching_all_samples_2_D28_binary_metadata$time_point != "2_D_minus_2",]
matching_all_samples_2_D28_binary_metadata <- matching_all_samples_2_D28_binary_metadata[matching_all_samples_2_D28_binary_metadata$time_point != "2_D2",]
matching_all_samples_2_D28_binary_metadata <- matching_all_samples_2_D28_binary_metadata[matching_all_samples_2_D28_binary_metadata$time_point != "2_D5",]
matching_all_samples_2_D28_binary_metadata <- matching_all_samples_2_D28_binary_metadata[matching_all_samples_2_D28_binary_metadata$time_point != "2_D8",]
matching_all_samples_2_D28_binary_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, matching_all_samples_2_D28_binary_metadata, contrast = c("2_D_minus_1", "2_D28"))

matching_all_samples_2_D8_binary_metadata <- matching_all_samples_second_period_metadata
matching_all_samples_2_D8_binary_metadata <- matching_all_samples_2_D8_binary_metadata[matching_all_samples_2_D8_binary_metadata$time_point != "2_D_minus_2",]
matching_all_samples_2_D8_binary_metadata <- matching_all_samples_2_D8_binary_metadata[matching_all_samples_2_D8_binary_metadata$time_point != "2_D2",]
matching_all_samples_2_D8_binary_metadata <- matching_all_samples_2_D8_binary_metadata[matching_all_samples_2_D8_binary_metadata$time_point != "2_D5",]
matching_all_samples_2_D8_binary_metadata <- matching_all_samples_2_D8_binary_metadata[matching_all_samples_2_D8_binary_metadata$time_point != "2_D28",]
matching_all_samples_2_D8_binary_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, matching_all_samples_2_D8_binary_metadata, contrast = c("2_D_minus_1", "2_D8"))

matching_all_samples_2_D_minus_2_binary_metadata <- placebo_metadata[placebo_metadata$period == "2",]
matching_all_samples_2_D_minus_2_binary_metadata <- rbind(matching_all_samples_2_D_minus_2_binary_metadata, vaccinated_metadata[vaccinated_metadata$period == "2",])
matching_all_samples_2_D_minus_2_binary_metadata <- matching_all_samples_2_D_minus_2_binary_metadata[matching_all_samples_2_D_minus_2_binary_metadata$time_point != "2_D8",]
matching_all_samples_2_D_minus_2_binary_metadata <- matching_all_samples_2_D_minus_2_binary_metadata[matching_all_samples_2_D_minus_2_binary_metadata$time_point != "2_D2",]
matching_all_samples_2_D_minus_2_binary_metadata <- matching_all_samples_2_D_minus_2_binary_metadata[matching_all_samples_2_D_minus_2_binary_metadata$time_point != "2_D5",]
matching_all_samples_2_D_minus_2_binary_metadata <- matching_all_samples_2_D_minus_2_binary_metadata[matching_all_samples_2_D_minus_2_binary_metadata$time_point != "2_D28",]
matching_all_samples_2_D_minus_2_binary_metadata <- matching_all_samples_2_D_minus_2_binary_metadata[matching_all_samples_2_D_minus_2_binary_metadata$subject_id  %in% names(table(matching_all_samples_2_D_minus_2_binary_metadata$subject_id)[table(matching_all_samples_2_D_minus_2_binary_metadata$subject_id) == 2]),]
matching_all_samples_2_D_minus_2_binary_wayne_classifier <- apply_wayne_classifier(gene_counts_normalized_without_scale, matching_all_samples_2_D_minus_2_binary_metadata, contrast = c("2_D_minus_1", "2_D_minus_2"))






# Other approaches


# HIGH PLACEBO
hvl_placebo_second_period_metadata <- hvl_placebo_metadata[hvl_placebo_metadata$period == "2",]
# hvl_placebo_second_period_metadata <- hvl_placebo_second_period_metadata[hvl_placebo_second_period_metadata$time_point != "2_D_minus_2",]
hvl_placebo_second_period_training <- create_classifier_for_bulk_data(gene_counts_normalized_without_scale, hvl_placebo_second_period_metadata)
hvl_placebo_second_period_training_results <- hvl_placebo_second_period_training[[2]]

hvl_placebo_second_period_training_results_2_D_minus_1 <- hvl_placebo_second_period_training_results[hvl_placebo_second_period_training_results$max_score == "2_D_minus_1",]
hvl_placebo_2_D_minus_1_metadata <- hvl_placebo_second_period_metadata[hvl_placebo_second_period_metadata$time_point == "2_D_minus_1",]
hvl_placebo_second_period_training_results_2_D_minus_1$AVAL <- hvl_placebo_2_D_minus_1_metadata$AVAL

colnames(hvl_placebo_second_period_training_results_2_D_minus_1) <- c("time_2_D2", "time_2_D_minus_2", "time_2_D8", "time_2_D5", "time_2_D_minus_1", "time_2_D28", "max_score", "tie_flag", "AVAL")

ggplot(data = hvl_placebo_second_period_training_results_2_D_minus_1, mapping = aes(x = time_2_D28, y = AVAL)) +
  geom_point(size = 2) +
  sm_statCorr() + xlab("Misclassification Prob") + ylab("Viral Load") + labs(title = "Day")

# LOW PLACEBO
lvl_placebo_second_period_metadata <- lvl_placebo_metadata[lvl_placebo_metadata$period == "2",]
# lvl_placebo_second_period_metadata <- lvl_placebo_second_period_metadata[lvl_placebo_second_period_metadata$time_point != "2_D_minus_2",]
lvl_placebo_second_period_training <- create_classifier_for_bulk_data(gene_counts_normalized, lvl_placebo_second_period_metadata)
lvl_placebo_second_period_training_results <- lvl_placebo_second_period_training[[2]]

lvl_placebo_second_period_training_results_2_D_minus_1 <- lvl_placebo_second_period_training_results[lvl_placebo_second_period_training_results$max_score == "2_D_minus_1",]
lvl_placebo_2_D_minus_1_metadata <- lvl_placebo_second_period_metadata[lvl_placebo_second_period_metadata$time_point == "2_D_minus_1",]
lvl_placebo_second_period_training_results_2_D_minus_1$AVAL <- lvl_placebo_2_D_minus_1_metadata$AVAL

colnames(lvl_placebo_second_period_training_results_2_D_minus_1) <- c("time_2_D28", "time_2_D8", "time_2_D_minus_1", "time_2_D_minus_2", "time_2_D5", "time_2_D2", "max_score", "tie_flag", "AVAL")

ggplot(data = lvl_placebo_second_period_training_results_2_D_minus_1, mapping = aes(x = time_2_D28, y = AVAL)) +
  geom_point(size = 2) +
  sm_statCorr() + xlab("Misclassification Prob") + ylab("Viral Load") + labs(title = "Day")

# BOTH (HIGH AND LOW) PLACEBO
both_placebo_second_period_metadata <- both_placebo_metadata[both_placebo_metadata$period == "2",]
# both_placebo_second_period_metadata <- both_placebo_second_period_metadata[both_placebo_second_period_metadata$time_point != "2_D_minus_2",]
both_placebo_second_period_training <- create_classifier_for_bulk_data(gene_counts_normalized, both_placebo_second_period_metadata)
both_placebo_second_period_training_results <- both_placebo_second_period_training[[2]]

both_placebo_aliquot_d_minus_1_subset <- rownames(both_placebo_second_period_metadata[both_placebo_second_period_metadata$time_point == "2_D_minus_1",])
both_placebo_second_period_training_results_2_D_minus_1 <- both_placebo_second_period_training_results[rownames(both_placebo_second_period_training_results) %in% both_placebo_aliquot_d_minus_1_subset,]
both_placebo_2_D_minus_1_metadata <- both_placebo_second_period_metadata[both_placebo_second_period_metadata$time_point == "2_D_minus_1",]
both_placebo_second_period_training_results_2_D_minus_1$AVAL <- both_placebo_2_D_minus_1_metadata$AVAL

colnames(both_placebo_second_period_training_results_2_D_minus_1) <- c("time_2_D28", "time_2_D8", "time_2_D2", "time_2_D_minus_2", "time_2_D_minus_1", "time_2_D5", "max_score", "tie_flag", "AVAL")

ggplot(data = both_placebo_second_period_training_results_2_D_minus_1, mapping = aes(x = time_2_D5, y = AVAL)) +
  geom_point(size = 2) +
  sm_statCorr() + xlab("Misclassification Score") + ylab("Viral Load") + labs(title = "2_D2")

# ALL PLACEBO
all_placebo_second_period_metadata <- placebo_metadata[placebo_metadata$period == "2",]
all_placebo_second_period_metadata <- all_placebo_second_period_metadata[all_placebo_second_period_metadata$time_point != "2_D_minus_2",]
all_placebo_second_period_training <- create_classifier_for_bulk_data(gene_counts_normalized, all_placebo_second_period_metadata)
all_placebo_second_period_training_results <- all_placebo_second_period_training[[2]]

all_placebo_aliquot_d_minus_1_subset <- rownames(all_placebo_second_period_metadata[all_placebo_second_period_metadata$time_point == "2_D_minus_1",])
all_placebo_second_period_training_results_2_D_minus_1 <- all_placebo_second_period_training_results[rownames(all_placebo_second_period_training_results) %in% all_placebo_aliquot_d_minus_1_subset,]
all_placebo_2_D_minus_1_metadata <- all_placebo_second_period_metadata[all_placebo_second_period_metadata$time_point == "2_D_minus_1",]
all_placebo_second_period_training_results_2_D_minus_1$AVAL <- all_placebo_2_D_minus_1_metadata$AVAL

colnames(all_placebo_second_period_training_results_2_D_minus_1) <- c("time_2_D_minus_1", "time_2_D28", "time_2_D8", "time_2_D2", "time_2_D5", "max_score", "tie_flag", "AVAL")

ggplot(data = all_placebo_second_period_training_results_2_D_minus_1, mapping = aes(x = time_2_D28, y = AVAL)) +
  geom_point(size = 2) +
  sm_statCorr() + xlab("Misclassification Score") + ylab("Viral Load") + labs(title = "Day")

# ALL PLACEBO BUT MATCHING TIME POINTS
matching_all_placebo_second_period_metadata <- placebo_metadata[placebo_metadata$period == "2",]
matching_all_placebo_second_period_metadata <- matching_all_placebo_second_period_metadata[matching_all_placebo_second_period_metadata$time_point != "2_D_minus_2",]
matching_all_placebo_second_period_metadata <- matching_all_placebo_second_period_metadata[matching_all_placebo_second_period_metadata$time_point != "2_D2",]
matching_all_placebo_second_period_metadata <- matching_all_placebo_second_period_metadata[matching_all_placebo_second_period_metadata$time_point != "2_D5",]
matching_all_placebo_second_period_training <- create_classifier_for_bulk_data(gene_counts_normalized, matching_all_placebo_second_period_metadata)
matching_all_placebo_second_period_training_results <- matching_all_placebo_second_period_training[[2]]

matching_all_placebo_aliquot_d_minus_1_subset <- rownames(matching_all_placebo_second_period_metadata[matching_all_placebo_second_period_metadata$time_point == "2_D_minus_1",])
matching_all_placebo_second_period_training_results_2_D_minus_1 <- matching_all_placebo_second_period_training_results[rownames(matching_all_placebo_second_period_training_results) %in% matching_all_placebo_aliquot_d_minus_1_subset,]
matching_all_placebo_2_D_minus_1_metadata <- matching_all_placebo_second_period_metadata[matching_all_placebo_second_period_metadata$time_point == "2_D_minus_1",]
matching_all_placebo_second_period_training_results_2_D_minus_1$AVAL <- matching_all_placebo_2_D_minus_1_metadata$AVAL

colnames(matching_all_placebo_second_period_training_results_2_D_minus_1) <- c("time_2_D_minus_1", "time_2_D28", "time_2_D8", "max_score", "tie_flag", "AVAL")

ggplot(data = matching_all_placebo_second_period_training_results_2_D_minus_1, mapping = aes(x = time_2_D8, y = AVAL)) +
  geom_point(size = 2) +
  sm_statCorr() + xlab("Misclassification Score") + ylab("Viral Load") + labs(title = "Day")

ggplot(data = matching_all_placebo_second_period_training_results_2_D_minus_1[matching_all_placebo_second_period_training_results_2_D_minus_1$AVAL >= hvl_threshold,], 
       mapping = aes(x = time_2_D28, y = AVAL)) +
  geom_point(size = 2) +
  sm_statCorr() + xlab("Misclassification Score") + ylab("Viral Load") + labs(title = "Day")

# RANDOM FOREST TEST

# Set up classifier object
matching_all_placebo_second_period_metadata <- placebo_metadata[placebo_metadata$period == "2",]
matching_all_placebo_second_period_metadata <- matching_all_placebo_second_period_metadata[matching_all_placebo_second_period_metadata$time_point != "2_D_minus_2",]
matching_all_placebo_second_period_metadata <- matching_all_placebo_second_period_metadata[matching_all_placebo_second_period_metadata$time_point != "2_D2",]
matching_all_placebo_second_period_metadata <- matching_all_placebo_second_period_metadata[matching_all_placebo_second_period_metadata$time_point != "2_D5",]
matching_all_placebo_second_period_counts <- gene_counts_normalized[,rownames(matching_all_placebo_second_period_metadata)]
rownames(matching_all_placebo_second_period_counts) <- gsub("-", "_", rownames(matching_all_placebo_second_period_counts))
time_points <- matching_all_placebo_second_period_metadata$time_point
classes <- as.vector(unique(time_points))

classifier_object <- ReadData(Data = matching_all_placebo_second_period_counts,
                              Labels = time_points,
                              verbose = TRUE)

# Sort genes
genes_RF <- sort_genes_RF(data_object = classifier_object,
                          # featureNo_altogether, it is better not to specify a number here
                          # featureNo_one_vs_rest, it is better not to specify a number here
                          rank_data = TRUE,
                          platform_wise = FALSE,
                          num.trees = 1000, # more features, more tress are recommended
                          seed=get_speedi_seed(), # for reproducibility
                          verbose = TRUE)
genes_RF # sorted genes object

summary_genes <- summary_genes_RF(sorted_genes_RF = genes_RF,
                                  genes_altogether = c(10,20,50,100,150,200),
                                  genes_one_vs_rest = c(10,20,50,100,150,200))

rules_RF <- sort_rules_RF(data_object = object, 
                          sorted_genes_RF = genes_RF,
                          genes_altogether = 50,
                          genes_one_vs_rest = 50, 
                          num.trees = 1000,# more rules, more tress are recommended 
                          seed=get_speedi_seed(),
                          verbose = TRUE)

parameters <- data.frame(
  gene_repetition=c(3,2,1),
  rules_one_vs_rest=c(2,3,10),
  rules_altogether=c(2,3,10),
  run_boruta=c(FALSE,"make_error",TRUE), # I want to produce error in the 2nd trial
  plot_boruta = FALSE,
  num.trees=c(100,200,300),
  stringsAsFactors = FALSE)

# parameters
# for overall and byclass possible options, check the help files
para_opt <- optimize_RF(data_object = classifier_object,
                        sorted_rules_RF = rules_RF,
                        parameters = parameters,
                        test_object = NULL,
                        overall = c("Accuracy","Kappa"), # wanted overall measurements 
                        byclass = c("F1"), # wanted measurements per class
                        verbose = TRUE,
                        seed = get_speedi_seed())

RF_classifier <- train_RF(data_object = classifier_object,
                          sorted_rules_RF = rules_RF,
                          gene_repetition = 1,
                          rules_altogether = 10,
                          rules_one_vs_rest = 10,
                          run_boruta = TRUE, 
                          plot_boruta = FALSE,
                          probability = TRUE,
                          num.trees = 300,
                          boruta_args = list(),
                          verbose = TRUE)

proximity_matrix_RF(object = classifier_object,
                    classifier = RF_classifier, 
                    plot = TRUE,
                    return_matrix = FALSE, # if we need to get the matrix itself
                    title = "Test",
                    cluster_cols = TRUE)

training_pred <- RF_classifier$RF_scheme$RF_classifier$predictions
if (is.factor(training_pred)) {
  x <- as.character(training_pred)
}

# if the classifier trained using probability   = TRUE
if (is.matrix(training_pred)) {
  x <- colnames(training_pred)[max.col(training_pred)]
}

# training accuracy
caret::confusionMatrix(data =factor(x),
                       reference = factor(classifier_object$data$Labels),
                       mode = "everything")

training_pred <- as.data.frame(training_pred)

training_pred$correct_time_point <- time_points
training_pred$AVAL <- matching_all_placebo_second_period_metadata$AVAL


training_pred_d_minus_1 <- training_pred[training_pred$correct_time_point == "2_D_minus_1",]
colnames(training_pred_d_minus_1) <- c("time_2_D_minus_1", "time_2_D28", "time_2_D8", "correct_time_point", "AVAL")

ggplot(data = training_pred_d_minus_1, mapping = aes(x = time_2_D28, y = AVAL)) +
  geom_point(size = 2) +
  sm_statCorr() + xlab("Misclassification Prob") + ylab("Viral Load") + labs(title = "Day")