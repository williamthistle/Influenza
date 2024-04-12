# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# Wayne approach
# TODO: Try binary D minus 1 vs D28

# High placebo
high_placebo_second_period_metadata <- high_placebo_metadata[high_placebo_metadata$period == "2",]
high_placebo_second_period_metadata <- high_placebo_second_period_metadata[high_placebo_second_period_metadata$time_point != "2_D_minus_2",]
high_placebo_second_period_counts <- gene_counts_normalized[,rownames(high_placebo_second_period_metadata)]
rownames(high_placebo_second_period_counts) <- gsub("-", "_", rownames(high_placebo_second_period_counts))
high_placebo_second_period_counts <- t(high_placebo_second_period_counts)
high_placebo_second_period_label <- as.vector(high_placebo_second_period_metadata$time_point)

contrast <- c("2_D_minus_1", "2_D2", "2_D5", "2_D8", "2_D28")
high_placebo_cv_report <- cv_report(high_placebo_second_period_counts, high_placebo_second_period_label, task = "classification_multi", contrast = contrast)
high_placebo_cv_report_predictions <- as.data.frame(high_placebo_cv_report[[1]])
high_placebo_cv_report_predictions$label <- high_placebo_second_period_label
high_placebo_cv_report_predictions$AVAL <- high_placebo_second_period_metadata$AVAL
colnames(high_placebo_cv_report_predictions) <- c("time_2_D_minus_1", "time_2_D2", "time_2_D5", "time_2_D8", "time_2_D28", "correct_label", "AVAL")

high_placebo_cv_report_predictions_d_minus_1 <- high_placebo_cv_report_predictions[high_placebo_cv_report_predictions$correct_label == "2_D_minus_1",]

ggplot(data = high_placebo_cv_report_predictions_d_minus_1, mapping = aes(x = time_2_D28, y = AVAL)) +
  geom_point(size = 2) +
  sm_statCorr() + xlab("Misclassification Prob") + ylab("Viral Load") + labs(title = "Day")


# All matched placebo
matching_all_placebo_second_period_metadata <- placebo_metadata[placebo_metadata$period == "2",]
matching_all_placebo_second_period_metadata <- matching_all_placebo_second_period_metadata[matching_all_placebo_second_period_metadata$time_point != "2_D_minus_2",]
matching_all_placebo_second_period_metadata <- matching_all_placebo_second_period_metadata[matching_all_placebo_second_period_metadata$time_point != "2_D2",]
matching_all_placebo_second_period_metadata <- matching_all_placebo_second_period_metadata[matching_all_placebo_second_period_metadata$time_point != "2_D5",]

matching_all_placebo_second_period_counts <- gene_counts_normalized[,rownames(matching_all_placebo_second_period_metadata)]
rownames(matching_all_placebo_second_period_counts) <- gsub("-", "_", rownames(matching_all_placebo_second_period_counts))
matching_all_placebo_second_period_counts <- t(matching_all_placebo_second_period_counts)
matching_all_placebo_second_period_label <- as.vector(matching_all_placebo_second_period_metadata$time_point)

contrast <- c("2_D_minus_1", "2_D8", "2_D28")
matching_all_placebo_cv_report <- cv_report(matching_all_placebo_second_period_counts, matching_all_placebo_second_period_label, task = "classification_multi", contrast = contrast)
matching_all_placebo_cv_report_predictions <- as.data.frame(matching_all_placebo_cv_report[[1]])
matching_all_placebo_cv_report_predictions$label <- matching_all_placebo_second_period_label
matching_all_placebo_cv_report_predictions$AVAL <- matching_all_placebo_second_period_metadata$AVAL
colnames(matching_all_placebo_cv_report_predictions) <- c("time_2_D_minus_1", "time_2_D8", "time_2_D28", "correct_label", "AVAL")

matching_all_placebo_cv_report_predictions_d_minus_1 <- matching_all_placebo_cv_report_predictions[matching_all_placebo_cv_report_predictions$correct_label == "2_D_minus_1",]

ggplot(data = matching_all_placebo_cv_report_predictions_d_minus_1, mapping = aes(x = time_2_D_minus_1, y = AVAL)) +
  geom_point(size = 2) +
  sm_statCorr() + xlab("Misclassification Prob") + ylab("Viral Load") + labs(title = "Day")

# All matched vaccine
matching_all_vaccinated_second_period_metadata <- vaccinated_metadata[vaccinated_metadata$period == "2",]
matching_all_vaccinated_second_period_metadata <- matching_all_vaccinated_second_period_metadata[matching_all_vaccinated_second_period_metadata$time_point != "2_D_minus_2",]
matching_all_vaccinated_second_period_metadata <- matching_all_vaccinated_second_period_metadata[matching_all_vaccinated_second_period_metadata$time_point != "2_D2",]
matching_all_vaccinated_second_period_metadata <- matching_all_vaccinated_second_period_metadata[matching_all_vaccinated_second_period_metadata$time_point != "2_D5",]

matching_all_vaccinated_second_period_counts <- gene_counts_normalized[,rownames(matching_all_vaccinated_second_period_metadata)]
rownames(matching_all_vaccinated_second_period_counts) <- gsub("-", "_", rownames(matching_all_vaccinated_second_period_counts))
matching_all_vaccinated_second_period_counts <- t(matching_all_vaccinated_second_period_counts)
matching_all_vaccinated_second_period_label <- as.vector(matching_all_vaccinated_second_period_metadata$time_point)

contrast <- c("2_D_minus_1", "2_D8", "2_D28")
matching_all_vaccinated_cv_report <- cv_report(matching_all_vaccinated_second_period_counts, matching_all_vaccinated_second_period_label, task = "classification_multi", contrast = contrast)
matching_all_vaccinated_cv_report_predictions <- as.data.frame(matching_all_vaccinated_cv_report[[1]])
matching_all_vaccinated_cv_report_predictions$label <- matching_all_vaccinated_second_period_label
matching_all_vaccinated_cv_report_predictions$AVAL <- matching_all_vaccinated_second_period_metadata$AVAL
colnames(matching_all_vaccinated_cv_report_predictions) <- c("time_2_D_minus_1", "time_2_D8", "time_2_D28", "correct_label", "AVAL")

matching_all_vaccinated_cv_report_predictions_d_minus_1 <- matching_all_vaccinated_cv_report_predictions[matching_all_vaccinated_cv_report_predictions$correct_label == "2_D_minus_1",]

ggplot(data = matching_all_vaccinated_cv_report_predictions_d_minus_1, mapping = aes(x = time_2_D28, y = AVAL)) +
  geom_point(size = 2) +
  sm_statCorr() + xlab("Misclassification Prob") + ylab("Viral Load") + labs(title = "Day")

# All samples
matching_all_samples_second_period_metadata <- rbind(matching_all_placebo_second_period_metadata, matching_all_vaccinated_second_period_metadata)
matching_all_samples_second_period_counts <- rbind(matching_all_placebo_second_period_counts, matching_all_vaccinated_second_period_counts)
matching_all_samples_second_period_label <- as.vector(matching_all_samples_second_period_metadata$time_point)

contrast <- c("2_D_minus_1", "2_D8", "2_D28")
matching_all_samples_cv_report <- cv_report(matching_all_samples_second_period_counts, matching_all_samples_second_period_label, task = "classification_multi", contrast = contrast)
matching_all_samples_cv_report_predictions <- as.data.frame(matching_all_samples_cv_report[[1]])
matching_all_samples_cv_report_predictions$label <- matching_all_samples_second_period_label
matching_all_samples_cv_report_predictions$AVAL <- matching_all_samples_second_period_metadata$AVAL
colnames(matching_all_samples_cv_report_predictions) <- c("time_2_D_minus_1", "time_2_D8", "time_2_D28", "correct_label", "AVAL")

matching_all_samples_cv_report_predictions_d_minus_1 <- matching_all_samples_cv_report_predictions[matching_all_samples_cv_report_predictions$correct_label == "2_D_minus_1",]

ggplot(data = matching_all_samples_cv_report_predictions_d_minus_1, mapping = aes(x = time_2_D28, y = AVAL)) +
  geom_point(size = 2) +
  sm_statCorr() + xlab("Misclassification Prob") + ylab("Viral Load") + labs(title = "Day")

# All matched placebo (binary)
matching_all_placebo_second_period_metadata <- placebo_metadata[placebo_metadata$period == "2",]
matching_all_placebo_second_period_metadata <- matching_all_placebo_second_period_metadata[matching_all_placebo_second_period_metadata$time_point != "2_D_minus_2",]
matching_all_placebo_second_period_metadata <- matching_all_placebo_second_period_metadata[matching_all_placebo_second_period_metadata$time_point != "2_D2",]
matching_all_placebo_second_period_metadata <- matching_all_placebo_second_period_metadata[matching_all_placebo_second_period_metadata$time_point != "2_D5",]
matching_all_placebo_second_period_metadata <- matching_all_placebo_second_period_metadata[matching_all_placebo_second_period_metadata$time_point != "2_D8",]

matching_all_placebo_second_period_counts <- gene_counts_normalized[,rownames(matching_all_placebo_second_period_metadata)]
rownames(matching_all_placebo_second_period_counts) <- gsub("-", "_", rownames(matching_all_placebo_second_period_counts))
matching_all_placebo_second_period_counts <- t(matching_all_placebo_second_period_counts)
matching_all_placebo_second_period_label <- as.vector(matching_all_placebo_second_period_metadata$time_point)

contrast <- c("2_D_minus_1", "2_D28")
matching_all_placebo_cv_report <- cv_report(matching_all_placebo_second_period_counts, matching_all_placebo_second_period_label, task = "classification", contrast = contrast)
matching_all_placebo_cv_report_predictions <- as.data.frame(matching_all_placebo_cv_report[[1]])
matching_all_placebo_cv_report_predictions$label <- matching_all_placebo_second_period_label
matching_all_placebo_cv_report_predictions$AVAL <- matching_all_placebo_second_period_metadata$AVAL
colnames(matching_all_placebo_cv_report_predictions) <- c("time_2_D_minus_1", "time_2_D8", "time_2_D28", "correct_label", "AVAL")

matching_all_placebo_cv_report_predictions_d_minus_1 <- matching_all_placebo_cv_report_predictions[matching_all_placebo_cv_report_predictions$correct_label == "2_D_minus_1",]

ggplot(data = matching_all_placebo_cv_report_predictions_d_minus_1, mapping = aes(x = time_2_D_minus_1, y = AVAL)) +
  geom_point(size = 2) +
  sm_statCorr() + xlab("Misclassification Prob") + ylab("Viral Load") + labs(title = "Day")

# Other approaches

create_classifier_for_bulk_data <- function(counts, metadata) {
  counts <- counts[,rownames(metadata)]
  rownames(counts) <- gsub("-", "_", rownames(counts))
  time_points <- metadata$time_point
  classes <- as.vector(unique(time_points))
  
  # Create classifier object
  classifier_object <- ReadData(Data = counts,
                                Labels = time_points,
                                verbose = TRUE)
  
  # Get list of filtered genes
  filtered_genes <- filter_genes_TSP(data_object = classifier_object,
                                     filter = "one_vs_one",
                                     platform_wise = FALSE,
                                     featureNo = 1000,
                                     UpDown = TRUE,
                                     verbose = TRUE)
  
  # Train classifier
  classifier <- train_one_vs_rest_TSP(data_object = classifier_object,
                                      filtered_genes = filtered_genes,
                                      k_range = 5:50,
                                      include_pivot = FALSE,
                                      one_vs_one_scores = TRUE,
                                      platform_wise_scores = FALSE,
                                      seed = get_speedi_seed(),
                                      verbose = TRUE)
  
  # Test classifier on training data
  results_train <- predict_one_vs_rest_TSP(classifier = classifier,
                                           Data = classifier_object,
                                           tolerate_missed_genes = TRUE,
                                           weighted_votes = TRUE,
                                           classes = classes,
                                           verbose = TRUE)
  
  # Check that predictions are accurate
  results_train$max_score == time_points
  
  caret::confusionMatrix(data = factor(results_train$max_score, 
                                       levels = unique(classifier_object$data$Labels)),
                         reference = factor(classifier_object$data$Labels, 
                                            levels = unique(classifier_object$data$Labels)),
                         mode="everything")
  return(list(classifier, results_train))
}

# HIGH PLACEBO
high_placebo_second_period_metadata <- high_placebo_metadata[high_placebo_metadata$period == "2",]
# high_placebo_second_period_metadata <- high_placebo_second_period_metadata[high_placebo_second_period_metadata$time_point != "2_D_minus_2",]
high_placebo_second_period_training <- create_classifier_for_bulk_data(gene_counts_normalized, high_placebo_second_period_metadata)
high_placebo_second_period_training_results <- high_placebo_second_period_training[[2]]

high_placebo_second_period_training_results_2_D_minus_1 <- high_placebo_second_period_training_results[high_placebo_second_period_training_results$max_score == "2_D_minus_1",]
high_placebo_2_D_minus_1_metadata <- high_placebo_second_period_metadata[high_placebo_second_period_metadata$time_point == "2_D_minus_1",]
high_placebo_second_period_training_results_2_D_minus_1$AVAL <- high_placebo_2_D_minus_1_metadata$AVAL

colnames(high_placebo_second_period_training_results_2_D_minus_1) <- c("time_2_D2", "time_2_D_minus_2", "time_2_D8", "time_2_D5", "time_2_D_minus_1", "time_2_D28", "max_score", "tie_flag", "AVAL")

ggplot(data = high_placebo_second_period_training_results_2_D_minus_1, mapping = aes(x = time_2_D28, y = AVAL)) +
  geom_point(size = 2) +
  sm_statCorr() + xlab("Misclassification Prob") + ylab("Viral Load") + labs(title = "Day")

# LOW PLACEBO
low_placebo_second_period_metadata <- low_placebo_metadata[low_placebo_metadata$period == "2",]
# low_placebo_second_period_metadata <- low_placebo_second_period_metadata[low_placebo_second_period_metadata$time_point != "2_D_minus_2",]
low_placebo_second_period_training <- create_classifier_for_bulk_data(gene_counts_normalized, low_placebo_second_period_metadata)
low_placebo_second_period_training_results <- low_placebo_second_period_training[[2]]

low_placebo_second_period_training_results_2_D_minus_1 <- low_placebo_second_period_training_results[low_placebo_second_period_training_results$max_score == "2_D_minus_1",]
low_placebo_2_D_minus_1_metadata <- low_placebo_second_period_metadata[low_placebo_second_period_metadata$time_point == "2_D_minus_1",]
low_placebo_second_period_training_results_2_D_minus_1$AVAL <- low_placebo_2_D_minus_1_metadata$AVAL

colnames(low_placebo_second_period_training_results_2_D_minus_1) <- c("time_2_D28", "time_2_D8", "time_2_D_minus_1", "time_2_D_minus_2", "time_2_D5", "time_2_D2", "max_score", "tie_flag", "AVAL")

ggplot(data = low_placebo_second_period_training_results_2_D_minus_1, mapping = aes(x = time_2_D28, y = AVAL)) +
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