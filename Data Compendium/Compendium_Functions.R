library(org.Hs.eg.db)

# This function takes a MetaIntegrator object and logical index vector indicating
# which samples to keep/remove, and returns a MetaIntegrator object containing 
# the indicated subset 
filterMIobj <- function(MIobj, index){
  tmpMIobj <- MIobj
  tmpMIobj$expr <- tmpMIobj$expr[, index] # filter expression
  tmpMIobj$pheno <- tmpMIobj$pheno[index, ] # filter pheno metadata
  tmpMIobj$class <- tmpMIobj$class[index] # filter class vector
  
  if(checkDataObject(tmpMIobj, 'Dataset')){ # verify that the resulting object is still a MetaIntegrator object
    return(tmpMIobj)
  } else {
    stop()
  }
}

getClassSum <- function(MIobj, input_class, minimum_N = 4){
  class_match <- grepl(pattern = input_class, x = MIobj$pheno$Class)
  return(ifelse(sum(class_match) > minimum_N, T, F))
}

getPathogenSum <- function(MIobj, input_class, minimum_N = 4){
  path_match <- grepl(pattern = input_class, x = MIobj$pheno$Pathogen)
  return(ifelse(sum(path_match) > minimum_N, T, F))
}

# for an input MIobj and  pathogen (character)
# this function returns a logical vector indicating whether a dataset 
# contains the pathogen of interest
findPathogen <- function(MIobj, pathogen){
  pathogen_vector <- MIobj$pheno$Pathogen
  pathogen_vector <- pathogen_vector[!is.na(pathogen_vector)]
  if(sum(pathogen_vector == pathogen) > 0){
    return(T)
  } else {
    return(F)
  }
}

# Create MetaIntegrator friendly object for processing
formatContrast <- function(MIobj, target_class) {
  keep_ix <- !is.na(MIobj$pheno$Class)
  if(sum(keep_ix) == 0){print(paste0(MIobj$formattedName, ' has no valid samples')); return(NULL)}
  MIobj <- filterMIobj(MIobj, keep_ix)
  keep_ix <- MIobj$pheno$Class %in% c("Healthy", "Convalescent", target_class)
  if(sum(keep_ix) == 0){print(paste0(MIobj$formattedName, ' has no valid samples')); return(NULL)}
  MIobj <- filterMIobj(MIobj, keep_ix)
  contrast_vector <- as.numeric(!MIobj$pheno$Class %in% c('Healthy', 'Convalescent'))
  names(contrast_vector) <- rownames(MIobj$pheno)
  MIobj$class <- contrast_vector
  return(MIobj)
}

# for an input MIobj and pathogen (character),
# this function returns the subset of subjects infected by the pathogen
# as well as healthy controls, with a modified $class vector for use with
# MetaIntegrator functions
formatPathogenContrast <- function(MIobj, pathogen){
  #print(MIobj$formattedName)
  keep_ix <- !is.na(MIobj$pheno$Pathogen) & !is.na(MIobj$pheno$Class)
  if(sum(keep_ix) == 0){print(paste0(MIobj$formattedName, ' has no valid samples')); return(NULL)}
  MIobj <- filterMIobj(MIobj, keep_ix)
  keep_ix <- MIobj$pheno$Pathogen == pathogen | MIobj$pheno$Class %in% c("Healthy", "Convalescent")
  if(sum(keep_ix) == 0){print(paste0(MIobj$formattedName, ' has no valid samples')); return(NULL)}
  MIobj <- filterMIobj(MIobj, keep_ix)
  contrast_vector <- as.numeric(!MIobj$pheno$Class %in% c('Healthy', 'Convalescent'))
  names(contrast_vector) <- rownames(MIobj$pheno)
  MIobj$class <- contrast_vector
  return(MIobj)
}

# return a logical indicating whether MIobj contains at least N cases and N controls
filterSampleSizes <- function(MIobj, N = 4){
  if(is.null(MIobj)){return(F)}
  keep_ix <- !is.na(MIobj$pheno$Class) & !is.na(MIobj$class)
  MIobj <- filterMIobj(MIobj, keep_ix)
  return(length(unique(MIobj$class)) == 2 & min(table(MIobj$pheno$Class)) >= N)
}

# return number of samples in a dataset
countSampleSizes <- function(MIobj){
  if(is.null(MIobj)){return(F)}
  keep_ix <- !is.na(MIobj$pheno$Class) & !is.na(MIobj$class)
  print(nrow(MIobj$pheno))
}




# for an input MIobj, pathogen (character), and target_class (character)
# this function returns the subset of subjects with infections of type target_class that are NOT infected by the pathogen
# as well as healthy controls, with a modified $class vector for use with
# MetaIntegrator functions
removePathogen <- function(MIobj, pathogen, target_class){
  #print(MIobj$formattedName)
  keep_ix <- !is.na(MIobj$pheno$Class)
  if(sum(keep_ix) == 0){print(paste0(MIobj$formattedName, ' has no valid samples')); return(NULL)}
  MIobj <- filterMIobj(MIobj, keep_ix)
  keep_ix <- MIobj$pheno$Pathogen != pathogen & MIobj$pheno$Class %in% c("Healthy", "Convalescent", target_class)
  if(sum(keep_ix) == 0){print(paste0(MIobj$formattedName, ' has no valid samples')); return(NULL)}
  MIobj <- filterMIobj(MIobj, keep_ix)
  contrast_vector <- as.numeric(!MIobj$pheno$Class %in% c('Healthy', 'Convalescent'))
  names(contrast_vector) <- rownames(MIobj$pheno)
  MIobj$class <- contrast_vector
  return(MIobj)
}

calculateAUROC <- function(MIobj, signature, method = 'geomMean', suppressMessages = F){
  labels <- MIobj$class
  scores <- calculateScoreRobust(filterObject = signature, datasetObject = MIobj, method = method, suppressMessages = suppressMessages)
  if(all(sapply(scores, is.na))){
    roc <- NA
  } else {
    roc <- calculateROC(labels = labels, predictions = scores, AUConly = T)
  }
  return(as.numeric(roc))
}

calculateScoreRobust <- function(filterObject, datasetObject, suppressMessages = T, method = 'geomMean', zScore = T){
  if(!method %in% c('original', 'geomMean')){
    stop('method must be original or geomMean')
  }
  
  if(method == 'original'){
    totalScore = calculateScore(filterObject = filterObject, datasetObject = datasetObject, suppressMessages = suppressMessages)
  } else {
    # method options: geomMean, geomMean_exp2, mean, mean_exp2
    # assumes data is being entered on a log2 scale
    
    datasetObjectmin <- min(datasetObject$expr, na.rm = TRUE)
    if (datasetObjectmin < 0) {
      datasetObject$expr <- datasetObject$expr + abs(datasetObjectmin) + 1
    }
    
    # pull list of genes present in data set
    expr_genes <- datasetObject$keys
    filterObject$posGeneNames <- setdiff(filterObject$posGeneNames,'')
    filterObject$negGeneNames <- setdiff(filterObject$negGeneNames,'')
    pos_genes = intersect(filterObject$posGeneNames, expr_genes)
    neg_genes = intersect(filterObject$negGeneNames, expr_genes)
    
    # for reporting fraction of genes used
    if (!suppressMessages){
      N_pos_str <- paste0(length(pos_genes), ' of ', length(filterObject$posGeneNames))
      N_neg_str <- paste0(length(neg_genes), ' of ', length(filterObject$negGeneNames))
      print(paste0('Used ', N_pos_str, ' pos genes and ', N_neg_str, ' neg genes'))
    }
    
    
    if(length(pos_genes) == 0 & length(neg_genes) == 0){
      if (!suppressMessages){
        print(filterObject$posGeneNames)
        print(filterObject$negGeneNames)
        print('no common genes between data and signature')
      }
      return(NA)
    }
    
    posScore = .calculate_scores(datasetObject, pos_genes, method = method)
    negScore = .calculate_scores(datasetObject, neg_genes, method = method)
    
    if(grepl(pattern = 'exp2', x = method)){
      totalScore = posScore/negScore
    } else {
      totalScore = posScore - negScore
    }
    
    
    if (sum(abs(totalScore), na.rm = T) != 0){
      if(zScore){
        totalScore = as.numeric(scale(totalScore))
      } else {
        totalScore = as.numeric(totalScore)
      }
      
    }
  }
  
  return(totalScore)
}


# implemented w input from Antonio Cappuccio!
.calculate_scores = function(datasetObject, genes, method = 'geomMean'){
  
  # find the genes to average
  genes_idx = which(datasetObject$keys %in% genes)
  
  # if there aren't any genes, return 0
  if (length(genes_idx) == 0){
    if(method == 'mean'){
      output_val = 1
    } else {
      output_val = 0
    }
    
    return(rep(output_val, ncol(datasetObject$expr)))
  }
  
  # filter expr to sig genes
  gene_expr = datasetObject$expr[genes_idx,] 

  #in case only one gene is present, take the vector itself
  if (is.null(nrow(gene_expr))){
    sample_scores = gene_expr
  }
  
  if (!is.null(nrow(gene_expr))){
    sample_scores = apply(gene_expr, 2, function(method, x){switch(method, 
                                                                   geomMean = geom_mean(x),
                                                                   geomMean_exp2 = geom_mean_exp2(x),
                                                                   mean = mean(x, na.rm = T),
                                                                   mean_exp2 = mean_exp2(x))
    }, method = method)}
  
  return(sample_scores)
}

# Generate AUCS for individual genes on a list of datasets
# I use source to declare the data source for the genes I'm testing (e.g., single cell, multiome) 
# and disease_tag to provide a little more detail on what I'm testing on (e.g., bulk_day_2)
test_individual_genes_on_datasets <- function(gene_list, data_list, source, disease_tag) {
  gene_aucs <- c()
  # Sometimes genes aren't found in our ENTREZ mapping database - for those, we just remove them
  final_gene_list <- c()
  for(gene in gene_list) {
    if(gene %in% keys(org.Hs.eg.db, keytype = "SYMBOL")) {
      final_gene_list <- c(final_gene_list, gene)
      sig <- list()
      # Convert gene name to ENTREZ ID
      gene_name <- mapIds(org.Hs.eg.db, c(gene), "ENTREZID", "SYMBOL")[1]
      print(gene_name)
      sig$posGeneNames <- gene_name
      sig$negGeneNames <- ''
      sig$filterDescription <- 'test'
      sig$FDRThresh <- 0
      sig$effectSizeThresh <- 0
      sig$numberStudiesThresh <- 1
      sig$isLeaveOneOut <- F
      sig$heterogeneityPvalThresh <- 0
      sig$timestamp <- Sys.time()
      # Test gene set signature on samples and print median AUROC
      flu_aucs <- sapply(X = data_list, FUN = calculateAUROC, signature = sig)
      gene_aucs <- c(gene_aucs, median(flu_aucs, na.rm = TRUE))
    } else {
      print(paste0("Gene ", gene, " not found in database"))
    }
  }
  gene_auc_name <- paste0(disease_tag, "_gene_auc")
  final_df <- data.frame("gene_name" = final_gene_list, temp_name = gene_aucs, "source" = rep(source, length(final_gene_list)))
  names(final_df)[names(final_df) == "temp_name"] <- gene_auc_name
  return(final_df)
}

# Create a signature based on positive and negative gene lists and find pooled AUC on a collection of datasets
plot_pooled_auc <- function(pos_gene_list, neg_gene_list, data_list, title) {
  meta_obj <- list()
  meta_obj$originalData <- data_list
  
  sig <- list()
  # Convert gene names to ENTREZ ID
  entrez_pos_gene_list <- mapIds(org.Hs.eg.db, pos_gene_list, "ENTREZID", "SYMBOL")
  sig$posGeneNames <- entrez_pos_gene_list
  if(length(neg_gene_list) > 0) {
    entrez_neg_gene_list <- mapIds(org.Hs.eg.db, neg_gene_list, "ENTREZID", "SYMBOL")
    sig$negGeneNames <- entrez_neg_gene_list
  } else {
    sig$negGeneNames <- ''
  }
  sig$filterDescription <- 'my_gene_list'
  sig$FDRThresh <- 0
  sig$effectSizeThresh <- 0
  sig$numberStudiesThresh <- 1
  sig$isLeaveOneOut <- F
  sig$heterogeneityPvalThresh <- 0
  sig$timestamp <- Sys.time()
  
  final_plot <- pooledROCPlot(metaObject = meta_obj, filterObject = sig, title = title)
  ggsave(filename = paste0("plots/pooled_auc/", title, ".png"), plot = final_plot, device = "png", dpi = 300)
}

geom_mean <- function(x) {
  return(exp(mean(log(x))))
}

# We will fix any datasets that need fixing
fix_data_list <- function(data_list) {
  current_dataset <- data_list[[127]]
  current_dataset$pheno[grepl("acute", current_dataset$pheno$title),]$Study <- "GSE97741_GPL10558"
  current_dataset$pheno[grepl("acute", current_dataset$pheno$title),]$Class <- "Virus"
  current_dataset$pheno[grepl("acute", current_dataset$pheno$title),]$Pathogen <- "Respiratory syncytial virus;Rhinovirus"
  current_dataset$pheno[grepl("discharge", current_dataset$pheno$title),]$Study <- "GSE97741_GPL10558"
  current_dataset$pheno[grepl("discharge", current_dataset$pheno$title),]$Class <- "Convalescent"
  current_dataset$pheno[grepl("discharge", current_dataset$pheno$title),]$Pathogen <- "Healthy"
  data_list[[127]] <- current_dataset
  return(data_list)
}

split_flu_data <- function(flu_list, study_metadata, daniel_paper_datasets, usage_token) {
  current_flu_list <- list()
  current_flu_dataset_names <- c()
  current_index <- 1
  daniel_paper_discovery_datasets <- daniel_paper_datasets[daniel_paper_datasets$Usage == usage_token,]
  for(current_flu_dataset in flu_list) {
    current_accession <- current_flu_dataset$formattedName
    for(daniel_accession_index in 1:nrow(daniel_paper_discovery_datasets)) {
      if(grepl(daniel_paper_discovery_datasets[daniel_accession_index,]$Accession, current_accession) & !(current_accession %in% current_flu_dataset_names)) {
        current_flu_dataset_names <- c(current_flu_dataset_names, current_accession)
        current_flu_list[[current_index]] <- current_flu_dataset
        current_index <- current_index + 1
      }
    }
  }
  current_flu_metadata <- study_metadata[study_metadata$Study %in% current_flu_dataset_names,]
  return(list(current_flu_list, current_flu_metadata))
}

# Need to filter time series data so it only contains most acute time point and maybe one healthy time point (convalescent / healthy)
fix_time_series_for_flu_discovery <- function(flu_discovery_list) {
  # Dataset 5
  kept_indices_5 <- c(1, 5, 10, 11, 12, 13, 16, 19, 22, 26, 27, 28, 29, 30, 31, 32, 33, 35, 36, 37, 38, 39, 40)
  kept_accessions_5 <- flu_discovery_list[[5]]$pheno$geo_accession[kept_indices_5]
  flu_discovery_list[[5]]$pheno <- flu_discovery_list[[5]]$pheno[flu_discovery_list[[5]]$pheno$geo_accession %in% kept_accessions_5,]
  flu_discovery_list[[5]]$expr <- flu_discovery_list[[5]]$expr[,colnames(flu_discovery_list[[5]]$expr) %in% kept_accessions_5]
  flu_discovery_list[[5]]$class <- flu_discovery_list[[5]]$class[names(flu_discovery_list[[5]]$class) %in% kept_accessions_5]
  # Dataset 6
  kept_indices_6 <- c(1, 5, 10, 15, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37)
  kept_accessions_6 <- flu_discovery_list[[6]]$pheno$geo_accession[kept_indices_6]
  flu_discovery_list[[6]]$pheno <- flu_discovery_list[[6]]$pheno[flu_discovery_list[[6]]$pheno$geo_accession %in% kept_accessions_6,]
  kept_accessions_6_alternate <- rownames(flu_discovery_list[[6]]$pheno)
  flu_discovery_list[[6]]$expr <- flu_discovery_list[[6]]$expr[,colnames(flu_discovery_list[[6]]$expr) %in% kept_accessions_6_alternate]
  flu_discovery_list[[6]]$class <- flu_discovery_list[[6]]$class[names(flu_discovery_list[[6]]$class) %in% kept_accessions_6_alternate]
  # Dataset 7
  flu_discovery_list[[7]]$pheno <- flu_discovery_list[[7]]$pheno[flu_discovery_list[[7]]$pheno$`Time.Point` != "late period",]
  kept_accessions_7 <- rownames(flu_discovery_list[[7]]$pheno)
  flu_discovery_list[[7]]$expr <- flu_discovery_list[[7]]$expr[,colnames(flu_discovery_list[[7]]$expr) %in% kept_accessions_7]
  flu_discovery_list[[7]]$class <- flu_discovery_list[[7]]$class[names(flu_discovery_list[[7]]$class) %in% kept_accessions_7]
  # Dataset 8
  flu_discovery_list[[8]]$pheno <- flu_discovery_list[[8]]$pheno[flu_discovery_list[[8]]$pheno$`Time.Point` == "Less than 72 hour after symptoms",]
  flu_discovery_list[[8]]$pheno <- flu_discovery_list[[8]]$pheno[flu_discovery_list[[8]]$pheno$Pathogen == "Influenza virus" | flu_discovery_list[[8]]$pheno$Pathogen == "Healthy",]
  kept_accessions_8 <- rownames(flu_discovery_list[[8]]$pheno)
  flu_discovery_list[[8]]$expr <- flu_discovery_list[[8]]$expr[,colnames(flu_discovery_list[[8]]$expr) %in% kept_accessions_8]
  flu_discovery_list[[8]]$class <- flu_discovery_list[[8]]$class[names(flu_discovery_list[[8]]$class) %in% kept_accessions_8]
  # Dataset 9
  flu_discovery_list[[9]]$pheno <- flu_discovery_list[[9]]$pheno[flu_discovery_list[[9]]$pheno$`Time.Point` == "Baseline" | flu_discovery_list[[9]]$pheno$`Time.Point` == "Day0",]
  kept_accessions_9 <- rownames(flu_discovery_list[[9]]$pheno)
  flu_discovery_list[[9]]$expr <- flu_discovery_list[[9]]$expr[,colnames(flu_discovery_list[[9]]$expr) %in% kept_accessions_9]
  flu_discovery_list[[9]]$class <- flu_discovery_list[[9]]$class[names(flu_discovery_list[[9]]$class) %in% kept_accessions_9]
  # Dataset 10
  flu_discovery_list[[10]]$pheno <- flu_discovery_list[[10]]$pheno[flu_discovery_list[[10]]$pheno$`Time.Point` == "Pre-Exposure" | flu_discovery_list[[10]]$pheno$`Time.Point` == "48 hours",]
  kept_accessions_10 <- rownames(flu_discovery_list[[10]]$pheno)
  flu_discovery_list[[10]]$expr <- flu_discovery_list[[10]]$expr[,colnames(flu_discovery_list[[10]]$expr) %in% kept_accessions_10]
  flu_discovery_list[[10]]$class <- flu_discovery_list[[10]]$class[names(flu_discovery_list[[10]]$class) %in% kept_accessions_10]
  return(flu_discovery_list)
}

# Need to filter time series data so it only contains most acute time point and maybe one healthy time point (convalescent / healthy)
fix_time_series_for_flu_validation <- function(flu_validation_list) {
  # Dataset 3
  flu_validation_list[[3]]$pheno <- flu_validation_list[[3]]$pheno[flu_validation_list[[3]]$pheno$`Time.Point` != "T2" & flu_validation_list[[3]]$pheno$`Time.Point` != "T3",]
  kept_accessions_3 <- rownames(flu_validation_list[[3]]$pheno)
  flu_validation_list[[3]]$expr <- flu_validation_list[[3]]$expr[,colnames(flu_validation_list[[3]]$expr) %in% kept_accessions_3]
  flu_validation_list[[3]]$class <- flu_validation_list[[3]]$class[names(flu_validation_list[[3]]$class) %in% kept_accessions_3]
  # Dataset 4
  flu_validation_list[[4]]$pheno <- flu_validation_list[[4]]$pheno[flu_validation_list[[4]]$pheno$`Time.Point` == 1,]
  kept_accessions_4 <- rownames(flu_validation_list[[4]]$pheno)
  flu_validation_list[[4]]$expr <- flu_validation_list[[4]]$expr[,colnames(flu_validation_list[[4]]$expr) %in% kept_accessions_4]
  flu_validation_list[[4]]$class <- flu_validation_list[[4]]$class[names(flu_validation_list[[4]]$class) %in% kept_accessions_4]
  return(flu_validation_list)
}

fix_time_series_for_non_flu_virus <- function(non_flu_virus_list) {
  # Dataset 10 (does it matter that this is also in flu_discovery? Maybe make sure it's in discovery set for non-influenza virus?)
  non_flu_virus_list[[10]]$pheno <- non_flu_virus_list[[10]]$pheno[non_flu_virus_list[[10]]$pheno$`Time.Point` == "Less than 72 hour after symptoms",]
  kept_accessions_10 <- rownames(non_flu_virus_list[[10]]$pheno)
  non_flu_virus_list[[10]]$expr <- non_flu_virus_list[[10]]$expr[,colnames(non_flu_virus_list[[10]]$expr) %in% kept_accessions_10]
  non_flu_virus_list[[10]]$class <- non_flu_virus_list[[10]]$class[names(non_flu_virus_list[[10]]$class) %in% kept_accessions_10]
  # Dataset 11 (does it matter that this is also in flu_discovery? Maybe make sure it's in discovery set for non-influenza virus?)
  non_flu_virus_list[[11]]$pheno <- non_flu_virus_list[[11]]$pheno[non_flu_virus_list[[11]]$pheno$`Time.Point` == "Baseline" | non_flu_virus_list[[11]]$pheno$`Time.Point` == "Day0",]
  non_flu_virus_list[[11]]$pheno <- non_flu_virus_list[[11]]$pheno[!(grepl("Influenza virus", non_flu_virus_list[[11]]$pheno$`Pathogen`)),]
  kept_accessions_11 <- rownames(non_flu_virus_list[[11]]$pheno)
  non_flu_virus_list[[11]]$expr <- non_flu_virus_list[[11]]$expr[,colnames(non_flu_virus_list[[11]]$expr) %in% kept_accessions_11]
  non_flu_virus_list[[11]]$class <- non_flu_virus_list[[11]]$class[names(non_flu_virus_list[[11]]$class) %in% kept_accessions_11]
  return(non_flu_virus_list)
}

# Append final sample size as column to metadata (since dataset sizes can change due to removing time series points, etc.)
append_final_sample_size <- function(current_metadata, current_list) {
  final_sample_sizes <- c()
  for(dataset_id in current_metadata$Study) {
    for(dataset in current_list) {
      if(dataset$formattedName == dataset_id) {
        final_sample_sizes <- c(final_sample_sizes, nrow(dataset$pheno))
      }
    }
  }
  current_metadata$final_sample_size <- final_sample_sizes
  return(current_metadata)
}

# Find overlap with a given set of metadata
append_overlap_with_other_metadata <- function(base_metadata, other_metadata) {
  overlapping <- c()
  for(dataset_id in base_metadata$Study) {
    if(dataset_id %in% other_metadata$Study) {
      overlapping <- c(overlapping, TRUE)
    } else {
      overlapping <- c(overlapping, FALSE)
    }
  }
  base_metadata$overlap <- overlapping
  return(base_metadata)
}

# Find discovery AUCs
get_discovery_aucs <- function(MAGICAL_single_cell_genes, MAGICAL_multiome_14_genes, MAGICAL_multiome_19_genes, flu_discovery_list, non_flu_virus_discovery_list, 
                               bacteria_discovery_list, noninfectious_discovery_list) {
  # Calculate flu (discovery) AUCs for individual genes
  sc_discovery_flu_aucs <- test_individual_genes_on_datasets(MAGICAL_single_cell_genes, flu_discovery_list, "Single_Cell_Paired", "flu_discovery")
  multiome_discovery_paired_flu_aucs <- test_individual_genes_on_datasets(MAGICAL_multiome_14_genes, flu_discovery_list, "Multiome_Paired", "flu_discovery")
  multiome_discovery_all_flu_aucs <- test_individual_genes_on_datasets(MAGICAL_multiome_19_genes, flu_discovery_list, "Multiome_All", "flu_discovery")
  # Calculate non-influenza virus (discovery) AUCs for individual genes
  sc_discovery_non_flu_virus_aucs <- test_individual_genes_on_datasets(MAGICAL_single_cell_genes, non_flu_virus_discovery_list, "Single_Cell_Paired", "non_flu_virus_discovery")
  multiome_discovery_paired_non_flu_virus_aucs <- test_individual_genes_on_datasets(MAGICAL_multiome_14_genes, non_flu_virus_discovery_list, "Multiome_Paired", "non_flu_virus_discovery")
  multiome_discovery_all_non_flu_virus_aucs <- test_individual_genes_on_datasets(MAGICAL_multiome_19_genes, non_flu_virus_discovery_list, "Multiome_All", "non_flu_virus_discovery")
  # Calculate bacteria (discovery) AUCs for individual genes
  sc_discovery_bacteria_aucs <- test_individual_genes_on_datasets(MAGICAL_single_cell_genes, bacteria_discovery_list, "Single_Cell_Paired", "bacteria_discovery")
  multiome_discovery_paired_bacteria_aucs <- test_individual_genes_on_datasets(MAGICAL_multiome_14_genes, bacteria_discovery_list, "Multiome_Paired", "bacteria_discovery")
  multiome_discovery_all_bacteria_aucs <- test_individual_genes_on_datasets(MAGICAL_multiome_19_genes, bacteria_discovery_list, "Multiome_All", "bacteria_discovery")
  # Calculate non-infectious (discovery) AUCs for individual genes
  sc_discovery_noninfectious_aucs <- test_individual_genes_on_datasets(MAGICAL_single_cell_genes, noninfectious_discovery_list, "Single_Cell_Paired", "noninfectious_discovery")
  multiome_discovery_paired_noninfectious_aucs <- test_individual_genes_on_datasets(MAGICAL_multiome_14_genes, noninfectious_discovery_list, "Multiome_Paired", "noninfectious_discovery")
  multiome_discovery_all_noninfectious_aucs <- test_individual_genes_on_datasets(MAGICAL_multiome_19_genes, noninfectious_discovery_list, "Multiome_All", "noninfectious_discovery")
  
  # Combine results from the different gene lists into comprehensive discovery dataframes
  all_discovery_flu_aucs <- rbind(sc_discovery_flu_aucs, multiome_discovery_paired_flu_aucs, multiome_discovery_all_flu_aucs)
  all_discovery_non_flu_virus_aucs <- rbind(sc_discovery_non_flu_virus_aucs, multiome_discovery_paired_non_flu_virus_aucs, multiome_discovery_all_non_flu_virus_aucs)
  all_discovery_bacteria_aucs <- rbind(sc_discovery_bacteria_aucs, multiome_discovery_paired_bacteria_aucs, multiome_discovery_all_bacteria_aucs)
  all_discovery_noninfectious_aucs <- rbind(sc_discovery_noninfectious_aucs, multiome_discovery_paired_noninfectious_aucs, multiome_discovery_all_noninfectious_aucs)
  
  # Create final, overall discovery dataframe
  all_discovery_aucs <- all_discovery_flu_aucs
  all_discovery_aucs$non_flu_virus_discovery_gene_auc <- all_discovery_non_flu_virus_aucs$non_flu_virus_discovery_gene_auc
  all_discovery_aucs$bacteria_discovery_gene_auc <- all_discovery_bacteria_aucs$bacteria_discovery_gene_auc
  all_discovery_aucs$noninfectious_discovery_gene_auc <- all_discovery_noninfectious_aucs$noninfectious_discovery_gene_auc
  all_discovery_aucs <- all_discovery_aucs[,c(1,2,4,5,6,3)]
  return(all_discovery_aucs)
}

capture_flu_pos_and_neg <- function(all_discovery_aucs) {
  all_discovery_aucs$diff <- all_discovery_aucs$flu_discovery_gene_auc - all_discovery_aucs$non_flu_virus_discovery_gene_auc
  flu_pos_genes <- c()
  flu_neg_genes <- c()
  # Plot non-influenza virus AUC vs influenza AUC
  # NOTE - check for overlap between different data types (e.g., noninfectious vs flu)
  # NOTE - one row is removed because we got an NA value - this is probably fine
  # No labels
  auc_plot <- ggplot(all_discovery_aucs, aes(x=flu_discovery_gene_auc, y=non_flu_virus_discovery_gene_auc, group = source)) + geom_point(aes(col = source)) + 
    coord_fixed(xlim = c(0,1), ylim = c(0,1)) + geom_vline(xintercept=0.3, linetype=2) + geom_vline(xintercept=0.7, linetype=2) + 
    geom_abline(slope = 1, intercept = 0.075, linetype=2) + geom_abline(slope = 1, intercept = -0.075, linetype=2) + xlab("Influenza Median AUROC") + ylab("Non-influenza Virus Median AUROC") +
    ggtitle("Influenza Detection") + theme(plot.title = element_text(hjust = 0.5))
  ggsave("plots/individual_auc/non_flu_virus_vs_flu_for_flu_auc.png", plot = auc_plot, device = "png", width = 10, height = 10, units = "in")
  # AUC > 0.7 genes
  current_df <- all_discovery_aucs[all_discovery_aucs$flu_discovery_gene_auc > 0.7,]
  for(current_row_index in 1:nrow(current_df)) {
    current_row <- current_df[current_row_index,]
    line_coord <- current_row$flu_discovery_gene_auc - 0.075
    if(current_row$non_flu_virus_discovery_gene_auc < line_coord) {
      flu_pos_genes <- c(flu_pos_genes, current_row$gene_name)
      print(paste0(current_row$gene_name, " - ", current_row$source))
    }
  }
  flu_pos_df <- current_df[order(current_df$diff),]
  print(flu_pos_df)
  
  # AUC < 0.3 genes
  all_discovery_aucs[all_discovery_aucs$flu_discovery_gene_auc < 0.3,]
  current_df <- all_discovery_aucs[all_discovery_aucs$flu_discovery_gene_auc < 0.3,]
  for(current_row_index in 1:nrow(current_df)) {
    current_row <- current_df[current_row_index,]
    line_coord <- current_row$flu_discovery_gene_auc + 0.075
    if(current_row$non_flu_virus_discovery_gene_auc > line_coord) {
      flu_neg_genes <- c(flu_neg_genes, current_row$gene_name)
      print(paste0(current_row$gene_name, " - ", current_row$source))
    }
  }
  flu_neg_df <- current_df[order(current_df$diff),]
  print(flu_neg_df)
  return(list(flu_pos_df, flu_neg_df, flu_pos_genes, flu_neg_genes))
}

capture_vir_pos_and_neg <- function(all_discovery_aucs) {
  all_discovery_aucs$diff <- all_discovery_aucs$non_flu_virus_discovery_gene_auc - all_discovery_aucs$flu_discovery_gene_auc
  vir_pos_genes <- c()
  vir_neg_genes <- c()
  # Plot non-influenza virus AUC vs influenza AUC for viral specificity
  auc_plot <- ggplot(all_discovery_aucs, aes(x=flu_discovery_gene_auc, y=non_flu_virus_discovery_gene_auc, group = source)) + geom_point(aes(col = source)) + 
    coord_fixed(xlim = c(0,1), ylim = c(0,1)) + geom_hline(yintercept=0.3, linetype=2) + geom_hline(yintercept=0.7, linetype=2) + 
    geom_abline(slope = 1, intercept = 0.075, linetype=2) + geom_abline(slope = 1, intercept = -0.075, linetype=2) + xlab("Influenza Median AUROC") + ylab("Non-influenza Virus Median AUROC") +
    ggtitle("Viral Specificity") + theme(plot.title = element_text(hjust = 0.5))
  ggsave("plots/individual_auc/non_flu_virus_vs_flu_for_non_flu_virus_auc.png", plot = auc_plot, device = "png", width = 10, height = 10, units = "in")
  current_df <- na.omit(all_discovery_aucs[all_discovery_aucs$non_flu_virus_discovery_gene_auc > 0.7,])
  for(current_row_index in 1:nrow(current_df)) {
    current_row <- current_df[current_row_index,]
    line_coord <- current_row$flu_discovery_gene_auc + 0.075
    if(current_row$non_flu_virus_discovery_gene_auc > line_coord) {
      vir_pos_genes <- c(vir_pos_genes, current_row$gene_name)
      print(paste0(current_row$gene_name, " - ", current_row$source))
    }
  }
  vir_pos_df <- current_df[order(current_df$diff),]
  print(vir_pos_df)
  
  # AUC < 0.3 genes
  current_df <- na.omit(all_discovery_aucs[all_discovery_aucs$non_flu_virus_discovery_gene_auc < 0.3,])
  for(current_row_index in 1:nrow(current_df)) {
    current_row <- current_df[current_row_index,]
    line_coord <- current_row$flu_discovery_gene_auc - 0.075
    if(current_row$non_flu_virus_discovery_gene_auc < line_coord) {
      vir_neg_genes <- c(vir_neg_genes, current_row$gene_name)
      print(paste0(current_row$gene_name, " - ", current_row$source))
    }
  }
  vir_neg_df <- current_df[order(current_df$diff),]
  print(vir_neg_df)
  return(list(vir_pos_df, vir_neg_df, vir_pos_genes, vir_neg_genes))
}

capture_bac_pos_and_neg <- function(all_discovery_aucs) {
  all_discovery_aucs$diff <- all_discovery_aucs$bacteria_discovery_gene_auc - all_discovery_aucs$flu_discovery_gene_auc
  bac_pos_genes <- c()
  bac_neg_genes <- c()
  auc_plot <- ggplot(all_discovery_aucs, aes(x=flu_discovery_gene_auc, y=bacteria_discovery_gene_auc , group = source)) + geom_point(aes(col = source)) + 
    coord_fixed(xlim = c(0,1), ylim = c(0,1)) + geom_hline(yintercept=0.3, linetype=2) + geom_hline(yintercept=0.7, linetype=2) + 
    geom_abline(slope = 1, intercept = 0.125, linetype=2) + geom_abline(slope = 1, intercept = -0.125, linetype=2) + xlab("Influenza Median AUROC") + ylab("Bacteria Median AUROC") +
    ggtitle("Bacterial Specificity") + theme(plot.title = element_text(hjust = 0.5))
  ggsave("plots/individual_auc/bacteria_vs_flu_auc.png", plot = auc_plot, device = "png", width = 10, height = 10, units = "in")
  
  # AUC > 0.7 genes
  current_df <- na.omit(all_discovery_aucs[all_discovery_aucs$bacteria_discovery_gene_auc > 0.7,])
  for(current_row_index in 1:nrow(current_df)) {
    current_row <- current_df[current_row_index,]
    line_coord <- current_row$flu_discovery_gene_auc + 0.125
    if(current_row$bacteria_discovery_gene_auc > line_coord) {
      bac_pos_genes <- c(bac_pos_genes, current_row$gene_name)
      print(paste0(current_row$gene_name, " - ", current_row$source))
    }
  }
  bac_pos_df <- current_df[order(current_df$diff),]
  print(bac_pos_df)

  # AUC < 0.3 genes
  current_df <- na.omit(all_discovery_aucs[all_discovery_aucs$bacteria_discovery_gene_auc < 0.3,])
  for(current_row_index in 1:nrow(current_df)) {
    current_row <- current_df[current_row_index,]
    line_coord <- current_row$flu_discovery_gene_auc - 0.125
    if(current_row$bacteria_discovery_gene_auc < line_coord) {
      bac_neg_genes <- c(bac_neg_genes, current_row$gene_name)
      print(paste0(current_row$gene_name, " - ", current_row$source))
    }
  }
  bac_neg_df <- current_df[order(current_df$diff),]
  print(bac_neg_df)
  return(list(bac_pos_df, bac_neg_df, bac_pos_genes, bac_neg_genes))
}

capture_nif_pos_and_neg <- function(all_discovery_aucs) {
  all_discovery_aucs$diff <- all_discovery_aucs$noninfectious_discovery_gene_auc - all_discovery_aucs$flu_discovery_gene_auc
  nif_pos_genes <- c()
  nif_neg_genes <- c()
  auc_plot <- ggplot(all_discovery_aucs, aes(x=flu_discovery_gene_auc, y=noninfectious_discovery_gene_auc, group = source)) + geom_point(aes(col = source)) + 
    coord_fixed(xlim = c(0,1), ylim = c(0,1)) + geom_hline(yintercept=0.3, linetype=2) + geom_hline(yintercept=0.7, linetype=2) + 
    geom_abline(slope = 1, intercept = 0.075, linetype=2) + geom_abline(slope = 1, intercept = -0.075, linetype=2) + xlab("Influenza Median AUROC") + ylab("Non-infectious Median AUROC") +
    ggtitle("Non-infectious Specificity") + theme(plot.title = element_text(hjust = 0.5))
  ggsave("plots/individual_auc/noninfectious_vs_flu_auc.png", plot = auc_plot, device = "png", width = 10, height = 10, units = "in")
  
  # AUC > 0.7 genes
  current_df <- na.omit(all_discovery_aucs[all_discovery_aucs$noninfectious_discovery_gene_auc > 0.7,])
  if(nrow(current_df) > 0) {
    for(current_row_index in 1:nrow(current_df)) {
      current_row <- current_df[current_row_index,]
      line_coord <- current_row$flu_discovery_gene_auc + 0.075
      if(current_row$noninfectious_discovery_gene_auc > line_coord) {
        nif_pos_genes <- c(nif_pos_genes, current_row$gene_name)
        print(paste0(current_row$gene_name, " - ", current_row$source))
      }
    }
    nif_pos_df <- current_df[order(current_df$diff),]
    print(nif_pos_df)
  } else {
    nif_pos_df <- NA
  }

  # AUC < 0.3 genes
  current_df <- na.omit(all_discovery_aucs[all_discovery_aucs$noninfectious_discovery_gene_auc < 0.3,])
  if(nrow(current_df) > 0) {
    for(current_row_index in 1:nrow(current_df)) {
      current_row <- current_df[current_row_index,]
      line_coord <- current_row$flu_discovery_gene_auc - 0.075
      if(current_row$noninfectious_discovery_gene_auc < line_coord) {
        nif_neg_genes <- c(nif_neg_genes, current_row$gene_name)
        print(paste0(current_row$gene_name, " - ", current_row$source))
      }
    }
    nif_neg_df <- current_df[order(current_df$diff),]
    print(nif_neg_df)
  } else {
    nif_neg_df <- NA
  }
  return(list(nif_pos_df, nif_neg_df, nif_pos_genes, nif_neg_genes))
}

create_final_gene_signatures <- function(flu_pos_genes, flu_neg_genes, vir_pos_genes, vir_neg_genes, bac_pos_genes, bac_neg_genes, nif_pos_genes, nif_neg_genes) {
  # Using difference in AUC to create a smaller sig
  #flu_pos_genes <- c("ELF1", "CAPN2", "RAB8B", "ARIH1", "MEF2A")
  #flu_neg_genes <- c("AUTS2", "USP36", "LPCAT1", "SECISBP2", "UQCR11", "ETS1", "HLA-DRA", "PRKCA", "HNRNPDL", "RASSF1", "CD247", "ZNF831", "BRD1")
  #vir_pos_genes <- c("IRAK3")
  #vir_neg_genes <- c("SLC38A1")
  bac_pos_genes <- c("ADGRE5", "TUBA1A", "LRRK2", "SIRPA", "HCAR3")
  bac_neg_genes <- c("CD69", "EZR", "CEBPZ", "SERBP1", "SLC38A1", "PPP3CC", "PRMT1", "NCL")
  
  # My signature
  sig_pos_genes <- c(flu_pos_genes, vir_neg_genes, bac_neg_genes, nif_neg_genes)
  sig_neg_genes <- c(flu_neg_genes, vir_pos_genes, bac_pos_genes, nif_pos_genes)
  
  # Daniel's signature
  old_flu_pos_genes <- c("CAPN2", "CITED2", "DUSP1", "ELF1", "IER2", "KLF6", "MAP3K8", "MOB1A", "RAB8B", "SAMHD1", "USP3", "USP8")
  old_flu_neg_genes <- c("HNRNPDL", "KRT10", "OAZ2", "OST4", "PRMT1", "TOMM6", "UQCR11")
  old_vir_pos_genes <- c("S100A4", "UBE2J1")
  old_vir_neg_genes <- c("JUN")
  old_bac_pos_genes <- c("APLP2", "HCAR3", "IL1RAP", "IQSEC1", "LRRK2", "VCAN")
  old_bac_neg_genes <- c("CD69", "CEBPZ", "SERBP1")
  old_nif_pos_genes <- c()
  old_nif_neg_genes <- c("CHURC1")
  
  old_sig_pos_genes <- c(old_flu_pos_genes, old_vir_neg_genes, old_bac_neg_genes, old_nif_neg_genes)
  old_sig_neg_genes <- c(old_flu_neg_genes, old_vir_pos_genes, old_bac_pos_genes, old_nif_pos_genes)
  
  # Signature 1 (28619954)
  other_sig_1_pos_genes <- c("IFI27")
  other_sig_1_neg_genes <- c()
  
  # Signature 2 (26682989)
  other_sig_2_pos_genes <- c("CD38", "HERC5", "HERC6", "IFI6", "IFIH1", "LGALS3BP", "LY6E", "MX1", "PARP12", "RTP4", "ZBP1")
  other_sig_2_neg_genes <- c()
  return(list(sig_pos_genes, sig_neg_genes, old_sig_pos_genes, old_sig_neg_genes, other_sig_1_pos_genes, other_sig_1_neg_genes,
              other_sig_2_pos_genes, other_sig_2_neg_genes))
}

test_gene_signatures <- function(sig_pos_genes, sig_neg_genes, old_sig_pos_genes, old_sig_neg_genes, other_sig_1_pos_genes, other_sig_1_neg_genes, 
                                 other_sig_2_pos_genes, other_sig_2_neg_genes) {
  # Flu
  plot_pooled_auc(sig_pos_genes, sig_neg_genes, flu_validation_list, "Influenza Pooled AUC (New)")
  plot_pooled_auc(old_sig_pos_genes, old_sig_neg_genes, flu_validation_list, "Influenza Pooled AUC (Old)")
  plot_pooled_auc(other_sig_1_pos_genes, other_sig_1_neg_genes, flu_validation_list, "Influenza Pooled AUC (Other Sig 1)")
  plot_pooled_auc(other_sig_2_pos_genes, other_sig_2_neg_genes, flu_validation_list, "Influenza Pooled AUC (Other Sig 2)")
  # Non-influenza virus
  plot_pooled_auc(sig_pos_genes, sig_neg_genes, non_flu_virus_validation_list, "Non-Influenza Virus Pooled AUC (New)")
  plot_pooled_auc(old_sig_pos_genes, old_sig_neg_genes, non_flu_virus_validation_list, "Non-Influenza Virus Pooled AUC (Old)")
  plot_pooled_auc(other_sig_1_pos_genes, other_sig_1_neg_genes, non_flu_virus_validation_list, "Non-Influenza Virus Pooled AUC (Other Sig 1)")
  plot_pooled_auc(other_sig_2_pos_genes, other_sig_2_neg_genes, non_flu_virus_validation_list, "Non-Influenza Virus Pooled AUC (Other Sig 2)")
  # Bacteria
  plot_pooled_auc(sig_pos_genes, sig_neg_genes, bacteria_validation_list, "Bacteria Pooled AUC (New)")
  plot_pooled_auc(old_sig_pos_genes, old_sig_neg_genes, bacteria_validation_list, "Bacteria Pooled AUC (Old)")
  plot_pooled_auc(other_sig_1_pos_genes, other_sig_1_neg_genes, bacteria_validation_list, "Bacteria Pooled AUC (Other Sig 1)")
  plot_pooled_auc(other_sig_2_pos_genes, other_sig_2_neg_genes, bacteria_validation_list, "Bacteria Pooled AUC (Other Sig 2)")
  # Noninfectious
  plot_pooled_auc(sig_pos_genes, sig_neg_genes, noninfectious_validation_list, "Noninfectious Pooled AUC (New)")
  plot_pooled_auc(old_sig_pos_genes, old_sig_neg_genes, noninfectious_validation_list, "Noninfectious Pooled AUC (Old)")
  plot_pooled_auc(other_sig_1_pos_genes, other_sig_1_neg_genes, noninfectious_validation_list, "Noninfectious Pooled AUC (Other Sig 1)")
  plot_pooled_auc(other_sig_2_pos_genes, other_sig_2_neg_genes, noninfectious_validation_list, "Noninfectious Pooled AUC (Other Sig 2)")
}