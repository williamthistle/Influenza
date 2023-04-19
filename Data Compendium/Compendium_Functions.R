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

test_individual_genes_on_datasets <- function(gene_list, data_list, source, disease_tag) {
  gene_aucs <- c()
  for(gene in gene_list) {
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
    # Test gene set signature on influenza samples and print median AUROC
    flu_aucs <- sapply(X = data_list, FUN = calculateAUROC, signature = sig)
    gene_aucs <- c(gene_aucs, median(flu_aucs, na.rm = TRUE))
  }
  gene_auc_name <- paste0(disease_tag, "_gene_auc")
  final_df <- data.frame("gene_name" = gene_list, temp_name = gene_aucs, "source" = rep(source, length(gene_list)))
  names(final_df)[names(final_df) == "temp_name"] <- gene_auc_name
  return(final_df)
}

geom_mean <- function(x) {
  return(exp(mean(log(x))))
}