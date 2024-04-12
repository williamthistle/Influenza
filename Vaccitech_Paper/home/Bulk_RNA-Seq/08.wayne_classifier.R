library(glmnet)

getBestIndex=function(gres, se=F){
  if(!se){
    which(gres$lambda==gres$lambda.min)
  }
  else{
    which(gres$lambda==gres$lambda.1se)
  }
}#getBestIndex

#parameter
########################################################
#input X and y
#X: sample-by-feature

#nfold <- 5
#se <- F
#alp <- 0.9
#n_bootstrap <- 100
#test_frac <- 0.1
#bootstrap_report <- T


cv_report <- function(X, y, nfold = 5, se = F, alp = 0.9, n_bootstrap = 100, test_frac = 0.1, bootstrap_report = T, task = "regression", contrast = c("Control", "LatePost")){
  #generate test, train
  ########################################################
  
  test.include <- list()
  train.include <- list()
  
  if (bootstrap_report){
    seeds <- 1:n_bootstrap
    for (i in 1:length(seeds)){
      set.seed(seeds[i])
      test.include[[i]] <- sample(1:length(y), length(y)*test_frac)
      train.include[[i]] <- setdiff(1:length(y), test.include[[i]])
    }#for i
  }else{
    #split by fold, everything appears once
    test.include <- split(1:length(y), 1: (1/test_frac))
    for (i in 1:lenggth(test.include)){
      train.include[[i]] <- setdiff(1:length(y), test.include[[i]])
    }#for i
  }#else
  
  
  #generate fold for inner loop
  set.seed(1)
  temp <- sample(length(train.include[[1]]))
  foldid.temp <- rep(1:nfold, each = floor(length(temp)/nfold)+1)[1:length(temp)]
  foldid <- rep(0, length(temp))
  for (i in 1:length(foldid)){
    foldid[temp[i]] <- foldid.temp[i]
  }#for i
  
  
  #bootstrap feature selection or fold wise
  ########################################################
  
  if (task!="classification_multi"){
    y.pred.rec <- c()    
    y.test.rec <- c()
  }else{
    y.pred.rec <- array(NA, dim = c(length(y), length(contrast), length(test.include)))
    y.test.rec <- array(NA, dim = c(length(y), length(contrast), length(test.include)))  
  }#else
  
  feature <- list()
  lambda.collection <- list()
  gres.list <- list()
  
  
  if (task == "regression"){
    print("regrgession")
  }else if (task == "classification"){
    print("classification")
    
    #convert y
    #############################
    for (i in 1:length(contrast)){
      y[y==contrast[i]] <- as.character(i)  
    }#for i
    
  }else if (task == "classification_multi"){
    print("classification mutli")
    
    #convert y
    #############################
    for (i in 1:length(contrast)){
      y[y==contrast[i]] <- as.character(i)  
    }#for i
    
  }#else if 
  
  
  for (ii in 1:length(test.include)){
    print(ii)
    
    train.index <- train.include[[ii]]
    test.index <- test.include[[ii]]
    X.train <- X[train.index,]
    X.test <- X[test.index,]
    y.train <- y[train.index]
    y.test <- y[test.index]
    
    
    if (task == "regression"){
      cvgres <- cv.glmnet(x=X.train, y=y.train, foldid = foldid, alpha=alp,family="gaussian", type.measure = "mse", intercept = T)
      gres <- glmnet(x=X.train, y=y.train,  family="gaussian", lambda = cvgres$lambda[getBestIndex(cvgres,se)], alpha=alp, intercept = T)  
      preds=predict(gres, newx=X.test, type="response")
      
      
    }else if (task == "classification"){
      cvgres <- cv.glmnet(x = X.train, y =y.train, foldid = foldid, alpha = alp, family = "binomial", type.measure = "auc", intercept = T)
      gres <- glmnet(x=X.train, y=y.train,  family="binomial", lambda = cvgres$lambda[getBestIndex(cvgres,se)], alpha=alp, intercept = T)  
      preds=predict(gres, newx=X.test, type="response")
      
      
    }else if (task == "classification_multi"){
      cvgres <- cv.glmnet(x=X.train, y=y.train, foldid = foldid, alpha=alp, family="multinomial", type.measure = "class", intercept = T)
      gres <- glmnet(x=X.train, y=y.train, family="multinomial", lambda = cvgres$lambda[getBestIndex(cvgres,se)], alpha=alp, intercept = T)
      preds=predict(gres, newx=X.test, type="response")
      
    }#else if 
    
    lambda.collection[[ii]] <- cvgres$lambda[getBestIndex(cvgres,se)]
    
    if (task!="classification_multi"){
      feature[[ii]] <- rownames(gres$beta)[which(gres$beta[,1]!=0)]
      gres.list[[ii]] <- gres
      
      temp <- rep(NA, length(y))
      temp[test.index] <- preds
      y.pred.rec <- cbind(y.pred.rec, temp)
      
      temp <- rep(NA, length(y))
      temp[test.index] <- y.test
      y.test.rec <- cbind(y.test.rec, temp)
      
    }else{
      tem <- c()
      for (jj in 1:length(gres$beta)){
        tem <- c(tem, rownames(gres$beta[[jj]])[which(gres$beta[[jj]][,1]!=0)])
      }#for jj
      feature[[ii]] <- unique(tem)
      gres.list[[ii]] <- gres
      
      
      temp <- matrix(NA, nrow = length(y), ncol = length(contrast))
      temp[test.index, ] <- preds
      y.pred.rec[,,ii] <- temp
      
      temp <- matrix(NA, nrow = length(y), ncol = length(contrast))
      for (jj in 1:length(test.index)){
        temp[test.index[jj], as.numeric(y.test)[jj]] <- 1
      }#for jj
      y.test.rec[,,ii] <- temp
    }#else
    
    
  }#for ii
  
  
  #check features
  ########################################################
  #sort(table(unlist(feature)),decreasing = T)[1:10]
  sig.report <- table(unlist(feature))
  sig <- names(sig.report)[sig.report>5]
  
  
  if (task == "regression"){
    y.pred <- c()
    y.sd <- c()
    #report the average and standard deviation
    for (ii in 1:nrow(y.pred.rec)){
      temp <- y.pred.rec[ii,which(!is.na(y.pred.rec[ii,]))]
      y.pred <- c(y.pred, mean(temp))
      y.sd <- c(y.sd, sd(temp))
    }#for ii
    
    return(list(y.pred = y.pred, y.sd = y.sd, y.pred.rec = y.pred.rec, sig = sig, sig.report = sig.report, feature = feature, lambda.collection = lambda.collection, y = y, y.test.rec = y.test.rec, gres.list = gres.list))
    
  }else if (task == "classification"){
    
    #average of the probability
    y.pred <- c()
    y.sd <- c()
    #report the average and standard deviation
    for (ii in 1:nrow(y.pred.rec)){
      temp <- y.pred.rec[ii,which(!is.na(y.pred.rec[ii,]))]
      y.pred <- c(y.pred, mean(temp))
      y.sd <- c(y.sd, sd(temp))
    }#for ii
    
    return(list(y.pred = y.pred, y.sd = y.sd, y.pred.rec = y.pred.rec, sig = sig, sig.report = sig.report, feature = feature, lambda.collection = lambda.collection, y = y, y.test.rec = y.test.rec, gres.list = gres.list))
    
    
  }else if (task == "classification_multi"){
    
    #also average of the probability, given the class that samples are supposed to be
    y.pred <- matrix(0, nrow = length(y), ncol = length(contrast))
    y.sd <- matrix(0, nrow = length(y), ncol = length(contrast))
    
    for (ii in 1:length(y)){
      for (jj in 1:length(contrast)){
        temp <- y.pred.rec[ii,jj, which(!is.na(y.pred.rec[ii,jj,]))]
        y.pred[ii,jj] <- mean(temp)
        y.sd[ii,jj] <- sd(temp)
      }#for jj
    }#for ii
    
    return(list(y.pred = y.pred, y.sd = y.sd, y.pred.rec = y.pred.rec, sig = sig, sig.report = sig.report, feature = feature, lambda.collection = lambda.collection, y = y, y.test.rec = y.test.rec, gres.list = gres.list))
    
  }#else if 
  
}#cv_report