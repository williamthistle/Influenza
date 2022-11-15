library(DESeq2)
library(limma)
library(sva)
library(glmnet)

funcByF=function (data, groups, func=mean, ...){
  meanstats = function(x) {
    tapply(t(x), groups, func, ...)
  }
  t(apply(data, 1, meanstats))
}


#map is a dataframe with a single column and rownames corresponding to rownames of data
collapseNamesSum=function (data, map) {
  data=as.matrix(data)
  cm=intersect(rownames(data), rownames(map))
  map=map[cm,,drop=F]
  unames = unique(map[,1])
  
  dataout = matrix(nrow = length(unames), ncol = ncol(data))
  
  rownames(dataout) = unames
  colnames(dataout) = colnames(data)
  
  for (i in 1:length(unames)) {
    ii = rownames(map)[(which(map[,1] == unames[i]))]
    dataout[i, ] = apply(data[ii, , drop=F],2,sum)
  }
  return(dataout)
}


#useful limma wrapper
fitLimma=function(data, grp, contr, extra=NULL, pathway=NULL, method="ls", plotHists=F, trend=F, subset=NULL){
  if(is.null(dim(grp))){
    mod=model.matrix(~0+grp)
    colnames(mod)=levels(grp)
  }
  else{
    mod=grp
  }
  if(!is.null(extra)){
    mod=cbind(mod,model.matrix(~0+extra))
  }
  colnames(mod)=make.names(colnames(mod))
  #show(colnames(mod))
  
  if(!is.null(subset)){
    mod=mod[subset,]
    data=data[, subset]
  }
  ne=nonEstimable(mod)
  if(!is.null(ne)){
    iirm=match(ne, colnames(mod))
    mod=mod[,-iirm]
  }
  fit=lmFit(data, mod, method=method)
  contrastMat=makeContrasts(contrasts=contr, levels=mod)
  fitC=contrasts.fit(fit, contrasts=contrastMat)
  fitC=eBayes(fitC, trend=trend, robust=T)
  fitC=addQV(fitC)
  if(plotHists){
    par(mfrow=c(2,3))
    for(i in 1:ncol(fitC$p.value)){
      hist(fitC$p.value[,i], main=contr[i])
    }
  }
  if(! is.null(pathway)){
    warning("running pathway analysis")
    pathwayout=list()
    np=names(pathway)
    for( i in 1:length(pathway)){
      warning(paste0("pathway ",np[i]))
      
      respath=wilcoxGMT(fitC$t[,1], pathway[[i]], simple=T)
      ii=which(respath[,4]<0.2)
      if(length(ii)>0){
        oo=order(respath[ii,2])
        pathwayout[[np[i]]]=respath[ii[oo],]
      }
      else{
        warning(paste0("min QV=",min(respath[,4], na.rm=T), " for pathway ", np[i]))
      } 
    }
    
    fitC$pathway=pathwayout
  }
  fitC$fit=fit
  fitC$mod=mod
  #show(colSums(fitC$q.value<0.2))
  fitC 
}


addQV=function(fit2){
  fit2$q.value=fit2$p.value
  for ( i in 1:ncol(fit2$p.value)){
    ii=which(!is.na(fit2$p.value[,i]))
    fit2$q.value[ii,i]=QV(fit2$p.value[ii,i])
  }
  return(fit2)
}

QV=function(pval){
  
  x=try(qvalue(pval))
  
  if(!is.list(x)){
    warning("Q-value error")
    # hist(pval)
    return(p.adjust(pval, method="BH"))
  }
  else{
    return(x$qvalue)
  }
}


commonRows=function(data1, data2){
  intersect(rownames(data1), rownames(data2))
}


myresid=function(dat, lab, useMean=T){
  if (is.null(dim(lab))){
    mod=model.matrix(~1+lab);
  }
  else{
    mod=lab
  }
  ne <- nonEstimable(mod)
  if (!is.null(ne)){ 
    cat("Coefficients not estimable:", paste(ne, collapse = " "), 
        "\n")
    mod=mod[, -match(ne, colnames(mod))]
  }
  
  n=dim(dat)[2]
  Id=diag(n)
  resid=dat %*% (Id - mod %*% solve(t(mod) %*% mod) %*% 
                   t(mod))
  colnames(resid)=colnames(dat)
  if (useMean){
    resid=sweep(resid,1,apply(dat,1,mean), "+")  
  }
  
  resid
}


mycoef=function(dat, lab, useMean=T){
  if (is.null(dim(lab))){
    mod=model.matrix(~1+lab);
  }
  else{
    mod=lab
  }
  ne <- nonEstimable(mod)
  if (!is.null(ne)){ 
    cat("Coefficients not estimable:", paste(ne, collapse = " "), 
        "\n")
    mod=mod[, -match(ne, colnames(mod))]
  }
  
  n=dim(dat)[2]
  Id=diag(n)
  coef=dat %*%  mod %*% solve(t(mod) %*% mod)
  coef
}


getBestIndex=function(gres, se=F, min.beta.size=2){
  if(!se){
    ii=which(gres$lambda==gres$lambda.min)
  }
  else{
    ii=which(gres$lambda==gres$lambda.1se)
  }
  beta.size=colSums(as.matrix(gres$glmnet.fit$beta)!=0)
  if(beta.size[ii]<min.beta.size){
    message("min.beta.size not reached")
    while(beta.size[ii]<min.beta.size){
      ii=ii+1
    }
  }
  ii
  
}


