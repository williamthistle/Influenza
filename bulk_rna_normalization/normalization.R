library(DESeq2)
library(ggplot2)
library(biomaRt)
library(limma)
library(HGNChelper)
library(data.table)

tscale=function(mat){
  t(scale(t(mat)))
}#tscale


setwd("C:/Users/wat2/Documents/GitHub/Influenza/bulk_rna_normalization")
source("normalizeDataHelper.R")
genemap <- readRDS("genemap.RDS")
gene <- readRDS("gene.RDS")
setwd("C:/Users/wat2/Documents/local_data_files/")



counts <- fread("rsem_genes_count.txt", header = T, sep = ",")
row.names <- as.character(counts$gene_id)
counts <- counts[,2:ncol(counts)]
counts <- as.matrix(counts)
rownames(counts) <- row.names



#convert to homo symbol
##############################################
iim=match(rownames(counts), rownames(gene))
notmatched=which(is.na(iim))
if(length(notmatched)>0){
  counts=counts[-notmatched,]
}#if

geneIDs=gene[iim, "name", drop=F]
rownames(geneIDs)=rownames(gene)[iim]

countsSymbol=collapseNamesSum(counts, geneIDs)



#rename the rownames of countsSymbol
gene.diff <- setdiff(rownames(countsSymbol), genemap$hgnc_symbol)

before_mapping <- c()
after_mapping <- c()
outlier <- c()
for (i in 1:length(gene.diff)){
  temp <- which(hgnc.table$Symbol==gene.diff[i])
  if (length(temp) > 1){
    print(gene.diff[i])
    if (gene.diff[i]=="SARS"){
      before_mapping <- c(before_mapping, "SARS")
      after_mapping <- c(after_mapping, "SARS1")
    }else if (gene.diff[i]=="SEPT2"){
      before_mapping <- c(before_mapping, "SEPT2")
      after_mapping <- c(after_mapping, "SEPTIN2")
    }else if (gene.diff[i]=="QARS"){
      before_mapping <- c(before_mapping, "QARS")
      after_mapping <- c(after_mapping, "QARS1")
    }else if (gene.diff[i]=="DEC1"){
      before_mapping <- c(before_mapping, "DEC1")
      after_mapping <- c(after_mapping, "DELEC1")
    }else if (gene.diff[i]=="LOR"){
      before_mapping <- c(before_mapping, "LOR")
      after_mapping <- c(after_mapping,"LORICRIN")
    }else if (gene.diff[i]=="H3F3AP1"){
      before_mapping <- c(before_mapping, "H3F3AP1")
      after_mapping <- c(after_mapping, "H3P38")
    }#else if 
  }else{
    if (length(temp)==1){
      before_mapping <- c(before_mapping, gene.diff[i])
      after_mapping <- c(after_mapping, hgnc.table$Approved.Symbol[temp])
    }else{
      outlier <- c(outlier, gene.diff[i])
    }#else
  }#else
}#for i

row.names <- rownames(countsSymbol)
row.names.index <- match(before_mapping, row.names)
#all(row.names[row.names.index]==before_mapping)
row.names[row.names.index] <- after_mapping




collapse_hgnc <- function(mat, gene.name){
  gene.names.unique <- names(table(gene.name))
  res <- matrix(0, nrow = length(gene.names.unique), ncol= ncol(mat))
  
  for (i in 1:length(gene.names.unique)){
    temp <- which(gene.name==gene.names.unique[i])
    if (length(temp) > 1){
      res[i,] <- apply(mat[temp,],2,sum)
    }else{
      res[i,] <- mat[temp,]
    }#else
  }#for i
  
  rownames(res) <- gene.names.unique
  colnames(res) <- colnames(mat)
  return(res)
}#collapse_hgnc


countsSymbol <- collapse_hgnc(countsSymbol, row.names)
# Note - need to add extra tab to beginning of header (TODO: fix this)
write.table(countsSymbol, "rsem_genes_count.processed.txt", quote = FALSE, sep = "\t")



#further normalization
#######################################################################################
#make a  biotype lookup for downstream filtering
iim=match(rownames(countsSymbol), genemap$hgnc_symbol)
#make a biotype label by symbol id
mybiotype=genemap$gene_biotype[iim]
names(mybiotype)=genemap$hgnc_symbol[iim]



#RPKM
#http://metagenomics.wiki/pdf/definition/rpkm-calculation
#######################################################################################
gene.length <- genemap$end_position-genemap$start_position
gene.length <- gene.length[iim]
names(gene.length) <- genemap$hgnc_symbol[iim]


#convert to RPKM
total.reads <- apply(countsSymbol,2,sum)

countsSymbol.RPKM <- countsSymbol
countsSymbol.RPKM <- sweep(countsSymbol.RPKM, 2, total.reads/1e6, "/")
countsSymbol.RPKM <- sweep(countsSymbol.RPKM, 1, gene.length/1e3, "/")



#look at the mean to pick a reasonable non-expressed cutoff
#######################################################################################
mm=apply(log2(1+countsSymbol),1,mean)
#hist(mm)
#it should be a bimodal distribution

#only >5 (expression levels) and protein_coding genes are included. You can adjust this accordingly
goodgenes=which(mm>5&mybiotype[rownames(countsSymbol)]=="protein_coding")
countsSymbolUse=countsSymbol.RPKM[goodgenes,]


sf=estimateSizeFactorsForMatrix(countsSymbolUse) #just a wrapper around DESeq2 normalization (see helper.R)
countsSymbolN=sweep(countsSymbolUse,2,sf, "/")
#compute an offset for log, smallest positive value
adj=min(countsSymbolN[countsSymbolN>0])

data=log2(adj+countsSymbolN)
boxplot(data, outline=F) #visually check the distribution for each sample and see if there are outliers



#quantile normalization
###########################################################################################
data <- normalizeBetweenArrays(data, method = "quantile")



#get the gc and length correlations
#we will get the average symbol value
###########################################################################################
gc_content <- genemap$percentage_gene_gc_content
gcval=t(funcByF(rbind(gc_content), genemap$hgnc_symbol, func = mean))

gene_length <- genemap$transcript_length
len=t(funcByF(rbind(gene_length), genemap$hgnc_symbol, func = mean))
#check
#all(rownames(gcval)==rownames(len))
gclen=cbind(gcval, len)


cm=commonRows(data , gclen)
gclen.cor=cor(data[cm,], gclen[cm,])




#normalize for gc and length
###########################################################################################
techMod1 <- model.matrix(~0+gclen.cor[,1]+gclen.cor[,2])

data.rpkm <- countsSymbolN
data1 <- myresid(data, techMod1, useMean = F)
data1.T <- myresid(data, techMod1, useMean = T)


#data1.T is usually used for the next step

#data.rpkm: RPKM
#data: RPKM+log+quantile normalization
#data1: data -> data + correct for gene length+gc content
#data1.T: data1 + add the mean value back
data1.T.scaled <- tscale(data1.T)
write.table(data1.T.scaled, "rsem_genes_count.processed.normalized.txt", quote = FALSE, sep = "\t")

save(data.rpkm , data1, data1.T, pData, file = paste0("normalization.RData"))




