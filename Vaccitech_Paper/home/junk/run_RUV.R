library("RUVSeq")

# RUV
current_analysis_run <- DESeq(dds)
res <- results(current_analysis_run, contrast = c("time_point", test_time, baseline_time), alpha = 0.05)
res <- res[order(res$padj),]
res <- subset(res, padj < 0.05)
sc_pseudobulk_genes <- unique(sc_pseudobulk_deg_table$Gene_Name)
length(intersect(rownames(res), sc_pseudobulk_genes))

res <- results(current_analysis_run, contrast = c("time_point", test_time, baseline_time), alpha = 0.05)
set <- newSeqExpressionSet(counts(dds))
idx  <- rowSums(counts(set) > 5) >= 2
set  <- set[idx, ]
set <- betweenLaneNormalization(set, which="upper")
not.sig <- rownames(res)[which(res$pvalue > .95)]
empirical <- rownames(set)[ rownames(set) %in% not.sig ]
set <- RUVg(set, empirical, k=2)

metadata_subset$W1 <- set$W_1
metadata_subset$W2 <- set$W_2
dds <- DESeqDataSetFromMatrix(countData = counts_subset, colData = metadata_subset, design = ~ subject_id + time_point + W1 + W2)

# Redo PCA
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = c( "subject_id", "time_point"), pcsToUse = c(1,2), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = subject_id, shape = time_point)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")

# Rerun analysis
current_analysis_run <- DESeq(dds)
res <- results(current_analysis_run, contrast = c("time_point", test_time, baseline_time), alpha = 0.05)
res <- res[order(res$padj),]
res <- subset(res, padj < 0.05)
