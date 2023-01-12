library(DESeq2)

set.seed(1)

#### SIMULATION ####
# Simulation to test whether much different dispersion parameter in random condition can create additional significant p-values in pairwise comparison (fake DEGs)
# First, compare A and B when they come from the same distributions - we should find 0 DEGs (or maybe 1 from noise)
first <- rnbinom(10, mu = 4, size = 1)
second <- rnbinom(10, mu = 4, size = 1)
third <- rnbinom(10, mu = 4, size = 1)
fourth <- rnbinom(10, mu = 4, size = 1)
simulation_counts <- data.frame(first, second, third, fourth)
condition <- c("A", "A", "B", "B")
simulation_metadata <- data.frame(condition)
rownames(simulation_metadata) <- c("first", "second", "third", "fourth")
simulation_analysis <- DESeqDataSetFromMatrix(countData = simulation_counts,
                                              colData = simulation_metadata,
                                              design = ~ condition)
simulation_analysis <- DESeq(simulation_analysis)
simulation_analysis_results <- results(simulation_analysis, alpha = 0.05)
simulation_analysis_results <- simulation_analysis_results[order(simulation_analysis_results$padj),]
simulation_analysis_results <- subset(simulation_analysis_results, padj < 0.05)
simulation_analysis_results
# We find 1 DEG (false positive)
# Next, we add C and it's very different (different mean and, more importantly, much different dispersion parameter)
# Since DESeq2 estimates an overall dispersion parameter for each gene, this could affect the accuracy of B vs A
fifth <- rnbinom(10, mu = 400, size = 10)
sixth <- rnbinom(10, mu = 400, size = 10)
simulation_counts <- cbind(simulation_counts, fifth)
simulation_counts <- cbind(simulation_counts, sixth)
simulation_metadata[nrow(simulation_metadata) + 1,] = c("C")
simulation_metadata[nrow(simulation_metadata) + 1,] = c("C")
rownames(simulation_metadata) <- c("first", "second", "third", "fourth", "fifth", "sixth")
simulation_analysis <- DESeqDataSetFromMatrix(countData = simulation_counts,
                                              colData = simulation_metadata,
                                              design = ~ condition)
simulation_analysis <- DESeq(simulation_analysis)
simulation_analysis_results <- results(simulation_analysis, contrast = c("condition", "B", "A"), alpha = 0.05)
simulation_analysis_results <- simulation_analysis_results[order(simulation_analysis_results$padj),]
simulation_analysis_results <- subset(simulation_analysis_results, padj < 0.05)
simulation_analysis_results
# We now have 4 DEGs for B vs A. So this answers our question - by having a third group with a different dispersion parameter,
# we found fake DEGs for B vs A.
# This means we should do pairwise comparisons for any groups of interest (alternatively, we could figure out which groups have very different variance and exclude those from our overall model
# when doing pairwise comparisons between all other conditions, but that sounds like unnecessary work)