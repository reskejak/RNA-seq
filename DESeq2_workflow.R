# Jake Reske
# Michigan State University, 2020
# reskejak@msu.edu

# Consider an RNA-seq experimental design of n=3 samples for a treatment condition ("Treat") and control condition ("Control").
# This is a quick-start workflow for DESeq2 differential gene expression analysis, beginning with a raw counts table/matrix of m genes by n samples.
###############################################
###############################################

# load dependencies
library("ggplot2")
library("IHW")
library("DESeq2")

###############################################
# set working directory
setwd("~/my_experiment")

# import raw counts data
raw.counts <- read.csv("raw_counts.csv")
# convert to numeric matrix
ids <- raw.counts$ID
raw.counts <- raw.counts[, -1] # remove gene ID column
raw.counts <- as.matrix(raw.counts)
rownames(raw.counts) <- ids
rm(ids)

###############################################

# prepare samples table containing identifiers
samples <- data.frame(c(colnames(raw.counts)))
# rename samples column
colnames(samples)[1] <- "sample"
# add sample condition
samples$condition <- factor(c("Control", "Control", "Control", "Treat", "Treat", "Treat"))
# add row names
rownames(samples) <- samples$sample

# relevel to set Control as reference, or another level if desired
samples$condition <- relevel(samples$condition, "Control")

###############################
# DESeq2
dds <- DESeqDataSetFromMatrix(countData = raw.counts,
                              colData = samples,
                              design = ~ condition)

# Pre-filter low count genes (minimum set to average of 1 count per sample)
keep <- rowSums(counts(dds)) >= ncol(counts(dds))
dds <- dds[keep, ]
rm(keep)
# e.g. 23,282 expressed genes

##################
# calculate normalized counts (median of ratios method)
dds <- estimateSizeFactors(dds)
normalized.counts <- counts(dds, normalized=TRUE)

# write normalized counts output to csv
# write.csv(normalized.counts, file="normalized_counts.csv")

#################
# regularized-logarithm transformation (rlog) for downstream PCA/MDS etc.
# note: blind=FALSE elicits slightly-supervised transformation
# "differences between [experimental design variables] will not contribute to expected variance-mean trend of data"
# blind=TRUE elicits fully unsupervised transformation
rld <- rlog(dds, blind = FALSE)
rlog.counts <- assay(rld)

# write rlog counts output to csv
# write.csv(rlog.counts, file="rlog_counts.csv")

##########################################

# plot PCA based on rlog-transformed counts
pcaData.rlog <- plotPCA(rld, ntop=500, returnData=TRUE)
percentVar.rlog <- round(100 *attr(pcaData.rlog, "percentVar"))

ggplot(pcaData.rlog, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar.rlog[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar.rlog[2], "% variance")) +
  coord_fixed() +
  theme_bw() +
  theme(text = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

###########################################
# DGE
dds <- DESeq(dds)

# for results etc., determine coef as output from resultsNames(dds)
# coef=1: "Intercept"
# coef=2: "condition_Treat_vs_Control"

res <- results(dds, name="condition_Treat_vs_Control") # DGE results
resLFC <- lfcShrink(dds, coef=2) # log2 fold change shrinkage
res.IHW <- results(dds, name="condition_Treat_vs_Control", filterFun=ihw) # independent hypothesis weighting for p-value correction

################
# combine analyses into one data.frame
resFinal.df <- data.frame(res)
resFinal.df$shrink.log2FoldChange <- resLFC$log2FoldChange
resFinal.df$shrink.lfcSE <- resLFC$lfcSE
resFinal.df$padj.IHW <- res.IHW$padj
resFinal.df$weight <- res.IHW$weight
# reorder columns
resFinal.df <- resFinal.df[c("baseMean",
                             "log2FoldChange",
                             "lfcSE",
                             "shrink.log2FoldChange",
                             "shrink.lfcSE",
                             "stat",
                             "pvalue",
                             "padj",
                             "padj.IHW",
                             "weight")]

# order genes by FDR (IHW-adjusted)
resFinal.df <- resFinal.df[order(resFinal.df$padj.IHW), ]

# write full output to csv
# write.csv(resFinal.df, file="DESeq2_Treat_vs_Control_DGE_results.csv")

