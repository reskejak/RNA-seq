# Downstream gene-level analysis of RNA-seq data
# April 2018
# Michigan State University
# reskejak@msu.edu

setwd("/mnt/usr/foo")

# experiment identifier for file labeling
prefix <- "foo_RNA"

#################################
# define raw count matrices for input
raw.counts.matrix <- "foo_RNA_raw_counts.csv"
#################################

# for converting Ensembl IDs to symbols
# use "hsapiens_gene_ensembl" for humans; use "mmusculus_gene_ensembl" for mouse
library("biomaRt")
human.ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

# import samples table containing identifiers i.e. fOO1
samples <- data.frame(c("foo_treat1", "foo_treat2", "foo_treatn", "foo_control1", "foo_control2", "foo_controln"))
# rename samples column
colnames(samples)[1] <- "sample"
# add sample condition (treat and control)
samples$condition <- factor(c(rep("treat", 3), 
                              rep("control", 3)))
# add row names
rownames(samples) <- samples$sample
# relevel to set control as reference
samples$condition <- relevel(samples$condition, "control")

# import raw counts table and convert to numeric matrix
raw.counts.df <- read.csv(raw.counts.matrix)
ids <- raw.counts.df$X
raw.counts <- raw.counts.df[, -1]
raw.counts <- as.matrix(raw.counts)
rownames(raw.counts) <- ids
colnames(raw.counts) <- samples$sample

# remove suffix from gene identifiers
rownames(raw.counts) <- sub("\\..*","", rownames(raw.counts))

rm(raw.counts.df)
rm(ids)

#################
# convert normalized counts table from Ensembl IDs to symbols
converted.raw.counts <- raw.counts
G_list <- getBM(filters= "ensembl_gene_id", 
                attributes= c("ensembl_gene_id","external_gene_name"),
                values=rownames(converted.raw.counts),
                mart=human.ensembl)
converted_unique_ids <- G_list[!duplicated(G_list$external_gene_name), ]

# subset raw counts for unique symbols
converted.raw.counts <- subset(converted.raw.counts, 
                                      rownames(converted.raw.counts) %in% converted_unique_ids$ensembl_gene_id,
                                      select = c(0:ncol(converted.raw.counts)))
rownames(converted.raw.counts) <- converted_unique_ids$external_gene_name

# order gene symbols alphabetically for final table
converted.raw.counts <- converted.raw.counts[order(rownames(converted.raw.counts)), ]

# write converted raw counts table to CSV
write.csv(converted.raw.counts, 
          file=paste(prefix, "DESeq2_converted_raw-counts.csv", sep="_"))

################
# DESeq2
library("digest")
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = converted.raw.counts,
                              colData = samples,
                              design = ~ condition)
                              
# Pre-filter low count genes (minimum set to average of 1 count per sample)
keep <- rowSums(counts(dds)) >= 6
dds <- dds[keep, ]
rm(keep)

##################
# calculate normalized counts (median of ratios method)
dds <- estimateSizeFactors(dds)
normalized.counts <- counts(dds, normalized=TRUE)
# order gene identifiers numerically
normalized.counts <- normalized.counts[order(rownames(normalized.counts)), ]

# write normalized counts table to CSV
write.csv(normalized.counts, 
          file=paste(prefix, "DESeq2_filtered_normalized-counts.csv", sep="_"))


#################
# regularized-logarithm transformation (rlog) for downstream PCA/MDS etc.
# note: blind=FALSE elicits slightly-supervised transformation
# "differences between [experimental design variables] will not contribute to expected variance-mean trend of data"
# blind=TRUE elicits fully unsupervised transformation
rld <- rlog(dds, blind = FALSE)
rlog.counts <- assay(rld)
# order gene identifiers numerically
rlog.counts <- rlog.counts[order(rownames(rlog.counts)), ]

#############
# determine variance distribution and plot
rld.df <- data.frame(assay(rld))
rld.var <- apply(rld.df, 1, var)
rld.df$var <- rld.var

var.plot <- ggplot(rld.df, aes(var)) + geom_density(aes(y=..scaled..)) + xlim(0, 2)

# export as PDF
pdf(paste(prefix, "DESeq2_rlog-transformation_variance_density-plot.pdf", sep="_"))
var.plot
dev.off()

# subset rld.df by variance for unsupervised clustering
# note: adjust var threshold to call around 500-2500 top variable genes
rld.df.subset <- subset(rld.df, 
                        var > 0.05, 
                        select=c(0:(ncol(rld.df)-1)))
rm(rld.df)
rm(rld.var)

#############     
# write rld$condition levels to treat and control for plotting purposes
rld$condition <- factor(c(rep("treat", 3),
                          rep("control", 3)))

# plot PCA based on rlog-transformed counts
pcaData.rlog <- plotPCA(rld, returnData=TRUE)
percentVar.rlog <- round(100 *attr(pcaData.rlog, "percentVar"))

pca.plot <- ggplot(pcaData.rlog, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar.rlog[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar.rlog[2], "% variance")) +
  coord_fixed() +
  theme_bw() +
  theme(text = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
		
# export as PDF
pdf(paste(prefix, "DESeq2_rlog-transformation_PCA.pdf", sep="_"))
pca.plot
dev.off()

rm(pcaData.rlog)
rm(percentVar.rlog)

#############
# DIFFERENTIAL GENE EXPRESSION
dds <- DESeq(dds)
res <- results(dds)
res

# log2 fold change shrinkage
# determine coef as output from resultsNames(dds)
resLFC <- lfcShrink(dds, coef=2)

# independent hypothesis weighting for p-value correction
library("IHW")
res.IHW <- results(dds, filterFun=ihw)

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

# order genes by IHW-adjusted p-value
resFinal.df.p.ordered <- resFinal.df[order(resFinal.df$padj.IHW), ]

# write full output to csv
write.csv(resFinal.df.p.ordered,
          file=paste(prefix, "DESeq2_treat-vs-control_results.csv", sep="_"))

##############
# filter by padj.IHW < 0.05 and corrected-|LFC| >= 1
# adjust these thresholds to match your own needs

resFinal.df.sig.p.ordered <- subset(resFinal.df.p.ordered, 
                                    (padj.IHW < 0.05) & (abs(shrink.log2FoldChange) >= 1), 
                                    select = c(0:ncol(resFinal.df.p.ordered)))

# write signature to csv
write.csv(resFinal.df.sig.p.ordered, 
          file=paste(prefix, "DESeq2_treat-vs-control_signature-genes.csv", sep="_"))

######################
# subset rlog.counts by signature genes for downstream plotting, etc.
rlog.counts.sig <- subset(rlog.counts, 
                          rownames(rlog.counts) %in% rownames(resFinal.df.sig.p.ordered), 
                          select=c(0:ncol(rlog.counts)))

# convert ensembl IDs to symbols for rlog.counts.sig
converted.rlog.counts.sig <- rlog.counts.sig

# convert ensembl IDs to symbols
rownames(converted.rlog.counts.sig) <- G_list_sig$external_gene_name

# order gene symbols alphabetically
converted.rlog.counts.sig <- converted.rlog.counts.sig[order(rownames(converted.rlog.counts.sig)), ]

#########
# Generate Heatmaps
#########
library("circlize")
library("ComplexHeatmap")

# UNSUPERVISED CLUSTERING
var.threshold = 0.05
# rld.df.subset = [var > threshold] genes
# adjust labeling to match selected variance threshold
centered.scaled.rld.df.subset <- t(scale(t(rld.df.subset)))
Unsup <- Heatmap(centered.scaled.rld.df.subset, 
                 clustering_distance_rows = "euclidean",
                 width= unit(8, "cm"),
                 column_title = paste("Unsupervised clustering (var > ", var.threshold, "; ", nrow(rld.df.subset), " genes)", sep=""),
                 name="scaled rlog(counts)",
                 show_row_names=FALSE)

pdf(paste(prefix, "DESeq2_rlog-counts", var.threshold, "Unsupervised_clustering.pdf", sep="_"))
Unsup
dev.off()

# TREAT VS CONTROL GENE SIGNATURE
# center and scale rlog counts from treat vs. control gene signature
centered.scaled.converted.rlog.counts.sig <- t(scale(t(converted.rlog.counts.sig)))
sig.ht1 <- Heatmap(centered.scaled.converted.rlog.counts.sig, 
                   clustering_distance_rows = "euclidean",
                   show_row_names=FALSE,
                   show_column_names=FALSE,
                   km = 2,
                   km_title = c("cluster %i"),
                   column_title = paste("Transcriptome signature (p < 0.05; |log2FC| > 1; ", nrow(centered.scaled.converted.rlog.counts.sig), " genes)", sep=""),
                   name="scaled rlog(counts)",
                   width = unit(8, "cm"))
                          
pdf(paste(prefix, "DESeq2_rlog-counts_treat-vs-control_signature-genes_clustering.pdf", sep="_"))
sig.ht1
dev.off()

#############
# label genes of interest with annotation WIP
genes.of.interest<- c("GAPDH", "ACTB", "TP53")
                      
positions <- which(rownames(centered.scaled.converted.rlog.counts.sig) %in% genes.of.interest, arr.ind = TRUE)
labels <- rownames(centered.scaled.converted.rlog.counts.sig[positions,])

sig.subset <- subset(centered.scaled.converted.rlog.counts.sig, 
                     rownames(centered.scaled.converted.rlog.counts.sig) %in% genes.of.interest, 
                     select=c(0:ncol(centered.scaled.converted.rlog.counts.sig)))
sig.labels <- rowAnnotation(link = row_anno_link(at = c(positions), 
                                                 labels = labels),
                            width = unit(2, "cm"))

pdf(paste(prefix, "DESeq2_rlog-counts_treat-vs-control_signature-genes_clustering_with-labels.pdf", sep="_"))
sig.ht1 + sig.labels
dev.off()
