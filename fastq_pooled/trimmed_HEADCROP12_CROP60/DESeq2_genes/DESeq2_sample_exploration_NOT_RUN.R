#!/applications/R/R-3.4.0/bin/Rscript

### Perform differential expression analysis using read counts generated by Salmon

# R version 3.4.0
# DESeq2 version 1.16.1
# Note that this R version or later is required for DESeq2 version 1.16 or later,
# in which "the log2 fold change shrinkage is no longer default for the DESeq and nbinomWaldTest functions".
# DESeq2 version 1.16 introduces "a separate function lfcShrink, which performs log2 fold change shrinkage
# for visualization and ranking of genes." (see https://support.bioconductor.org/p/95695/ and
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#changes)

library(DESeq2)
print(packageVersion("DESeq2"))
#[1] ‘1.16.1’

inDir <- "/home/ajt200/analysis/180928_Pallas_RNAseq_series1/fastq_pooled/trimmed_HEADCROP12_CROP60/DESeq2_genes/"
plotDir <- paste0(inDir, "sample_exploration_plots/")
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

load(paste0(inDir, "tximport.RData"))

sampleTable <- data.frame(sample = c("wt RNA-seq Rep1",
                                     "hta6 RNA-seq Rep1",
                                     "hta7 RNA-seq Rep1",
                                     "cmt3 RNA-seq Rep1",
                                     "hta6 hta7 RNA-seq Rep1",
                                     "cmt3 hta7 RNA-seq Rep1",
                                     "cmt3 hta6 hta7 RNA-seq Rep1"),
                          condition = factor(rep.int(c("wt", "hta6", "hta7", "cmt3",
                                                       "hta6 hta7", "cmt3 hta7", "cmt3 hta6 hta7"),
                                                     times = c(1, 1, 1, 1, 1, 1, 1))))
rownames(sampleTable) <- colnames(txi$counts)
print(sampleTable)
#                             sample      condition
#sample1             wt RNA-seq Rep1             wt
#sample2           hta6 RNA-seq Rep1           hta6
#sample3           hta7 RNA-seq Rep1           hta7
#sample4           cmt3 RNA-seq Rep1           cmt3
#sample5      hta6 hta7 RNA-seq Rep1      hta6 hta7
#sample6      cmt3 hta7 RNA-seq Rep1      cmt3 hta7
#sample7 cmt3 hta6 hta7 RNA-seq Rep1 cmt3_hta6_hta7

dds <- DESeqDataSetFromTximport(txi = txi,
                                colData = sampleTable,
                                design = ~condition)

# Pre-filter the dataset
print(nrow(dds))
#[1] 27586 
# Retain only rows that have more than a single count across all samples
dds <- dds[rowSums(counts(dds)) > 1,]
print(nrow(dds))
#[1] 27437

## The rlog and variance stabilizing transformations
# see http://www.bioconductor.org/help/workflows/rnaseqGene/#the-rlog-and-variance-stabilizing-transformations

rld <- rlog(dds, blind = TRUE)
print(head(assay(rld), 3))

vsd <- vst(dds, blind = TRUE)
print(head(assay(vsd), 3))

# Visualise the effect of transformation
library(dplyr)
library(ggplot2)
library(hexbin)

# For the log2 approach, estimate size factors to account for sequencing depth
# Sequencing-depth correction is done automatically for rlog and vst
dds <- estimateSizeFactors(dds)

df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized = T)[,1:2]+1)) %>%
    mutate(transformation = "log2(normalized counts + 1)"),
  as_data_frame(assay(rld)[,1:2]) %>% mutate(transformation = "rlog"),
  as_data_frame(assay(vsd)[,1:2]) %>% mutate(transformation = "vst"))

colnames(df)[1:2] <- c("Sample 1 transformed counts",
                       "Sample 2 transformed counts")

plot_transformed_counts <- ggplot(df,
                                  aes(x = `Sample 1 transformed counts`,
                                      y = `Sample 2 transformed counts`)) +
                           geom_hex(bins = 80) +
                           coord_fixed() +
                           facet_grid(. ~ transformation) +
                           labs(fill = "Occurrences") +
                           theme_classic()
                           #theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
ggsave(plot_transformed_counts,
       file = paste0(plotDir,
                     "Sample1_vs_Sample2_transformed_counts_log2countsPlus1_rlog_vst_genes.pdf"))


## Sample distances

# Sample distances using the rlog-transformed counts
sampleDists <- dist(t(assay(rld)))
print(sampleDists)

library(pheatmap)
library(RColorBrewer)

# Heatmap of sample distances using the rlog-transformed counts
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- as.vector(sampleTable$sample)
colnames(sampleDistMatrix) <- NULL
mycols <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pdf(paste0(plotDir, "sample_distances_heatmap_rlog_genes.pdf"),
    height = 5, width = 7.5, onefile = F)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = mycols)
dev.off()

# Sample distances using the Poisson Distance
library(PoiClaClu)
poisd <- PoissonDistance(t(counts(dds)))
print(poisd)
#Value of alpha used to transform data:  0.5757143
#This type of normalization was performed: mle
#Dissimilarity computed for  7  observations.

# Heatmap
samplePoisDistMatrix <- as.matrix(poisd$dd)
rownames(samplePoisDistMatrix) <- as.vector(sampleTable$sample)
colnames(samplePoisDistMatrix) <- NULL
pdf(paste0(plotDir, "sample_distances_heatmap_Poisson_genes.pdf"), height = 5, width = 7.5, onefile = F)
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = mycols)
dev.off()

# PCA plots for visualising sample-to-sample distances
PCAplot_rlog <- DESeq2::plotPCA(rld, intgroup = c("condition", "sample")) +
                  theme(plot.margin = grid::unit(c(0,0,0,0), "mm")) +
                  theme_classic() +
                  theme(panel.border = element_rect(colour = "black",
                                                    fill = NA, 
                                                    size = 1)) +
                  coord_fixed()
ggsave(PCAplot_rlog,
       file = paste0(plotDir, "PCAplot_rlog_genes.pdf"),
       width = 20, height = 20, units = "cm")

PCAplot_rlog_data <- DESeq2::plotPCA(rld, intgroup = c("condition", "sample"),
                                     returnData = T)
print(PCAplot_rlog_data)
# Obtain percentage variance explained by PC1 and PC2 for plotting using ggplot2
percentVar <- round(100 * attr(PCAplot_rlog_data, "percentVar"))

PCAggplot_rlog <- ggplot(PCAplot_rlog_data,
                         aes(x = PC1, y = PC2,
                             colour = sample)) +
                  geom_point(size = 3, shape = c(0, 1, 2, 3, 7, 8, 11)) +
                  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
                  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
                  theme(plot.margin = grid::unit(c(0,0,0,0), "mm")) +
                  theme_classic() +
                  theme(panel.border = element_rect(colour = "black",
                                                    fill = NA,
                                                    size = 1)) +
                  coord_fixed()
ggsave(PCAggplot_rlog, 
       file = paste0(plotDir, "PCAggplot_rlog_genes.pdf"),
       height = 20, width = 20, units = "cm")


# MDS (multi-dimensional scaling) plots for visualising sample-to-sample distances
# rlog-transformed counts
mds <- as.data.frame(colData(rld)) %>%
         cbind(cmdscale(sampleDistMatrix))
MDSplot_rlog <- ggplot(mds, aes(x = `1`, y = `2`,
                                #shape = condition,
                                colour = sample)) +
                  geom_point(size = 3, shape = c(0, 1, 2, 3, 7, 8, 11)) +
                  theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
                  theme_classic() +
                  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
                  coord_fixed()
ggsave(MDSplot_rlog,
       file = paste0(plotDir, "MDSplot_rlog_genes.pdf"), height = 20, width = 20, units = "cm")

# Poisson distance
mdsPois <- as.data.frame(colData(dds)) %>%
             cbind(cmdscale(samplePoisDistMatrix))
MDSplot_Pois <- ggplot(mdsPois, aes(x = `1`, y = `2`,
                                #shape = condition,
                                colour = sample)) +
                  geom_point(size = 3, shape = c(0, 1, 2, 3, 7, 8, 11)) +
                  theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
                  theme_classic() +
                  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
                  coord_fixed()
ggsave(MDSplot_Pois,
       file = paste0(plotDir, "MDSplot_Pois_genes.pdf"), height = 20, width = 20, units = "cm")

# gene clustering
# Clustering of 20 genes with greatest rlog-transformed variance across samples
library(genefilter)
topVargenes <- head(order(rowVars(assay(rld)), decreasing = T), 20)
mat <- assay(rld)[topVargenes, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[, c("sample", "condition")])
pdf(paste0(plotDir, "gene_clustering_rld_topVar20_genes.pdf"), height = 5, width = 7.5, onefile = F)
pheatmap(mat, annotation_col = anno)
dev.off()
