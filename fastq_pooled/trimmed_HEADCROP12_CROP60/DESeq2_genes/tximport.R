#!/applications/R/R-3.4.0/bin/Rscript
### Convert gene-level counts into format for use with DESeq2

# R version 3.4.0

library(tximport)
print(packageVersion("tximport"))
#[1] ‘1.4.0’

inDir <- "/home/ajt200/analysis/180928_Pallas_RNAseq_series1/fastq_pooled/trimmed_HEADCROP12_CROP60"
outDir <- "/home/ajt200/analysis/180928_Pallas_RNAseq_series1/fastq_pooled/trimmed_HEADCROP12_CROP60/DESeq2_genes/"

# Read in table of sample IDs that will be used to specify paths to count files
samples <- read.table(file.path(inDir, "/samples_DESeq2_genes.txt"), header = T)
print(samples)
#        genotype           directory                                  sample
#1             WT salmon_quants_genes             WT_RNAseq_Pallas_Rep1_quant
#2           hta6 salmon_quants_genes           hta6_RNAseq_Pallas_Rep1_quant
#3           hta7 salmon_quants_genes           hta7_RNAseq_Pallas_Rep1_quant
#4           cmt3 salmon_quants_genes           cmt3_RNAseq_Pallas_Rep1_quant
#5      hta6_hta7 salmon_quants_genes      hta6_hta7_RNAseq_Pallas_Rep1_quant
#6      cmt3_hta7 salmon_quants_genes      cmt3_hta7_RNAseq_Pallas_Rep1_quant
#7 cmt3_hta6_hta7 salmon_quants_genes cmt3_hta6_hta7_RNAseq_Pallas_Rep1_quant

# Specify paths to count files
files <- file.path(inDir, samples$genotype, samples$directory, samples$sample, "quant.sf")
# Set 1:#_samples
names(files) <- c("wt",
                  "hta6",
                  "hta7",
                  "cmt3",
                  "hta6_hta7",
                  "cmt3_hta7",
                  "cmt3_hta6_hta7")
all(file.exists(files))

# Create a dataframe of transcript IDs and corresponding gene IDs
transID <- read.table(files[1], colClasses = c(NA, rep("NULL", 4)), header = T)
tx2gene <- data.frame(cbind(as.vector(transID[,1]), substr(transID[,1], 1, 9)))
colnames(tx2gene) <- c("TXNAME", "GENEID")

# Import transcript-level counts, summarised at gene level
# (reads that map to transcript IDs with a common parent gene ID are pooled)
library(readr)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
print(names(txi))
#[1] "abundance"           "counts"              "length"             
#[4] "countsFromAbundance"

# Import transcript-level counts, summarised at transcript level with "txOut = TRUE"
txi.tx <- tximport(files, type = "salmon", txOut = TRUE, tx2gene = tx2gene)
# Then summarise to gene level
txi.sum <- summarizeToGene(txi.tx, tx2gene)
# These two approaches should produce identical results
print(all.equal(txi$counts, txi.sum$counts))
#[1] TRUE

print(head(txi$counts))
#            wt hta6      hta7      cmt3 hta6_hta7 cmt3_hta7 cmt3_hta6_hta7
#AT1G01010  105   53  131.0000  103.0000   83.0000    85.000       141.0000
#AT1G01020  143  152  103.2074  147.0000  141.0000   147.000       181.0000
#AT1G01030  197  194  231.0000  251.0000  191.0000   201.000       234.0000
#AT1G01040 1429 1394 1527.0000 2217.0000 1950.0000  2201.000      1931.0000
#AT1G01050  244  207  177.2382  139.0000  191.0000   156.000       114.0000
#AT1G01060  224   67  225.0000  321.3197  197.0741   264.384       194.3265

save(txi,
     file = paste0(outDir, "tximport.RData"))

