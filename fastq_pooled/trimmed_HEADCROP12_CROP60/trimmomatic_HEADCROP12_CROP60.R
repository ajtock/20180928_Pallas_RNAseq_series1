#!/applications/R/R-3.4.0/bin/Rscript

# Remove Illumina adapter sequences, PCR primer fragments and low-quality bases using
# Trimmomatic version 0.36 (Bolger et al., 2014) in paired-end and palindrome modes,
# which enable reliable detection of adapter read-through.
# Sequence quality metrics for raw and trimmed reads were evaluated using FastQC version 0.11.4 (Andrews, 2010).

# Usage:
# csmit -m 20G -c 8 "./trimmomatic_HEADCROP12_CROP60.R Pallas_RNAseq_Col_series1"

args <- commandArgs(trailingOnly = T)
libName <- args[1]

progDir <- "/home/ajt200/tools/Trimmomatic-0.36/"
fastqDir <- "/home/ajt200/analysis/180928_Pallas_RNAseq_series1/fastq_pooled/"
outDir <- paste0(fastqDir, "trimmed_HEADCROP12_CROP60/")
adapterDir <- "/home/ajt200/tools/Trimmomatic-0.36/adapters/"

# Paired-end fastq files with "R1" and "R2" naming convention must exactly match "*R1_001.fastq" and "*R2_001.fastq"
# for correct handling by trimmomatic -basein option
#fastq <- system(paste0("ls ", fastqDir, "*_R1_001.fastq.gz | xargs -n 1 basename"), intern = T)
fastq <- system(paste0("ls ", fastqDir, libName, "_R1_001.fastq.gz"), intern = T)
fastqOut <- sub("R1", "trimmed", fastq)

system(paste0("java -jar ", progDir, "trimmomatic-0.36.jar PE -threads 8 -phred33 -basein ",
              fastq,
              " -baseout ",
              fastqOut,
              " ILLUMINACLIP:",
              adapterDir, "cat_all_and_TruSeq_Single_Indexes.fa:2:30:10:1:true",
              " HEADCROP:12 CROP:60 SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:36"))
print(paste0(fastq, " trimming complete"))
system(paste0("mv ", fastqDir, "*trimmed_001*.fastq.gz ", outDir))

