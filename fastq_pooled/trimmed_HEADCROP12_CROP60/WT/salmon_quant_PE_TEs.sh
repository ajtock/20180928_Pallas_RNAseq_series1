#!/bin/bash
# Quantify gene expression levels using salmon 0.9.1

# Example usage via condor submission system on hydrogen node7
# csmit -m 1G -c 24 "./salmon_quant_PE_TEs.sh WT_RNAseq_Pallas"

prefix=$1

for i in $prefix"_Rep"{1..1};
do
  samp=`basename ${i}`
  echo "Processing sample ${samp}"
  salmon quant --index /projects/ajt200/TAIR10/salmon_TAIR10_Buisine_TEs_strandAware_index \
               --libType A \
               --mates1 ${i}_R1.fastq.gz \
               --mates2 ${i}_R2.fastq.gz \
               --threads 24 \
               --output ./salmon_quants_TEs_strandAware/${samp}_quant
done
