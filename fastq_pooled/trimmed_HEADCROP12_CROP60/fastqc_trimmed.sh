#!/bin/bash

# Run fastqc on each trimmed fastq.gz file in directory
# and move output files to fastqc directory

[ -d fastqc ] || mkdir fastqc
[ -d fastqc/trimmed ] || mkdir fastqc/trimmed
for f in *P.fastq.gz
do
( echo "Processing $f"
  fastqc $f
  echo "$f processing complete" ) &
done
wait
mv *_fastqc.* fastqc/trimmed

