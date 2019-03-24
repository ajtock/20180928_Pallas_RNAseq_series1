#!/bin/bash

# Run fastqc on each fastq.gz file in directory
# and move output files to fastqc directory

[ -d fastqc ] || mkdir fastqc
[ -d fastqc/raw ] || mkdir fastqc/raw
for f in *.fastq.gz
do
( echo "Processing $f"
  fastqc $f
  echo "$f processing complete" ) &
done
wait
mv *_fastqc.* fastqc/raw

