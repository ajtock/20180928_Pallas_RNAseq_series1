#!/bin/bash

source activate RNAseq_mapping
snakemake -p --cores 48
conda deactivate
