#!/bin/bash

mkdir -p microbiome

# extract microbiome data from raw reads
for i in {1..3};do
	d=S${i}
	zgrep -A 3 "microbiome" fastq_raw/${d}_1.fastq.gz | grep -v "^--$" > microbiome/${d}_1.fq
	zgrep -A 3 "microbiome" fastq_raw/${d}_2.fastq.gz | grep -v "^--$" > microbiome/${d}_2.fq
done