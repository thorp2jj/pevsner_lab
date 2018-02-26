#!/bin/bash

#bam -> sorted by read group ID -> fastq



#samtools sort -n /work-zfs/jpevsne1/wgs_fg/bams_rg/98222_1031.rg.bam | bedtools bamtofastq -i /dev/stdin -fq1 98222_1031.rg.1.fastq -fq2 98222_1031.rg.2.fastq

#gzip 98222_1031.rg.1.fastq
#gzip 98222_1031.rg.2.fastq


samtools sort -m 1G -O bam -T $2 -n $1 | bedtools bamtofastq -i /dev/stdin -fq ${2}_1.fastq -fq2 ${2}_2.fastq

gzip ${2}_1.fastq
gzip ${2}_2.fastq
