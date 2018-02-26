#!/bin/bash

#loop through bam to fastq conversion bash script

for bam in $(ls /mnt/data/bpd/90_trios/*.bam | head -n 290 | tail -n +271) 
do
	echo ${bam}
	basename=${bam##*/}
	echo $basename
	outname=${basename%%.bam}
	echo $outname
    sbatch -N 1 -n 1 -c 1 --mem 6144 /mnt/data/bpd/90_trios/fastqs/scripts/bam_to_fastq_2.sh $bam /mnt/data/bpd/90_trios/fastqs/${outname}
	sleep 10
done

# |tail -n +86 | head -n 1
