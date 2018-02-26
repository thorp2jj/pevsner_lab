#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 5

for bam in /mnt/data/jeremy/projects/sws/collab_samples/cleaned_bams/*.bam
do
remove_extension=${bam%.*}
echo $remove_extension
newname=${remove_extension##*/}
echo $newname
samtools index $bam $newname
done
