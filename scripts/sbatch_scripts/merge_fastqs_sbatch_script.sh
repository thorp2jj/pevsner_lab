#!/bin/bash
#
#SBATCH --job-name merge_fastqs
#SBATCH --ntasks 1
#SBATCH --time 06:59:00
#SBATCH --cpus-per-task 10
#SBATCH --nodes 1
#SBATCH --nodelist plab-node04
#

cat 80kbp_mean_S1_L001_R2_001.fastq.gz 80kbp_mean_S1_L002_R2_001.fastq.gz > 80kbpcat_S1_L001_R2_001.fastq.gz

