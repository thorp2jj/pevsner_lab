#!/bin/bash
#
#SBATCH --job-name common_experiment_bcl_to_fastq
#SBATCH --output /mnt/data/10X/common_experiment_bcl/output.txt
#SBATCH --error /mnt/data/10X/common_experiment_bcl/error.txt
#SBATCH --ntasks 1
#SBATCH --time 06:59:00
#SBATCH --cpus-per-task 32
#

longranger demux --run=/mnt/data/10X/1611UNHX-0003/L006 --localcore=32 --localmem=128

#longranger demux --run=/mnt/data/10X/1611UNHX-0003/L006 --localcores=32 --localmem=128
