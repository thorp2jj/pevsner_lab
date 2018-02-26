#!/usr/bin/env bash
#SBATCH -c 1
#SBATCH -J parse_vcf
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 7-00:05


module add bcftools/1.3

for file in *.vcf*; do
      for sample in `bcftools query -l $file`; do
              bcftools view -c1 -Oz -s $sample -o ${sample}.vcf.gz $file
                done
            done
