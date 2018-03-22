#!/usr/bin/env bash
#SBATCH -c 1
#SBATCH -J parse_vcf
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 7-00:05



for file in /mnt/data/jeremy/projects/smri/Data_Files/VCF/combined_filt.vcf; do
  for sample in `bcftools query -l $file`; do
      if [ "$sample" = "TAC-94" ]; then
        bcftools view -c1 -Oz -s $sample -o ${file/.vcf*/.$sample.vcf.gz} $file
      fi
    done
done
