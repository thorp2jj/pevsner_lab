#!/usr/bin/env bash
#SBATCH -c 1
#SBATCH -J parse_vcf
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 7-00:05


module add bcftools/1.3
int_samp="TAC-94"
for file in /mnt/data/jeremy/projects/smri/Data_Files/VCF/combined_filt.vcf; do
    for sample in `bcftools query -l $file`; do
        if [ "$sample" == "$int_samp" ]; then
              bcftools view -c1 -Ov -s $sample -o ${sample}.vcf $file
        else
            continue
        fi      
    done
done
