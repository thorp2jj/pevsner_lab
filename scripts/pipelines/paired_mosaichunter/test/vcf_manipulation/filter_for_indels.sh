#!/usr/bin/env bash
#SBATCH -c 1
#SBATCH -J parse_vcf
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 7-00:05
#SBATCH -w plab-node04


for file in /mnt/data/jeremy/scripts/pipelines/mosaichunter/test/vcf_manipulation/*.vcf;
do
    echo $file
    basename=${file##*/combined_filt.}
    echo $basename
    vcfname=${basename%%.gz}
    outname=${basename%%.vcf.gz}
    echo $outname
    vcftools --gzvcf $file --out /mnt/data/jeremy/scripts/pipelines/mosaichunter/test/vcf_manipulation/$outname.indels --keep-only-indels --recode
    #bgzip /mnt/data/jeremy/scripts/pipelines/mosaichunter/test/vcf_manipulation/$outname.indels.recode.vcf
done
