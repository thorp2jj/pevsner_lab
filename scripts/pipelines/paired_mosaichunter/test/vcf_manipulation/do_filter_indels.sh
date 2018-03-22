#!/usr/bin/env bash
#SBATCH -c 1
#SBATCH -J parse_vcf
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 7-00:05
#SBATCH -w plab-node04


for file in /mnt/data/jeremy/projects/sws/collab_samples/gatk_output/Data_Files/VCF/indels_only_vcfs/*.vcf.gz; 
do
    echo $file
    basename=${file##*/}
    echo $basename
    vcfname=${basename%%.gz}
    outname=${basename%%.recode.vcf.gz}
    echo $outname
    gunzip $file
    
    python3 /mnt/data/jeremy/projects/sws/collab_samples/gatk_output/Data_Files/VCF/indels_only_vcfs/filter_indels.py \
    /mnt/data/jeremy/projects/sws/collab_samples/gatk_output/Data_Files/VCF/indels_only_vcfs/$vcfname > /mnt/data/jeremy/projects/sws/collab_samples/gatk_output/Data_Files/VCF/indels_only_vcfs/$outname.indels.5bpflanks.bed
    
    bgzip /mnt/data/jeremy/projects/sws/collab_samples/gatk_output/Data_Files/VCF/indels_only_vcfs/$vcfname
done
