#!/usr/bin/env bash
#SBATCH --job-name mosaic_hunter
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 7-00:05

for file in /mnt/data/jeremy/projects/sws/collab_samples/cleaned_bams/*treat_recal_cleaned.bam; do
echo $file
basename=${file##*/}
echo $basename
outname=${basename%%_recal_cleaned.bam}
echo $outname
cntrl_name=${outname%%treat}Ctrl_recal_cleaned.bam
echo $cntrl_name
mkdir /mnt/data/jeremy/projects/sws/collab_samples/mosaic_hunter/output/${outname}
    for individual in $(ls /mnt/data/jeremy/projects/sws/collab_samples/gatk_output/Data_Files/VCF/individual_vcfs/${outname}.vcf.gz); do
           sex=$(bcftools +vcf2sex $individual | cut -f2)
           #echo $individual
           #echo $sex
    done
    
sbatch -N 1 -n 1 -c 11 /mnt/data/jeremy/scripts/sbatch_scripts/general.sh "java -jar /mnt/data/jeremy/programs/MosaicHunter/build/mosaichunter.jar genome \
    -P input_file=$file \
    -P reference_file=/mnt/data/reference/hs37d5.fa \
    -P mosaic_filter.sex=$sex \
    -P mosaic_filter.dbsnp_file=/mnt/data/jeremy/programs/MosaicHunter/resources/dbsnp_137.b37.tsv \
    -P repetitive_region_filter.bed_file=/mnt/data/jeremy/programs/MosaicHunter/resources/all_repeats.b37.bed \
    -P indel_region_filter.bed_file=/mnt/data/jeremy/projects/sws/collab_samples/error_prone_regions_beds/beds/${outname}.indel.cnvs.bed \
    -P common_site_filter.bed_file=/mnt/data/jeremy/programs/MosaicHunter/resources/WES_Agilent_50M.error_prone.b37.bed \
    -P output_dir=/mnt/data/jeremy/projects/sws/collab_samples/mosaic_hunter/output/${outname} \
    -P mode=paired_naive \
    -P control_bam_file=/mnt/data/jeremy/projects/sws/collab_samples/cleaned_bams/${cntrl_name}"

sleep 10
done







