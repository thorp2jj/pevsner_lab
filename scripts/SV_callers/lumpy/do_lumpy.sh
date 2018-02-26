#!/usr/bin/env bash

for i in /mnt/data/jeremy/projects/kk_patients/700_series/analysis/Data_Files/BAM/*recal.bam
do
out1=${i%_recal.bam}
#echo $out1
out2=${out1##*/}
#echo $out2
out3=$out2.lumpy.vcf
echo $out3
sbatch -N 1 -n 1 -c 31 /mnt/data/jeremy/scripts/sbatch_scripts/general.sh "lumpyexpress -v -B $i -o $out3"
sleep 3
done 
