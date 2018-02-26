#!/bin/bash

#Pull filenames for snakemake pipeline

#Terminal command: sh sample_names_sm.sh | awk 'NR % 2 == 0' - > test.txt

for gz in `ls /work-zfs/jpevsne1/BP_Wgs_Data/Data_Files/snakemake_test_bipolar_subset/*.fastq.gz`
do
    #echo ${gz}
    basename=${gz##*/}
    #echo $basename
    outname=${basename%%??.fastq.gz}
    echo $outname


#    sbatch -o /work-zfs/jpevsne1/BP_Wgs_Data/Logs/bam_to_fastq_%j.out \
#        -e /work-zfs/jpevsne1/BP_Wgs_Data/Logs/bam_to_fastq_%j.err \
#        -t 6-23 \
#        --mem 3072 \
#        -N 1 \
#        -n 1 \
#        -c 1 \
#        /home-2/jthorpe6@jhu.edu/scratch/bash_scripts/bam_fastq_conversion.sh \
#        $bam /work-zfs/jpevsne1/BP_Wgs_Data/Data_Files/Fastq/${outname}
#    sleep 10
done

