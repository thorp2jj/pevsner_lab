#!/user/bin/env bash

for bam in `cat /mnt/data/jeremy/projects/smri/Data_Files/BAM/100X/100x_samplelist.txt`
do
    remove_suffix=${bam%_recal.bam}
    sample_name=${remove_suffix##*/}
    #echo ${remove_prefix}
    for fastq in /mnt/data/DNASeq/BPD_BSMN/stanley_collection/split_trimmed_fastqs/*.fastq.gz
    do
        name=${fastq##*/}
        remove_s=${name%%_*}
        #echo ${name}
        #echo ${remove_s}
        if [[ ${sample_name} == ${remove_s} ]];
        then
            ln -s ${fastq} .
        else
            continue
        fi
    done
done
