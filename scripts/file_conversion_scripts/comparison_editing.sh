#!/usr/bin/env bash

#Editing BAM files

for bam in `ls /work-zfs/jpevsne1/wgs_fg/bams_rg/*bam`
do
    basename=${bam##*/}
#    echo $basename
    outname=${basename%%.rg.bam}
    echo $outname
done > test1.txt

for fastq in `ls /work-zfs/jpevsne1/BP_Wgs_Data/Data_Files/Fastq/*_1.fastq.gz`
do
    fbase=${fastq##*/}
#    echo $fbase
    fout1=${fbase%%_[1].fastq.gz}
    echo $fout1
done > test2.txt

diff <(cat test1.txt | sort) <(cat test2.txt | sort)

#diff <( echo "$outname" ) <( echo "$fout1" )
