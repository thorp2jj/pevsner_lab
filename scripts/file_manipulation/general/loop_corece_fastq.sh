#!/user/bin/env bash

for fq in ./PWS10CTRL_ACTTGA_L003_R1.fastq.gz
do
    python3 /mnt/data/jeremy/projects/sws/collab_samples/gatk_output/coerce_fastq.py $fq -o /mnt/data/jeremy/projects/sws/collab_samples/data/test
done
