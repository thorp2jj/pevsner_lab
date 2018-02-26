#!/user/bin/env bash

for i in `cat /mnt/data/jeremy/projects/smri/Data_Files/BAM/100X/100x_samplelist_bai.txt`
do
    ln $i .
done

#find /mnt/data/raw_data/1704KHF-0126 -type f -iname "*.fastq.gz" -exec ln -s {} . \;
