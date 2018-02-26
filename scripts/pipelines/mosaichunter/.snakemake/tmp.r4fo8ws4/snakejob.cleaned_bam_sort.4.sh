#!/bin/sh
# properties = {"rule": "cleaned_bam_sort", "local": false, "input": ["/mnt/data/jeremy/scripts/pipelines/mosaichunter/Data_Files/tmp/TAC-94_cleaned.bam"], "output": ["/mnt/data/jeremy/scripts/pipelines/mosaichunter/Data_Files/tmp/TAC-94_cleaned_sorted.bam"], "wildcards": ["TAC-94"], "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 4, "cluster": {}}
cd /mnt/data/jeremy/scripts/pipelines/mosaichunter && \
/cm/shared/plab/apps/anaconda/3-5.0.0/bin/python \
-m snakemake /mnt/data/jeremy/scripts/pipelines/mosaichunter/Data_Files/tmp/TAC-94_cleaned_sorted.bam --snakefile /mnt/data/jeremy/scripts/pipelines/mosaichunter/smk/wes_mosaic_hunter_200x.smk \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /mnt/data/jeremy/scripts/pipelines/mosaichunter/.snakemake/tmp.r4fo8ws4 /mnt/data/jeremy/scripts/pipelines/mosaichunter/Data_Files/tmp/TAC-94_cleaned.bam --latency-wait 30 \
--benchmark-repeats 1 --attempt 1 \
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
 --configfile /mnt/data/jeremy/scripts/pipelines/mosaichunter/test/mosaichunter_config.json -p --nocolor \
--notemp --no-hooks --nolock --timestamp  --printshellcmds  --force-use-threads  --allowed-rules cleaned_bam_sort  && touch "/mnt/data/jeremy/scripts/pipelines/mosaichunter/.snakemake/tmp.r4fo8ws4/4.jobfinished" || (touch "/mnt/data/jeremy/scripts/pipelines/mosaichunter/.snakemake/tmp.r4fo8ws4/4.jobfailed"; exit 1)

