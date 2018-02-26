#!/bin/sh
# properties = {"rule": "prepare_reads", "local": false, "input": ["/mnt/data/jeremy/scripts/pipelines/mosaichunter/test/TAC-94_recal.bam"], "output": ["/mnt/data/jeremy/scripts/pipelines/mosaichunter/Data_Files/tmp/TAC-94_cleaned.bam"], "wildcards": ["TAC-94"], "params": {"sample": "TAC-94"}, "log": [], "threads": 1, "resources": {}, "jobid": 6, "cluster": {}}
cd /mnt/data/jeremy/scripts/pipelines/mosaichunter && \
/cm/shared/plab/apps/anaconda/3-5.0.0/bin/python \
-m snakemake /mnt/data/jeremy/scripts/pipelines/mosaichunter/Data_Files/tmp/TAC-94_cleaned.bam --snakefile /mnt/data/jeremy/scripts/pipelines/mosaichunter/smk/wes_mosaic_hunter_200x.smk \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /mnt/data/jeremy/scripts/pipelines/mosaichunter/.snakemake/tmp.3x17h8sa /mnt/data/jeremy/scripts/pipelines/mosaichunter/test/TAC-94_recal.bam --latency-wait 60 \
--benchmark-repeats 1 --attempt 1 \
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
 --configfile /mnt/data/jeremy/scripts/pipelines/mosaichunter/test/mosaichunter_config.json  --nocolor \
--notemp --no-hooks --nolock --timestamp  --force-use-threads  --allowed-rules prepare_reads  && touch "/mnt/data/jeremy/scripts/pipelines/mosaichunter/.snakemake/tmp.3x17h8sa/6.jobfinished" || (touch "/mnt/data/jeremy/scripts/pipelines/mosaichunter/.snakemake/tmp.3x17h8sa/6.jobfailed"; exit 1)

