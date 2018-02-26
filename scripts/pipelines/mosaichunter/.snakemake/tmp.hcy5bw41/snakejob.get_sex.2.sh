#!/bin/sh
# properties = {"rule": "get_sex", "local": false, "input": ["/mnt/data/jeremy/scripts/pipelines/mosaichunter/test/TAC-94.vcf.gz"], "output": ["/mnt/data/jeremy/scripts/pipelines/mosaichunter/Data_Files/tmp/TAC-94_sex"], "wildcards": ["TAC-94"], "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 2, "cluster": {}}
cd /mnt/data/jeremy/scripts/pipelines/mosaichunter && \
/cm/shared/plab/apps/anaconda/3-5.0.0/bin/python \
-m snakemake /mnt/data/jeremy/scripts/pipelines/mosaichunter/Data_Files/tmp/TAC-94_sex --snakefile /mnt/data/jeremy/projects/smri/mosaic_hunter/MH_pipe/smk/200x/wes_mosaic_hunter_200x.smk \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /mnt/data/jeremy/scripts/pipelines/mosaichunter/.snakemake/tmp.hcy5bw41 /mnt/data/jeremy/scripts/pipelines/mosaichunter/test/TAC-94.vcf.gz --latency-wait 60 \
--benchmark-repeats 1 --attempt 1 \
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
 --configfile /mnt/data/jeremy/scripts/pipelines/mosaichunter/test/mosaichunter_config.json  --nocolor \
--notemp --no-hooks --nolock --timestamp  --force-use-threads  --allowed-rules get_sex  && touch "/mnt/data/jeremy/scripts/pipelines/mosaichunter/.snakemake/tmp.hcy5bw41/2.jobfinished" || (touch "/mnt/data/jeremy/scripts/pipelines/mosaichunter/.snakemake/tmp.hcy5bw41/2.jobfailed"; exit 1)

