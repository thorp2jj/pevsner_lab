#!/bin/sh
# properties = {"rule": "params_MH", "local": false, "input": ["/mnt/data/jeremy/scripts/pipelines/mosaichunter/Data_Files/tmp/TAC-94_cleaned_sorted.bam", "/mnt/data/jeremy/scripts/pipelines/mosaichunter/Data_Files/tmp/TAC-94_sex", "/mnt/data/jeremy/programs/MosaicHunter/resources/dbsnp_137.b37.tsv", "/mnt/data/jeremy/programs/MosaicHunter/resources/WES_Agilent_50M.error_prone.b37.bed", "/mnt/data/jeremy/programs/MosaicHunter/resources/all_repeats.b37.bed", "/mnt/data/reference/hs37d5.fa", "/mnt/data/jeremy/scripts/pipelines/mosaichunter/Data_Files/tmp/TAC-94_cleaned_sorted.bai", "/mnt/data/jeremy/scripts/pipelines/mosaichunter/test/TAC-94_indels_CNVs.bed"], "output": ["/mnt/data/jeremy/scripts/pipelines/mosaichunter/Data_Files/tmp/TAC-94_exomeparams"], "wildcards": ["TAC-94"], "params": {"max_depth_filt": 400}, "log": [], "threads": 1, "resources": {}, "jobid": 3, "cluster": {}}
cd /mnt/data/jeremy/scripts/pipelines/mosaichunter && \
/cm/shared/plab/apps/anaconda/3-5.0.0/bin/python \
-m snakemake /mnt/data/jeremy/scripts/pipelines/mosaichunter/Data_Files/tmp/TAC-94_exomeparams --snakefile /mnt/data/jeremy/projects/smri/mosaic_hunter/MH_pipe/smk/200x/wes_mosaic_hunter_200x.smk \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /mnt/data/jeremy/scripts/pipelines/mosaichunter/.snakemake/tmp.k97y6uly /mnt/data/jeremy/scripts/pipelines/mosaichunter/Data_Files/tmp/TAC-94_cleaned_sorted.bam /mnt/data/jeremy/scripts/pipelines/mosaichunter/Data_Files/tmp/TAC-94_sex /mnt/data/jeremy/programs/MosaicHunter/resources/dbsnp_137.b37.tsv /mnt/data/jeremy/programs/MosaicHunter/resources/WES_Agilent_50M.error_prone.b37.bed /mnt/data/jeremy/programs/MosaicHunter/resources/all_repeats.b37.bed /mnt/data/reference/hs37d5.fa /mnt/data/jeremy/scripts/pipelines/mosaichunter/Data_Files/tmp/TAC-94_cleaned_sorted.bai /mnt/data/jeremy/scripts/pipelines/mosaichunter/test/TAC-94_indels_CNVs.bed --latency-wait 60 \
--benchmark-repeats 1 --attempt 1 \
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
 --configfile /mnt/data/jeremy/scripts/pipelines/mosaichunter/test/mosaichunter_config.json  --nocolor \
--notemp --no-hooks --nolock --timestamp  --force-use-threads  --allowed-rules params_MH  && touch "/mnt/data/jeremy/scripts/pipelines/mosaichunter/.snakemake/tmp.k97y6uly/3.jobfinished" || (touch "/mnt/data/jeremy/scripts/pipelines/mosaichunter/.snakemake/tmp.k97y6uly/3.jobfailed"; exit 1)
