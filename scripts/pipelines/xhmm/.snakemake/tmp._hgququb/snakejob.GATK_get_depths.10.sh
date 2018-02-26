#!/bin/sh
# properties = {"rule": "GATK_get_depths", "local": false, "input": ["/mnt/data/jeremy/scripts/pipelines/xhmm/test/TAC-94_recal.bam"], "output": ["/mnt/data/jeremy/scripts/pipelines/xhmm/test/run/tmp/TAC-94.sample_interval_statistics"], "wildcards": ["TAC-94"], "params": {"sample": "TAC-94", "intervals": "/mnt/data/jeremy/projects/smri/sureselect_v5_capture_targets/S04380110_Covered.interval_list", "reference": "/mnt/data/reference/hs37d5.fa", "outname": "/mnt/data/jeremy/scripts/pipelines/xhmm/test/run/tmp/TAC-94"}, "log": [], "threads": 1, "resources": {}, "jobid": 10, "cluster": {}}
cd /mnt/data/jeremy/scripts/pipelines/xhmm && \
/cm/shared/plab/apps/anaconda/3-5.0.0/bin/python \
-m snakemake /mnt/data/jeremy/scripts/pipelines/xhmm/test/run/tmp/TAC-94.sample_interval_statistics --snakefile /mnt/data/jeremy/scripts/pipelines/xhmm/smk/wes_xhmm.smk \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /mnt/data/jeremy/scripts/pipelines/xhmm/.snakemake/tmp._hgququb /mnt/data/jeremy/scripts/pipelines/xhmm/test/TAC-94_recal.bam --latency-wait 30 \
--benchmark-repeats 1 --attempt 1 \
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
 --configfile /mnt/data/jeremy/scripts/pipelines/xhmm/xhmm_config.json -p --nocolor \
--notemp --no-hooks --nolock --timestamp  --printshellcmds  --force-use-threads  --allowed-rules GATK_get_depths  && touch "/mnt/data/jeremy/scripts/pipelines/xhmm/.snakemake/tmp._hgququb/10.jobfinished" || (touch "/mnt/data/jeremy/scripts/pipelines/xhmm/.snakemake/tmp._hgququb/10.jobfailed"; exit 1)

