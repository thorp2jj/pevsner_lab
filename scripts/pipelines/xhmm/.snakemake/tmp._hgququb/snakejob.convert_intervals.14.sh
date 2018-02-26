#!/bin/sh
# properties = {"rule": "convert_intervals", "local": false, "input": ["/mnt/data/jeremy/projects/smri/sureselect_v5_capture_targets/S04380110_Covered.interval_list"], "output": ["/mnt/data/jeremy/scripts/pipelines/xhmm/test/run/DATA.EXOME.targets.reg"], "wildcards": [], "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 14, "cluster": {}}
cd /mnt/data/jeremy/scripts/pipelines/xhmm && \
/cm/shared/plab/apps/anaconda/3-5.0.0/bin/python \
-m snakemake /mnt/data/jeremy/scripts/pipelines/xhmm/test/run/DATA.EXOME.targets.reg --snakefile /mnt/data/jeremy/scripts/pipelines/xhmm/smk/wes_xhmm.smk \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /mnt/data/jeremy/scripts/pipelines/xhmm/.snakemake/tmp._hgququb /mnt/data/jeremy/projects/smri/sureselect_v5_capture_targets/S04380110_Covered.interval_list --latency-wait 30 \
--benchmark-repeats 1 --attempt 1 \
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
 --configfile /mnt/data/jeremy/scripts/pipelines/xhmm/xhmm_config.json -p --nocolor \
--notemp --no-hooks --nolock --timestamp  --printshellcmds  --force-use-threads  --allowed-rules convert_intervals  && touch "/mnt/data/jeremy/scripts/pipelines/xhmm/.snakemake/tmp._hgququb/14.jobfinished" || (touch "/mnt/data/jeremy/scripts/pipelines/xhmm/.snakemake/tmp._hgququb/14.jobfailed"; exit 1)

