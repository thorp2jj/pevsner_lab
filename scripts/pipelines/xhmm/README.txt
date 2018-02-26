## XHMM Pipeline

## Requirements

Whole exome sequencing aligned, sorted, indexed BAMs
Exome capture interval list (format: chr:start-stop)
Human reference genome (BAM aligned reference)


##Software

module add xhmm/1.0 GATK/3.5-0 pseq/0.10 


##Create config file

```
mkdir -p Configs
python3 ./src/make_config.py \
    --bam /mnt/data/jeremy/scripts/pipelines/xhmm/test/*.bam \
    --exome_interval /mnt/data/jeremy/projects/smri/sureselect_v5_capture_targets/S04380110_Covered.interval_list \
    --out_dir /mnt/data/jeremy/scripts/pipelines/xhmm/test/run \
    --config xhmm_config.json

```


##Run pipeline

```
snakemake -w 30 --configfile xhmm_config.json -p -j 12 --nocolor --verbose -s ./smk/wes_xhmm.smk --cluster "sbatch -N 1 -n 1 -c {threads} --nodelist plab-node0[1,4] --nice=200 -e Logs/xhmm_{rule}_%j.log -o Logs/xhmm_{rule}_%j.out"
```





