### MosaicHunter Pipeline

### Requirements

Reference file (.fa & .fai)
Exome: selected regions file (.bed)
Repeat region file (.bed)
Error prone region file (.bed)
Aligned sample reads (.bam)
Called sample INDELS (.bed)
Called sample CNVs (.bed)
Called sample SNP/INDELS (bgzipped .vcf.gz)
Indexed sample SNP/INDELS (bcftools index $sample, .vcf.gz.csi)

### Software
'''
module add samtools/1.3.1 mosaichunter/1.0 bcftools/1.3 tabix_bgzip/1.3.1 
export BCFTOOLS_PLUGINS=/cm/shared/plab/apps/bcftools/bcftools-1.3/plugins/
'''
### Quickstart

Here are quickstart instructions to use this pipeline: 


#### Create config file

```
mkdir -p Configs
python3 /src/make_config.py --bam path/*.bam \ 
    --exome_region path/regions.bed \
    --config path/config.json \
    MH --repeats path/repeats.bed \
    --common path/error_prone.bed \
    --vcf_dir path/out_dir \ 
    --indels path/cnv_indel.bed 

#EXAMPLE
python3 ./src/make_config.py --bam /mnt/data/jeremy/scripts/pipelines/mosaichunter/test/*.bam \
    --exome_region /mnt/data/DNASeq/BPD_BSMN/stanley_collection/sureselect_v5_capture_targets/S04380110_Covered.bed \
    --config /mnt/data/jeremy/scripts/pipelines/mosaichunter/test/mosaichunter_config.json MH \
    --repeats /mnt/data/jeremy/programs/MosaicHunter/resources/all_repeats.b37.bed \
    --common /mnt/data/jeremy/programs/MosaicHunter/resources/WES_Agilent_50M.error_prone.b37.bed \
    --vcf_dir /mnt/data/jeremy/scripts/pipelines/mosaichunter/test --indels /mnt/data/jeremy/scripts/pipelines/mosaichunter/test/TAC-94_indels_CNVs.bed / 
    --max_depth=400
```

#### Run pipeline through snakemake

```
mkdir -p ./Logs
snakemake -npr -w 60 --configfile ../config.json -j 30 --nocolor --verbose -s wes_mosaic_hunter.smk --cluster "sbatch -N 1 -n 1 --error /Logs/MH_{rule}_%j.log --output /Logs/MH_{rule}_%j.log"

#EXAMPLE
snakemake -npr -w 60 --configfile ../test/mosaichunter_config.json -j 30 --nocolor --verbose -s wes_mosaic_hunter.smk --cluster "sbatch -N 1 -n 1 --error /Logs/MH_%j.log --output /Logs/MH_%j.log"
```
