Mutect2 Snakemake Analysis

Require modules

module add JDK/1.8.0_11 picard/2.09.0 GATK/4.0.1.0 samtools/1.3.1


Create panel of normals
Save a list of BAM files, used as "normal" samples, with absolute paths to file (sample_vcfs.list)

srun java -jar $GATK CreateSomaticPanelOfNormals -vcfs samples_vcfs.list -O samples_pon.vcf.gz


Create config file

python3 ./src/make_config.py --bam /mnt/data/jeremy/projects/smri/100X/Data_Files/BAM/CM*.bam 
--exome_region /mnt/data/jeremy/projects/smri/sureselect_v5_capture_targets/S04380110_Padded.bed
--config CM-T_100X.json MT2 --pon /mnt/data/jeremy/projects/smri/mutect2/panel_of_normals/100X/TAC_norms/TAC_100X_pon.vcf.gz
--matched_index /mnt/data/DNASeq/BPD_BSMN/stanley_collection/Matched_Sample_Info/Corrected_File_Matching/Matched_Stanley_IDs_FLIPFLOP_corrected.tsv 
--matched_bam_dir /mnt/data/jeremy/projects/smri/100X/Data_Files/BAM


Run Mutect2 snakemake

snakemake -npr -w 60 --configfile CM-T_100X.json -j 30 --nocolor --verbose -s smk/wes_mutect2_norm_var_removed.smk 
--cluster "sbatch -N 1 -n 1 --error /Logs/MT2_{rule}_%j.log --output /Logs/MT2_{rule}%j.log"








