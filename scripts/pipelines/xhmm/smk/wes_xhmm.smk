'''
This pipeline runs XHMM analysis on Whole Exome Sequencing.
1) import modules
2) create a config file
3) perform a dry-run with snakemake -npr ...
4) run the pipeline
'''

#include: '../rules/make_directories.smk'

rule all:
    input:
        cnvs=("{out}/DATA.xcnv".format(out=config["dir"]["out"]))

'''
rule genotype_cnvs:
    input:
    
    output:
    
    shell:
'''

rule discover_cnvs:
    input:
        "{out}/DATA.PCA_normalized.filtered.sample_zscores.RD.txt".format(out=config["dir"]["out"]),
        "{out}/DATA.sample_filtered.RD.txt".format(out=config["dir"]["out"]) 

    params:
        out_directory=config["dir"]["out"]
    
    output:
        "{out}/DATA.xcnv".format(out=config["dir"]["out"]),
        "{out}/DATA.aux_xcnv".format(out=config["dir"]["out"])

    shell:
        "xhmm --discover -p /cm/shared/plab/apps/xhmm/1.0/params.txt "
        "-r {input[0]} "
        "-R {input[1]} "
        "-c {output[0]} "
        "-a {output[1]} "
        "-s {params.out_directory}"


rule filter_read_depth:
    input:
        "{out}/DATA.RD.txt".format(out=config["dir"]["out"]),
        "{out}/DATA.filtered_centered.RD.txt.filtered_targets.txt".format(out=config["dir"]["out"]),
        "{out}/DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt".format(out=config["dir"]["out"]),
        "{out}/DATA.filtered_centered.RD.txt.filtered_samples.txt".format(out=config["dir"]["out"]),
        "{out}/DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt".format(out=config["dir"]["out"]) 

    output:
        "{out}/DATA.sample_filtered.RD.txt".format(out=config["dir"]["out"])

    shell:
        "xhmm --matrix -r {input[0]} "
        "--excludeTargets {input[1]} "
        "--excludeTargets {input[2]} "
        "--excludeSamples {input[3]} "
        "--excludeSamples {input[4]} "
        "-o {output}"


rule filter_zscore_PCA_norm:
    input:
        "{out}/DATA.PCA_normalized.txt".format(out=config["dir"]["out"])

    output:
        "{out}/DATA.PCA_normalized.filtered.sample_zscores.RD.txt".format(out=config["dir"]["out"]),
        "{out}/DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt".format(out=config["dir"]["out"]),
        "{out}/DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt".format(out=config["dir"]["out"])

    shell:
        "xhmm --matrix -r {input} "
        "--centerData --centerType sample --zScoreData "
        "-o {output[0]} --outputExcludedTargets {output[1]} --outputExcludedSamples {output[2]} "
        "--maxSdTargetRD 30"

rule norm_PCA:
    input:
        "{out}/DATA.filtered_centered.RD.txt".format(out=config["dir"]["out"]),
        "{out}/DATA.RD_PCA.PC.txt".format(out=config["dir"]["out"]) 

    params:
        "{out}/DATA.RD_PCA".format(out=config["dir"]["out"])

    output:
        "{out}/DATA.PCA_normalized.txt".format(out=config["dir"]["out"])

    shell:
        "xhmm --normalize -r {input[0]} "
        "--PCAfiles {params[0]} "
        "--normalizeOutput {output} "
        "--PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7"


rule PCA:
    input:
        "{out}/DATA.filtered_centered.RD.txt".format(out=config["dir"]["out"])

    params:
        "{out}/DATA.RD_PCA".format(out=config["dir"]["out"])

    output:
        "{out}/DATA.RD_PCA.PC.txt".format(out=config["dir"]["out"])

    shell:
        "xhmm --PCA -r {input[0]} --PCAfiles {params[0]}"


rule filter_samp_targ_meancenter:
    input:
        "{out}/DATA.RD.txt".format(out=config["dir"]["out"]),
        "{out}/DATA.low_complexity_targets.txt".format(out=config["dir"]["out"]),
        "{out}/extreme_gc_targets.txt".format(out=config["dir"]["out"])

    output:
        "{out}/DATA.filtered_centered.RD.txt".format(out=config["dir"]["out"]),
        "{out}/DATA.filtered_centered.RD.txt.filtered_targets.txt".format(out=config["dir"]["out"]),
        "{out}/DATA.filtered_centered.RD.txt.filtered_samples.txt".format(out=config["dir"]["out"])

    shell:
        "xhmm --matrix -r {input[0]} --centerData --centerType target "
        "-o {output[0]} --outputExcludedTargets {output[1]} --outputExcludedSamples {output[2]} "
        "--excludeTargets {input[2]} --excludeTargets {input[1]} "
        "--minTargetSize 10 --maxTargetSize 10000 "
        "--minMeanTargetRD 10 --maxMeanTargetRD 500 "
        "--minMeanSampleRD 25 --maxMeanSampleRD 200 "
        "--maxSdSampleRD 150"


rule generate_low_complex_targets:
    input:
        "{out}/DATA.locus_complexity.txt".format(out=config["dir"]["out"])

    output:
        "{out}/DATA.low_complexity_targets.txt".format(out=config["dir"]["out"])

    shell:
        "cat {input} | awk '{{if ($2 > 0.25) print $1}}' > {output}"


rule pseq_loc_stats:
    input:
        "{out}/DATA.EXOME.targets.LOCDB".format(out=config["dir"]["out"])

    params:
        intervals=config["exome_region"]        

    output:
        "{out}/DATA.locus_complexity.txt".format(out=config["dir"]["out"])

    shell:
        "pseq . loc-stats --locdb {input} --group targets --seqdb /mnt/data/jeremy/resources/pseq/seqdb.hg19 | awk '{{if (NR > 1) print $_}}' | sort -k1 -g | awk '{{print $10}}' | paste {params.intervals} - | awk '{{print $1\"\\t\"$2}}' > {output}"


rule pseq_loc_load:
    input:
        "{out}/DATA.EXOME.targets.reg".format(out=config["dir"]["out"])

    params:
        "{out}/DATA.EXOME.targets.LOCDB.loc-load".format(out=config["dir"]["out"])

    output:
        "{out}/DATA.EXOME.targets.LOCDB".format(out=config["dir"]["out"])

    shell:
        "pseq . loc-load --locdb {output[0]} --file {input[0]} --group targets --out {params[0]}"


rule convert_intervals:
    input:
        intervals=config["exome_region"]
                
    output:
        "{out}/DATA.EXOME.targets.reg".format(out=config["dir"]["out"])

    shell:
        "interval_list_to_pseq_reg {input.intervals} > {output}"


rule extreme_GC:
    input:
        "{out}/DATA.locus_GC.txt".format(out=config["dir"]["out"])

    output:
        "{out}/extreme_gc_targets.txt".format(out=config["dir"]["out"])        

    shell:
        "cat {input} | awk '{{if ($2 < 0.1 || $2 > 0.9) print $1}}' > {output}"


rule per_target_GC:
    input:
        "{out}/DATA.RD.txt".format(out=config["dir"]["out"])

    params:
        intervals=config["exome_region"],
        reference=config["reference"]

    output:
        "{out}/DATA.locus_GC.txt".format(out=config["dir"]["out"])

    shell:
        "java -jar -Xmx4g $GATK "
        "-T GCContentByInterval "
        "-L {params.intervals} "
        "-R {params.reference} "
        "-o {output}"


rule GATK_merge_depths:
    input:
         ("{out}/{tmp}/{sample}.sample_interval_summary".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], sample=sample) for sample in config["samples"].keys())

    output:
        "{out}/DATA.RD.txt".format(out=config["dir"]["out"])

    run:
        arg=""
        for samp in input:
             arg = arg + "--GATKdepths " + samp
        shell(
            "xhmm --mergeGATKdepths -o {output} {arg}"
             )



rule GATK_get_depths:
    input:
        lambda wildcards : config["samples"][wildcards.sample][int(0)]["bam"]

    output:
        #"{out}/{tmp}/{sample}".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], sample="{sample}"),
        "{out}/{tmp}/{sample}.sample_interval_summary".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], sample="{sample}")

    params:
        sample=lambda wildcards : config["samples"][wildcards.sample][int(0)]["sample"],
        intervals=config["exome_region"],
        reference=config["reference"],
        outname="{out}/{tmp}/{sample}".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], sample="{sample}")

    run:
        shell(
            "java -jar -Xmx8g -Djava.io.tmpdir=/mnt/data/jeremy/scratch/tmp $GATK "
            "-T DepthOfCoverage "
            "-R {params.reference} "
            "-o {params.outname} "
            "-I {input[0]} "
            "-dt BY_SAMPLE -dcov 5000 -l INFO --omitDepthOutputAtEachBase "
            "--minBaseQuality 0 --minMappingQuality 20 --start 1 --stop 5000 --nBins 200 "
            "--omitLocusTable --includeRefNSites "
            "--countType COUNT_FRAGMENTS "
            "-L {params.intervals}"
             )




