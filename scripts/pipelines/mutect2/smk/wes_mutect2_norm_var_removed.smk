'''
This pipeline runs Mutect2 analysis on Whole Exome Sequencing.
1) import modules
2) create panel of normals
3) create config file
4) perform a dry-run with snakemake -npr ...
5) run the pipeline
'''

include: '../rules/make_directories.smk'

rule all:
    input:
        mosaics=("{out}/{tmp}/{sample}_complete.txt".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], sample=sample) for sample in config["samples"].keys())


rule completeness:
     input:
        "{out}/{samp_dir}/{sample}_filtered_oxogfilt.vcf.gz".format(out=config["dir"]["out"], samp_dir="{sample}", sample="{sample}")

     output:
        "{out}/{tmp}/{sample}_complete.txt".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], sample="{sample}")

     shell:
        "touch {output}"


rule filter_orientation_bias:
    input:
        "{out}/{samp_dir}/{sample}_filtered.vcf.gz".format(out=config["dir"]["out"], samp_dir="{sample}", sample="{sample}"),
        "{out}/{tmp}/{sample}_artifact_metrics.pre_adapter_detail_metrics".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], sample="{sample}") 

    output:
        "{out}/{samp_dir}/{sample}_filtered_oxogfilt.vcf.gz".format(out=config["dir"]["out"], samp_dir="{sample}", sample="{sample}")

    shell:
        "java -jar $GATK FilterByOrientationBias --artifact-modes 'G/T' -V {input[0]} -P {input[1]} --output {output}"


rule generate_bamouts:
    input:
        lambda wildcards : config["samples"][wildcards.sample][int(0)]["bam"],
        lambda wildcards : config["samples"][wildcards.sample][int(0)]["matched_normal_bam"]

    params:
        reference=config["reference"],
        exome_interval=config["exome_region"],
        germline_resource=config["resource"]["germline_resource"],
        tumor_sample=lambda wildcards : config["samples"][wildcards.sample][int(0)]["sample"],
        normal_sample=lambda wildcards : config["samples"][wildcards.sample][int(0)]["matched_normal_sample"],
        panel_of_normals=config["resource"]["panel_of_normals"] 

    output:
        "{out}/{samp_dir}/{sample}_BAMOUT.bam".format(out=config["dir"]["out"], samp_dir="{sample}", sample="{sample}"),
        "{out}/{samp_dir}/{sample}_BAMOUT.vcf.gz".format(out=config["dir"]["out"], samp_dir="{sample}", sample="{sample}")

    shell:
        "java -jar $GATK Mutect2 -R {params.reference} -I {input[0]} -tumor {params.tumor_sample} -I {input[1]} -normal {params.normal_sample} "
        "-L {params.exome_interval} --germline-resource {params.germline_resource} --panel-of-normals {params.panel_of_normals} "
        "-bamout {output[0]} --dont-trim-active-regions --active-probability-threshold 0.000 -ip 300 -O {output[1]}"


rule filter_calls:
    input:
        "{out}/{samp_dir}/{sample}.vcf.gz".format(out=config["dir"]["out"], samp_dir="{sample}", sample="{sample}"),
        "{out}/{tmp}/{sample}_contamination.table".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], sample="{sample}")

    output:
        "{out}/{samp_dir}/{sample}_filtered.vcf.gz".format(out=config["dir"]["out"], samp_dir="{sample}", sample="{sample}")

    shell:
        "java -jar $GATK FilterMutectCalls -V {input[0]} -contaminationTable {input[1]} -O {output}"


rule sequencing_artifacts:
    input:
        lambda wildcards : config["samples"][wildcards.sample][int(0)]["bam"]

    params:
        reference=config["reference"],
        outname = "{out}/{tmp}/{sample}_artifact_metrics".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], sample="{sample}")

    output:
        "{out}/{tmp}/{sample}_artifact_metrics.pre_adapter_detail_metrics".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], sample="{sample}")

    shell:
        "java -jar -Djava.io.tmpdir=/mnt/data/jeremy/scratch/tmp /cm/shared/plab/apps/picard/picard-tools-2.09.0/picard.jar CollectSequencingArtifactMetrics I={input[0]} "
        "O={params.outname} R={params.reference}"


rule estimate_contamination:
    input:
        "{out}/{tmp}/{sample}_pileupsummary.table".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], sample="{sample}")

    output:
         "{out}/{tmp}/{sample}_contamination.table".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], sample="{sample}")

    shell:
        "java -jar $GATK CalculateContamination -I {input[0]} -O {output}" 


rule getpileupsummaries:
    input:
        lambda wildcards : config["samples"][wildcards.sample][int(0)]["bam"]

    params:
        reference=config["reference"]

    output:
        "{out}/{tmp}/{sample}_pileupsummary.table".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], sample="{sample}")

    shell:
        "java -jar $GATK -R {params.reference} -I {input} -V /mnt/data/jeremy/resources/gnomad/small_exac_common_3.vcf -O {output}"    


rule call_somatic_vars:
    input:
        lambda wildcards : config["samples"][wildcards.sample][int(0)]["bam"],
        lambda wildcards : config["samples"][wildcards.sample][int(0)]["matched_normal_bam"]

    params:
        reference=config["reference"],
        exome_interval=config["exome_region"],
        germline_resource=config["resource"]["germline_resource"],
        panel_of_normals=config["resource"]["panel_of_normals"],
        tumor_sample = lambda wildcards : config["samples"][wildcards.sample][int(0)]["sample"],
        normal_sample = lambda wildcards : config["samples"][wildcards.sample][int(0)]["matched_normal_sample"]

    output:
        "{out}/{samp_dir}/{sample}.vcf.gz".format(out=config["dir"]["out"], samp_dir="{sample}", sample="{sample}")

    shell:
        "java -jar $GATK Mutect2 -I {input[0]} -tumor {params.tumor_sample} -I {input[1]} -normal {params.normal_sample} -O {output} -R {params.reference} -L {params.exome_interval} "
        "--germline-resource {params.germline_resource} --panel-of-normals {params.panel_of_normals}"


