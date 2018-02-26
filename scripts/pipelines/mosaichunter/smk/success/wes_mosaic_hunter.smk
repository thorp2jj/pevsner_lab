'''
This pipeline runs MosaicHunter analysis on Whole Exome Sequencing.
1) import modules
2) create a config file
3) perform a dry-run with snakemake -npr ...
4) run the pipeline
'''

#include: '../rules/make_directories.smk'

rule all:
    input:
        mosaics=("{out}/{sample}/final.passed.tsv".format(out=config["dir"]["out"], sample=sample) for sample in config["samples"].keys())

rule MH:
    input:
        "{out}/{tmp}/{sample}_cleaned_sorted.bam".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], sample="{sample}"),
        sex="{out}/{tmp}/{sample}_sex".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], sample="{sample}"),
        dbsnp=config["resource"]["dbsnp"],
        common=config["resource"]["common_sites"],
        repetitive_regions=config["resource"]["repeats"],
        ref=config["reference"],
        bai="{out}/{tmp}/{sample}_cleaned_sorted.bai".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], sample="{sample}"),
        indels= lambda wildcards : config["samples"][wildcards.sample][int(0)]["indels"]        
    
    params:
        samp_out="{out}/{sample}".format(out=config["dir"]["out"], sample="{sample}") 

    output:
        "{out}/{sample}/final.passed.tsv".format(out=config["dir"]["out"], sample="{sample}")
    
    run:
        #sample_sex=str({input.sex})
        with open(input.sex, "r") as f:
            gender=f.readline()[0]
        f.close()
 
        shell(
        "java -jar /mnt/data/jeremy/programs/MosaicHunter/build/mosaichunter.jar exome "
        "-P input_file={input[0]} "
        "-P reference_file={input.ref} "
        "-P mosaic_filter.sex={gender} "
        "-P mosaic_filter.dbsnp_file={input.dbsnp} "
        "-P repetitive_region_filter.bed_file={input.repetitive_regions} "
        "-P indel_region_filter.bed_file={input.indels} "
        "-P common_site_filter.bed_file={input.common} "
        "-P output_dir={params.samp_out}") 


rule get_sex:
    input:
        lambda wildcards : config["samples"][wildcards.sample][int(0)]["vcf"]

    output:
        "{out}/{tmp}/{sample}_sex".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], sample="{sample}")

    shell:
        "bcftools +vcf2sex {input[0]} | cut -f 2 > {output}"

rule cleaned_bam_index:
    input:
        "{out}/{tmp}/{sample}_cleaned_sorted.bam".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], sample="{sample}")

    output:
        "{out}/{tmp}/{sample}_cleaned_sorted.bai".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], sample="{sample}")

    shell:
        "samtools index {input[0]} {output}" 

rule cleaned_bam_sort:
    input:
        "{out}/{tmp}/{sample}_cleaned.bam".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], sample="{sample}")

    output:
        "{out}/{tmp}/{sample}_cleaned_sorted.bam".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], sample="{sample}")

    shell:
        "samtools sort -o {output} -O bam {input[0]}"    

'''
rule get_depth:
    input:
        "{out}/{sample}_cleaned.bam".format(out=config["dir"]["out"], sample="{sample}")

    output:
        "{out}/{tmp}/{sample}_depth".format(out=config["dir"]["out"], sample="{sample}")

    shell:
        "samtools depth {input.vcf}  |  awk '{sum+=$3} END { print sum/NR}' > {output}"
'''

rule prepare_reads:
    input:
        lambda wildcards : config["samples"][wildcards.sample][int(0)]["bam"]

    params:
        sample=lambda wildcards : config["samples"][wildcards.sample][int(0)]["sample"]
    output:
        "{out}/{tmp}/{sample}_cleaned.bam".format(out=config["dir"]["out"], tmp=config["dir"]["tmp"], sample="{sample}")

    shell:
        "samtools view -h -f 0x2 {input[0]} | perl -ne 'print if (/^@/||(/NM:i:(\d+)/&&$1<=4))' | samtools view -Sb - >{output}"






