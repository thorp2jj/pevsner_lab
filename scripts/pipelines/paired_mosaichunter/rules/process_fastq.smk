
rule all:
    input:
        "{out}/finished.txt".format(out=config["out_dir"])

def get_all_output(wildcards):
    infiles = []
    out_dir = config["out_dir"]
    for sample in config["Samples"].keys():
        for library in config["Samples"][sample].keys():
            for i, fq1 in enumerate(config["Samples"][sample][library][0]):
                infiles.append("{}/{}_{}_{}_1_finished.txt".format(out_dir, sample, library, i))
            for i, fq2 in enumerate(config["Samples"][sample][library][1]):
                infiles.append("{}/{}_{}_{}_2_finished.txt".format(out_dir, sample, library, i))
    return infiles

rule collector:
    input:
        infiles = get_all_output
    output:
        "{out}/finished.txt".format(out=config["out_dir"])
    run:
        infiles = " ".join(input.infiles)
        shell("touch {output}")
        shell("rm {infiles}")
        shell("sleep 5")

def get_fastq(wildcards):
    ori = int(wildcards.ori) - 1
    idx = int(wildcards.idx)
    return config["Samples"][wildcards.sample][wildcards.library][ori][idx]

rule split_fastq_1:
    input:
        fastq = get_fastq,
        split_fastq = config["split_fastq"]
    output:
        "{out}/{sample}_{library}_{idx}_{ori}_finished.txt".format(out=config["out_dir"], sample="{sample}", library="{library}", idx="{idx}", ori="{ori}")
    run:
        sample = wildcards.sample
        library = wildcards.library
        out_dir = config["out_dir"]
        shell("python3 {input.split_fastq} {sample} {library} {input.fastq} {out_dir}")
        shell("touch {output}")
        shell("sleep 5")
