#!/usr/bin/env python3

import argparse
import sys
import shutil
import subprocess

male_min_y_fraction = 0.003
female_max_y_fraction = 0.0011

male_max_x_fraction = 0.035
female_min_x_fraction = 0.04

def process_args():
    parser = argparse.ArgumentParser(description="Determine the sex of samples based on the fraction of reads aligned to the X and Y chromosomes")
    parser.add_argument("outfile", help="The output text file of sample sex")
    parser.add_argument("infiles", nargs='+', help="The input BAM files")
    samtools = shutil.which("samtools")
    if not samtools:
        sys.exit("The samtools excutable is not in the PATH")
    return parser.parse_args()

def main(args):
    if not args:
        args = process_args()

    output = []
    for bam_file in args.infiles:
        sample_name = subprocess.check_output("""samtools view -H {bam} | grep "^@RG" | sed 's/.*SM://' | sed 's/\s.*//' | uniq""".format(bam=bam_file), shell=True, universal_newlines=True).rstrip()
        idxstats_out = subprocess.check_output("samtools idxstats {bam}".format(bam=bam_file), shell=True, universal_newlines=True)
        x_reads, y_reads, total_reads = 0, 0, 0
        for line in idxstats_out.split('\n'):
            line = line.rstrip().split('\t')
            if not line or not line[0]:
                continue
            reads = int(line[2])
            total_reads += reads
            if line[0] == "X" or line[0] == "chrX":
                x_reads = reads
            elif line[0] == "Y" or line[0] == "chrY":
                y_reads = reads
        
        y_fraction = y_reads / total_reads
        x_fraction = x_reads / total_reads
        
        if y_fraction > male_min_y_fraction and x_fraction < male_max_x_fraction:
            output.append(sample_name + '\t' + 'M')
        elif y_fraction < female_max_y_fraction and x_fraction > female_min_x_fraction:
            output.append(sample_name + '\t' + 'F')
        else:
            sys.exit("Sample {sample} in file {file} has an unknown sex".format(sample=sample_name, file=bam_file))
    with open(args.outfile, 'w') as f:
        print('\n'.join(output), file=f)

if __name__ == "__main__":
    main(None)
