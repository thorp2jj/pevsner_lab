#!/usr/bin/env python3

import sys
version = sys.version_info
if version[0] < 3 or version[1] < 3:
    sys.exit("This script requires at least python 3.3")

import argparse
import json
import string
import random
import os.path
import shutil
import subprocess
import os
import re

def test_env(commands):
    for cmd in commands:
        if not shutil.which(cmd):
            sys.exit("Error: {cmd} not found in the PATH".format(cmd=cmd))

def mt2_test_env():
    test_env(["samtools", "picard"])

def parse_args():
    parser = argparse.ArgumentParser(description="Create a json config file for snakemake",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--bam", type=str, nargs='+', help="The input BAM file")
    parser.add_argument("--exome_region", required=True, help="BED file with capture regions used in WES.")
    parser.add_argument("--reference", type=str, default="/mnt/data/reference/hs37d5.fa", help="The reference aligned to.")
    parser.add_argument("--out_dir", type=str, default="./Data_Files", help="The output directory for the files")
    parser.add_argument("--config", help="The config file to write to")
    subparsers = parser.add_subparsers(dest="subparser_name")

    #  MT2 parser #
    mt2_parser = subparsers.add_parser("MT2")
    mt2_parser.add_argument("--tmp_dir", type=str, default="tmp", help="A subdirectory for temporary files")
    mt2_parser.add_argument("--germline_resource", type=str, default="/mnt/data/jeremy/resources/gnomad/af-only-gnomad.raw.sites.b37.vcf.gz", help="VCF file of GNOMAD common variant sites")
    mt2_parser.add_argument("--pon", type=str, required=True, help="A VCF file of Panel of Normal samples")
    mt2_parser.add_argument("--matched_index", type=str, required=True, help="File listing matched tumor-normal sample pairs")
    mt2_parser.add_argument("--matched_bam_dir", type=str, required = True, help="Directory containing matched sample BAM files")

    return parser.parse_args()

def main_MT2(args):
    mt2_test_env()
    
    args.directories.update( {
            "tmp": args.tmp_dir } )

    resource = {
        "germline_resource": args.germline_resource,
        "panel_of_normals": args.pon }

    args.configs.update( {
            "resource": resource,
            "dir": args.directories } )


if __name__ == "__main__":
    args = parse_args()

    samples = {}
    for i, sample_bam_file in enumerate(args.bam):
        
        bam = os.path.abspath(sample_bam_file)

        m = re.match(r".*/([\da-zA-Z-]+_.*)", sample_bam_file)
        sample = m.group(1)[:-10]

        matched_normal_sample = ""
        with open(args.matched_index, 'r') as y:
            for line in y:
                line = line.rstrip().split("\t")
                if line[0] == sample:
                    matched_normal_sample = line[1]
                else:
                    continue
        y.close()

        matched_normal_bam = args.matched_bam_dir + "/" + matched_normal_sample + "_recal.bam"

        if not sample in samples:
            samples[sample] = []

        cur_sample = {
            "sample": sample,
            "bam": bam,
            "matched_normal_sample": matched_normal_sample,
            "matched_normal_bam": matched_normal_bam }

        samples[sample].append(cur_sample)

    args.directories = {"out": os.path.abspath(args.out_dir)}

    args.configs = { "reference": args.reference,
                     "exome_region": os.path.abspath(args.exome_region),
                     "samples": samples}

    if args.subparser_name == "MT2":
        main_MT2(args)
    else:
        sys.exit("Unrecognized command")
    
    with open(args.config, 'w') as f:
        json.dump(args.configs, f, indent=2)
