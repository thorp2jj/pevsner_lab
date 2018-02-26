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

def try_python2_import(modules):
    for mod in modules:
        try:
            subprocess.check_call("python2 -c \"import imp; imp.find_module('{mod}')\"".format(mod=mod), shell=True)
        except subprocess.CalledProcessError:
            sys.exit("Error: requires python2 with the {mod} module".format(mod=mod))

def mh_test_env():
    test_env(["samtools", "bcftools"])

def parse_args():
    parser = argparse.ArgumentParser(description="Create a json config file for snakemake",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--bam", type=str, nargs='+', help="The input BAM file")
    parser.add_argument("--exome_region", required=True, help="BED file with capture regions used in WES.")
    parser.add_argument("--reference", type=str, default="/mnt/data/reference/hs37d5.fa", help="The reference aligned to.")
    parser.add_argument("--out_dir", type=str, default="./Data_Files", help="The output directory for the files")
    parser.add_argument("--config", help="The config file to write to")
    subparsers = parser.add_subparsers(dest="subparser_name")

    # The MH parser #
    mh_parser = subparsers.add_parser("MH")
    mh_parser.add_argument("--tmp_dir", type=str, default="tmp", help="A subdirectory for temporary files")
    mh_parser.add_argument("--dbsnp", type=str, default="/mnt/data/jeremy/programs/MosaicHunter/resources/dbsnp_137.b37.tsv", help="BED file of dbSNP common variant sites")
    mh_parser.add_argument("--repeats", type=str, required=True, help="A BED file of genomic repetitive regions")
    mh_parser.add_argument("--common", type=str, required=True, help="A BED file of error prone genomic regions")
    mh_parser.add_argument("--vcf_dir", type=str, required=True, help="VCF directory to identify sample sex")
    mh_parser.add_argument("--indels", type=str, required=True, help="BED file of Sample INDELs")
    mh_parser.add_argument("--max_depth", type=str, default="400", required=True, help="2X average sample sequencing depth")

    return parser.parse_args()

def main_MH(args):
    mh_test_env()
    
    args.directories.update( {
            "tmp": args.tmp_dir } )

    resource = {
        "dbsnp": args.dbsnp,
        "repeats": args.repeats,
        "common_sites": args.common }

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

        sample_vcf = args.vcf_dir + "/" + sample + ".vcf" 
        
        #indels = args.vcf_dir + "/" + sample + ".indels.recode.bed"
        indels = args.indels

        if not sample in samples:
            samples[sample] = []

        cur_sample = {
            "sample": sample,
            "bam": bam,
            "vcf": sample_vcf,
            "indels": indels }

        samples[sample].append(cur_sample)

    args.directories = {"out": os.path.abspath(args.out_dir)}

    args.configs = { "reference": args.reference,
                     "max_depth":args.max_depth,
                     "exome_region": os.path.abspath(args.exome_region),
                     "samples": samples}

    if args.subparser_name == "MH":
        main_MH(args)
    else:
        sys.exit("Unrecognized command")
    
    with open(args.config, 'w') as f:
        json.dump(args.configs, f, indent=2)
