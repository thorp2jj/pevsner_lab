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

def xhmm_test_env():
    test_env(["xhmm", "pseq"])

def parse_args():
    parser = argparse.ArgumentParser(description="Create a json config file for snakemake",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--bam", type=str, nargs='+', help="The input BAM file")
    parser.add_argument("--exome_interval", required=True, help="Interval file with capture regions used in WES.")
    parser.add_argument("--reference", type=str, default="/mnt/data/reference/hs37d5.fa", help="The reference aligned to.")
    parser.add_argument("--out_dir", type=str, default="./Data_Files", help="The output directory for the files")
    parser.add_argument("--config", type=str, default="xhmm_config.json", help="The config file to write to")
    parser.add_argument("--tmp_dir", type=str, default="tmp", help="A subdirectory for temporary files")

    return parser.parse_args()

if __name__ == "__main__":
    xhmm_test_env()
    
    args = parse_args()

    samples = {}
    for i, sample_bam_file in enumerate(args.bam):
        
        bam = os.path.abspath(sample_bam_file)

        m = re.match(r".*/([\da-zA-Z-]+_.*)", sample_bam_file)
        sample = m.group(1)[:-10]

        if not sample in samples:
            samples[sample] = []

        cur_sample = {
            "sample": sample,
            "bam": bam}

        samples[sample].append(cur_sample)

    args.directories = {"out": os.path.abspath(args.out_dir)}
    
    args.directories.update({"tmp": args.tmp_dir})

    args.configs = { "reference": args.reference,
                     "exome_region": os.path.abspath(args.exome_interval),
                     "dir": args.directories,
                     "samples": samples}

    with open(args.config, 'w') as f:
        json.dump(args.configs, f, indent=2)
