#!/usr/bin/env python3

import argparse
import sys
import math


parser = argparse.ArgumentParser(description='Parse variants from a VCF file to BED format')
parser.add_argument('invcf', metavar='.vcf', help='A VCF format file.')
#parser.add_argument('outfile', metavar='.bed', help='A file to output variant coordinates in BED format')
args = parser.parse_args()

with open(args.invcf, 'r') as f:
    for line in f:
        if line[0:2] == '##': #Skip header
            continue
        elif line[0:4] == '#CHR': #Populate sample list
            continue
        else:
            line = line.rstrip().split()
            chr = line[0]
            orig_pos = line[1]
            new_start = str(int(orig_pos) - 5)
            new_end = str(int(orig_pos) + 5)
            print (chr,"\t",new_start,"\t",new_end)

