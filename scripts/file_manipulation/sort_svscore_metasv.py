#!/usr/bin/env python

from __future__ import print_function

import argparse

parser = argparse.ArgumentParser(description='Sort SVScore output')
parser.add_argument('infile', help='input svscore output file to sort')
parser.add_argument('outfile', help='file to write sorted svscore output to')

args = parser.parse_args()

header = []
vcfline = []
outdata = []

with open(args.infile) as f:
    for line in f:
        if line.startswith("#"):
            header.append(line.strip())
        else:
#       print(str(outdata))
            splitline = line.split()
#            print(splitline)
            info = splitline[7]
            info = info.split(";")
            for i in info:
                if i.startswith("SVSCORETOP10"):
                    score = float(i[13:])
                    vcfline.append([line, score])              
#                    print(vcfline)
            else:
                    continue
vcfline.sort(key=lambda x:x[1], reverse=True)
outfile = open(args.outfile, 'w')

for i in header:
    print(i, file=outfile)

for i in vcfline:
    print(i[0].strip(), file=outfile)


