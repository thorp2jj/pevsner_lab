#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='Sort Xhmm CNV output')
parser.add_argument('infile', help='input xcnv output file to sort')


args = parser.parse_args()

samples = {}

with open(args.infile) as f:
    for line in f:
        if line.startswith("SAMPLE"):        
            continue
        else:
            split_line=line.split("\t")
            interval_split=split_line[2].split(":")
            interval_chr=interval_split[0]
            
            interval_start=interval_split[1].split("-")[0]
            
            interval_end=interval_split[1].split("-")[1]
            #print(interval_chr, interval_start, interval_end)
            sample = split_line[0]
            if  sample  not in samples:
                samples[sample] = [[interval_chr, interval_start, interval_end]]
#                print(samples) 
            else:
                samples[sample].append([interval_chr, interval_start, interval_end])
#print(samples)

for key, value in samples.items():
    outfile_name = key + ".bed"
    outfile = open(outfile_name, 'w')
    for i in value:
        line = ""
        for a in i:
            line+= a + "\t" 
        outfile.write(line + "\n")
    outfile.close()

