#!/usr/bin/env python3

import argparse

chroms = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
          '13', '14', '15', '16', '17', '18', '19', '20', '21', '22',
           'X', 'Y', 'MT', 'M'}

def process_args():
    parser = argparse.ArgumentParser(description="Make an interval list of human autosomes and sex chromosomes from a reference")
    parser.add_argument("infile", help="The reference fasta index (.fai)")
    parser.add_argument("outfile", help="The list of cononcial chromsomes")
    return parser.parse_args()

def main(args):
    if not args:
        args = process_args()

    outfile = open(args.outfile, 'w')

    with open(args.infile, 'r') as f:
        for line in f:
            chrom = line.rstrip().split('\t')[0]
            if chrom in chroms or "chr" + chrom in chroms:
                print(chrom, file=outfile)
    
    outfile.close()

if __name__ == "__main__":
    main(None)
