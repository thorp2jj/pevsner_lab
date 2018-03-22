#!/usr/bin/env python3

import argparse
import gzip
import os.path
import re

barcode_end_regex = re.compile(b"[ \t\n\r\f\v]")

def process_args():
    parser = argparse.ArgumentParser(description="Split a fastq file into a single file for each run/lane")
    parser.add_argument("sample", help="The sample name")
    parser.add_argument("library", help="The sample library")
    parser.add_argument("infile", help="The input gzipped fastq file")
    parser.add_argument("out_dir", help="The output directory")
    return parser.parse_args()

def main(args):
    if not args:
        args = process_args()

    out_suffix = args.infile[-11:]
    out_filehandles = {}
    outfile = ''
    lane_pos = 0
    past_header = ''

    with gzip.open(args.infile, 'r') as f:
        for line in f:
            if not lane_pos or line[:lane_pos + 1] != past_header:
                machine_pos = line.index(b':')
                run_pos = line.index(b':', machine_pos + 1)
                flowcell_pos = line.index(b':', run_pos + 1)
                lane_pos = line.index(b':', flowcell_pos + 1)
                past_header = line[:lane_pos + 1]

                instrument = line[1:machine_pos]
                flowcell = line[1 + run_pos:flowcell_pos]
                lane = line[1 + flowcell_pos:lane_pos]

                # Get the barcode #
                barcode_pos = line.index(b' ')
                barcode_pos = line.index(b':', barcode_pos + 1)
                barcode_pos = line.index(b':', barcode_pos + 1)
                barcode_pos = line.index(b':', barcode_pos + 1)
                barcode_end = barcode_end_regex.search(line, barcode_pos + 1).start()
                barcode = line[barcode_pos + 1 : barcode_end]

                flowcell_barcode = flowcell + barcode

                out_name = '_'.join([args.sample, args.library] + [x.decode("UTF-8") for x in [instrument, flowcell_barcode, lane]])
                if not out_name in out_filehandles:
                    out_filehandles[out_name] = gzip.open(args.out_dir + '/' + out_name + '_' + out_suffix, 'wb')
                outfile = out_filehandles[out_name]
            
            outfile.write(line)
            outfile.write(f.readline())
            outfile.write(f.readline())
            outfile.write(f.readline())

    for outfile in out_filehandles.values():
        outfile.close()

if __name__ == "__main__":
    main(None)
