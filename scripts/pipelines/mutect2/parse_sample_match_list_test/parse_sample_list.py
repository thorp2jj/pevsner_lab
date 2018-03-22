#!/usr/bin/env python


with open("/mnt/data/DNASeq/BPD_BSMN/stanley_collection/Matched_Sample_Info/Corrected_File_Matching//mnt/data/DNASeq/BPD_BSMN/stanley_collection/Matched_Sample_Info/Corrected_File_Matching/Matched_Stanley_IDs_FLIPFLOP_corrected.tsv", 'r') as y:
    for line in y:
        line = line.rstrip().split("\t")
        if line[0] == :
            matched_normal_sample = line[1]
        else:
            continue

print("hello")
print(matched_normal_sample)
y.close()





