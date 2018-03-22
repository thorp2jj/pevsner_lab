#!/usr/bin/env python

with open("TAC-94_exomeparams.txt", "r") as y:
    for line in y:
        line = line.strip()
        if line.startswith('alpha'):
            alpha = line[7:]
            print(alpha)
        elif line.startswith('beta'):
            beta = line[6:]
            print(beta)
        elif line.startswith('average depth'):
            depth = line[15:]
            print(depth)
        else:
            continue

y.close()






