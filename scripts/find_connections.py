#!/usr/bin/env python

import sys
import os
import argparse
from datetime import datetime

def is_file(filename):
    """Checks if a file exists"""
    if not os.path.isfile(filename):
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename


def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='find connections using a SNP distance cutoff')
    parser.add_argument("differences", help="output of SNP differences from snp dist",
        type=is_file)
    parser.add_argument("cutoff", help='SNP difference cutoff', type=int)
    return parser.parse_args()

args = get_args()
transmissions = {}
with open(args.differences, "r") as infile:
    for i,line in enumerate(infile):
        if i == 0:
            strains = line.strip().split("\t")[1:]
        else:
            line = line.strip().split("\t")
            strain1 = line[0]
            connected_strains = [strains[j] for j,x in enumerate(line[1:]) if int(x) <= args.cutoff]
            connected_strains.remove(strain1)
            transmissions[strain1] = connected_strains

connected_isolates = []
nodes = set(transmissions.keys())
while nodes:
    n = nodes.pop()
    group = {n}
    queue = [n]
    while queue:
        n = queue.pop(0)
        neighbors = set(transmissions[n])
        neighbors.difference_update(group)
        nodes.difference_update(neighbors)
        group.update(neighbors)
        queue.extend(neighbors)
    connected_isolates.append(group)

with open("{0}_{1}_clusters.txt".format(datetime.now().strftime("%Y-%m-%d"), args.differences[:-4]), "w") as outfile:
    for g in connected_isolates:
        outfile.write("\t".join(list(g)) + "\n")



