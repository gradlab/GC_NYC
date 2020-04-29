#!/usr/bin/env python

import argparse
import sys
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator
from collections import defaultdict


# This script assigns genogroups based on ngmaster output
# Genogroup is defined as identical tbpB allele and >99% similarity in porB
# Note that assemblies with "___" (3 underscores) in name will cause a problem

def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Assign genogroups from NG-MAST typing')
    parser.add_argument("ngmast", help="ngmast designation file")
    parser.add_argument("fasta", help="fasta of porB and tbpB alleles")
    return parser.parse_args()


def read_typing(ngmast):
    ngmast_dict = defaultdict(list)
    ngmast_count = defaultdict(int)
    isolates_dict = defaultdict(list)
    with open(ngmast, "r") as infile:
        for i,line in enumerate(infile):
            if i > 0:
                line = line.strip().split()
                if "multiple" in line:
                    print("{0} has multiple alleles. Will not be typed".format(line[0]))
                    continue
                if "new" in line:
                    line = [line[0] if x == "new" else x for x in line]
                if line[2] == "-":
                    print("{0} has no porB allele. Will not be typed".format(line[0]))
                    continue
                if line[3] == "-":
                    print("{0} has no tbpB allele. Will not be typed".format(line[0]))
                    continue
                isolates_dict[line[0]] = [line[1], line[2], line[3]]
                ngmast_dict[line[1]] = [line[2], line[3]]
                ngmast_count[line[1]] += 1
    return (ngmast_dict, ngmast_count, isolates_dict)

def allele_alignments(fasta, ngmast_dict):
    porB = []
    porB_length = []
    porB_types = []
    tbpB = []
    tbpB_length = []
    tbpB_types = []
    for rec in SeqIO.parse(fasta, "fasta"):
        if "POR" in rec.description:
            porB_length.append(len(rec.seq))
            try:
                if ngmast_dict[rec.id][1] not in porB_types:
                    porB.append(rec)
                    porB_types.append(ngmast_dict[rec.id][1])
            except IndexError:
                print("{0} has no types assigned.".format(rec.id))
                continue
        elif "TBPB" in rec.description:
            tbpB_length.append(len(rec.seq))
            try:
                if ngmast_dict[rec.id][2] not in tbpB_types:
                    tbpB.append(rec)
                    tbpB_types.append(ngmast_dict[rec.id][2])
            except IndexError:
                print("{0} has no types assigned.".format(rec.id))
                continue
        else:
            print("Warning: {0} did not have porB/tbpB designation".format(seq.id))
    if len(set(porB_length)) != 1:
        print("porB alleles are unequal lengths")
        sys.exit(1)
    if len(set(tbpB_length)) != 1:
        print("tbpB alleles are unequal lengths")
    porB_aln = MultipleSeqAlignment(porB)
    tbpB_aln = MultipleSeqAlignment(tbpB)
    return (porB_types, porB_aln, tbpB_types, tbpB_aln)

def allele_distances(porB_aln, tbpB_aln):
    calculator = DistanceCalculator('identity')
    porB_dm = calculator.get_distance(porB_aln)
    tbpB_dm = calculator.get_distance(tbpB_aln)
    return (porB_dm, tbpB_dm)


args = get_args()
ng_dict, ng_count, isolates_dict = read_typing(args.ngmast)
porB_types, porB_aln, tbpB_types, tbpB_aln = allele_alignments(args.fasta, isolates_dict)
porB_dm, tbpB_dm = allele_distances(porB_aln, tbpB_aln)


isolates = set([k for k in isolates_dict.keys() if len(isolates_dict[k]) == 3])
isolates_known = set([x for x in list(isolates) if isolates_dict[x][0].isdigit()])

genogroups_dict = defaultdict(list)

ng_count.pop('-')

# assign genogroups for isolates with known NG-MAST types

while len(isolates_known) > 0:
    # calculate unassigned NG-MAST type with the most isolates
    max_ngmast = max(ng_count, key=lambda key: ng_count[key])
    # assign new genogroup with max NG-MAST type
    max_ngmast_isolates = [y for y in list(isolates) if isolates_dict[y][0] == max_ngmast]
    genogroups_dict[max_ngmast].extend(max_ngmast_isolates)
    isolates = isolates - set(max_ngmast_isolates)
    isolates_known = isolates_known - set(max_ngmast_isolates)
    # find other isolates with identical tbpB
    same_tbpB = [z for z in list(isolates) if isolates_dict[z][2] == ng_dict[max_ngmast][1]]
    # check if porB is > 99% identical
    genogroup_isolates = []
    for s in same_tbpB:
        max_ngmast_porB = ng_dict[max_ngmast][0]
        isolate_porB = isolates_dict[s][1]
        distance = porB_dm[porB_types.index(max_ngmast_porB),porB_types.index(isolate_porB)]
        if distance < 0.01:
            genogroup_isolates.append(s)
    genogroups_dict[max_ngmast].extend(genogroup_isolates)
    isolates = isolates - set(genogroup_isolates)
    isolates_known = isolates_known - set(genogroup_isolates)
    ng_count.pop(max_ngmast)


# assign genogroups for isolates with novel porB/tbpB combinations

tbpB_count = defaultdict(int)
for t in isolates:
    tbpB_count[isolates_dict[t][1] + "___" + isolates_dict[t][2]] += 1

while len(isolates) > 0:
    # calculate unassigned tbpB allele with the most isolates
    max_new_ngMast = max(tbpB_count, key=lambda key: tbpB_count[key])
    max_porB, max_tbpB = max_new_ngMast.split("___")
    max_new_ngMast_isolates = [j for j in list(isolates) if isolates_dict[j][1] == max_porB and isolates_dict[j][2] == max_tbpB]
    genogroups_dict[max_new_ngMast].extend(max_new_ngMast_isolates)
    isolates = isolates - set(max_new_ngMast_isolates)
    max_tbpB_isolates = [m for m in list(isolates) if isolates_dict[m][2] == max_tbpB]
    # check if porB is > 99% identical
    known_tbpB_genogroup_isolates = []
    for n in max_tbpB_isolates:
        isolate_new_porB = isolates_dict[n][1]
        d = porB_dm[porB_types.index(max_porB), porB_types.index(isolate_new_porB)]
        if d < 0.01:
            known_tbpB_genogroup_isolates.append(n)
    genogroups_dict[max_new_ngMast].extend(known_tbpB_genogroup_isolates)
    isolates = isolates - set(known_tbpB_genogroup_isolates)
    tbpB_count.pop(max_new_ngMast)

with open("genogroups.txt", "w") as outfile:
    for gen in genogroups_dict:
        for genome in genogroups_dict[gen]:
            outfile.write("{0}\t{1}\n".format(genome, gen))
