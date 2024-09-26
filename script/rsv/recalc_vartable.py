#!/usr/bin/env python3
import argparse

##Create a variant table from a vcf and the corresponding reference fasta and gff3
#v0.0.1

def list_intersect(containedl, containingl, mode, threshold):
    if mode == 'match':
        return [ x for x in containedl if x in containingl ]
    elif mode == 'missing':
        return [ x for x in containedl if x not in containingl ]

def count_commented(file):
    lines = open(file, 'r').read().replace('\r\n','\n').rstrip('\n').split('\n')
    count = 0
    for line in lines:
        if line[0] == "#":
            count += 1
    return count

def list_flatten(toflat_list):
    toflat_list = [item for ilist in toflat_list for item in ilist]
    return toflat_list

def list_bgapos(bga):
    bgapos = []
    for i in range(len(bga)):
        if bga[i][0] not in [ x[0] for x in bgapos ]:
            bgapos.append([bga[i][0], []])
            loop_index = [ x[0] for x in bgapos ].index(bga[i][0])
        else:
            loop_index = [ x[0] for x in bgapos ].index(bga[i][0])
        bgapos[loop_index][1].append([int(bga[i][3]) for x in range(int(bga[i][1]),int(bga[i][2]))])
    return bgapos

def list_bedpos(bed):
    bedpos = []
    for i in range(len(bed)):
        if bed[i][0] not in [ x[0] for x in bedpos ]:
            bedpos.append([bed[i][0], []])
            loop_index = [ x[0] for x in bedpos ].index(bed[i][0])
        else:
            loop_index = [ x[0] for x in bedpos ].index(bed[i][0])
        bedpos[loop_index][1].append([int(x) for x in range(int(bed[i][1]),int(bed[i][2]))])
    return bedpos

def vcf_parse(vcf, fields):
    if len(vcf) == 0:
        return []
    fields_header = [ x.split('=')[0] for x in vcf[0].split('\t')[7].split(';') if '=' in x ]
    fields_index = [ fields_header.index(x) for x in fields.split(',') ]
    var = [ [x.split('\t')[0], int(x.split('\t')[1])-1, x.split('\t')[3], x.split('\t')[4]] for x in vcf]
    for i in range(len(var)):
        for j in range(len(fields_index)):
            var[i].append(vcf[i].split('\t')[7].split(';')[fields_index[j]].split('=')[1])
        if var[i][3][0] == '-':
            var_temp = var[i][2]
            var[i][2] = var[i][2] + var[i][3][1:]
            var[i][3] = var_temp
        if var[i][3][0] == '+':
            var[i][3] = var[i][2] + var[i][3][1:]
    return var

def getseq(positions, refseq):
    '''Create a list of nucleotides with complemented base when negative'''
    seq = []
    for pos in positions:
        if pos < 0:
            seq.append(complement[refseq[abs(pos)]])
        else:
            seq.append(refseq[pos])
    return seq

parser = argparse.ArgumentParser(description='Count the number of minor variants in a target vcf reported as major variant in a reference vcf')
debugmode = parser.add_mutually_exclusive_group()
debugmode.add_argument('-v', '--verbose', action='store_true')
debugmode.add_argument('-q', '--quiet', action='store_true')
parser.add_argument('--version', action='version', version='0.0.1')
parser.add_argument('-i', '--input', help='the base')
parser.add_argument('-t', '--target', help='the base')
parser.add_argument('-d', '--depth', help='the base')
parser.add_argument('--min_depth', type=int, default=100, help='the base')
parser.add_argument('-o', '--output', help='output file', default='./')

if __name__ == '__main__':
    args = parser.parse_args()

inp = [x.split('\t') for x in open(args.input, 'r').read().replace('\r\n','\n').rstrip('\n').split('\n')]

bga = [x.split('\t') for x in open(args.depth, 'r').read().replace('\r\n','\n').rstrip('\n').split('\n')]

depth = [ [x[0], list_flatten(x[1])] for x in list_bgapos(bga) ]
depth_chrom = [y[0] for y in depth]

target_index = inp[0].index(args.target)

for i in range(1,len(inp)):
    chrom = inp[i][0].split('|')[0] + '|' + inp[i][0].split('|')[1]
    pos = int(inp[i][0].split('|')[3])-1
    if int(inp[i][target_index]) < 50 and depth[depth_chrom.index(chrom)][1][pos]<args.min_depth:
        inp[i][target_index] = 'NA'

w = open(args.output, 'w')

for i in range(len(inp)):
    w.write('\t'.join(inp[i]) + '\n')

w.close()
