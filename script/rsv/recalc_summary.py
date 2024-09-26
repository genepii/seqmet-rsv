#!/usr/bin/env python3
import argparse

##
#v0.0.1

def list_intersect(containedl, containingl, mode, threshold):
    if mode == 'match':
        return [ x for x in containedl if x in containingl ]
    elif mode == 'missing':
        return [ x for x in containedl if x not in containingl ]

def list_flatten(toflat_list):
    toflat_list = [item for ilist in toflat_list for item in ilist]
    return toflat_list

def count_commented(file):
    lines = open(file, 'r').read().replace('\r\n','\n').rstrip('\n').split('\n')
    count = 0
    for line in lines:
        if line[0] == "#":
            count += 1
    return count

parser = argparse.ArgumentParser(description='')
debugmode = parser.add_mutually_exclusive_group()
debugmode.add_argument('-v', '--verbose', help='the base', action='store_true')
debugmode.add_argument('-q', '--quiet', help='the base', action='store_false')
parser.add_argument('--version', action='version', version='0.0.1')
parser.add_argument('-i', '--input', help='the base')
parser.add_argument('-x', '--header', help='the base', action='store_false')
parser.add_argument('-k', '--keys', help='the base', default='sample_id,fastq_readcount,bam_readcount,run_id')
parser.add_argument('-o', '--output', help='output file', default='./')

if __name__ == '__main__':
    args = parser.parse_args()

tsv = [ x.split('\t') for x in open(args.input, 'r').read().replace('\r\n','\n').rstrip('\n').split('\n')[count_commented(args.input):] ]

if args.header:
    tsv_header = tsv[0]
    tsv_lines = tsv[1:]
    keys_index = [ tsv_header.index(x) for x in args.keys.split(',') ]
else:
    tsv_lines = tsv[0:]
    keys_index = [ int(x) for x in args.keys.split(',') ]

keys_lines = []
for i in range(len(tsv_lines)):
    keys_lines.append([])
    for j in range(len(keys_index)):
        keys_lines[i].append(tsv_lines[i][keys_index[j]])

uniq_lines = [ x if keys_lines.count(x) == 1 else [] for x in keys_lines ]

dup_lines = [ x if keys_lines.count(x) > 1 else [] for x in keys_lines ]

dup_lines_work = []
for i in range(len(dup_lines)):
    dup_lines_work.append([])
    for j in range(len(tsv_header)):
        if j in keys_index:
            dup_lines_work[i].append(tsv_lines[i][j])
        elif tsv_header[j] in ['reference_id']:
            dup_lines_work[i].append('NA')
        else:
            dup_lines_work[i].append('0')

dup_lines_out = [ list(x) for x in set(tuple(x) for x in dup_lines_work if x != []) ]

w = open(args.output, 'w')

if args.header:
    w.write('\t'.join(tsv_header) + "\n")

uniq_lines_sampleid = [ x[0] for x in uniq_lines if x != [] ]

for i in range(len(tsv_lines)):
    if uniq_lines[i] != [] and uniq_lines_sampleid.count(tsv_lines[i][tsv_header.index('sample_id')]) == 1:
        w.write('\t'.join(tsv_lines[i]) + "\n")
    elif uniq_lines[i] != [] and tsv_lines[i][tsv_header.index('bam_readcount')] == 'NA':
        continue
    elif uniq_lines[i] != [] and int(tsv_lines[i][tsv_header.index('bam_readcount')]) > 0:
        w.write('\t'.join(tsv_lines[i]) + "\n")

for i in range(len(dup_lines_out)):
    if dup_lines_out[i][0] not in uniq_lines_sampleid:
        w.write('\t'.join(dup_lines_out[i]) + "\n")

w.close()
