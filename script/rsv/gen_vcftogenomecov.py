#!/usr/bin/env python3
import argparse
import copy

##Create a genomecov bed file from a comprehensive vcf file
#v0.0.1

def count_commented(file):
    lines = open(file, 'r').read().replace('\r\n','\n').rstrip('\n').split('\n')
    count = 0
    for line in lines:
        if line[0] == "#":
            count += 1
    return count

def vcf_parse(vcf, fields):
    if len(vcf) == 0:
        return []
    var = []
    newvcf = []
    for i in range(len(vcf)):
        if [vcf[i].split('\t')[0], int(vcf[i].split('\t')[1])-1] not in var:
            var.append([vcf[i].split('\t')[0], int(vcf[i].split('\t')[1])-1])
            newvcf.append(vcf[i])
    for i in range(len(var)):
        fields_header = [ x.split('=')[0] for x in newvcf[i].split('\t')[7].split(';') ]
        fields_index = [ fields_header.index(x) for x in fields.split(',') ]
        for j in range(len(fields_index)):
            var[i].append(newvcf[i].split('\t')[7].split(';')[fields_index[j]].split('=')[1])
    return var

parser = argparse.ArgumentParser(description='Create a genomecov bed file from a comprehensive vcf file')
debugmode = parser.add_mutually_exclusive_group()
debugmode.add_argument('--verbose', action='store_true')
debugmode.add_argument('--quiet', action='store_true')
parser.add_argument('--version', action='version', version='0.0.1')
parser.add_argument('-v', '--vcf', help='vcf')
parser.add_argument('-r', '--reference', help='fasta')
parser.add_argument('-m', '--mode', help='bed,bga', default='bga')
parser.add_argument('-o', '--output', help='output file', default='output.bed')

if __name__ == '__main__':
    args = parser.parse_args()

tvcf = open(args.vcf, 'r').read().replace('\r\n','\n').rstrip('\n').split('\n')[count_commented(args.vcf):]
ref = [ [x.split('\n')[0], ''.join(x.split('\n')[1:])] for x in open(args.reference, 'r').read().replace('\r\n', '\n').rstrip('\n').split('>')[1:]]

#chrom, pos + INFOFIELD
vcf = vcf_parse(tvcf, 'DP')

chrom_list = [ x[0] for x in ref ]
bedpos = [ [ 0 for x in range(len(y[1])) ] for y in ref ]

for i in range(len(vcf)):
    chrom = chrom_list.index(vcf[i][0])
    pos = vcf[i][1]
    bedpos[chrom][pos] = int(vcf[i][2])

w = open(args.output, 'w')

if args.mode == "bga":
    for j in range(len(bedpos)):
        i = 0
        reflength = len(bedpos[j])
        while i < reflength:
            curpos = copy.deepcopy(i)
            curposdp = bedpos[j][i]
            endpos = -1
            while endpos < 0:
                if i+1 == reflength:
                    endpos = i
                elif curposdp != bedpos[j][i+1]:
                    endpos = i
                i += 1
            w.write(chrom_list[j] + '\t' + str(curpos) + '\t' + str(endpos+1) + '\t' + str(curposdp) + '\n')

elif args.mode == "bed":
    for j in range(len(bedpos)):
        for i in range(len(bedpos[j])):
            curposdp = bedpos[j][i]
            w.write(chrom_list[j] + '\t' + str(i) + '\t' + str(i+1) + '\t' + str(curposdp) + '\n')

w.close()
