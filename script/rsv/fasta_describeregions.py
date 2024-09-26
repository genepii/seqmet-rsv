#!/usr/bin/env python3
import argparse
import copy

##
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

def getseq(positions, refseq):
    '''Create a list of nucleotides with complemented base when negative'''
    seq = []
    for pos in positions:
        if pos < 0:
            seq.append(complement[refseq[abs(pos)]])
        else:
            seq.append(refseq[pos])
    return seq

parser = argparse.ArgumentParser(description='')
debugmode = parser.add_mutually_exclusive_group()
debugmode.add_argument('-v', '--verbose', action='store_true')
debugmode.add_argument('-q', '--quiet', action='store_true')
parser.add_argument('--version', action='version', version='0.0.1')
parser.add_argument('-r', '--reference', help='the base')
parser.add_argument('-f', '--fasta', help='the base')
parser.add_argument('-g', '--gff', help='the base')
parser.add_argument('-m', '--mode', help='the base')
parser.add_argument('-o', '--output', help='output file', default='./')
parser.add_argument('-x', '--faoutput', help='output file', default='./')

if __name__ == '__main__':
    args = parser.parse_args()

stoplist = ['TAA', 'TGA', 'TAG']

cdt = {'TTT': 'F', 'CTT': 'L', 'ATT': 'I', 'GTT': 'V',
       'TTC': 'F', 'CTC': 'L', 'ATC': 'I', 'GTC': 'V',
       'TTA': 'L', 'CTA': 'L', 'ATA': 'I', 'GTA': 'V',
       'TTG': 'L', 'CTG': 'L', 'ATG': 'M', 'GTG': 'V',
       'TCT': 'S', 'CCT': 'P', 'ACT': 'T', 'GCT': 'A',
       'TCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A',
       'TCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A',
       'TCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A',
       'TAT': 'Y', 'CAT': 'H', 'AAT': 'N', 'GAT': 'D',
       'TAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D',
       'TAA': '*', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E',
       'TAG': '*', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E',
       'TGT': 'C', 'CGT': 'R', 'AGT': 'S', 'GGT': 'G',
       'TGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G',
       'TGA': '*', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G',
       'TGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'}

complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'X': 'X', 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H', 'V': 'B', 'H': 'D'}

fasta = [ [x.split('\n')[0], ''.join(x.split('\n')[1:])] for x in open(args.fasta, 'r').read().replace('\r\n', '\n').rstrip('\n').split('>')[1:]]

msapos_nocov = [ 0 for x in range(len(fasta[0][1])) ]
msapos_gap = [ 0 for x in range(len(fasta[0][1])) ]

reflength = len(fasta[0][1])
msapos_ref = [ "" for x in range(len(fasta[0][1])) ]

msapos_gff = [ "" for x in range(len(fasta[0][1])) ]

for i in range(len(fasta)):
    for j in range(len(fasta[0][1])):
        if fasta[i][1][j] == '-':
            msapos_gap[j]+=1
        if fasta[i][1][j] == 'N':
            msapos_nocov[j]+=1

msapos_range = []
range_current = ""
range_current_type = "todo"
for i in range(len(msapos_gap)):
    if i == 0:
        pass
    elif range_current_type == "todo" and (msapos_gap[i] > 0 or msapos_nocov[i] > 0):
        pass
    elif range_current_type == "ok" and msapos_gap[i] == 0 and msapos_nocov[i] == 0:
        pass
    else:
        msapos_range.append(range_current)
        range_current = ""
        
    if msapos_gap[i] > 0:
        range_current += "-"
    elif msapos_nocov[i] > 0:
        range_current += "N"
    else:
        range_current += "X"
    
    if range_current[-1] in ["-", "N"]:
        range_current_type = "todo"
    else:
        range_current_type = "ok"
    
    if i == len(msapos_gap)-1:
        msapos_range.append(range_current)
        range_current = ""

x=1
msapos_range_info = []
for i in range(len(msapos_range)):
    startpos = x
    endpos = x + len(msapos_range[i]) - 1
    if msapos_range[i][0] in ["-", "N"]:
        cmd = "realign --maxiterate 100 --nwildcard --anysymbol --op 5.0 --leavegappyregion --auto"
    else:
        cmd = "preserve"
    seq = msapos_range[i]
    msapos_range_info.append([startpos, endpos, cmd, seq])
    x = x + len(msapos_range[i])

loop_todo = len(msapos_range_info)
i = 0
while i < loop_todo:
    if i < 2:
        i+=1
        continue
    if msapos_range_info[i][0] - msapos_range_info[i-2][1] < 100 and msapos_range_info[i][2] == msapos_range_info[i-2][2] != "preserve" and msapos_range_info[i-1][2] == "preserve" and msapos_range_info[i][2] != "preserve":
        msapos_range_info[i-2][1] = msapos_range_info[i][1]
        msapos_range_info[i-2][3] = msapos_range_info[i-2][1] + msapos_range_info[i-1][1] + msapos_range_info[i][1]
        msapos_range_info.remove(msapos_range_info[i])
        msapos_range_info.remove(msapos_range_info[i-1])
        i=0
        loop_todo = len(msapos_range_info)
    else:
        i+=1

for i in range(len(fasta)):
    for j in range(len(msapos_range_info)):
        #print(len(fasta[i][1]))
        #print(msapos_range_info[j][0]-1)
        #print(msapos_range_info[j][1]-1)
        #print(fasta[i][1][msapos_range_info[j][0]-1:msapos_range_info[j][1]-1])
        #print(msapos_range_info[j])
        if fasta[i][1][msapos_range_info[j][0]-1] in ["-", "N"] and fasta[i][1][msapos_range_info[j][1]-1] in ["-", "N"] and fasta[i][1][msapos_range_info[j][0]-1:msapos_range_info[j][1]].count('N')+fasta[i][1][msapos_range_info[j][0]-1:msapos_range_info[j][1]].count('-') > 0.05*len(fasta[i][1][msapos_range_info[j][0]-1:msapos_range_info[j][1]]):
            #print(len(fasta[i][1]))
            fasta[i][1] = fasta[i][1][0:msapos_range_info[j][0]-1] + 'N'*len(fasta[i][1][msapos_range_info[j][0]-1:msapos_range_info[j][1]]) + fasta[i][1][msapos_range_info[j][1]:]
            #print(len(fasta[i][1]))

w = open(args.output, 'w')
for i in range(len(msapos_range_info)):
    w.write('\t'.join([ str(x) for x in msapos_range_info[i]][0:3]) + '\n')
w.close()

w = open(args.faoutput, 'w')
for i in range(len(fasta)):
    w.write('>' + fasta[i][0] + '\n' + fasta[i][1] + '\n')
w.close()
