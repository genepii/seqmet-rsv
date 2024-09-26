#!/usr/bin/env python3
import argparse
import copy

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
    var = [ [x.split('\t')[0], int(x.split('\t')[1]), x.split('\t')[3], x.split('\t')[4]] for x in vcf]
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
parser.add_argument('-t', '--target', help='the base')
parser.add_argument('-r', '--reference', help='the base')
parser.add_argument('-g', '--gff', help='the base')
parser.add_argument('-m', '--mode', help='the base')
parser.add_argument('-o', '--output', help='output file', default='./')

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

tvcf = open(args.target, 'r').read().replace('\r\n','\n').rstrip('\n').split('\n')[count_commented(args.target):]

#chrom, pos, ref, alt, af, ao, dp, type
tvar = vcf_parse(tvcf, 'AF,AO,DP,TYPE')

ichrom = 0
ipos = 1
iref = 2
ivar = 3
iaf = 4
icount = 5
idepth = 6
itype = 7

ref = [ [x.split('\n')[0], ''.join(x.split('\n')[1:])] for x in open(args.reference, 'r').read().replace('\r\n', '\n').rstrip('\n').split('>')[1:]]
con = copy.deepcopy(ref)
qsp = copy.deepcopy(con)

for var in tvar:
    refindex = [ x[0] for x in ref ].index(var[ichrom])
    if var[itype] == 'snp' and len(var[ivar]) == 1 and float(var[iaf]) >= 0.5:
        con[refindex][1] = con[refindex][1][0:int(var[ipos])-1] + var[ivar] + con[refindex][1][int(var[ipos]):]
    elif float(var[iaf]) >= 0.5:
        con[refindex][1] = con[refindex][1][0:int(var[ipos])-1] + 'X' + con[refindex][1][int(var[ipos]):]

for var in tvar:
    refindex = [ x[0] for x in ref ].index(var[ichrom])
    if var[itype] == 'snp' and len(var[ivar]) == 1 and float(var[iaf]) < 0.5:
        qsp[refindex][1] = qsp[refindex][1][0:int(var[ipos])-1] + var[ivar] + qsp[refindex][1][int(var[ipos]):]
    elif float(var[iaf]) < 0.5:
        qsp[refindex][1] = qsp[refindex][1][0:int(var[ipos])-1] + 'X' + qsp[refindex][1][int(var[ipos]):]

gff_pnname = []
gff_ref = []
gff_strand = []
gff_position = []

refpn_seq = []
refpn_chromname = []
refpn_pnname = []
refpn_pos = []
refpn_fna = []
refpn_faa = []

conpn_seq = []
conpn_chromname = []
conpn_pnname = []
conpn_pos = []
conpn_fna = []
conpn_faa = []

qsppn_seq = []
qsppn_chromname = []
qsppn_pnname = []
qsppn_pos = []
qsppn_fna = []
qsppn_faa = []

for item in open(args.gff, 'r').read().replace("\r", "").rstrip('\n').split('\n'):
    if item.split('\t')[2] == 'CDS':
        pn = item.split('\t')
        pnname_index = [ x.split('=')[0] for x in pn[8].split(';') ].index("gene")
        pnname = [ x.split('=')[1] for x in pn[8].split(';') ][pnname_index]

        if pnname not in gff_pnname:
            gff_pnname.append(pnname)
            gff_ref.append(pn[0])
            gff_strand.append(pn[6])
            gff_position.append([])
        
        if gff_strand[gff_pnname.index(pnname)] == '+':
            gff_position[gff_pnname.index(pnname)] += [ x for x in range(int(pn[3])-1,int(pn[4])) ]
        elif gff_strand[gff_pnname.index(pnname)] == '-':
            gff_position[gff_pnname.index(pnname)] += [ -x for x in range(int(pn[4])-1,int(pn[3])-2,-1) ]


for i in range(len(gff_pnname)):
    refpn_seq.append(ref[[x[0] for x in ref].index(gff_ref[i])][1])
    conpn_seq.append(con[[x[0] for x in con].index(gff_ref[i])][1])
    qsppn_seq.append(qsp[[x[0] for x in qsp].index(gff_ref[i])][1])
    refpn_chromname += [gff_ref[i]]
    refpn_pnname += [gff_pnname[i]]
    refpn_pos.append([ abs(x) for x in gff_position[i]])
    refpn_fna.append(getseq(gff_position[i], refpn_seq[-1]))
    refpn_faa.append([cdt[''.join(refpn_fna[-1][x:x+3])] if ''.join(refpn_fna[-1][x:x+3]) in cdt else 'X' for x in range(0,len(refpn_fna[-1]),3)])
    conpn_chromname += [gff_ref[i]]
    conpn_pnname += [gff_pnname[i]]
    conpn_pos.append([ abs(x) for x in gff_position[i]])
    conpn_fna.append(getseq(gff_position[i], conpn_seq[-1]))
    conpn_faa.append([cdt[''.join(conpn_fna[-1][x:x+3])] if ''.join(conpn_fna[-1][x:x+3]) in cdt else 'X' for x in range(0,len(conpn_fna[-1]),3)])
    qsppn_chromname += [gff_ref[i]]
    qsppn_pnname += [gff_pnname[i]]
    qsppn_pos.append([ abs(x) for x in gff_position[i]])
    qsppn_fna.append(getseq(gff_position[i], qsppn_seq[-1]))
    qsppn_faa.append([cdt[''.join(qsppn_fna[-1][x:x+3])] if ''.join(qsppn_fna[-1][x:x+3]) in cdt else 'X' for x in range(0,len(qsppn_fna[-1]),3)])

w = open(args.output, 'w')

for var in tvar:
    varoup = []
    refindex = [ x[0] for x in ref ].index(var[ichrom])
    pnname = []
    pnname = [ x for x in refpn_pnname if var[ichrom] == refpn_chromname[refpn_pnname.index(x)] and (int(var[ipos])-1) in refpn_pos[refpn_pnname.index(x)] ]
    varoup.append(var[ichrom].split('|')[0])
    varoup.append(var[ichrom].split('|')[1])
    varoup.append('CODING') if len(pnname) > 0 else varoup.append('NON_CODING')
    varoup.append(var[ipos])
    varoup.append(','.join(pnname)) if len(pnname) > 0 else varoup.append('__')
    varoup.append(','.join([ str(int(((refpn_pos[refpn_pnname.index(x)].index(int(var[ipos])-1))+3)/3)) for x in pnname ])) if len(pnname) > 0 else varoup.append('__')
    if var[itype] == 'snp':
        if int(float(var[iaf])*100) >= 50:
            varoup.append(','.join([ refpn_faa[refpn_pnname.index(x)][int(((refpn_pos[refpn_pnname.index(x)].index(int(var[ipos])-1))+3)/3)-1] for x in pnname ]) + '>' + ','.join([ conpn_faa[conpn_pnname.index(x)][int(((conpn_pos[conpn_pnname.index(x)].index(int(var[ipos])-1))+3)/3)-1] for x in pnname ])) if len(pnname) > 0 else varoup.append('__')
            varoup.append( ref[refindex][1][int(var[ipos])-1] + '>' + var[ivar] )
        elif int(float(var[iaf])*100) < 50:
            varoup.append(','.join([ refpn_faa[refpn_pnname.index(x)][int(((refpn_pos[refpn_pnname.index(x)].index(int(var[ipos])-1))+3)/3)-1] for x in pnname ]) + '>' + ','.join([ qsppn_faa[qsppn_pnname.index(x)][int(((qsppn_pos[qsppn_pnname.index(x)].index(int(var[ipos])-1))+3)/3)-1] for x in pnname ])) if len(pnname) > 0 else varoup.append('__')
            varoup.append( ref[refindex][1][int(var[ipos])-1] + '>' + qsp[refindex][1][int(var[ipos])-1] )
    elif var[itype] == 'ins':
        varoup.append(','.join([ 'INS' for x in pnname ])) if len(pnname) > 0 else varoup.append('INS')
        varoup.append( ref[refindex][1][int(var[ipos])-1] + '>' + var[ivar] )
    elif var[itype] == 'del':
        varoup.append(','.join([ 'DEL' for x in pnname ])) if len(pnname) > 0 else varoup.append('DEL')
        varoup.append( ref[refindex][1][int(var[ipos])-1] + var[iref][1:] + '>' + var[ivar] )
    elif var[itype] == 'mnp':
        varoup.append(','.join([ 'MNP' for x in pnname ])) if len(pnname) > 0 else varoup.append('MNP')
        varoup.append( var[iref] + '>' + var[ivar] )
    else:
        varoup.append(','.join([ 'COMPLEX' for x in pnname ])) if len(pnname) > 0 else varoup.append('COMPLEX')
        varoup.append( var[iref] + '>' + var[ivar] )
    varoup.append('MAJOR') if int(float(var[iaf])*100) >= 50 else varoup.append('MINOR')
    val = ''
    if args.mode == 'freq':
        val = str(int(float(var[iaf])*100))
    elif args.mode == 'abs':
        val = '1'
    elif args.mode == 'norm':
        val = str(int(float(var[iaf])*100) * int(var[idepth]))
    w.write('|'.join([str(x) for x in varoup]) + '\t' + val + '\n')

w.close()
