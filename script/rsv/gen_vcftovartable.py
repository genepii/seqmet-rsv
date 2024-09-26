#!/usr/bin/env python3
import argparse
import math
import copy

##Create a variant table from a vcf file
#v0.0.4

parser = argparse.ArgumentParser(description='Create a variant table from a vcf file')
debugmode = parser.add_mutually_exclusive_group()
debugmode.add_argument('--verbose', action='store_true')
debugmode.add_argument('--quiet', action='store_true')
parser.add_argument('--version', action='version', version='0.0.4')
parser.add_argument('-v', '--vcf', help='vcf')
parser.add_argument('-r', '--reference', help='fasta')
parser.add_argument('-g', '--gff', help='gff')
parser.add_argument('-d', '--datatype', help='freq,abs,norm', default='freq')
parser.add_argument('-o', '--output', help='output file', default='output.tsv')

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

def getseq(positions, refseq):
    '''Create a list of nucleotides with complemented base when negative'''
    seq = []
    for pos in positions:
        if pos < 0:
            seq.append(complement[refseq[abs(pos)]])
        else:
            seq.append(refseq[pos])
    return seq

def count_commented(file):
    lines = open(file, 'r').read().rstrip('\n').split('\n')
    count = 0
    for line in lines:
        if line[0] == "#":
            count += 1
    return count

def iaf(line):
    info = line[iinfo]
    headers = [ x.split('=')[0] for x in info.split(';') ]
    data = [ x.split('=')[-1] for x in info.split(';') ]
    field_index = headers.index('VAF')
    return data[field_index]

def idepth(line):
    info = line[iinfo]
    headers = [ x.split('=')[0] for x in info.split(';') ]
    data = [ x.split('=')[-1] for x in info.split(';') ]
    field_index = headers.index('DP')
    return data[field_index]

headers = open(args.vcf, 'r').read().replace('\r\n', '\n').rstrip('\n').split('\n')[count_commented(args.vcf)-1].split('\t')
ichrom = headers.index("#CHROM")
ipos = headers.index("POS")
iref = headers.index("REF")
ivar = headers.index("ALT")
iinfo = headers.index("INFO")
iformat = headers.index("FORMAT")

sampleid = args.vcf.split('/')[-1].split('.')[0]

ref = [ [x.split('\n')[0], ''.join(x.split('\n')[1:])] for x in open(args.reference, 'r').read().replace("\r", "").rstrip('\n').split('>')[1:]]
con = copy.deepcopy(ref)
vcf_file = open(args.vcf, 'r').read().replace('\r\n', '\n').rstrip('\n').split('\n')[count_commented(args.vcf):]

for item in vcf_file:
    var = item.split('\t')
    refindex = [ x[0] for x in ref ].index(var[ichrom])
    if len(var[iref]) == 1 and len(var[iref]) == len(var[ivar]) and len(var[ivar]) == 1 and float(iaf(var)) >= 0.5:
        con[refindex][1] = con[refindex][1][0:int(var[ipos])-1] + var[ivar] + con[refindex][1][int(var[ipos]):]
    elif float(iaf(var)) >= 0.5:
        con[refindex][1] = con[refindex][1][0:int(var[ipos])-1] + 'X' + con[refindex][1][int(var[ipos]):]

qsp = copy.deepcopy(con)

for item in vcf_file:
    var = item.split('\t')
    refindex = [ x[0] for x in ref ].index(var[ichrom])
    if len(var[iref]) == 1 and len(var[iref]) == len(var[ivar]) and len(var[ivar]) == 1 and float(iaf(var)) < 0.5:
        qsp[refindex][1] = qsp[refindex][1][0:int(var[ipos])-1] + var[ivar] + qsp[refindex][1][int(var[ipos]):]
    elif float(iaf(var)) < 0.5:
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
        pnname = [ x.split('=')[-1] for x in pn[8].split(';') ][pnname_index]

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

for item in vcf_file:
    var = item.split('\t')
    varoup = []
    refindex = [ x[0] for x in ref ].index(var[ichrom])
    pnname = []
    pnname = [ x for x in refpn_pnname if var[ichrom] == refpn_chromname[refpn_pnname.index(x)] and (int(var[ipos])-1) in refpn_pos[refpn_pnname.index(x)] ]
    varoup.append(var[ichrom].split('|')[0])
    varoup.append(var[ichrom].split('|')[1])
    varoup.append('CODING') if len(pnname) > 0 else varoup.append('NON_CODING')
    varoup.append(var[ipos])
    varoup.append(','.join(pnname)) if len(pnname) > 0 else varoup.append('__')
    varoup.append(','.join([ str((((refpn_pos[refpn_pnname.index(x)].index(int(var[ipos])-1))+3)/3)) for x in pnname ])) if len(pnname) > 0 else varoup.append('__')
    if len(var[iref]) == 1 and len(var[iref]) == len(var[ivar]):
        print(var)
        print(refpn_faa)
        print(refpn_pos)
        print(refpn_pnname)
        print(pnname)
        print(conpn_faa)
        print(conpn_pos)
        print(conpn_pnname)
        print(var[ipos])
        if int(float(iaf(var))*100) >= 50:
            varoup.append(','.join([ refpn_faa[refpn_pnname.index(x)][(((refpn_pos[refpn_pnname.index(x)].index(int(var[ipos])-1))+3)/3)-1] for x in pnname ]) + '>' + ','.join([ conpn_faa[conpn_pnname.index(x)][(((conpn_pos[conpn_pnname.index(x)].index(int(var[ipos])-1))+3)/3)-1] for x in pnname ])) if len(pnname) > 0 else varoup.append('__')
            varoup.append( ref[refindex][1][int(var[ipos])-1] + '>' + var[ivar] )
        elif int(float(iaf(var))*100) < 50:
            varoup.append(','.join([ refpn_faa[refpn_pnname.index(x)][(((refpn_pos[refpn_pnname.index(x)].index(int(var[ipos])-1))+3)/3)-1] for x in pnname ]) + '>' + ','.join([ qsppn_faa[qsppn_pnname.index(x)][(((qsppn_pos[qsppn_pnname.index(x)].index(int(var[ipos])-1))+3)/3)-1] for x in pnname ])) if len(pnname) > 0 else varoup.append('__')
            varoup.append( ref[refindex][1][int(var[ipos])-1] + '>' + qsp[refindex][1][int(var[ipos])-1] )
    elif len(var[iref]) < len(var[ivar]):
        varoup.append(','.join([ 'INS' for x in pnname ])) if len(pnname) > 0 else varoup.append('INS')
        varoup.append( ref[refindex][1][int(var[ipos])-1] + '>' + var[ivar] )
    elif len(var[iref]) > len(var[ivar]):
        varoup.append(','.join([ 'DEL' for x in pnname ])) if len(pnname) > 0 else varoup.append('DEL')
        varoup.append( ref[refindex][1][int(var[ipos])-1] + var[iref][1:] + '>' + var[ivar] )
    elif len(var[iref]) > 1 and len(var[iref]) == len(var[ivar]):
        if int(float(iaf(var))*100) >= 50:
            varoup.append(','.join([ 'MNP' for x in pnname ])) if len(pnname) > 0 else varoup.append('MNP')
            varoup.append( ''.join([ ref[refindex][1][int(var[ipos])-1+x] for x in range(len(var[ivar])) ]) + '>' + var[ivar] )
        elif int(float(iaf(var))*100) < 50:
            varoup.append(','.join([ 'MNP' for x in pnname ])) if len(pnname) > 0 else varoup.append('MNP')
            varoup.append( ''.join([ ref[refindex][1][int(var[ipos])-1+x] for x in range(len(var[ivar])) ]) + '>' + qsp[refindex][1][int(var[ipos])-1] )
    else:
        varoup.append(','.join([ 'COMPLEX' for x in pnname ])) if len(pnname) > 0 else varoup.append('COMPLEX')
        varoup.append( ref[refindex][1][int(var[ipos])-1] + var[iref][1:] + '>' + var[ivar] )
    varoup.append('MAJOR') if int(float(iaf(var))*100) >= 50 else varoup.append('MINOR')
    val = ''
    if args.datatype == 'freq':
        val = str(int(float(iaf(var))*100))
    elif args.datatype == 'abs':
        val = '1'
    elif args.datatype == 'norm':
        val = str(int(float(iaf(var))*100) * int(idepth(var)))
    w.write('|'.join(varoup) + '\t' + val + '\n')

w.close()
w = open("refpn.faa", "w")
for i in range(len(refpn_faa)):
    w.write('>' + sampleid + '|' + refpn_pnname[i] + '\n' + ''.join(refpn_faa[i]) + '\n')
w.close()

w.close()
w = open("conpn.faa", "w")
for i in range(len(conpn_faa)):
    w.write('>' + sampleid + '|' + conpn_pnname[i] + '\n' + ''.join(conpn_faa[i]) + '\n')
w.close()

w.close()
w = open("qsppn.faa", "w")
for i in range(len(qsppn_faa)):
    w.write('>' + sampleid + '|' + qsppn_pnname[i] + '\n' + ''.join(qsppn_faa[i]) + '\n')
w.close()

w.close()
w = open("refpn.fna", "w")
for i in range(len(refpn_fna)):
    w.write('>' + sampleid + '|' + refpn_pnname[i] + '\n' + ''.join(refpn_fna[i]) + '\n')
w.close()

w.close()
w = open("conpn.fna", "w")
for i in range(len(conpn_fna)):
    w.write('>' + sampleid + '|' + conpn_pnname[i] + '\n' + ''.join(conpn_fna[i]) + '\n')
w.close()

w.close()
w = open("qsppn.fna", "w")
for i in range(len(qsppn_fna)):
    w.write('>' + sampleid + '|' + qsppn_pnname[i] + '\n' + ''.join(qsppn_fna[i]) + '\n')
w.close()
