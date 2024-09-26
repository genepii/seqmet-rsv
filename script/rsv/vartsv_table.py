#!/usr/bin/env python2.7
import os
import sys
import getopt
import math
import copy

##Create a table of annotated variants from a vcf and gff3 file
#v0.1.0

def main(argv):
    global vcf
    global ref
    global gff
    global oup
    global datatype
    vcf = ''
    ref = ''
    gff = ''
    oup = ''
    datatype = 'freq'
    try:
        opts, args = getopt.getopt(argv, 'hv:r:g:o:d:', ['--help', '--vcf', '--reference', '--gff', '--output', '--datatype'])
        for opt, arg in opts:
            if opt in ('-h', '--help'):
                usage()
                sys.exit()
            elif opt in ('-v', '--vcf'):
                vcf = arg
            elif opt in ('-r', '--reference'):
                ref = arg
            elif opt in ('-g', '--gff'):
                gff = arg
            elif opt in ('-o', '--output'):
                oup = arg
            elif opt in ('-d', '--datatype'):
                datatype = arg
        if vcf == '':
            usage()
            sys.exit()
        if ref == '':
            usage()
            sys.exit()
        if gff == '':
            usage()
            sys.exit()
        if oup == '':
            oup = vcf.split('/')[-1].split('.')[0]
    except getopt.GetoptError:
        usage()
        sys.exit(2)

def usage():
    print 'usage: ', sys.argv[0], '-h --help -v --vcf -r --ref [fasta] -g --gff [gff3] -o --output [table] -d --datatype [freq,abs,norm]'

if __name__ == '__main__':
    main(sys.argv[1:])

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

bases = ['A', 'C', 'G', 'T']

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
    field_index = headers.index('AF')
    return data[field_index]

def icount(line):
    info = line[iinfo]
    headers = [ x.split('=')[0] for x in info.split(';') ]
    data = [ x.split('=')[-1] for x in info.split(';') ]
    field_index = headers.index('AO')
    return data[field_index]

def idepth(line):
    info = line[iinfo]
    headers = [ x.split('=')[0] for x in info.split(';') ]
    data = [ x.split('=')[-1] for x in info.split(';') ]
    field_index = headers.index('DP')
    return data[field_index]

def itype(line):
    info = line[iinfo]
    headers = [ x.split('=')[0] for x in info.split(';') ]
    data = [ x.split('=')[-1] for x in info.split(';') ]
    field_index = headers.index('TYPE')
    return data[field_index]

headers = open(vcf, 'r').read().replace('\r\n', '\n').rstrip('\n').split('\n')[count_commented(vcf)-1].split('\t')
ichrom = headers.index("#CHROM")
ipos = headers.index("POS")
iref = headers.index("REF")
ivar = headers.index("ALT")
iinfo = headers.index("INFO")
iformat = headers.index("FORMAT")

sampleid = vcf.split('/')[-1].split('.')[0]

ref = [ [x.split('\n')[0], ''.join(x.split('\n')[1:])] for x in open(ref, 'r').read().replace("\r", "").rstrip('\n').split('>')[1:]]
con = copy.deepcopy(ref)
vcf_file = open(vcf, 'r').read().replace('\r\n', '\n').rstrip('\n').split('\n')[count_commented(vcf):]

for item in vcf_file:
    var = item.split('\t')
    refindex = [ x[0] for x in ref ].index(var[ichrom])
    if itype(var) == 'snp' and float(iaf(var)) >= 0.5:
        con[refindex][1] = con[refindex][1][0:int(var[ipos])-1] + var[ivar] + con[refindex][1][int(var[ipos]):]
    elif float(iaf(var)) >= 0.5:
        con[refindex][1] = con[refindex][1][0:int(var[ipos])-1] + 'X' + con[refindex][1][int(var[ipos]):]

qsp = copy.deepcopy(con)

for item in vcf_file:
    var = item.split('\t')
    refindex = [ x[0] for x in ref ].index(var[ichrom])
    if itype(var) == 'snp' and float(iaf(var)) < 0.5:
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
qsppn_posframe = []

for item in open(gff, 'r').read().replace("\r", "").rstrip('\n').split('\n')[count_commented(gff):]:
    if item.split('\t')[2] == 'CDS':
        pn = item.split('\t')
        if "gene" in [ x.split('=')[0] for x in pn[8].split(';') ]:
            pnname_index = [ x.split('=')[0] for x in pn[8].split(';') ].index("gene")
        elif "Parent" in [ x.split('=')[0] for x in pn[8].split(';') ]:
            pnname_index = [ x.split('=')[0] for x in pn[8].split(';') ].index("Parent")
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
    refpn_nestedlist_loop = [ [refpn_fna[-1][x]+refpn_fna[-1][x+1]+refpn_fna[-1][x+2], refpn_fna[-1][x]+refpn_fna[-1][x+1]+refpn_fna[-1][x+2], refpn_fna[-1][x]+refpn_fna[-1][x+1]+refpn_fna[-1][x+2]] if 'X' not in [refpn_fna[-1][x], refpn_fna[-1][x+1], refpn_fna[-1][x+2]] else ['XXX', 'XXX', 'XXX'] for x in range(0,len(refpn_fna[-1]),3)]
    refpn_faa.append([item for sublist in refpn_nestedlist_loop for item in sublist])
    conpn_chromname += [gff_ref[i]]
    conpn_pnname += [gff_pnname[i]]
    conpn_pos.append([ abs(x) for x in gff_position[i]])
    conpn_fna.append(getseq(gff_position[i], conpn_seq[-1]))
    conpn_nestedlist_loop = [ [conpn_fna[-1][x]+conpn_fna[-1][x+1]+conpn_fna[-1][x+2], conpn_fna[-1][x]+conpn_fna[-1][x+1]+conpn_fna[-1][x+2], conpn_fna[-1][x]+conpn_fna[-1][x+1]+conpn_fna[-1][x+2]] if 'X' not in [conpn_fna[-1][x], conpn_fna[-1][x+1], conpn_fna[-1][x+2]] else ['XXX', 'XXX', 'XXX'] for x in range(0,len(conpn_fna[-1]),3)]
    conpn_faa.append([item for sublist in conpn_nestedlist_loop for item in sublist])
    qsppn_chromname += [gff_ref[i]]
    qsppn_pnname += [gff_pnname[i]]
    qsppn_pos.append([ abs(x) for x in gff_position[i]])
    qsppn_nestedlist_loop = [ [0, 1, 2] for x in range(0,len(qsppn_pos[-1]),3)]
    qsppn_posframe.append([item for sublist in qsppn_nestedlist_loop for item in sublist])

w = open(oup, 'w')

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
    if itype(var) == 'snp':
        if int(float(iaf(var))*100) >= 50:
            refpn_loop = []
            conpn_loop = []

            if len(pnname) == 0:
                varoup.append('__')
                varoup.append( ref[refindex][1][int(var[ipos])-1] + '>' + var[ivar] )

            for i in range(len(pnname)):
                x = pnname[i]
                pos_loop = refpn_pos[refpn_pnname.index(x)].index(int(var[ipos])-1)
                pnname_loop = refpn_pnname.index(x)
                if refpn_faa[pnname_loop][pos_loop] in cdt:
                    refpn_loop.append(cdt[refpn_faa[pnname_loop][pos_loop]])
                else:
                    refpn_loop.append('X')
                if conpn_faa[pnname_loop][pos_loop] in cdt:
                    conpn_loop.append(cdt[conpn_faa[pnname_loop][pos_loop]])
                else:
                    conpn_loop.append('X')

            if len(pnname) > 0:
                varoup.append(','.join(refpn_loop) + '>' + ','.join(conpn_loop))
                varoup.append( ref[refindex][1][int(var[ipos])-1] + '>' + var[ivar] )
            
        elif int(float(iaf(var))*100) < 50:

            refpn_loop = []
            qsppn_loop = []

            if len(pnname) == 0:
                varoup.append('__')
                varoup.append( ref[refindex][1][int(var[ipos])-1] + '>' + var[ivar] )

            for i in range(len(pnname)):
                x = pnname[i]
                pos_loop = refpn_pos[refpn_pnname.index(x)].index(int(var[ipos])-1)
                if (int(var[ipos])-3) in refpn_pos[refpn_pnname.index(x)]:
                    pos_loop_2l = refpn_pos[refpn_pnname.index(x)].index(int(var[ipos])-3)
                else:
                    pos_loop_2l = 0
                if (int(var[ipos])-2) in refpn_pos[refpn_pnname.index(x)]:
                    pos_loop_1l = refpn_pos[refpn_pnname.index(x)].index(int(var[ipos])-2)
                else:
                    pos_loop_1l = 0
                if (int(var[ipos])) in refpn_pos[refpn_pnname.index(x)]:
                    pos_loop_1 = refpn_pos[refpn_pnname.index(x)].index(int(var[ipos]))
                else:
                    pos_loop_1 = 0
                if (int(var[ipos])+1) in refpn_pos[refpn_pnname.index(x)]:
                    pos_loop_2 = refpn_pos[refpn_pnname.index(x)].index(int(var[ipos])+1)
                else:
                    pos_loop_2 = 0

                pnname_loop = refpn_pnname.index(x)
                if refpn_faa[pnname_loop][pos_loop] in cdt:
                    refpn_loop.append(cdt[refpn_faa[pnname_loop][pos_loop]])
                else:
                    refpn_loop.append('X')
                if qsppn_posframe[pnname_loop][pos_loop] == 0 and (var[ivar] + conpn_fna[pnname_loop][pos_loop_1] + conpn_fna[pnname_loop][pos_loop_2]) in cdt:
                    qsppn_loop.append(cdt[var[ivar] + conpn_fna[pnname_loop][pos_loop_1] + conpn_fna[pnname_loop][pos_loop_2]])
                elif qsppn_posframe[pnname_loop][pos_loop] == 1 and (conpn_fna[pnname_loop][pos_loop_1l] + var[ivar] + conpn_fna[pnname_loop][pos_loop_1]) in cdt:
                    qsppn_loop.append(cdt[conpn_fna[pnname_loop][pos_loop_1l] + var[ivar] + conpn_fna[pnname_loop][pos_loop_1]])
                elif qsppn_posframe[pnname_loop][pos_loop] == 2 and (conpn_fna[pnname_loop][pos_loop_2l] + conpn_fna[pnname_loop][pos_loop_1l] + var[ivar]) in cdt:
                    qsppn_loop.append(cdt[conpn_fna[pnname_loop][pos_loop_2l] + conpn_fna[pnname_loop][pos_loop_1l] + var[ivar]])
                else:
                    qsppn_loop.append('X')
            
            if len(pnname) > 0:
                varoup.append(','.join(refpn_loop) + '>' + ','.join(qsppn_loop))
                varoup.append( ref[refindex][1][int(var[ipos])-1] + '>' + var[ivar] )

    elif itype(var) == 'ins':
        varoup.append(','.join([ 'INS' for x in pnname ])) if len(pnname) > 0 else varoup.append('INS')
        varoup.append( ref[refindex][1][int(var[ipos])-1] + '>' + var[ivar] )
    elif itype(var) == 'del':
        varoup.append(','.join([ 'DEL' for x in pnname ])) if len(pnname) > 0 else varoup.append('DEL')
        varoup.append( ref[refindex][1][int(var[ipos])-1] + var[iref][1:] + '>' + var[ivar] )
    elif itype(var) == 'mnp':
        if int(float(iaf(var))*100) >= 50:
            varoup.append(','.join([ 'MNP' for x in pnname ])) if len(pnname) > 0 else varoup.append('MNP')
            varoup.append( ''.join([ ref[refindex][1][int(var[ipos])-1+x] for x in range(len(var[ivar])) ]) + '>' + var[ivar] )
        elif int(float(iaf(var))*100) < 50:
            varoup.append(','.join([ 'MNP' for x in pnname ])) if len(pnname) > 0 else varoup.append('MNP')
            varoup.append( ''.join([ ref[refindex][1][int(var[ipos])-1+x] for x in range(len(var[ivar])) ]) + '>' + var[ivar] )
    else:
        varoup.append(','.join([ 'COMPLEX' for x in pnname ])) if len(pnname) > 0 else varoup.append('COMPLEX')
        varoup.append( ref[refindex][1][int(var[ipos])-1] + var[iref][1:] + '>' + var[ivar] )
    if int(float(iaf(var))*100) >= 50:
        varoup.append('MAJOR')
    elif 'REVERSED=t' in var[iinfo]:
        varoup.append('REFASMINOR')
    else:
        varoup.append('MINOR')
    val = ''
    if datatype == 'freq':
        val = str(int(float(iaf(var))*100))
    elif datatype == 'abs':
        val = '1'
    elif datatype == 'norm':
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
w = open("refpn.fna", "w")
for i in range(len(refpn_fna)):
    w.write('>' + sampleid + '|' + refpn_pnname[i] + '\n' + ''.join(refpn_fna[i]) + '\n')
w.close()

w.close()
w = open("conpn.fna", "w")
for i in range(len(conpn_fna)):
    w.write('>' + sampleid + '|' + conpn_pnname[i] + '\n' + ''.join(conpn_fna[i]) + '\n')
w.close()
