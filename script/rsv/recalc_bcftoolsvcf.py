#!/usr/bin/env python3
import argparse
import copy

##Add reference as a minor variant on a new line and add custom tags on a bcftools vcf
#v0.0.1

def count_commented(file):
    lines = open(file, 'r').read().replace('\r\n','\n').rstrip('\n').split('\n')
    count = 0
    for line in lines:
        if line[0] == "#":
            count += 1
    return count

def gen_type(ref, alt):
    if len(ref) == len(alt) and len(ref) == 1:
        return 'snp'
    elif len(ref) == len(alt):
        return 'mnp'
    elif len(ref) > len(alt) and len(alt) == 1:
        return 'del'
    elif len(ref) < len(alt) and len(ref) == 1:
        return 'ins'
    else:
        return 'complex'

parser = argparse.ArgumentParser(description='Add reference as a minor variant on a new line and add custom tags on a bcftools vcf')
debugmode = parser.add_mutually_exclusive_group()
debugmode.add_argument('--verbose', action='store_true')
debugmode.add_argument('--quiet', action='store_true')
parser.add_argument('--version', action='version', version='0.0.1')
parser.add_argument('-v', '--vcf', help='vcf')
parser.add_argument('-o', '--output', help='output file', default='output.vcf')

if __name__ == '__main__':
    args = parser.parse_args()

tvcf = open(args.vcf, 'r').read().replace('\r\n','\n').rstrip('\n').split('\n')[count_commented(args.vcf):]

vcf = [ [x.split('\t')[0], x.split('\t')[1], x.split('\t')[2], x.split('\t')[3], x.split('\t')[4], x.split('\t')[5], x.split('\t')[6], x.split('\t')[7].split(';'), x.split('\t')[8].split(':'), x.split('\t')[9].split(':')] for x in tvcf]

vcf_newentry = []
vcf_newentry_simple = []

w = open(args.output, 'w')

comment = open(args.vcf, 'r').read().replace('\r\n','\n').rstrip('\n').split('\n')[:count_commented(args.vcf)]
w.write('\n'.join(comment[:-1]) + '\n')
w.write('##recalc_bcftoolsvcf="Reference as minor variant, INFO/VAF and INFO/TYPE added"\n')
w.write('##INFO=<ID=TYPE,Number=A,Type=String,Description="The type of allele, either snp, mnp, ins, del, or complex.">\n')
w.write('##INFO=<ID=AF,Number=A,Type=Float,Description="The fraction of reads with alternate allele (nALT/nSumAll)">\n')
w.write('##INFO=<ID=REVERSED,Number=A,Type=String,Description="Indicates that the reference is found as a minor variant.">\n')
w.write(comment[-1] + '\n')

##INFO=<ID=AD,Number=R,Type=Integer,Description="Total allelic depths (high-quality bases)">

for i in range(len(vcf)):

    info_header = [ x.split('=')[0] for x in vcf[i][7] ]
    format_header = vcf[i][8]
    
    if float(vcf[i][7][info_header.index('RAF')].split('=')[1]) == 0.5:
        vcf[i][7][info_header.index('RAF')] = 'RAF=0.499'

    if float(vcf[i][7][info_header.index('RAF')].split('=')[1]) > 0 and float(vcf[i][7][info_header.index('RAF')].split('=')[1]) < 0.5 and vcf[i][0:2] not in vcf_newentry_simple:
        
        vcf_newentry_simple.append(vcf[i][0:2])
        newentry = copy.deepcopy(vcf[i])
        
        newentry[7][info_header.index('ADF')] = 'ADF=' + newentry[7][info_header.index('ADF')].split('=')[1].split(',')[0] + ',' + newentry[7][info_header.index('ADF')].split('=')[1].split(',')[0]
        newentry[7][info_header.index('ADR')] = 'ADR=' + newentry[7][info_header.index('ADR')].split('=')[1].split(',')[0] + ',' + newentry[7][info_header.index('ADR')].split('=')[1].split(',')[0]
        newentry[7][info_header.index('AD')] = 'AD=' + newentry[7][info_header.index('AD')].split('=')[1].split(',')[0] + ',' + newentry[7][info_header.index('AD')].split('=')[1].split(',')[0]
        newentry[7][info_header.index('DP4')] = 'DP4=' + newentry[7][info_header.index('DP4')].split('=')[1].split(',')[0] + ',' + newentry[7][info_header.index('DP4')].split('=')[1].split(',')[1] + ',' + newentry[7][info_header.index('DP4')].split('=')[1].split(',')[0] + ',' + newentry[7][info_header.index('DP4')].split('=')[1].split(',')[1]
        
        if 'RPBZ' in info_header:
            newentry[7][info_header.index('RPBZ')] = 'RPBZ=' + str(-float(newentry[7][info_header.index('RPBZ')].split('=')[1]))
        if 'MQBZ' in info_header:
            newentry[7][info_header.index('MQBZ')] = 'MQBZ=' + str(-float(newentry[7][info_header.index('MQBZ')].split('=')[1]))
        if 'BQBZ' in info_header:
            newentry[7][info_header.index('BQBZ')] = 'BQBZ=' + str(-float(newentry[7][info_header.index('BQBZ')].split('=')[1]))
        if 'MQSBZ' in info_header:
            newentry[7][info_header.index('MQSBZ')] = 'MQSBZ=' + str(-float(newentry[7][info_header.index('MQSBZ')].split('=')[1]))
        if 'SCBZ' in info_header:
            newentry[7][info_header.index('SCBZ')] = 'SCBZ=' + str(-float(newentry[7][info_header.index('SCBZ')].split('=')[1]))

        if float(newentry[7][info_header.index('RAF')].split('=')[1]) >= 0.5:
            newentry[7][info_header.index('AC')] = 'AC=1'
            newentry[7][info_header.index('AN')] = 'AN=1'
        else:
            newentry[7][info_header.index('AC')] = 'AC=0'
            newentry[7][info_header.index('AN')] = 'AN=0'
            newentry[9][format_header.index('GT')] = '.'
            newentry[9][format_header.index('PL')] = '255,255'
        
        newentry[9][format_header.index('ADF')] = newentry[9][format_header.index('ADF')].split(',')[0] + ',' + newentry[9][format_header.index('ADF')].split(',')[0]
        newentry[9][format_header.index('ADR')] = newentry[9][format_header.index('ADR')].split(',')[0] + ',' + newentry[9][format_header.index('ADR')].split(',')[0]
        newentry[9][format_header.index('AD')] = newentry[9][format_header.index('AD')].split(',')[0] + ',' + newentry[9][format_header.index('AD')].split(',')[0]
        newentry[9][format_header.index('VAF')] = newentry[7][info_header.index('RAF')].split('=')[1]
        newentry[9][format_header.index('VAF1')] = newentry[7][info_header.index('RAF')].split('=')[1]

        w.write('\t'.join(newentry[0:3]) + '\t' + newentry[4] + '\t' + newentry[3] + '\t' + '\t'.join(newentry[5:7]) + '\tTYPE=' + gen_type(newentry[4],newentry[3]) + ';AF=' + newentry[9][format_header.index('VAF')] + ';REVERSED=t;' + ';'.join(newentry[7]) + '\t' + ':'.join(newentry[8]) + '\t' + ':'.join(newentry[9]) + '\n')

    if 'VAF' in format_header:
        w.write('\t'.join(vcf[i][0:7]) + '\tTYPE=' + gen_type(vcf[i][3],vcf[i][4]) + ';AF=' + vcf[i][9][format_header.index('VAF')] + ';' + ';'.join(vcf[i][7]) + '\t' + ':'.join(vcf[i][8]) + '\t' + ':'.join(vcf[i][9]) + '\n')
    else:
        w.write('\t'.join(vcf[i][0:7]) + '\tTYPE=NA;AF=0;' + ';'.join(vcf[i][7]) + '\t' + ':'.join(vcf[i][8]) + '\t' + ':'.join(vcf[i][9]) + '\n')

w.close()
