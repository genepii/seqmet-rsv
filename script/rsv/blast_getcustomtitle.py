from __future__ import print_function
import os
import sys
import getopt

##Get title of an accession id for blasts results on custom db
#v0.0.1

def main(argv):
    global inp
    global oup
    global db
    inp = ''
    oup = ''
    db = ''
    try:
        opts, args = getopt.getopt(argv, 'hi:o:d:', ['help', 'input', 'output', 'db'])
        for opt, arg in opts:
            if opt in ('-h', '--help'):
                usage()
                sys.exit()
            elif opt in ('-i', '--input'):
                inp = arg
            elif opt in ('-o', '--output'):
                oup = arg
            elif opt in ('-d', '--db'):
                db = arg
        if oup == '':
            usage()
            sys.exit()
    except getopt.GetoptError:
        usage()
        sys.exit(2)

def usage():
    print('usage: ' + sys.argv[0] + ' -h --help -o --output [path]')

if __name__ == '__main__':
    main(sys.argv[1:])

blast = [ x.split('\t') for x in open(inp, 'r').read().replace('\r\n','\n').rstrip('\n').split('\n')]
acc2title = [ x.split('\t') for x in open(db, 'r').read().replace('\r\n','\n').rstrip('\n').split('\n')]
acc = [ x[0] for x in acc2title ]
foundacc_acc = []
foundacc_title = []

w = open(oup, 'a')

for i in range(len(blast)):
    if blast[i][8] in foundacc_acc:
        w.write( '\t'.join(blast[i][0:7]) + '\t' + foundacc_title[foundacc_acc.index(blast[i][8])] + '\t' + '\t'.join(blast[i][8:]) + '\n')
        continue
    elif blast[i][8] in acc:
        foundacc_acc.append(acc2title[acc.index(blast[i][8])][0])
        foundacc_title.append(acc2title[acc.index(blast[i][8])][1])
        w.write( '\t'.join(blast[i][0:7]) + '\t' + foundacc_title[-1] + '\t' + '\t'.join(blast[i][8:]) + '\n')
        continue
    else:
        w.write( '\t'.join(blast[i][0:7]) + '\t' + 'NA' + '\t' + '\t'.join(blast[i][8:]) + '\n')

w.close()
