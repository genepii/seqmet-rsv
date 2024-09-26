##Convert a gff file to bed format
#v0.0.1

o = [x.replace(' ', '_').replace(';', '\t').replace('=', '\t').replace('%', '-').split('\t') for x in open('ecoli_O157H7.gff', 'r').read().rstrip('\n').split('\n') if x[0] != '#']
w = open('ecoli_O157H7.bed', 'w')
for i in o:
	if i[2] == 'CDS':
		if 'product' in i:
			w.write(i[0] + '\t' + str(int(i[3])-1) + '\t' + i[4] + '\t' + i[i.index('product')+1] + '\n')
		elif 'protein_id' in i:
			w.write(i[0] + '\t' + str(int(i[3])-1) + '\t' + i[4] + '\t' + i[i.index('protein_id')+1] + '\n')
w.close()
