##Retrieve region of sequences delimited by a forward and reverse pattern
#v0.0.1

f = [ [x.split('\n')[0], ''.join(x.split('\n')[1:])] for x in open('gg-13-8-99.fasta', 'r').read().rstrip('\n').split('>') ][1:]
fw = ['AGAGTTTGATCMTGGCTCAG']
rv = ['GWATTACCGCGGCKGCTG']
mtol = 1
endstr = 3
w = open('ggmb.fa', 'w')
nl = []

print('Forward:' + fw[0])
print('Reverse:' + rv[0])
print('Mismatches tolerated:' + str(mtol))
print('3\' end stringency length:' + str(endstr))

def degenerate(primer):
	r = []
	for it in primer:
		p = list(it)
		for i in range(len(p)):
			if p[i] == 'R':
				r.append(''.join([p[x] if x!=i else 'A' for x in range(len(p))]))
				r.append(''.join([p[x] if x!=i else 'G' for x in range(len(p))]))
				break
			if p[i] == 'Y':
				r.append(''.join([p[x] if x!=i else 'C' for x in range(len(p))]))
				r.append(''.join([p[x] if x!=i else 'T' for x in range(len(p))]))
				break
			if p[i] == 'S':
				r.append(''.join([p[x] if x!=i else 'G' for x in range(len(p))]))
				r.append(''.join([p[x] if x!=i else 'C' for x in range(len(p))]))
				break
			if p[i] == 'W':
				r.append(''.join([p[x] if x!=i else 'A' for x in range(len(p))]))
				r.append(''.join([p[x] if x!=i else 'T' for x in range(len(p))]))
				break
			if p[i] == 'K':
				r.append(''.join([p[x] if x!=i else 'G' for x in range(len(p))]))
				r.append(''.join([p[x] if x!=i else 'T' for x in range(len(p))]))
				break
			if p[i] == 'M':
				r.append(''.join([p[x] if x!=i else 'A' for x in range(len(p))]))
				r.append(''.join([p[x] if x!=i else 'C' for x in range(len(p))]))
				break
			if p[i] == 'B':
				r.append(''.join([p[x] if x!=i else 'C' for x in range(len(p))]))
				r.append(''.join([p[x] if x!=i else 'G' for x in range(len(p))]))
				r.append(''.join([p[x] if x!=i else 'T' for x in range(len(p))]))
				break
			if p[i] == 'D':
				r.append(''.join([p[x] if x!=i else 'A' for x in range(len(p))]))
				r.append(''.join([p[x] if x!=i else 'G' for x in range(len(p))]))
				r.append(''.join([p[x] if x!=i else 'T' for x in range(len(p))]))
				break
			if p[i] == 'H':
				r.append(''.join([p[x] if x!=i else 'A' for x in range(len(p))]))
				r.append(''.join([p[x] if x!=i else 'C' for x in range(len(p))]))
				r.append(''.join([p[x] if x!=i else 'T' for x in range(len(p))]))
				break
			if p[i] == 'V':
				r.append(''.join([p[x] if x!=i else 'A' for x in range(len(p))]))
				r.append(''.join([p[x] if x!=i else 'C' for x in range(len(p))]))
				r.append(''.join([p[x] if x!=i else 'G' for x in range(len(p))]))
				break
			if p[i] == 'N':
				r.append(''.join([p[x] if x!=i else 'A' for x in range(len(p))]))
				r.append(''.join([p[x] if x!=i else 'C' for x in range(len(p))]))
				r.append(''.join([p[x] if x!=i else 'G' for x in range(len(p))]))
				r.append(''.join([p[x] if x!=i else 'T' for x in range(len(p))]))
				break
	return r

def mismatch(primer):
	r = []
	for it in primer:
		p = list(it)
		for i in range(len(p) - endstr):
			r.append(''.join([p[x] if x!=i else 'A' for x in range(len(p))]))
			r.append(''.join([p[x] if x!=i else 'C' for x in range(len(p))]))
			r.append(''.join([p[x] if x!=i else 'G' for x in range(len(p))]))
			r.append(''.join([p[x] if x!=i else 'T' for x in range(len(p))]))
	return r
a = ['A','C','G','T']
fwndeg = len([x for x in list(fw[0])if x not in a])
rvndeg = len([x for x in list(rv[0])if x not in a])

for i in range(fwndeg):
	fw = degenerate(fw)

for i in range(rvndeg):
	rv = degenerate(rv)

for i in range(mtol):
	fw = mismatch(fw)
	rv = mismatch(rv)

rv = [ x[::-1].replace('A','t').replace('C','g').replace('G','c').replace('T','a').upper() for x in rv ]

for it in f:
	rt = ''
	cut = 0
	for i in range(len(it[1])-len(fw[0])):
		if it[1][i:i+len(fw[0])] in fw:
			rt = it[1][i:]
			cut += 0
			break
	for i in range(len(rt)-len(rv[0])):
		if rt[i:i+len(rv[0])] in rv:
			rt = rt[:i+len(rv[0])]
			cut += 2
			break
	if cut == 2:
		nl.append([it[0], rt])

lnl = [len(x[1]) for x in nl]

print('Matches:' + str(len(nl))+ '/' + str(len(f)))
print('Min:' + str(min(lnl)))
print('Max:' + str(max(lnl)))
print('Mean:' + str(sum(lnl) / len(nl)))
print('Matches shorter than 300:' + str(len([x for x in lnl if x < 300])))
print('Matches longer than 580:' + str(len([x for x in lnl if x > 580])))
w.write('\n'.join(['>' + '\n'.join(x) for x in nl]))
w.close()