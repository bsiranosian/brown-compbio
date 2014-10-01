# Haplotype phasing with Clark's method 
# Ben Siranosian | CSCI2951
import os
import itertools
import random
import sys

def main():
	# get arguments 
	if len(sys.argv) != 2:
		sys.exit('USAGE: python clark_phasing.py genotype_file \n Takes as input a list of genotypes (defined by 0,1,2), one per line. \n Computes the clark phasing of the input with rules defined in the accompanying document. \n By default runs 1000 iterations of the algorithm and picks the best solution \n (fewest orphans and fewest explaining haplotypes).')
	inFile=sys.argv[1]
	# need to get input file of genotypes, output file for solution, posssibly some parameters about how to run
	# os.chdir('C://Users/Admin/Documents/GitHub/brown-compbio/CSCI2951/hw1/')
	# inFile='AACB2014-dataset1-genotypes.txt'

	with open(inFile,'r') as inf:
		lines = [map(int, list(l.strip())) for l in inf.readlines()]

	# our vector of genotypes by removing duplicates
	# note that this is a choice and could not be done in the furture
	# there's something useful about keeping the numbers of things
	# could do some sort of ML approach where previously seen is used 
	# to resolve ambiguities.
	lines.sort()
	gen = list(l for l,_ in itertools.groupby(lines))

	answer = returnBest(gen, 1000)
	print "Number of orphans remaining: " + str(len(answer[1]))
	if len(answer[1]) > 0:
		print "Orphan genotypes are"
		for a in answer[1]:
			print '\t' + str(a)
	print "Number of haplotypes explaining resolved genotypes: " + str(len(answer[0]))
	print "Haplotypes explaining genotypes"
	for a in answer[0]:
		print '\t' + str(a)

## HELPER FUNCTIONS
# canExplain(genotype, haplotype): returns True if a haplotype can 
# explain part of a genotype, False otherwise
def canExplain(g, h): 
	for i in range(len(g)):
		if g[i] == 0 or g[i] ==1:
			if g[i] != h[i]:
				return False
	# if we make it through to the end without failing
	# h can explain g.  
	return True

# resolveGenotype(ambiguous, resolved): resolves an ambiguous genotype given the chosen 
# resoled haplotype. Returns the explaination for the genotype opposite of resolved.
def resolveGenotype(a, r):
	o = []
	# iterate over positions in vector
	for i, j in zip(a,r):
		# if same, oppsite is already defined
		if i==j:
			o.append(i)
		# if i is ambiguous, take opposite of j
		if i==2 and j==0:
			o.append(1)
		if i==2 and j==1:
			o.append(0)
	return o 

## STARTING THE ALGORITHM ##
# resolveHomozygotes(genotypes, resolved=None): finds homozygote genotpyes,
# removes them from the genotype set and adds to resolved set. 
# if resovled==None (the case in starting the algorithm) makes a new resolved list 
# returns [updated genotypes, resolved list]
def resolveHomozygotes(g, r=None):
	# make a resolved list if we need to 
	if r==None:
		r=[]
	g2 = g[:]
	r2 = r[:]
	for genotype in g2:
		# check for homozygotes
		if 2 not in genotype:
			g2.remove(genotype)
			# check if genotype is in the resolved list already
			# as said before, this duplicates method can be switched around
			if genotype not in r2:
				r2.append(genotype)
	return [g2,r2]

# resolveHeterozygotes(genotypes, resolved): finds heterozygotes with 
# a single ambiguous site, resolves them and adds to the resolved list
# returns [updated genotypes, resolved list]
def resolveHeterozygotes(g, r):
	g2 = g[:]
	r2 = r[:]
	for genotype in g2:
		# look for a single ambiguuous site
		if genotype.count(2) == 1:
			o1=[]
			o2=[]
			for i in genotype:
				if i==2:
					o1.append(0)
					o2.append(1)
				else:
					o1.append(i)
					o2.append(i)
			# append the new unique haplotypes
			if o1 not in r2: r2.append(o1)
			if o2 not in r2: r2.append(o2)
			# remove from genotype list
			g2.remove(genotype)
	return[g2,r2]

# inferenceRule(genotypes, resolved): applys the inference Rule once
# to the list of genotypes. If a genotype can be resolved by the resolved list, 
# append the new haplotype to resolved and remove from genotypes
def inferenceRule(g, r):
	g2 = g[:]
	r2 = r[:]
	counter=0

	for genotype in g2:
		#print genotype
		notExplained=True
		for possible in r2:
			if canExplain(genotype,possible):
				g2.remove(genotype)
				resolved =resolveGenotype(genotype,possible)
				if resolved not in r2:
					r2.append(resolved)
				counter +=1
				break
		else:
			continue
		break
	#print 'One iteration resolved ' +str(counter) +' genotypes'
	#print str(len(g2)) + ' genotypes remain'
	return [g2,r2]

# doPhasing(g): runs the whole phasing algorithm. Takes as input a set of genotypes. 
# return [list of explanations, orphan genotypes]
def doPhasing(g):
	g2 = g[:]
	# resolve homozygotes
	res1 = resolveHomozygotes(g2)
	g2=res1[0]
	r2=res1[1]
	# resolve heterozygotes
	res2=resolveHeterozygotes(g2, r2)
	g2=res2[0]
	r2=res2[1]

	#apply the inference rule until we're stuck with orphans
	oldGLen = None
	while len(g2) != oldGLen:
		oldGLen = len(g2)
		res3=inferenceRule(g2,r2)
		g2=res3[0]
		r2=res3[1]

	# print "haplotypes explaining original genotypes: "
	# for l in r2:
	# 	print l

	# if len(g2) > 0:
	# 	print "orphans remaining: " + str(len(g))
	# 	for l in g2:
	# 		print l
	# return list of haplotypes
	return [r2, g2]

#returnBest(g, n): repeats the Clark Phasing algorithm n number of times on the input genotypes,
# randomizng the list before each iteration. 
def returnBest(g, n):
	g2 = g[:]
	# run n times on randomized data
	resList = []
	for i in range(n):
		random.shuffle(g2)
		resList.append(doPhasing(g2))

	# find best solution 
	# min number of orphans
	minO = min(map(len,[res[1] for res in resList]))
	# print "minO: " + str(minO)
	# subset on this 
	resListSubset1 = [res for res in resList if len(res[1])==minO]

	#min number of explanations
	minE = min(map(len,[res[0] for res in resListSubset1])) 
	# print "minE: " + str(minE)

	# subset on this 
	resListSubset2 = [res for res in resList if len(res[0])==minE]

	#randomize and return
	random.shuffle(resListSubset2)
	return resListSubset2[0]

if __name__ == '__main__':
	main()