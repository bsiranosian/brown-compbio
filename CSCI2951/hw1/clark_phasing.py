# Haplotype phasing with Clark's method 
# Ben Siranosian | CSCI2951
import argparse
import os
import itertools
import random

# argparse stuff here
# need to get input file of genotypes, output file for solution, posssibly some parameters about how to run
os.chdir('C://Users/Admin/Documents/GitHub/brown-compbio/CSCI2951/hw1/')
inFile='AACB2014-dataset1-genotypes.txt'

with open(inFile,'r') as inf:
	lines = [map(int, list(l.strip())) for l in inf.readlines()]

# our vector of genotypes by removing duplicates
# note that this is a choice and could not be done in the furture
# there's something useful about keeping the numbers of things
# could do some sort of ML approach where previously seen is used 
# to resolve ambiguities.
lines.sort()
gen = list(l for l,_ in itertools.groupby(lines))

# ensure we get a new result each time 
# disable for testing
#random.shuffle(g)

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
	for genotype in g:
		# check for homozygotes
		if 2 not in genotype:
			g.remove(genotype)
			# check if genotype is in the resolved list already
			# as said before, this duplicates method can be switched around
			if genotype not in r:
				r.append(genotype)
	return [g,r]

# resolveHeterozygotes(genotypes, resolved): finds heterozygotes with 
# a single ambiguous site, resolves them and adds to the resolved list
# returns [updated genotypes, resolved list]
def resolveHeterozygotes(g, r):
	for genotype in g:
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
			if o1 not in r: r.append(o1)
			if o2 not in r: r.append(o2)
			# remove from genotype list
			g.remove(genotype)
	return[g,r]

# inferenceRule(genotypes, resolved): applys the inference Rule once
# to the list of genotypes. If a genotype can be resolved by the resolved list, 
# append the new haplotype to resolved and remove from genotypes
def inferenceRule(g, r):
	counter=0
	for genotype in g:
		notExplained=True
		for possible in r:
			if canExplain(genotype,possible):
				g.remove(genotype)
				r.append(resolveGenotype(genotype,possible))
				counter +=1
				break
		else:
			continue
		break
	print 'One iteration resolved ' +str(counter) +' genotypes'
	print str(len(g)) + ' genotypes remain'
	return [g,r]

#schematic for doing inference
res1 = resolveHomozygotes(gen)
g1=res1[0]
r1=res1[1]
for l in g1: print l
print
for l in r1: print l
res2=resolveHeterozygotes(g1, r1)
g2=res2[0]
r2=res2[1]
for l in g2: print l
print
for l in r2: print l
res3=inferenceRule(g2,r2)
g3=res3[0]
r3=res3[1]
for l in g3: print l
print
for l in r3: print l

for i in range(20):
	res4=inferenceRule(g3,r3)
	g3=res4[0]
	r3=res4[1]
	for l in g3: print l
	print
	for l in r3: print l
	