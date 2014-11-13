# DOMINATING SET and SET COVER algorithms 
# Ben Siranosian | CSCI2951N | hw3 | 2014-11-07
# This file will contain the functions necessary to compute the dominating set and set cover
# for a given input

import itertools
import os
import numpy as np
import scipy.stats
import math
import random

os.chdir('Documents/GitHub/brown-compbio/CSCI2951/hw3')

# dominatingSet(adjMat): computes the smallest dominating set for the adjacency matrix given. 
# The graph to compute dominating set on should be represented by an adjacency matrix, size |V|*|V|
# where v1,v2 is 1 if there is an edge between v1 and v2 in the original graph, and 0 otherwise. 
# Diagonal of the adjacency matrix should be all 1
# adjacency matrix should be defined as a list of lists. 
def dominatingSet(adjMat):
	n = len(adjMat)
	# check all lengths of subsets, starting with the smallest to speed things up. 
	for m in range(1,len(adjMat)+1):
		subsets = itertools.combinations(range(n),m)
		#print "subsets of size: " + str(m)
		for s in subsets:
			if dominating(adjMat, s):
				#print "Found a dominating set of size " + str(m)
				return s 

# dominating(adjMat, s): returns true if a subset sub is dominating for a graph defined by 
# adjacency matrix adjMat, false otherwise
def dominating(adjMat, sub):
	# select colums in adjMat from s
	cols = []
	for a in adjMat:
		row = []
		for s in sub:
			row.append(a[s])
		cols.append(row)

	# check if rowsums are 1
	return sum([sum(c)>0 for c in cols])==len(adjMat)

# testing dominatingSet
# Example presented in class
ex1 = [[1,1,0,0,1,0],
	   [1,1,1,0,0,0],
	   [0,1,1,1,0,1],
	   [0,0,1,1,0,0],
	   [1,0,0,0,1,1],
	   [0,0,1,0,1,1]]
# dominatingSet(ex1)
# solution [0,2] is valid!

# setCover(subsets): Given a set of subsets defined as a list of lists, return the subsets
# that represent a set cover of the union of all sets.
def setCover(subsets):
	# find set of universe 
	universe = []
	for s in subsets:
		for i in s:
			if i not in universe: universe.append(i)

	# solve this through conversion to dominating set problem
	# make an adjacency matrix 
	totalSize = len(subsets)+ len(universe)
	adj = [[0 for i in range(totalSize)] for j in range(totalSize)]

	# make diagonal 1
	for i in range(totalSize):
		adj[i][i] =1

	# Fill in as desired  - first the clique for conections between index sets.
	for i in range(len(subsets)):
		for j in range(len(subsets)):
			adj[i][j] =1

	# Fill in connections between subsets and index sets
	for i in range(len(subsets)):
		for j in range(len(universe)):
			if universe[j] in subsets[i]:
				# connection between subset i and element j
				adj[i][j+len(subsets)] = 1
				adj[j+len(subsets)][i] = 1

	for a in adj:
		print a 

	return dominatingSet(adj)	

# example from wikipedia
ex2 = [[1,2,3],
	   [1,2],
	   [2,3,4],
	   [3,4,5]]
# setCover(ex2)
# solution [0,3] is valid!

# hapMapToGenotype(inFile, numLines): Reads a file of genotypes for a single chromosome 
# downloaded from HapMap. Converts to a matrix where rows are SNP sites and columns are individuals.
# If the genotype data contains Ns, that SNP is discarded. 
# only reads numLines of the input file to reduce the size of the output. 
# Reads information in second column to get information on the two alleles at a given site
# returns [n*m matrix of genotypes, each entry is 0 / 1 for homozygote or 2 for heterozygote, 
# list of length n of SNP names and major/minor alleles, list of length,
# list of length m of individual names
# if discardFixed is True, doesn't keep alleles with all heterozygotes
def hapMapToGenotype(inFile, numLines, discardFixed=True):
	# make sure we have the file
	assert os.path.isfile(inFile)

	with open(inFile, 'r') as inf:
		lines = []
		# read first line of headers
		firstLine = inf.readline().strip().split(' ')
		# read in numLines of the input
		for l in range(numLines):
			lines.append(inf.readline().strip().split(' '))

	# get names of individuals - fist 9 items are header info
	individuals = firstLine[9:]

	# names of SNPs | numpy array so we can subset later
	SNPnames = np.array([l[0] for l in lines])

	# alleles at each position
	alleles = np.array([l[1].split('/') for l in lines])
	# data for each row
	calls = np.array([l[9:] for l in lines])
	# Keep rows that dont have any Ns
	noN = np.array([''.join(d).count('N')==0 for d in calls])

	# subset all data to noN
	SNPnames = SNPnames[noN]
	alleles = alleles[noN]
	calls = calls[noN]

	# rows that are all one symbol 
	allHomo = np.array([''.join(d).count(d[0][0])==(len(d)*2) for d in calls])
	# if discardedFixed, only keep SNPs that aren't all homozygotes
	if discardFixed:
		SNPnames = SNPnames[~allHomo]
		alleles = alleles[~allHomo]
		calls = calls[~allHomo]

	# print alleles
	# print SNPnames
	# print calls
	#convert calls to genotype 0/1/2
	genotypes = np.zeros((len(SNPnames), len(individuals)))
	for row in range(len(SNPnames)):
		for col in range(len(individuals)):
			# switch genotpyes to major/minor
			#print calls[row]
			#print ''.join(calls[row]).count(alleles[row][1]) 
			#print ''.join(calls[row]).count(alleles[row][0])
			if ''.join(calls[row]).count(alleles[row][1]) > ''.join(calls[row]).count(alleles[row][0]):
				temp = alleles[row][1]
				alleles[row][1] = alleles[row][0]
				alleles[row][0] = temp
			genotypes[row, col] = stringToGenotype(calls[row,col],alleles[row])

	# return all the data we got!
	return [genotypes, SNPnames, individuals]

# stringToGenotype(call, alleles): helper function to call 0/1/2 genotypes
def stringToGenotype(call, alleles): 
	if call[0]==alleles[0] and call[1]==alleles[0]:
		return 2
	elif call[0]==alleles[1] and call[1]==alleles[1]:
		return 1
	else:
		return 0

# calculateR2(genotypes): calculates pairwise r^2 between alleles, given the output
# of hapMapToGenotype, an n*m matrix of 0/1/2 genotype calls
def calculateR2(genotypes):
	# make pairwise array
	r2array = np.zeros((np.shape(genotypes)[0],np.shape(genotypes)[0]))

	# go over each element - pairwise interaction
	for i in range(0, len(r2array)):
		for j in range(i+1, len(r2array)):
			# calculate r^2 for the interaction between i and j
			# I dont have phased data! 
			# use an approximation from Rogers and Huff (2009) where r = correlation(Y,Z)

			# get data for row
			gi = genotypes[i] 
			gj = genotypes[j] 

			# correlation between genotpye data
			# if two vectors are the same, correlation value is 0
			if np.array_equal(gi,gj):
				cor=0
			else:
				cor = scipy.stats.pearsonr(gi,gj)[0]
				if np.isnan(cor): cor =0
			# square this!
			r2array[i,j] = math.pow(cor,2)
			r2array[j,i] = math.pow(cor,2)

	# return correlations
	return r2array


# r2ToAdjacency(r2array, threshold): creates an adjacency matrix from a given r2array.
# Nodes are individual SNPs, edges exist between two alleles if their r2 value is > threshold
# if zeroDiagonal, doesn't add one to the diagonal of the matrix.
def r2ToAdjacency(r2array, threshold, zeroDiagonal=False):
	# make values greater than threshold 1, 0 otherwise
	toReturn = (r2array > threshold) * 1
	# set diagonal to 1
	if not zeroDiagonal:
		for i in range(len(toReturn)):
			toReturn[i,i] =1

	return toReturn

# code to produce the figures I used
h = hapMapToGenotype('genotypes_chr1_subset.txt', 3386)
r= calculateR2(h[0])
import matplotlib.pyplot as plt
adj = r2ToAdjacency(r, 0.15)

#figure 1
plt.matshow(r); plt.colorbar(); plt.show()

#figure 2a
plt.matshow(r2ToAdjacency(r, 0.75), cmap='Blues'); plt.show()
#figure 2b
plt.matshow(r2ToAdjacency(r, 0.50), cmap='Blues'); plt.show()
#figure 2c
plt.matshow(r2ToAdjacency(r, 0.15), cmap='Blues'); plt.show()

# testing relationship between informativeness and r2 in LD-select
reps = 100
r2values = [float(f)/100 for f in range(100)]
sizes = []
for val in r2values:
	print val
	adj = adj = r2ToAdjacency(r, val) 
	# for a given r2 value
	setSize=[]
	for i in range(reps):
		#print "Replicate " + str(i)
		# pick a random subset of 10 adjacent SNPs
		start = random.choice(range(600))
		end = start + 13
		sub = adj[start:end, start:end]
		if np.sum(sub)==13:
			setSize.append(13)
		else:
			dom = dominatingSet(sub)
		setSize.append(len(dom))
	sizes.append(setSize)

# figure 3a 	
plt.plot(r2values, [np.mean(s)/13.0 for s in sizes])
plt.xlabel('r^2')
plt.ylabel('Mean informativeness')
plt.show()
# figure 3b
plt.plot(r2values, [s.count(13) for s in sizes])
plt.xlabel('r^2')
plt.ylabel('Sets of size 13')
plt.show()