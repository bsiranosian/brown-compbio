# calclating haplotype frequencies by the EM method
# Ben Siranosian | CSCI2951
import os
import itertools
import random
import sys
import math

def main():
	# get arguments 
	if len(sys.argv) != 2:
		sys.exit('USAGE: python em_haplotype.py genotype_file \n Takes as input a list of genotypes (defined by 0,1,2), one per line. \n Computes the expectation maximization phasing of the input, returning the frequencies of the haplotypes')
	inFile=sys.argv[1]
	# need to get input file of genotypes, output file for solution, posssibly some parameters about how to run
	# os.chdir('C://Users/Admin/Documents/GitHub/brown-compbio/CSCI2951/hw1/')
	# inFile='AACB2014-dataset1-genotypes.txt'

	with open(inFile,'r') as inf:
		# store genotypes as strings
		genotypeList = [l.strip() for l in inf.readlines()]

	# do some basic calcs 
	# genotype to data dictionary [list of ints, number, frequency, probability, {haplotypes that can explain: [prob, list of ints]}]
	G = dict()

	for g in genotypeList:
		if g not in G:
			exHap = explainingHaplotypes(g)
			G[g] =  [map(int, list(g)), 1, 0.0, 0.0, \
			{''.join([str(t) for t in l[0]])+'/'+''.join([str(t) for t in l[1]]):[0, l] for l in exHap}]
		else: 
			G[g][1]+=1

	# set frequencies
	for g in G.keys():
		G[g][2] = len(genotypeList)/float(G[g][1])

	# get list of possible haplotypes
	haps = getHaplotypes([l[0] for l in G.values()])
	# initial frequencies: assume all are equally likely 
	# haplotype data. dictionary haplotype: [probability, list of ints]
	H = {''.join([str(t) for t in l]):[float(1)/len(haps), l] for l in haps}
	
# eStep(gDict, hDict): takes in the dictionaries of genotypes and associated data, haplotypes and assoicated data
# and computes one iteration of the expectation step
def eStep(gDict, hDict):
	# going to calculate the probability of genotype j
	for j in gDict.keys():
		Pj = 0
		# sum over all explanations...
		explanations = gDict[j][4]
		for pair in explainations.keys():
			if explanations[pair][1][0] == explanations[pair][1][1]:
				# add square of frequency to sum
				prob = hDict[''.join([str(t) for t in explanations[pair][1][0]])] * hDict[''.join([str(t) for t in explanations[pair][1][0]])]
				Pj += prob
				gDict[j][4][pair][0]=prob 
			else:
				# add 2 * pk * pl
				prob = 2*hDict[''.join([str(t) for t in explanations[pair][1][0]])] * hDict[''.join([str(t) for t in explanations[pair][1][1]])]
				Pj += prob
				gDict[j][4][pair][0]=prob 

		# total probability for genotype J
		gDict[j][3] = Pj

		# update individual explanation probabilities
		for pair in gDict[j][4].keys():
			gDict[j][4][pair][0] = gDict[j][4][pair][0] / Pj 

	# calculate likelihood
	likelihood = 1
	for j in gDict.keys():
		likelihood = likelihood * math.pow(gDict[j][3],gDict[j][1])
	print "log lik: " + str(math.log(likelihood))

	# return updated data
	return [gDict, hDict]

# mStep(gDict, hDict): takes in the dictionaries of genotypes and associated data, haplotypes and assoicated data
# and computes one iteration of the maximization step
def mStep(gDict, hDict):
	# update haplotype frequencies now!
	for r in hDict.keys():
		#haplotpype as list
		hap  = hDict[r][1]
		prob = 0.5
		# for every genotype
		for j in gDict.keys():
			# genotype as list
			gen = gDict[i][0]
			exProb = 0
			# for every explanation
			for i in gDict[j][4].keys():
				# multiply by 0,1,or 2
				delta = gDict[j][4][i][1].count(hap)
				exProb += delta * gDict[j][4][i][0]
			prob += exProb
		# update r with final prob
		hDict[r][0] = prob

# explainingHaplotypes(genotype): given a genotype as a list of ints,
#nreturns a list of haplotype pairs that can explain that genotype
def explainingHaplotypes(genotype): 
	explainList = []
	# get possible haplotypes for the individual genotype
	pHaps = getHaplotypes([genotype])
	# solve this the naive way so I don't have to rewrite getHaplotypes...
	for hap in pHaps:
		pair = [hap, resolveGenotype(genotype, hap)]
		if (pair not in explainList) and (pair[::-1] not in explainList):
			explainList.append(pair)


# getHaplotypes(gList): returns a list of all possible haplotypes that can explain the list of genotypes. 
# gList must be a list of lists representing genotpes as integers in {0,1,2} 
def getHaplotypes(gList):
	# remove duplicates 
	gList.sort()
	uList = list(l for l,_ in itertools.groupby(gList))
	# haplotype list to add to 
	hList = []

	for u in uList:
		# if no ambiguities 
		if (2 not in u) and (u not in hList):
			hList.append(u)
		# if ambiguities 
		else:
			#temporary list to store
			temp = []
			for pos in u:
				if pos == 2:
					if len(temp) ==0:
						temp.append([0])
						temp.append([1])
					else:
						temp0 = [l+[0] for l in temp]
						temp1 = [l+[1] for l in temp]
						temp = temp0 + temp1
				else:
					if len(temp) ==0:
						temp.append(pos)
					else:
						temp = [l+[pos] for l in temp]

			# add unique elements
			for t in temp:
				if t not in hList: hList.append(t)

	return hList

## CODE REUSED FROM HW1
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