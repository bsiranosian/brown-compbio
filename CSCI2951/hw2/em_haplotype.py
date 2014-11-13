# calclating haplotype frequencies by the EM method
# Ben Siranosian | CSCI2951
import os
import itertools
import random
import sys
import math

def main():
	# get arguments 
	if (len(sys.argv) < 3) or (len(sys.argv) > 4):
		sys.exit('USAGE: python em_haplotype.py genotype_file output_file [iterations] \n Takes as input a list of genotypes (defined by 0,1,2), one per line. \n Computes the expectation maximization phasing of the input, and writes most liekly explanations to output_file \n By default does 50 rounds of EM, this can be changed with the third parameter.')
	inFile=sys.argv[1]
	outFile=sys.argv[2]
	if (len(sys.argv) == 4): 
		counter= int(sys.argv[3])
	else: 
		counter = 50

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
			exHap = explainingHaplotypes([int(a) for a in g])
			G[g] =  [map(int, list(g)), 1, 0.0, 0.0, \
			{''.join([str(t) for t in l[0]])+'/'+''.join([str(t) for t in l[1]]):[0, l] for l in exHap}]
		else: 
			G[g][1]+=1

	# set frequencies
	for g in G.keys():
		G[g][2] = float(G[g][1])/len(genotypeList)

	# get list of possible haplotypes
	haps = getHaplotypes([l[0] for l in G.values()])
	# initial frequencies: assume all are equally likely 
	# haplotype data. dictionary haplotype: [probability, list of ints]
	H = {''.join([str(t) for t in l]):[float(1)/len(haps), l] for l in haps}
	
	# do the EM algorithm 
	# right now, do 100 steps
	#counter = 15
	#print G
	# print H
	while counter > 0:
		#print '-----------------' + str(counter) 
		result = eStep(G, H)
		G = result[0]
		H = mStep(G, 	H)
		counter -= 1
		# for g in G.items():
		# 	print g
		# for h in H.items():
		# 	print h

	# select the most likely explanation for each genotype and write to output file
	finalList = []
	for g in genotypeList:
		ex = G[g][4].items()
		maxProb = 0
		maxItem = ex[0][0]
		for e in ex:
			if e[1][0] > maxProb:
				maxProb = e[1][0]
				maxItem = e[0]
		finalList.append(maxItem)

	with open(outFile, 'w') as of:
		for  f in finalList:
			of.write(f.split('/')[0] + '\n')
			of.write(f.split('/')[1] + '\n')


# eStep(gDict, hDict): takes in the dictionaries of genotypes and associated data, haplotypes and assoicated data
# and computes one iteration of the expectation step
# returns [updated gDict, log-likelihood]
def eStep(gDict, hDict):
	# going to calculate the probability of genotype j
	sumPj = 0
	for j in gDict.keys():
		Pj = 0
		# sum over all explanations...
		explanations = gDict[j][4]
		for pair in explanations.keys():
			if explanations[pair][1][0] == explanations[pair][1][1]:
				# add square of frequency to sum
				prob = hDict[''.join([str(t) for t in explanations[pair][1][0]])][0] * hDict[''.join([str(t) for t in explanations[pair][1][0]])][0]
				Pj += prob
				gDict[j][4][pair][0]=prob 
				#print 'equal ' + ''.join([str(t) for t in explanations[pair][1][0]]) + ' : ' + str(prob)
			else:
				# add 2 * pk * pl
				prob = 2*hDict[''.join([str(t) for t in explanations[pair][1][0]])][0] * hDict[''.join([str(t) for t in explanations[pair][1][1]])][0]
				Pj += prob
				gDict[j][4][pair][0]=prob 
				#print 'unequal ' + ''.join([str(t) for t in explanations[pair][1][0]]) + ' , ' + ''.join([str(t) for t in explanations[pair][1][1]]) +' : ' + str(prob)


		# total probability for genotype J
		gDict[j][3] = Pj
		sumPj += Pj

		# update individual explanation probabilities
	#for j in gDict.keys():
		for pair in gDict[j][4].keys():
			gDict[j][4][pair][0] = gDict[j][4][pair][0] / Pj 

	# calculate likelihood
	likelihood = 1
	for j in gDict.keys():
		likelihood = likelihood * math.pow(gDict[j][3],gDict[j][1])
	#print "log lik: " + str(math.log(likelihood))

	# return updated data
	return [gDict, math.log(likelihood)]

# mStep(gDict, hDict): takes in the dictionaries of genotypes and associated data, haplotypes and assoicated data
# and computes one iteration of the maximization step
# returns updated hDict
def mStep(gDict, hDict):
	# update haplotype frequencies now!
	for r in hDict.keys():
		#print 'haplotpype: ' + r
		#haplotpype as list
		hap  = hDict[r][1]
		prob = 0
		# for every genotype
		for j in gDict.keys():
			#print 'genotype: ' + j
			# genotype as list
			gen = gDict[j][0]
			exProb = 0
			# for every explanation
			for i in gDict[j][4].keys():
				#print 'explanation: ' + i
				# multiply by 0,1,or 2
				delta = gDict[j][4][i][1].count(hap)
				#print 'delta: ' + str(delta)
				exProb += (delta * gDict[j][4][i][0])
				#print 'exProb: ' + str(exProb)
			prob += exProb
		# update r with final prob
		hDict[r][0] = prob * 0.5

	return hDict

# explainingHaplotypes(genotype): given a genotype as a list of ints,
# returns a list of haplotype pairs that can explain that genotype
def explainingHaplotypes(genotype): 
	explainList = []
	# get possible haplotypes for the individual genotype
	pHaps = getHaplotypes([genotype])
	# solve this the naive way so I don't have to rewrite getHaplotypes...
	for hap in pHaps:
		#print hap
		pair = [hap, resolveGenotype(genotype, hap)]
		if (pair not in explainList) and (pair[::-1] not in explainList):
			explainList.append(pair)
	return explainList

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
			temp=[]
			for pos in u:
				if pos == 2:
					if len(temp) ==0:
						temp=[[0],[1]]
					else:
						temp0 = [l+[0] for l in temp]
						temp1 = [l+[1] for l in temp]
						temp = temp0 + temp1
				else:
					if len(temp) ==0:
						temp = [[pos]]
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

if __name__ == '__main__':
	main()