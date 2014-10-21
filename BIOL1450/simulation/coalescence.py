import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import sys
import random
import time

# algorithm for coalescent simulations 
# • Set k = n, T=0, V = {V1, V2, …, Vn} (vertices or leaves) and A1 = A2
# = … = An = 0 (ages).
# • Draw t from an exponential distribution with mean 2/(k*(k-1))
# • Draw a and b uniformly, with . These are the indices of two 
# nodes in the tree to coalesce.
# • Set k = k-1, T = T+t
# • Connect Va and Vb to a new node V2n-k and set A2n-k = T
# • Remove Va and Vb from V and add V2n-k
# • Stop when k = 1.

# check usage - eventually use argparse
if len(sys.argv) != 2:
	sys.exit('USAGE: coalescence.py k \n Where k is the starting population size')

# get input and set initial parameters 
k=int(sys.argv[1])
assert k >1
n=k
T=0
# list of vertices
V = range(1,k+1)
# graph to represent the tree for later
G=nx.Graph()
# add nodes to represent initial population
G.add_nodes_from(V)
# dictionary of ages
A = {i:0 for i in V}

# set initial positions for drawing
pos = {a: (a-1,0) for a in V}

#interactive plotting
plt.ion()

# while we have more than one node left:
while len(V) > 1:
	# draw t from exponential
	t = np.random.exponential(scale=float(2)/(k*(k-1)), size=1)[0]
	print 't: ' + str(t)

	# draw a and b uniformly 
	sample = random.sample(V,2)
	a=sample[0]
	b=sample[1]
	print 'nodes to coalesce: ' + str(a) + ' , ' + str(b)

	k=k-1
	T=T+t
	print 'Total time at coalescence: ' + str(T)

	# add new node from coalescent event
	V.append(2*n-k)
	A[2*n-k] = T

	#remove old nodes
	V.remove(a)
	V.remove(b)

	#update the graph with new vertix and time
	G.add_edge(a, 2*n-k, length=T - A[a])
	G.add_edge(b, 2*n-k, length=T - A[b])
	# update position of new node
	pos[2*n-k] = (float(pos[a][0]+pos[b][0])/2,T)
	
	nx.draw(G,pos)
	nx.draw_networkx_labels(G,pos,labels={g:g for g in G.nodes()})
	plt.draw()
	time.sleep(0.25)	


# Notes from meeting
What else can we add to the simulation?
Not always constant population size
	simulate selection 
	growth rate

Put mutations on the tree afterwards so we can see how they appear 
Migration

With recombination:
	Layout? Graph based?

run multliple simulations, disibuation of time across them 