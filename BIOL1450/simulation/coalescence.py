import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import sys
import random
import time

# algorithm for coalescent simulations 
# Set k = n, T=0, V = {V1, V2, , Vn} (vertices or leaves) and A1 = A2
# =  = An = 0 (ages).
# Draw t from an exponential distribution with mean 2/(k*(k-1))
# Draw a and b uniformly, with . These are the indices of two 
# nodes in the tree to coalesce.
# Set k = k-1, T = T+t
# Connect Va and Vb to a new node V2n-k and set A2n-k = T
# Remove Va and Vb from V and add V2n-k
# Stop when k = 1.

#### SOME HELPER FUCTIONS

# need to do a lot of work to make the graph layout look nice
# layout(graph, Gpos): applies a layout algorithm to the coalescence tree. returns a dictionary 
# of x positions
def Glayout(graph, root, k):
	x = {n:0 for n in graph.nodes()}
	x[root] = k/2.0

	a = nx.dfs_preorder_nodes(graph, root)
	n = next(a, None)
	while n != None:
		#print n
		if graph[n] != {}:
			left = graph[n].keys()[0]
			right = graph[n].keys()[1]
			val = subgraphX(graph, n, x[n])
			x[left] = val[0] 
			x[right] = val[1] 
		n = next(a, None)
	return x

# gets the X position for a subgraph
def subgraphX(graph, node, prev):
	left = prev - subgraphLeaves(graph, graph[node].keys()[0]) /2.0
	right = prev + subgraphLeaves(graph, graph[node].keys()[1]) /2.0
	return [left, right]

# subgraphLeaves: computes the number of leaves in a subgraph
def subgraphLeaves(graph, node):
	#print node
	if graph[node] == {}:
		return 1
	else: 
		return subgraphLeaves(graph, graph[node].keys()[0]) + subgraphLeaves(graph, graph[node].keys()[1])

# returns the nodes in a connected subgraph, following directed edges. 
def getSubgraph(graph, node):
	a = nx.dfs_successors(graph, node)
	if a == {}:
		return [node]
	else:
		return [node] + getSubgraph(graph, a[node][0]) + getSubgraph(graph, a[node][1])

# computing with changes in time. 
# need to stretch or shrink time locally depending if time at the current is greater 
# or less than time at the start. 
# popGrowth(growth, factor, k): returns a list of relative population sizes, given
# the growht model and growth factor
def popGrowth(growth, factor, k):
	factor=1.0/factor
	# linear growth: increase by percent a of initial every time unit
	if growth == 'linear':
		populations = [1.0]
		for i in range(k-1):
			populations.append(populations[-1]+factor)


	if growth == 'exponential':
	# exponential growth: with parameter b
		populations = [1.0]
		for i in range(k-1):
			populations.append(populations[-1]*factor)

	return populations


def main(k, growth= None, factor = None):	
	# Set initial parameters, list of vertices
	currentNodes = range(1,k+1)
	remaining = k
	eventTimes = []
	coalescePairs = []
	draw = 0
	# if we want to simulate population change, calculate that here. 
	if growth != None:	
		populations = popGrowth(growth, factor, k)
		populations = populations[::-1]
	# do the randomized drawing
	# while we have more than one node left:
	while remaining > 1:
		# draw t from exponential. scale by population change factor if we need to. 
		if growth !=None:
			eventTimes.append(np.random.exponential(scale=float(2)/(remaining*(remaining-1)), size=1)[0] * populations[-draw])
		else:
			eventTimes.append(np.random.exponential(scale=float(2)/(remaining*(remaining-1)), size=1)[0])
		#print 'time for this event: ' + str(eventTimes[-1])
		draw +=1
		# draw a and b uniformly 
		sample = random.sample(currentNodes,2)
		a=sample[0]
		b=sample[1]
		#print 'nodes to coalesce: ' + str(a) + ' , ' + str(b)
		coalescePairs.append([a,b])

		# subtract one from remaining nodes
		remaining += -1

		# add new node from coalescent event
		currentNodes.append(k+(k-remaining))

		#remove old nodes
		currentNodes.remove(a)
		currentNodes.remove(b)
		#print currentNodes

	# at the end, calculate cumulative times to coalescence 
	totalTimes = np.cumsum(eventTimes)

	# graph to represent the tree for later	
	G=nx.DiGraph()
	# add nodes to represent initial population
	G.add_nodes_from(range(1,k+1))

	# add new edges for the coalesced pairs
	# update the graph with new vertix and time
	for i in range(len(coalescePairs)):
		a = coalescePairs[i][0]
		b = coalescePairs[i][1]
		nodeNum = k+i+1
		#intermediate time
		t =  eventTimes[i]

		G.add_edge(nodeNum, a, length=t)
		G.add_edge(nodeNum, b, length=t)

	# calculate layout for graph
	xpos = Glayout(G, k+k-1, k)
	# set initial positions for drawing
	Gpos = {a: np.array([xpos[a], 0.0]) for a in range(1,k*2)}

	for i in range(len(coalescePairs)):
		# update y position of new node, from total coalescece time
		nodeNum = k+i+1
		Gpos[nodeNum][1] = totalTimes[i]

	#interactive plotting
	plt.ion()
	# start figure. configure subplots
	fig = plt.figure()
	fig.show()
	# 3 panels if modeling growth
	if growth != None:
		ax = fig.add_subplot(111)    # The big subplot
		ax1 = fig.add_subplot(131)
		ax1.set_title('Coalescence Graph')
		ax2 = fig.add_subplot(132)
		ax2.set_title('Time to Coalescence')
		ax3 = fig.add_subplot(133)
		ax3.set_title('Population size')
	# 2 panels otherwise
	else:
		ax = fig.add_subplot(111)    # The big subplot
		ax1 = fig.add_subplot(121)
		ax1.set_title('Coalescence Graph')
		ax2 = fig.add_subplot(122)
		ax2.set_title('Time to Coalescence')
	# make things look (slightly more) pretty
	plt.tight_layout()

	# plot initial nodes at the start	
	nx.draw(G.subgraph(range(1,k+1)),pos=Gpos, with_labels=True, node_size=300, ax =ax1)
	
	# plot future coalescence events
	for i in range(k+1, 2*k):
		# plot new subgraphs
		nx.draw(G.subgraph(getSubgraph(G, i)),pos=Gpos, with_labels=True, node_size=300, arrows=False, ax=ax1)
		# plot time to coalescnce
		ax2.plot(range(i-k), eventTimes[:i-k], color='blue')
		if growth != None:
			# plot population growh
			ax3.plot(range(i-k), populations[:i-k], color='red')
		
		# push update to figure
		plt.draw()

		# wait between drawing
		time.sleep(0.25)	

	# keep it open
	plt.show()

# check usage - eventually use argparse (nope no time for that)
if __name__ == '__main__':
	if not ((len(sys.argv) == 2) or (len(sys.argv)== 4)):
		sys.exit('USAGE: coalescence.py k growth_type growth_factor \n \
		 k is the starting population size and the only mandatory argument \n \
		 Changes in population size can be modeled by specifying \'exponential\' or \n \
		 \'linear\' for growth_type. Then you must specify a growth_factor.')

	#get input and set initial parameters 
	k=int(sys.argv[1])
	assert k >1

	# set growth if we need to
	if len(sys.argv) == 4:
		growth = sys.argv[2]
		factor = float(sys.argv[3])
		# run it!
		main(k, growth, factor)

	else:
		# run it! (without growth)
		main(k)
