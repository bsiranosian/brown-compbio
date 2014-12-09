The program coalescence.py can be used to simulate coalescent events
looking back in time. It can take a population size and model of growth
and produce a tree of coalescence events where the height of a node is relative
to the time of coalescence for that node. 

Two or three plots are produced. The first is always a graph of the coalescence 
tree. The second is the time to coalescence for each event, looking back in time.
The third only shows up if you simulate population growth. It is the relative population
size looking FORWARD in time. 

USAGE: python coalescence.py k [growth_type] [growth_rate]
MANDATORY ARGUMENTS: specify a population of alleles of size k
OPTIONAL ARGUMENTS: growth_type can be one of 'exponential' or 'linear'
					if growth type is specified, growth_rate must also be specified
					growth_rate serves to regulate how the population changes looking
					forward in time. Specify exponental and a rate of 2 to double every
					generation, etc. With linear growth, the value for rate will increase 
					the population by that fraction of the original value in every genereation
					So a rate of 0.05 will increase population by 5% of original every time. s
DEPENDS: numpy, matplotlib, networkx.
		 These packages can be downloaded from the internet.