README for em_haplotype.py
Ben Siranosian | CSCI2951 | Fall 2014

USAGE: python em_haplotype.py genotype_file output_file [iterations]

Takes as input a list of genotypes (defined by 0,1,2), one per line. Computes the expectation maximization phasing of the input, and writes most liekly explanations to output_file. By default does 50 rounds of EM, this can be changed with the third parameter. This implementation assumes that all haplotypes are equally likely at the start of the algorithm. 