import sys

# read from infile, output vcf to outfile
def main(infile, outfile_base):
	#open file and read
	with open(infile,'r') as inf:
		lines = [l.strip().split('\t') for l in inf.readlines()]
		#print lines[0]

	#change Chr to chr
	for i in range(len(lines)):
		lines[i][0] =  lines[i][0].lower()

	# Sep into snp, insert, deletions
	snp = [l for l in lines if int(l[2]) ==0]
	insertions = [l for l in lines if int(l[2]) > 0]
	deletions = [l for l in lines if int(l[2]) <	 0]

	# write SNPS
	with open(outfile_base+'_SNP.bed', 'w') as of:
		for l in snp:
			of.write(l[0]+'\t'+l[1]+'\t'+str(int(l[1])+1)+'\t'+l[3]+'/'+l[4]+'\n')
	# write insertions
	with open(outfile_base+'_insertions.bed', 'w') as of:
		for l in insertions:
			of.write(l[0]+'\t'+l[1]+'\t'+str(int(l[1])+1)+'\t'+l[3]+'/'+l[4]+'\n')
	# write deletions
	with open(outfile_base+'_deletions.bed', 'w') as of:
		for l in deletions:
			of.write(l[0]+'\t'+l[1]+'\t'+str(int(l[1])+1)+'\t'+l[3]+'/'+l[4]+'\n')
		
		

if __name__ == '__main__':
	if len(sys.argv) != 3:
		sys.exit("USAGE: python austin_sdi.py infile outfile_base  \n Creates three files by splitting the SDI file up into insertions, deletions and snps. 3 files are nammed according to the base name specified.")
	main(sys.argv[1],sys.argv[2])