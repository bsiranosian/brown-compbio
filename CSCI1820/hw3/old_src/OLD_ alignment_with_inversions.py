def alignment_score(seq1, seq2, gap_penalty):
	n = len(seq1)
	m = len(seq2)
	score_matrix=[[0 for i in range(m+1)] for j in range(n+1)]

	#set initial values of matrix
	for k in range(n+1):
		score_matrix[k][0] = gap_penalty *k
	for k in range(m+1):
		score_matrix[0][k] = gap_penalty * k 	

	for i in range(1, n+1):
		for j in range(1, m+1):

			gap1 =[0 for k in range(i+1)]
			for k in range(i+1):
				gap1[k] = score_matrix[i-1-k][j-1] + (gap_penalty * (k+1))

			gap2 = [0 for k in range(j+1)]
			for k in range(j+1):
				gap2[k] = score_matrix[i-1][j-1-k] + (gap_penalty * (k+1))

			print gap1
			print gap2

			score_matrix[i][j] = max((score_matrix[i-1][j-1] + single_score(seq1[i-1], seq2[j-1])),
				(max(gap1)),
				(max(gap2)))


			for t in  score_matrix:
				print t

	return score_matrix[n][m]

#single_score(char1, char2): defines the alignment and mismatch scores for a single base
def single_score(char1, char2):
	if char1==char2:
		return 2
	else:
		return -2

# TEST CASES
# alignment_score('','', -1) == 0
# alignment_score('A','A', -1) == 2
# alignment_score('A','C', -1) == 0
# alignment_score('ATCG','ATCG', -1) == 8
# alignment_score('ATTCG','AAACG', -1) == 8



