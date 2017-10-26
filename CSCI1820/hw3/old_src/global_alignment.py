#global_alignment(seq1, seq2, match_score, mismatch_score, gap_penalty, return_score): returns the optimal local alignment of two sequences by using the given parameters.
# uses the needleman-wunsch algorithm for alignment. Depends on helper function single_score. 
# if return_score is true, simply return the highest score for the alignment. otherwise, return a list of the two aligned strings. 
def global_alignment(seq1, seq2, match_score, mismatch_score, gap_penalty, return_score):
	n = len(seq1)
	m = len(seq2)
	score_matrix=[[0 for i in range(m+1)] for j in range(n+1)]
	trace_matrix = [[0 for i in range(m+1)] for j in range(n+1)]
	for i in range(1, n+1):
		for j in range(1, m+1):
			#arguments to smith-waterman maximization 
			args=[(score_matrix[i-1][j-1] + single_score(seq1[i-1], seq2[j-1], match_score, mismatch_score)),
				(score_matrix[i-1][j] + gap_penalty),
				(score_matrix[i][j-1] + gap_penalty)]

			#pick max
			score_matrix[i][j] = max(args)
			#argmax for traceback
			trace_matrix[i][j] = args.index(max(args))

	# for row in score_matrix:
	# 	print row
	# print " --- "
	# for row in trace_matrix:
	# 	print row
	# print " ****** "

	# For GLOBAL alignment we want to start at the LAST score in the matrix.
	# Traceback
	al1 = ''
	al2 = ''
	cur_i = n
	cur_j = m

	while (cur_i!=0) or (cur_j!=0):
		if cur_i ==0:
			al1 += "-"
			al2 += seq2[cur_j-1]
			cur_j+=-1
		elif cur_j ==0:
			al2 += "-"
			al1 += seq1[cur_i-1]
			cur_i+=-1
		#if match or mismatch
		elif trace_matrix[cur_i][cur_j] == 0:
			al1 += seq1[cur_i-1]
			al2 += seq2[cur_j-1]
			cur_i += -1
			cur_j += -1
			
		#if gap in i 
		elif trace_matrix[cur_i][cur_j] == 1:
			al1 += seq1[cur_i-1]
			al2 += "-"
			cur_i += -1
		#if gap in j
		elif trace_matrix[cur_i][cur_j] == 2:
			al1 += "-"
			al2 += seq2[cur_j-1]
			cur_j += -1

	if return_score:
		return score_matrix[n][m]
	else:
		return[al1[::-1],al2[::-1]]

#single_score(char1, char2): defines the match and mismatch scores for a single base
def single_score(char1, char2, match_score, mismatch_score):
	if char1==char2:
		return match_score
	else:
		return mismatch_score

global_alignment("CATG","ATC",10,-10,-5,False)
global_alignment("MISMATCHINGATGATG","ATGCATG",1,-1,-1,False)