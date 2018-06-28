import os, sys, subprocess
import numpy as np, matplotlib.pyplot as plt, pandas as pd
from skbio import Protein
from skbio.alignment import local_pairwise_align, local_pairwise_align_ssw
from sklearn import metrics
from numba import jit

def read_seq(f_a,f_b):
	'''Read in two sequences for alignment'''

	# Read the first sequence, skipping header
	with open(f_a, 'r') as f1:
		next(f1) 
		seq1 = ''.join(line.strip() for line in f1)
		seq1=seq1.strip()
		f1.close()
	# Read the second sequence, skipping header
	with open(f_b, 'r') as f2:
		next(f2) 
		seq2 = ''.join(line.strip() for line in f2)
		seq2=seq2.strip() 
		f2.close()

	return seq1, seq2

def read_matrix (filepath):
	with open(filepath) as mat:
		matrix = mat.read()
		lines = matrix.strip().split('\n')
		score_matrix = []
		for line in lines:
			if not line.startswith("#"):
				score_matrix.append(line)
	aa_column = pd.DataFrame(score_matrix[0].strip().split("  "))
	scores = []
	for i in range(len(score_matrix)-1):
		ls = score_matrix[i+1].strip()
		new = ls.split(" ")
		sub = []
		for n in new:
			if n != "":
				sub.append(int(n))
		scores.append(sub)
		score_df = pd.DataFrame(scores)
		score_df.columns = score_matrix[0].strip().split("  ")
	df_score = aa_column.join(score_df)
	df_score = df_score.set_index([0])
	return df_score

@jit
def align(seq1, seq2, go, ge):
	''' Perform alignment using scikit-bio for any two given sequences, gap penalties, and score matrix. '''
	
	a, b = read_seq(seq1, seq2)
	scoreMatrix = read_matrix(sys.argv[2])
	# alignment, score, start_end_positions = local_pairwise_align_protein(Protein(a, lowercase = True), Protein(b, lowercase = True), gap_open_penalty=go, gap_extend_penalty=ge, substitution_matrix=None)
	alignment, score, start_end_positions = local_pairwise_align_ssw(Protein(a, lowercase = True), Protein(b, lowercase = True), gap_open_penalty=go, gap_extend_penalty=ge, substitution_matrix=scoreMatrix) 
	print ("\nScore:", score)
	return score

def new_align(seq1, seq2, gap_open_penalty, gap_extend_penalty, substitution_matrix):
	alignment, score, start_end_positions = local_pairwise_align(Protein(seq1, lowercase = True), Protein(seq2, lowercase = True), gap_open_penalty, gap_extend_penalty, substitution_matrix)
	return score

@jit
def optimize_gap_penalties():
	''' Run alignment with different gap penalties to find best ones. '''

	print ("Running main alignment function using data from %s" %sys.argv[1])
	pwd = os.getcwd()
	folder = pwd + "/" + sys.argv[1]
	l = subprocess.check_output(["ls", "%s" % folder], universal_newlines=True).split()
	results = np.zeros((len(l), len(l)))
	normres = np.zeros((len(l), len(l)))
	
	for y in range(1,21):
		go = y
		for z in range(1,6):
			ge = z
			print (go, ge)
			# Run alignment for given sequences
			for i, fa in enumerate(l):
				for j in range(len(l)):
					print (fa, l[j])
					# results[i][j] = align(folder + fa, folder + l[j], go, ge)
					results[i][j] = new_align(folder + fa, folder + l[j], go, ge, read_matrix(sys.argv[2]))
			maxres = np.ndarray.max(results)
			for index, r in np.ndenumerate(results):
				normres[index] = r / maxres
			print(normres)
			
	# print ("Gap optimization:\n",temp)

optimize_gap_penalties()

''' skbio.alignment.local_pairwise_align_protein(seq1, seq2, gap_open_penalty=11, gap_extend_penalty=1, substitution_matrix=None) '''
