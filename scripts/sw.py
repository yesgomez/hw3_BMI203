import os, sys, subprocess
import numpy as np, matplotlib.pyplot as plt, pandas as pd
from skbio import Protein
from skbio.alignment import local_pairwise_align, local_pairwise_align_ssw
from sklearn import metrics
from numba import jit


def read_seq(f_a):
	'''Read in two sequences for alignment'''

	# Read the sequence, skipping header
	with open(f_a, 'r') as f1:
		next(f1) 
		seq1 = ''.join(line.strip() for line in f1)
		seq1=seq1.strip()
		print (seq1[0:5])
		f1.close()
	return seq1


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


def import_pairs():
	#Loop through sequence files in folder
	seq1_negpairs = []
	seq2_negpairs = []
	seq1_pospairs = []
	seq2_pospairs = []
	negpairs = []
	pospairs = []

	f = open('Negpairs.txt', 'rt')
	for line in f:
		seq1_negpairs.append(str.split(line)[0])
		seq2_negpairs.append(str.split(line)[1])
	f.close()

	f = open('Pospairs.txt', 'rt')
	for line in f:
		seq1_pospairs.append(str.split(line)[0])
		seq2_pospairs.append(str.split(line)[1])
	f.close()
	
	# Return pairs as single lists
	return seq1_negpairs, seq2_negpairs, seq1_pospairs, seq2_pospairs


@jit
def new_align(seq1, seq2, go, ge, scoreMatrix):
	''' Perform alignment using scikit-bio for any two given sequences, gap penalties, and score matrix. '''
	a, b = read_seq(seq1), read_seq(seq2)
	alignment, score, start_end_positions = local_pairwise_align_ssw(Protein(a, lowercase = True), Protein(b, lowercase = True), gap_open_penalty=go, gap_extend_penalty=ge, substitution_matrix = scoreMatrix)
	print ("\nScore:", score)
	return score


@jit
def optimize_gap_penalties(list_a, list_b):
	''' Run alignment with different gap penalties to find best ones. '''

	print ("Running gap penalty optimization.\n")
	# folder = pwd + "/" + sys.argv[1]
	results = np.zeros((len(list_a), len(list_b)))
	# # normres = np.zeros((len(list_a), len(list_b)))
	# i = j = 0
	for y in range(1,21):
		go = y
		for z in range(1,6):
			ge = z
			print (go, ge)
			# Run alignment for given sequences
			smatrix = read_matrix(sys.argv[2])
			for i in range(list_a):
				print (list_a[i], list_b[i])
				results[y-1][z-1] = new_align(list_a[i], list_b[i], go, ge, smatrix)
			# maxres = np.ndarray.max(results)
			# for index, r in np.ndenumerate(results):
				# normres[index] = r / maxres
			# print(normres)
	return results

	
def temporary():
	print ("Importing data from %s" %sys.argv[1])
	pwd = os.getcwd()
	folder = pwd + "/" + sys.argv[1]
	l = subprocess.check_output(["ls", "%s" % folder], universal_newlines=True).split()
	newl = l
	return l, newl


