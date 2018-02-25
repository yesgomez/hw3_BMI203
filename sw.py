import os, sys, subprocess
import numpy as np, matplotlib.pyplot as plt
from skbio import Protein
from skbio.alignment import local_pairwise_align_protein
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


def read_matrix(mfile):
	''' Reads in score matrix file and turns it into lists of mapped objects '''
	matrix = {}
	i = 0
	with open(mfile, 'r') as f:
		for line in f:
			if line[0] != '#': # remove header 
				if line[1] == 'A': 
					keys = line.split()
				if line[1] != 'A':
					line = line.split()
					for l in line: int(l)
					matrix[keys[i]] = dict(zip(keys, line))
					i += 1
	# print (matrix)
	return matrix

@jit
def align(seq1, seq2, go, ge):
	''' Perform alignment using scikit-bio for any two given sequences, gap penalties, and score matrix. '''
	a, b = read_seq(seq1, seq2)
	# scoreMatrix = read_matrix(sys.argv[1])
	alignment, score, start_end_positions = local_pairwise_align_protein(Protein(a, lowercase = True), Protein(b, lowercase = True), gap_open_penalty=go, gap_extend_penalty=ge, substitution_matrix=None)
	print ("\nScore:", score)
	return score

# @jit
def main():
	''' Run alignment with different gap penalties to find best ones. '''

	pwd = os.getcwd()
	folder = pwd + "/sequences/"
	l = subprocess.check_output(["ls", "%s" % folder], universal_newlines=True).split()
	results = np.zeros((len(l), len(l)))
	# temp = np.zeros((21, 6))

	mini_list = l[0:19]
	
	for y in range(1,21):
		go = y
		for z in range(1,6):
			ge = z
			print (go, ge)
			# Run alignment for all sequences.
			# for i, fa in enumerate(l):
			# 	for j in range(len(l)):
			# 		print (fa, l[j])
			# 		results[i][j] = align(folder + fa, folder + l[j], go, ge)
			# print(results)

			for i, fa in enumerate(mini_list):
				for j in range(len(mini_list)):
					print (fa, mini_list[j])
					results[i][j] = align(folder + fa, folder + mini_list[j], go, ge)
			print(results, np.matrix.sum(results))
			
			# temp[y][z] = align(a, b, go, ge)

	# temp[0] = range(6)
	# print ("Gap optimization:\n",temp)

main()

''' skbio.alignment.local_pairwise_align_protein(seq1, seq2, gap_open_penalty=11, gap_extend_penalty=1, substitution_matrix=None) '''