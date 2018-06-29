import os, sys
from .sw import optimize_gap_penalties, import_pairs, temporary
from .hybrid import init_files, penalty_optimize, all_score_mats, read_fasta, get_sequences, read_matrix, align, get_negpair_seq, get_pospair_seq, penalty_optimize

# Set working directory and sequence folder location
pwd = os.getcwd()
directory = pwd + "/" + sys.argv[1]

# Import positive and negative pairs, and scoring matrices
negpairlist_filename, pospairlist_filename, negpairs1, negpairs2, pospairs1, pospairs2, names, options = init_files(pwd, directory)

# Run optimization using BLOSUM50 and pair lists
gap, extension = optimize_gap_penalties(options[0], pospairs1, pospairs2, negpairs1, negpairs2)
all_score_mats(options, gap, extension, pospairs1, pospairs2, negpairs1, negpairs2)

# Normalize best results by seq length
	# For best resulting matrix
	# Return all val as matrix and all lengths as matrix, then divide element-wise
	# Voila - new matrix!