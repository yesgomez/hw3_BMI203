import os, sys
from .sw import *
#from .sw import init_files, optimize_gap_penalties, all_score_mats, read_fasta, get_sequences, read_matrix, align, get_negpair_seq, get_pospair_seq, make_roc_curve, single_scoring, norm_score

# Set working directory and sequence folder location
pwd = os.getcwd()
directory = pwd + "/" + sys.argv[1]

# Import positive and negative pairs, and scoring matrices
negpairlist_filename, pospairlist_filename, negpairs1, negpairs2, pospairs1, pospairs2, names, options = init_files(pwd, directory)

## Run optimization using BLOSUM50 and +/- pair lists
# gap, extension = optimize_gap_penalties(options[0], pospairs1, pospairs2, negpairs1, negpairs2)

## Find best scoring matrix by false positive rate
# all_score_mats(options, gap, extension, pospairs1, pospairs2, negpairs1, negpairs2)

## Normalize best results by seq length
# for i, each in enumerate(options):
	# fpr, pos_align_score, neg_align_score = single_scoring(each, gap, extension, pospairs1, pospairs2, negpairs1, negpairs2)
	# make_roc_curve(pos_align_score, neg_align_score, each, names[i])
	
	# pos_norm_score = norm_score(pos_align_score, pospairs1, pospairs2)
	# neg_norm_score = norm_score(neg_align_score, negpairs1, negpairs2)
	# make_roc_curve(pos_norm_score, neg_norm_score, each, names[i]+'norm')

fpr, pos_align_score, neg_align_score = single_scoring(options[0], 4, 3, pospairs1, pospairs2, negpairs1, negpairs2)
	
