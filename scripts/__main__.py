import os, sys
from .sw import optimize_gap_penalties, import_pairs, temporary
from .hybrid import init_files, optimize_gap_penalties, all_score_mats, read_fasta, get_sequences, read_matrix, align, get_negpair_seq, get_pospair_seq, make_roc_curve, single_scoring

# Set working directory and sequence folder location
pwd = os.getcwd()
directory = pwd + "/" + sys.argv[1]

# Import positive and negative pairs, and scoring matrices
negpairlist_filename, pospairlist_filename, negpairs1, negpairs2, pospairs1, pospairs2, names, options = init_files(pwd, directory)

# Run optimization using BLOSUM50 and +/- pair lists
# gap, extension = optimize_gap_penalties(options[0], pospairs1, pospairs2, negpairs1, negpairs2)

# Find best scoring matrix by false positive rate
# all_score_mats(options, gap, extension, pospairs1, pospairs2, negpairs1, negpairs2)
## all_score_mats(options, 4, 3, pospairs1, pospairs2, negpairs1, negpairs2, names) # hardcoded penalties for testing

# Normalize best results by seq length
for i, each in enumerate(options):
	## fpr, pos_align_score, neg_align_score = single_scoring(each, gap, extension, pospairs1, pospairs2, negpairs1, negpairs2)
	fpr, pos_align_score, neg_align_score = single_scoring(each, 4, 3, pospairs1, pospairs2, negpairs1, negpairs2)
	make_roc_curve(pos_align_score, neg_align_score, each, names[i])
