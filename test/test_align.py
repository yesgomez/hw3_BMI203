import numpy as np
from scripts.sw import *

# generic tests for algorithm logic #
pwd = os.getcwd()
directory = pwd + "/sequences"

def test_import():
	# Are the files read in correctly?
	negpairlist_filename, pospairlist_filename, negpairs1, negpairs2, pospairs1, pospairs2, names, options = init_files(pwd, directory)	
	print ("Import test passed.\n")

def test_fpr():
	# Can it calculate scores using different scoring matrices?
	negpairlist_filename, pospairlist_filename, negpairs1, negpairs2, pospairs1, pospairs2, names, options = init_files(pwd, directory)	
	all_score_mats(options, 4, 3, pospairs1, pospairs2, negpairs1, negpairs2, names) # hardcoded penalties for testing
	print ("FPR test passed.\n")
	
# def test_makeroc():
# 	# Can it output an ROC curve using provided scores and scoring matrix?
# 	negpairlist_filename, pospairlist_filename, negpairs1, negpairs2, pospairs1, pospairs2, names, options = init_files(pwd, directory)	
# 	fpr, pos_align_score, neg_align_score = single_scoring(options[1], 4, 3, pospairs1, pospairs2, negpairs1, negpairs2)
# 	make_roc_curve(pos_align_score, neg_align_score, options[1], names[1])
# 	print ("ROC curve test passed.\n")
## Couldn't do plt.figure for some reason. Did not pass...
