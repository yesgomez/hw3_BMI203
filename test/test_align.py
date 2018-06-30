import numpy as np
from scripts.hybrid import *

# generic tests for algorithm logic #
pwd = os.getcwd()
directory = pwd + "/" + sys.argv[1]

def test_import():
	# Are the files read in correctly?
	negpairlist_filename, pospairlist_filename, negpairs1, negpairs2, pospairs1, pospairs2, names, options = init_files(pwd, directory)	
	print ("Import test passed.\n")

def test_fpr():
	# Can it calculate scores using different scoring matrices?
	all_score_mats(options, 4, 3, pospairs1, pospairs2, negpairs1, negpairs2, names) # hardcoded penalties for testing
	print ("FPR test passed.\n")
	
def test_makeroc():
	# Can it output an ROC curve using provided scores and scoring matrix?
	fpr, pos_align_score, neg_align_score = single_scoring(each, 4, 3, pospairs1, pospairs2, negpairs1, negpairs2)
	make_roc_curve(pos_align_score, neg_align_score, each, names[i])
	print ("ROC curve test passed.\n")

# def odd_even():