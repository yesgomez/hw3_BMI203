# coding: utf-8
import numpy as np
import os
import re
import glob
from skbio import Protein
from skbio.alignment import local_pairwise_align, TabularMSA
import itertools 
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
import pandas as pd


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

def get_sequences(dir):
    files = glob.glob(dir + '/*.fa')
    sequences = []
    for filepath in glob.iglob(os.path.join(dir, "*.fa")):
        with open(filepath) as fp:
            for name, seq in read_fasta(fp):
                sequences.append(seq)
    print("Read in %d fasta files"%len(files))
    return sequences


def read_score_matrix(filepath):
    with open(filepath) as mat:
        matrix = mat.read()
        lines = matrix.strip().split('\n')
        score_matrix = []
        for line in lines:
            if not line.startswith("#"):
                score_matrix.append(line)
    aas_column = pd.DataFrame(score_matrix[0].strip().split("  "))
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
    df_score = aas_column.join(score_df)
    df_score = df_score.set_index([0])
    return df_score


def align(seq1, seq2, open_penalty, extend_penalty, substitution_matrix):
    '''
    Given two sequences to algin substitution matrix, gap open penalty and
    extend penalty, output the alignment score matrix.

    '''
    # Initialization of score matrix and state state matrix by making an empty matrix.
    # Score Matrix: Matrix with alignment score calulated in each cell
    # State Matrix: Matrix with info on where maximum score is calculated

    mat_score = [[None]*len(seq2) for _ in range(0,len(seq1))]
    mat_state = [[None]*(len(seq2)+1) for _ in range(0,len(seq1)+1)]
    mat_score.insert(0,[0] * len(seq2))
    # Set the first column equal to zero
    for row in mat_score:
        row.insert(0,0)
	
	# Start scoring matrix by first initializing max score as zero.     
    for i in range(1,len(seq1)+1):
        for j in range(1, len(seq2)+1):        
            max_score = 0
            state = ""
	# Calculate left well deletion score by determining whether it is gap open or extension 
            if mat_state[i][j-1] == "del":
                score = mat_score[i][j-1] + extend_penalty
            else:
                score = mat_score[i][j-1] + open_penalty

            if score > max_score:
                max_score = score
                state = "del"
	# Calculate up well insertion score by determining whether it is gap open or extension         
            if mat_state[i-1][j] == "ins":
                score = mat_score[i-1][j] + extend_penalty
            else:
                score = mat_score[i-1][j] + open_penalty

            if score > max_score:
                max_score = score
                state = "ins"
	# Calculate diagonal well alignment score by reffering to substitution matrix 
            aa_i = seq1[i-1]
            aa_j = seq2[j-1]
            align_score = substitution_matrix.loc[str(aa_i),str(aa_j)]
            score = mat_score[i-1][j-1] + align_score
            if score > max_score:
                max_score = score
                state = "align"

            mat_score[i][j] = max_score
            mat_state[i][j] = state

    return np.asarray(mat_score), np.asanyarray(mat_state), max_score, state


def get_negpair_seq(filepath, negpairlist_filename):
    '''
    Get sequences for negative pairs.
    Input: Directory for all the sequences
    Output: List of negative pairs sequences
    '''
    files = glob.glob(filepath + '/*.fa')
    negpair_list = open(negpairlist_filename).read().splitlines()
    neg_files, neg_files2 = [],[]
    negpair_sequences, negpair2_sequences  = [], []
    for seq_file in files:
        for negpair in negpair_list:
            if negpair[10:22] in seq_file[-12:]:
                neg_files.append(seq_file)
            if negpair[33:] in seq_file[-12:]:
                neg_files2.append(seq_file)
    for neg_file in neg_files:
        with open(neg_file) as neg:
            for name, neg_seq in read_fasta(neg):
                negpair_sequences.append(neg_seq)
    for neg_file in neg_files2:
        with open(neg_file) as neg:
            for name, neg2_seq in read_fasta(neg):
                negpair2_sequences.append(neg2_seq)
    #Some of the sequences has "x" as a gap. I'm replacing those 'x's with '-'
    #so it works with skbio package. Not really a good way but it works.
    without_x = [s for s in negpair2_sequences if 'x' not in s]
    with_x= [s for s in negpair2_sequences if 'x' in s]
    with_x =' '.join(with_x)
    with_x = with_x.replace('x', '*')
    without_x.append(with_x)
    negpair2_sequences = without_x
    
    print("Read in %d negative pair fasta files" %len(negpair_sequences))
    print("Read in %d negative pair fasta files" %len(negpair2_sequences))
    return negpair_sequences, negpair2_sequences

def get_pospair_seq(filepath,pospairlist_filename):
    '''
    Get sequences for positive pairs.
    Input: Directory for all the sequences
    Output: List of positive pairs sequences
    '''
    files = glob.glob(filepath + '/*.fa')
    pospair_list = open(pospairlist_filename).read().splitlines()
    pos_files,pos2_files = [],[]
    pospair_sequences,pospair2_sequences = [],[]
    for seq_file in files:
        for pospair in pospair_list:
            if pospair[10:22] in seq_file[-12:]:
                pos_files.append(seq_file)
            if pospair[33:] in seq_file[-12:]:
                pos2_files.append(seq_file)
    for pos_file in pos_files:
        with open(pos_file) as pos:
            for name, pos_seq in read_fasta(pos):
                pospair_sequences.append(pos_seq)
    for pos_file in pos2_files:
        with open(pos_file) as pos:
            for name, pos2_seq in read_fasta(pos):
                pospair2_sequences.append(pos2_seq)    
    print("Read in %d positive pair fasta files" %len(pospair_sequences))
    print("Read in %d positive pair fasta files" %len(pospair2_sequences))
    return pospair_sequences, pospair2_sequences


directory = '/Users/student/Winter2018/Algorithms/hw3_BMI203/sequences/'
negpairlist_filename = '/Users/student/Winter2018/Algorithms/hw3_BMI203/Negpairs.txt'
pospairlist_filename = '/Users/student/Winter2018/Algorithms/hw3_BMI203/Pospairs.txt'
negpairs1, negpairs2 = get_negpair_seq(directory, negpairlist_filename)
pospairs1, pospairs2 = get_pospair_seq(directory, pospairlist_filename)
blosum50 = read_score_matrix('/Users/student/Winter2018/Algorithms/hw3_BMI203/matrices/BLOSUM50')
blosum62 = read_score_matrix('/Users/student/Winter2018/Algorithms/hw3_BMI203/matrices/BLOSUM62')
matio = read_score_matrix('/Users/student/Winter2018/Algorithms/hw3_BMI203/matrices/MATIO')
pam100 = read_score_matrix('/Users/student/Winter2018/Algorithms/hw3_BMI203/matrices/PAM100')
pam250 = read_score_matrix('/Users/student/Winter2018/Algorithms/hw3_BMI203/matrices/PAM250')


possible_gap_extend_penalties = list(range(1,6))
possible_gap_open_penalties = list(range(1,21))
best_fpr = float('inf')
best_gap = None
best_extension = None

for open_penalty in possible_gap_open_penalties: 
    for extend_penalty in possible_gap_extend_penalties:
        print ("At open penalty score " + str(open_penalty) + " and extend penalty "+ str(extend_penalty))
        pos_align_score = []
        neg_align_score = []
        for i in range(0,len(pospairs1)):
            mat_score_p,mat_state_p, score_p, state_p = align(pospairs1[i],pospairs2[i], -open_penalty,-extend_penalty, blosum50)
            print ("In " + str(i) + " Pospair: score_p is " + str(score_p))
            mat_score_n,mat_state_n, score_n, state_n = align(negpairs1[i],negpairs2[i], -open_penalty,-extend_penalty, blosum50)
            print ("In " + str(i) + " Negpair: score_n is " + str(score_n))
            pos_align_score.append(score_p)
            neg_align_score.append(score_n)
        pos_align_score.sort()
        neg_align_score.sort()
        tpr = 0.7 #tpr = #>thresh / total scores
        cutoff_index = int((1-tpr)*len(pos_align_score))# find the index that makes it so tpr*lenscores is > threshold
        threshold = pos_align_score[cutoff_index]
        neg_score_np = np.array(neg_align_score)
        fpr = len(neg_score_np[neg_score_np > threshold])/len(neg_score_np) #False Positive Rate
        if fpr < best_fpr:
            best_fpr = fpr
            best_gap = open_penalty
            best_extension = extend_penalty
            
            print ("Best gap now is " + str(best_gap))
            print ("Best extend now is "+ str(best_extension))
    
print('Best Open and extension gap penalty scores are %d and %d.' % (best_gap,best_extension))
print('Best False positive rate is %d.' %(best_fpr))

