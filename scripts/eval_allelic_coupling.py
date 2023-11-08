#! /usr/bin/env python

from scmail_libs.sequence_lib import read_sequences,read_charMtrx
from scmail_libs.eval_lib import *
from treeswift import *
from eval_tree_coupling import read_groundtruth

if __name__ == "__main__":
    from sys import argv
    
    groundtruthFile = argv[1]
    testFile = argv[2]
    tag = argv[3]

    true_charMtrx, true_tree = read_groundtruth(groundtruthFile)
    k = len(true_charMtrx['I0'])

    fin = open(testFile,'r')
    test_charMtrx,_ = read_charMtrx(fin,replace_mchar='d',convert_to_int=False)
    fin.close()

    cells = [cell for cell in test_charMtrx.keys()]

    true_answer = tree_coupling(true_tree,cells,true_charMtrx)
    test_answer = allelic_coupling(test_charMtrx,cells,masked_symbol='d')

    for (x,y) in true_answer:
        s_true, d_true = true_answer[(x,y)]
        s_test, d_test = test_answer[(x,y)]
        d_error = abs(d_true-d_test)/k
        s_score = score_seq(s_true,s_test,soft_assignment=False,masked_symbol='d')
        print(tag,x,y,d_true,d_test,d_error,s_score)
