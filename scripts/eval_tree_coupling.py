#! /usr/bin/env python

from problin_libs.sequence_lib import read_sequences,read_charMtrx
from problin_libs.eval_lib import *
from treeswift import *

if __name__ == "__main__":
    from sys import argv
    
    groundtruthFile = argv[1]
    testFile = argv[2]
    tag = argv[3]

    true_charMtrx, true_tree = read_groundtruth(groundtruthFile)
    k = len(true_charMtrx['I0'])

    test_charMtrx, test_tree = read_annotation(testFile)

    cells = [node.label for node in test_tree.traverse_leaves()]

    true_answer = tree_coupling(true_tree,cells,true_charMtrx)
    test_answer = tree_coupling(test_tree,cells,test_charMtrx)

    for (x,y) in true_answer:
        s_true, d_true = true_answer[(x,y)]
        s_test, d_test = test_answer[(x,y)]
        d_error = abs(d_true-d_test)/k
        s_score = score_seq(s_true,s_test,soft_assignment=False)
        print(tag,x,y,d_true,d_test,d_error,s_score)
