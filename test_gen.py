from laml_libs.IO_handler.sequence_lib import read_sequences
from laml_libs.PMM_original.EM_solver import EM_solver
from math import *

trees = ["((a:1,b:1)ab:1,c:1)abc:1;","(((a:0.47,b:1.3)e:1.1,c:0.14)f:0.8,d:1.1)r:0.2;"];
Q = [{1:1}]
charMatrices = [{'a':[1],'b':[1],'c':[1]},{'a':['?'],'b':['?'],'c':[1],'d':[0]}]
nus = [0,0.22]
phis = [0,0.2]

test_id = 1

for T,charMtrx,nu,phi in zip(trees,charMatrices,nus,phis):
    with open("laml_unit_tests/test_data/test_Count_model/test_PMMN/test_posterior/test_" + str(test_id) + ".txt",'w') as fout:
        mySolver = EM_solver([T],{'charMtrx':charMtrx},{'Q':Q},{'phi':phi,'nu':nu})
        mySolver.az_partition()
        mySolver.Estep()
        for node in mySolver.trees[0].traverse_preorder():
            if node.post0[0] != -float("inf"):
                fout.write(node.label + ' post_node_z '+ str(node.post0[0]) + "\n")
            if node.post1[0] != -float("inf"):
                fout.write(node.label + ' post_node_s '+ str(node.post1[0]) + "\n") 
            if node.S0[0] != 0:
                fout.write(node.label + ' post_edge_zz ' + str(log(node.S0[0])) + "\n") 
            if node.S1[0] != 0:
                fout.write(node.label + ' post_edge_za ' + str(log(node.S1[0])) + "\n")
            if node.S2[0] != 0:
                fout.write(node.label + ' post_edge_zs ' + str(log(node.S2[0])) + "\n")
            if node.S3[0] != 0:
                fout.write(node.label + ' post_edge_aa ' + str(log(node.S3[0])) + "\n")
            if node.S4[0] != 0:
                fout.write(node.label + ' post_edge_as ' + str(log(node.S4[0])) + "\n")
    test_id += 1        
