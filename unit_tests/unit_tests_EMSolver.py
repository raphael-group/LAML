import unittest
from scmail_libs.sequence_lib import read_sequences
from scmail_libs.EM_solver import EM_solver
from scmail_libs.ML_solver import ML_solver
from treeswift import *
from math import log
from random import random
from scmail_libs import DEFAULT_STRATEGY
from copy import deepcopy

class EMTest(unittest.TestCase):
    # test likelihood computation
    def test_1(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[1],'b':[1],'c':[1]}
        true_nllh = 0.20665578828621584

        mySolver = EM_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0,'nu':0})
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="EMTest: test_1 failed.")
    
    def test_2(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[1],'b':[1],'c':[0]}
        true_nllh = 2.2495946917551692 

        mySolver = EM_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0,'nu':0})
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="EMTest: test_2 failed.")
    
    def test_3(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[1],'b':[0],'c':[1]}
        true_nllh = 3.917350291274164 

        mySolver = EM_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0,'nu':0})
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="EMTest: test_3 failed.")
    
    def test_4(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[0],'b':[1],'c':[1]}
        true_nllh = 3.917350291274164 

        mySolver = EM_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0,'nu':0})
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="EMTest: test_4 failed.")
    
    def test_5(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[1],'b':[0],'c':[0]}
        true_nllh = 4.4586751457870815

        mySolver = EM_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0,'nu':0})
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="EMTest: test_5 failed.")
    
    def test_6(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[0],'b':[1],'c':[0]}
        true_nllh = 4.4586751457870815

        mySolver = EM_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0,'nu':0})
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="EMTest: test_6 failed.")
    
    def test_7(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[0],'b':[0],'c':[1]}
        true_nllh = 4.4586751457870815

        mySolver = EM_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0,'nu':0})
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="EMTest: test_7 failed.")
    
    def test_8(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[0],'b':[0],'c':[0]}
        true_nllh = 5.0

        
        mySolver = EM_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0,'nu':0})
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="EMTest: test_8 failed.")
    
    def test_9(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[0],'b':[0],'c':['?']}
        true_nllh = 6.513306124309698

        
        mySolver = EM_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':0})
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="EMTest: test_9 failed.")
    
    def test_10(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[0],'b':['?'],'c':[0]}
        true_nllh = 6.513306124309698

        
        mySolver = EM_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':0})
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="EMTest: test_10 failed.")
    
    def test_11(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':['?'],'b':[0],'c':[0]}
        true_nllh = 6.513306124309698

        
        mySolver = EM_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':0})
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="EMTest: test_11 failed.")
    
    def test_12(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[0],'b':[1],'c':['?']}
        true_nllh = 5.97198126969678

        
        mySolver = EM_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':0})
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="EMTest: test_12 failed.")
    
    def test_13(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[0],'b':['?'],'c':[1]}
        true_nllh = 5.97198126969678

        
        mySolver = EM_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':0})
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="EMTest: test_13 failed.")
    
    def test_14(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':['?'],'b':[0],'c':[1]}
        true_nllh = 5.97198126969678

        
        mySolver = EM_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':0})
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="EMTest: test_14 failed.")
    
    def test_15(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[1],'b':['?'],'c':[0]}
        true_nllh = 4.658719582178557

        
        mySolver = EM_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':0})
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="EMTest: test_15 failed.")
    
    def test_16(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':['?'],'b':[1],'c':[0]}
        true_nllh = 4.658719582178557

        
        mySolver = EM_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':0})
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="EMTest: test_16 failed.")
    
    def test_17(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[1],'b':[1],'c':['?']}
        true_nllh = 2.5980566021648364

        
        mySolver = EM_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':0})
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="EMTest: test_17 failed.")
    
    def test_18(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[1],'b':['?'],'c':[1]}
        true_nllh = 2.695795750497349

        
        mySolver = EM_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':0})
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="EMTest: test_18 failed.")
    
    def test_19(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':['?'],'b':[1],'c':[1]}
        true_nllh = 2.695795750497349

        
        mySolver = EM_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':0})
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="EMTest: test_19 failed.")
    
    def test_20(self): 
        Q = [{1:0.5,2:0.5}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[1],'b':[1],'c':[1]}
        true_nllh = 1.0297894223949402

        
        mySolver = EM_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0,'nu':0})
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="EMTest: test_20 failed.")
    # test Estep_out_llh
    def __get_reduced_trees__(self,tree_str):
        tree_obj = read_tree_newick(tree_str)
        tree_reduced = {}
        for node in tree_obj.traverse_leaves():
            tree_reduced[node.label] = tree_str
        node_list = [node for node in tree_obj.traverse_postorder() if not node.is_root() and not node.is_leaf()]
        for node in node_list:
            u,v = node.children
            u_len = u.edge_length
            v_len = v.edge_length
            node.remove_child(u)
            node.remove_child(v)
            tree_reduced[node.label] = tree_obj.newick()
            node.add_child(u)
            node.add_child(v)
            u.edge_length = u_len
            v.edge_length = v_len
        return tree_reduced    

    # computing the outside likelihood
    def __test_outllh__(self,T,Q,msa,phi,nu,test_no,give_label=False):
        # generic function to test the computed out0 and out1 after calling Estep_out_llh
        mySolver = EM_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':phi,'nu':nu})
        if give_label:
            currIdx = 0
            for node in mySolver.trees[0].traverse_preorder():
                if not node.is_leaf():
                    node.label = "I" + str(currIdx)
                    currIdx += 1
        mySolver.az_partition()
        mySolver.Estep_in_llh()
        mySolver.Estep_out_llh()
        out0 = {} # mapping node label to node.out0
        out1 = {} # mapping node label to node.out1
        
        for node in mySolver.trees[0].traverse_postorder():
            out0[node.label] = node.out0
            out1[node.label] = node.out1
        tree_reduced = self.__get_reduced_trees__(mySolver.trees[0].newick())
        for x in tree_reduced:    
            # test out0
            msa0 = {y:msa[y] for y in msa}
            msa0[x]= [0]*mySolver.numsites
            tree_str = tree_reduced[x]
            mySolver0 = EM_solver([tree_str],{'charMtrx':msa0},{'Q':Q},{'phi':phi,'nu':nu})
            mySolver0.az_partition()
            mySolver0.Estep_in_llh()
            for true,est in zip(mySolver0.trees[0].root.L0,out0[x]):
                self.assertAlmostEqual(true,est+log(1-phi),places=5,msg="EMTest: test_" + str(test_no) + " failed.")
            # test out1            
            msa1 = {y:msa[y] for y in msa}
            msa1[x] = [-1]*mySolver.numsites
            mySolver1 = EM_solver([tree_str],{'charMtrx':msa1},{'Q':Q},{'phi':phi,'nu':nu})
            mySolver1.az_partition()
            mySolver1.Estep_in_llh()
            for true,est in zip(mySolver1.trees[0].root.L0,out1[x]):
                self.assertAlmostEqual(true,est,places=5,msg="EMTest: test_" + str(test_no) + " failed.")

    def test_21(self):
        T = "((a:1,b:1)e:1,(c:1,d:1)f:1)r:1;"
        Q = [{1:1.0}]
        msa = {'a':[1],'b':[1],'c':[1],'d':[1]}
        phi = 0
        nu = 0.5
        self.__test_outllh__(T,Q,msa,phi,nu,21)
    
    def test_22(self):
        T = "((a:1,b:1)e:1,(c:1,d:1)f:1)r:1;"
        Q = [{1:1.0}]
        msa = {'a':[1],'b':[1],'c':[1],'d':[1]}
        nu = 0.1
        phi = 0
        self.__test_outllh__(T,Q,msa,phi,nu,22)

    def test_23(self):
        T = "((a:1,b:1)e:1,(c:1,d:1)f:1)r:1;"
        Q = [{1:1.0}]
        msa = {'a':[0],'b':[0],'c':[0],'d':[0]}
        nu = 0.25
        phi = 0
        self.__test_outllh__(T,Q,msa,phi,nu,23)
    
    def test_24(self):
        T = "((a:0.1,b:1)e:1,(c:0.1,d:1)f:1)r:1;"
        Q = [{1:1.0}]
        msa = {'a':[0],'b':[1],'c':[0],'d':[1]}
        nu = 0.15
        phi = 0
        self.__test_outllh__(T,Q,msa,phi,nu,24)
    
    def test_25(self):
        T = "((a:0.1,b:1)e:1,(c:0.1,d:1)f:1)r:1;"
        Q = [{1:1.0}]
        msa = {'a':[1],'b':[0],'c':[1],'d':[0]}
        nu = 0.1
        phi = 0
        self.__test_outllh__(T,Q,msa,phi,nu,25)
        
    def test_26(self):
        T = "((a:0.1,b:1)e:1,(c:0.1,d:1)f:0.7)r:0.3;"
        Q = [{1:1.0}]
        msa = {'a':[1],'b':[0],'c':[1],'d':[0]}
        nu = 0.19
        phi = 0
        self.__test_outllh__(T,Q,msa,phi,nu,26)
    
    def test_27(self):
        T = "((a:0.1,b:1)e:1,(c:0.1,d:1)f:0.7)r:0.3;"
        Q = [{1:1.0}]
        msa = {'a':['?'],'b':[0],'c':[1],'d':[0]}
        nu = 0.39
        phi = 0
        self.__test_outllh__(T,Q,msa,phi,nu,27)
    
    def test_28(self):
        T = "((a:0.1,b:1)e:1,(c:0.1,d:1)f:0.7)r:0.3;"
        Q = [{1:1.0}]
        msa = {'a':['?'],'b':[0],'c':[1],'d':['?']}
        nu = 0.39
        phi = 0
        self.__test_outllh__(T,Q,msa,phi,nu,28)
    
    def test_29(self):
        T = "((a:0.1,b:1)e:1,(c:0.1,d:1)f:0.7)r:0.3;"
        Q = [{1:0.5,2:0.5}]
        msa = {'a':['?'],'b':[2],'c':[1],'d':['?']}
        nu = 0.3
        phi = 0
        self.__test_outllh__(T,Q,msa,phi,nu,29)
    
    def test_30(self):
        T = "((a:0.47,b:1.3)e:1.1,(c:0.14,d:1.1)f:0.72)r:0.39;"
        Q = [{1:0.5,2:0.3,3:0.2}]
        msa = {'a':['?'],'b':[2],'c':[1],'d':['?']}
        nu = 0.22
        phi = 0.3
        self.__test_outllh__(T,Q,msa,phi,nu,30)
    
    def test_31(self):
        T = "(((a:0.47,b:1.3)e:1.1,c:0.14)f:0.8,d:1.1)r:0.2;"
        Q = [{1:0.5,2:0.3,3:0.2}]
        msa = {'a':[1],'b':[1],'c':[1],'d':[1]}
        nu = 0.22
        phi = 0.01
        self.__test_outllh__(T,Q,msa,phi,nu,31)

    def test_32(self):
        T = "(((a:0.47,b:1.3)e:1.1,c:0.14)f:0.8,d:1.1)r:0.2;"
        Q = [{1:0.5,2:0.3,3:0.2}]
        msa = {'a':['?'],'b':['?'],'c':[1],'d':[0]}
        nu = 0.22
        phi = 0.2
        self.__test_outllh__(T,Q,msa,phi,nu,32)
    
    def test_32(self):
        T = "(((a:0.47,b:1.3)e:1.1,c:0.14)f:0.8,d:1.1)r:0.2;"
        Q = [{1:0.5,2:0.3,3:0.2}]
        msa = {'a':['?'],'b':['?'],'c':['?'],'d':[0]}
        nu = 0.22
        phi = 0.5
        self.__test_outllh__(T,Q,msa,phi,nu,32)
    
    def test_33(self):
        T = "(((a:0.47,b:1.3)e:1.1,c:0.14)f:0.8,d:1.1)r:0.2;"
        Q = [{1:0.5,2:0.3,3:0.2}]
        msa = {'a':['?'],'b':['?'],'c':['?'],'d':[1]}
        nu = 0.22
        phi = 0.01
        self.__test_outllh__(T,Q,msa,phi,nu,33)
    
    def test_34(self):
        T = "(((a:0.47,b:1.3)e:1.1,c:0.14)f:0.8,d:1.1)r:0.2;"
        Q = [{1:0.5,2:0.3,3:0.2}]
        msa = {'a':['?'],'b':['?'],'c':['?'],'d':['?']}
        nu = 0.22
        phi = 0.9
        self.__test_outllh__(T,Q,msa,phi,nu,34)
    
    def test_35(self):
        T = "(((a:0.47,b:1.3)e:1.1,c:0.14)f:0.8,d:1.1)r:0.2;"
        Q = [{1:0.5,2:0.3,3:0.2}]
        msa = {'a':[1],'b':[2],'c':[3],'d':['?']}
        nu = 0.22
        phi = 0.03
        self.__test_outllh__(T,Q,msa,phi,nu,35)
    
    def test_36(self):
        T = "(((a:0.47,b:1.3)f:1.1,c:0.14)g:0.8,(d:1.1,e:0.2)h:0.2)r:0.01;"
        Q = [{1:0.5,2:0.3,3:0.2}]
        msa = {'a':[1],'b':[2],'c':[3],'d':['?'],'e':['?']}
        nu = 0.22
        phi = 0.5
        self.__test_outllh__(T,Q,msa,phi,nu,36)
    
    def test_37(self):
        T = "(((a:0.47,b:1.3)f:1.1,c:0.14)g:0.8,(d:1.1,e:0.2)h:0.2)r:0.01;"
        Q = [{1:0.5,2:0.3,3:0.2}]
        msa = {'a':['?'],'b':['?'],'c':['?'],'d':['?'],'e':['?']}
        nu = 0.26
        phi = 0.7
        self.__test_outllh__(T,Q,msa,phi,nu,37)
    
    def test_38(self):
        T = "(((a:0.47,b:1.3)f:1.1,c:0.14)g:0.8,(d:1.1,e:0.2)h:0.2)r:0.01;"
        Q = [{1:0.2,2:0.6,3:0.2}]
        msa = {'a':['?'],'b':['?'],'c':['?'],'d':[2],'e':[2]}
        nu = 0.26
        phi = 0.1
        self.__test_outllh__(T,Q,msa,phi,nu,38)
    
    def test_39(self):
        T = "(((a:0.47,b:1.3)f:1.1,c:0.14)g:0.8,(d:1.1,e:0.2)h:0.2)r:0.01;"
        Q = [{1:0.2,2:0.6,3:0.2}]
        msa = {'a':['?'],'b':['?'],'c':[2],'d':[2],'e':[2]}
        nu = 0.52
        phi = 0.8
        self.__test_outllh__(T,Q,msa,phi,nu,39)
    
    def test_40(self):
        T = "(((a:0.47,b:1.3)f:1.1,c:0.14)g:0.8,(d:1.1,e:0.2)h:0.2)r:0.01;"
        Q = [{1:0.2,2:0.6,3:0.2}]
        msa = {'a':['?'],'b':['?'],'c':[2],'d':[1],'e':[2]}
        nu = 0.52
        phi = 0.001
        self.__test_outllh__(T,Q,msa,phi,nu,40)
    
    def test_41(self):
        T = "((((a:0.47,b:1.3)f:1.1,c:0.14)g:0.8,d:1.1)h:0.2,e:0.2)r:0.2;"
        Q = [{1:0.2,2:0.6,3:0.2}]
        msa = {'a':['?'],'b':['?'],'c':[2],'d':['?'],'e':[2]}
        nu = 0.2
        phi = 0.82
        self.__test_outllh__(T,Q,msa,phi,nu,41)
    
    def test_42(self):
        T = "((((a:0.47,b:1.3)f:1.1,c:0.14)g:0.8,d:1.1)h:0.2,e:0.2)r:0.2;"
        Q = [{1:0.2,2:0.6,3:0.2}]
        msa = {'a':['?'],'b':['?'],'c':[2],'d':['?'],'e':[0]}
        nu = 0.2
        phi = 0.1
        self.__test_outllh__(T,Q,msa,phi,nu,42)
    
    def test_43(self):
        T = "((((a:0.47,b:1.3)f:1.1,c:0.14)g:0.8,d:1.1)h:0.2,e:0.2)r:0.2;"
        Q = [{1:0.2,2:0.6,3:0.2}]
        msa = {'a':['?'],'b':[0],'c':[2],'d':['?'],'e':[0]}
        nu = 0.2
        phi = 0.3
        self.__test_outllh__(T,Q,msa,phi,nu,43)
    
    def test_44(self):
        T = "((((a:0.47,b:1.3)f:1.1,c:0.14)g:0.8,d:1.1)h:0.2,e:0.2)r:0.2;"
        Q = [{1:0.2,2:0.6,3:0.2}]
        msa = {'a':[-1],'b':[0],'c':[2],'d':[-1],'e':[1]}
        nu = 0.2
        phi = 0.3
        self.__test_outllh__(T,Q,msa,phi,nu,44)
    
    def test_45(self):
        T = read_tree_newick("unit_tests/test_data/test_EM/test1.tre")
        phi = 0.05231954386883335
        nu = 0.15877477685098262
        msa,_ = read_sequences("unit_tests/test_data/test_EM/test1_charMtrx.txt",filetype="charMtrx",delimiter=",",masked_symbol='-',suppress_warnings=True)
        Q = []
        k = 60
        for i in range(k):
            M_i = set(msa[x][i] for x in msa if msa[x][i] not in [0,"?"])
            m_i = len(M_i)
            q = {x:1/m_i for x in M_i}
            q[0] = 0
            Q.append(q)
        self.__test_outllh__(T,Q,msa,phi,nu,45,give_label=True)
    
    def test_46(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):0,c:1):1;"
        msa = {'a':[1],'b':[1],'c':[1]}
        true_nllh = 0.3215288449416738

        mySolver = EM_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0,'nu':0})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="EMTest: test_46 failed.")

    # test optimize_one
    def test_47(self):
        T = read_tree_newick("unit_tests/test_data/test_EM/test4.tre")
        msa,_ = read_sequences("unit_tests/test_data/test_EM/test4_charMtrx.txt",filetype="charMtrx",delimiter=",",masked_symbol='?')
        k = len(msa[next(iter(msa.keys()))])
        Q = [{0:0} for i in range(k)]
        seen_sites = set()
        with open("unit_tests/test_data/test_EM/test4_prior.csv",'r') as fin:
            lines = fin.readlines()
            for line in lines:
                site_idx,char_state,prob = line.strip().split(',')
                site_idx = int(site_idx)
                if site_idx not in seen_sites:
                    seen_sites.add(site_idx)
                char_state = int(char_state)
                prob = float(prob)
                Q[len(seen_sites) - 1][char_state] = prob
        true_nllh = 1007.009158767873 
        true_phi = 0
        true_nu = 0

        mySolver = EM_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0,'nu':0})
        randseed = 1221
        nllh,status = mySolver.optimize(initials=1,random_seeds=randseed,verbose=-1,ultra_constr=False)
        phi = mySolver.params.phi
        nu = mySolver.params.nu
        self.assertAlmostEqual(0,abs(true_nllh-nllh)/true_nllh,places=4,msg="EMTest: test_47 failed.")
        self.assertAlmostEqual(true_phi,phi,places=4,msg="EMTest: test_47 failed.")
        self.assertAlmostEqual(true_nu,nu,places=4,msg="EMTest: test_47 failed.")

    # test score_tree
    def test_48(self):
        T = read_tree_newick("unit_tests/test_data/test_EM/test4.tre")
        msa,_ = read_sequences("unit_tests/test_data/test_EM/test4_charMtrx.txt",filetype="charMtrx",delimiter=",",masked_symbol='?')
        k = len(msa[next(iter(msa.keys()))])
        Q = [{0:0} for i in range(k)]
        seen_sites = set()
        with open("unit_tests/test_data/test_EM/test4_prior.csv",'r') as fin:
            lines = fin.readlines()
            for line in lines:
                site_idx,char_state,prob = line.strip().split(',')
                site_idx = int(site_idx)
                if site_idx not in seen_sites:
                    seen_sites.add(site_idx)
                char_state = int(char_state)
                prob = float(prob)
                Q[len(seen_sites) - 1][char_state] = prob
        true_nllh = 1007.009158767873
        true_phi = 0
        true_nu = 0

        mySolver = EM_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0,'nu':0})
        randseed = 1221
        my_strategy = deepcopy(DEFAULT_STRATEGY) 
        my_strategy['fixed_brlen'] = {}
        score,status = mySolver.score_tree(strategy=my_strategy)
        nllh = -score
        phi = mySolver.params.phi
        nu = mySolver.params.nu
        self.assertAlmostEqual(abs(true_nllh-nllh)/true_nllh,0,places=4,msg="EMTest: test_48 failed.")
        self.assertAlmostEqual(true_phi,phi,places=4,msg="EMTest: test_48 failed.")
        self.assertAlmostEqual(true_nu,nu,places=4,msg="EMTest: test_48 failed.")
