import unittest
from laml_libs.IO_handler.sequence_lib import read_sequences
from laml_libs.Count_model.PMM_base import PMM_model, Alphabet,AlleleTable
from treeswift import *
from math import log
from random import random
from laml_libs import DEFAULT_STRATEGY
from copy import deepcopy
import pkg_resources

class PMMTest2(unittest.TestCase):
    # test Estep_out_llh
    def __get_reduced_trees(self,tree_str):
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

    def __countgen(self,alphabet,chosen_state,maxcount=1000):
        M = len(alphabet)
        counts = [0]*M
        is_missing = True
        for s in chosen_state:
            if s != '?':
                is_missing = False
                break
        if not is_missing:
            for i in range(M):
                counts[i] = int(random()*maxcount)
            m = max(counts) + int(maxcount/M)
        C = {}    
        for i,a in enumerate(alphabet):    
            C[a] = counts[i]
            if not is_missing and a == chosen_state:
                C[a] = m
        return C        
    
    def __charMtrx_2_alleleTable(self,charMtrx,alphabet):
        K = alphabet.K
        J = alphabet.J
        data_struct = {}
        for cell_name in charMtrx:
            counts = [{}]*K
            for k in range(K):
                counts[k] = self.__countgen(alphabet.get_cassette_alphabet(k),tuple([charMtrx[cell_name][k]]))
            data_struct[cell_name] = counts
        allele_table = AlleleTable(K,J,data_struct,alphabet)
        return allele_table

    # computing the outside likelihood
    def __test_outllh(self,T,Q,charMtrx,phi,nu,test_no,give_label=False):
        # generic function to test the computed out0 and out1 after calling Estep_out_llh
        K = len(list(charMtrx.values())[0])
        J = 1
        alphabet = Alphabet(K,J,[[[0,-1]+list(Q[k][0].keys())] for k in range(K)])
        allele_table = self.__charMtrx_2_alleleTable(charMtrx,alphabet)
        myModel = PMM_model([T],{'alleleTable':allele_table},{'Q':Q},1,nu,phi)
        if give_label:
            currIdx = 0
            for node in myModel.trees[0].traverse_preorder():
                if not node.is_leaf():
                    node.label = "I" + str(currIdx)
                    currIdx += 1
        myModel.Estep_in_llh()
        myModel.Estep_out_llh()
        out_llh = {}
        for node in myModel.trees[0].traverse_postorder():
            out_llh[node.label] = node.out_llh
        tree_reduced = self.__get_reduced_trees(myModel.trees[0].newick())
        for x in tree_reduced:    
            charMtrx0 = {y:charMtrx[y] for y in charMtrx}
            charMtrx0[x]= [0]*K
            allele_table0 = self.__charMtrx_2_alleleTable(charMtrx0,alphabet)
            tree_str = tree_reduced[x]
            myModel0 = PMM_model([tree_str],{'alleleTable':allele_table0},{'Q':Q},1,nu,phi)
            myModel0.Estep_in_llh()
            for k in range(K):
                true = myModel0.trees[0].root.in_llh[k][tuple([0])]
                est = out_llh[x][k][tuple([0])]
                self.assertAlmostEqual(true,est+log(1-phi),places=5,msg="PMMTest out llh: test_" + str(test_no) + " failed.")                
    def test_1(self):
        T = "(((a:1,b:1)e:1,(c:1,d:1)f:1)g:1)r;"
        Q = [[{1:1.0}]]
        charMtrx = {'a':[1],'b':[1],'c':[1],'d':[1]}
        phi = 0
        nu = 0.5
        self.__test_outllh(T,Q,charMtrx,phi,nu,1)
   
    def test_2(self):
        T = "(((a:1,b:1)e:1,(c:1,d:1)f:1)g:1)r;"
        Q = [[{1:1.0}]]
        charMtrx = {'a':[1],'b':[1],'c':[1],'d':[1]}
        nu = 0.1
        phi = 0
        self.__test_outllh(T,Q,charMtrx,phi,nu,2)
    
    def test_3(self):
        T = "(((a:1,b:1)e:1,(c:1,d:1)f:1)g:1)r;"
        Q = [[{1:1.0}]]
        charMtrx = {'a':[0],'b':[0],'c':[0],'d':[0]}
        nu = 0.25
        phi = 0
        self.__test_outllh(T,Q,charMtrx,phi,nu,3)
    
    def test_4(self):
        T = "(((a:0.1,b:1)e:1,(c:0.1,d:1)f:1)g:1)r;"
        Q = [[{1:1.0}]]
        charMtrx = {'a':[0],'b':[1],'c':[0],'d':[1]}
        nu = 0.15
        phi = 0
        self.__test_outllh(T,Q,charMtrx,phi,nu,4)
    
    def test_5(self):
        T = "(((a:0.1,b:1)e:1,(c:0.1,d:1)f:1)g:1)r;"
        Q = [[{1:1.0}]]
        charMtrx = {'a':[1],'b':[0],'c':[1],'d':[0]}
        nu = 0.1
        phi = 0
        self.__test_outllh(T,Q,charMtrx,phi,nu,5)
        
    def test_6(self):
        T = "(((a:0.1,b:1)e:1,(c:0.1,d:1)f:0.7)g:0.3)r;"
        Q = [[{1:1.0}]]
        charMtrx = {'a':[1],'b':[0],'c':[1],'d':[0]}
        nu = 0.19
        phi = 0
        self.__test_outllh(T,Q,charMtrx,phi,nu,6)
    
    def test_7(self):
        T = "(((a:0.1,b:1)e:1,(c:0.1,d:1)f:0.7)g:0.3)r;"
        Q = [[{1:1.0}]]
        charMtrx = {'a':['?'],'b':[0],'c':[1],'d':[0]}
        nu = 0.39
        phi = 0
        self.__test_outllh(T,Q,charMtrx,phi,nu,7)
    
    def test_8(self):
        T = "(((a:0.1,b:1)e:1,(c:0.1,d:1)f:0.7)g:0.3)r;"
        Q = [[{1:1.0}]]
        charMtrx = {'a':['?'],'b':[0],'c':[1],'d':['?']}
        nu = 0.39
        phi = 0
        self.__test_outllh(T,Q,charMtrx,phi,nu,8)
    
    def test_9(self):
        T = "(((a:0.1,b:1)e:1,(c:0.1,d:1)f:0.7)g:0.3)r;"
        Q = [[{1:0.5,2:0.5}]]
        charMtrx = {'a':['?'],'b':[2],'c':[1],'d':['?']}
        nu = 0.3
        phi = 0
        self.__test_outllh(T,Q,charMtrx,phi,nu,9)
    
    def test_10(self):
        T = "(((a:0.47,b:1.3)e:1.1,(c:0.14,d:1.1)f:0.72)g:0.39)r;"
        Q = [[{1:0.5,2:0.3,3:0.2}]]
        charMtrx = {'a':['?'],'b':[2],'c':[1],'d':['?']}
        nu = 0.22
        phi = 0.3
        self.__test_outllh(T,Q,charMtrx,phi,nu,30)
    
    def test_11(self):
        T = "((((a:0.47,b:1.3)e:1.1,c:0.14)f:0.8,d:1.1)g:0.2)r;"
        Q = [[{1:0.5,2:0.3,3:0.2}]]
        charMtrx = {'a':[1],'b':[1],'c':[1],'d':[1]}
        nu = 0.22
        phi = 0.01
        self.__test_outllh(T,Q,charMtrx,phi,nu,31)

    def test_12(self):
        T = "((((a:0.47,b:1.3)e:1.1,c:0.14)f:0.8,d:1.1)g:0.2)r;"
        Q = [[{1:0.5,2:0.3,3:0.2}]]
        charMtrx = {'a':['?'],'b':['?'],'c':[1],'d':[0]}
        nu = 0.22
        phi = 0.2
        self.__test_outllh(T,Q,charMtrx,phi,nu,32)
    
    def test_13(self):
        T = "((((a:0.47,b:1.3)e:1.1,c:0.14)f:0.8,d:1.1)g:0.2)r;"
        Q = [[{1:0.5,2:0.3,3:0.2}]]
        charMtrx = {'a':['?'],'b':['?'],'c':['?'],'d':[0]}
        nu = 0.22
        phi = 0.5
        self.__test_outllh(T,Q,charMtrx,phi,nu,32)
    
    def test_14(self):
        T = "((((a:0.47,b:1.3)e:1.1,c:0.14)f:0.8,d:1.1)g:0.2)r;"
        Q = [[{1:0.5,2:0.3,3:0.2}]]
        charMtrx = {'a':['?'],'b':['?'],'c':['?'],'d':[1]}
        nu = 0.22
        phi = 0.01
        self.__test_outllh(T,Q,charMtrx,phi,nu,33)
    
    def test_15(self):
        T = "((((a:0.47,b:1.3)e:1.1,c:0.14)f:0.8,d:1.1)g:0.2)r;"
        Q = [[{1:0.5,2:0.3,3:0.2}]]
        charMtrx = {'a':['?'],'b':['?'],'c':['?'],'d':['?']}
        nu = 0.22
        phi = 0.9
        self.__test_outllh(T,Q,charMtrx,phi,nu,34)
    
    def test_16(self):
        T = "((((a:0.47,b:1.3)e:1.1,c:0.14)f:0.8,d:1.1)g:0.2)r;"
        Q = [[{1:0.5,2:0.3,3:0.2}]]
        charMtrx = {'a':[1],'b':[2],'c':[3],'d':['?']}
        nu = 0.22
        phi = 0.03
        self.__test_outllh(T,Q,charMtrx,phi,nu,35)
    
    def test_17(self):
        T = "((((a:0.47,b:1.3)f:1.1,c:0.14)g:0.8,(d:1.1,e:0.2)h:0.2)g:0.01)r;"
        Q = [[{1:0.5,2:0.3,3:0.2}]]
        charMtrx = {'a':[1],'b':[2],'c':[3],'d':['?'],'e':['?']}
        nu = 0.22
        phi = 0.5
        self.__test_outllh(T,Q,charMtrx,phi,nu,36)
    
    def test_18(self):
        T = "((((a:0.47,b:1.3)f:1.1,c:0.14)g:0.8,(d:1.1,e:0.2)h:0.2)g:0.01)r;"
        Q = [[{1:0.5,2:0.3,3:0.2}]]
        charMtrx = {'a':['?'],'b':['?'],'c':['?'],'d':['?'],'e':['?']}
        nu = 0.26
        phi = 0.7
        self.__test_outllh(T,Q,charMtrx,phi,nu,37)
    
    def test_19(self):
        T = "((((a:0.47,b:1.3)f:1.1,c:0.14)g:0.8,(d:1.1,e:0.2)h:0.2)g:0.01)r;"
        Q = [[{1:0.2,2:0.6,3:0.2}]]
        charMtrx = {'a':['?'],'b':['?'],'c':['?'],'d':[2],'e':[2]}
        nu = 0.26
        phi = 0.1
        self.__test_outllh(T,Q,charMtrx,phi,nu,38)
    
    def test_20(self):
        T = "((((a:0.47,b:1.3)f:1.1,c:0.14)g:0.8,(d:1.1,e:0.2)h:0.2)g:0.01)r;"
        Q = [[{1:0.2,2:0.6,3:0.2}]]
        charMtrx = {'a':['?'],'b':['?'],'c':[2],'d':[2],'e':[2]}
        nu = 0.52
        phi = 0.8
        self.__test_outllh(T,Q,charMtrx,phi,nu,39)
    
    def test_21(self):
        T = "((((a:0.47,b:1.3)f:1.1,c:0.14)g:0.8,(d:1.1,e:0.2)h:0.2)g:0.01)r;"
        Q = [[{1:0.2,2:0.6,3:0.2}]]
        charMtrx = {'a':['?'],'b':['?'],'c':[2],'d':[1],'e':[2]}
        nu = 0.52
        phi = 0.001
        self.__test_outllh(T,Q,charMtrx,phi,nu,40)
    
    def test_22(self):
        T = "(((((a:0.47,b:1.3)f:1.1,c:0.14)g:0.8,d:1.1)h:0.2,e:0.2)g:0.2)r;"
        Q = [[{1:0.2,2:0.6,3:0.2}]]
        charMtrx = {'a':['?'],'b':['?'],'c':[2],'d':['?'],'e':[2]}
        nu = 0.2
        phi = 0.82
        self.__test_outllh(T,Q,charMtrx,phi,nu,41)
    
    def test_23(self):
        T = "(((((a:0.47,b:1.3)f:1.1,c:0.14)g:0.8,d:1.1)h:0.2,e:0.2)g:0.2)r;"
        Q = [[{1:0.2,2:0.6,3:0.2}]]
        charMtrx = {'a':['?'],'b':['?'],'c':[2],'d':['?'],'e':[0]}
        nu = 0.2
        phi = 0.1
        self.__test_outllh(T,Q,charMtrx,phi,nu,42)
    
    def test_24(self):
        T = "(((((a:0.47,b:1.3)f:1.1,c:0.14)g:0.8,d:1.1)h:0.2,e:0.2)g:0.2)r;"
        Q = [[{1:0.2,2:0.6,3:0.2}]]
        charMtrx = {'a':['?'],'b':[0],'c':[2],'d':['?'],'e':[0]}
        nu = 0.2
        phi = 0.3
        self.__test_outllh(T,Q,charMtrx,phi,nu,43)
    
    def test_25(self):
        T = "(((((a:0.47,b:1.3)f:1.1,c:0.14)g:0.8,d:1.1)h:0.2,e:0.2)g:0.2)r;"
        Q = [[{1:0.2,2:0.6,3:0.2}]]
        charMtrx = {'a':['?'],'b':[0],'c':[2],'d':['?'],'e':[1]}
        nu = 0.2
        phi = 0.3
        self.__test_outllh(T,Q,charMtrx,phi,nu,44)
    
    def test_26(self):
        #treedata_path = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_EM/test1.tre')
        treedata_path = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_Count_model/test_PMM_base/test1_n25.tre')
        #charMtrx_path = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_EM/test1_charMtrx.txt')
        charMtrx_path = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_Count_model/test_PMM_base/test1_charMtrx.txt')
        T = read_tree_newick(treedata_path)
        phi = 0.05231954386883335
        nu = 0.15877477685098262
        charMtrx,_ = read_sequences(charMtrx_path,filetype="charMtrx",delimiter=",",masked_symbol='?',suppress_warnings=True)
        k = 60
        Q = [[] for _ in range(k)]
        for i in range(k):
            M_i = set(charMtrx[x][i] for x in charMtrx if charMtrx[x][i] not in [0,"?"])
            m_i = len(M_i)
            q = {x:1/m_i for x in M_i}
            #q[0] = 0
            Q[i].append(q)
        self.__test_outllh(T,Q,charMtrx,phi,nu,45,give_label=True)
