import os 
import unittest
from laml_libs.Count_model.PMM_base import *
from treeswift import *
from laml_libs.IO_handler.sequence_lib import read_sequences
from random import random

class PMMTest(unittest.TestCase):
    def __countgen(self,alphabet,chosen_state,maxcount=1000):
        M = len(alphabet)
        counts = [0]*M
        for i in range(M):
            counts[i] = int(random()*maxcount)
        m = max(counts) + int(maxcount/M)
        C = {}    
        for i,a in enumerate(alphabet):    
            C[a] = counts[i]
            if a == chosen_state:
                C[a] = m
        return C        

    def test_1(self): 
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1]]])
        #allele_table = AlleleTable(K,J,{'a':[{(0,):0,(1,):1}],'b':[{(0,):0,(1,):1}],'c':[{(0,):0,(1,):1}]},alphabet)
        counts_a = self.__countgen([(0,),(1,)],(1,))
        counts_b = self.__countgen([(0,),(1,)],(1,))
        counts_c = self.__countgen([(0,),(1,)],(1,))
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b],'c':[counts_c]},alphabet)

        true_nllh = 0.20665578828621584

        myModel = PMM_model([T],{'alleleTable':allele_table},{'Q':Q},1,0,0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest: test_1 failed.")
    
    def test_2(self): 
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1]]])
        #allele_table = AlleleTable(K,J,{'a':[{(0,):0,(1,):1}],'b':[{(0,):0,(1,):1}],'c':[{(0,):1,(1,):0}]},alphabet)
        counts_a = self.__countgen([(0,),(1,)],(1,))
        counts_b = self.__countgen([(0,),(1,)],(1,))
        counts_c = self.__countgen([(0,),(1,)],(0,))        
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b],'c':[counts_c]},alphabet)

        true_nllh = 2.2495946917551692 

        myModel = PMM_model([T],{'alleleTable':allele_table},{'Q':Q},1,0,0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest: test_2 failed.")
   
    def test_3(self): 
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1]]])
        #allele_table = AlleleTable(K,J,{'a':[{(0,):0,(1,):1}],'b':[{(0,):1,(1,):0}],'c':[{(0,):0,(1,):1}]},alphabet)
        counts_a = self.__countgen([(0,),(1,)],(1,))
        counts_b = self.__countgen([(0,),(1,)],(0,))
        counts_c = self.__countgen([(0,),(1,)],(1,))        
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b],'c':[counts_c]},alphabet)
        true_nllh = 3.917350291274164 

        myModel = PMM_model([T],{'alleleTable':allele_table},{'Q':Q},1,0,0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest: test_3 failed.")
    
    def test_4(self): 
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1]]])
        #allele_table = AlleleTable(K,J,{'a':[{(0,):0,(1,):1}],'b':[{(0,):1,(1,):0}],'c':[{(0,):0,(1,):1}]},alphabet)
        counts_a = self.__countgen([(0,),(1,)],(0,))
        counts_b = self.__countgen([(0,),(1,)],(1,))
        counts_c = self.__countgen([(0,),(1,)],(1,))        
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b],'c':[counts_c]},alphabet)
        true_nllh = 3.917350291274164 

        myModel = PMM_model([T],{'alleleTable':allele_table},{'Q':Q},1,0,0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest: test_4 failed.")
    
    '''def test_5(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[1],'b':[0],'c':[0]}
        true_nllh = 4.4586751457870815

        #mySolver = ML_solver(msa,Q,T,phi=0,nu=0)
        mySolver = ML_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0,'nu':0})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_5 failed.")
    
    def test_6(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[0],'b':[1],'c':[0]}
        true_nllh = 4.4586751457870815

        #mySolver = ML_solver(msa,Q,T,phi=0,nu=0)
        mySolver = ML_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0,'nu':0})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_6 failed.")
    
    def test_7(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[0],'b':[0],'c':[1]}
        true_nllh = 4.4586751457870815

        #mySolver = ML_solver(msa,Q,T,phi=0,nu=0)
        mySolver = ML_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0,'nu':0})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_7 failed.")
    
    def test_8(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[0],'b':[0],'c':[0]}
        true_nllh = 5.0

        #mySolver = ML_solver(msa,Q,T,phi=0,nu=0)
        mySolver = ML_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0,'nu':0})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_8 failed.")
    
    def test_9(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[0],'b':[0],'c':['?']}
        true_nllh = 6.513306124309698

        #mySolver = ML_solver(msa,Q,T,phi=0.1,nu=0)
        mySolver = ML_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':0})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_9 failed.")
    
    def test_10(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[0],'b':['?'],'c':[0]}
        true_nllh = 6.513306124309698

        #mySolver = ML_solver(msa,Q,T,phi=0.1,nu=0)
        mySolver = ML_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':0})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_10 failed.")
    
    def test_11(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':['?'],'b':[0],'c':[0]}
        true_nllh = 6.513306124309698

        #mySolver = ML_solver(msa,Q,T,phi=0.1,nu=0)
        mySolver = ML_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':0})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_11 failed.")
    
    def test_12(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[0],'b':[1],'c':['?']}
        true_nllh = 5.97198126969678

        #mySolver = ML_solver(msa,Q,T,phi=0.1,nu=0)
        mySolver = ML_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':0})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_12 failed.")
    
    def test_13(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[0],'b':['?'],'c':[1]}
        true_nllh = 5.97198126969678

        #mySolver = ML_solver(msa,Q,T,phi=0.1,nu=0)
        mySolver = ML_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':0})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_13 failed.")
    
    def test_14(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':['?'],'b':[0],'c':[1]}
        true_nllh = 5.97198126969678

        #mySolver = ML_solver(msa,Q,T,phi=0.1,nu=0)
        mySolver = ML_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':0})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_14 failed.")
    
    def test_15(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[1],'b':['?'],'c':[0]}
        true_nllh = 4.658719582178557

        #mySolver = ML_solver(msa,Q,T,phi=0.1,nu=0)
        mySolver = ML_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':0})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_15 failed.")
    
    def test_16(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':['?'],'b':[1],'c':[0]}
        true_nllh = 4.658719582178557

        #mySolver = ML_solver(msa,Q,T,phi=0.1,nu=0)
        mySolver = ML_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':0})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_16 failed.")
    
    def test_17(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[1],'b':[1],'c':['?']}
        true_nllh = 2.5980566021648364

        #mySolver = ML_solver(msa,Q,T,phi=0.1,nu=0)
        mySolver = ML_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':0})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_17 failed.")
    
    def test_18(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[1],'b':['?'],'c':[1]}
        true_nllh = 2.695795750497349

        #mySolver = ML_solver(msa,Q,T,phi=0.1,nu=0)
        mySolver = ML_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':0})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_18 failed.")
    
    def test_19(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':['?'],'b':[1],'c':[1]}
        true_nllh = 2.695795750497349

        #mySolver = ML_solver(msa,Q,T,phi=0.1,nu=0)
        mySolver = ML_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':0})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_19 failed.")
    
    def test_20(self): 
        Q = [{1:0.5,2:0.5}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[1],'b':[1],'c':[1]}
        true_nllh = 1.0297894223949402

        mySolver = ML_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0,'nu':0})
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_20 failed.")

    def test_21(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):0,c:1):1;"
        msa = {'a':[1],'b':[1],'c':[1]}
        true_nllh = 0.3215288449416738

        mySolver = ML_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0,'nu':0})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_21 failed.") '''
