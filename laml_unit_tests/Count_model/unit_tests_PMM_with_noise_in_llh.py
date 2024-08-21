import os 
import unittest
from laml_libs.Count_model.PMM_with_noise import *
from treeswift import *
from laml_libs.IO_handler.sequence_lib import read_sequences
from random import random
from laml_libs.Count_model.utils import *
from .virtual_unit_tests import VirtualUnitTest

class PMMNTest_in_llh(VirtualUnitTest):
    # test in_llh computation
    def test_1(self): 
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1]]])
        counts_a = countgen([(0,),(1,)],(1,))
        counts_b = countgen([(0,),(1,)],(1,))
        counts_c = countgen([(0,),(1,)],(1,))
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b],'c':[counts_c]},alphabet)

        true_nllh = 0.20665578828621584

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0,phi=0,eta=0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_1 failed.")
    
    def test_2(self): 
        # evaluate when cell C's argmax cassette state is not the true cassette state
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1]]])
        counts_a = countgen([(0,),(1,)],(1,))
        counts_b = countgen([(0,),(1,)],(1,))
        counts_c = countgen([(0,),(1,)],(0,))        
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b],'c':[counts_c]},alphabet)

        true_nllh = 2.2495946917551692 

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0,phi=0,eta=0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_2 failed.")
   
    def test_3(self): 
        # evaluate when cell B's argmax cassette state is not the true cassette state
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1]]])
        counts_a = countgen([(0,),(1,)],(1,))
        counts_b = countgen([(0,),(1,)],(0,))
        counts_c = countgen([(0,),(1,)],(1,))        
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b],'c':[counts_c]},alphabet)
        true_nllh = 3.917350291274164 

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0,phi=0,eta=0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_3 failed.")
    
    def test_4(self): 
        # evaluate when cell A's argmax cassette state is not the true cassette state
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1]]])
        counts_a = countgen([(0,),(1,)],(0,))
        counts_b = countgen([(0,),(1,)],(1,))
        counts_c = countgen([(0,),(1,)],(1,))        
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b],'c':[counts_c]},alphabet)
        true_nllh = 3.917350291274164 

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0,phi=0,eta=0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_4 failed.")
    
    def test_5(self): 
        # evaluate when cell B and C's argmax cassette state is not the true cassette state
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1]]])
        counts_a = countgen([(0,),(1,)],(1,))
        counts_b = countgen([(0,),(1,)],(0,))
        counts_c = countgen([(0,),(1,)],(0,))        
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b],'c':[counts_c]},alphabet)
        true_nllh = 4.4586751457870815

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0,phi=0,eta=0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_5 failed.")
    
    def test_6(self): 
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1]]])
        counts_a = countgen([(0,),(1,)],(0,))
        counts_b = countgen([(0,),(1,)],(1,))
        counts_c = countgen([(0,),(1,)],(0,))        
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b],'c':[counts_c]},alphabet)
        true_nllh = 4.4586751457870815

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0,phi=0,eta=0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_6 failed.")
    
    def test_7(self): 
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1]]])
        counts_a = countgen([(0,),(1,)],(0,))
        counts_b = countgen([(0,),(1,)],(0,))
        counts_c = countgen([(0,),(1,)],(1,))        
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b],'c':[counts_c]},alphabet)
        true_nllh = 4.4586751457870815

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0,phi=0,eta=0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_7 failed.")
    
    
    def test_8(self): 
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1]]])
        counts_a = countgen([(0,),(1,)],(0,))
        counts_b = countgen([(0,),(1,)],(0,))
        counts_c = countgen([(0,),(1,)],(0,))        
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b],'c':[counts_c]},alphabet)
        true_nllh = 5.0

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0,phi=0,eta=0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_8 failed.")
    
    
    def test_9(self): 
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1]]])
        counts_a = countgen([(0,),(1,)],(0,))
        counts_b = countgen([(0,),(1,)],(0,))
        counts_c = countgen([(0,),(1,)],('?',))
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b],'c':[counts_c]},alphabet)
        true_nllh = 6.513306124309698

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0,phi=0.1,eta=0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_9 failed.")
    
    def test_10(self): 
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1]]])
        counts_a = countgen([(0,),(1,)],(0,))
        counts_b = countgen([(0,),(1,)],('?',))
        counts_c = countgen([(0,),(1,)],(0,))
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b],'c':[counts_c]},alphabet)
        true_nllh = 6.513306124309698

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0,phi=0.1,eta=0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_10 failed.")
    
    
    def test_11(self): 
        true_nllh = 6.513306124309698
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1]]])
        counts_a = countgen([(0,),(1,)],('?',))
        counts_b = countgen([(0,),(1,)],(0,))
        counts_c = countgen([(0,),(1,)],(0,))
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b],'c':[counts_c]},alphabet)

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0,phi=0.1,eta=0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_11 failed.")
    
    def test_12(self): 
        true_nllh = 5.97198126969678
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1]]])
        counts_a = countgen([(0,),(1,)],(0,))
        counts_b = countgen([(0,),(1,)],(1,))
        counts_c = countgen([(0,),(1,)],('?',))
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b],'c':[counts_c]},alphabet)

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0,phi=0.1,eta=0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_12 failed.")
    
    def test_13(self): 
        true_nllh = 5.97198126969678
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1]]])
        counts_a = countgen([(0,),(1,)],(0,))
        counts_b = countgen([(0,),(1,)],('?',))
        counts_c = countgen([(0,),(1,)],(1,))
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b],'c':[counts_c]},alphabet)

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0,phi=0.1,eta=0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_13 failed.")
    
    def test_14(self): 
        true_nllh = 5.97198126969678
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1]]])
        counts_a = countgen([(0,),(1,)],('?',))
        counts_b = countgen([(0,),(1,)],(0,))
        counts_c = countgen([(0,),(1,)],(1,))
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b],'c':[counts_c]},alphabet)

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0,phi=0.1,eta=0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_14 failed.")
    
    def test_15(self): 
        true_nllh = 4.658719582178557
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1]]])
        counts_a = countgen([(0,),(1,)],(1,))
        counts_b = countgen([(0,),(1,)],('?',))
        counts_c = countgen([(0,),(1,)],(0,))
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b],'c':[counts_c]},alphabet)

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0,phi=0.1,eta=0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_15 failed.")
    
    def test_16(self): 
        true_nllh = 4.658719582178557
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1]]])
        counts_a = countgen([(0,),(1,)],('?',))
        counts_b = countgen([(0,),(1,)],(1,))
        counts_c = countgen([(0,),(1,)],(0,))
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b],'c':[counts_c]},alphabet)

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0,phi=0.1,eta=0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_16 failed.")
    
    def test_17(self): 
        true_nllh = 2.5980566021648364
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1]]])
        counts_a = countgen([(0,),(1,)],(1,))
        counts_b = countgen([(0,),(1,)],(1,))
        counts_c = countgen([(0,),(1,)],('?',))
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b],'c':[counts_c]},alphabet)

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0,phi=0.1,eta=0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_17 failed.")
    
    def test_18(self): 
        true_nllh = 2.695795750497349
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1]]])
        counts_a = countgen([(0,),(1,)],(1,))
        counts_b = countgen([(0,),(1,)],('?',))
        counts_c = countgen([(0,),(1,)],(1,))
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b],'c':[counts_c]},alphabet)

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0,phi=0.1,eta=0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_18 failed.")
    
    def test_19(self): 
        true_nllh = 2.695795750497349
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1]]])
        counts_a = countgen([(0,),(1,)],('?',))
        counts_b = countgen([(0,),(1,)],(1,))
        counts_c = countgen([(0,),(1,)],(1,))
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b],'c':[counts_c]},alphabet)

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0,phi=0.1,eta=0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_19 failed.")
    
    def test_20(self): 
        Q = [[{1:0.5,2:0.5}]]
        true_nllh = 1.0297894223949402
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1]]])
        counts_a = countgen([(0,),(1,)],(1,))
        counts_b = countgen([(0,),(1,)],(1,))
        counts_c = countgen([(0,),(1,)],(1,))
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b],'c':[counts_c]},alphabet)

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0,phi=0,eta=0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_20 failed.")

    def test_21(self): 
        Q = [[{1:1}]]
        true_nllh = 0.3215288449416738
        T = "(((a:1,b:1)ab:0,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1]]])
        counts_a = countgen([(0,),(1,)],(1,))
        counts_b = countgen([(0,),(1,)],(1,))
        counts_c = countgen([(0,),(1,)],(1,))
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b],'c':[counts_c]},alphabet)

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0,phi=0,eta=0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_21 failed.")
    
    def test_22(self): 
        # Testing when I introduce UMI count ties explicity.
        Q = [[{1:1}]]
        true_nllh = 0.14541345786885906 # with tie
        T = "((a:1,b:1)ab:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1]]])
        counts_a = countgen([(0,),(1,)],(1,))
        counts_a[(0,)] = counts_a[(1,)]
        counts_b = countgen([(0,),(1,)],(1,))
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b]},alphabet)

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0,phi=0,eta=0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_22 failed.")
    
    def test_23(self): 
        # Testing with silencing rate and no missing data.
        Q = [[{1:1}]]
        true_nllh = 0.3995946911551692
        T = "((a:1,b:1)ab:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]]) # only add -1 if we have silencing rate > 0 
        counts_a = countgen([(0,),(1,)],(1,),silencing=True)
        counts_b = countgen([(0,),(1,)],(1,),silencing=True)
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b]},alphabet)

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0.05,phi=0,eta=0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_23 failed.")

    def test_24(self):
        # Testing with silencing and missing data.
        Q = [[{1:1}]]
        true_nllh = 3.266041566926236
        T = "((a:1,b:1)ab:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]]) # only add -1 if we have silencing rate > 0
        counts_a = countgen([(0,),(1,)],(1,),silencing=True)
        counts_b = countgen([(0,),(1,)],('?',),silencing=True)
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b]},alphabet)

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0.05,phi=0,eta=0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_24 failed.")

    def test_25(self):
        # Testing with multistate.
        Q = [[{1:0.5,2:0.5}]]
        true_nllh = 1.0418276933439998
        T = "((a:1,b:1)ab:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,2]]])
        counts_a = countgen([(0,),(1,),(2,)],(1,))
        counts_b = countgen([(0,),(1,),(2,)],(1,))
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b]},alphabet)

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0,phi=0,eta=0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_25 failed.")

    def test_26(self):
        # Testing with phi and missing data.
        Q = [[{1:1}]]
        true_nllh = 3.1924390258104007
        T = "((a:1,b:1)ab:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1]]])
        counts_a = countgen([(0,),(1,)],(1,))
        counts_b = countgen([(0,),(1,)],('?',))
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b]},alphabet)

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0,phi=0.05,eta=0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_26 failed.")

    def test_27(self):
        # Testing with phi but no missing data.
        Q = [[{1:1}]]
        true_nllh = 0.35218127993027026
        T = "((a:1,b:1)ab:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1]]])
        counts_a = countgen([(0,),(1,)],(1,))
        counts_b = countgen([(0,),(1,)],(1,))
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b]},alphabet)

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0,phi=0.05,eta=0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_27 failed.")

    def test_28(self):
        # Testing with no missing data and no missing parameters.
        Q = [[{1:1}]]
        true_nllh = 0.24959469115516922
        # lower is better, no missing rate is better model fit
        T = "((a:1,b:1)ab:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1]]])
        counts_a = countgen([(0,),(1,)],(1,))
        counts_b = countgen([(0,),(1,)],(1,))
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b]},alphabet)

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0,phi=0,eta=0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_28 failed.")

    def test_29(self):
        # Testing with varying J and K
        Q = [[{1:1},{1:1}]]
        true_nllh = 3.249594691155169
        T = "((a:1,b:1)ab:1)r;"
        K = 1
        J = 2
        alphabet = Alphabet(K,J,[[[0,1],[0,1]]])
        counts_a = countgen([(0,0),(1,0),(0,1),(1,1)],(1,0))
        counts_b = countgen([(0,0),(1,0),(0,1),(1,1)],(1,0))
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b]},alphabet)

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0,phi=0,eta=0)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_29 failed.")

    def test_30(self):
        # Testing with confusion probability (ambient noise).
        Q = [[{1:1}]]
        true_nllh = 0.43543848121500794
        T = "((a:1,b:1)ab:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1]]])
        counts_a = countgen([(0,),(1,)],(1,))
        counts_b = countgen([(0,),(1,)],(1,))
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b]},alphabet)

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0,phi=0,eta=0.1)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_30 failed.")

    def test_31(self):
        # Test a tie with noise.
        Q = [[{1:1}]]
        true_nllh = 0.2335326144652978
        T = "((a:1,b:1)ab:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1]]])
        counts_a = countgen([(0,),(1,)],(1,))
        counts_a[(0,)] = counts_a[(1,)]
        counts_b = countgen([(0,),(1,)],(1,))
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b]},alphabet)

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0,phi=0,eta=0.1)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_31 failed.")

    def test_32(self):
        # Testing multistate and confusion rate.
        Q = [[{1:0.5,2:0.5}]]
        true_nllh = 1.2236554273339757
        T = "((a:1,b:1)ab:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,2]]])
        counts_a = countgen([(0,),(1,),(2,)],(1,))
        counts_b = countgen([(0,),(1,),(2,)],(1,))
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b]},alphabet)

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0,phi=0,eta=0.05)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_32 failed.")

    def test_33(self):
        # Test a tie bewteen two out of three cassettes.
        Q = [[{1:0.5,2:0.5}]]
        true_nllh = 1.052567463339014
        T = "((a:1,b:1)ab:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,2]]])
        counts_a = countgen([(0,),(1,),(2,)],(1,))
        counts_a[(0,)] = counts_a[(1,)]
        counts_b = countgen([(0,),(1,),(2,)],(1,))
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b]},alphabet)

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0,phi=0,eta=0.05)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_33 failed.")

    def test_34(self):
        # Test confusion probability with silenced states, but no missing data.
        Q = [[{1:1}]]
        true_nllh = 0.5720634603438768
        T = "((a:1,b:1)ab:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        counts_a = countgen([(0,),(1,)],(1,),silencing=True)
        counts_b = countgen([(0,),(1,)],(1,),silencing=True)
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b]},alphabet)

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0.05,phi=0,eta=0.1)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_34 failed.")

    def test_35(self):
        # Test confusion probability with silenced states, with missing data.
        Q = [[{1:1}]]
        true_nllh = 3.2178264980657496
        T = "((a:1,b:1)ab:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        counts_a = countgen([(0,),(1,)],(1,),silencing=True)
        counts_b = countgen([(0,),(1,)],('?',),silencing=True)
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b]},alphabet)

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0.05,phi=0,eta=0.1)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_35 failed.")

    def test_36(self):
        # Test confusion probability with silenced states and dropout.
        Q = [[{1:1}]]
        true_nllh = 2.6511966785803462
        T = "((a:1,b:1)ab:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        counts_a = countgen([(0,),(1,)],(1,),silencing=True)
        counts_b = countgen([(0,),(1,)],('?',),silencing=True)
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b]},alphabet)

        myModel = PMMN_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0.05,phi=0.05,eta=0.1)
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMTest in llh: test_36 failed.")

