import os 
import unittest
from laml_libs.Count_model.PMMC_model import *
from treeswift import *
from laml_libs.IO_handler.sequence_lib import read_sequences
from random import random
from laml_libs.Count_model.utils import *
from laml_libs.Count_model.AlleleTable import AlleleTable
from .virtual_unit_tests import VirtualUnitTest

class PMMCTest_in_llh(VirtualUnitTest):
    # test in_llh computation 
    def test_1(self): 
        alleleTable = {'a':[{1:1,0:0}],'b':[{1:1,0:0}],'c':[{1:1,0:0}]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = AlleleTable(alleleTable,alphabet)

        true_nllh = 0.20665578828621584

        myModel = PMMC_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMCTest in llh: test_1 failed.")
    
    def test_2(self): 
        alleleTable = {'a':[{1:1,0:0}],'b':[{1:1,0:0}],'c':[{1:0,0:1}]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = AlleleTable(alleleTable,alphabet)

        true_nllh = 2.2495946917551692 

        myModel = PMMC_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMCTest in llh: test_2 failed.")
   
    def test_3(self): 
        alleleTable = {'a':[{1:1,0:0}],'b':[{1:0,0:1}],'c':[{1:1,0:0}]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = AlleleTable(alleleTable,alphabet)
        
        true_nllh = 3.917350291274164 

        myModel = PMMC_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMCTest in llh: test_3 failed.")
    
    def test_4(self): 
        # evaluate when cell A's argmax cassette state is not the true cassette state
        alleleTable = {'a':[{1:0,0:1}],'b':[{1:1,0:0}],'c':[{1:1,0:0}]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = AlleleTable(alleleTable,alphabet)
        
        true_nllh = 3.917350291274164 

        myModel = PMMC_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMCTest in llh: test_4 failed.")
    
    def test_5(self): 
        # evaluate when cell B and C's argmax cassette state is not the true cassette state
        alleleTable = {'a':[{1:1,0:0}],'b':[{1:0,0:1}],'c':[{1:0,0:1}]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = AlleleTable(alleleTable,alphabet)
        
        true_nllh = 4.4586751457870815

        myModel = PMMC_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMCTest in llh: test_5 failed.")
    
    def test_6(self): 
        alleleTable = {'a':[{1:0,0:1}],'b':[{1:1,0:0}],'c':[{1:0,0:1}]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = AlleleTable(alleleTable,alphabet)
        
        true_nllh = 4.4586751457870815

        myModel = PMMC_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMCTest in llh: test_6 failed.")
    
    def test_7(self): 
        alleleTable = {'a':[{1:0,0:1}],'b':[{1:0,0:1}],'c':[{1:1,0:0}]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = AlleleTable(alleleTable,alphabet)
        
        true_nllh = 4.4586751457870815

        myModel = PMMC_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMCTest in llh: test_7 failed.")
    
    
    def test_8(self): 
        alleleTable = {'a':[{1:0,0:1}],'b':[{1:0,0:1}],'c':[{1:0,0:1}]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = AlleleTable(alleleTable,alphabet)
        
        true_nllh = 5.0

        myModel = PMMC_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMCTest in llh: test_8 failed.")
    
    
    def test_9(self): 
        alleleTable = {'a':[{1:0,0:1}],'b':[{1:0,0:1}],'c':[{}]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = AlleleTable(alleleTable,alphabet)
        
        true_nllh = 6.513306124309698

        myModel = PMMC_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0.1,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMCTest in llh: test_9 failed.")
    
    def test_10(self): 
        alleleTable = {'a':[{1:0,0:1}],'b':[{}],'c':[{1:0,0:1}]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = AlleleTable(alleleTable,alphabet)
        
        true_nllh = 6.513306124309698

        myModel = PMMC_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0.1,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMCTest in llh: test_10 failed.")
    
    def test_11(self): 
        alleleTable = {'a':[{}],'b':[{1:0,0:1}],'c':[{1:0,0:1}]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = AlleleTable(alleleTable,alphabet)
        
        true_nllh = 6.513306124309698

        myModel = PMMC_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0.1,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMCTest in llh: test_11 failed.")
    
    def test_12(self): 
        alleleTable = {'a':[{1:0,0:1}],'b':[{1:1,0:0}],'c':[{}]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = AlleleTable(alleleTable,alphabet)
        
        true_nllh = 5.97198126969678

        myModel = PMMC_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0.1,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMCTest in llh: test_12 failed.")
    
    def test_13(self): 
        alleleTable = {'a':[{1:0,0:1}],'b':[{}],'c':[{1:1,0:0}]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = AlleleTable(alleleTable,alphabet)
        
        true_nllh = 5.97198126969678

        myModel = PMMC_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0.1,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMCTest in llh: test_13 failed.")
    
    def test_14(self): 
        alleleTable = {'a':[{}],'b':[{1:0,0:1}],'c':[{1:1,0:0}]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = AlleleTable(alleleTable,alphabet)
        
        true_nllh = 5.97198126969678

        myModel = PMMC_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0.1,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMCTest in llh: test_14 failed.")
    
    def test_15(self): 
        alleleTable = {'a':[{1:1,0:0}],'b':[{}],'c':[{1:0,0:1}]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = AlleleTable(alleleTable,alphabet)
        
        true_nllh = 4.658719582178557

        myModel = PMMC_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0.1,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMCTest in llh: test_15 failed.")
    
    def test_16(self): 
        alleleTable = {'a':[{}],'b':[{1:1,0:0}],'c':[{1:0,0:1}]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = AlleleTable(alleleTable,alphabet)
        
        true_nllh = 4.658719582178557

        myModel = PMMC_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0.1,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMCTest in llh: test_16 failed.")
    
    def test_17(self): 
        alleleTable = {'a':[{1:1,0:0}],'b':[{1:1,0:0}],'c':[{}]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = AlleleTable(alleleTable,alphabet)
        
        true_nllh = 2.5980566021648364

        myModel = PMMC_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0.1,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMCTest in llh: test_17 failed.")
    
    def test_18(self): 
        alleleTable = {'a':[{1:1,0:0}],'b':[{}],'c':[{1:1,0:0}]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = AlleleTable(alleleTable,alphabet)
        
        true_nllh = 2.695795750497349

        myModel = PMMC_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0.1,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMCTest in llh: test_18 failed.")
    
    def test_19(self): 
        alleleTable = {'a':[{}],'b':[{1:1,0:0}],'c':[{1:1,0:0}]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = AlleleTable(alleleTable,alphabet)
        
        true_nllh = 2.695795750497349

        myModel = PMMC_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0.1,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMCTest in llh: test_19 failed.")
    
    def test_20(self): 
        alleleTable = {'a':[{1:1,0:0}],'b':[{1:1,0:0}],'c':[{1:1,0:0}]}
        Q = [[{1:0.5,2:0.5}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1,2]]])
        DLT_data = AlleleTable(alleleTable,alphabet)
        
        true_nllh = 1.0297894223949402

        myModel = PMMC_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMCTest in llh: test_20 failed.")
    
    def test_21(self): 
        alleleTable = {'a':[{1:1,0:0}],'b':[{1:1,0:0}],'c':[{1:1,0:0}]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:0,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = AlleleTable(alleleTable,alphabet)
        
        true_nllh = 0.3215288449416738

        myModel = PMMC_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMCTest in llh: test_21 failed.")
    
    def test_22(self): 
        # Testing with varying J and K
        alleleTable = {'a':[{(1,0):1}],'b':[{(1,0):1}]}
        Q = [[{1:1},{1:1}]]
        T = "((a:1,b:1)ab:1)r;"
        K = 1
        J = 2
        alphabet = Alphabet(K,J,[[[0,1,-1],[0,1,-1]]])
        DLT_data = AlleleTable(alleleTable,alphabet)
        
        true_nllh = 3.249594691155169

        myModel = PMMC_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':11,'nu':0,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMCTest in llh: test_22 failed.")
    
    def test_23(self): 
        # Testing with silencing rate and no missing data.
        alleleTable = {'a':[{(1,):1,(0,):0}],'b':[{(1,):1,(0,):0}]}
        Q = [[{1:1}]]
        T = "((a:1,b:1)ab:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = AlleleTable(alleleTable,alphabet)
        
        true_nllh = 0.3995946911551692

        myModel = PMMC_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':11,'nu':0.05,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMCTest in llh: test_23 failed.")
    
    def test_24(self): 
        # Testing with silencing and missing data.
        #alleleTable = {'a':[{1:1,0:0}],'b':[{}]}
        alleleTable = {'a':[{1:1}],'b':[{}]}
        Q = [[{1:1}]]
        T = "((a:1,b:1)ab:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = AlleleTable(alleleTable,alphabet)

        true_nllh = 3.266041566926236

        myModel = PMMC_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':11,'nu':0.05,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMCTest in llh: test_24 failed.")
    
    def test_25(self): 
        # Testing with multistate.
        alleleTable = {'a':[{1:1}],'b':[{1:1}]}
        Q = [[{1:0.5,2:0.5}]]
        T = "((a:1,b:1)ab:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,2,-1]]])
        DLT_data = AlleleTable(alleleTable,alphabet)
        
        true_nllh = 1.0418276933439998

        myModel = PMMC_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':11,'nu':0,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMCTest in llh: test_25 failed.")
    
    def test_26(self): 
        # Testing with phi and missing data.
        alleleTable = {'a':[{1:1}],'b':[{}]}
        Q = [[{1:1}]]
        T = "((a:1,b:1)ab:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = AlleleTable(alleleTable,alphabet)
        
        true_nllh = 3.1924390258104007

        myModel = PMMC_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':11,'nu':0,'phi':0.05,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMCTest in llh: test_26 failed.")
    
    def test_27(self): 
        # Testing with phi and no missing data.
        alleleTable = {'a':[{1:1}],'b':[{1:1}]}
        Q = [[{1:1}]]
        T = "((a:1,b:1)ab:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = AlleleTable(alleleTable,alphabet)
        
        true_nllh = 0.35218127993027026

        myModel = PMMC_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':11,'nu':0,'phi':0.05,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMCTest in llh: test_27 failed.")
    
    def test_28(self): 
        # Testing with no missing data and no missing parameters.
        alleleTable = {'a':[{1:1}],'b':[{1:1}]}
        Q = [[{1:1}]]
        T = "((a:1,b:1)ab:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = AlleleTable(alleleTable,alphabet)
        
        true_nllh = 0.24959469115516922

        myModel = PMMC_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':11,'nu':0,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMCTest in llh: test_28 failed.")
