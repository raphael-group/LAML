import os 
import unittest
from laml_libs.Count_model.PMMN_model import *
from treeswift import *
from laml_libs.IO_handler.sequence_lib import read_sequences
from random import random
from laml_libs.Count_model.utils import *
from laml_libs.Count_model.CharMtrx import CharMtrx
from .virtual_unit_tests import VirtualUnitTest

class PMMNTest_in_llh(VirtualUnitTest):
    # test in_llh computation 
    def test_1a(self): 
        charMtrx = {'a':[1],'b':[1],'c':[1]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        
        K = 1
        J = 1
        silence_mechanism = 'convolve'
        alphabet = Alphabet(K,J,[[[0,1,-1]]],silence_mechanism=silence_mechanism)
        DLT_data = CharMtrx(charMtrx,alphabet)
        priors = {'Q':Q,'silence_mechanism':silence_mechanism}

        true_nllh = 0.20665578828621584

        myModel = PMMN_model([T],{'DLT_data':DLT_data},priors,{'mu':1,'nu':0,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_1a failed.")
    
    def test_1b(self): 
        charMtrx = {'a':[1],'b':[1],'c':[1]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        
        K = 1
        J = 1
        silence_mechanism = 'separated'
        alphabet = Alphabet(K,J,[[[0,1]]],silence_mechanism=silence_mechanism)
        DLT_data = CharMtrx(charMtrx,alphabet)
        priors = {'Q':Q,'silence_mechanism':silence_mechanism}

        true_nllh = 0.20665578828621584

        myModel = PMMN_model([T],{'DLT_data':DLT_data},priors,{'mu':1,'nu':0,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_1b failed.")
    
    def test_2a(self): 
        charMtrx = {'a':[1],'b':[1],'c':[0]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        silence_mechanism = 'convolve'
        alphabet = Alphabet(K,J,[[[0,1,-1]]],silence_mechanism=silence_mechanism)
        DLT_data = CharMtrx(charMtrx,alphabet)
        priors = {'Q':Q,'silence_mechanism':silence_mechanism}

        true_nllh = 2.2495946917551692 

        myModel = PMMN_model([T],{'DLT_data':DLT_data},priors,{'mu':1,'nu':0,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_2a failed.")
    
    def test_2b(self): 
        charMtrx = {'a':[1],'b':[1],'c':[0]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        silence_mechanism = 'separated'
        alphabet = Alphabet(K,J,[[[0,1]]],silence_mechanism=silence_mechanism)
        DLT_data = CharMtrx(charMtrx,alphabet)
        priors = {'Q':Q,'silence_mechanism':silence_mechanism}

        true_nllh = 2.2495946917551692 

        myModel = PMMN_model([T],{'DLT_data':DLT_data},priors,{'mu':1,'nu':0,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_2b failed.")
   
    def test_3(self): 
        charMtrx = {'a':[1],'b':[0],'c':[1]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = CharMtrx(charMtrx,alphabet)
        
        true_nllh = 3.917350291274164 

        myModel = PMMN_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_3 failed.")
    
    def test_4(self): 
        # evaluate when cell A's argmax cassette state is not the true cassette state
        charMtrx = {'a':[0],'b':[1],'c':[1]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = CharMtrx(charMtrx,alphabet)
        
        true_nllh = 3.917350291274164 

        myModel = PMMN_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_4 failed.")
    
    def test_5(self): 
        # evaluate when cell B and C's argmax cassette state is not the true cassette state
        charMtrx = {'a':[1],'b':[0],'c':[0]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = CharMtrx(charMtrx,alphabet)
        
        true_nllh = 4.4586751457870815

        myModel = PMMN_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_5 failed.")
    
    def test_6(self): 
        charMtrx = {'a':[0],'b':[1],'c':[0]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = CharMtrx(charMtrx,alphabet)
        
        true_nllh = 4.4586751457870815

        myModel = PMMN_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_6 failed.")
    
    def test_7(self): 
        charMtrx = {'a':[0],'b':[0],'c':[1]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = CharMtrx(charMtrx,alphabet)
        
        true_nllh = 4.4586751457870815

        myModel = PMMN_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_7 failed.")
    
    
    def test_8(self): 
        charMtrx = {'a':[0],'b':[0],'c':[0]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = CharMtrx(charMtrx,alphabet)
        
        true_nllh = 5.0

        myModel = PMMN_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_8 failed.")
    
    
    def test_9a(self): 
        charMtrx = {'a':[0],'b':[0],'c':['?']}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        silence_mechanism = 'convolve'
        alphabet = Alphabet(K,J,[[[0,1,-1]]],silence_mechanism=silence_mechanism)
        DLT_data = CharMtrx(charMtrx,alphabet)
        priors = {'Q':Q,'silence_mechanism':silence_mechanism}
        
        true_nllh = 6.513306124309698

        myModel = PMMN_model([T],{'DLT_data':DLT_data},priors,{'mu':1,'nu':0,'phi':0.1,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_9 failed.")
    
    def test_9b(self): 
        charMtrx = {'a':[0],'b':[0],'c':['?']}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        silence_mechanism = 'separated'
        alphabet = Alphabet(K,J,[[[0,1]]],silence_mechanism=silence_mechanism)
        DLT_data = CharMtrx(charMtrx,alphabet)
        priors = {'Q':Q,'silence_mechanism':silence_mechanism}
        
        true_nllh = 6.513306124309698

        myModel = PMMN_model([T],{'DLT_data':DLT_data},priors,{'mu':1,'nu':0,'phi':0.1,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_9 failed.")
    
    def test_10a(self): 
        charMtrx = {'a':[0],'b':['?'],'c':[0]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        silence_mechanism = 'convolve'
        alphabet = Alphabet(K,J,[[[0,1,-1]]],silence_mechanism=silence_mechanism)
        DLT_data = CharMtrx(charMtrx,alphabet)
        priors = {'Q':Q,'silence_mechanism':silence_mechanism}
        
        true_nllh = 6.513306124309698

        myModel = PMMN_model([T],{'DLT_data':DLT_data},priors,{'mu':1,'nu':0,'phi':0.1,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_10 failed.")
    
    def test_10b(self): 
        charMtrx = {'a':[0],'b':['?'],'c':[0]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        silence_mechanism = 'separated'
        alphabet = Alphabet(K,J,[[[0,1]]],silence_mechanism=silence_mechanism)
        DLT_data = CharMtrx(charMtrx,alphabet)
        priors = {'Q':Q,'silence_mechanism':silence_mechanism}
        
        true_nllh = 6.513306124309698

        myModel = PMMN_model([T],{'DLT_data':DLT_data},priors,{'mu':1,'nu':0,'phi':0.1,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_10 failed.")
    
    def test_11(self): 
        charMtrx = {'a':['?'],'b':[0],'c':[0]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = CharMtrx(charMtrx,alphabet)
        
        true_nllh = 6.513306124309698

        myModel = PMMN_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0.1,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_11 failed.")
    
    def test_12(self): 
        charMtrx = {'a':[0],'b':[1],'c':['?']}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = CharMtrx(charMtrx,alphabet)
        
        true_nllh = 5.97198126969678

        myModel = PMMN_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0.1,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_12 failed.")
    
    def test_13(self): 
        charMtrx = {'a':[0],'b':['?'],'c':[1]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = CharMtrx(charMtrx,alphabet)
        
        true_nllh = 5.97198126969678

        myModel = PMMN_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0.1,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_13 failed.")
    
    def test_14(self): 
        charMtrx = {'a':['?'],'b':[0],'c':[1]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = CharMtrx(charMtrx,alphabet)
        
        true_nllh = 5.97198126969678

        myModel = PMMN_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0.1,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_14 failed.")
    
    def test_15(self): 
        charMtrx = {'a':[1],'b':['?'],'c':[0]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = CharMtrx(charMtrx,alphabet)
        
        true_nllh = 4.658719582178557

        myModel = PMMN_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0.1,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_15 failed.")
    
    def test_16(self): 
        charMtrx = {'a':['?'],'b':[1],'c':[0]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = CharMtrx(charMtrx,alphabet)
        
        true_nllh = 4.658719582178557

        myModel = PMMN_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0.1,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_16 failed.")
    
    def test_17(self): 
        charMtrx = {'a':[1],'b':[1],'c':['?']}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = CharMtrx(charMtrx,alphabet)
        
        true_nllh = 2.5980566021648364

        myModel = PMMN_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0.1,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_17 failed.")
    
    def test_18(self): 
        charMtrx = {'a':[1],'b':['?'],'c':[1]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = CharMtrx(charMtrx,alphabet)
        
        true_nllh = 2.695795750497349

        myModel = PMMN_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0.1,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_18 failed.")
    
    def test_19(self): 
        charMtrx = {'a':['?'],'b':[1],'c':[1]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = CharMtrx(charMtrx,alphabet)
        
        true_nllh = 2.695795750497349

        myModel = PMMN_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0.1,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_19 failed.")
    
    def test_20(self): 
        charMtrx = {'a':[1],'b':[1],'c':[1]}
        Q = [[{1:0.5,2:0.5}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1,2]]])
        DLT_data = CharMtrx(charMtrx,alphabet)
        
        true_nllh = 1.0297894223949402

        myModel = PMMN_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_20 failed.")
    
    def test_21(self): 
        charMtrx = {'a':[1],'b':[1],'c':[1]}
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:0,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = CharMtrx(charMtrx,alphabet)
        
        true_nllh = 0.3215288449416738

        myModel = PMMN_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_21 failed.")
    
    def test_22(self): 
        # Testing with varying J and K
        charMtrx = {'a':[(1,0)],'b':[(1,0)]}
        Q = [[{1:1},{1:1}]]
        T = "((a:1,b:1)ab:1)r;"
        K = 1
        J = 2
        alphabet = Alphabet(K,J,[[[0,1,-1],[0,1,-1]]])
        DLT_data = CharMtrx(charMtrx,alphabet)
        
        true_nllh = 3.249594691155169

        myModel = PMMN_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':11,'nu':0,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_22 failed.")
    
    def test_23a(self): 
        # Testing with silencing rate and no missing data.
        charMtrx = {'a':[(1,)],'b':[(1,)]}
        Q = [[{1:1}]]
        T = "((a:1,b:1)ab:1)r;"
        K = 1
        J = 1
        silence_mechanism = 'convolve'
        alphabet = Alphabet(K,J,[[[0,1,-1]]],silence_mechanism=silence_mechanism)
        DLT_data = CharMtrx(charMtrx,alphabet)
        priors = {'Q':Q,'silence_mechanism':silence_mechanism}
        
        true_nllh = 0.3995946911551692

        myModel = PMMN_model([T],{'DLT_data':DLT_data},priors,{'mu':11,'nu':0.05,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_23a failed.")
    
    def test_23b(self): 
        # Testing with silencing rate and no missing data.
        charMtrx = {'a':[(1,)],'b':[(1,)]}
        Q = [[{1:1}]]
        T = "((a:1,b:1)ab:1)r;"
        K = 1
        J = 1
        silence_mechanism = 'separated'
        alphabet = Alphabet(K,J,[[[0,1]]],silence_mechanism=silence_mechanism)
        DLT_data = CharMtrx(charMtrx,alphabet)
        priors = {'Q':Q,'silence_mechanism':silence_mechanism}
        
        true_nllh = 0.3995946911551692

        myModel = PMMN_model([T],{'DLT_data':DLT_data},priors,{'mu':11,'nu':0.05,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_23b failed.")
    
    def test_24a(self): 
        # Testing with silencing and missing data.
        charMtrx = {'a':[(1,)],'b':[('?',)]}
        Q = [[{1:1}]]
        T = "((a:1,b:1)ab:1)r;"
        K = 1
        J = 1
        silence_mechanism = 'convolve'
        alphabet = Alphabet(K,J,[[[0,1,-1]]],silence_mechanism=silence_mechanism)
        DLT_data = CharMtrx(charMtrx,alphabet)
        priors = {'Q':Q,'silence_mechanism':silence_mechanism}
        
        true_nllh = 3.266041566926236

        myModel = PMMN_model([T],{'DLT_data':DLT_data},priors,{'mu':11,'nu':0.05,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_24a failed.")
    
    def test_24b(self): 
        # Testing with silencing and missing data.
        charMtrx = {'a':[(1,)],'b':[('?',)]}
        Q = [[{1:1}]]
        T = "((a:1,b:1)ab:1)r;"
        K = 1
        J = 1
        silence_mechanism = 'separated'
        alphabet = Alphabet(K,J,[[[0,1]]],silence_mechanism=silence_mechanism)
        DLT_data = CharMtrx(charMtrx,alphabet)
        priors = {'Q':Q,'silence_mechanism':silence_mechanism}
        
        true_nllh = 3.266041566926236

        myModel = PMMN_model([T],{'DLT_data':DLT_data},priors,{'mu':11,'nu':0.05,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_24b failed.")
    
    def test_25a(self): 
        # Testing with multistate.
        charMtrx = {'a':[(1,)],'b':[(1,)]}
        Q = [[{1:0.5,2:0.5}]]
        T = "((a:1,b:1)ab:1)r;"
        K = 1
        J = 1
        silence_mechanism = 'convolve'
        alphabet = Alphabet(K,J,[[[0,1,2,-1]]],silence_mechanism=silence_mechanism)
        DLT_data = CharMtrx(charMtrx,alphabet)
        priors = {'Q':Q,'silence_mechanism':silence_mechanism}

        true_nllh = 1.0418276933439998

        myModel = PMMN_model([T],{'DLT_data':DLT_data},priors,{'mu':11,'nu':0,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_25a failed.")
    
    def test_25b(self): 
        # Testing with multistate.
        charMtrx = {'a':[(1,)],'b':[(1,)]}
        Q = [[{1:0.5,2:0.5}]]
        T = "((a:1,b:1)ab:1)r;"
        K = 1
        J = 1
        silence_mechanism = 'separated'
        alphabet = Alphabet(K,J,[[[0,1,2]]],silence_mechanism=silence_mechanism)
        DLT_data = CharMtrx(charMtrx,alphabet)
        priors = {'Q':Q,'silence_mechanism':silence_mechanism}

        true_nllh = 1.0418276933439998

        myModel = PMMN_model([T],{'DLT_data':DLT_data},priors,{'mu':11,'nu':0,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_25b failed.")
    
    def test_26(self): 
        # Testing with phi and missing data.
        charMtrx = {'a':[(1,)],'b':[('?',)]}
        Q = [[{1:1}]]
        T = "((a:1,b:1)ab:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = CharMtrx(charMtrx,alphabet)
        
        true_nllh = 3.1924390258104007

        myModel = PMMN_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':11,'nu':0,'phi':0.05,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_26 failed.")
    
    def test_27(self): 
        # Testing with phi and no missing data.
        charMtrx = {'a':[(1,)],'b':[(1,)]}
        Q = [[{1:1}]]
        T = "((a:1,b:1)ab:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = CharMtrx(charMtrx,alphabet)
        
        true_nllh = 0.35218127993027026

        myModel = PMMN_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':11,'nu':0,'phi':0.05,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_27 failed.")
    
    def test_28(self): 
        # Testing with no missing data and no missing parameters.
        charMtrx = {'a':[(1,)],'b':[(1,)]}
        Q = [[{1:1}]]
        T = "((a:1,b:1)ab:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = CharMtrx(charMtrx,alphabet)
        
        true_nllh = 0.24959469115516922

        myModel = PMMN_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':11,'nu':0,'phi':0,'rho':1})
        my_nllh = myModel.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="PMMNTest in llh: test_28 failed.")
