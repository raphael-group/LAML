import os 
import unittest
from laml_libs.PMM_original.ML_solver import ML_solver
from treeswift import *
from laml_libs.IO_handler.sequence_lib import read_sequences

class MLTest(unittest.TestCase):
    def test_1(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[1],'b':[1],'c':[1]}
        true_nllh = 0.20665578828621584

        mySolver = ML_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0,'nu':0})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_1 failed.")
    
    def test_2(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[1],'b':[1],'c':[0]}
        true_nllh = 2.2495946917551692 

        mySolver = ML_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0,'nu':0})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_2 failed.")
    
    def test_3(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[1],'b':[0],'c':[1]}
        true_nllh = 3.917350291274164 

        #mySolver = ML_solver(msa,Q,T,phi=0,nu=0)
        mySolver = ML_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0,'nu':0})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_3 failed.")
    
    def test_4(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[0],'b':[1],'c':[1]}
        true_nllh = 3.917350291274164 

        #mySolver = ML_solver(msa,Q,T,phi=0,nu=0)
        mySolver = ML_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0,'nu':0})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_4 failed.")
    
    def test_5(self): 
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
        T = "((a:1,b:1)e:1,(c:1,d:1)f:1)r:1;"
        Q = [{1:1.0}]
        msa = {'a':[1],'b':[1],'c':[1],'d':[1]}
        phi = 0
        nu = 0.5
        
        true_nllh = 3.656149065448911

        mySolver = ML_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':phi,'nu':nu})
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_21 failed.")    
    
    def test_22(self):
        T = "((a:1,b:1)e:1,(c:1,d:1)f:1)r:1;"
        Q = [{1:1.0}]
        msa = {'a':[1],'b':[1],'c':[1],'d':[1]}
        nu = 0.1
        phi = 0
        
        true_nllh = 0.8561490654489112

        mySolver = ML_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':phi,'nu':nu})
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_22 failed.")    

    def test_23(self):
        T = "((a:1,b:1)e:1,(c:1,d:1)f:1)r:1;"
        Q = [{1:1.0}]
        msa = {'a':[0],'b':[0],'c':[0],'d':[0]}
        nu = 0.25
        phi = 0
    
        true_nllh = 8.75 

        mySolver = ML_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':phi,'nu':nu})
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_23 failed.")    

    def test_24(self):
        T = "((a:0.1,b:1)e:1,(c:0.1,d:1)f:1)r:1;"
        Q = [{1:1.0}]
        msa = {'a':[0],'b':[1],'c':[0],'d':[1]}
        nu = 0.15
        phi = 0
        
        true_nllh = 4.897350290774163

        mySolver = ML_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':phi,'nu':nu})
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_24 failed.")    
    
    def test_25(self):
        T = "((a:0.1,b:1)e:1,(c:0.1,d:1)f:1)r:1;"
        Q = [{1:1.0}]
        msa = {'a':[1],'b':[0],'c':[1],'d':[0]}
        nu = 0.1
        phi = 0

        true_nllh = 10.22433692208818

        mySolver = ML_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':phi,'nu':nu})
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_25 failed.")    
    
    def test_26(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):0,c:1):1;"
        msa = {'a':[1],'b':[1],'c':[1]}
        true_nllh = 0.3215288449416738

        mySolver = ML_solver([T],{'charMtrx':msa},{'Q':Q},{'phi':0,'nu':0})
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_26 failed.")   
    
    # test optimization
    def test_27(self):
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[1, 1, 0, 0, 0], 'b':[1, 1, 1, 0, 0], 'c':[0, 0, 0, 1, 0], 'd':[0, 0, 0, 1, 0]}
        data = {'charMtrx':msa}
        prior = {'Q':Q}
         
        tree_list = ['(((b,c),d),a);', '((b,c),(d,a));', '(d,((b,c),a));', '(d,(b,(c,a)));', '(d,(c,(b,a)));', '(((b,d),c),a);', '((b,d),(c,a));', '(c,((b,d),a));', '(c,(b,(d,a)));', '(c,(d,(b,a)));', '(((c,d),b),a);', '((c,d),(b,a));', '(b,((c,d),a));', '(b,(c,(d,a)));', '(b,(d,(c,a)));']  
        true_nllh_list = [11.809140931727208, 11.809141253319298, 11.80914097392261, 11.809140974410134, 10.338626804278578, 11.809141098008908, 11.809140926336006, 11.809141047672332, 11.809141363494154, 10.33862520405755, 9.322029697571756, 7.851513459377366, 9.32202949252424, 11.809140996738018, 11.809140939639644]
        
        for T,true_nllh in zip(tree_list,true_nllh_list):
            mySolver = ML_solver([T],data,prior)
            test_nllh,_ = mySolver.optimize(initials=1,verbose=-1,ultra_constr=False)
            self.assertAlmostEqual(true_nllh,test_nllh,places=4,msg="MLTest: test_27 failed.")
    
    def test_28(self):
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[0, 1, 1, 1, 1], 'b':[1, 0, 0, 0, 0], 'c':[1, 0, 0, 0, 0], 'd':[0, 1, 1, 1, 1]}
        data = {'charMtrx':msa}
        prior = {'Q':Q}
       
        tree_list = ['(((b,c),d),a);', '((b,c),(d,a));', '(d,((b,c),a));', '(d,(b,(c,a)));', '(d,(c,(b,a)));', '(((b,d),c),a);', '((b,d),(c,a));', '(c,((b,d),a));', '(c,(b,(d,a)));', '(c,(d,(b,a)));', '(((c,d),b),a);', '((c,d),(b,a));', '(b,((c,d),a));', '(b,(c,(d,a)));', '(b,(d,(c,a)));']
        true_nllh_list = [7.595936888280069, 5.078900745505063, 7.595936937640891, 10.0830486451717, 10.083048625619265, 10.08305235692635, 10.083048544279704, 10.083048478217442, 7.566011711889361, 10.083048503399493, 10.083048490449752, 10.083048551149702, 10.083048664742668, 7.566011487530715, 10.083048480708316]

        for T,true_nllh in zip(tree_list,true_nllh_list):
            mySolver = ML_solver([T],data,prior)
            test_nllh,_ = mySolver.optimize(initials=1,verbose=-1,ultra_constr=False)
            self.assertAlmostEqual(true_nllh,test_nllh,places=4,msg="MLTest: test_28 failed.")
