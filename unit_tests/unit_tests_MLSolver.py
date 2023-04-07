import os 
import unittest
from problin_libs.ML_solver import ML_solver

class MLTest(unittest.TestCase):
    def test_1(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[1],'b':[1],'c':[1]}
        true_nllh = 0.20665578828621584

        mySolver = ML_solver(T,{'charMtrx':msa},{'Q':Q},{'phi':1e-10,'nu':1e-10})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_1 failed.")
    
    def test_2(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[1],'b':[1],'c':[0]}
        true_nllh = 2.2495946917551692 

        mySolver = ML_solver(T,{'charMtrx':msa},{'Q':Q},{'phi':1e-10,'nu':1e-10})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_2 failed.")
    
    def test_3(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[1],'b':[0],'c':[1]}
        true_nllh = 3.917350291274164 

        #mySolver = ML_solver(msa,Q,T,phi=1e-10,nu=1e-10)
        mySolver = ML_solver(T,{'charMtrx':msa},{'Q':Q},{'phi':1e-10,'nu':1e-10})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_3 failed.")
    
    def test_4(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[0],'b':[1],'c':[1]}
        true_nllh = 3.917350291274164 

        #mySolver = ML_solver(msa,Q,T,phi=1e-10,nu=1e-10)
        mySolver = ML_solver(T,{'charMtrx':msa},{'Q':Q},{'phi':1e-10,'nu':1e-10})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_4 failed.")
    
    def test_5(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[1],'b':[0],'c':[0]}
        true_nllh = 4.4586751457870815

        #mySolver = ML_solver(msa,Q,T,phi=1e-10,nu=1e-10)
        mySolver = ML_solver(T,{'charMtrx':msa},{'Q':Q},{'phi':1e-10,'nu':1e-10})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_5 failed.")
    
    def test_6(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[0],'b':[1],'c':[0]}
        true_nllh = 4.4586751457870815

        #mySolver = ML_solver(msa,Q,T,phi=1e-10,nu=1e-10)
        mySolver = ML_solver(T,{'charMtrx':msa},{'Q':Q},{'phi':1e-10,'nu':1e-10})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_6 failed.")
    
    def test_7(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[0],'b':[0],'c':[1]}
        true_nllh = 4.4586751457870815

        #mySolver = ML_solver(msa,Q,T,phi=1e-10,nu=1e-10)
        mySolver = ML_solver(T,{'charMtrx':msa},{'Q':Q},{'phi':1e-10,'nu':1e-10})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_7 failed.")
    
    def test_8(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[0],'b':[0],'c':[0]}
        true_nllh = 5.0

        #mySolver = ML_solver(msa,Q,T,phi=1e-10,nu=1e-10)
        mySolver = ML_solver(T,{'charMtrx':msa},{'Q':Q},{'phi':1e-10,'nu':1e-10})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_8 failed.")
    
    def test_9(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[0],'b':[0],'c':['?']}
        true_nllh = 6.513306124309698

        #mySolver = ML_solver(msa,Q,T,phi=0.1,nu=1e-10)
        mySolver = ML_solver(T,{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':1e-10})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_9 failed.")
    
    def test_10(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[0],'b':['?'],'c':[0]}
        true_nllh = 6.513306124309698

        #mySolver = ML_solver(msa,Q,T,phi=0.1,nu=1e-10)
        mySolver = ML_solver(T,{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':1e-10})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_10 failed.")
    
    def test_11(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':['?'],'b':[0],'c':[0]}
        true_nllh = 6.513306124309698

        #mySolver = ML_solver(msa,Q,T,phi=0.1,nu=1e-10)
        mySolver = ML_solver(T,{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':1e-10})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_11 failed.")
    
    def test_12(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[0],'b':[1],'c':['?']}
        true_nllh = 5.97198126969678

        #mySolver = ML_solver(msa,Q,T,phi=0.1,nu=1e-10)
        mySolver = ML_solver(T,{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':1e-10})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_12 failed.")
    
    def test_13(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[0],'b':['?'],'c':[1]}
        true_nllh = 5.97198126969678

        #mySolver = ML_solver(msa,Q,T,phi=0.1,nu=1e-10)
        mySolver = ML_solver(T,{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':1e-10})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_13 failed.")
    
    def test_14(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':['?'],'b':[0],'c':[1]}
        true_nllh = 5.97198126969678

        #mySolver = ML_solver(msa,Q,T,phi=0.1,nu=1e-10)
        mySolver = ML_solver(T,{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':1e-10})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_14 failed.")
    
    def test_15(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[1],'b':['?'],'c':[0]}
        true_nllh = 4.658719582178557

        #mySolver = ML_solver(msa,Q,T,phi=0.1,nu=1e-10)
        mySolver = ML_solver(T,{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':1e-10})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_15 failed.")
    
    def test_16(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':['?'],'b':[1],'c':[0]}
        true_nllh = 4.658719582178557

        #mySolver = ML_solver(msa,Q,T,phi=0.1,nu=1e-10)
        mySolver = ML_solver(T,{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':1e-10})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_16 failed.")
    
    def test_17(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[1],'b':[1],'c':['?']}
        true_nllh = 2.5980566021648364

        #mySolver = ML_solver(msa,Q,T,phi=0.1,nu=1e-10)
        mySolver = ML_solver(T,{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':1e-10})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_17 failed.")
    
    def test_18(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[1],'b':['?'],'c':[1]}
        true_nllh = 2.695795750497349

        #mySolver = ML_solver(msa,Q,T,phi=0.1,nu=1e-10)
        mySolver = ML_solver(T,{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':1e-10})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_18 failed.")
    
    def test_19(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':['?'],'b':[1],'c':[1]}
        true_nllh = 2.695795750497349

        #mySolver = ML_solver(msa,Q,T,phi=0.1,nu=1e-10)
        mySolver = ML_solver(T,{'charMtrx':msa},{'Q':Q},{'phi':0.1,'nu':1e-10})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_19 failed.")
    
    def test_20(self): 
        Q = [{1:0.5,2:0.5}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[1],'b':[1],'c':[1]}
        true_nllh = 1.0297894223949402

        #mySolver = ML_solver(msa,Q,T,phi=1e-10,nu=1e-10)
        mySolver = ML_solver(T,{'charMtrx':msa},{'Q':Q},{'phi':1e-10,'nu':1e-10})
        #mySolver.az_partition(mySolver.params)
        mySolver.az_partition()
        my_nllh = mySolver.negative_llh()
        self.assertAlmostEqual(true_nllh,my_nllh,places=5,msg="MLTest: test_20 failed.")

    def test_21(self):
        T = "((a:1,b:1)e:2,(c:1,d:1)f:2)g:1;"
        Q = [{0:0, 1:0.5, 2:0.5}, {0:0, 1:0.5, 2:0.5}]
        msa = {'a':[0, 0], 'b':[1, 1], 'c':[1, 2], 'd':[1, 2]}
        correct_branches = set(['b', 'e'])
        
        #mySolver = ML_solver(msa,Q,T)
        mySolver = ML_solver(T,{'charMtrx':msa},{'Q':Q})#,{'phi':1e-10,'nu':1e-10})
        branches = mySolver.score_branches()
        branch = max(branches, key=lambda item:item[1])[0].label

        self.assertIn(branch, correct_branches, msg="MLTest: test_21 failed.")
    
    ''' 
    def test_22(self):
        T = "[&R] ((a:0.5,b:1)e:2,(c:1,d:0.5)f:1)g:1;"
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[1, 1, 1, 1, 1], 'b':[0, 0, 0, 0, 0], 'c':[0, 0, 0, 0, 0], 'd':[1, 1, 1, 1, 1]}
        correct_topology = "[&R] (((a:0.5,d:0.5)f:1,c:1)e:2,b:1)g:1;"
        
        mySolver = ML_solver(msa,Q,T)
        mySolver.topology_search(maxiter=10, verbose=False, strategy="vanilla", trynextbranch=False, outdir="./unit_tests", conv=0.2)
        opt_topo = mySolver.params.tree.newick()
        self.assertEqual(opt_topo, correct_topology, msg="MLTest: test_22 failed.")
        os.remove("./unit_tests/results_nni_topo_llh.txt")
        os.remove("./unit_tests/results_nni_topo_progress.nwk") '''
