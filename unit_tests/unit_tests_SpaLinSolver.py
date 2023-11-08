import unittest
from scmail_libs.SpaLin_solver import SpaLin_solver

class SpaLinTest(unittest.TestCase):
    def test_1(self): 
        Q = [{1:1}]
        T = "((a:1,b:1):1,c:1):1;"
        msa = {'a':[1],'b':[1],'c':[1]}
        true_llh = -0.20665578828621584
        locations = {'a':(0,0),'b':(0,0),'c':(0,0),'d':(0,0)}
        sigma = 0

        data = {'charMtrx':msa,'locations':locations}
        prior = {'Q':Q}
        params = {'phi':0,'nu':0,'sigma':sigma}

        mySolver = SpaLin_solver([T],data,prior,params)
        mySolver.az_partition()
        my_llh = -mySolver.negative_llh()
        self.assertAlmostEqual(true_llh,my_llh,places=5,msg="SpaLinTest: test_1 failed.")
