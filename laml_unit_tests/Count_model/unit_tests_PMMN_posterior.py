import os 
import pkg_resources
import unittest
from laml_libs.Count_model.PMMN_model import *
from treeswift import *
from laml_libs.IO_handler.sequence_lib import read_sequences
from random import random
from laml_libs.Count_model.utils import *
from .virtual_unit_tests import VirtualUnitTest
from math import *
from laml_libs.Count_model.CharMtrx import CharMtrx

class PMMNTest_posterior(VirtualUnitTest):
    # test posterior computation
    def check_posterior(self,charMtrx,Q,T,phi,nu,rho,test_no): ##### only works if K = 1, J = 1, and alphabet = {0,1,-1} #####
        ref_file = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_Count_model/test_PMMN/test_posterior/test_'+str(test_no)+'.txt')
        value_types = ['post_node_z','post_node_s','post_edge_zz','post_edge_za','post_edge_zs','post_edge_aa','post_edge_as']
        ref_values = {x:{} for x in value_types}
        with open(ref_file,'r') as fin:
            for line in fin:
                tokens = line.strip().split()
                node_label = tokens[0]
                value_type = tokens[1]
                value = float(tokens[2])
                if value_type in ref_values:
                    ref_values[value_type][node_label] = value

        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = CharMtrx(charMtrx,alphabet)
        
        myModel = PMMN_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':nu,'phi':phi,'rho':rho})
        myModel.Estep()

        #print(ref_values)
        for node in myModel.trees[0].traverse_postorder():
            #print(node.label)
            # test post_node_z
            if node.label in ref_values['post_node_z']:
                true = ref_values['post_node_z'][node.label]
                est = node.log_node_posterior[0][(0,)]
                self.assertAlmostEqual(true,est,places=5,msg="PMMNTest posterior: test_" + str(test_no) +  " failed for post_node_z.")
            # test post_node_s
            if node.label in ref_values['post_node_s']:
                true = ref_values['post_node_s'][node.label]
                #print(node.label,node.out_llh[0][(0,)],node.out_llh[0][(-1,)])
                est = node.log_node_posterior[0][(-1,)]
                self.assertAlmostEqual(true,est,places=5,msg="PMMNTest posterior: test_" + str(test_no) +  " failed for post_node_s.")
            # test post_edge_zz
            if node.label in ref_values['post_edge_zz']:
                true = ref_values['post_edge_zz'][node.label]
                est = node.log_edge_posterior[0][((0,),(0,))]
                self.assertAlmostEqual(true,est,places=5,msg="PMMNTest posterior: test_" + str(test_no) + " failed for post_edge_zz.")
            # test post_edge_za
            if node.label in ref_values['post_edge_za']:
                true = ref_values['post_edge_za'][node.label]
                est = node.log_edge_posterior[0][((0,),(1,))]
                self.assertAlmostEqual(true,est,places=5,msg="PMMNTest posterior: test_" + str(test_no) + " failed for post_edge_za.")
            # test post_edge_zs
            if node.label in ref_values['post_edge_zs']:
                true = ref_values['post_edge_zs'][node.label]
                est = node.log_edge_posterior[0][((0,),(-1,))]
                self.assertAlmostEqual(true,est,places=5,msg="PMMNTest posterior: test_" + str(test_no) + " failed for post_edge_zs.")
            # test post_edge_aa
            if node.label in ref_values['post_edge_aa']:
                true = ref_values['post_edge_aa'][node.label]
                est = node.log_edge_posterior[0][((1,),(1,))]
                self.assertAlmostEqual(true,est,places=5,msg="PMMNTest posterior: test_" + str(test_no) + " failed for post_edge_aa.")
            # test post_edge_as
            if node.label in ref_values['post_edge_as']:
                true = ref_values['post_edge_as'][node.label]
                est = node.log_edge_posterior[0][((1,),(-1,))]
                self.assertAlmostEqual(true,est,places=5,msg="PMMNTest posterior: test_" + str(test_no) + " failed for post_edge_as.")

    def test_1(self): 
        Q = [[{1:1}]]
        charMtrx = {'a':[1],'b':[1],'c':[1]}
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        phi = 0
        nu = 0
        rho = 1
        test_no = 1
        self.check_posterior(charMtrx,Q,T,phi,nu,rho,test_no)
    
    def test_2(self): 
        Q = [[{1:1}]]
        charMtrx = {'a':['?'],'b':['?'],'c':[1],'d':[0]}
        T = "(((a:0.47,b:1.3)e:1.1,c:0.14)f:0.8,d:1.1)r:0.2;"
        phi = 0.2
        nu = 0.22
        rho = 1
        test_no = 2
        self.check_posterior(charMtrx,Q,T,phi,nu,rho,test_no)
