import unittest
from laml_libs.IO_handler.sequence_lib import read_sequences
from laml_libs.Count_model.PMM_base import PMM_model, Alphabet,AlleleTable
from treeswift import *
from math import log
from random import random
from laml_libs import DEFAULT_STRATEGY
from copy import deepcopy
import pkg_resources
from .utils import *
from .virtual_unit_tests import VirtualUnitTest

class PMM_Test_out_llh(VirtualUnitTest):
    def test_1(self):
        T = "(((a:1,b:1)e:1,(c:1,d:1)f:1)g:1)r;"
        Q = [[{1:1.0}]]
        charMtrx = {'a':[1],'b':[1],'c':[1],'d':[1]}
        phi = 0
        nu = 0.5
        self.check_outllh(T,Q,charMtrx,mu=1,phi=phi,nu=nu,test_no=1)
   
    def test_2(self):
        T = "(((a:1,b:1)e:1,(c:1,d:1)f:1)g:1)r;"
        Q = [[{1:1.0}]]
        charMtrx = {'a':[1],'b':[1],'c':[1],'d':[1]}
        nu = 0.1
        phi = 0
        self.check_outllh(T,Q,charMtrx,mu=1,phi=phi,nu=nu,test_no=2)
    
    def test_3(self):
        T = "(((a:1,b:1)e:1,(c:1,d:1)f:1)g:1)r;"
        Q = [[{1:1.0}]]
        charMtrx = {'a':[0],'b':[0],'c':[0],'d':[0]}
        nu = 0.25
        phi = 0
        self.check_outllh(T,Q,charMtrx,mu=1,phi=phi,nu=nu,test_no=1)
    
    def test_4(self):
        T = "(((a:0.1,b:1)e:1,(c:0.1,d:1)f:1)g:1)r;"
        Q = [[{1:1.0}]]
        charMtrx = {'a':[0],'b':[1],'c':[0],'d':[1]}
        nu = 0.15
        phi = 0
        self.check_outllh(T,Q,charMtrx,mu=1,phi=phi,nu=nu,test_no=4)
    
    def test_5(self):
        T = "(((a:0.1,b:1)e:1,(c:0.1,d:1)f:1)g:1)r;"
        Q = [[{1:1.0}]]
        charMtrx = {'a':[1],'b':[0],'c':[1],'d':[0]}
        nu = 0.1
        phi = 0
        self.check_outllh(T,Q,charMtrx,mu=1,phi=phi,nu=nu,test_no=5)
        
    def test_6(self):
        T = "(((a:0.1,b:1)e:1,(c:0.1,d:1)f:0.7)g:0.3)r;"
        Q = [[{1:1.0}]]
        charMtrx = {'a':[1],'b':[0],'c':[1],'d':[0]}
        nu = 0.19
        phi = 0
        self.check_outllh(T,Q,charMtrx,mu=1,phi=phi,nu=nu,test_no=6)
    
    def test_7(self):
        T = "(((a:0.1,b:1)e:1,(c:0.1,d:1)f:0.7)g:0.3)r;"
        Q = [[{1:1.0}]]
        charMtrx = {'a':['?'],'b':[0],'c':[1],'d':[0]}
        nu = 0.39
        phi = 0
        self.check_outllh(T,Q,charMtrx,mu=1,phi=phi,nu=nu,test_no=7)
    
    def test_8(self):
        T = "(((a:0.1,b:1)e:1,(c:0.1,d:1)f:0.7)g:0.3)r;"
        Q = [[{1:1.0}]]
        charMtrx = {'a':['?'],'b':[0],'c':[1],'d':['?']}
        nu = 0.39
        phi = 0
        self.check_outllh(T,Q,charMtrx,mu=1,phi=phi,nu=nu,test_no=8)
    
    def test_9(self):
        T = "(((a:0.1,b:1)e:1,(c:0.1,d:1)f:0.7)g:0.3)r;"
        Q = [[{1:0.5,2:0.5}]]
        charMtrx = {'a':['?'],'b':[2],'c':[1],'d':['?']}
        nu = 0.3
        phi = 0
        self.check_outllh(T,Q,charMtrx,mu=1,phi=phi,nu=nu,test_no=9)
    
    def test_10(self):
        T = "(((a:0.47,b:1.3)e:1.1,(c:0.14,d:1.1)f:0.72)g:0.39)r;"
        Q = [[{1:0.5,2:0.3,3:0.2}]]
        charMtrx = {'a':['?'],'b':[2],'c':[1],'d':['?']}
        nu = 0.22
        phi = 0.3
        self.check_outllh(T,Q,charMtrx,mu=1,phi=phi,nu=nu,test_no=10)
    
    def test_11(self):
        T = "((((a:0.47,b:1.3)e:1.1,c:0.14)f:0.8,d:1.1)g:0.2)r;"
        Q = [[{1:0.5,2:0.3,3:0.2}]]
        charMtrx = {'a':[1],'b':[1],'c':[1],'d':[1]}
        nu = 0.22
        phi = 0.01
        self.check_outllh(T,Q,charMtrx,mu=1,phi=phi,nu=nu,test_no=11)

    def test_12(self):
        T = "((((a:0.47,b:1.3)e:1.1,c:0.14)f:0.8,d:1.1)g:0.2)r;"
        Q = [[{1:0.5,2:0.3,3:0.2}]]
        charMtrx = {'a':['?'],'b':['?'],'c':[1],'d':[0]}
        nu = 0.22
        phi = 0.2
        self.check_outllh(T,Q,charMtrx,mu=1,phi=phi,nu=nu,test_no=12)
    
    def test_13(self):
        T = "((((a:0.47,b:1.3)e:1.1,c:0.14)f:0.8,d:1.1)g:0.2)r;"
        Q = [[{1:0.5,2:0.3,3:0.2}]]
        charMtrx = {'a':['?'],'b':['?'],'c':['?'],'d':[0]}
        nu = 0.22
        phi = 0.5
        self.check_outllh(T,Q,charMtrx,mu=1,phi=phi,nu=nu,test_no=13)
    
    def test_14(self):
        T = "((((a:0.47,b:1.3)e:1.1,c:0.14)f:0.8,d:1.1)g:0.2)r;"
        Q = [[{1:0.5,2:0.3,3:0.2}]]
        charMtrx = {'a':['?'],'b':['?'],'c':['?'],'d':[1]}
        nu = 0.22
        phi = 0.01
        self.check_outllh(T,Q,charMtrx,mu=1,phi=phi,nu=nu,test_no=14)
    
    def test_15(self):
        T = "((((a:0.47,b:1.3)e:1.1,c:0.14)f:0.8,d:1.1)g:0.2)r;"
        Q = [[{1:0.5,2:0.3,3:0.2}]]
        charMtrx = {'a':['?'],'b':['?'],'c':['?'],'d':['?']}
        nu = 0.22
        phi = 0.9
        self.check_outllh(T,Q,charMtrx,mu=1,phi=phi,nu=nu,test_no=15)
    
    def test_16(self):
        T = "((((a:0.47,b:1.3)e:1.1,c:0.14)f:0.8,d:1.1)g:0.2)r;"
        Q = [[{1:0.5,2:0.3,3:0.2}]]
        charMtrx = {'a':[1],'b':[2],'c':[3],'d':['?']}
        nu = 0.22
        phi = 0.03
        self.check_outllh(T,Q,charMtrx,mu=1,phi=phi,nu=nu,test_no=16)
    
    def test_17(self):
        T = "((((a:0.47,b:1.3)f:1.1,c:0.14)g:0.8,(d:1.1,e:0.2)h:0.2)g:0.01)r;"
        Q = [[{1:0.5,2:0.3,3:0.2}]]
        charMtrx = {'a':[1],'b':[2],'c':[3],'d':['?'],'e':['?']}
        nu = 0.22
        phi = 0.5
        self.check_outllh(T,Q,charMtrx,mu=1,phi=phi,nu=nu,test_no=17)
    
    def test_18(self):
        T = "((((a:0.47,b:1.3)f:1.1,c:0.14)g:0.8,(d:1.1,e:0.2)h:0.2)g:0.01)r;"
        Q = [[{1:0.5,2:0.3,3:0.2}]]
        charMtrx = {'a':['?'],'b':['?'],'c':['?'],'d':['?'],'e':['?']}
        nu = 0.26
        phi = 0.7
        self.check_outllh(T,Q,charMtrx,mu=1,phi=phi,nu=nu,test_no=18)
    
    def test_19(self):
        T = "((((a:0.47,b:1.3)f:1.1,c:0.14)g:0.8,(d:1.1,e:0.2)h:0.2)g:0.01)r;"
        Q = [[{1:0.2,2:0.6,3:0.2}]]
        charMtrx = {'a':['?'],'b':['?'],'c':['?'],'d':[2],'e':[2]}
        nu = 0.26
        phi = 0.1
        self.check_outllh(T,Q,charMtrx,mu=1,phi=phi,nu=nu,test_no=19)
    
    def test_20(self):
        T = "((((a:0.47,b:1.3)f:1.1,c:0.14)g:0.8,(d:1.1,e:0.2)h:0.2)g:0.01)r;"
        Q = [[{1:0.2,2:0.6,3:0.2}]]
        charMtrx = {'a':['?'],'b':['?'],'c':[2],'d':[2],'e':[2]}
        nu = 0.52
        phi = 0.8
        self.check_outllh(T,Q,charMtrx,mu=1,phi=phi,nu=nu,test_no=20)
    
    def test_21(self):
        T = "((((a:0.47,b:1.3)f:1.1,c:0.14)g:0.8,(d:1.1,e:0.2)h:0.2)g:0.01)r;"
        Q = [[{1:0.2,2:0.6,3:0.2}]]
        charMtrx = {'a':['?'],'b':['?'],'c':[2],'d':[1],'e':[2]}
        nu = 0.52
        phi = 0.001
        self.check_outllh(T,Q,charMtrx,mu=1,phi=phi,nu=nu,test_no=21)
    
    def test_22(self):
        T = "(((((a:0.47,b:1.3)f:1.1,c:0.14)g:0.8,d:1.1)h:0.2,e:0.2)g:0.2)r;"
        Q = [[{1:0.2,2:0.6,3:0.2}]]
        charMtrx = {'a':['?'],'b':['?'],'c':[2],'d':['?'],'e':[2]}
        nu = 0.2
        phi = 0.82
        self.check_outllh(T,Q,charMtrx,mu=1,phi=phi,nu=nu,test_no=22)
    
    def test_23(self):
        T = "(((((a:0.47,b:1.3)f:1.1,c:0.14)g:0.8,d:1.1)h:0.2,e:0.2)g:0.2)r;"
        Q = [[{1:0.2,2:0.6,3:0.2}]]
        charMtrx = {'a':['?'],'b':['?'],'c':[2],'d':['?'],'e':[0]}
        nu = 0.2
        phi = 0.1
        self.check_outllh(T,Q,charMtrx,mu=1,phi=phi,nu=nu,test_no=23)
    
    def test_24(self):
        T = "(((((a:0.47,b:1.3)f:1.1,c:0.14)g:0.8,d:1.1)h:0.2,e:0.2)g:0.2)r;"
        Q = [[{1:0.2,2:0.6,3:0.2}]]
        charMtrx = {'a':['?'],'b':[0],'c':[2],'d':['?'],'e':[0]}
        nu = 0.2
        phi = 0.3
        self.check_outllh(T,Q,charMtrx,mu=1,phi=phi,nu=nu,test_no=24)
    
    def test_25(self):
        T = "(((((a:0.47,b:1.3)f:1.1,c:0.14)g:0.8,d:1.1)h:0.2,e:0.2)g:0.2)r;"
        Q = [[{1:0.2,2:0.6,3:0.2}]]
        charMtrx = {'a':['?'],'b':[0],'c':[2],'d':['?'],'e':[1]}
        nu = 0.2
        phi = 0.3
        self.check_outllh(T,Q,charMtrx,mu=1,phi=phi,nu=nu,test_no=25)
    
    def test_26(self):
        treedata_path = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_Count_model/test_PMM_base/test1_n25.tre')
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
            Q[i].append(q)
        self.check_outllh(T,Q,charMtrx,mu=1,phi=phi,nu=nu,test_no=26,give_label=True)
