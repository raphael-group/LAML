#! /usr/bin/env python
from laml_unit_tests.unit_tests_MLSolver import *
from laml_unit_tests.unit_tests_mosek import *
from laml_unit_tests.unit_tests_TopoSearch import *
from laml_unit_tests.unit_tests_EMSolver import *
from laml_unit_tests.unit_tests_fastEMSolver import *
from laml_unit_tests.unit_tests_TopoSearch_wfastEMSolver import *

#from laml_unit_tests.utils import *

##from unit_tests.unit_tests_Simulator import *
##from unit_tests.unit_tests_SpaLinSolver import *

##from laml_unit_tests.unit_tests_TopoSearchParallel import *

import sys
import os

#def main():
    



if __name__ == '__main__':
    #print(os.getcwd())
    sys.path.append(os.path.dirname(__file__))
    sys.path.append(os.path.join(os.path.dirname(__file__), 'helpers'))
    print("Running tests for LAML...")
    unittest.main()

