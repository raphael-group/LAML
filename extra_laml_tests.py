#! /usr/bin/env python
##from laml_unit_tests.unit_tests_MLSolver import *
##from laml_unit_tests.unit_tests_Simulator import *
#from laml_unit_tests.unit_tests_TopoSearch import *
#from laml_unit_tests.unit_tests_TopoSearchParallel import *
from laml_unit_tests.unit_tests_io import *
import sys
import os

if __name__ == '__main__':

    sys.path.append(os.path.dirname(__file__))
    sys.path.append(os.path.join(os.path.dirname(__file__), 'helpers'))
    print("Running extra tests for LAML...")
    unittest.main()

