#! /usr/bin/env python
from laml_unit_tests.PMM_original.unit_tests_MLSolver import *
from laml_unit_tests.TopoSearch.unit_tests_TopoSearch import *
from laml_unit_tests.PMM_original.unit_tests_EMSolver import *
from laml_unit_tests.ThirdParty_packages.unit_tests_mosek import *
import sys
import os

if __name__ == '__main__':
    sys.path.append(os.path.dirname(__file__))
    sys.path.append(os.path.join(os.path.dirname(__file__), 'helpers'))
    print("Running tests for LAML...")
    unittest.main()

