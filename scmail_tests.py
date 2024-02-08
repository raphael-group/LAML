#! /usr/bin/env python
from scmail_unit_tests.unit_tests_MLSolver import *
#from unit_tests.unit_tests_Simulator import *
#from unit_tests.unit_tests_SpaLinSolver import *
from scmail_unit_tests.unit_tests_TopoSearch import *
#from scmail_unit_tests.unit_tests_TopoSearchParallel import *
from scmail_unit_tests.unit_tests_EMSolver import *
from scmail_unit_tests.utils import *
import sys
import os

#def main():


if __name__ == '__main__':
    #print(os.getcwd())
    sys.path.append(os.path.dirname(__file__))
    sys.path.append(os.path.join(os.path.dirname(__file__), 'helpers'))
    print("Running tests for sc-MAIL...")
    unittest.main()
