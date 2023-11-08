from unit_tests.unit_tests_MLSolver import *
from unit_tests.unit_tests_Simulator import *
#from unit_tests.unit_tests_SpaLinSolver import *
from unit_tests.unit_tests_TopoSearch import *
from unit_tests.unit_tests_TopoSearchParallel import *
from unit_tests.unit_tests_EMSolver import *
from unit_tests.utils import *
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), 'helpers'))
if __name__ == '__main__':
    unittest.main()
