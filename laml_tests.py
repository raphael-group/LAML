#! /usr/bin/env python
from laml_unit_tests.unit_tests_MLSolver import *
from laml_unit_tests.unit_tests_mosek import *
from laml_unit_tests.unit_tests_TopoSearch import *
from laml_unit_tests.unit_tests_EMSolver import *
from laml_unit_tests.unit_tests_fastEMSolver import *
from laml_unit_tests.unit_tests_TopoSearch_wfastEMSolver import *
from laml_unit_tests.utils import *

##from unit_tests.unit_tests_Simulator import *
##from unit_tests.unit_tests_SpaLinSolver import *
##from laml_unit_tests.unit_tests_TopoSearchParallel import *

import sys
import os

def main():
    sys.path.append(os.path.dirname(__file__))
    sys.path.append(os.path.join(os.path.dirname(__file__), 'helpers'))
    print("Running tests for LAML...")

    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    # Import your test modules
    import laml_unit_tests.unit_tests_MLSolver as ml
    import laml_unit_tests.unit_tests_mosek as mosek
    import laml_unit_tests.unit_tests_TopoSearch as topo
    import laml_unit_tests.unit_tests_EMSolver as em
    import laml_unit_tests.unit_tests_fastEMSolver as fastem
    import laml_unit_tests.unit_tests_TopoSearch_wfastEMSolver as topo_fastem

    # Add all tests from these modules
    for module in [ml, mosek, topo, em, fastem, topo_fastem]:
        suite.addTests(loader.loadTestsFromModule(module))

    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)


if __name__ == '__main__':
    #print(os.getcwd())
    main()



