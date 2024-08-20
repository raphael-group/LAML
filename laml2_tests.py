#! /usr/bin/env python
from laml_unit_tests.Count_model.unit_tests_PMM_base_in_llh import *
from laml_unit_tests.Count_model.unit_tests_PMM_base_out_llh import *
from laml_unit_tests.Count_model.unit_tests_PMM_base_posterior import *
from laml_unit_tests.Count_model.unit_tests_PMM_base_scipy_opt import *
from laml_unit_tests.Count_model.unit_tests_PMM_base_EM_opt import *
from laml_unit_tests.Count_model.unit_tests_PMM_with_noise_in_llh import *
from laml_unit_tests.Count_model.unit_tests_PMM_with_noise_out_llh import *
import sys
import os

if __name__ == '__main__':
    sys.path.append(os.path.dirname(__file__))
    sys.path.append(os.path.join(os.path.dirname(__file__), 'helpers'))
    print("Running tests for LAML2...")
    unittest.main()

