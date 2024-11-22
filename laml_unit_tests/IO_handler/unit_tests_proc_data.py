import sys
import os 
import unittest
import cvxpy as cp
import pkg_resources
import subprocess
import re 
from random import randint 
from laml_libs.IO_handler.DLT_parser import *
import shutil
from math import *
import ast

# Functional tests
class ProcData(unittest.TestCase):

    def test_1(self):
        ## Test unified command character matrix as input, process into json.

        charmtrx_file = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_proc_data/small_character_matrix.csv')
        delimiter = ","
        missing_state = "?"
        outputfile = "laml_unit_tests/test_data/test_proc_data/small_character_matrix.json"

        if os.path.exists(outputfile):
            os.remove(outputfile)

        tmp = DLT_parser(datafile=charmtrx_file, delimiter=delimiter, missing_state=missing_state, outputfile=outputfile)
        
        # check that the outputfile was generated
        self.assertTrue(os.path.exists(outputfile), f"File {outputfile} was not created")
        self.assertEqual(tmp.K, 6)
        self.assertEqual(tmp.J, 1)
        self.assertEqual(len(tmp.data), 4)

        alphabet_ds = tmp.set_alphabet()
        self.assertEqual(tmp.alphabet.M, [4, 2, 3, 3, 2, 3])
        self.assertEqual(tmp.alphabet.get_site_alphabet(0,0), {0, 1, 2, 3, -1})

        # check that it added a 0 in
        self.assertEqual(tmp.alphabet.get_site_alphabet(4,0), {0, 1, -1})

        #teardown
        if os.path.exists(outputfile):
            os.remove(outputfile)
    
    def test_2(self):
        ## Test default outputfile path, with unified command character matrix as input, process into json.

        charmtrx_file = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_proc_data/small_character_matrix.csv')
        delimiter = ","
        missing_state = "?"

        outputfile = "laml_unit_tests/test_data/test_proc_data/small_character_matrix.json"
        if os.path.exists(outputfile):
            os.remove(outputfile)

        tmp = DLT_parser(datafile=charmtrx_file, delimiter=delimiter, missing_state=missing_state)
       
        # check that the default outputfile was generated
        self.assertTrue(os.path.exists(outputfile), f"File {outputfile} was not created")
        self.assertEqual(tmp.K, 6)
        self.assertEqual(tmp.J, 1)
        self.assertEqual(len(tmp.data), 4)

        alphabet_ds = tmp.set_alphabet()
        self.assertEqual(tmp.alphabet.M, [4, 2, 3, 3, 2, 3])
        self.assertEqual(tmp.alphabet.get_site_alphabet(0,0), {0, 1, 2, 3, -1})

        # check that it added a 0 in
        self.assertEqual(tmp.alphabet.get_site_alphabet(4,0), {0, 1, -1})

        #teardown
        if os.path.exists(outputfile):
            os.remove(outputfile)

    def test_3(self):
        ## Test default outputfile path, with unified command json character matrix as input, process into json.

        charmtrx_file = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_proc_data/example_small_character_matrix.json')
        tmp = DLT_parser(datafile=charmtrx_file)
       
        self.assertTrue(tmp.outputfile is None)
        self.assertEqual(tmp.K, 6)
        self.assertEqual(tmp.J, 1)
        self.assertEqual(len(tmp.data), 4)

        # check that the alphabet is read correctly
        alphabet_ds = tmp.set_alphabet()
        self.assertEqual(tmp.alphabet.M, [4, 2, 3, 3, 2, 3])
        self.assertEqual(tmp.alphabet.get_site_alphabet(0,0), {0, 1, 2, 3, -1})

        # check that it added a 0 in
        self.assertEqual(tmp.alphabet.get_site_alphabet(4,0), {0, 1, -1})

        
    def test_4(self):
        ## Test single command character matrix as input, process into json.

        charmtrx_file = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_proc_data/character_matrix.csv')
        delimiter = ","
        missing_state = "?"
        outputfile = "laml_unit_tests/test_data/test_proc_data/example1_test.json"

        if os.path.exists(outputfile):
            os.remove(outputfile)

        tmp = DLT_parser(datafile=charmtrx_file, delimiter=delimiter, missing_state=missing_state, outputfile=outputfile)
        
        # check that the outputfile was generated
        self.assertTrue(os.path.exists(outputfile), f"File {outputfile} was not created")
        self.assertEqual(tmp.K, 30)
        self.assertEqual(tmp.J, 1)
        self.assertEqual(len(tmp.data), 250)

        #teardown
        if os.path.exists(outputfile):
            os.remove(outputfile)

    def test_5(self):
        ## Test character matrix as input with individual commands, process into json.

        charmtrx_file = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_proc_data/character_matrix.csv')
        delimiter = ","
        missing_state = "?"
        outputfile = "laml_unit_tests/test_data/test_proc_data/example1_test.json"

        tmp = DLT_parser(delimiter=delimiter, missing_state=missing_state)
        
        if os.path.exists(outputfile):
            os.remove(outputfile)

        tmp.process_datafile(datafile=charmtrx_file, delimiter=delimiter, missing_state=missing_state, outputfile=outputfile, max_allele_per_cassette=None)
        
        # check that the outputfile was generated
        self.assertTrue(os.path.exists(outputfile), f"File {outputfile} was not created")
        self.assertEqual(tmp.K, 30)
        self.assertEqual(tmp.J, 1)
        self.assertEqual(len(tmp.data), 250)
        
        #teardown
        if os.path.exists(outputfile):
            os.remove(outputfile)

    def test_6(self):
        ## Test small json allele table as input.
        alleletable_file= pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_proc_data/example_small_alleletable.json')

        tmp = DLT_parser(datafile=alleletable_file)
       
        self.assertTrue(tmp.outputfile is None)
        self.assertEqual(tmp.K, 2)
        self.assertEqual(tmp.J, 3)
        self.assertEqual(tmp.num_cells, 4)

        # check that the alphabet is read correctly
        alphabet_ds = tmp.set_alphabet()
        self.assertEqual(tmp.alphabet.M, [prod([4, 2, 3]), prod([3, 2, 2])])
        self.assertEqual(tmp.alphabet.get_site_alphabet(0,0), {-1, 0, 1, 2, 3})
        self.assertEqual(tmp.alphabet.get_site_alphabet(1,1), {-1, 0, 1})  
        
    def test_7(self):
        ## Test big json allele table as input.

        alleleTableFile = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_proc_data/TLS_Bar8.json')
        delimiter = ","
        missing_state = "?"

        tmp = DLT_parser(datafile=alleleTableFile, delimiter=delimiter, missing_state=missing_state)

        self.assertEqual(tmp.K, 7)
        self.assertEqual(tmp.J, 3)
    
    def test_8(self):
        ## Test parse_prior for character matrix in csv format
        charmtrx_file = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_proc_data/small_character_matrix.csv')
        priorfile = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_proc_data/small_priors.csv')
        delimiter = ","
        missing_state = "?"
        outputfile = "laml_unit_tests/test_data/test_proc_data/small_character_matrix.json"

        if os.path.exists(outputfile):
            os.remove(outputfile)

        tmp = DLT_parser(datafile=charmtrx_file, delimiter=delimiter, missing_state=missing_state, outputfile=outputfile)
        Q = tmp.parse_prior(priorfile,tmp.J)

        # I introduced a new state in the prior file at site index 3 that isn't observed in the json, should be in the alphabet
        self.assertEqual(tmp.alphabet.M, [4, 2, 3, 4, 2, 3]) # alphabet is made according to priorfile
        self.assertEqual(tmp.alphabet.get_site_alphabet(0,0), {0, 1, 2, 3, -1})
        self.assertEqual(tmp.alphabet.get_site_alphabet(4,0), {0, 1, -1})

        # check the priors
        self.assertEqual(Q[0], [{1: 0.33, 2: 0.33, 3: 0.33, 0: 0}])
        self.assertEqual(Q[3], [{1: 0.5, 2: 0.4, 3: 0.1, 0: 0}])

        # teardown
        if os.path.exists(outputfile):
            os.remove(outputfile)
    
    def test_9(self):
        ## Test parse_prior for character matrix in json format
        charmtrx_file = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_proc_data/example_small_character_matrix.json')
        priorfile = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_proc_data/small_priors.csv')
        tmp = DLT_parser(datafile=charmtrx_file)
        Q = tmp.parse_prior(priorfile,tmp.J)
        
        # I introduced a new state in the prior file at site index 3 that isn't observed in the json, should be in the alphabet
        self.assertEqual(tmp.alphabet.M, [4, 2, 3, 4, 2, 3]) # alphabet is made according to priorfile
        self.assertEqual(tmp.alphabet.get_site_alphabet(0,0), {0, 1, 2, 3, -1})
        self.assertEqual(tmp.alphabet.get_site_alphabet(4,0), {0, 1, -1})
        
        # check the priors
        self.assertEqual(Q[0], [{1: 0.33, 2: 0.33, 3: 0.33, 0: 0}])
        self.assertEqual(Q[0][0], {1: 0.33, 2: 0.33, 3: 0.33, 0: 0})
        self.assertEqual(Q[3][0], {1: 0.5, 2: 0.4, 3: 0.1, 0: 0})

    
    def test_10(self):
        ## Test parse_prior for allele table in json format
        alleletable_file = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_proc_data/example_small_alleletable.json')
        priorfile = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_proc_data/small_priors.csv')
        tmp = DLT_parser(datafile=alleletable_file)

        self.assertEqual(tmp.J, 3)
        Q = tmp.parse_prior(priorfile,tmp.J)

        # I introduced a new state in the prior file at site index 3 that isn't observed in the json, should be in the alphabet
        self.assertEqual(tmp.alphabet.M, [24, 24]) # [ 4 * 2 * 3 ], ...
        self.assertEqual(tmp.alphabet.get_site_alphabet(0,0), {-1, 0, 1, 2, 3})
        self.assertEqual(tmp.alphabet.get_site_alphabet(1,0), {-1, 0, 1, 2, 3})
        self.assertEqual(tmp.alphabet.get_site_alphabet(1,1), {-1, 0, 1})

        # check the priors, split into their cassettes
        self.assertEqual(Q[0][0], {1: 0.33, 2: 0.33, 3: 0.33, 0: 0})
        self.assertEqual(Q[1][0], {1: 0.5, 2: 0.4, 3: 0.1, 0: 0}) # site idx 3
    
    def test_11(self):
        ## Test get_from_path for allele table in json format with uniform priors
        alleletable_file = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_proc_data/example_small_alleletable.json')
        tmp = DLT_parser()
        DLT_data, Q = tmp.get_from_path(datafile=alleletable_file)

        ### Uniform priors should drop the -1
        
        # I introduced a new state in the prior file at site index 3 that isn't observed in the json, should be in the alphabet
        self.assertEqual(tmp.alphabet.M, [24, 12]) # the extra entry at site index 3 is not there in the uniform priors
        self.assertEqual(tmp.alphabet.get_site_alphabet(0,0), {-1, 0, 1, 2, 3})
        self.assertEqual(tmp.alphabet.get_site_alphabet(1,1), {-1, 0, 1})

        # check the priors, split into their cassettes
        self.assertEqual(Q[0][1], {1: 1.0, 0: 0})
        self.assertEqual(Q[0][2], {1: 0.5, 2: 0.5, 0: 0})
        self.assertEqual(Q[1][0], {1: 0.5, 2: 0.5, 0: 0})
        
    def test_12(self):
        ## Test get_from_path for allele table in json format with prior file
        alleletable_file = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_proc_data/example_small_alleletable.json')
        priorfile = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_proc_data/small_priors.csv')
        tmp = DLT_parser()
        DLT_data, Q = tmp.get_from_path(datafile=alleletable_file, priorfile=priorfile)

        # I introduced a new state in the prior file at site index 3 that isn't observed in the json, should be in the alphabet
        self.assertEqual(tmp.alphabet.M, [24, 24])
        self.assertEqual(tmp.alphabet.get_site_alphabet(1,0), {-1, 0, 1, 2, 3})
        self.assertEqual(tmp.alphabet.get_site_alphabet(0,0), {-1, 0, 1, 2, 3})
        self.assertEqual(tmp.alphabet.get_site_alphabet(1,1), {-1, 0, 1})

    def test_13(self):
        ## Test max_allele_per_cassette. 
        alleletable_file = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_proc_data/example_small_alleletable.json')
        priorfile = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_proc_data/small_priors.csv')
        tmp = DLT_parser()
        self.assertEqual(tmp.max_allele_per_cassette, None)
        tmp = DLT_parser(max_allele_per_cassette=2)
        self.assertEqual(tmp.max_allele_per_cassette, 2)

        DLT_data, Q = tmp.get_from_path(datafile=alleletable_file, priorfile=priorfile, max_allele_per_cassette=1)

        # I introduced a new state in the prior file at site index 3 that isn't observed in the json, should be in the alphabet
        self.assertEqual(tmp.alphabet.M, [24, 24])
        self.assertEqual(tmp.max_allele_per_cassette, 1)
        self.assertEqual(tmp.alphabet.get_site_alphabet(0,0), {-1, 0, 1, 2, 3})
        self.assertEqual(tmp.alphabet.get_site_alphabet(1,0), {-1, 0, 1, 2, 3})
        self.assertEqual(tmp.alphabet.get_site_alphabet(1,1), {-1, 0, 1})
    
    def test_14(self):
        ## Test example1 data structs.
        charMtrx_file = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_proc_data/example1/character_matrix.csv')
        priorfile = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_proc_data/priors.csv')
        tmp = DLT_parser(datafile=charMtrx_file, priorfile=priorfile)

        # compare to structs in old_readin_structs
        
        true_alphabet_data_struct_file = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_proc_data/example1/old_readin_structs/alphabet_data_struct.csv')
        true_charMtrx_data_struct_file = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_proc_data/example1/old_readin_structs/charMtrx_data_struct.csv')
        true_priors_file = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_proc_data/example1/old_readin_structs/priors.json')
        true_cassette_alphabet_file = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_proc_data/example1/old_readin_structs/cassette_alphabet.csv')

        with open(true_alphabet_data_struct_file, "r") as f:
            true_alphabet_data_struct = json.load(f)
        with open(true_charMtrx_data_struct_file, "r") as f:
            true_charMtrx_data_struct = json.load(f)
        with open(true_priors_file, "r") as f:
            true_priors = json.load(f)

        # 1. tmp.DLT_data.data_struct and true_charMtrx_data_struct
        for cell_name in true_charMtrx_data_struct:
            self.assertIn(cell_name, tmp.DLT_data.data_struct)
            self.assertEqual([(x,) for x in true_charMtrx_data_struct[cell_name]], tmp.DLT_data.data_struct[cell_name])
        # 2. tmp.priors and true_priors
        true_priors = true_priors['Q']
        for i, list_of_dicts in enumerate(true_priors):
            q_dict = list_of_dicts[0]
            for key in q_dict:
                self.assertEqual(q_dict[key], tmp.priors[i][0][int(key)])
        # 3. tmp.CharMtrx.alphabet and true_alphabet_data_struct
        for idx, cassette_alphabet in enumerate(true_alphabet_data_struct):
            set1 = set(cassette_alphabet[0])
            set2 = set(tmp.DLT_data.alphabet.data_struct[idx][0])
            self.assertEqual(set1, set2)
        # 4. tmp.CharMtrx.alphabet.get_cassette_alphabet(k) for k in range(tmp.CharMtrx.K) and true_alphabet_data_struct
        with open(true_cassette_alphabet_file, "r") as r:
            lines = r.readlines()
            for line in lines:
                k, true_cassette_alphabet = ast.literal_eval(line)
                set1 = set(true_cassette_alphabet)
                set2 = set(tmp.DLT_data.alphabet.get_cassette_alphabet(k))
                self.assertEqual(set1, set2)

    def test_15(self):
        """"
        ## Test reading in PMMI. 
        alleletable_file = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_proc_data/example_small_alleletable.json')
        priorfile = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_proc_data/small_priors.csv')
        tmp = DLT_parser()
        self.assertEqual(tmp.max_allele_per_cassette, None)
        tmp = DLT_parser(max_allele_per_cassette=2)
        self.assertEqual(tmp.max_allele_per_cassette, 2)

        DLT_data, Q = tmp.get_from_path(datafile=alleletable_file, priorfile=priorfile, max_allele_per_cassette=1)
        """

        # I introduced a new state in the prior file at site index 3 that isn't observed in the json, should be in the alphabet
        #self.assertEqual(tmp.alphabet.M, [60, 60])
        #self.assertEqual(tmp.max_allele_per_cassette, 1)
        #self.assertEqual(tmp.alphabet.get_site_alphabet(0,0), {-1, 0, 1, 2, 3})
        #self.assertEqual(tmp.alphabet.get_site_alphabet(1,0), {-1, 0, 1, 2, 3})
        #self.assertEqual(tmp.alphabet.get_site_alphabet(1,1), {-1, 0, 1})
        pass
