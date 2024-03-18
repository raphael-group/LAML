import sys
import os 
import unittest
import cvxpy as cp
import pkg_resources
import subprocess
import re 
from random import randint
import shutil

# Functional tests
class InputOutputTest(unittest.TestCase):

    def test_1(self):

        out = subprocess.run(["run_laml", "--help"], capture_output=True, text=True)
        self.assertEqual(out.returncode, 0)

    def test_2(self):
        ### run with character matrix with column of missing 
        tmp_output_dir = f'tmp_{randint(10**2, 10**3-1)}'
        os.mkdir(tmp_output_dir)
        # python run_laml.py -c /Users/gc3045/scmail_v1/laml/laml_unit_tests/test_data/test_inputs/test1_charMtrx.txt -t /Users/gc3045/scmail_v1/laml/laml_unit_tests/test_data/test_inputs/test1.tre
        treedata_path = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_inputs/test1.tre')
        msa_path = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_inputs/test1_charMtrx.txt')
       
        rand_file_name = os.path.join(tmp_output_dir, f'LAML_output_{randint(10**2, 10**3-1)}') 
        # waits for the process to end
        out = subprocess.run(["run_laml", "-c", msa_path, "-t", treedata_path, "-o", rand_file_name], capture_output=True, text=True) 
        
        self.assertEqual(bool(re.search('uniform', out.stdout)), True)
        self.assertEqual(out.returncode, 0)

        # check that the four files are generated
        param_file = f'{rand_file_name}_params.txt'
        os.path.isfile(param_file)
        tree_file = f'{rand_file_name}_trees.nwk'
        os.path.isfile(tree_file)
        log_file = f'{rand_file_name}.log'
        os.path.isfile(log_file)
        annotations_file = f'{rand_file_name}_annotations.txt'
        os.path.isfile(annotations_file)

        shutil.rmtree(tmp_output_dir)
       
    def test_3(self):
        ### check --noSilence
        tmp_output_dir = f'tmp_{randint(10**2, 10**3-1)}'
        os.mkdir(tmp_output_dir)

        # run_laml -c examples/example1/character_matrix.csv -t examples/example1/starting.tree -p examples/example1/priors.csv -o example1 --nInitials 1 --noSilence
        treedata_path = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_inputs/example1/starting.tree')
        msa_path = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_inputs/example1/character_matrix.csv')
        priors_path = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_inputs/example1/priors.csv')
       
        rand_file_name = os.path.join(tmp_output_dir, f'LAML_output_{randint(10**2, 10**3-1)}')
        out = subprocess.run(["run_laml", "-c", msa_path, "-t", treedata_path, "-p", priors_path, "-o", rand_file_name, "--nInitials", "1", "--noSilence"], capture_output=True, text=True) 
       
        # check tree is scaled to 1
        pattern = "Tree height after scaling:(.*), mutation"
        matches = re.findall(pattern, out.stdout, re.DOTALL)
        tree_height = float(matches[0])
        self.assertAlmostEqual(tree_height,1.0,places=4,msg="InputOutputTest: test_3 failed.")

        # check silencing rate is close to 0
        pattern = "Optimal nu:(.*). Optimal nllh"
        matches = re.findall(pattern, out.stdout, re.DOTALL)
        silencing_rate = float(matches[-1])
        self.assertAlmostEqual(silencing_rate,0.0,places=4,msg="InputOutputTest: test_3 failed.")
        
        # check dropout rate is >0
        pattern = "Optimal phi:(.*). Optimal nu"
        matches = re.findall(pattern, out.stdout, re.DOTALL)
        dropout_rate = float(matches[-1])
        self.assertGreater(dropout_rate,0.0,msg="InputOutputTest: test_3 failed.")
        
        shutil.rmtree(tmp_output_dir)

    def test_4(self):
        ### check --noSilence --timescale 10
        tmp_output_dir = f'tmp_{randint(10**2, 10**3-1)}'
        os.mkdir(tmp_output_dir)

        treedata_path = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_inputs/example1/starting.tree')
        msa_path = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_inputs/example1/character_matrix.csv')
        priors_path = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_inputs/example1/priors.csv')
       
        rand_file_name = os.path.join(tmp_output_dir, f'LAML_output_{randint(10**2, 10**3-1)}')
        out = subprocess.run(["run_laml", "-c", msa_path, "-t", treedata_path, "-p", priors_path, "-o", rand_file_name, "--nInitials", "1", "--noSilence", "--timescale", "10"], capture_output=True, text=True) 
       
        # check tree is scaled to 10
        pattern = "Tree height after scaling:(.*), mutation"
        matches = re.findall(pattern, out.stdout, re.DOTALL)
        tree_height = float(matches[0])
        self.assertAlmostEqual(tree_height,10,places=4,msg="InputOutputTest: test_4 failed.")

        # check silencing rate is close to 0
        pattern = "Optimal nu:(.*). Optimal nllh"
        matches = re.findall(pattern, out.stdout, re.DOTALL)
        silencing_rate = float(matches[-1])
        self.assertAlmostEqual(silencing_rate,0.0,places=4,msg="InputOutputTest: test_4 failed.")
        
        # check dropout rate is >0
        pattern = "Optimal phi:(.*). Optimal nu"
        matches = re.findall(pattern, out.stdout, re.DOTALL)
        dropout_rate = float(matches[-1])
        self.assertGreater(dropout_rate,0.0,msg="InputOutputTest: test_4 failed.")
        
        shutil.rmtree(tmp_output_dir)

    def test_5(self):
        ### check --timescale 10
        tmp_output_dir = f'tmp_{randint(10**2, 10**3-1)}'
        os.mkdir(tmp_output_dir)

        treedata_path = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_inputs/example1/starting.tree')
        msa_path = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_inputs/example1/character_matrix.csv')
        priors_path = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_inputs/example1/priors.csv')
       
        rand_file_name = os.path.join(tmp_output_dir, f'LAML_output_{randint(10**2, 10**3-1)}')
        out = subprocess.run(["run_laml", "-c", msa_path, "-t", treedata_path, "-p", priors_path, "-o", rand_file_name, "--nInitials", "1", "--timescale", "10"], capture_output=True, text=True) 
       
        # check tree is scaled to 10
        pattern = "Tree height after scaling:(.*), mutation"
        matches = re.findall(pattern, out.stdout, re.DOTALL)
        tree_height = float(matches[0])
        self.assertAlmostEqual(tree_height,10,places=4,msg="InputOutputTest: test_5 failed.")

        # check silencing rate is >0
        pattern = "Optimal nu:(.*). Optimal nllh"
        matches = re.findall(pattern, out.stdout, re.DOTALL)
        silencing_rate = float(matches[-1])
        self.assertGreater(silencing_rate,0.0,msg="InputOutputTest: test_5 failed.")
        
        # check dropout rate is >0
        pattern = "Optimal phi:(.*). Optimal nu"
        matches = re.findall(pattern, out.stdout, re.DOTALL)
        dropout_rate = float(matches[-1])
        self.assertGreater(dropout_rate,0.0,msg="InputOutputTest: test_5 failed.")
        
        shutil.rmtree(tmp_output_dir)

    def test_6(self):
        ### check --noDropout
        tmp_output_dir = f'tmp_{randint(10**2, 10**3-1)}'
        os.mkdir(tmp_output_dir)

        treedata_path = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_inputs/example1/starting.tree')
        msa_path = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_inputs/example1/character_matrix.csv')
        priors_path = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_inputs/example1/priors.csv')
       
        rand_file_name = os.path.join(tmp_output_dir, f'LAML_output_{randint(10**2, 10**3-1)}')
        out = subprocess.run(["run_laml", "-c", msa_path, "-t", treedata_path, "-p", priors_path, "-o", rand_file_name, "--nInitials", "1", "--noDropout"], capture_output=True, text=True) 
       
        # check tree is scaled to 1
        pattern = "Tree height after scaling:(.*), mutation"
        matches = re.findall(pattern, out.stdout, re.DOTALL)
        tree_height = float(matches[0])
        self.assertAlmostEqual(tree_height,1.0,places=4,msg="InputOutputTest: test_6 failed.")
        
        # check silencing rate is >0
        pattern = "Optimal nu:(.*). Optimal nllh"
        matches = re.findall(pattern, out.stdout, re.DOTALL)
        silencing_rate = float(matches[-1])
        self.assertGreater(silencing_rate,0.0,msg="InputOutputTest: test_6 failed.")

        # check dropout rate is 0
        pattern = "Optimal phi:(.*). Optimal nu"
        matches = re.findall(pattern, out.stdout, re.DOTALL)
        dropout_rate = float(matches[-1])
        self.assertAlmostEqual(dropout_rate,0.0,places=4,msg="InputOutputTest: test_6 failed.")
        
        shutil.rmtree(tmp_output_dir)
                

        



