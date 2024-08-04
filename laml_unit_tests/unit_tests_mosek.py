import os 
import unittest
from treeswift import *
from laml_libs.PMM_model.EM_solver import EM_solver
from laml_libs.IO_handler.sequence_lib import read_sequences
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
import cvxpy as cp

class MosekTest(unittest.TestCase):
    def test_1(self): 
        np.random.seed(1)
        x = cp.Variable(pos=True)
        y = cp.Variable(pos=True)
        z = cp.Variable(pos=True)

        objective_fn = x * y * z
        constraints = [
          4 * x * y * z + 2 * x * z <= 10, x <= 2*y, y <= 2*x, z >= 1]
        problem = cp.Problem(cp.Maximize(objective_fn), constraints)
        problem.solve(gp=True, solver=cp.MOSEK)

