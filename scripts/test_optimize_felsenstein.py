import numpy as np
from scipy import optimize
import math
from math import log, exp

import sys
sys.path.append("/Users/gillianchu/raphael/repos/problin/problin_libs")
# from distance_based_lib import ML_pairwise_estimate
from ml import wrapper_felsenstein

T = '[&R] ((0,1),(2,3));'

msa = np.array([[1,0,2], # species 1
                [1,0,0], # species 2
                [0,1,2], # species 3
                [0,1,0]] # species 4
            )

# prob of transitioning 0 -> 0, 0 -> 1, 0 -> 2
Q = [{0: 0.2, 1:0.3, 2:0.5}, # site 1
     {0: 0.2, 1:0.3, 2:0.5}, # site 2
     {0: 0.2, 1:0.3, 2:0.5}, # site 3
    ]
wrapper_felsenstein(T, Q, msa, 3, root_edge_len=0.2)