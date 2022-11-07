#! /usr/bin/env python

import numpy as np
import sys


M = np.loadtxt(sys.argv[1], usecols=[1,2,3,4,5,6,7])
print(M.mean(0))
