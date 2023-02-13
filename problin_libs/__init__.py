# below we define the constants that are shared among all libraries

min_llh = -1000 # the minimum log-likelihood value
eps = 1e-10 # epsilon value, usually used as the lower bound for non-zero values
conv_eps = 1e-8 # convergence threshold (to stop a search algorithm such as EM)

PROGRAM_NAME = "Problin"
PROGRAM_AUTHOR = ["Uyen Mai","Gillian Chu","Ben Raphael"]
PROGRAM_LICENSE = "GNU General Public License, version 3"
PROGRAM_VERSION = "0.1"
PROGRAM_YEAR = "2023"
PROGRAM_INSTITUTE = "Computer Science Department, Princeton University"
PROGRAM_DESCRIPTION = "Probabilistic Lineage Tracing with Mixed-type Missing Data"
nni_conv_eps = 1e-15 # additional NNI convergence threshold since conv_eps has several purposes
