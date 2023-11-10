# below we define the constants that are shared among all libraries


PROGRAM_NAME = "scmail" # sc-MAIL"
PROGRAM_AUTHOR = ["Uyen Mai","Gillian Chu","Ben Raphael"]
PROGRAM_LICENSE = "GNU General Public License, version 3"
PROGRAM_VERSION = "1.0.0"
PROGRAM_YEAR = "2023"
PROGRAM_INSTITUTE = "Computer Science Department, Princeton University"
PROGRAM_DESCRIPTION = "sc-MAIL: single-cell Maximum-likelihood Ancestral Inference for Lineage-tracing"


###### Hyper/default parameters #########

min_llh = -float("inf") # the minimum log-likelihood value
eps = 1e-10 # epsilon value, usually used as the lower bound for non-zero values
conv_eps = 1e-8 # convergence threshold (to stop a search algorithm such as EM)
nni_conv_eps = 1e-15 # additional NNI convergence threshold since conv_eps has several purposes
DEFAULT_STRATEGY={'resolve_search_only':False,'only_marked':False,'ultra_constr':False,'fixed_phi':None,'fixed_nu':None,'local_brlen_opt':True}
