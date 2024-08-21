# below we define the constants that are shared among all libraries

PROGRAM_NAME = "LAML" 
PROGRAM_AUTHOR = ["Uyen Mai","Gillian Chu","Ben Raphael"]
PROGRAM_LICENSE = "GNU General Public License, version 3"
PROGRAM_VERSION = "0.0.4"
PROGRAM_YEAR = "2023"
PROGRAM_INSTITUTE = "Computer Science Department, Princeton University"
PROGRAM_DESCRIPTION = "LAML: Lineage Analysis via Maximum Likelihood"


###### Hyper/default parameters #########
min_llh = -float("inf") # the minimum log-likelihood value
eps = 1e-10 # epsilon value, usually used as the lower bound for non-zero values
DEFAULT_conv_eps = 1e-8 # convergence threshold (to stop a search algorithm such as EM)
DEFAULT_dmin = 0.005
DEFAULT_dmax = 10
DEFAULT_max_mu = 10
DEFAULT_max_nu = 10
DEFAULT_scipy_options = {'iprint':3,'maxiter':1000}
DEFAULT_EM_options = {'conv_eps':DEFAULT_conv_eps,'max_iter':1000}
DEFAULT_STRATEGY={'resolve_search_only':False,'only_marked':False,'ultra_constr':False,'fixed_phi':None,'fixed_nu':None,'local_brlen_opt':True}
nni_conv_eps = 1e-15 # additional NNI convergence threshold since conv_eps has several purposes
chkpt_freq = 10
