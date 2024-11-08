from . import Count_model
from . import IO_handler
from . import PMM_original
from . import TopoSearch
from . import Utils
from . import mixins
from . import preprocess

# below we define the constants that are shared among all libraries

PROGRAM_NAME = "LAML" 
PROGRAM_AUTHOR = ["Uyen Mai","Gillian Chu","Ben Raphael"]
PROGRAM_LICENSE = "GNU General Public License, version 3"
PROGRAM_VERSION = "2.0.0"
PROGRAM_YEAR = "2023"
PROGRAM_INSTITUTE = "Computer Science Department, Princeton University"
PROGRAM_DESCRIPTION = "LAML: Lineage Analysis via Maximum Likelihood"


###### Hyper/default parameters #########
min_llh = -float("inf") # the minimum log-likelihood value
eps = 1e-10 # epsilon value, usually used as the lower bound for non-zero values
DEFAULT_conv_eps = 1e-6 # convergence threshold (to stop a search algorithm such as EM)
DEFAULT_dmin = 0.005
DEFAULT_dmax = 10
DEFAULT_max_mu = 10
DEFAULT_max_nu = 10
DEFAULT_scipy_options = {'iprint':3,'maxiter':1000}
DEFAULT_EM_options = {'conv_eps':DEFAULT_conv_eps,'max_iter':1000}
#DEFAULT_STRATEGY={'resolve_search_only':False,'only_marked':False,'ultra_constr':False,'fixed_phi':None,'fixed_nu':None,'local_brlen_opt':True,'fixed_params':{}}
DEFAULT_STRATEGY={'resolve_search_only':False,'only_marked':False,'ultra_constr':False,'local_brlen_opt':True,'fixed_params':{},'compute_cache':None}
nni_conv_eps = 1e-15 # additional NNI convergence threshold since conv_eps has several purposes
#chkpt_freq = 10
chkpt_freq = 1
