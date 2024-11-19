from .Base_model import Base_model
from .AlleleTable import AlleleTable
from .Param import Param
from math import *

DEFAULT_max_mu = 10
DEFAULT_max_nu = 1
DEFAULT_min_rho = 0.5

class PMM_base_model(Base_model):
    """
    The base class for all models that use PMM as a generative process for the dynamic lineage tracing (DLT) data
    It inherits all attributes and methods from Base_model and only overrides Psi
    """
    def Psi(self,c_node,k,j,alpha,beta):
        """
        Compute Layer 1 transition probabilities  
        This function overrides the Base_model class
        """    
        delta = c_node.edge_length
        nu = self.params.get_value('nu')
        if alpha == 0:
            if beta == 0:
                p = exp(-delta*(1+nu)) if self.silence_mechanism == 'convolve' else exp(-delta)
            elif beta == -1:
                p = 1-exp(-delta*nu)
            else:
                q = self.Q[k][j][beta]
                p = q*exp(-delta*nu)*(1-exp(-delta)) if self.silence_mechanism == 'convolve' else q*(1-exp(-delta)) 
        else:
            if alpha == -1 and beta == -1:
                p = 1
            elif alpha == beta:
                p = exp(-delta*nu) if self.silence_mechanism == 'convolve' else 1
            elif beta == -1:
                p = 1-exp(-delta*nu) # shouldn't happen if self.silence_mechanism == 'separated'
            else:
                p = 0        
        return p
