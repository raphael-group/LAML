from .PMM_base_model import PMM_base_model
from .AlleleTable import AlleleTable
from .Param import Param
from math import *

DEFAULT_max_mu = 10
DEFAULT_max_nu = 1
DEFAULT_min_rho = 0.5

class PMMN_model(PMM_base_model):
    """
    This class represents the PMMN model as a generative process for the dynamic lineage tracing (DLT) data.
    The DLT data of this model must be an instance of the CharMtrx class
    This class inherits all attributes and methods from PMM_base_model. 
    It has the following parameters for layer 1:
        mu: genome edit rate
        nu: silencing rate
        phi: dropout probability
    Layer 2 is either parameterized by    
        rho: the accuracy of readout in layer 2        
    or a prior emission matrix. 
    NOTE: rho is a parameter while emission matrix governs hyperparameters. Only one of them should be specified.
    This class overrides logGamma of PMM_base_model
    """
    # PMMN = probabilistic mixed-type missing with noise
    def __init__(self,treeList,data,prior,kw_params={}):
    # data is a dictionary; must have 'DLT_data' and data['DLT_data'] must be an instance of the CharMtrx class
        mu = kw_params['mu'] if 'mu' in kw_params else 1 # arbitrarily chosen to set as a default value
        nu = kw_params['nu'] if 'nu' in kw_params else 0 # arbitrarily chosen to set as a default value
        phi = kw_params['phi'] if 'phi' in kw_params else 0.1 # arbitrarily chosen to set as a default value
        
        if 'emission_mtrx' in prior:
            self.emission_mtrx = prior['emission_mtrx']
            params = Param(['mu','nu','phi'],[mu,nu,phi],[0,0,0],[DEFAULT_max_mu,DEFAULT_max_nu,1])
        else:  
            self.emission_mtrx = None  
            rho = kw_params['rho'] if 'rho' in kw_params else 1 # arbitrarily chose to set as a default value
            params = Param(['mu','nu','phi','rho'],[mu,nu,phi,rho],[0,0,0,DEFAULT_min_rho],[DEFAULT_max_mu,DEFAULT_max_nu,1,1])
        super(PMMN_model,self).__init__(treeList,data,prior,params)

    def logGamma(self,k,x,c):
        """
        Compute the emission probability in Layer 2
        This method overrides that of the Base_model class
            x is a cassette state of cassette k (data type: tuple of length J)
            c is a cassette state of cassette k (data type: tuple of length J)
        """    
        J = self.data['DLT_data'].J
        M = self.data['DLT_data'].alphabet.get_M(k)
        if self.silence_mechanism == 'separated':
            x_is_silenced = (x[-1] == -1)
            x = x[:-1] # remove the silence flag
        else: # self.silence_mechanism is 'convolve'
            x_is_silenced = False
            for x_j in x:
                if x_j == -1:
                    x_is_silenced = True
                    break   
        phi = self.params.get_value('phi')
        rho = self.params.get_value('rho') # will be None if 'rho' not in self.params
        missing_state = tuple(['?']*J)   

        if x_is_silenced:
            p = 1 if c == missing_state else 0
        else:
            if c == missing_state:
                p = phi
            else:
                if self.emission_mtrx is not None:
                    p = self.emission_mtrx[x][c]
                else:
                    p = (1-phi)*rho if c==x else (1-phi)*(1-rho)/(M-1) 
        return log(p) if p>0 else None

    def set_closed_form_optimal(self,fixed_params={},verbose=1):
        """    
        For every param that has a closed-form M-step optimal, 
        if that param is in fixed_params, set it to the specified fixed value;
        otherwise, compute the closed-form solution and set that to the param's value
        This method overrides that of the Base_model class
        This class (i.e. the PMMN model) has two params with closed-form M-step, which are phi and rho
        NOTE: rho will only be optimized if the prior emission matrix is not provided. 
        Otherwise, rho must have been set to None in __init__ and will still be None
        """    
        def __is_silenced(x):
            # "separated" silencing mechanism, we only need to check the flag
            if self.silence_mechanism == 'separated': 
                    return x[-1] == -1
            # "convolved" silencing mechanism, return True if any of the site is silenced
            for x_j in x:
                if x_j == -1:
                    return True
                return False        
        A = 0
        B = 0 
        C = 0
        
        K = self.data['DLT_data'].K
        J = self.data['DLT_data'].J
        #silenced_state = tuple([-1]*J)
        missing_state = tuple(['?']*J)

        for tree in self.trees:
            for v in tree.traverse_leaves():
                for k in range(K):
                    c = self.data['DLT_data'].get(v.label,k)
                    if c == missing_state:
                        p_A = p_B = 0
                        #p_silenced = 0
                        #for x in v.log_node_posterior[k]:
                        #    if __is_silenced(x):
                        #        p_silenced += exp(v.log_node_posterior[k][x])
                        #p_C = 1-p_silenced        
                        p_C = 0
                        for x in v.log_node_posterior[k]:
                            if not __is_silenced(x):
                                p_C += exp(v.log_node_posterior[k][x])
                    else:
                        x = c if self.silence_mechanism == 'convolve' else tuple(list(c)+[0])
                        p_A = exp(v.log_node_posterior[k][x])
                        p_B = 1-p_A
                        if self.silence_mechanism == 'separated':
                            for y in v.log_node_posterior[k]:
                                if __is_silenced(y):
                                    p_B -= exp(v.log_node_posterior[k][y])
                        p_C = 0    
                    A += p_A
                    B += p_B
                    C += p_C

        if 'phi' in fixed_params:
            phi_star = fixed_params['phi']
            if verbose > 0:
                print("Fixed phi to " + str(phi_star))
        else:
            phi_star = C/(A+B+C)
            if verbose > 0:
                print("Current optimal phi: " + str(phi_star))
        # set the value of phi to phi_star
        self.params.set_value('phi',phi_star)
        
        if self.emission_mtrx is None: # only optimize rho if the emission matrix is not provided
            if 'rho' in fixed_params:
                rho_star = fixed_params['rho']
                if verbose > 0:
                    print("Fixed rho to " + str(rho_star))
            else:
                rho_star = A/(A+B)
                if verbose > 0:
                    print("Current optimal rho: " + str(rho_star))
            # set the value of rho to rho_star
            self.params.set_value('rho',rho_star)
