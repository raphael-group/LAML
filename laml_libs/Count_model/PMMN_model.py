from .Base_model import Base_model
from .AlleleTable import AlleleTable
from .Param import Param
from math import *

DEFAULT_max_mu = 10
DEFAULT_max_nu = 10

class PMMN_model(Base_model):
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
            params = Param(['mu','nu','phi','rho'],[mu,nu,phi,rho],[0,0,0,0],[DEFAULT_max_mu,DEFAULT_max_nu,1,1])
        super(PMMN_model,self).__init__(treeList,data,prior,params)
    
    def Psi(self,c_node,k,j,alpha,brho):
        # Layer 1: transition probabilities  
        # override the Base_model class
        delta = c_node.edge_length
        nu = self.params.get_value('nu')
        if alpha == 0:
            if brho == 0:
                p = exp(-delta*(1+nu))
            elif brho == -1:
                p = 1-exp(-delta*nu)
            else:
                q = self.Q[k][j][brho]
                p = q*exp(-delta*nu)*(1-exp(-delta))   
        else:
            if alpha == -1 and brho == -1:
                p = 1
            elif alpha == brho:
                p = exp(-delta*nu)
            elif brho == -1:
                p = 1-exp(-delta*nu)
            else:
                p = 0        
        return p

    def Gamma(self,k,x,c):
        # Layer 2: emission probabilities  
        # override the Base_model class
        # x is a cassette state of cassette k (data type: tuple of length J)
        # c is a cassette state of cassette k (data type: tuple of length J)
        J = self.data['DLT_data'].J
        M = self.data['DLT_data'].alphabet.get_M(k)
        x_is_silenced = (x == tuple([-1]*J))
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
                    p = (1-phi)*rho if c==x else (1-phi)*(1-rho)/(M-2) 
        return p

    def set_closed_form_optimal(self,fixed_params={},verbose=1):
        # For every param that has a closed-form M-step optimal, 
        # if that param is in fixed_params, set it to the specified fixed value;
        # otherwise, compute the closed-form solution and set that to the param's value
        # Override the Base_model class. 
        # This class (i.e. the PMM model) only has one param with closed-form M-step, which is phi
        if 'phi' in fixed_params:
            phi_star = fixed_params['phi']
            if verbose > 0:
                print("Fixing phi to " + str(phi_star))
        else:
            if verbose > 0:
                print("Optimizing phi")
            R = []
            R_tilde = []
            K = self.data['DLT_data'].K
            J = self.data['DLT_data'].J
            silenced_state = tuple([-1]*J)

            for tree in self.trees:
                for v in tree.traverse_leaves():
                    R.append(sum([not self.data['DLT_data'].is_missing(v.label,k) for k in range(K)]))
                    R_tilde.append(sum([1-exp(v.log_node_posterior[silenced_state]) for k in range(K) if self.data['DLT_data'].is_missing(v.label,k)])) 
            phi_star = sum(R_tilde)/(sum(R)+sum(R_tilde))
            if verbose > 0:
                print("Current phi: " + str(phi_star))
        # set the value of phi to phi_star
        self.params.set_value('phi',phi_star)
