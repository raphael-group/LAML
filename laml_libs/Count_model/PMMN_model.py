from .PMM_base_model import PMM_base_model
from .AlleleTable import AlleleTable
from .Param import Param
from math import *

DEFAULT_max_mu = 10
DEFAULT_max_nu = 1
DEFAULT_min_rho = 0.5

class PMMN_model(PMM_base_model):
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
    
    def Gamma(self,k,x,c,node_label):
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
                    # add node_label, etc..
                else:
                    p = (1-phi)*rho if c==x else (1-phi)*(1-rho)/(M-2) 
        return p

    def set_closed_form_optimal(self,fixed_params={},verbose=1):
        # For every param that has a closed-form M-step optimal, 
        # if that param is in fixed_params, set it to the specified fixed value;
        # otherwise, compute the closed-form solution and set that to the param's value
        # Override the Base_model class. 
        # This class (i.e. the PMMN model) has two params with closed-form M-step, which are phi and rho
        A = 0
        B = 0 
        C = 0
        
        K = self.data['DLT_data'].K
        J = self.data['DLT_data'].J
        silenced_state = tuple([-1]*J)
        missing_state = tuple(['?']*J)

        for tree in self.trees:
            for v in tree.traverse_leaves():
                for k in range(K):
                    observed_state = self.data['DLT_data'].get(v.label,k) 
                    if observed_state == missing_state: # p_A = p_B = 0
                        p_A = p_B = 0
                        if silenced_state in v.log_node_posterior[k]:
                            p_C = 1-exp(v.log_node_posterior[k][silenced_state])
                        else:
                            p_C = 1    
                    else: # p_C = 0; p_A + p_B = 1
                        p_A = exp(v.log_node_posterior[k][observed_state])
                        p_B = 1-p_A
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
