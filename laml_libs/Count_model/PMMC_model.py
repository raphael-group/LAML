from .PMM_base_model import PMM_base_model
from .AlleleTable import AlleleTable
from .Param import Param
from math import *
from scipy.stats import multinomial
import time

DEFAULT_max_mu = 10
DEFAULT_max_nu = 1
DEFAULT_min_rho = 0.5

class PMMC_model(PMM_base_model):
    """
        This class represents the PMMC model as a generative process for the dynamic lineage tracing (DLT) data.
        The DLT data of this model must be an instance of the AlleleTable class
        This class inherits all attributes and methods from PMM_base_model. 
        It has the following parameters for layer 1:
            mu: genome edit rate
            nu: silencing rate
            phi: dropout probability
        Layer 2 is either parameterized by    
            rho: the accuracy of readout in layer 2        
        This class overrides logGamma of PMM_base_model
    """
    # PMMC = probabilistic mixed-type missing with counts
    def __init__(self,treeList,data,prior,kw_params={}):
    # data is a dictionary; must have 'DLT_data' and data['DLT_data'] must be an instance of the AlleleTable class
        mu = kw_params['mu'] if 'mu' in kw_params else 1 # arbitrarily chosen to set as a default value
        nu = kw_params['nu'] if 'nu' in kw_params else 0 # arbitrarily chosen to set as a default value
        phi = kw_params['phi'] if 'phi' in kw_params else 0.1 # arbitrarily chosen to set as a default value
        rho = kw_params['rho'] if 'rho' in kw_params else 1 # arbitrarily chose to set as a default value
        params = Param(['mu','nu','phi','rho'],[mu,nu,phi,rho],[0,0,0,DEFAULT_min_rho],[DEFAULT_max_mu,DEFAULT_max_nu,1,1])
        super(PMMC_model,self).__init__(treeList,data,prior,params)
    
    def logGamma(self,k,x,c):
        """
            Compute the emission probability in Layer 2
            This method overrides that of the Base_model class
                x is a cassette state of cassette k (data type: tuple of length J)
                c is a mapping (cassette state -> count); its data type is a dictionary whose keys are tuples of length J and whose values are integers
        """    
        #J = self.data['DLT_data'].J # self.data['DLT_data'] is an instance of AlleleTable
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
        #x_is_silenced = (x == tuple([-1]*J))
        phi = self.params.get_value('phi')
        rho = self.params.get_value('rho') # will be None if 'rho' not in self.params
        c_is_missing = True
        for y in c:
            if c[y] != 0:
                c_is_missing = False
                break   
        if x_is_silenced:
            log_p_ans = 0 if c_is_missing else None
        else:
            if c_is_missing:
                log_p_ans = log(phi) if phi>0 else None
            elif phi == 1:
                log_p_ans = None
            else:
                P = [rho]
                if x in c:
                    X = [c[x]] 
                    N = c[x]
                else:
                    X = [0]
                    N = 0    
                M_c = 0
                for y in c:
                    if y != x:
                        N += c[y]
                        X.append(c[y])
                        P.append((1-rho)/(M-1))        
                        M_c += 1
                # the cassette states that are not included in c have counts 0
                M_missing =  M-M_c
                while M_missing > 0:
                    X.append(0)
                    P.append((1-rho)/(M-1))
                    M_missing -= 1
                log_p_ans = log(1-phi) + log(factorial(N))
                for X_i,P_i in zip(X,P):
                    if P_i == 0 and X_i != 0:
                        log_p_ans = None
                        break
                    if P_i != 0:
                        log_p_ans += X_i*log(P_i)
                    log_p_ans -= log(factorial(X_i))
        return log_p_ans
    
    def set_closed_form_optimal(self,fixed_params={},verbose=1):
        """
            For every param that has a closed-form M-step optimal, 
            if that param is in fixed_params, set it to the specified fixed value;
            otherwise, compute the closed-form solution and set that to the param's value
            Override the Base_model class. 
            This class (i.e. the PMMC model) has two params with closed-form M-step, which are phi and rho 
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
        A = 0 # x is not silenced, c is missing
        B = 0 # x is not silenced, c is not missing 
        C = 0 # x is not silenced, c is not missing, y = x
        D = 0 # x is not silenced, c is not missing, y != x
        
        K = self.data['DLT_data'].K
        #J = self.data['DLT_data'].J
        
        for tree in self.trees:
            for v in tree.traverse_leaves():
                for k in range(K):
                    p_A = p_B = p_C = p_D = 0
                    p_silenced = 0
                    for x in v.log_node_posterior[k]:
                        if __is_silenced(x):
                            p_silenced += exp(v.log_node_posterior[k][x])
                    if self.data['DLT_data'].is_missing(v.label,k): # c is missing
                        p_A = 1-p_silenced  # p_A is P(c is missing, x is not silenced | D; T,\Theta)
                    else: 
                        c = self.data['DLT_data'].get(v.label,k) # NOTE: c is not missing
                        p_B = 1-p_silenced # p_B is P(c is not missing, x is not silenced | D; T,\Theta)
                        for x in v.log_node_posterior[k]:
                            if __is_silenced(x):
                                continue
                            w = exp(v.log_node_posterior[k][x])
                            x_no_flag = x if self.silence_mechanism == 'convolve' else x[:-1]
                            for y in c:
                                if y == x_no_flag:
                                    p_C += w*c[y] # w*c[y] is P(c is not missing, x is not silenced, y=x)
                                else:    
                                    p_D += w*c[y] # w*c[y] is P(c is not missing, x is not silenced, y!=x)
                    A += p_A
                    B += p_B
                    C += p_C
                    D += p_D

        if 'phi' in fixed_params:
            phi_star = fixed_params['phi']
            if verbose > 0:
                print("Fixed phi to " + str(phi_star))
        else:
            # optimizaztion problem: max_{\phi}Alog\phi + Blog(1-\phi)
            phi_star = A/(A+B)
            if verbose > 0:
                print("Current optimal phi: " + str(phi_star))
        # set the value of phi to phi_star
        self.params.set_value('phi',phi_star)
        
        if 'rho' in fixed_params:
            rho_star = fixed_params['rho']
            if verbose > 0:
                print("Fixed rho to " + str(rho_star))
        else:
            # optimization problem: max_{\rho}C\log\rho + Dlog(1-\rho)
            rho_star = C/(C+D)
            if verbose > 0:
                print("Current optimal rho: " + str(rho_star))
        # set the value of rho to rho_star
        self.params.set_value('rho',rho_star)
