from .Base_model import Base_model
from .AlleleTable import AlleleTable
from .Param import Param
from math import *
from scipy.stats import multinomial
import time

DEFAULT_max_mu = 10
DEFAULT_max_nu = 1
DEFAULT_min_rho = 0.5

class PMMC_model(Base_model):
    # PMMC = probabilistic mixed-type missing with counts
    def __init__(self,treeList,data,prior,kw_params={}):
    # data is a dictionary; must have 'DLT_data' and data['DLT_data'] must be an instance of the AlleleTable class
        mu = kw_params['mu'] if 'mu' in kw_params else 1 # arbitrarily chosen to set as a default value
        nu = kw_params['nu'] if 'nu' in kw_params else 0 # arbitrarily chosen to set as a default value
        phi = kw_params['phi'] if 'phi' in kw_params else 0.1 # arbitrarily chosen to set as a default value
        rho = kw_params['rho'] if 'rho' in kw_params else 1 # arbitrarily chose to set as a default value
        params = Param(['mu','nu','phi','rho'],[mu,nu,phi,rho],[0,0,0,DEFAULT_min_rho],[DEFAULT_max_mu,DEFAULT_max_nu,1,1])
        super(PMMC_model,self).__init__(treeList,data,prior,params)
    
    def Psi(self,c_node,k,j,alpha,beta):
        # Layer 1: transition probabilities  
        # override the Base_model class
        delta = c_node.edge_length
        nu = self.params.get_value('nu')
        if alpha == 0:
            if beta == 0:
                p = exp(-delta*(1+nu))
            elif beta == -1:
                p = 1-exp(-delta*nu)
            else:
                q = self.Q[k][j][beta]
                p = q*exp(-delta*nu)*(1-exp(-delta))   
        else:
            if alpha == -1 and beta == -1:
                p = 1
            elif alpha == beta:
                p = exp(-delta*nu)
            elif beta == -1:
                p = 1-exp(-delta*nu)
            else:
                p = 0        
        return p

    def Gamma(self,k,x,c):
        # Layer 2: emission probabilities  
        # override the Base_model class
        # x is a cassette state of cassette k (data type: tuple of length J)
        # c is a mapping (cassette state -> count) (data type: a dictionary, where keys are tuples of length J, values are integers)
        # NOTE: assume without checking that x is in c.keys()
        J = self.data['DLT_data'].J # self.data['DLT_data'] is an instance of AlleleTable
        M = self.data['DLT_data'].alphabet.get_M(k)
        x_is_silenced = (x == tuple([-1]*J))
        phi = self.params.get_value('phi')
        rho = self.params.get_value('rho') # will be None if 'rho' not in self.params
        c_is_missing = True
        for y in c:
            if c[y] != 0:
                c_is_missing = False
                break   
        if x_is_silenced:
            p_ans = 1 if c_is_missing else 0
        else:
            if c_is_missing:
                p_ans = phi
            else:
                #start = time.time()
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
                #end = time.time()
                #print("setup: ",end-start)
                #p_ans = (1-phi)*multinomial.pmf(X,p=P,n=N)
                #start = time.time()
                log_p_ans = log(1-phi) + log(factorial(N))
                for X_i,P_i in zip(X,P):
                    if P_i == 0 and X_i != 0:
                        log_p_ans = -float("inf")
                        break
                    if P_i != 0:
                        log_p_ans += X_i*log(P_i)
                    log_p_ans -= log(factorial(X_i))
                p_ans = exp(log_p_ans)
                #end = time.time()
                #print("compute: ",end-start)
        return p_ans
    
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
        
        for tree in self.trees:
            for v in tree.traverse_leaves():
                for k in range(K):
                    if self.data['DLT_data'].is_missing(v.label,k): # p_A = p_B = 0
                        p_A = p_B = 0
                        if silenced_state in v.log_node_posterior[k]:
                            p_C = 1-exp(v.log_node_posterior[k][silenced_state])
                        else:
                            p_C = 1    
                    else: # p_C = 0; p_A + p_B = 1
                        c = self.data['DLT_data'].get(v.label,k) 
                        for x in v.log_node_posterior[k]:
                            w = exp(v.log_node_posterior[k][x])
                            p_A = w*c[x]
                            p_B = w*sum([c[y] for y in c if y != x])
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
