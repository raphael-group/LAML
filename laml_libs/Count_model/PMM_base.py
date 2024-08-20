from .Base_model import *
from math import *

DEFAULT_max_mu = 10
DEFAULT_max_nu = 10

class PMM_model(Count_base_model):
    # PMM = probabilistic mixed-type missing
    def __init__(self,treeList,data,prior,**params):
    # NOTE: **params must have 'mu', 'nu', 'phi'
        params = Param(['mu','nu','phi'],[params['mu'],params['nu'],params['phi']],[0,0,0],[DEFAULT_max_mu,DEFAULT_max_nu,1])
        super(PMM_model,self).__init__(treeList,data,prior,params)
    
    def Psi(self,c_node,k,j,alpha,beta):
        # Layer 1 transition probabilities  
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
        # Layer 2 transition probabilities  
        # override the Base_model class
        # x is a cassette state of length J corresponding to cassette k
        J = self.data['alleleTable'].J
        x_is_silenced = (x == tuple([-1]*J))
        phi = self.params.get_value('phi')
        max_count = 0
        for y in c:
            if c[y] > max_count:
                max_count = c[y]
        if max_count == 0: # all counts in c are 0, dropout
            p = 1 if x_is_silenced else phi
        else:
            p = 1-phi if (c[x] == max_count and not x_is_silenced) else 0
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
            K = self.data['alleleTable'].K
            J = self.data['alleleTable'].J
            silenced_state = tuple([-1]*J)

            for tree in self.trees:
                for v in tree.traverse_leaves():
                    R.append(sum([not self.data['alleleTable'].is_missing(v.label,k) for k in range(K)]))
                    R_tilde.append(sum([1-exp(v.log_node_posterior[silenced_state]) for k in range(K) if self.data['alleleTable'].is_missing(v.label,k)])) 
            phi_star = sum(R_tilde)/(sum(R)+sum(R_tilde))
            if verbose > 0:
                print("Current phi: " + str(phi_star))
        # set the value of phi to phi_star
        self.params.set_value('phi',phi_star)
