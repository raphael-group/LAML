import numpy as np
from math import *
import cvxpy as cp
from copy import deepcopy

eps_nu = 1e-5

class Mstep_PMM_base:
    def __init__(self,model_obj,ultra_constr_cache=None,local_brlen_opt=True,nIters=1):
        self.model_obj = model_obj
        self.ultra_constr_cache = ultra_constr_cache
        self.local_brlen_opt = local_brlen_opt
        self.nIters = nIters

        N = model_obj.num_edges
        N_free = N
        if self.local_brlen_opt:
            for tree in model_obj.trees:
                N_free -= len([node for node in tree.traverse_postorder() if node.mark_fixed])
        self.N = N
        self.N_free = N_free
        self.compute_counts()
    
    def compute_counts(self):    
        self.C_z2z = np.zeros(self.N)
        self.C_z2a = np.zeros(self.N)
        self.C_a2a = np.zeros(self.N)
        self.C_za2s = np.zeros(self.N)

        self.d_ini = np.zeros(self.N)
        self._j2i = np.zeros(self.N_free,dtype='int')
        i = 0
        j = 0
        K = self.model_obj.data['DLT_data'].K
        
        for tree in self.model_obj.trees:
            for v in tree.traverse_postorder():
                if not v.is_root(): #and not (v.mark_fixed and self.local_brlen_opt):   
                    for k in range(K):    
                        for x,y in v.log_edge_posterior[k]:
                            w = exp(v.log_edge_posterior[k][(x,y)])
                            for x_j,y_j in zip(x,y):
                                self.C_z2z[i] += w*(x_j == 0 and y_j == 0)
                                self.C_z2a[i] += w*(x_j == 0 and y_j != 0 and y_j != -1)
                                self.C_a2a[i] += w*(x_j == y_j and x_j != 0 and x_j != -1)
                                self.C_za2s[i] += w*(y_j == -1 and x_j != -1)
                    self.d_ini[i] = v.edge_length
                    if not (v.mark_fixed and self.local_brlen_opt):
                        self._j2i[j] = i
                        j += 1
                    i += 1
    
    def optimize_brlen(self,nu,verbose=False): # nu is a single number
    # will be defined in a derived class, which should return an array of self.N_free elements
         pass 
    
    def optimize_nu(self,d,verbose=False): # d is a vector of all branch lengths
        pass # will be defined in a derived class
   
    def _free_to_full_list(self,d_free):
        d_full = deepcopy(self.d_ini)
        for j,x in enumerate(d_free):
            d_full[self._j2i[j]] = x
        return d_full

    def solve(self,fixed_params={},verbose=1):
        # optimize the params that have a closed-form solution
        self.model_obj.set_closed_form_optimal(fixed_params=fixed_params,verbose=verbose)
        nu_star = self.model_obj.params.get_value('nu')
        for r in range(self.nIters):
            if self.N_free > 0:
                if verbose > 0:
                    print("Optimizing branch lengths. Current nu:" + str(nu_star), flush=True)
                #try:
                d_star,status_d = self.optimize_brlen(nu_star,verbose=False)
                #except:
                #    d_star = self.d_ini
                #    status_d = "failure"
                #if status_d == "infeasible": # [Perhaps it's okay to remove this line?] should only happen with local EM 
                #    return False,"d_infeasible"
            else:
                if verbose > 0:
                    print("Fixed branch lengths to the given inputs")
                d_star,status_d = [],"fixed"        
            if 'nu' in fixed_params:
                if verbose > 0:
                    print("Fixed nu to " + str(fixed_params['nu']))
                nu_star = fixed_params['nu']
                status_nu = "optimal"    
            else:    
                d_star_full = self._free_to_full_list(d_star)
                if verbose > 0:
                    print("Optimizing nu", flush=True)
                #try:
                nu_star,status_nu = self.optimize_nu(d_star_full,verbose=False)                 
                #except:
                #    status_nu = "failure"
        return d_star,status_d,nu_star,status_nu

class Mstep_PMMconv(Mstep_PMM_base):
    def optimize_brlen(self,nu,verbose=False): # nu is a single number
        var_d = cp.Variable(self.N_free,nonneg=True) # the branch length variables
        SA = -(nu+1)*self.C_z2z[self._j2i].T @ var_d
        SB = -nu*(self.C_z2a[self._j2i]+self.C_a2a[self._j2i]).T @ var_d
        SC = (self.C_z2a[self._j2i]).T @ cp.log(1-cp.exp(-var_d))
        SD = (self.C_za2s[self._j2i]).T @ cp.log(1-cp.exp(-nu*var_d)) if sum(self.C_za2s) > 0 and nu > eps_nu else 0
        
        objective = cp.Maximize(SA+SB+SC+SD)
        constraints = [np.zeros(self.N_free)+self.model_obj.dmin <= var_d, var_d <= np.zeros(self.N_free)+self.model_obj.dmax]             
        if self.ultra_constr_cache is not None:
            M,b = self.ultra_constr_cache
            constraints += [np.array(M) @ var_d == np.array(b)]
        prob = cp.Problem(objective,constraints)
        prob.solve(verbose=verbose,solver=cp.MOSEK)
        return var_d.value,prob.status
   
    def optimize_nu(self,d,verbose=False): # d is a vector of all branch lengths
        var_nu = cp.Variable(1,nonneg=True) # the nu variable
        SA = -(var_nu+1)*(self.C_z2z).T @ d
        SB = -var_nu*(self.C_z2a+self.C_a2a).T @ d
        SD = self.C_za2s.T @ cp.log(1-cp.exp(-var_nu*d)) if sum(self.C_za2s) > 0 else 0
        
        objective = cp.Maximize(SA+SB+SD)
        prob = cp.Problem(objective)
        prob.solve(verbose=verbose,solver=cp.MOSEK)
        return var_nu.value[0],prob.status

class Mstep_PMMsep(Mstep_PMM_base):
    def compute_counts(self): # override the base class   
        self.C_z2z = np.zeros(self.N)
        self.C_z2a = np.zeros(self.N)
        self.C_omega_z2z = np.zeros(self.N)
        self.C_omega_z2s = np.zeros(self.N)

        self.d_ini = np.zeros(self.N)
        self._j2i = np.zeros(self.N_free,dtype='int')
        i = 0
        j = 0
        K = self.model_obj.data['DLT_data'].K
        for tree in self.model_obj.trees:
            for v in tree.traverse_postorder():
                if not v.is_root(): #and not (v.mark_fixed and self.local_brlen_opt):   
                    for k in range(K):    
                        for x,y in v.log_edge_posterior[k]:
                            w = exp(v.log_edge_posterior[k][(x,y)])
                            self.C_omega_z2z[i] += w*(x[-1]==0 and y[-1]==0)
                            self.C_omega_z2s[i] += w*(x[-1]==0 and y[-1]==-1)
                            for x_j,y_j in zip(x[:-1],y[:-1]):
                                self.C_z2z[i] += w*(x_j == 0 and y_j == 0)
                                self.C_z2a[i] += w*(x_j == 0 and y_j != 0 and y_j != -1)
                    self.d_ini[i] = v.edge_length
                    if not (v.mark_fixed and self.local_brlen_opt):
                        self._j2i[j] = i
                        j += 1
                    i += 1
    def optimize_brlen(self,nu,verbose=False): # nu is a single number
        # override Mstep_PMM_base
        var_d = cp.Variable(self.N,nonneg=True) # the branch length variables
        SA = -self.C_z2z[self._j2i].T @ var_d # [Gillian]: possibly causing errors
        SB = (self.C_z2a[self._j2i]).T @ cp.log(1-cp.exp(-var_d))
        SC = -nu*self.C_omega_z2z[self._j2i].T @ var_d
        SD = (self.C_omega_z2s[self._j2i]).T @ cp.log(1-cp.exp(-nu*var_d)) #if sum(self.C_omega_z2s[self._j2i]) > 0 and nu > eps_nu else 0
        
        objective = cp.Maximize(SA+SB+SC+SD)
        constraints = [np.zeros(self.N)+self.model_obj.dmin <= var_d, var_d <= np.zeros(self.N)+self.model_obj.dmax]             
        if self.ultra_constr_cache is not None:
            M,b = self.ultra_constr_cache
            constraints += [np.array(M) @ var_d == np.array(b)]
        prob = cp.Problem(objective,constraints)
        prob.solve(verbose=verbose,solver=cp.MOSEK)
        return var_d.value,prob.status
   
    def optimize_nu(self,d,verbose=False): # d is a vector of all branch lengths
        var_nu = cp.Variable(1,nonneg=True) # the nu variable
        SC = -var_nu*self.C_omega_z2z.T @ d
        SD = self.C_omega_z2s.T @ cp.log(1-cp.exp(-var_nu*d)) #if sum(self.C_omega_z2s) > 0 else 0
        
        objective = cp.Maximize(SC+SD)
        prob = cp.Problem(objective)
        prob.solve(verbose=verbose,solver=cp.MOSEK)
        return var_nu.value[0],prob.status
