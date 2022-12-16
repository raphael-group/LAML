from treeswift import *
from math import log,exp,sqrt
from random import random, seed
from scipy import optimize
import warnings
import numpy as np

min_llh = -1000
eps = 1e-10

class Params:
    def __init__(self,nwkTree,nu=eps,phi=eps,sigma=None):
    # nwkTree: a newick tree (i.e. a string) with branch lengths
    # nu, phi: positive float numbers
        self.tree = read_tree_newick(nwkTree) # store a treewfit object in self.tree
        self.nu = nu
        self.phi = phi
        self.sigma = sigma

class ML_solver:
    # at this stage, the tree topology must be given. Only branch lengths
    # and other parameters can be optimized
    def __init__(self,charMtrx,Q,nwkTree,nu=eps,phi=eps):
        self.charMtrx = charMtrx
        self.Q = Q
        self.params = Params(nwkTree,nu=nu,phi=phi)
        self.numsites = len(self.charMtrx[next(iter(self.charMtrx.keys()))])
        self.num_edges = len(list(self.params.tree.traverse_postorder()))
   
    def az_partition(self,params):
    # partition the tree into edge-distjoint alpha-clades and z-branches
    # there is a different partition for each target-site
        numsites = len(self.charMtrx[next(iter(self.charMtrx.keys()))])
        # a-z decomposition
        # z-branches are given tag 'z'
        # each of other branches is given a tag 
        # alpha where alpha is the alpha-tree it belongs to
        for node in params.tree.traverse_postorder():
            if node.is_leaf():
                node.alpha = [self.charMtrx[node.label][site] if self.charMtrx[node.label][site] != 0 else 'z' for site in range(numsites)]
            else:
                C = node.children
                node.alpha = [None]*numsites
                for site in range(numsites):
                    S = set(c.alpha[site] for c in C)
                    R = S-set(['z','?'])
                    if 'z' in S or len(R)>1:
                        node.alpha[site] = 'z'
                    elif len(R) == 1:
                        node.alpha[site] = list(R)[0]    
                    else:
                        node.alpha[site] = "?"
    
    def lineage_llh(self,params):
        # assume az_partition has been performed so
        # each node has the attribute node.alpha
        numsites = len(self.charMtrx[next(iter(self.charMtrx.keys()))])
        phi = params.phi
        nu = params.nu
        llh = [0]*numsites
        for node in params.tree.traverse_postorder():
            p = exp(-node.edge_length)
            node.L0 = [0]*numsites # L0 and L1 are stored in log-scale
            node.L1 = [0]*numsites
            for site in range(numsites):    
                if node.alpha[site] != 'z':
                    q = self.Q[site][node.alpha[site]] if node.alpha[site] != "?" else 1.0
                    if node.is_leaf():
                        if node.alpha[site] == "?":         
                            #print("masked clade",node.label,site)                  
                            masked_llh = log(1-(1-phi)*p**nu) #if (phi>0 or nu>0) else min_llh
                            node.L0[site] = node.L1[site] = masked_llh
                        else:    
                            node.L0[site] = nu*(-node.edge_length) + log(1-p) + log(q) + log(1-phi)
                            node.L1[site] = nu*(-node.edge_length) + log(1-phi)
                    else:
                        C = node.children
                        l0 = l1 = 0
                        for c in C:
                            l0 += c.L0[site]
                            l1 += c.L1[site]
                        L0 = exp(l0+(nu+1)*(-node.edge_length)) + exp(l1 + log(1-p)+log(q) + nu*(-node.edge_length)) + (1-p**nu)*int(node.alpha[site]=="?")   
                        L1 = exp(l1+nu*(-node.edge_length)) + (1-p**nu)*int(node.alpha[site]=="?")
                        node.L0[site] = min_llh if L0==0 else log(L0)
                        node.L1[site] = min_llh if L1==0 else log(L1)
                            
                    if node.is_root() or node.parent.alpha[site] == 'z':
                        llh[site] += node.L0[site]
                else:
                    llh[site] += (-node.edge_length + int(node.is_leaf())*log(1-phi))
        return sum(llh)         

    def ini_brlens(self):
        dmax = -log(1/self.numsites)*2
        dmin = -log(1-1/self.numsites)/2
        return [random() * (dmax/2 - 2*dmin) + 2*dmin for i in range(self.num_edges)]        

    def ini_nu(self,fixed_nu=None):
        return random()*0.99 if fixed_nu is None else fixed_nu
        
    def ini_phi(self,fixed_phi=None):
        return random()*0.99 if fixed_phi is None else fixed_phi   

    def ini_all(self,fixed_phi=None,fixed_nu=None):
        return self.ini_brlens() + [self.ini_nu(fixed_nu=fixed_nu),self.ini_phi(fixed_phi=fixed_phi)]

    def bound_nu(self,fixed_nu=None):
        return (eps,10) if fixed_nu is None else (fixed_nu-eps,fixed_nu+eps)
    
    def bound_phi(self,fixed_phi=None):
        return (eps,0.99) if fixed_phi is None else (fixed_phi-eps,fixed_phi+eps)

    def bound_brlen(self):        
        dmax = -log(1/self.numsites)*2
        dmin = -log(1-1/self.numsites)/2
        return [dmin]*self.num_edges,[dmax]*self.num_edges
        
    def get_bound(self,keep_feasible=False,fixed_phi=None,fixed_nu=None):
        br_lower,br_upper = self.bound_brlen()  
        phi_lower,phi_upper = self.bound_phi(fixed_phi=fixed_phi)
        nu_lower,nu_upper = self.bound_nu(fixed_nu=fixed_nu)
        bounds = optimize.Bounds(br_lower+[nu_lower,phi_lower],br_upper+[nu_upper,phi_upper],keep_feasible=keep_feasible)
        return bounds

    def x2brlen(self,x):
        for i, node in enumerate(self.params.tree.traverse_postorder()):
            node.edge_length = x[i]

    def x2nu(self,x,fixed_nu=None):
        self.params.nu = x[self.num_edges] if fixed_nu is None else fixed_nu
    
    def x2phi(self,x,fixed_phi=None):
        self.params.phi = x[self.num_edges+1] if fixed_phi is None else fixed_phi

    def x2params(self,x,fixed_nu=None,fixed_phi=None):
        self.x2brlen(x)
        self.x2nu(x,fixed_nu=fixed_nu)
        self.x2phi(x,fixed_phi=fixed_phi)

    def __llh__(self):
        return self.lineage_llh(self.params)

    def negative_llh(self):
        self.az_partition(self.params)
        return -self.__llh__()

    def show_params(self):
        print("tree: " + self.params.tree.newick())
        print("nu: " + str(self.params.nu))
        print("phi: " + str(self.params.phi))
        print("negative-llh: " + str(self.negative_llh()))

    def optimize(self,initials=20,fixed_phi=None,fixed_nu=None,verbose=True,max_trials=100):
    # optimize tree branch lengths and nu and phi 
        self.az_partition(self.params)
        warnings.filterwarnings("ignore")
        def nllh(x): 
            self.x2params(x,fixed_nu=fixed_nu,fixed_phi=fixed_phi)            
            return -self.__llh__()
        
        bounds = self.get_bound(fixed_phi=fixed_phi,fixed_nu=fixed_nu)
        
        x_star = None
        f_star = float("inf")

        x0 = []
        all_failed = True
        all_trials = 0
        while all_failed and all_trials < max_trials:
            for i in range(initials):
                randseed = int(random()*10000)
                print("Initial point " + str(i+1) + ". Random seed: " + str(randseed))
                seed(a=randseed)
                x0 = self.ini_all(fixed_phi=fixed_phi,fixed_nu=fixed_nu)
                out = optimize.minimize(nllh, x0, method="SLSQP", options={'disp':verbose,'iprint':3,'maxiter':1000}, bounds=bounds)
                if out.success:
                    all_failed = False
                    print("Optimal point found for initial " + str(i+1))
                    self.x2params(out.x,fixed_phi=fixed_phi,fixed_nu=fixed_nu)
                    self.show_params()
                    if out.fun < f_star:
                        x_star = out.x
                        f_star = out.fun
                else:
                    print("Failed to optimize using initial point " + str(i+1)) 
            all_trials += initials
        
        # store the optimal values in x_star to self.params
        self.x2params(x_star,fixed_phi=fixed_phi,fixed_nu=fixed_nu)    
        
        return self.negative_llh()     

class SpaLin_solver(ML_solver):
    # at this stage, the tree topology and sig,a must be given. Only branch lengths
    # and other parameters can be optimized
    def __init__(self,charMtrx,Q,nwkTree,locations,sigma,nu=eps,phi=eps):
        super(SpaLin_solver,self).__init__(charMtrx,Q,nwkTree,nu=nu,phi=phi)
        self.given_locations = locations
        self.params.sigma = sigma
        self.inferred_locations = {}
        for x in self.given_locations:
            self.inferred_locations[x] = self.given_locations[x]
  
    def spatial_llh(self,locations):
        llh = 0
        for node in self.params.tree.traverse_preorder():
            if node.is_root() or not node.label in locations or not node.parent.label in locations:
                continue
            d = node.edge_length
            curr_sigma = self.params.sigma*sqrt(d)
            x,y = locations[node.label]
            x_par,y_par = locations[node.parent.label]
            llh -= (0.5*((x-x_par)/curr_sigma)**2 + log(curr_sigma))
            llh -= (0.5*((y-y_par)/curr_sigma)**2 + log(curr_sigma))
        return llh 

    def __llh__(self):
        return self.lineage_llh(self.params) + self.spatial_llh(self.inferred_locations)
    
    def ini_all(self,fixed_phi=None,fixed_nu=None):
        x_lin = self.ini_brlens() + [self.ini_nu(fixed_nu=fixed_nu),self.ini_phi(fixed_phi=fixed_phi)]
        x_spa = []
        for node in self.params.tree.traverse_postorder():
            if not node.label in self.given_locations:
                x_spa += [random(),random()]
        x_sigma = 22 # hard code for now        
        return x_lin + x_spa + [x_sigma]
    
    def x2params(self,x,fixed_nu=None,fixed_phi=None):
        self.x2brlen(x)
        self.x2nu(x,fixed_nu=fixed_nu)
        self.x2phi(x,fixed_phi=fixed_phi)
        i = self.num_edges + 2
        for node in self.params.tree.traverse_postorder():
            if not node.label in self.given_locations:
                self.inferred_locations[node.label] = (x[i],x[i+1])
                i += 2
        self.params.sigma = x[-1]        
               
    def bound_locations(self,lower=-np.inf,upper=np.inf):
        N = 2*len([node for node in self.params.tree.traverse_postorder() if not node.label in self.given_locations])    
        return [lower]*N,[upper]*N

    def bound_sigma(self):
        return (eps,np.inf)    

    def get_bound(self,keep_feasible=False,fixed_phi=None,fixed_nu=None):
        br_lower,br_upper = self.bound_brlen()  
        phi_lower,phi_upper = self.bound_phi(fixed_phi=fixed_phi)
        nu_lower,nu_upper = self.bound_nu(fixed_nu=fixed_nu)
        spa_lower,spa_upper = self.bound_locations()
        sigma_lower,sigma_upper = self.bound_sigma()
        bounds = optimize.Bounds(br_lower+[nu_lower,phi_lower]+spa_lower+[sigma_lower],br_upper+[nu_upper,phi_upper]+spa_upper+[sigma_upper],keep_feasible=keep_feasible)
        return bounds

def main(): 
    from sequence_lib import read_sequences
    from ml_log import wrapper_felsenstein as wf_log
    
    k = 50
    m = 10
    Q = []
    for i in range(k):
        q = {j+1:1/m for j in range(m)}
        q[0] = 0
        Q.append(q)
    T = "((a:0.0360971597765934,b:3.339535381892265)e:0.0360971597765934,(c:0.0360971597765934,d:3.339535381892265)f:0.0360971597765934)r:0.0;"
    #T = "((a,b)e,(c,d)f)r;"
    #T = "((a:1,b:1):1,c:1):1;"
    S = read_sequences("../Experiments/MP_inconsistent/seqs_m10_k" + str(k) + ".txt",filetype="fasta")
    msa = S[6]
    #msa['d'][0] = '?'
    #msa['b'][0] = '?'
    #msa['c'][0] = '?'
    #msa['a'][0] = '?'
    #msa = {'a':[0],'b':[0],'c':['?']}
    #print(wf_log(T, Q, msa, optimize_branchlengths=True,initials=1))

    mySolver = ML_solver(msa,Q,T,nu=eps,phi=eps)
    mySolver.az_partition(mySolver.params)
    print(mySolver.negative_llh())
    #print(mySolver.optimize(initials=1,fixed_phi=eps,fixed_nu=eps,verbose=True))
    #print(mySolver.params.phi,mySolver.params.nu)
    #print(mySolver.params.tree.newick())

if __name__ == "__main__":
    main()        
