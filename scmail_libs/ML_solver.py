from treeswift import *
from math import log,exp,sqrt, isclose
from random import random, seed, choice
from scipy import optimize
import warnings
import numpy as np
from scmail_libs import min_llh, eps, nni_conv_eps
from scmail_libs.Virtual_solver import Virtual_solver
from scipy.sparse import csr_matrix
from copy import deepcopy
from scmail_libs.lca_lib import find_LCAs

class Params:
    def __init__(self,nu,phi):
        self.nu = nu
        self.phi = phi

class ML_solver(Virtual_solver):
    def __init__(self,treeList,data,prior,params={'nu':eps,'phi':eps}):
        charMtrx = data['charMtrx']
        Q = prior['Q']
        nu = params['nu']
        phi = params['phi']
        self.charMtrx = charMtrx
        self.trees = []
        self.num_edges = 0
        for tree in treeList:
            tree_obj = read_tree_newick(tree)
            tree_obj.suppress_unifurcations()
            self.num_edges += len(list(tree_obj.traverse_postorder()))
            self.trees.append(tree_obj)
        # normalize Q
        self.Q = []
        for Q_i in Q:
            s = sum([Q_i[x] for x in Q_i])
            Q_i_norm = {x:Q_i[x]/s for x in Q_i}
            self.Q.append(Q_i_norm)        
        # setup params
        self.params = Params(nu,phi)        
        # compute numsites, num_edges, dmin, and dmax 
        self.numsites = len(self.charMtrx[next(iter(self.charMtrx.keys()))])
        self.dmin = 0.005
        zerocount = sum([self.charMtrx[e].count(0) for e in self.charMtrx])
        totalcount = self.numsites * len(self.charMtrx)
        zeroprop = zerocount/totalcount
        self.dmax = -log(zeroprop) if zeroprop != 0 else 10

    def get_tree_newick(self):
        return [tree.newick() for tree in self.trees]

    def get_params(self):
        return {'phi':self.params.phi,'nu':self.params.nu}

    def ultrametric_constr(self):
        N = len(self.ini_all())
        M = []
        idx = 0
        for tree in self.trees: 
            for node in tree.traverse_postorder():
                if node.is_leaf():
                    node.constraint = [0.]*N
                else:
                    c1,c2 = node.children
                    m = [x-y for (x,y) in zip(c1.constraint,c2.constraint)]
                    M.append(m)
                    node.constraint = c1.constraint
                node.constraint[idx] = 1
                idx += 1
        for tree in self.trees[1:]:
            m = [x-y for (x,y) in zip(self.trees[0].root.constraint,tree.root.constraint)]
            M.append(m)
        return M

    def score_tree(self,strategy={'ultra_constr':False,'fixed_phi':None,'fixed_nu':None,'fixed_brlen':{}}):
        ultra_constr = strategy['ultra_constr']
        fixed_phi = strategy['fixed_phi']
        fixed_nu = strategy['fixed_nu']
        fixed_brlen = strategy['fixed_brlen']
        nllh,status = self.optimize(initials=1,verbose=-1,ultra_constr=ultra_constr,fixed_phi=fixed_phi,fixed_nu=fixed_nu,fixed_brlen=fixed_brlen)
        score = None if nllh is None else -nllh
        #if score is None:
        #    print("Fatal error: failed to score tree " + self.get_tree_newick() + ". Optimization status: " + status)
        return score,status

    def az_partition(self):
    # Purpose: partition the tree into edge-distjoint alpha-clades and z-branches
    # Note: there is a different partition for each target-site
    # Output: annotate each node of the tree by node.alpha
        # z-branches are given tag 'z'
        # each of other branches is given a tag 
        # alpha where alpha is the alpha-tree it belongs to
        for tree in self.trees:
            for node in tree.traverse_postorder():
                if node.is_leaf():
                    node.alpha = [None]*self.numsites
                    for site in range(self.numsites):
                        if self.charMtrx[node.label][site] == 0:
                            node.alpha[site] = 'z'
                        elif self.charMtrx[node.label][site] == -1:   
                            node.alpha[site] = '?' 
                        else:
                            node.alpha[site] = self.charMtrx[node.label][site]
                else:
                    C = node.children
                    node.alpha = [None]*self.numsites
                    for site in range(self.numsites):
                        S = set(c.alpha[site] for c in C)
                        R = S-set(['z','?',-1])
                        if 'z' in S or len(R)>1:
                            node.alpha[site] = 'z'
                        elif len(R) == 1:
                            node.alpha[site] = list(R)[0]    
                        else:
                            node.alpha[site] = "?"
    
    def lineage_llh(self):
        # assume az_partition has been performed so
        # each node has the attribute node.alpha
        phi = self.params.phi
        nu = self.params.nu
        llh = [0]*self.numsites
        for tree in self.trees:
            for node in tree.traverse_postorder():
                p = exp(-node.edge_length)
                node.L0 = [0]*self.numsites # L0 and L1 are stored in log-scale
                node.L1 = [0]*self.numsites
                for site in range(self.numsites):    
                    if node.alpha[site] != 'z':
                        q = self.Q[site][node.alpha[site]] if node.alpha[site] != "?" else 1.0
                        if node.is_leaf():
                            if node.alpha[site] == "?":         
                                masked_llh = log(1-(1-phi)*p**nu) if (1-(1-phi)*p**nu)>0 else min_llh
                                node.L0[site] = node.L1[site] = masked_llh
                            else:    
                                node.L0[site] = nu*(-node.edge_length) + log(1-p) + log(q) + log(1-phi) if (1-p)*q*(1-phi)>0 else min_llh
                                node.L1[site] = nu*(-node.edge_length) + log(1-phi)
                        else:
                            C = node.children
                            l0 = l1 = 0
                            for c in C:
                                l0 += c.L0[site]
                                l1 += c.L1[site]
                            #L0 = exp(l0+(nu+1)*(-node.edge_length)) + exp(l1 + log(1-p)+log(q) + nu*(-node.edge_length)) + (1-p**nu)*int(node.alpha[site]=="?")   
                            L0 = exp(l0+(nu+1)*(-node.edge_length)) + q*(1-p)*exp(l1 + nu*(-node.edge_length)) + (1-p**nu)*int(node.alpha[site]=="?")   
                            L1 = exp(l1+nu*(-node.edge_length)) + (1-p**nu)*int(node.alpha[site]=="?")
                            node.L0[site] = min_llh if L0==0 else log(L0)
                            node.L1[site] = min_llh if L1==0 else log(L1)
                                
                        if node.is_root() or node.parent.alpha[site] == 'z':
                            llh[site] += node.L0[site]
                    else:
                        llh[site] += (-node.edge_length*(1+nu) + int(node.is_leaf())*log(1-phi))
        return sum(llh)         

    def ini_brlens(self):
        x = [random() * (self.dmax/2 - 2*self.dmin) + 2*self.dmin for i in range(self.num_edges)]        
        idx = 0
        for tree in self.trees:
            for node in tree.traverse_postorder():
                if node.edge_length is not None:
                    x[idx] = max(2*self.dmin,node.edge_length)
                idx += 1     
        return x    

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
        return [self.dmin]*self.num_edges,[self.dmax]*self.num_edges
        
    def get_bound(self,keep_feasible=False,fixed_phi=None,fixed_nu=None):
        br_lower,br_upper = self.bound_brlen()  
        phi_lower,phi_upper = self.bound_phi(fixed_phi=fixed_phi)
        nu_lower,nu_upper = self.bound_nu(fixed_nu=fixed_nu)
        bounds = optimize.Bounds(br_lower+[nu_lower,phi_lower],br_upper+[nu_upper,phi_upper],keep_feasible=keep_feasible)
        return bounds

    def x2brlen(self,x):
        i = 0
        for tree in self.trees:
            for node in tree.traverse_postorder():
                node.edge_length = x[i]
                i += 1

    def x2nu(self,x,fixed_nu=None):
        self.params.nu = x[self.num_edges] if fixed_nu is None else fixed_nu
    
    def x2phi(self,x,fixed_phi=None):
        self.params.phi = x[self.num_edges+1] if fixed_phi is None else fixed_phi

    def x2params(self,x,fixed_nu=None,fixed_phi=None):
        self.x2brlen(x)
        self.x2nu(x,fixed_nu=fixed_nu)
        self.x2phi(x,fixed_phi=fixed_phi)

    def __llh__(self):
        return self.lineage_llh()

    def negative_llh(self):
        self.az_partition()
        return -self.__llh__()

    def optimize(self,initials=20,fixed_phi=None,fixed_nu=None,fixed_brlen={},verbose=1,max_trials=100,random_seeds=None,ultra_constr=False):
    # random_seeds can either be a single number or a list of intergers where len(random_seeds) = initials
    # verbose level: 1 --> show all messages; 0 --> show minimal messages; -1 --> completely silent
    # fixed_brlen is a dictionary that maps a tuple (a,b) to a number. Each pair a, b is a tuple of two leaf nodes whose LCA define the node for the branch above it to be fixed.
        results = []
        all_failed = True
        all_trials = 0
        if random_seeds is None:
            rseeds = [int(random()*10000) for i in range(initials)]
        elif type(random_seeds) == int:
            if verbose >= 0:
                print("Global random seed: " + str(random_seeds))
            seed(a=random_seeds)
            rseeds = [int(random()*10000) for i in range(initials)]
        elif type(random_seeds) == list:
            if len(random_seeds) < initials:
                if verbose >= 0:
                    print("Fatal: the number of random seeds is smaller than the number of initials!")
                return None
            elif len(random_seeds) > initials:
                if verbose >= 0:
                    print("Warning: the number of random seeds is larger than the number of initials. Ignoring the last " + str(len(random_seeds)-initials) + " seeds")
            rseeds = random_seeds[:initials]    
        else:
            if verbose >= 0:
                print("Fatal: incorrect random_seeds type provided")        
            return None
        while all_failed and all_trials < max_trials:
            if verbose > 0:
                print("Optimization start with " + str(initials) + " initials")
            for rep in range(initials):
                randseed = rseeds[rep] + all_trials
                if verbose >= 0:
                    print("Initial point " + str(rep+1) + ". Random seed: " + str(randseed))
                if verbose >= 0:
                    if ultra_constr:
                        print("Numerical optimization started with ultrametric constraint")
                    else:      
                        print("Numerical optimization started without ultrametric constraint")
                # read in fixed_brlen and mark the tree nodes
                for tree in self.trees:
                    for node in tree.traverse_postorder():
                        node.mark_fixed=False        
                    fixed_nodes = find_LCAs(tree,list(fixed_brlen.keys()))        
                    for i,(a,b) in enumerate(fixed_brlen):
                        u = fixed_nodes[i]
                        u.edge_length = fixed_brlen[(a,b)]
                        u.mark_fixed = True
                nllh,status = self.optimize_one(randseed,fixed_phi=fixed_phi,fixed_nu=fixed_nu,verbose=verbose,ultra_constr=ultra_constr)
                
                if nllh is not None:
                    all_failed = False
                    if verbose >= 0:
                        print("Optimal point found for initial point " + str(rep+1))
                        #self.show_params()
                    # remove zero-length branches
                    processed_trees = []
                    for tree in self.trees:
                        tree_copy = read_tree_newick(tree.newick())
                        tree_copy.collapse_short_branches(self.dmin*0.01)
                        processed_trees.append(tree_copy.newick())
                    results.append((nllh,rep,deepcopy(self.params),processed_trees,status))
                elif verbose >= 0:
                    print("Fatal: failed to optimize using initial point " + str(rep+1))    
            all_trials += initials    
        if all_failed:
            if verbose >= 0:
                print("Fatal: Optimization failed on more than " + str(max_trials) + " retries. Aborting ...")
            return None
        else:    
            results.sort()
            best_nllh,_,best_params,best_trees,status = results[0]
            self.trees = []
            for tree in best_trees:
                self.trees.append(read_tree_newick(tree))
            self.params = best_params
            if verbose >= 0:
                print("Numerical optimization finished successfully")
            return results[0][0],status

    def optimize_one(self,randseed,fixed_phi=None,fixed_nu=None,verbose=1,ultra_constr=False):
        # optimize using a specific initial point identified by the input randseed
        # verbose level: 1 --> show all messages; 0 --> show minimal messages; -1 --> completely silent
        warnings.filterwarnings("ignore")
        def nllh(x): 
            self.x2params(x,fixed_nu=fixed_nu,fixed_phi=fixed_phi)            
            return -self.__llh__()
        
        seed(a=randseed)
        x0 = self.ini_all(fixed_phi=fixed_phi,fixed_nu=fixed_nu)
        self.az_partition()
        bounds = self.get_bound(fixed_phi=fixed_phi,fixed_nu=fixed_nu)
        constraints = []    

        A = []
        b = []
        idx = 0
        for tree in self.trees:
            for node in tree.traverse_postorder():
                if node.mark_fixed:
                    a = [0]*len(x0)
                    a[idx] = 1
                    A.append(a)
                    b.append(node.edge_length)
                idx += 1   
        if len(A) > 0:     
            constraints.append(optimize.LinearConstraint(csr_matrix(A),b,b,keep_feasible=False))
        if ultra_constr:
            M = self.ultrametric_constr()
            constraints.append(optimize.LinearConstraint(csr_matrix(M),[0]*len(M),[0]*len(M),keep_feasible=False))
        disp = (verbose > 0)
        out = optimize.minimize(nllh, x0, method="SLSQP", options={'disp':disp,'iprint':3,'maxiter':1000}, bounds=bounds,constraints=constraints)
        if out.success:
            self.x2params(out.x,fixed_phi=fixed_phi,fixed_nu=fixed_nu)
            params = self.params
            f = out.fun
        else:
            f,params = None,None
        status = "optimal" if out.success else out.message
        return f,status
