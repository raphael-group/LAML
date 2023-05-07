from treeswift import *
from math import log,exp,sqrt, isclose
from random import random, seed, choice
from scipy import optimize
import warnings
import numpy as np
from problin_libs import min_llh, eps, nni_conv_eps
from problin_libs.Virtual_solver import Virtual_solver
from scipy.sparse import csr_matrix
from copy import deepcopy
from problin_libs.lca_lib import find_LCAs

class Params:
    def __init__(self,nu,phi):
        self.nu = nu
        self.phi = phi

class ML_solver(Virtual_solver):
    def __init__(self,treeTopo,data,prior,params={'nu':eps,'phi':eps}):
        charMtrx = data['charMtrx']
        Q = prior['Q']
        nu = params['nu']
        phi = params['phi']
        self.charMtrx = charMtrx
        self.tree = read_tree_newick(treeTopo)  
        self.tree.suppress_unifurcations()      
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
        self.num_edges = len(list(self.tree.traverse_postorder()))
        self.dmin = 0.005
        self.dmax = 10
        zerocount = sum([self.charMtrx[e].count(0) for e in self.charMtrx])
        totalcount = self.numsites * len(self.charMtrx)
        zeroprop = zerocount/totalcount
        self.dmax = -log(zeroprop) if zeroprop != 0 else 10

    def get_tree_newick(self):
        return self.tree.newick()

    def get_params(self):
        return {'phi':self.params.phi,'nu':self.params.nu}

    def ultrametric_constr(self):
        #N = self.num_edges + 2 # add 2: phi and nu
        N = len(self.ini_all())
        M = []
        idx = 0 
        for node in self.tree.traverse_postorder():
            if node.is_leaf():
                node.constraint = [0.]*N
            else:
                c1,c2 = node.children
                m = [x-y for (x,y) in zip(c1.constraint,c2.constraint)]
                M.append(m)
                node.constraint = c1.constraint
            node.constraint[idx] = 1
            idx += 1        
        return M
    
    def compare_tags(self, tags1, tags2):
        # Purpose: compute similarity score from az-partition tags
        total, same = 0.0, 0.0

        assert len(tags1) == len(tags2)
        for t1, t2 in zip(tags1, tags2):
            # consider locations where neither is ?
            if t1 == '?' or t2 == '?':
                continue
            else:
                total += 1
                if t1 == t2 and t1 != 'z':
                    # NOTE: maybe this isn't the best way to handle comparison with z
                    same += 1
        return same/total

    def similarity_score(self, a, b, c, strat):
        # Purpose: score the branch according to the maximum similarity 
        d_ab = self.compare_tags(a.alpha, b.alpha)
        d_ac = self.compare_tags(a.alpha, c.alpha)
        d_bc = self.compare_tags(b.alpha, c.alpha)

        if strat == "vanilla":
            return max(d_ab, d_ac)
        elif strat == "shouldchange":
            return max(d_ab - d_bc, d_ac - d_bc)
        else:
            return 1    

    def score_terminal_branch(self, u, strat):
        v = u.get_parent()
        gp = v.get_parent()
        uncle = [w for w in gp.child_nodes() if w is not v][0]
        sister = [w for w in v.child_nodes() if w is not u][0]
        
        d_cu = self.compare_tags(uncle.alpha, u.alpha)
        d_su = self.compare_tags(sister.alpha, u.alpha)
        
        if strat == "vanilla":
            return d_cu
        elif strat == "shouldchange":
            return d_cu - d_su

    def score_internal_branch(self, u, strat):
        v = u.get_parent()
        cladeA, cladeB = [w for w in u.child_nodes()]
        cladeC = [w for w in v.child_nodes() if w is not u][0]
        return self.similarity_score(cladeC, cladeB, cladeA, strat)

    def resolve_keybranches(self, keybranches):
        num_internal = 0
        kb = []
        l2n = self.tree.label_to_node(selection="all")
        for nlabel in keybranches:
            node = l2n[nlabel]
            if not node.is_leaf():
                num_internal += 1
                kb.append(node)
        return num_internal, kb

    def score_branches(self, strategy="vanilla", keybranches=[]):
        #if self.params.tree.num_nodes(internal=True, leaves=False) <= 2:
        if self.tree.num_nodes(internal=True, leaves=False) <= 2:
            print("Provided tree does not have enough internal branches to perform a nearest neighbor interchange operation.")
            return None

        # Purpose: Score all branches before returning one to consider nnis around
        self.az_partition()
        branches = []

        if keybranches != []:
            for node in keybranches:
                s = self.score_internal_branch(node, strategy)
                branches.append((node, s))
        else:
            #for node in self.params.tree.traverse_postorder():
            for node in self.tree.traverse_postorder():
                if node.is_root():
                    continue
                if not node.is_leaf():
                    # consider moving it inside the tree
                    if strategy == "random":
                        branches.append((node, 1))
                    else:
                        s = self.score_internal_branch(node, strategy)
                        branches.append((node, s))
        return branches 

    def score_tree(self,strategy={'ultra_constr':False,'fixed_phi':None,'fixed_nu':None,'fixed_brlen':{}}):
        ultra_constr = strategy['ultra_constr']
        fixed_phi = strategy['fixed_phi']
        fixed_nu = strategy['fixed_nu']
        fixed_brlen = strategy['fixed_brlen']
        nllh,status = self.optimize(initials=1,verbose=-1,ultra_constr=ultra_constr,fixed_phi=fixed_phi,fixed_nu=fixed_nu,fixed_brlen=fixed_brlen)
        score = None if nllh is None else -nllh
        if score is None:
            print("Fatal error: failed to score tree " + self.get_tree_newick() + ". Optimization status: " + status)
        return score,status

    def tree_copy(self):
        tree = self.tree
        return tree.extract_subtree(tree.root)

    def num_internal_branches(self):
        nib = 0
        #tree = self.params.tree
        #for node in self.params.tree.traverse_postorder():
        for node in self.tree.traverse_postorder():
            if not node.is_leaf():
                nib += 1
        return float(nib) 
    
    def az_partition(self):
    # Purpose: partition the tree into edge-distjoint alpha-clades and z-branches
    # Note: there is a different partition for each target-site
    # Output: annotate each node of the tree by node.alpha
        # z-branches are given tag 'z'
        # each of other branches is given a tag 
        # alpha where alpha is the alpha-tree it belongs to
        for node in self.tree.traverse_postorder():
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
        #for node in params.tree.traverse_postorder():
        for node in self.tree.traverse_postorder():
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
                            #if p == 1:
                            #    node.L0[site] = -float("inf")
                            #else:    
                            #    node.L0[site] = nu*(-node.edge_length) + log(1-p) + log(q) + log(1-phi)
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
                    #llh[site] += (-node.edge_length*(1+params.nu) + int(node.is_leaf())*log(1-phi))
                    llh[site] += (-node.edge_length*(1+nu) + int(node.is_leaf())*log(1-phi))
        return sum(llh)         

    def ini_brlens(self):
        x = [random() * (self.dmax/2 - 2*self.dmin) + 2*self.dmin for i in range(self.num_edges)]        
        for i,node in enumerate(self.tree.traverse_postorder()):
            #if node.mark_fixed:
            if node.edge_length is not None:
                x[i] = node.edge_length
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
        #for i, node in enumerate(self.params.tree.traverse_postorder()):
        for i, node in enumerate(self.tree.traverse_postorder()):
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
        return self.lineage_llh()

    def negative_llh(self):
        self.az_partition()
        return -self.__llh__()

    def show_params(self):
        print("tree: " + self.tree.newick())
        print("nu: " + str(self.params.nu))
        print("phi: " + str(self.params.phi))
        print("negative-llh: " + str(self.negative_llh()))
    
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
                randseed = rseeds[rep]
                if verbose >= 0:
                    print("Initial point " + str(rep+1) + ". Random seed: " + str(randseed))
                if verbose >= 0:
                    if ultra_constr:
                        print("Solving ML with ultrametric constraint")
                    else:      
                        print("Solving ML without ultrametric constraint")
                for node in self.tree.traverse_postorder():
                    node.mark_fixed=False        
                fixed_nodes = find_LCAs(self.tree,list(fixed_brlen.keys()))        
                for i,(a,b) in enumerate(fixed_brlen):
                    u = fixed_nodes[i]
                    u.edge_length = fixed_brlen[(a,b)]
                    u.mark_fixed = True
                nllh,status = self.optimize_one(randseed,fixed_phi=fixed_phi,fixed_nu=fixed_nu,verbose=verbose,ultra_constr=ultra_constr)
                
                if nllh is not None:
                    all_failed = False
                    if verbose >= 0:
                        print("Optimal point found for initial point " + str(rep+1))
                        self.show_params()
                    # remove zero-length branches
                    tree_copy = read_tree_newick(self.tree.newick())
                    tree_copy.collapse_short_branches(self.dmin*0.01)
                    results.append((nllh,rep,deepcopy(self.params),tree_copy.newick(),status))
                elif verbose >= 0:
                    print("Fatal: failed to optimize using initial point " + str(rep+1))    
            all_trials += initials    
        if all_failed:
            if verbose >= 0:
                print("Fatal: Optimization failed on more than " + str(max_trials) + " retries")
            return None
        else:    
            results.sort()
            best_nllh,_,best_params,best_tree,status = results[0]
            self.tree = read_tree_newick(best_tree)
            self.params = best_params
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
        for node in self.tree.traverse_postorder():
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
