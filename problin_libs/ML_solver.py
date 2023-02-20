from treeswift import *
from math import log,exp,sqrt
from random import random, seed, choice
from scipy import optimize
import warnings
import numpy as np
from problin_libs import min_llh, eps, nni_conv_eps

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
    # beta_prior is a tuple (alpha,beta) that are the parameters of the beta prior of phi
    # if beta_prior is set to 'auto', run self.compute_beta_prior to estimate alpha and beta
    def __init__(self,charMtrx,Q,nwkTree,nu=eps,phi=eps):
    #def __init__(self,charMtrx,Q,nwkTree,nu=eps,phi=eps,beta_prior=(1,1)):
        self.charMtrx = charMtrx
        # normalize Q
        self.Q = []
        for Q_i in Q:
            s = sum([Q_i[x] for x in Q_i])
            Q_i_norm = {x:Q_i[x]/s for x in Q_i}
            self.Q.append(Q_i_norm)
        self.params = Params(nwkTree,nu=nu,phi=phi)
        self.numsites = len(self.charMtrx[next(iter(self.charMtrx.keys()))])
        self.num_edges = len(list(self.params.tree.traverse_postorder()))
        # TODO put in dmax and dmin here, remove everywhere else!!
        #self.dmin
        #self.dmax
    
    def compute_beta_prior(self):
        msa = self.charMtrx
        D = []
        alphabet = list(msa.keys())
        for i,x in enumerate(alphabet):
            for y in alphabet[i+1:]:
                A = 0
                B = 0
                for (a,b) in zip(msa[x],msa[y]):
                    A += (a != '?' and b != '?')
                    B += (a != '?' or b != '?')   
                D.append((B-A)/2/B)
        N = len(D)
        mu = sum(D)/N
        var = sum([(x-mu)**2 for x in D])/N
        self.alpha = (mu*(1-mu)/var-1)*mu
        self.beta = self.alpha*(1-mu)/mu

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

        #print("a.alpha:", a.alpha)
        #print("b.alpha:", b.alpha)
        #print("c.alpha:", c.alpha)
        if strat == "vanilla":
            return max(d_ab, d_ac)
        elif strat == "shouldchange":
            return max(d_ab - d_bc, d_ac - d_bc)

    def score_terminal_branch(self, u, strat):
        v = u.get_parent()
        gp = v.get_parent()
        uncle = [w for w in gp.child_nodes() if w is not v][0]
        sister = [w for w in v.child_nodes() if w is not u][0]
        
        d_cu = self.compare_tags(uncle.alpha, u.alpha)
        d_su = self.compare_tags(sister.alpha, u.alpha)
        #print("uncle.alpha:", uncle.alpha)
        #print("u.alpha:", u.alpha)
        #print("sister.alpha", sister.alpha)
        
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
        l2n = self.params.tree.label_to_node(selection="all")
        for nlabel in keybranches:
            node = l2n[nlabel]
            if not node.is_leaf():
                num_internal += 1
                kb.append(node)
        return num_internal, kb


    def score_branches(self, strategy="vanilla", keybranches=[]):
        if self.params.tree.num_nodes(internal=True, leaves=False) <= 2:
            print("Provided tree does not have enough internal branches to perform a nearest neighbor interchange operation.")
            return None

        # Purpose: Score all branches before returning one to consider nnis around
        self.az_partition(self.params)
        branches = []

        if keybranches != []:
            for node in keybranches:
                s = self.score_internal_branch(node, strategy)
                branches.append((node, s))
        else:
            for node in self.params.tree.traverse_postorder():
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

    def score_tree(self):
        self.az_partition(self.params)
        return self.lineage_llh(self.params)

    def apply_nni(self, u, verbose):
        # apply nni [DESTRUCTIVE FUNCTION! Changes tree inside this function.]
        v = u.get_parent()
        u_edges = [w for w in u.child_nodes()]
        v_edges = [w for w in v.child_nodes() if w is not u]
        nni_moves = []

        a, b = u_edges
        c = v_edges[0]
        d_ab = self.compare_tags(a.alpha, b.alpha)
        d_ac = self.compare_tags(a.alpha, c.alpha)
        d_bc = self.compare_tags(b.alpha, c.alpha)

        w = v_edges[0] 
        pre_llh = self.score_tree()
        #if verbose:
        #    print("pre_llh", pre_llh)
            #print("pre_llh", self.params.tree.newick(), pre_llh)

        # explore in order of importance
        if d_bc > d_ac:
            # move a out
            u_children = [a, b] 
        else:
            # move b out
            u_children = [b, a] 
        
        for u_child in u_children:

            u_child.set_parent(v)
            u.remove_child(u_child)
            v.add_child(u_child)

            w.set_parent(u)
            v.remove_child(w)
            u.add_child(w)
            
            new_llh = self.score_tree()
            #if verbose:
            #    print("new_llh", new_llh)
                #print("new_llh", self.params.tree.newick(), new_llh)

            if new_llh > pre_llh:
                # log likelihood improved
                return True
            elif new_llh == pre_llh:
                #if verbose:
                #    print("same log likelihood", new_llh)
                    #print("same log likelihood", self.params.tree.newick(), new_llh)
                return True
            else:
                # REVERSE IF LIKELIHOOD IS NOT BETTER
                #if verbose:
                #    print("reversing...")
                u_child.set_parent(u)
                v.remove_child(u_child)
                u.add_child(u_child)
                
                w.set_parent(v)
                u.remove_child(w)
                v.add_child(w)
                
                new_llh = self.score_tree()
                #if verbose:
                #    print("new_llh", new_llh)
                    #print("new_llh", self.params.tree.newick(), new_llh)
        #if verbose:
        #    print(new_llh)
            #print(new_llh, self.params.tree.newick())
        return False

    def single_nni(self, verbose, trynextbranch=True, strategy="vanilla", keybranches=[]):
        branches = self.score_branches(strategy, keybranches)
        took = False
        bidx = 0
        while not took:
            if verbose:
                print("Branch Attempt:", bidx)

            if len(branches) == 0:
                break
            elif strategy == "random":
                m = choice(branches)
                u, u_score = m
            else:
                m = max(branches, key=lambda item:item[1])
                u, u_score = m
            
            midx = branches.index(m)

            branches.pop(midx)
            took = self.apply_nni(u, verbose)
            bidx += 1
            if not trynextbranch:
                took = True 
        if verbose:
            print(bidx, " branch attempts.")
        llh = self.score_tree()
        return llh

    def tree_copy(self):
        tree = self.params.tree
        return tree.extract_subtree(tree.root)

    def num_internal_branches(self):
        nib = 0
        tree = self.params.tree
        for node in self.params.tree.traverse_postorder():
            if not node.is_leaf():
                nib += 1
        return float(nib) 
    def topology_search(self, maxiter=100, verbose=False, prefix="results_nni", trynextbranch=False, strategy="vanilla", keybranches=[], nreps=1, outdir="", conv=0.2):
        nib = self.num_internal_branches()
        t = round(0.2 * nib) + 1
        k = -int(log(conv)/log(t) * nib)
        if verbose:
            print("Running topology search for", k, "iterations.")
<<<<<<< HEAD
=======

>>>>>>> fceae3f (Fixed the unit tests, deleted scripts from problin_libs.)
        resolve_polytomies = False
        if verbose:
            print("Starting topology search.")
        if keybranches != []:
            resolve_polytomies = True
            if verbose:
                print("Doing topology search on polytomies.")
            nib, keybranches = self.resolve_keybranches(keybranches)
        elif verbose:
            print("Doing topology search on polytomy-resolved tree.")

        nni_replicates = dict()
        starting_tree = self.tree_copy()
        for i in range(nreps):

            topo_dict = {}
            seen = set()
            self.params.tree = starting_tree
            nni_iter = 0
            same = 0
            pre_llh = self.score_tree()
            
            while 1:
                if verbose:
                    print("NNI Iter:", nni_iter)
                opt_score = self.single_nni(verbose, trynextbranch=trynextbranch, strategy=strategy, keybranches=keybranches)
                
                tstr = self.params.tree.newick()
                topo_dict[nni_iter] = (tstr, opt_score)
                
                seen.add(tstr)
                new_llh = self.score_tree()

                if new_llh == pre_llh:
                    same += 1
                else:
                    same = 0
                
                if (new_llh - pre_llh < nni_conv_eps and tstr in seen and same > k) or nni_iter > maxiter:
                    break
                #if (nni_iter / nib) > conv:
                #    break
                pre_llh = new_llh
                nni_iter += 1

            nni_replicates[i] = (new_llh, topo_dict)
       
        #m = max(nni_replicates, key=lambda item:nni_replicates[item][0])
        #llh, topo_dict = nni_replicates[m]
        
        #if verbose:
        #    print("Recording the best result by ending llh.")
        if resolve_polytomies:
            out1 = outdir + "/" + prefix + "_topo_llh_resolvingpolytomies.txt"
            out2 = outdir + "/" + prefix + "_topo_progress_resolvingpolytomies.nwk"
        else:
            out1 = outdir + "/" + prefix + "_topo_llh.txt"
            out2 = outdir + "/" + prefix + "_topo_progress.nwk"

        with open(out1, "w+") as w:
            w.write("RandomRep\tnniIter\tNLLH\n")
            for rep in nni_replicates:
                llh, topo_dict = nni_replicates[rep]
                for nni_iter in topo_dict:
                    w.write(str(rep) + "\t" + str(nni_iter) + "\t" + str(-topo_dict[int(nni_iter)][1]) + "\n")
        with open(out2, "w+") as w:
            w.write("RandomRep\tnniIter\tTopology\n")
            for rep in nni_replicates:
                llh, topo_dict = nni_replicates[rep]
                for nni_iter in topo_dict:
                    w.write(str(rep) + "\t" + str(nni_iter) + "\t" + topo_dict[int(nni_iter)][0] + "\n") 
        #if resolve_polytomies:
        #    self.topology_search(maxiter, verbose, prefix, trynextbranch, strategy, [], nreps, outdir, conv)
    def az_partition(self,params):
    # Purpose: partition the tree into edge-distjoint alpha-clades and z-branches
    # Note: there is a different partition for each target-site
    # Output: annotate each node of the tree by node.alpha
        # z-branches are given tag 'z'
        # each of other branches is given a tag 
        # alpha where alpha is the alpha-tree it belongs to
        for node in params.tree.traverse_postorder():
            if node.is_leaf():
                node.alpha = [None]*self.numsites
                for site in range(self.numsites):
                    if self.charMtrx[node.label][site] == 0:
                        node.alpha[site] = 'z'
                    elif self.charMtrx[node.label][site] == -1:   
                        node.alpha[site] = '?' 
                    else:
                        node.alpha[site] = self.charMtrx[node.label][site]
                #node.alpha = [self.charMtrx[node.label][site] if self.charMtrx[node.label][site] != 0 else 'z' for site in range(self.numsites)]
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
    
    def lineage_llh(self,params):
        # assume az_partition has been performed so
        # each node has the attribute node.alpha
        phi = params.phi
        nu = params.nu
        llh = [0]*self.numsites
        for node in params.tree.traverse_postorder():
            p = exp(-node.edge_length)
            node.L0 = [0]*self.numsites # L0 and L1 are stored in log-scale
            node.L1 = [0]*self.numsites
            for site in range(self.numsites):    
                if node.alpha[site] != 'z':
                    q = self.Q[site][node.alpha[site]] if node.alpha[site] != "?" else 1.0
                    if node.is_leaf():
                        if node.alpha[site] == "?":         
                            masked_llh = log(1-(1-phi)*p**nu)
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
                    llh[site] += (-node.edge_length*(1+params.nu) + int(node.is_leaf())*log(1-phi))
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
        return self.lineage_llh(self.params) #+ (self.alpha-1)*log(self.params.phi) + (self.beta-1)*log(1-self.params.phi)


    def negative_llh(self):
        self.az_partition(self.params)
        return -self.__llh__()

    def show_params(self):
        print("tree: " + self.params.tree.newick())
        print("nu: " + str(self.params.nu))
        print("phi: " + str(self.params.phi))
        print("negative-llh: " + str(self.negative_llh()))
    
    def optimize(self,initials=20,fixed_phi=None,fixed_nu=None,verbose=True,max_trials=100,random_seeds=None):
    # random_seeds can either be a single number or a list of intergers where len(random_seeds) = initials
        results = []
        all_failed = True
        all_trials = 0
        if random_seeds is None:
            rseeds = [int(random()*10000) for i in range(initials)]
        elif type(random_seeds) == int:
            print("Global random seed: " + str(random_seeds))
            seed(a=random_seeds)
            rseeds = [int(random()*10000) for i in range(initials)]
        elif type(random_seeds) == list:
            if len(random_seeds) < initials:
                print("Fatal: the number of random seeds is smaller than the number of initials!")
                return None
            elif len(random_seeds) > initials:
                print("Warning: the number of random seeds is larger than the number of initials. Ignoring the last " + str(len(random_seeds)-initials) + " seeds")
            rseeds = random_seeds[:initials]    
        else:
            print("Fatal: incorrect random_seeds type provided")        
            return None
        while all_failed and all_trials < max_trials:
            if verbose:
                print("Starting with initials: ", initials)
            for rep in range(initials):
                randseed = rseeds[rep]
                print("Initial point " + str(rep+1) + ". Random seed: " + str(randseed))
                #print("Running EM with initial point " + str(rep+1))
                nllh,params = self.optimize_one(randseed,fixed_phi=fixed_phi,fixed_nu=fixed_nu,verbose=verbose)
                if nllh is not None:
                    all_failed = False
                    print("Optimal point found for initial point " + str(rep+1))
                    print("Optimal phi: " + str(params.phi))
                    print("Optimal nu: " + str(params.nu))
                    print("Optimal tree: " + params.tree.newick())
                    print("Optimal nllh: " + str(nllh))
                    results.append((nllh,params))
                else:
                    print("Fatal: failed to optimize using initial point " + str(rep+1))    
            all_trials += initials    
        results.sort()
        best_nllh,best_params = results[0]
        self.params = best_params
        return results[0][0]

    def optimize_one(self,randseed,fixed_phi=None,fixed_nu=None,verbose=True):
        # optimize using a specific initial point identified by the input randseed
        warnings.filterwarnings("ignore")
        def nllh(x): 
            self.x2params(x,fixed_nu=fixed_nu,fixed_phi=fixed_phi)            
            return -self.__llh__()
        
        seed(a=randseed)
        x0 = self.ini_all(fixed_phi=fixed_phi,fixed_nu=fixed_nu)
        self.az_partition(self.params)
        bounds = self.get_bound(fixed_phi=fixed_phi,fixed_nu=fixed_nu)
        out = optimize.minimize(nllh, x0, method="SLSQP", options={'disp':verbose,'iprint':3,'maxiter':1000}, bounds=bounds)
        if out.success:
            self.x2params(out.x,fixed_phi=fixed_phi,fixed_nu=fixed_nu)
            params = self.params
            f = out.fun
        else:
            f,params = None,None
        return f,params

    def optimize_old(self,initials=20,fixed_phi=None,fixed_nu=None,verbose=True,max_trials=100,alpha=1,beta=1):
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
                #nllh,params = self.optimize_one(randseed,fixed_phi=fixed_phi,fixed_nu=fixed_nu,verbose=verbose)
                if out.success:
                #if type(nllh) == int:
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
