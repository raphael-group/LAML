from treeswift import *
from math import log,exp,sqrt, isclose
from random import random, seed, choice
from scipy import optimize
import warnings
import numpy as np
from laml_libs import *
from .Virtual_solver import Virtual_solver
from .Mstep_solver import *
from scipy.sparse import csr_matrix
from copy import deepcopy
from laml_libs.Utils.lca_lib import find_LCAs
from .Param import Param
from .Alphabet import Alphabet
from .AlleleTable import AlleleTable
from .utils import * 
import time
import cvxpy as cp

class Base_model(Virtual_solver):
    def __init__(self,treeList,data,prior,params):
    # `params` is an instance of the Param class
    # `data` is a dictionary of multiple data modules; it MUST have 'DLT_data'
    # this base class only uses `DLT_data` of `data`, but a derived class can use more attributes for joint likelihood computation
    # `prior` contains hyper-parameters and information about model variations (e.g. silencing mechanism)
    # this base class uses `Q` and `silence_mechanism` of `prior`, but a derived class can use more depending on the model
        self.data = data
        self.num_cassettes = data['DLT_data'].K
        self.site_per_cassette = data['DLT_data'].J
        self.params = params
        self.trees = []
        self.num_edges = 0
        for tree in treeList:
            tree_obj = read_tree_newick(tree)
            #tree_obj.suppress_unifurcations()
            if len(tree_obj.root.children) > 1:
                new_root = Node()
                new_root.add_child(tree_obj.root)
                tree_obj.root = new_root
            self.num_edges += len(list(tree_obj.traverse_postorder()))-1
            for node in tree_obj.traverse_postorder():
                node.mark_recompute = True
            self.trees.append(tree_obj)

        # get silence_mechanism; default to 'convolve'
        self.silence_mechanism = 'convolve'
        if 'silence_mechanism' in prior:
            self.silence_mechanism = prior['silence_mechanism'] # should be one of {'convolve','separated'}

        # normalize Q
        Q = prior['Q']
        self.Q = []
        for Q_k in Q:
            Q_k_norm = []
            for Q_kj in Q_k:
                s = sum([Q_kj[x] for x in Q_kj])
                Q_kj_norm = {x:Q_kj[x]/s for x in Q_kj} 
                Q_k_norm.append(Q_kj_norm)                   
            self.Q.append(Q_k_norm)        
            
        # set dmin and dmax 
        self.dmin = DEFAULT_dmin
        self.dmax = DEFAULT_dmax

    def get_tree_newick(self):
    # override Virtual solver
        return [tree.newick() for tree in self.trees]

    def get_params(self):
    # override Virtual solver
        return self.params.get_name2value_dict()

    def score_tree(self,strategy={'ultra_constr':False,'fixed_params':None,'fixed_brlen':None,'compute_cache':None}):
    # override Virtual solver
        ultra_constr = strategy['ultra_constr']
        fixed_params = strategy['fixed_params']
        fixed_brlen = strategy['fixed_brlen']
        compute_cache = strategy['compute_cache']
        nllh,status = self.optimize('EM',initials=1,verbose=-1,ultra_constr=ultra_constr,fixed_params=fixed_params,
                                        fixed_brlen=fixed_brlen,compute_cache=compute_cache)
        score = None if nllh is None else -nllh
        return score,status
    
    def get_compute_cache(self):
    # override Virtual solver
        # run Estep with all node recomputed to make sure all values are up-to-date 
        for tree in self.trees:
            for node in tree.traverse_postorder():
                node.mark_recompute = True
        self.Estep()
        # pull values to cache
        cache_attrs = ['in_llh','out_llh','in_llh_edge','log_node_posterior','log_edge_posterior']
        #cache_attrs = ['log_node_posterior']
        full_cache = []
        for tree in self.trees:
            this_cache = {}
            for node in tree.traverse_postorder():
                if node.is_root():
                    continue
                if node.is_leaf():
                    node.anchors = (node.label,node.label)
                else:
                    C = node.child_nodes()
                    a = C[0].anchors[0]
                    b = C[-1].anchors[0]
                    node.anchors = (a,b)
                this_cache[node.anchors] = {}    
                node_dict = node.__dict__
                for attr in cache_attrs:  
                    if attr in node_dict:
                        this_cache[node.anchors][attr] = deepcopy(node_dict[attr])
                        #for x in node_dict[attr]:
                        #    this_cache[node.anchors][attr].append(deepcopy(x))
            full_cache.append(this_cache)
        #for x in 'a','b','c','d':
        #    print('send',x,full_cache[0][(x,x)]['log_node_posterior'])
        return full_cache        
    
    def ultrametric_constr(self,local_brlen_opt=True):
        N = self.num_edges#-self.num_polytomy_mark
        if local_brlen_opt:
            for tree in self.trees:
                N -= len([node for node in tree.traverse_postorder() if node.mark_fixed])
        constrs = {}        
        idx = 0
        for tree in self.trees:
            for node in tree.traverse_postorder():
                if node.is_leaf():
                    node.constraint = [0.]*N
                    node.constant = node.edge_length if node.mark_fixed else 0
                elif not node.is_root(): #####HACKING#####
                    if len(node.children) == 2: #####HACKING#####
                        c1,c2 = node.children
                        m = tuple(x-y for (x,y) in zip(c1.constraint,c2.constraint))
                        m_compl = tuple(-x for x in m)
                        c = c2.constant-c1.constant
                        if sum([x!=0 for x in m]) > 0 and not (m in constrs or m_compl in constrs):
                            constrs[m] = c
                    else: #####HACKING#####
                        c1 = node.children[0]
                    node.constraint = c1.constraint
                    if node.mark_fixed:
                        node.constant = c1.constant + node.edge_length
                    else:
                        node.constant = c1.constant
                #if not node.polytomy_mark and not (node.mark_fixed and local_brlen_opt):    
                if not node.is_root() and not (node.mark_fixed and local_brlen_opt): #####HACKING#####
                    node.constraint[idx] = 1
                    idx += 1
        for tree in self.trees[1:]:
            m = tuple(x-y for (x,y) in zip(self.trees[0].root.constraint,tree.root.constraint))
            c = tree.root.constant-self.trees[0].root.constant
            constrs[m] = c
        M = []
        b = []
        for m in constrs:
            M.append(list(m))
            b.append(constrs[m])
        return M,b

    def ini_brlens(self):
        x = [random() * (self.dmax/2 - 2*self.dmin) + 2*self.dmin for i in range(self.num_edges)]        
        idx = 0
        for tree in self.trees:
            for node in tree.traverse_postorder():
                if node.edge_length is not None:
                    x[idx] = node.edge_length if node.mark_fixed else max(2*self.dmin,node.edge_length)
                idx += 1 
        return x

    def ini_one_param(self,pname):
        # here we set a simple initialization scheme for all params
        # in most cases this method should be overrided in a class that inhirits from this base class
        lower_bound,upper_bound = self.params.get_bound(pname)
        return random()*(upper_bound-lower_bound)+lower_bound

    def ini_all_params(self,fixed_params={}):
        # initialize branch lengths
        x = self.ini_brlens()
        # initialize the params specified in self.params
        for p in self.params.get_names():
            x += [fixed_params[p]] if p in fixed_params else [self.ini_one_param(p)]
        return x    

    def bound_brlen(self):        
        return [self.dmin]*self.num_edges,[self.dmax]*self.num_edges
    
    def get_bound(self,keep_feasible=False,fixed_params={}):
        lower_bounds,upper_bounds = self.bound_brlen()  
        for i,(l,u) in enumerate(zip(self.params.lower_bounds,self.params.upper_bounds)):
            pname = self.params.names[i]
            if pname in fixed_params:
                l = fixed_params[pname]-eps
                u = fixed_params[pname]+eps
            lower_bounds.append(l)
            upper_bounds.append(u)                
        bounds = optimize.Bounds(lower_bounds,upper_bounds,keep_feasible=keep_feasible)
        return bounds

    def x2params(self,x,fixed_params={}):
        # x2brlens
        i = 0
        for tree in self.trees:
            for node in tree.traverse_postorder():
                if not node.is_root(): #and not node.mark_fixed:
                    node.edge_length = x[i]
                    i += 1
        # other x2params                
        for j,p in enumerate(self.params.names):
            self.params.values[j] = fixed_params[p] if p in fixed_params else x[i]
            i += 1
    
    def Psi(self,c_node,k,j,x,y):
        # Layer 1 transition probabilities
        # This is a placeholder (i.e. non-informative model) for this function in the base class
        # MUST be overrided in any derived class!
        return 1

    def log_Psi_cassette(self,c_node,k,x,y):
        # compute the log-transformed transition probability of cassette k  
        # return None if the probability is 0
        log_trans_p = 0
        for j in range(self.data['DLT_data'].J):
            p = self.Psi(c_node,k,j,x[j],y[j])
            if p > 0:
                log_trans_p += log(p) #####*****#####
            else:    
                log_trans_p = None
                break
        return log_trans_p

    def Gamma(self,k,x,c,node_name):
        # Layer 2 transition probabilities
        # A placeholder (i.e. non-informative model) for this function in the base class
        # MUST be overrided in any derived class!
        return 1

    def llh_DLT_data(self):
        # compute the log-likelihood of the dynamic lineage tracing (DLT) data (e.g. character matrix, allele table)
        K = self.data['DLT_data'].K
        DLT_data = self.data['DLT_data']
        root_state = tuple([0]*self.data['DLT_data'].J)
        for tree in self.trees:
            for node in tree.traverse_postorder():
                if node.mark_recompute:
                    node.in_llh = [{} for _ in range(K)] # a list of dictionaries
        
        total_llh = 0
        for k in range(K):
            allele_list = self.data['DLT_data'].alphabet.get_cassette_alphabet(k)
            for tree in self.trees:
                for node in tree.traverse_postorder():
                    if not node.mark_recompute:
                        continue
                    if node.is_leaf():
                        c = DLT_data.get(node.label,k)
                        for x in allele_list: 
                            trans_p = self.Gamma(k,x,c,node.label) # transition probability
                            if trans_p>0:
                                node.in_llh[k][x] = log(trans_p)
                    else:   
                        for x in allele_list:
                            llh = 0
                            for c_node in node.children:
                                llh_list = []                                    
                                for y in c_node.in_llh[k]:
                                    log_trans_p = self.log_Psi_cassette(c_node,k,x,y)
                                    if log_trans_p is not None:
                                        llh_list.append(log_trans_p + c_node.in_llh[k][y])
                                if llh_list: # if the list is not empty
                                    llh += log_sum_exp(llh_list)
                                else:
                                    llh = None
                                    break   
                            if llh is not None:
                                node.in_llh[k][x] = llh
                total_llh += tree.root.in_llh[k][root_state]
        return total_llh

    def llh_DLT_data_edge(self,node,k,x):
        # assume in_llh has been computed for every node
        llh_list = []
        for y in node.in_llh[k]:
            log_trans_p = self.log_Psi_cassette(node,k,x,y)
            if log_trans_p is not None: 
                llh_list.append(node.in_llh[k][y] + log_trans_p)        
        return log_sum_exp(llh_list) if llh_list else None       

    def negative_llh(self):
        # compute the negative log-likelihood of all data modules
        # in the base class, only the allele table is included
        # if a derived class has data modules other than the allele table, 
        # this function MUST be overrided to compute the joint-llh 
        # of all available data modules
        return -self.llh_DLT_data()

    def optimize(self,solver,initials=20,fixed_brlen=None,compute_cache=None,fixed_params={},verbose=1,max_trials=100,random_seeds=None,ultra_constr=False,**solver_opts):
    # solver can either be "Scipy" or "EM"
    # random_seeds can either be a single number or a list of intergers where len(random_seeds) = initials
    # verbose level: 1 --> show all messages; 0 --> show minimal messages; -1 --> completely silent
    # fixed_brlen can be one of the followings:
        # 1. None: don't fix any branch length
        # 2. (str) "All": fix all branch lengths to the current value. Assume without checking that solver.trees have (valid) branch lengths
        # 3. a list of t dictionaries where t is the number of trees in self.trees. Each dictionay fixed_brlen[t] maps a tuple (a,b) to a number. Each pair a, b is a tuple of two leaf nodes whose LCA define the node for the branch above it to be fixed
    # fixed_params is a dictionary mapping a param's name to the value we wish to fix it to
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

                # read in compute_cache
                if compute_cache is not None:
                    for t,tree in enumerate(self.trees):
                        cache_nodes = find_LCAs(tree,list(compute_cache[t].keys()))
                        for i,(a,b) in enumerate(compute_cache[t]):
                            u = cache_nodes[i]
                            for attr in compute_cache[t][(a,b)]:
                                setattr(u, attr, compute_cache[t][(a,b)][attr])
                
                # read in fixed_brlen and mark the tree nodes
                for t,tree in enumerate(self.trees):
                    # initialization
                    for node in tree.traverse_postorder():
                        node.mark_fixed = (fixed_brlen == 'All' and not node.is_root())
                    if fixed_brlen == 'All': # the above piece of code should have set all mark_fixed flags to True
                        continue    
                    if fixed_brlen is not None:
                        fixed_nodes = find_LCAs(tree,list(fixed_brlen[t].keys()))        
                        for i,(a,b) in enumerate(fixed_brlen[t]):
                            u = fixed_nodes[i]
                            u.edge_length = fixed_brlen[t][(a,b)]
                            u.mark_fixed = True
                    # mark the nodes that will need to be recomputed
                    for node in tree.traverse_postorder():
                        node.mark_recompute = (not node.mark_fixed) or (compute_cache is None)
                        for c_node in node.children:
                            node.mark_recompute = node.mark_recompute or c_node.mark_recompute
               
                # NOTE: the solvers will use the flags "mark_fixed" and "mark_recompute" decorated on the tree nodes
                # Be careful: don't modify these flags while doing the computation/optimization inside the solver!
                if solver == 'Scipy':
                    scipy_options = solver_opts if solver_opts else DEFAULT_scipy_options
                    scipy_options['disp'] = (verbose>0)
                    nllh,status = self.scipy_optimization(randseed,fixed_params=fixed_params,ultra_constr=ultra_constr,scipy_options=scipy_options)
                else:
                    EM_options = {}
                    for opts in DEFAULT_EM_options:
                        EM_options[opts] = DEFAULT_EM_options[opts]
                    for opts in solver_opts:
                        EM_options[opts] = solver_opts[opts]    
                    nllh,status = self.EM_optimization(randseed,verbose=verbose,fixed_params=fixed_params,ultra_constr=ultra_constr,EM_options=EM_options)
                
                if nllh is not None:
                    all_failed = False
                    if verbose >= 0:
                        print("Optimal point found for initial point " + str(rep+1))
                    # remove zero-length branches
                    processed_trees = []
                    for tree in self.trees:
                        tree_copy = read_tree_newick(tree.newick())
                        tree_copy.collapse_short_branches(self.dmin*0.01) ##### HACKING: using a hard-code for now #####
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

    def Estep_in_llh(self):
        self.llh_DLT_data()

    def Estep_out_llh(self):
        K = self.data['DLT_data'].K
        DLT_data = self.data['DLT_data']
        root_state = tuple([0]*self.data['DLT_data'].J)
        
        for tree in self.trees:
            for node in tree.traverse_preorder():
                if node.mark_recompute:
                    node.out_llh = [{} for _ in range(K)] # a list of dictionaries
                    node.in_llh_edge = [{} for _ in range(K)] # a list of dictionaries

        for k in range(K):
            allele_list = self.data['DLT_data'].alphabet.get_cassette_alphabet(k)
            for tree in self.trees:
                for node in tree.traverse_preorder():
                    if not node.mark_recompute:
                        continue
                    if node.is_root():
                        node.out_llh[k][root_state] = 0
                    elif node.parent.is_root():
                        for x in allele_list:
                            log_trans_p = self.log_Psi_cassette(node,k,root_state,x) 
                            if log_trans_p is not None:
                                node.out_llh[k][x] = log_trans_p
                    else:
                        for x in allele_list:
                            llh_list = []                                    
                            par_state_lists = []
                            for j,x_j in enumerate(x):
                                if x_j == -1:
                                    L = self.data['DLT_data'].alphabet.get_site_alphabet(k,j)
                                else:
                                    L = set([0,x_j])
                                par_state_lists.append(L)
                            candidate_par_states = join_lists(par_state_lists)
                            for y in candidate_par_states:
                            #for y in node.parent.out_llh[k]:
                                if y not in node.parent.out_llh[k]:
                                    continue
                                log_trans_p = self.log_Psi_cassette(node,k,y,x)
                                if log_trans_p is None:
                                    continue
                                par_out_llh = node.parent.out_llh[k][y]
                                sum_sisters_in_llh = 0
                                for w in node.parent.children:
                                    if w is node:
                                        continue
                                    if y not in w.in_llh_edge[k]:
                                        curr_llh = self.llh_DLT_data_edge(w,k,y)    
                                        w.in_llh_edge[k][y] = curr_llh
                                    else:
                                        curr_llh = w.in_llh_edge[k][y]    
                                    #curr_llh = self.llh_DLT_data_edge(w,k,y)
                                    if curr_llh is not None:
                                        sum_sisters_in_llh += curr_llh
                                    else:
                                        sum_sisters_in_llh = None
                                        break
                                if sum_sisters_in_llh is not None:
                                    llh_list.append(log_trans_p+par_out_llh+sum_sisters_in_llh)
                            if llh_list: # if the list is not empty
                                node.out_llh[k][x] = log_sum_exp(llh_list)

    def Estep_posterior(self):
    # compute the log-transformed of all posterior probabilities
        K = self.data['DLT_data'].K
        DLT_data = self.data['DLT_data']
        root_state = tuple([0]*self.data['DLT_data'].J)
        
        for tree in self.trees:
            for node in tree.traverse_preorder():
                if node.mark_recompute:
                    node.log_node_posterior = [{} for _ in range(K)] # a list of dictionaries
                    node.log_edge_posterior = [{} for _ in range(K)] # a list of dictionaries

        for k in range(K):
            allele_list = self.data['DLT_data'].alphabet.get_cassette_alphabet(k)
            for tree in self.trees:
                for node in tree.traverse_preorder():
                    if not node.mark_recompute:
                        continue
                    if node.is_root():
                        node.log_node_posterior[k][root_state] = 0
                    else:
                        for x in allele_list:
                            in_llh = node.in_llh[k][x] if x in node.in_llh[k] else None
                            out_llh = node.out_llh[k][x] if x in node.out_llh[k] else None
                            total_llh = tree.root.in_llh[k][root_state]
                            if (in_llh is not None) and (out_llh is not None):
                                node.log_node_posterior[k][x] = in_llh + out_llh - total_llh
                        for x in node.parent.log_node_posterior[k]:
                            for y in node.in_llh[k]:
                                log_p_trans = self.log_Psi_cassette(node,k,x,y)
                                if log_p_trans is None:
                                    continue
                                log_u_post = node.parent.log_node_posterior[k][x]
                                log_v_in = node.in_llh[k][y]
                                log_llh_edge = self.llh_DLT_data_edge(node,k,x)
                                node.log_edge_posterior[k][(x,y)] = log_u_post + log_v_in + log_p_trans - log_llh_edge

    def Estep(self,run_in_llh=True):
        if run_in_llh:
            #start = time.time()
            self.Estep_in_llh()
            #end = time.time()
            #print("Estep in-llh:",end-start)
        #start = time.time()
        self.Estep_out_llh()
        #end = time.time()
        #print("Estep out-llh:",end-start)
        
        #start = time.time()
        self.Estep_posterior()
        #end = time.time()
        #print("Estep posterior:",end-start)

    def set_closed_form_optimal(self,fixed_params={},verbose=1):
        # For every param that has a closed-form M-step optimal, 
        # if that param is in fixed_params, set it to the specified fixed value;
        # otherwise, compute the closed-form solution and set that to the param's value
        # IMPORTANT: this is a placeholder for this function in the base class
        # MUST be overrided in any derived class!
        return False # this function should never be called, so this line should never be reached!

    def Mstep(self,fixed_params={},verbose=1,local_brlen_opt=True,ultra_constr_cache=None):
        # assume that Estep has been performed so that all nodes have S0-S4 attributes
        # output: optimize all parameters: branch lengths, phi, and nu
        # verbose level: 1 --> show all messages; 0 --> show minimal messages; -1 --> completely silent        
        
        # setup and call a Mstep_solver
        selected_solver = Mstep_PMMconv if self.silence_mechanism == "convolve" else Mstep_PMMsep
        my_Mstep_solver = selected_solver(self,ultra_constr_cache=ultra_constr_cache,local_brlen_opt=local_brlen_opt,nIters=1)
        d_star,status_d,nu_star,status_nu = my_Mstep_solver.solve(fixed_params=fixed_params,verbose=verbose)

        # place the optimal value back to params
        if status_nu in ["optimal","UNKNOWN"]:
            self.params.set_value('nu',nu_star)
        if status_d in ["optimal","UNKNOWN","fixed"]:
            i = 0
            for tree in self.trees:
                for node in tree.traverse_postorder():
                    if not node.is_root() and not (node.mark_fixed and local_brlen_opt):   
                        node.edge_length = d_star[i]
                        i += 1
        success = (status_d in ["fixed","optimal","UNKNOWN"] and status_nu in ["optimal","UNKNOWN"])
        if success:
            status = "optimal"
        else:
            status = ""
            if status_d != "optimal":
                status += "failed_d"
            if status_nu != "optimal":
                status = ",failed_nu"
        return success, status

    def Mstep_old(self,fixed_params={},verbose=1,local_brlen_opt=True,ultra_constr_cache=None,eps_nu=1e-5,eps_s=1e-6):
        # assume that Estep has been performed so that all nodes have S0-S4 attributes
        # output: optimize all parameters: branch lengths, phi, and nu
        # verbose level: 1 --> show all messages; 0 --> show minimal messages; -1 --> completely silent        
       
        # optimize the params that have a closed-form solution
        self.set_closed_form_optimal(fixed_params=fixed_params,verbose=verbose)

        # optimize nu and all branch lengths
        K = self.data['DLT_data'].K
        N = self.num_edges
        if local_brlen_opt:
            for tree in self.trees:
                N -= len([node for node in tree.traverse_postorder() if node.mark_fixed])
        C_z2z = np.zeros(N)
        C_z2a = np.zeros(N)
        C_a2a = np.zeros(N)
        C_za2s = np.zeros(N)
        
        d_ini = np.zeros(N)
        i = 0
        for tree in self.trees:
            for v in tree.traverse_postorder():
                if not v.is_root() and not (v.mark_fixed and local_brlen_opt):   
                    for k in range(K):    
                        for x,y in v.log_edge_posterior[k]:
                            w = exp(v.log_edge_posterior[k][(x,y)])
                            for x_j,y_j in zip(x,y):
                                C_z2z[i] += w*(x_j == 0 and y_j == 0)
                                C_z2a[i] += w*(x_j == 0 and y_j != 0 and y_j != -1)
                                C_a2a[i] += w*(x_j == y_j and x_j != 0 and x_j != -1)
                                C_za2s[i] += w*(y_j == -1 and x_j != -1)
                    d_ini[i] = v.edge_length
                    i += 1
        def __optimize_brlen(nu,verbose=False): # nu is a single number
            var_d = cp.Variable(N,nonneg=True) # the branch length variables
            SA = -(nu+1)*C_z2z.T @ var_d
            SB = -nu*(C_z2a+C_a2a).T @ var_d
            SC = (C_z2a).T @ cp.log(1-cp.exp(-var_d))
            SD = (C_za2s).T @ cp.log(1-cp.exp(-nu*var_d)) if sum(C_za2s) > 0 and nu > eps_nu else 0
            
            objective = cp.Maximize(SA+SB+SC+SD)
            constraints = [np.zeros(N)+self.dmin <= var_d, var_d <= np.zeros(N)+self.dmax]             
            if ultra_constr_cache is not None:
                M,b = ultra_constr_cache
                constraints += [np.array(M) @ var_d == np.array(b)]
            prob = cp.Problem(objective,constraints)
            prob.solve(verbose=verbose,solver=cp.MOSEK)
            return var_d.value,prob.status
       
        def __optimize_nu(d,verbose=False): # d is a vector of all branch lengths
            var_nu = cp.Variable(1,nonneg=True) # the nu variable
            SA = -(var_nu+1)*(C_z2z).T @ d
            SB = -var_nu*(C_z2a+C_a2a).T @ d
            SD = C_za2s.T @ cp.log(1-cp.exp(-var_nu*d)) if sum(C_za2s) > 0 else 0
            
            objective = cp.Maximize(SA+SB+SD)
            prob = cp.Problem(objective)
            prob.solve(verbose=verbose,solver=cp.MOSEK)
            return var_nu.value[0],prob.status

        nIters = 1
        nu_star = self.params.get_value('nu')
        for r in range(nIters):
            if verbose > 0:
                print("Optimizing branch lengths. Current nu:" + str(nu_star))
            try:
                d_star,status_d = __optimize_brlen(nu_star,verbose=False)
            except:
                d_star = d_ini
                status_d = "failure"
            if status_d == "infeasible": # should only happen with local EM 
                return False,"d_infeasible"
            if 'nu' in fixed_params:
                if verbose > 0:
                    print("Fixed nu to " + str(fixed_params['nu']))
                nu_star = fixed_params['nu']
                status_nu = "optimal"    
            else:    
                if verbose > 0:
                    print("Optimizing nu")
                try:
                    nu_star,status_nu = __optimize_nu(d_star,verbose=False)                 
                except:
                    status_nu = "failure"
        
        # place the optimal value back to params
        self.params.set_value('nu',nu_star)
        i = 0
        for tree in self.trees:
            for node in tree.traverse_postorder():
                if not node.is_root() and not (node.mark_fixed and local_brlen_opt):   
                    node.edge_length = d_star[i]
                    i += 1
        success = (status_d == "optimal" or status_d == "UNKNOWN") and (status_nu == "optimal" or status_nu == "UNKNOWN")
        if success:
            status = "optimal"
        else:
            status = ""
            if status_d != "optimal":
                status += "failed_d"
            if status_nu != "optimal":
                status = ",failed_nu"
        return success, status
    
    def EM_optimization(self,randseed,verbose=1,fixed_params={},ultra_constr=True,EM_options=DEFAULT_EM_options):
        # IMPORTANT: this EM algorithm ONLY optimizes the likelihood of the allele table !
        # optimize the following parameters: branch lengths and nu
        # caution: this function will modify params and branch lengths in place!
        # verbose level: 1 --> show all messages; 0 --> show minimal messages; -1 --> completely silent
        seed(a=randseed)
        maxIter = EM_options['max_iter']
        x0 = self.ini_all_params(fixed_params=fixed_params)
        self.x2params(x0,fixed_params=fixed_params)
        pre_llh = self.llh_DLT_data()
        if verbose >= 0:
            print("Initial parameter values: " + self.params.show_values())
            print("Initial nllh: " + str(-pre_llh))
        em_iter = 1
        converged = False
        if ultra_constr:
            ultra_constr_cache = self.ultrametric_constr(local_brlen_opt=True) #####*****##### 
        else:
            ultra_constr_cache = None        
        while em_iter <= maxIter:
            if verbose > 0:
                print("Starting EM iter: " + str(em_iter))
                print("Estep")
            estep_start = time.time()
            self.Estep(run_in_llh=False) #####*****#####
            estep_end = time.time()
            if verbose > 0:
                print(f"Estep runtime (s): {estep_end - estep_start}")
                print("Mstep")
            mstep_start = time.time()
            m_success,status=self.Mstep(fixed_params=fixed_params,verbose=verbose,local_brlen_opt=True,ultra_constr_cache=ultra_constr_cache)
            mstep_end = time.time()
            if verbose > 0:
                print(f"Mstep runtime (s): {mstep_end - mstep_start}")
            if not m_success:
                if status == "d_infeasible": # should only happen with local EM
                    if verbose >= 0:
                        print("Warning: EM failed to optimize parameters in one Mstep due to infeasible constraints") 
                    return -pre_llh,status
                elif verbose >= 0:    
                    print("Warning: EM failed to optimize parameters in one Mstep.")                
            curr_llh = self.llh_DLT_data()
            if verbose > 0:
                print("Finished EM iter: " + str(em_iter) + ". Current nllh: " + str(-curr_llh))
            if abs((curr_llh - pre_llh)/pre_llh) < DEFAULT_conv_eps:
                converged = True
                break
            pre_llh = curr_llh
            em_iter += 1
        if not converged and verbose >= 0:
            print("Warning: exceeded maximum number of EM iterations (" + str(maxIter) + " iters)!")
        return -curr_llh,status    

    def scipy_optimization(self,randseed,fixed_params={},ultra_constr=True,scipy_options=DEFAULT_scipy_options):
        # optimize using a specific initial point identified by the input randseed
        warnings.filterwarnings("ignore")
        def nllh(x): 
            self.x2params(x,fixed_params=fixed_params)            
            return self.negative_llh()
        
        seed(a=randseed)
        x0 = self.ini_all_params(fixed_params=fixed_params)
        bounds = self.get_bound(fixed_params=fixed_params)
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
            M,b = self.ultrametric_constr(local_brlen_opt=False) # temporary solution: this version only works with local_brlen_opt=False
            #constraints.append(optimize.LinearConstraint(csr_matrix(M),[0]*len(M),[0]*len(M),keep_feasible=False))
            constraints.append(optimize.LinearConstraint(csr_matrix(M),b,b,keep_feasible=False))
        out = optimize.minimize(nllh, x0, method="SLSQP", options=scipy_options, bounds=bounds,constraints=constraints)
        if out.success:
            self.x2params(out.x,fixed_params=fixed_params)
            params = self.params
            f = out.fun
        else:
            f,params = None,None
        status = "optimal" if out.success else out.message
        return f,status
