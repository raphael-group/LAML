from treeswift import *
from math import log,exp,sqrt, isclose
from random import random, seed, choice
from scipy import optimize
import warnings
import numpy as np
from laml_libs import *
from .Virtual_solver import Virtual_solver
from scipy.sparse import csr_matrix
from copy import deepcopy
from laml_libs.Utils.lca_lib import find_LCAs

def log_sum_exp(numlist):
    # using log-trick to compute log(sum(exp(x) for x in numlist))
    # mitigate the problem of underflow
    maxx = max(numlist)
    result = maxx + log(sum([exp(x-maxx) for x in numlist]))
    return result

def pseudo_log(x):
    return log(x) if x>0 else min_llh

class Alphabet():
    def __init__(self,K,J,data_struct):
        self.K = K
        self.J = J
        self.data_struct = data_struct # data_struct is a 3-way nested list
    
    def get_cassette_alphabet(self,k):
        a_list = self.data_struct[k]   
        def __get_alphabet(j):
            if j == 0:
                return [[x] for x in a_list[j]]
            else:
                prev = __get_alphabet(j-1)
                curr = []
                for x in prev:
                    for y in a_list[j]:        
                        curr.append(x+[y])
                return curr
        return [tuple(x) for x in __get_alphabet(self.J-1)]

class AlleleTable():
    def __init__(self,K,J,data_struct,alphabet):
    # K: the number of cassettes
    # J: the number of sites per cassette
    # data_struct: a mapping: cell_name -> (cassette -> (cassette_state -> count))
    # alphabet: an instance of class Alphabet; must have the same K and J
        self.K = K
        self.J = J
        self.data_struct = data_struct
        self.alphabet = alphabet
        
    def get(self,w,k,x):
        # w: a cell name/label
        # k: a cassette index
        # x: a state of the specified cassette; must be a tuple of length J
        assert len(x) == self.J
        return self.data_struct[w][k][x]
    
    def get_all_counts(self,w,k):
        # w: a cell name/label
        # k: a cassette index
        return self.data_struct[w][k]

class Param():
    def __init__(self,names,values,lower_bounds,upper_bounds):
        self.names = names
        self.values = values
        self.lower_bounds = lower_bounds
        self.upper_bounds = upper_bounds
        self.name2bound = {n:(l,u) for (n,l,u) in zip(names,lower_bounds,upper_bounds)}

    def get_value(self,pname):
        for i,(n,v) in enumerate(zip(self.names,self.values)):
            if n == pname:
                return v
        return None # only get here if pname is not found in self.names

class Count_base_model(Virtual_solver):
    def __init__(self,treeList,data,prior,params):
    # params in an instance of the Param class
    # data is a dictionary of multiple data modules; it MUST have 'alleleTable'
    # this base class only uses allele_table, but any derived class can add more attributes for joint likelihood computation
        self.data = data
        self.num_cassettes = data['alleleTable'].K
        self.site_per_cassette = data['alleleTable'].J
        self.params = params
        self.trees = []
        self.num_edges = 0
        for tree in treeList:
            tree_obj = read_tree_newick(tree)
            #tree_obj.suppress_unifurcations()
            self.num_edges += len(list(tree_obj.traverse_postorder()))
            self.trees.append(tree_obj)

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
        return self.params

    def score_tree(self,strategy={'ultra_constr':False,'fixed_params':None,'fixed_brlen':None}):
    # override Virtual solver
        ultra_constr = strategy['ultra_constr']
        fixed_params = strategy['fixed_params']
        fixed_brlen = strategy['fixed_brlen']
        nllh,status = self.optimize(initials=1,verbose=-1,ultra_constr=ultra_constr,fixed_params=fixed_params,fixed_brlen=fixed_brlen)
        score = None if nllh is None else -nllh
        return score,status
    
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

    def ini_brlens(self):
        x = [random() * (self.dmax/2 - 2*self.dmin) + 2*self.dmin for i in range(self.num_edges)]        
        idx = 0
        for tree in self.trees:
            for node in tree.traverse_postorder():
                if node.edge_length is not None:
                    x[idx] = max(2*self.dmin,node.edge_length)
                idx += 1 
        return x

    def ini_one_param(self,pname):
        # here we set a simple initialization scheme for all params
        # in most cases this method should be overrided in a class that inhirits from this base class
        lower_bound,upper_bound = self.params.name2bounds[pname]
        return random()*(upper_bound-lower_bound)+lower_bound

    def ini_all_params(self,fixed_params={}):
        # initialize branch lengths
        x = self.ini_brlens()
        # initialize the params specified in self.params
        for p in self.params.names:
            x += [fixed_params[p]] if p in fixed_params else [ini_one_param]
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
                node.edge_length = x[i]
                i += 1
        # other x2params                
        for j,p in enumerate(self.params.names):
            self.params.values[j] = fixed_params[p] if p in fixed_params else x[i]
            i += 1
    
    def Psi(self,c_node,k,j,alpha,beta):
        # Layer 1 transition probabilities
        # This is a placeholder (i.e. non-informative model) for this function in the base class
        # MUST be overrided in any derived class!
        return 1

    def log_Psi_cassette(self,c_node,k,alpha,beta):
        # compute the log-transformed transition probability of cassette k  
        # return None if the probability is 0
        log_trans_p = 0
        for j in range(self.data['alleleTable'].J):
            p = self.Psi(c_node,k,j,alpha[j],beta[j])
            if p > 0:
                log_trans_p += log(p) #####*****#####
            else:    
                log_trans_p = None
                break
        return log_trans_p

    def Gamma(self,k,x,c):
        # Layer 2 transition probabilities
        # A placeholder (i.e. non-informative model) for this function in the base class
        # MUST be overrided in any derived class!
        return 1

    def llh_alleleTable(self):
        # compute the log-likelihood of the allele table
        K = self.data['alleleTable'].K
        allele_table = self.data['alleleTable']
        root_state = tuple([0]*self.data['alleleTable'].J)
        for tree in self.trees:
            for node in tree.traverse_postorder():
                node.in_llh = [{} for _ in range(K)] # a list of dictionaries
        
        total_llh = 0
        for k in range(K):
            allele_list = self.data['alleleTable'].alphabet.get_cassette_alphabet(k)
            for tree in self.trees:
                for node in tree.traverse_postorder():
                    #node.in_llh = [{} for _ in range(K)] # a list of dictionaries
                    if node.is_leaf():
                        c = allele_table.get_all_counts(node.label,k) # c is a dictionary of allele -> count
                        for alpha in allele_list: 
                            trans_p = self.Gamma(k,alpha,c) # transition probability
                            if trans_p>0:
                                node.in_llh[k][alpha] = log(trans_p)
                    else:   
                        for alpha in allele_list:
                            llh = 0
                            for c_node in node.children:
                                llh_list = []                                    
                                for beta in c_node.in_llh[k]:
                                    log_trans_p = self.log_Psi_cassette(c_node,k,alpha,beta)
                                    if log_trans_p is not None:
                                        llh_list.append(log_trans_p + c_node.in_llh[k][beta])
                                if llh_list: # if the list is not empty
                                    llh += log_sum_exp(llh_list)
                                else:
                                    llh = None
                                    break   
                            if llh is not None:
                                node.in_llh[k][alpha] = llh
                #total_llh += sum([tree.root.in_llh[k][root_state] for k in range(self.data['alleleTable'].K)])
                total_llh += tree.root.in_llh[k][root_state]
        return total_llh

    def llh_alleleTable_edge(self,node,k,alpha):
        # assume in_llh has been computed for every node
        llh_list = []
        for beta in node.in_llh[k]:
            log_trans_p = self.log_Psi_cassette(node,k,alpha,beta)
            if log_trans_p is not None: 
                llh_list.append(node.in_llh[k][beta] + log_trans_p)        
        return log_sum_exp(llh_list) if llh_list else None       

    def negative_llh(self):
        # compute the negative log-likelihood of all data modules
        # in the base class, only the allele table is included
        # if a derived class has data modules other than the allele table, 
        # this function MUST be overrided to compute the joint-llh 
        # of all available data modules
        return -self.llh_alleleTable()

    def optimize(self,solver,initials=20,fixed_brlen=None,fixed_params={},verbose=1,max_trials=100,random_seeds=None,ultra_constr=False,**solver_opts):
    # random_seeds can either be a single number or a list of intergers where len(random_seeds) = initials
    # verbose level: 1 --> show all messages; 0 --> show minimal messages; -1 --> completely silent
    # fixed_brlen is a list of t dictionaries, where t is the number of trees in self.trees, each maps a tuple (a,b) to a number. Each pair a, b is a tuple of two leaf nodes whose LCA define the node for the branch above it to be fixed
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
                if verbose >= 0:
                    if ultra_constr:
                        print("Numerical optimization started with ultrametric constraint (default)")
                    else:      
                        print("Numerical optimization started without ultrametric constraint [deprecated]")
                # read in fixed_brlen and mark the tree nodes
                for t,tree in enumerate(self.trees):
                    for node in tree.traverse_postorder():
                        node.mark_fixed=False
                    if fixed_brlen is None:
                        continue
                    fixed_nodes = find_LCAs(tree,list(fixed_brlen[t].keys()))        
                    for i,(a,b) in enumerate(fixed_brlen[t]):
                        u = fixed_nodes[i]
                        u.edge_length = fixed_brlen[t][(a,b)]
                        u.mark_fixed = True
                if solver == 'Scipy':
                    scipy_options = solver_opts if solver_opts else DEFAULT_scipy_options
                    scipy_options['disp'] = (verbose>0)
                    nllh,status = self.scipy_optimization(randseed,fixed_params=fixed_params,ultra_constr=ultra_constr,scipy_options=scipy_options)
                else:
                    EM_options = solver_opts if solver_opts else DEFAULT_EM_opts
                    nllh,status = self.EM_optimization(randseed,fixed_params=fixed_params,ultra_constr=ultra_constr,EM_options=EM_options)
                
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

    def Estep_in_llh(self):
        self.llh_alleleTable()

    def Estep_out_llh(self):
        K = self.data['alleleTable'].K
        allele_table = self.data['alleleTable']
        root_state = tuple([0]*self.data['alleleTable'].J)
        
        for tree in self.trees:
            for node in tree.traverse_preorder():
                node.out_llh = [{} for _ in range(K)] # a list of dictionaries

        for k in range(K):
            allele_list = self.data['alleleTable'].alphabet.get_cassette_alphabet(k)
            for tree in self.trees:
                for node in tree.traverse_preorder():
                    if node.is_root():
                        node.out_llh[k][root_state] = 0
                    elif node.parent.is_root():
                        for alpha in allele_list:
                            log_trans_p = self.log_Psi_cassette(node,k,root_state,alpha) 
                            if log_trans_p is not None:
                                node.out_llh[k][alpha] = log_trans_p
                    else:
                        for alpha in allele_list:
                            llh_list = []                                    
                            for beta in node.parent.out_llh[k]:
                                log_trans_p = self.log_Psi_cassette(node,k,beta,alpha)
                                if log_trans_p is None:
                                    continue
                                par_out_llh = node.parent.out_llh[k][beta]
                                sum_sisters_in_llh = 0
                                for w in node.parent.children:
                                    if w is node:
                                        continue
                                    curr_llh = self.llh_alleleTable_edge(w,k,beta)    
                                    if curr_llh is not None:
                                        sum_sisters_in_llh += curr_llh
                                    else:
                                        sum_sisters_in_llh = None
                                        break
                                if sum_sisters_in_llh is not None:
                                    llh_list.append(log_trans_p+par_out_llh+sum_sisters_in_llh)
                            if llh_list: # if the list is not empty
                                node.out_llh[k][alpha] = log_sum_exp(llh_list)

    def Estep_posterior(self):
    # compute the log-transformed of all posterior probabilities
    # TODO: compute S0-4 as in EM_solver of PMM_original
        K = self.data['alleleTable'].K
        allele_table = self.data['alleleTable']
        root_state = tuple([0]*self.data['alleleTable'].J)
        
        for tree in self.trees:
            for node in tree.traverse_preorder():
                node.log_posterior = [{} for _ in range(K)] # a list of dictionaries

        for k in range(K):
            allele_list = self.data['alleleTable'].alphabet.get_cassette_alphabet(k)
            for tree in self.trees:
                for node in tree.traverse_preorder():
                    if node.is_root():
                        node.log_posterior[k][root_state] = 1
                    else:
                        for alpha in allele_list:
                            in_llh = node.in_llh[k][alpha] if alpha in node.in_llh[k] else None
                            out_llh = node.out_llh[k][alpha] if alpha in node.out_llh[k] else None
                            total_llh = tree.root.in_llh[k][root_state]
                            if (in_llh is not None) and (out_llh is not None):
                                node.log_posterior[k][alpha] = in_llh + out_llh - total_llh

    def Estep(self):
        self.Estep_in_llh()
        self.Estep_out_llh()
        self.Estep_posterior()
   
    def Closed_formed_optimal(param_name):
        # A placeholder for this function in the base class
        # MUST be overrided in any derived class!
        return 0

    def Mstep(self,fixed_params={},verbose=1,local_brlen_opt=True,ultra_constr_cache=None,eps_nu=1e-5,eps_s=1e-6):
        # assume that Estep has been performed so that all nodes have S0-S4 attributes
        # output: optimize all parameters: branch lengths, phi, and nu
        # verbose level: 1 --> show all messages; 0 --> show minimal messages; -1 --> completely silent        
       
        ################## IN PROGRESS! ######################
        raise("I AM UNDERDEVELOPMENT. PLEASE DON'T CALL ME!!!")

        # optimize nu and all branch lengths
        N = self.num_edges-self.num_polytomy_mark
        if local_brlen_opt:
            for tree in self.trees:
                N -= len([node for node in tree.traverse_postorder() if node.mark_fixed])
        S0 = np.zeros(N)
        S1 = np.zeros(N)
        S2 = np.zeros(N)
        S3 = np.zeros(N)
        S4 = np.zeros(N)
        d_ini = np.zeros(N)
        i = 0
        for tree in self.trees:
            for v in tree.traverse_postorder():
                if not v.polytomy_mark and not (v.mark_fixed and local_brlen_opt):    
                    s = [sum(v.S0),sum(v.S1),sum(v.S2),sum(v.S3),sum(v.S4)]
                    s = [max(eps_s,x) for x in s]
                    s = [x/sum(s)*self.numsites for x in s]
                    S0[i],S1[i],S2[i],S3[i],S4[i] = s
                    d_ini[i] = v.edge_length
                    i += 1
        
        def __optimize_brlen(nu,verbose=False): # nu is a single number
            var_d = cp.Variable(N,nonneg=True) # the branch length variables
            C0 = -(nu+1)*S0.T @ var_d
            C1 = -nu*S1.T @ var_d + S1.T @ cp.log(1-cp.exp(-var_d)) 
            C2 = S2.T @ cp.log(1-cp.exp(-nu*var_d)) if sum(S2) > 0 and nu > eps_nu else 0 
            C3 = -nu*S3.T @ var_d
            C4 = S4.T @ cp.log(1-cp.exp(-nu*var_d)) if sum(S4) > 0 and nu > eps_nu else 0

            objective = cp.Maximize(C0+C1+C2+C3+C4)
            constraints = [np.zeros(N)+self.dmin <= var_d, var_d <= np.zeros(N)+self.dmax]             
            if ultra_constr_cache is not None:
                M,b = ultra_constr_cache
                constraints += [np.array(M) @ var_d == np.array(b)]
            prob = cp.Problem(objective,constraints)
            prob.solve(verbose=False,solver=cp.MOSEK)
            return var_d.value,prob.status
       
        def __optimize_nu(d): # d is a vector of all branch lengths
            var_nu = cp.Variable(1,nonneg=True) # the nu variable
            C0 = -(var_nu+1)*S0.T @ d
            C1 = -var_nu*S1.T @ d + S1.T @ cp.log(1-cp.exp(-d))
            C2 = S2.T @ cp.log(1-cp.exp(-var_nu*d)) if sum(S2) > 0 else 0
            C3 = -var_nu*S3.T @ d
            C4 = S4.T @ cp.log(1-cp.exp(-var_nu*d)) if sum(S4) > 0 else 0
            objective = cp.Maximize(C0+C1+C2+C3+C4)
            prob = cp.Problem(objective)
            prob.solve(verbose=False,solver=cp.MOSEK)
            return var_nu.value[0],prob.status

        nIters = 1
        nu_star = self.params.nu
        for r in range(nIters):
            if verbose > 0:
                print("Optimizing branch lengths. Current phi: " + str(phi_star) + ". Current nu:" + str(nu_star))
            try:
                d_star,status_d = __optimize_brlen(nu_star,verbose=False)
            except:
                d_star = d_ini
                status_d = "failure"
            if status_d == "infeasible": # should only happen with local EM 
                return False,"d_infeasible"
            if not optimize_nu:
                if verbose > 0:
                    print("Fixing nu to " + str(self.params.nu))
                nu_star = self.params.nu
                status_nu = "optimal"    
            else:    
                if verbose > 0:
                    print("Optimizing nu")
                try:
                    nu_star,status_nu = __optimize_nu__(d_star)                 
                except:
                    status_nu = "failure"
        # place the optimal value back to params
        self.params.phi = phi_star
        self.params.nu = nu_star
        i = 0
        for tree in self.trees:
            for node in tree.traverse_postorder():
                if not node.polytomy_mark and not (node.mark_fixed and local_brlen_opt):    
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
    
    def EM_optimization(self,verbose=1,fixed_params={},ultra_constr=True,EM_options=DEFAULT_EM_options):
        # IMPORTANT: this EM algorithm ONLY optimizes the likelihood of the allele table !
        # assume that az_partition has been performed
        # optimize all parameters: branch lengths, phi, and nu
        # if optimize_phi is False, it is fixed to the original value in params.phi
        # the same for optimize_nu
        # caution: this function will modify params in place!
        # verbose level: 1 --> show all messages; 0 --> show minimal messages; -1 --> completely silent
        maxIter = EM_options['max_iter']
        pre_llh = self.llh_alleleTable()
        if verbose >= 0:
            print("Initial nllh: " + str(pre_llh))
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
            self.Estep() #####*****#####
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
                    return -pre_llh, em_iter,status
                elif verbose >= 0:    
                    print("Warning: EM failed to optimize parameters in one Mstep.")                
            curr_llh = self.llh_alleleTable()
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
        x0 = self.ini_all(fixed_params=fixed_params)
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
            M = self.ultrametric_constr()
            constraints.append(optimize.LinearConstraint(csr_matrix(M),[0]*len(M),[0]*len(M),keep_feasible=False))
        out = optimize.minimize(nllh, x0, method="SLSQP", options=scipy_options, bounds=bounds,constraints=constraints)
        if out.success:
            self.x2params(out.x,fixed_params=fixed_params)
            params = self.params
            f = out.fun
        else:
            f,params = None,None
        status = "optimal" if out.success else out.message
        return f,status
