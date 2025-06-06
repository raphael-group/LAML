from math import log,isclose,exp
import timeit
from random import choice, shuffle, random
from laml_libs import *
from treeswift import *
from laml_libs.EM_solver import EM_solver
from copy import deepcopy
from laml_libs.lca_lib import find_LCAs

class Topology_search:
    def __init__(self,treeTopoList,solver,data={},prior={},params={},T_cooldown=20,alpha_cooldown=0.9):
        self.treeTopoList = treeTopoList # treeTopoList is a newick string
        self.solver = solver     # solver is a solver definition
        self.params = params   
        self.data = data
        self.prior = prior        
        self.__renew_treeList_obj__()
        # identify polytomies
        self.has_polytomy = False
        for tree in self.treeList_obj:
            for node in tree.traverse_preorder():
                if len(node.children) > 2:
                    self.has_polytomy = True
                    break        
        self.treeTopoList = []
        for tree in self.treeList_obj:
            self.treeTopoList.append(tree.newick())
        self.T_cooldown = T_cooldown
        self.alpha_cooldown = alpha_cooldown
        self.b = 1/(1-1/self.alpha_cooldown**self.T_cooldown)
        self.a = -self.b/self.alpha_cooldown**self.T_cooldown

    def __renew_treeList_obj__(self):
        self.treeList_obj = []
        for treeTopo in self.treeTopoList:
            tree = read_tree_newick(treeTopo)
            tree.suppress_unifurcations()
            self.treeList_obj.append(tree)
    
    def get_solver(self):
        return self.solver(self.treeTopoList,self.data,self.prior,self.params)
                
    def update_from_solver(self,mySolver):
        self.treeTopoList = mySolver.get_tree_newick()
        self.params = mySolver.get_params()
        if self.get_solver().solver_name == "fastEM_solver":
            self.__renew_treeList_obj__() # added Feb 28, 2025: I would rather pass in all the updated branch lengths
            # currently does not maintain marked nodes, so that it cannot be run to resolve polytomies!

    def __mark_polytomies__(self,eps_len=1e-3):
        # mark and resolve all polytomies in self.treeList_obj
        self.has_polytomy = False
        for tree in self.treeList_obj:
            for node in tree.traverse_preorder():
                node.mark = False
                if len(node.children) > 2:
                    self.has_polytomy = True
            tree.resolve_polytomies()
        
        for tree in self.treeList_obj:
            for node in tree.traverse_preorder():
                if not hasattr(node,'mark'):
                    node.mark = True
                    node.edge_length = eps_len  
        self.treeTopoList = [tree.newick() for tree in self.treeList_obj]

    def __accept_proposal__(self,curr_score,new_score,t):
        #print(f"curr_score: {curr_score}, new_score: {new_score}, relative improvement (positive): {new_score-curr_score}")
        if new_score > curr_score:
            #print(f"curr_score: {curr_score}, new_score: {new_score}, relative improvement (positive): {new_score-curr_score}")
            return True
        else:
            T = max(1e-12,self.a*self.alpha_cooldown**t + self.b)
            #T = max(1e-4,self.a*self.alpha_cooldown**t + self.b)
            p = min(exp((new_score-curr_score-1e-12)/T),1)
            #print(f"Acceptance probability: {p}")
            prob_accept = random() < p
            #if prob_accept:
            #    print(f"accept w.p. {p, prob_accept}")
            return prob_accept

        #relative_improvement = (curr_score - new_score)/abs(curr_score) # bigger is better, should be negative
        #print(f"curr_score: {curr_score}, new_score: {new_score}")
        #print(f"relative improvement: {relative_improvement}")
        #if relative_improvement < -0.001 and curr_score < new_score:
        #    return True
        #else:
        #    eps = 0.001
        #    T = max(eps,self.a*self.alpha_cooldown**t + self.b) # T decreases over time
        #    p = min(exp((new_score-curr_score-eps)/T),1) # p decreases over time
        #    print(f"Acceptance probability: {p}")
        #    return random() < p
        """"
        print(f"curr_score: {curr_score}, new_score: {new_score}")
        relative_improvement = (curr_score - new_score)/abs(curr_score)
        print(f"relative improvement: {relative_improvement}")
        if relative_improvement < 0.001 and curr_score < new_score: 
        #if new_score > curr_score:
            return True
        else:
            return False
        T = max(1e-12,self.a*self.alpha_cooldown**t + self.b)
        p = min(exp((new_score-curr_score-1e-12)/T),1)
        print(f"Acceptance probability: {p}")
        return random() < p
        """

    def search(self,resolve_polytomies=True,maxiter=100,verbose=False,nreps=1,strategy=DEFAULT_STRATEGY,checkpoint_file=None):
        original_topos = self.treeTopoList
        original_params = self.params
        #nni_replicates = [(None,None)]*nreps
        best_trees = None
        best_score = -float("inf")
        best_params = None
        for i in range(nreps):
            # resolve polytomies
            if verbose:
                print("Performing nni-search " + str(i+1), flush=True)
            self.treeTopoList = original_topos
            self.params = original_params
            self.__renew_treeList_obj__()
            if resolve_polytomies:
                self.__mark_polytomies__()
            if strategy['resolve_search_only']:
                if verbose:
                    print("Only perform local nni moves to resolve polytomies")
                if self.has_polytomy:
                    trees,score,params = self.__search_one__(strategy,maxiter=maxiter,verbose=verbose,only_marked=True,checkpoint_file=checkpoint_file)
                else: # score this tree topology (optimize all numerical params)
                    if verbose:
                        print("Found no polytomy to resolve. Optimizing numerical parameters without further topology search")
                    mySolver = self.get_solver()
                    score_tree_strategy = deepcopy(strategy)
                    score_tree_strategy['fixed_brlen'] = None
                    score_tree_strategy['nodes_to_recompute'] = None
                    score,status = mySolver.score_tree(strategy=score_tree_strategy)
                    self.update_from_solver(mySolver)
                    trees = self.treeTopoList
                    params = self.params
            else:    
                if verbose:
                    print("Perform nni moves for full topology search", flush=True)
                    if self.has_polytomy and resolve_polytomies:                        
                        print("Found polytomies in the input tree(s). Arbitrarily resolving them to obtain fully resolved initial tree(s).") 
                trees,score,params = self.__search_one__(strategy,maxiter=maxiter,verbose=verbose,only_marked=False,checkpoint_file=checkpoint_file)
            # The final optimization of parameters
            if verbose:
                print("Optimal topology found. Re-optimizing other parameters ...", flush=True)
            self.treeTopoList = trees
            self.params = params
            self.__renew_treeList_obj__()
            mySolver = self.get_solver()
            score_tree_strategy = deepcopy(strategy)
            score_tree_strategy['fixed_brlen'] = None
            score_tree_strategy['final_optimization'] = True
            score_tree_strategy['nodes_to_recompute'] = None
            score,status = mySolver.score_tree(strategy=score_tree_strategy)
            self.update_from_solver(mySolver)
            trees = self.treeTopoList
            params = self.params
            if verbose:
                print("Optimal score for this search: " + str(score), flush=True)
            # compare to the best_score of previous searches
            if score > best_score:
                best_score = score
                best_trees = trees    
                best_params = params

        # synchronization        
        self.treeTopoList = best_trees
        self.params = best_params
        
        return best_trees,best_score,best_params
    
    def __search_one__(self,strategy,maxiter=100,verbose=False,only_marked=False, checkpoint_file=None):
        # optimize branch lengths and other parameters for the starting tree
        mySolver = self.get_solver()
        score_tree_strategy = deepcopy(strategy)
        score_tree_strategy['fixed_brlen'] = None
        score_tree_strategy['nodes_to_recompute'] = None
        curr_score,status = mySolver.score_tree(strategy=score_tree_strategy)
        if verbose:
            if self.has_polytomy:
                print("Initial score (polytomies were arbitrarily resolved): " + str(curr_score), flush=True)
            else:
                print("Initial score: " + str(curr_score), flush=True)
        self.update_from_solver(mySolver)
        best_score = curr_score
        best_trees = self.treeTopoList
        best_params = self.params 
        # perform nni search
        for nni_iter in range(maxiter):
            if verbose:
                print("NNI Iter:", nni_iter, flush=True)
                start_time = timeit.default_timer()
            new_score,n_attempts,success = self.single_nni(curr_score,nni_iter,strategy,only_marked=only_marked)
            if verbose: 
                print(f"Number of trees checked: {n_attempts}", flush=True)
            if not success:
                break
            curr_score = new_score
            if curr_score > best_score:
                best_score = curr_score
                best_trees = self.treeTopoList
                best_params = self.params
            if verbose:
                print("Current score: " + str(curr_score), flush=True)
                stop_time = timeit.default_timer()
                print("Runtime (s):", stop_time - start_time, flush=True)
            if nni_iter % chkpt_freq == 0 and checkpoint_file is not None:
                with open(checkpoint_file, "a") as fout:
                    fout.write(f"NNI Iteration: {nni_iter}\n")
                    fout.write(f"Current newick tree: {best_trees}\n")
                    fout.write(f"Current negative-llh: {best_score}\n")
                    fout.write(f"Current dropout rate: {best_params['phi']}\n")
                    fout.write(f"Current silencing rate: {best_params['nu']}\n")
        if verbose:
            print("Best score for this search: " + str(best_score), flush=True)
        return best_trees,best_score,best_params 
    
    def single_nni(self,curr_score,nni_iter,strategy,only_marked=False,verbose=False):
        branches = []

        for tree in self.treeList_obj:
            for node in tree.traverse_preorder():
                if node.is_leaf() or node.is_root():
                    continue
                if not only_marked or node.mark:
                    branches.append(node)
        # branch ordering: random 
        shuffle(branches)

        took = False
        score = -float("inf")
        n_attempts = 0
        while not took and branches:
            u = branches.pop()
            took,score = self.apply_nni(u,curr_score,nni_iter,strategy)
            n_attempts += 2
        return score,n_attempts,took
    
    def apply_nni(self,u,curr_score,nni_iter,strategy):
        # apply nni [DESTRUCTIVE FUNCTION! Changes tree inside this function.]
        v = u.get_parent()
        for node in v.child_nodes():
            if node != u:
                w = node
                break                
        u_children = u.child_nodes()
        # shuffle the order of the nni moves
        shuffle(u_children)
        score_tree_strategy = deepcopy(strategy)
        score_tree_strategy['fixed_brlen'] = None
        score_tree_strategy['nodes_to_recompute'] = None

        if strategy['local_brlen_opt']:
            score_tree_strategy['fixed_nu'] = self.params['nu'] 
            score_tree_strategy['fixed_phi'] = self.params['phi'] 
            free_branches = set(u.child_nodes() + v.child_nodes() + [v])
            fixed_branches_anchors = [[] for _ in range(len(self.treeList_obj))]
            nodes_to_recompute = [[] for _ in range(len(self.treeList_obj))]
            for t,tree in enumerate(self.treeList_obj):
                for node in tree.traverse_postorder():
                    if node.is_leaf():
                        node.anchors = (node.label,node.label)
                    else:
                        C = node.child_nodes()
                        a = C[0].anchors[0]
                        b = C[-1].anchors[0]
                        node.anchors = (a,b)
                    if not node in free_branches:
                        fixed_branches_anchors[t].append(node.anchors)
                    else:
                        nodes_to_recompute[t].append(node.anchors)
            fixed_branches = [[] for _ in range(len(self.treeList_obj))]
            for t,treeTopo in enumerate(self.treeTopoList):
                tree = read_tree_newick(treeTopo)
                fixed_branches[t] += find_LCAs(tree,fixed_branches_anchors[t])
            fixed_brlen = [{} for _ in range(len(self.treeList_obj))]
            for t,B in enumerate(fixed_branches):
                for i,node in enumerate(B):
                    fixed_brlen[t][fixed_branches_anchors[t][i]] = node.edge_length

            # possibility 1: update the branch lengths in self.treeList_obj too

            score_tree_strategy['fixed_brlen'] = fixed_brlen
            score_tree_strategy['nodes_to_recompute'] = nodes_to_recompute

        for u_child in u_children:
            u_child.set_parent(v)
            u.remove_child(u_child)
            v.add_child(u_child)

            w.set_parent(u)
            v.remove_child(w)
            u.add_child(w)

            #print("before score tree (treeList_obj):", [tree.newick() for tree in self.treeList_obj])
            mySolver = self.solver([tree.newick() for tree in self.treeList_obj],self.data,self.prior,self.params)            
            new_score,status = mySolver.score_tree(strategy=score_tree_strategy)
            #print("after score tree:", [tree.newick() for tree in self.treeList_obj])
            if status != "optimal" and strategy['local_brlen_opt']:
                score_tree_strategy['fixed_brlen'] = None
                score_tree_strategy['nodes_to_recompute'] = None
                mySolver = self.solver([tree.newick() for tree in self.treeList_obj],self.data,self.prior,self.params)            
                new_score,status = mySolver.score_tree(strategy=score_tree_strategy)
            if self.__accept_proposal__(curr_score,new_score,nni_iter): # accept the new tree and params                
                self.update_from_solver(mySolver)
                #print("after update_from_solver:", [tree.newick() for tree in self.treeList_obj])
                return True,new_score
            
            # Score doesn't improve --> reverse to the previous state
            u_child.set_parent(u)
            v.remove_child(u_child)
            u.add_child(u_child)
            
            w.set_parent(v)
            u.remove_child(w)
            v.add_child(w)            

        # no move accepted
        return False,curr_score
