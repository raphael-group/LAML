problin_libs/Topology_search.pyfrom math import log,isclose,exp
from random import choice, shuffle, random
from problin_libs import *
from treeswift import *
from problin_libs.EM_solver import EM_solver


class Topology_search:
    def __init__(self,treeTopo,solver,data={},prior={},params={},T_cooldown=20,alpha_cooldown=0.9):
        self.treeTopo = treeTopo # treeTopo is a newick string
        self.solver = solver     # solver is a solver definition
        self.params = params   
        self.data = data
        self.prior = prior        
        self.__renew_tree_obj__()
        # identify polytomies
        self.has_polytomy = False
        for node in self.tree_obj.traverse_preorder():
            if len(node.children) > 2:
                self.has_polytomy = True
                break        
        self.treeTopo = self.tree_obj.newick()
        self.T_cooldown = T_cooldown
        self.alpha_cooldown = alpha_cooldown
        self.b = 1/(1-1/self.alpha_cooldown**self.T_cooldown)
        self.a = -self.b/self.alpha_cooldown**self.T_cooldown

    def __renew_tree_obj__(self):
        self.tree_obj = read_tree_newick(self.treeTopo)
        self.tree_obj.suppress_unifurcations()
        #self.__mark_polytomies__()
    
    def get_solver(self):
        return self.solver(self.treeTopo,self.data,self.prior,self.params)
                
    def update_from_solver(self,mySolver):
        self.treeTopo = mySolver.get_tree_newick()
        self.params = mySolver.get_params()

    def __mark_polytomies__(self,eps_len=1e-3):
        # mark and resolve all polytomies in self.tree_obj
        self.has_polytomy = False
        for node in self.tree_obj.traverse_preorder():
            node.mark = False
            if len(node.children) > 2:
                #for c in node.children:
                #    c.mark = True
                self.has_polytomy = True
        self.tree_obj.resolve_polytomies()
        for node in self.tree_obj.traverse_preorder():
            if not hasattr(node,'mark'):
                node.mark = True
                node.edge_length = eps_len  
                #self.has_polytomy = True              
        self.treeTopo = self.tree_obj.newick()        

    def __accept_proposal__(self,curr_score,new_score,t):
        if new_score > curr_score:
            return True
        T = max(1e-6,self.a*self.alpha_cooldown**t + self.b)
        p = min(exp((new_score-curr_score)/T),1)
        return random() < p

    def search(self,maxiter=100,verbose=False,nreps=1,strategy=DEFAULT_STRATEGY):
        original_topo = self.treeTopo
        original_params = self.params
        nni_replicates = [(None,None)]*nreps
        for i in range(nreps):
            # resolve polytomies
            if verbose:
                print("Performing nni-search " + str(i+1))
            self.treeTopo = original_topo
            self.params = original_params
            self.__renew_tree_obj__()
            self.__mark_polytomies__()
            #topo_list1 = []
            #topo_list2 = []
            if strategy['resolve_search_only']:
                if verbose:
                    print("Only perform local nni moves to resolve polytomies")
                if self.has_polytomy:
                    topo_list,best_score = self.__search_one__(strategy,maxiter=maxiter,verbose=verbose,only_marked=True)
                else: # score this tree topology (optimize all numerical params)
                    mySolver = self.get_solver()
                    best_score = mySolver.score_tree(strategy=strategy)
                    self.update_from_solver(mySolver)
                    topo_list = [(self.treeTopo,best_score)]            
            #if not strategy['only_marked']:    
            else:    
                if verbose:
                    print("Perform nni moves for full topology search")
                topo_list,best_score = self.__search_one__(strategy,maxiter=maxiter,verbose=verbose,only_marked=False)
            #topo_list = [(x,y,'resolve_polytomies') for x,y in topo_list1] + [(x,y,'full_search') for x,y in topo_list2]
            nni_replicates[i] = (best_score,topo_list)
        return nni_replicates
    
    def __search_one__(self,strategy,maxiter=100,verbose=False,only_marked=False):
        # optimize branch lengths and other parameters for the starting tree
        mySolver = self.get_solver()
        #strategy_copy = {x:strategy[x] for x in strategy}
        #strategy_copy['optimize'] = True
        curr_score = mySolver.score_tree(strategy=strategy) 
        if verbose:
            print("Initial score: " + str(curr_score))
        self.update_from_solver(mySolver)
        topo_list = [(self.treeTopo,curr_score)]            
        # perform nni search
        for nni_iter in range(maxiter):
            if verbose:
                print("NNI Iter:", nni_iter)
            new_score,n_attempts,success = self.single_nni(curr_score,nni_iter,strategy,only_marked=only_marked)
            if not success:
                break
            curr_score = new_score
            topo_list.append((self.treeTopo,curr_score))
            if verbose:
                print("Current score: " + str(curr_score))
        final_score = curr_score
        if verbose:
            print("Final score: " + str(final_score))
        return topo_list,final_score 
    
    def single_nni(self,curr_score,nni_iter,strategy,only_marked=False):
        branches = []
        for node in self.tree_obj.traverse_preorder():
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
            n_attempts += 1
        return score,n_attempts,took
    
    def apply_nni(self,u,curr_score,nni_iter,strategy):
        # apply nni [DESTRUCTIVE FUNCTION! Changes tree inside this function.]
        v = u.get_parent()
        for node in v.child_nodes():
            if node != u:
                w = node
                break                
        #mySolver = self.get_solver()        
        #curr_score = mySolver.score_tree(strategy=strategy)
        u_children = u.child_nodes()
        # shuffle the order of the nni moves
        shuffle(u_children)
        nni_moves = []

        for u_child in u_children:
            u_child.set_parent(v)
            u.remove_child(u_child)
            v.add_child(u_child)

            w.set_parent(u)
            v.remove_child(w)
            u.add_child(w)

            mySolver = self.solver(self.tree_obj.newick(),self.data,self.prior,self.params)
            new_score = mySolver.score_tree(strategy=strategy)

            #if new_score > curr_score or isclose(new_score,curr_score,rel_tol=1e-3): # accept the new tree and params
            if self.__accept_proposal__(curr_score,new_score,nni_iter): # accept the new tree and params
                #print(curr_score,new_score,"accept")
                self.update_from_solver(mySolver)
                return True,new_score
            
            #print(curr_score,new_score,"reject")
            
            # Score doesn't improve --> reverse to the previous state
            u_child.set_parent(u)
            v.remove_child(u_child)
            u.add_child(u_child)
            
            w.set_parent(v)
            u.remove_child(w)
            v.add_child(w)            

        # no move accepted
        return False,curr_score