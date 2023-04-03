from math import log,isclose
from random import choice, shuffle
from problin_libs import min_llh, eps, nni_conv_eps
from treeswift import *
from problin_libs.EM_solver import EM_solver

class Topology_search:
    def __init__(self,treeTopo,solver,data={},prior={},params={}):
        self.treeTopo = treeTopo # treeTopo is a newick string
        self.solver = solver     # solver is a solver definition
        self.params = params   
        self.data = data
        self.prior = prior        
        self.tree_obj = read_tree_newick(self.treeTopo)
    
    def __renew_tree_obj__(self):
        self.tree_obj = read_tree_newick(self.treeTopo)
    
    def get_solver(self):
        return self.solver(self.treeTopo,self.data,self.prior,self.params)

    def search(self,maxiter=100,verbose=False,nreps=1,strategy={}):
        nni_replicates = [(None,None)]*nreps
        for i in range(nreps):
            mySolver = self.get_solver()
            curr_score = mySolver.score_tree()
            topo_list = []
            
            for nni_iter in range(maxiter):
                if verbose:
                    print("NNI Iter:", nni_iter)
                new_score,n_attempts,success = self.single_nni()
                if not success:
                    break
                tstr = self.treeTopo
                topo_list.append((tstr, new_score))
                curr_score = new_score
            nni_replicates[i] = (curr_score,topo_list)    
        return nni_replicates    
    
    def single_nni(self):
        branches = [node for node in self.tree_obj.traverse_preorder() if (not node.is_leaf() and not node.is_root())]
        # random strategy 
        shuffle(branches)

        took = False
        n_attempts = 0
        while not took and branches:
            u = branches.pop()
            took,score = self.apply_nni(u)
            n_attempts += 1
        return score,n_attempts,took
    
    def apply_nni(self,u):
        # apply nni [DESTRUCTIVE FUNCTION! Changes tree inside this function.]
        v = u.get_parent()
        for node in v.child_nodes():
            if node != u:
                w = node
                break
        mySolver = self.get_solver()        
        curr_score = mySolver.score_tree()
        u_children = u.child_nodes()
        nni_moves = []

        for u_child in u_children:
            u_child.set_parent(v)
            u.remove_child(u_child)
            v.add_child(u_child)

            w.set_parent(u)
            v.remove_child(w)
            u.add_child(w)

            mySolver = self.solver(self.tree_obj.newick(),self.data,self.prior,self.params)
            new_score = mySolver.score_tree()

            if new_score >= curr_score: # accept the new tree and params
                self.treeTopo = mySolver.get_tree_newick()
                self.params = mySolver.get_params()
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
