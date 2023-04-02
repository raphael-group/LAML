from math import log,isclose
from random import choice, shuffle
from problin_libs import min_llh, eps, nni_conv_eps
from treeswift import *

class Topology_search:
    def __init__(self,solver):
        self.solver = solver        
    
    def search(self,maxiter=100,verbose=False,trynextbranch=False,strategy="vanilla",branch_list='all',nreps=1,conv=0.2):        
        if branch_list == 'all':
            #branch_list = [node.label for node in self.solver.params.tree.traverse_postorder() if (not node.is_root() and not node.is_leaf())]
            branch_list = [node.label for node in self.solver.tree.traverse_postorder() if (not node.is_root() and not node.is_leaf())]
        nni_replicates = dict()
        original_tree = self.solver.tree.newick()

        for i in range(nreps):
            self.solver.tree = read_tree_newick(original_tree)
            nib = 0
            key_branches = []

            l2n = self.solver.tree.label_to_node(selection="all")
            for nlabel in branch_list:
                node = l2n[nlabel]
                if not node.is_leaf():
                    nib += 1
                    key_branches.append(node)

            topo_dict = {}
            seen = set()
            nni_iter = 0
            same = 0
            pre_llh = self.solver.score_tree()
            
            for j in range(maxiter):
                if verbose:
                    print("NNI Iter:", nni_iter)
                opt_score,n_attempts = self.single_nni(verbose, trynextbranch=trynextbranch, strategy=strategy, keybranches=key_branches)
                
                tstr = self.solver.tree.newick()
                topo_dict[nni_iter] = (tstr, opt_score)
                
                new_llh = self.solver.score_tree()

                if isclose(new_llh, pre_llh, rel_tol=1e-9, abs_tol=0.0):
                    same += 1
                else:
                    same = 0
               
                if new_llh-pre_llh < nni_conv_eps or n_attempts >= nib:
                    break
                pre_llh = new_llh
                nni_iter += 1
            nni_replicates[i] = (new_llh, topo_dict)
        return nni_replicates    
    
    def single_nni(self, verbose, trynextbranch=True, strategy="vanilla", keybranches=[]):
        branches = self.solver.score_branches(strategy, keybranches)
        if strategy == "random":
            shuffle(branches) 
        else:    
            branches = sorted(branches,key=lambda item:item[1])

        for x,y in branches:
            print(x.label,y)    
        took = False
        bidx = 0
        while not took and branches:
            u,u_score = branches.pop()
            took = self.apply_nni(u, verbose)
            bidx += 1
            if not trynextbranch:
                took = True 
        if verbose:
            print(bidx, " branch attempts.")
        llh = self.solver.score_tree()
        return llh, bidx 
    
    def apply_nni(self, u, verbose):
        # apply nni [DESTRUCTIVE FUNCTION! Changes tree inside this function.]
        v = u.get_parent()
        for node in v.child_nodes():
            if node != u:
                w = node
                break
        pre_llh = self.solver.score_tree()
        u_children = u.child_nodes()
        nni_moves = []

        for u_child in u_children:
            u_child.set_parent(v)
            u.remove_child(u_child)
            v.add_child(u_child)

            w.set_parent(u)
            v.remove_child(w)
            u.add_child(w)
            
            new_llh = self.solver.score_tree()

            if new_llh > pre_llh:
                # log likelihood improved
                return True
            elif new_llh == pre_llh:
                return True
            else:
                # REVERSE IF LIKELIHOOD IS NOT BETTER
                u_child.set_parent(u)
                v.remove_child(u_child)
                u.add_child(u_child)
                
                w.set_parent(v)
                u.remove_child(w)
                v.add_child(w)
                
                new_llh = self.solver.score_tree()
        return False
