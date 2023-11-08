from math import log,isclose,exp
import timeit
from random import choice, shuffle, random
from scmail_libs import *
from treeswift import *
from scmail_libs.EM_solver import EM_solver
from scmail_libs.Topology_search import Topology_search
from copy import deepcopy
from scmail_libs.lca_lib import find_LCAs
from multiprocessing import Pool

class Topology_search_parallel(Topology_search):
    def single_nni(self,curr_score,nni_iter,strategy,only_marked=False):
        all_nni_moves = self.list_all_nni(strategy,only_marked=only_marked)
        N = len(all_nni_moves)
        batch_size = 16
        curr_start_idx = 0
        curr_end_idx = min(batch_size,N)
        retry_nni_moves = []
        checked_all = False
        took = False
        curr_queue = all_nni_moves
        while True:
            subset_nni_moves = curr_queue[curr_start_idx:curr_end_idx]
            with Pool() as pool:
                nni_results = pool.map(self.apply_nni,subset_nni_moves) 
            for i,nni_result in enumerate(nni_results):
                if nni_result['status'] == "optimal":
                    new_score = nni_result['score']
                    if self.__accept_proposal__(curr_score,new_score,nni_iter): # accept the new tree and params               
                        u,v,u_child,w = nni_result['cache']
                        u_child.set_parent(v)
                        u.remove_child(u_child)
                        v.add_child(u_child)
                        w.set_parent(u)
                        v.remove_child(w)
                        u.add_child(w)
                        self.update_from_solver(nni_result['mySolver'])
                        self.treeList_obj = nni_result['treeList_obj']
                        took = True
                        break
                elif not checked_all:
                    nwk_str,score_tree_strategy,(u,v,u_child,w) = subset_nni_moves[i]
                    score_tree_strategy['fixed_brlen'] = {}
                    retry_move = (nwk_str,score_tree_strategy,(u,v,u_child,w)) 
                    retry_nni_moves.append(retry_move)
            if checked_all or took:
                break    
            if not checked_all and curr_end_idx == N: # reach the end of the main queue  
                checked_all = True
                curr_queue = retry_nni_moves # switch queue
                curr_start_idx = 0 # restart index
            else:    
                curr_start_idx += batch_size
            curr_end_idx = min(curr_start_idx+batch_size,len(curr_queue))
        return new_score,curr_end_idx,took    
   
    def apply_nni(self,arguments):
        treeTopoList,score_tree_strategy,cache = arguments
        mySolver = self.solver(treeTopoList,self.data,self.prior,self.params)            
        score,status = mySolver.score_tree(strategy=score_tree_strategy)
        nni_result = {'mySolver':mySolver,'score':score,'status':status,'cache':cache,'treeList_obj':self.treeList_obj}
        return nni_result

    def list_all_nni(self,strategy,only_marked=False):    
        branches = []
        for tree in self.treeList_obj:
            for node in tree.traverse_preorder():
                if node.is_leaf() or node.is_root():
                    continue
                if not only_marked or node.mark:
                    branches.append(node)
        shuffle(branches)        
        all_nni_moves = []
        for u in branches:        
            v = u.get_parent()
            for node in v.child_nodes():
                if node != u:
                    w = node
                    break                
            u_children = u.child_nodes()
            score_tree_strategy = deepcopy(strategy)
            score_tree_strategy['fixed_brlen'] = {}
            if strategy['local_brlen_opt']:
                score_tree_strategy['fixed_nu'] = self.params['nu'] 
                score_tree_strategy['fixed_phi'] = self.params['phi'] 
                free_branches = set(u.child_nodes() + v.child_nodes() + [v])
                fixed_branches_anchors = []
                for tree in self.treeList_obj: 
                    for node in tree.traverse_postorder():
                        if node.is_leaf():
                            node.anchors = (node.label,node.label)
                        else:
                            C = node.child_nodes()
                            a = C[0].anchors[0]
                            b = C[-1].anchors[0]
                            node.anchors = (a,b)
                        if not node in free_branches:
                            fixed_branches_anchors.append(node.anchors)
                fixed_branches = []            
                for treeTopo in self.treeTopoList:
                    tree = read_tree_newick(treeTopo)
                    fixed_branches += find_LCAs(tree,fixed_branches_anchors)
                fixed_brlen = {}
                for i,node in enumerate(fixed_branches):
                    fixed_brlen[fixed_branches_anchors[i]] = node.edge_length
                score_tree_strategy['fixed_brlen'] = fixed_brlen
            shuffle(u_children)
            for u_child in u_children:
                # apply the nni move
                u_child.set_parent(v)
                u.remove_child(u_child)
                v.add_child(u_child)

                w.set_parent(u)
                v.remove_child(w)
                u.add_child(w)
                # get the new trees' strings
                nwk_strs = [tree.newick() for tree in self.treeList_obj]
                all_nni_moves.append((nwk_strs,score_tree_strategy,(u,v,u_child,w)))
                # turn back the move
                u_child.set_parent(u)
                v.remove_child(u_child)
                u.add_child(u_child)
                
                w.set_parent(v)
                u.remove_child(w)
                v.add_child(w)            
        return all_nni_moves             
