from treeswift import *
from math import log,exp,sqrt, isclose
from random import random, seed, choice

class Virtual_solver:
    def __init__(self,treeList,data,prior,params):
    # virtual method. Should never be called!
        pass
    def get_tree_newick(self):
    # virtual method. Should never be called!
        pass    
    def get_params(self):
    # virtual method, should be overried in most cases
        return dict()   
    def score_tree(self,strategy={}):     
    # virtual method. Should never be called!
        pass    
