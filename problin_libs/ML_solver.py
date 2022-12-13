import dendropy
from treeswift import *
from math import log,exp
from random import random, seed
from scipy import optimize
import warnings

min_llh = -800

class Params:
    def __init__(self,nwkTree,nu,phi):
    # nwkTree: a newick tree (i.e. a string) with branch lengths
    # nu, phi: positive float numbers
        self.tree = read_tree_newick(nwkTree) # store a treewfit object in self.tree
        self.nu = nu
        self.phi = phi

class ML_solver:
    # at this stage, the tree topology must be given. Only branch lengths
    # and other parameters can be optimized
    def __init__(self,charMtrx,Q,nwkTree,nu=None,phi=None):
        self.charMtrx = charMtrx
        self.Q = Q
        self.params = Params(nwkTree,nu,phi)
    
    def az_partition(self,params):
    # partition the tree into edge-distjoint alpha-clades and z-branches
    # there is a different partition for each target-site
        numsites = len(self.charMtrx[next(iter(self.charMtrx.keys()))])
        # a-z decomposition
        # z-branches are given tag 'z'
        # each of other branches is given a tag 
        # alpha where alpha is the alpha-tree it belongs to
        for node in params.tree.traverse_postorder():
            if node.is_leaf():
                node.alpha = [self.charMtrx[node.label][site] if self.charMtrx[node.label][site] != 0 else 'z' for site in range(numsites)]
            else:
                C = node.children
                node.alpha = [None]*numsites
                for site in range(numsites):
                    S = set(c.alpha[site] for c in C)
                    R = S-set(['z','?'])
                    if 'z' in S or len(R)>1:
                        node.alpha[site] = 'z'
                    elif len(R) == 1:
                        node.alpha[site] = list(R)[0]    
                    else:
                        node.alpha[site] = "?"
    
    def compute_llh(self,params):
        # assume az_partition has been performed so
        # each node has the attribute node.alpha
        numsites = len(self.charMtrx[next(iter(self.charMtrx.keys()))])
        phi = params.phi
        nu = params.nu
        llh = [0]*numsites
        for node in params.tree.traverse_postorder():
            p = exp(-node.edge_length)
            node.L0 = [0]*numsites # L0 and L1 are stored in log-scale
            node.L1 = [0]*numsites
            for site in range(numsites):    
                if node.alpha[site] != 'z':
                    q = self.Q[site][node.alpha[site]] if node.alpha[site] != "?" else 1.0
                    if node.is_leaf():
                        #try:
                        node.L0[site] = nu*(-node.edge_length) + log(1-p) + log(q) + log(1-phi) if node.alpha[site] != "?" else log(1-(1-phi)*p**nu)
                        node.L1[site] = nu*(-node.edge_length) + log(1-phi) if node.alpha[site] != "?" else log(1-(1-phi)*p**nu)
                        #except:
                        #    print(1-p,q,1-phi,1-(1-phi)*p**nu)    
                    else:
                        C = node.children
                        L0 = (nu+1)*(-node.edge_length)
                        L1 = log(1-p)+log(q) + nu*(-node.edge_length)
                        for c in C:
                            L0 += c.L0[site]
                            L1 += c.L1[site]
                        L0 = exp(L0)
                        L1 = exp(L1)
                        if L1 == 0 and node.alpha[site]!="?":    
                            node.L0[site] = min_llh if L0 == 0 else log(L0)
                            node.L1[site] = min_llh
                        else:        
                            node.L0[site] = log(L0 + L1 + (1-p**nu)*int(node.alpha[site]=="?"))
                            node.L1[site] = log(L1 + (1-p**nu)*int(node.alpha[site]=="?"))
                            
                    if node.is_root() or node.parent.alpha[site] == 'z':
                        llh[site] += node.L0[site] 
                else:
                    llh[site] += (-node.edge_length + int(node.is_leaf())*log(1-phi))
        return sum(llh)         

    def optimize(self,initials=20,fixed_phi=None,fixed_nu=None,verbose=True,max_trials=100):
    # optimize tree branch lengths    
        warnings.filterwarnings("ignore")
        nwkt = self.params.tree
        num_edges = len(list(nwkt.traverse_postorder()))
        self.az_partition(self.params)

        def nllh(x): 
            for i, node in enumerate(nwkt.traverse_postorder()):
                node.edge_length = x[i]
            self.params.phi = x[-1] if fixed_phi is None  else fixed_phi
            self.params.nu = x[-2]  if fixed_nu is None else fixed_nu      
            return -self.compute_llh(self.params)
        
        numsites = len(self.charMtrx[next(iter(self.charMtrx.keys()))])
        x_star = []
        dmax = -log(1/numsites)*2
        dmin = -log(1-1/numsites)/2
        bounds = optimize.Bounds([dmin]*num_edges+[1e-10,1e-10],[dmax]*num_edges+[10,0.99],keep_feasible=True)
        
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
                x0 = [random() * (dmax - dmin) + dmin for i in range(num_edges)] + [random()*0.99,random()*0.99]
                out = optimize.minimize(nllh, x0, method="SLSQP", options={'disp':verbose,'iprint':3,'maxiter':1000}, bounds=bounds)
                if out.success:
                    all_failed = False
                    print("Optimal point found for initial " + str(i+1))
                    print("Optimal x: ", out.x)
                    print("Optimal f: ", out.fun)
                    if out.fun < f_star:
                        x_star = out.x
                        f_star = out.fun
                else:
                    print("Failed to optimize using initial point " + str(i+1)) 
            all_trials += initials
        
        # place results onto the tree
        for i,node in enumerate(nwkt.traverse_postorder()):
            node.edge_length = x_star[i]
        self.params.phi = x_star[-1] if fixed_phi is None  else fixed_phi
        self.params.nu = x_star[-2]  if fixed_nu is None else fixed_nu  
        return self.compute_llh(self.params)     

def main(): 
    from sequence_lib import read_sequences
    from ml_log import wrapper_felsenstein as wf_log
    
    k = 500
    m = 10
    Q = []
    for i in range(k):
        q = {j+1:1/m for j in range(m)}
        q[0] = 0
        Q.append(q)
    #T = "((a:0.0360971597765934,b:3.339535381892265)e:0.0360971597765934,(c:0.0360971597765934,d:3.339535381892265)f:0.0360971597765934)r:0.0;"
    #T = "((a,b)e,(c,d)f)r;"
    T = "((a:1,b:1):1,c:1):1;"
    #S = read_sequences("../MP_inconsistent/seqs_m10_k" + str(k) + ".txt")
    #msa = S[231]
    msa = dict()
    msa['a'] = [1, 1]
    msa['b'] = [1, 1]
    msa['c'] = [1, 1]
    #msa['d'][1] = '?'
    #msa['b'][0] = '?'
    ##msa['c'][0] = '?'
    ##msa['b'][1] = '?'
    ##msa = {'a':[1],'b':[1],'c':[1]}
    #print(wf_log(T, Q, msa, optimize_branchlengths=True))

    mySolver = ML_solver(msa,Q,T,nu=0)
    #mySolver.az_partition(mySolver.params)
    #print(mySolver.compute_llh(mySolver.params))
    print(mySolver.optimize(initials=3))
    print(mySolver.params.phi,mySolver.params.nu)
    print(mySolver.params.tree.newick())

if __name__ == "__main__":
    main()        
