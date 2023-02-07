from problin_libs.ML_solver import *
from math import exp,log

em_eps = 1e-6 # threshold to stop the EM algorithm 

class EM_solver(ML_solver):
    def Estep_in_llh(self,params):
        # assume az_partition has been performed so each node has the attribute node.alpha
        # compute the inside llh, store in L0 and L1 of each node
        phi = params.phi
        nu = params.nu
        for node in params.tree.traverse_postorder():
            p = exp(-node.edge_length)
            node.L0 = [0]*self.numsites # L0 and L1 are stored in log-scale
            node.L1 = [0]*self.numsites
            for site in range(self.numsites):    
                q = self.Q[site][node.alpha[site]] if node.alpha[site] not in ['?','z'] else 1.0
                if node.is_leaf():
                    if node.alpha[site] == "?":         
                        masked_llh = log(1-(1-phi)*p**nu)
                        node.L0[site] = node.L1[site] = masked_llh
                    elif node.alpha[site] == 'z':
                        node.L0[site] = (nu+1)*(-node.edge_length) + log(1-phi)
                        node.L1[site] = min_llh
                    else:    
                        node.L0[site] = nu*(-node.edge_length) + log(1-p) + log(q) + log(1-phi)
                        node.L1[site] = nu*(-node.edge_length) + log(1-phi)
                else:
                    C = node.children
                    l0 = l1 = 0
                    for c in C:
                        l0 += c.L0[site]
                        l1 += c.L1[site]
                    L0 = exp(l0+(nu+1)*(-node.edge_length))
                    if node.alpha[site] != 'z':
                        L0 += exp(l1 + log(1-p)+log(q) + nu*(-node.edge_length)) 
                    if node.alpha[site] == '?':
                        L0 += (1-p**nu)                       
                    L1 = 0 if node.alpha[site] == 'z' else (exp(l1+nu*(-node.edge_length)) + (1-p**nu)*int(node.alpha[site]=="?"))
                    node.L0[site] = min_llh if L0==0 else log(L0)
                    node.L1[site] = min_llh if L1==0 else log(L1)
    
    def lineage_llh(self,params):
        self.Estep_in_llh(params)
        return sum(params.tree.root.L0)
    
    def Estep_out_llh(self,params):
        # assume binary tree
        # assume az_parition and Estep_in_llh have been performed 
        # so that all nodes have `alpha`, `L0` and `L1` attribues
        # output: add the attributes `out0` and `out1` to each node
        # where v.out0 = P(~D_v,v=0) and v.out1 = P(~D_v,v=-1) 
        for v in params.tree.traverse_preorder():
            if v.is_root(): # base case
                # Auxiliary
                v.A = [0]*self.numsites
                #v.B = [min_llh]*self.numsites
                v.X = [-params.nu*v.edge_length + log(1-exp(-v.edge_length))]*self.numsites
                v.out_alpha = [{} for i in range(self.numsites)]
                #v.sis_tag = [None]*self.numsites
                # Main components    
                v.out0 = [-(1+params.nu)*v.edge_length]*self.numsites
                v.out1 = [log(1-exp(-v.edge_length*params.nu))]*self.numsites
            else:
                u = v.parent
                w = None 
                # get the sister
                for x in u.children:
                    if x is not v:
                        w = x
                        break
                if w is None:
                    print("w is none", [x.label for x in u.traverse_leaves()])
                # Auxiliary components
                v.A = [None]*self.numsites
                #v.B = [None]*self.numsites
                v.X = [None]*self.numsites
                v.out_alpha = [{} for i in range(self.numsites)]
                #v.sis_tag = [None]*self.numsites
                # Main components
                v.out0 = [None]*self.numsites
                v.out1 = [None]*self.numsites
                for site in range(self.numsites): 
                    # compute out0: P(~D_v,v=0)
                    v.out0[site] = u.out0[site] + w.L0[site] - (1+params.nu)*v.edge_length
                    v.A[site] = u.out0[site] + w.L0[site] 
                    v.X[site] = -params.nu*v.edge_length + log(1-exp(-v.edge_length)) + v.A[site]
                    # compute out1: P(~D_v,v=-1)
                    if w.alpha[site] == 'z': # z-branch
                        v.out1[site] = log(1-exp(-v.edge_length*params.nu)) + v.A[site]
                        #v.B[site] = min_llh
                        #v.sis_tag[site] = 'z'
                    elif w.alpha[site] == '?': # masked branch
                        #alpha0 = w.alpha[site] if w.alpha[site] != '?' else u.sis_tag[site]
                        #v.sis_tag[site] = alpha0
                        #u_out_alpha = u.A[site] - params.nu*u.edge_length + log(1-exp(-u.edge_length)) + log(self.Q[site][alpha0]) # log P(~D_u,u != 0 and u != -1)  
                        #if u.sis_tag[site] is None or alpha0 == u.sis_tag[site]:
                        #    u_out_alpha = log(exp(u_out_alpha)+exp(u.B[site]))
                        #v.B[site] = -params.nu*v.edge_length + w.L1[site]  + u_out_alpha
                        v.X[site] = log(exp(v.X[site])+exp(u.X[site]+w.L1[site]-params.nu*v.edge_length))
                        v.out1[site] = log((1-exp(-v.edge_length*params.nu))*(exp(v.A[site])+exp(u.X[site]+w.L1[site])) + exp(u.out1[site]+1))
                    else:
                        alpha0 = w.alpha[site] 
                        if alpha0 not in u.out_alpha[site]:
                            self.__out_alpha_up__(u,site,alpha0,params)
                        B = u.out_alpha[site][alpha0] + params.nu*(-u.edge_length) + w.L1[site]  
                        C = v.A[site] - params.nu*v.edge_length + log(1-exp(-v.edge_length)) + log(self.Q[site][alpha0])
                        v.out_alpha[site][alpha0] = log(exp(C) + exp(B))
                        v.X[site] = log(exp(v.X[site])+exp(B))
                        try:
                            v.out1[site] = log(1-exp(-v.edge_length*params.nu)) + log(exp(v.A[site]) + exp(w.L1[site]+u.out_alpha[site][alpha0]))
                        except:
                            print(v.edge_length,params.nu)    

    def __out_alpha_up__(self,node,site,alpha0,params):
        v = node
        path = []
        while not v.is_root():
            u = v.parent
            for x in u.children:
                if x is not v:
                    w = x
                    break
            if w.alpha[site] != '?' and w.alpha[site] != alpha0: # the branch above u is a z-branch
                v.out_alpha[site][alpha0] = v.A[site] - params.nu*v.edge_length + log(1-exp(-v.edge_length)) + log(self.Q[site][alpha0])
                break
            path.append(v)
            v = u
        if v.is_root(): # no z-branch found along the way
            v.out_alpha[site][alpha0] = -params.nu*v.edge_length + log(1-exp(-v.edge_length)) + log(self.Q[site][alpha0])             
        while path:
            v = path.pop()
            u = v.parent
            B = u.out_alpha[site][alpha0] + params.nu*(-u.edge_length) + w.L1[site]  
            C = v.A[site] - params.nu*v.edge_length + log(1-exp(-v.edge_length)) + log(self.Q[site][alpha0])
            v.out_alpha[site][alpha0] = log(exp(C) + exp(B))

    def Estep_out_llh_old(self,params):
        # assume binary tree
        # assume az_parition and Estep_in_llh have been performed 
        # so that all nodes have `alpha`, `L0` and `L1` attribues
        # output: add the attributes `out0` and `out1` to each node
        # where v.out0 = P(~D_v,v=0) and v.out1 = P(~D_v,v=-1) 
        for v in params.tree.traverse_preorder():
            if v.is_root(): # base case
                # Auxiliary
                v.A = [0]*self.numsites
                v.B = [min_llh]*self.numsites
                v.X = [-params.nu*v.edge_length + log(1-exp(-v.edge_length))]*self.numsites
                v.sis_tag = [None]*self.numsites
                # Main components    
                v.out0 = [-(1+params.nu)*v.edge_length]*self.numsites
                v.out1 = [log(1-exp(-v.edge_length*params.nu))]*self.numsites
            else:
                u = v.parent
                # get the sister
                for x in u.children:
                    if x is not v:
                        w = x
                        break
                # Auxiliary components
                v.A = [None]*self.numsites
                v.B = [None]*self.numsites
                v.X = [None]*self.numsites
                v.sis_tag = [None]*self.numsites
                # Main components
                v.out0 = [None]*self.numsites
                v.out1 = [None]*self.numsites
                for site in range(self.numsites): 
                    # compute out0: P(~D_v,v=0)
                    v.out0[site] = u.out0[site] + w.L0[site] - (1+params.nu)*v.edge_length
                    v.A[site] = u.out0[site] + w.L0[site] 
                    v.X[site] = -params.nu*v.edge_length + log(1-exp(-v.edge_length)) + v.A[site]
                    # compute out1: P(~D_v,v=-1)
                    if w.alpha[site] == 'z': # z-branch
                        v.out1[site] = log(1-exp(-v.edge_length*params.nu)) + v.A[site]
                        v.B[site] = min_llh
                        v.sis_tag[site] = 'z'
                    else:
                        alpha0 = w.alpha[site] if w.alpha[site] != '?' else u.sis_tag[site]
                        v.sis_tag[site] = alpha0
                        u_out_alpha = u.A[site] - params.nu*u.edge_length + log(1-exp(-u.edge_length)) + log(self.Q[site][alpha0]) # log P(~D_u,u=alpha_0)  
                        if u.sis_tag[site] is None or alpha0 == u.sis_tag[site]:
                            u_out_alpha = log(exp(u_out_alpha)+exp(u.B[site]))
                        v.B[site] = -params.nu*v.edge_length + w.L1[site]  + u_out_alpha
                        if w.alpha[site] == '?': # masked branch
                            v.X[site] = log(exp(v.X[site])+exp(u.X[site]+w.L1[site]-params.nu*v.edge_length))
                            v.out1[site] = log((1-exp(-v.edge_length*params.nu))*(exp(v.A[site])+exp(u.X[site]+w.L1[site])) + exp(u.out1[site]+1))
                        else:
                            v.X[site] = log(exp(v.X[site])+exp(v.B[site]))
                            v.out1[site] = log(1-exp(-v.edge_length*params.nu)) + log(exp(v.A[site]) + exp(w.L1[site]+u_out_alpha))
    def Estep_posterior(self,params):
        # assume binary tree
        # assume az_parition, Estep_in_llh, and Estep_out_llh have been performed 
        # so that all nodes have `alpha`, `L0`, `L1`, `out0`, and `out1` attribues
        # output: add the attributes `S0-4` (refer to the paper for definitions)
        # and `post0` and `post1` to each node where v.post0 = P(v=0|D) and v.post1 = P(v=-1|D) 
        full_llh = params.tree.root.L0
        for v in params.tree.traverse_preorder():
            v.post0 = [None]*self.numsites
            v.post1 = [None]*self.numsites
            v.S0 = [None]*self.numsites
            v.S1 = [None]*self.numsites
            v.S2 = [None]*self.numsites
            v.S3 = [None]*self.numsites
            v.S4 = [None]*self.numsites
            for site in range(self.numsites):
                # compute auxiliary values: v_in1 = log P(D_v|v=-1),v_in0 = log P(D_v|v=0)
                # v_in_alpha = log P(D_v|v=alpha0)
                v_in1 = 0 if v.alpha[site] == '?' else None
                #alpha_0 = v.alpha[site]
                if v.is_leaf():
                        c = self.charMtrx[v.label][site]
                        if c == 0:
                            v_in0 = log(1-params.phi)
                            #v_in_alpha = min_llh # z-branch
                        else:
                            v_in0 = log(params.phi) if c == '?' else min_llh    
                            #v_in_alpha = log(params.phi) if c == '?' else log(1-params.phi)
                else:    
                    v1,v2 = v.children
                    v_in0 = v1.L0[site] + v2.L0[site]                 
                # compute posterior
                v.post0[site] = v_in0 + v.out0[site] - full_llh[site]
                v.post1[site] = v_in1 + v.out1[site] - full_llh[site] if v_in1 is not None else min_llh                
                # compute S (note that all S values are NOT in log-scale)
                if v.alpha[site] == 'z': # masked branch
                    v.S0[site] = 1.0
                    v.S1[site] = v.S2[site] = v.S3[site] = v.S4[site] = 0.0
                elif v.is_root():
                    v.S0[site] = exp(v_in0 + (1.0+params.nu)*(-v.edge_length) - v.L0[site])
                    v.S2[site] = 0.0 if v.alpha[site] != '?' else (1.0-exp(-params.nu*v.edge_length))/exp(v.L0[site])
                    v.S1[site] = 1.0-v.S0[site]-v.S2[site]
                    v.S3[site] = v.S4[site] = 0.0
                else:
                    u = v.parent
                    v.S0[site] = exp(u.post0[site] + v_in0 + (1.0+params.nu)*(-v.edge_length) - v.L0[site])
                    if v.alpha[site] != '?':
                        v.S2[site] = 0.0
                        v.S4[site] = 0.0
                    else:
                        v.S2[site] = exp(u.post0[site]-v.L0[site])*(1.0-exp(-params.nu*v.edge_length))
                        v.S4[site] = (1.0-exp(u.post0[site])-exp(u.post1[site]))*(1.0-exp(-params.nu*v.edge_length))/exp(v.L1[site])
                    v.S1[site] = exp(u.post0[site]) - v.S0[site] - v.S2[site] 
                    v.S3[site] = 1.0-v.S0[site]-v.S1[site]-v.S2[site]-v.S4[site]

    def __Mstep_nu0__(self,params):
    # assume without checking that params.nu is negligible
    # should only be called by a function that has fixed_nu=eps
    # assume that Estep have been performed so that all nodes have S0-S4 attributes
    # output: optimize branch lengths
        for v in params.tree.traverse_preorder():
            S0 = sum(v.S0)
            S1 = sum(v.S1)
            p = S0/(S0+S1)
            dmax = -log(1/self.numsites)*2
            dmin = -log(1-1/self.numsites)/2
            if p <= exp(-dmax):
                d = dmax
            elif p >= exp(-dmin):
                d = dmin
            else:
                d = -log(p)         
            v.edge_length = d

    def __EM_nu0__(self,params,verbose=True):
    # assume without checking that params.nu is eps and params.phi is optimal
    # assume that az_partition has been performed
    # optimize other params while fixing nu to eps
    # caution: this function will modify params in place!
        pre_llh = self.lineage_llh(params)
        if verbose:
            print("Initial nllh: " + str(-pre_llh))
        em_iter = 1
        while 1:
            self.Estep_in_llh(params)
            self.Estep_out_llh(params)
            self.Estep_posterior(params)
            self.__Mstep_nu0__(params)
            curr_llh = self.lineage_llh(params)
            if verbose:
                print("EM iter: " + str(em_iter) + ". Current llh: " + str(-curr_llh))
            if curr_llh - pre_llh < em_eps:
                break
            pre_llh = curr_llh
            em_iter += 1
        return -curr_llh    

    def optimize(self,initials=20,fixed_phi=None,fixed_nu=None,verbose=True,max_trials=100):
        if fixed_nu is not None and fixed_nu <= eps:
            if fixed_phi is not None:
                phi_star = fixed_phi
            else:
                total = 0
                missing = 0
                for node in self.params.tree.traverse_leaves():
                    x = node.label
                    total += len(self.charMtrx[x])
                    missing += sum([y=='?' for y in self.charMtrx[x]])
                phi_star = missing/total    
            if verbose:
                print("Optimal phi: " + str(phi_star))
            results = []
            for rep in range(initials):
                if verbose:
                    print("Running EM with initial point " + str(rep+1))
                x0 = self.ini_all(fixed_phi=phi_star,fixed_nu=fixed_nu)
                self.x2params(x0,fixed_phi=phi_star,fixed_nu=fixed_nu)
                params = self.params
                self.az_partition(params)
                nllh = self.__EM_nu0__(params,verbose=verbose)
                results.append((nllh,params))
            results.sort()
            best_nllh,best_params = results[0]
            self.params = best_params
            return results[0][0]
        else:
            print("COMING SOON")     

def main(): 
    from sequence_lib import read_sequences
    from ml_log import wrapper_felsenstein as wf_log
    
    k = 50
    m = 10
    Q = []
    for i in range(k):
        q = {j+1:1/m for j in range(m)}
        q[0] = 0
        Q.append(q)
    T = "((a:0.0360971597765934,b:3.339535381892265)e:0.0360971597765934,(c:0.0360971597765934,d:3.339535381892265)f:0.0360971597765934)r:0.0001;"
    #T = "((a:1,b:1)e:1,(c:1,d:1)f:1)r:1;"
    #S = read_sequences("../tests/seqs_m10_k" + str(k) + ".txt",filetype="fasta")
    #msa = S[6]
    #msa['c'][0] = '?'

    msa = {'a':[0],'b':[0],'c':[1],'d':[1]}

    mySolver = EM_solver(msa,Q,T,phi=1e-10,nu=1e-10)
    mySolver.az_partition(mySolver.params)
    mySolver.Estep_in_llh(mySolver.params)
    #mySolver.lineage_llh(mySolver.params)
    mySolver.Estep_out_llh(mySolver.params)
    mySolver.Estep_posterior(mySolver.params)
    for node in mySolver.params.tree.traverse_postorder():
        #print(node.label,node.alpha,node.out0,node.out1)
        print(node.label,node.alpha,node.S0,node.S1,node.S2,node.S3,node.S4)
    #print(mySolver.optimize(initials=1,verbose=True,fixed_nu=None,fixed_phi=None))
    #print("phi", mySolver.params.phi,"nu", mySolver.params.nu)
    #print(mySolver.params.tree.newick()) 

if __name__ == "__main__":
    main()        
