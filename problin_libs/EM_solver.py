from problin_libs.ML_solver import *
from math import exp,log
import cvxpy as cp
from problin_libs import min_llh, conv_eps, eps

def log_sum_exp(numlist):
    # using log-trick to compute log(sum(exp(x) for x in numlist))
    # mitigate the problem of underflow
    maxx = max(numlist)
    result = maxx + log(sum([exp(x-maxx) for x in numlist]))
    return result

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
                    l0_z = [l0 + (nu+1)*(-node.edge_length)]
                    l0_alpha = [l1 + log(1-p)+log(q) + nu*(-node.edge_length)] if node.alpha[site] != 'z' else []
                    l0_masked = [log(1-p**nu)] if (node.alpha[site] == '?' and nu > 0) else []
                    node.L0[site] = log_sum_exp(l0_z + l0_alpha + l0_masked)
                    node.L1[site] = min_llh if node.alpha[site] == 'z' else log((exp(l1+nu*(-node.edge_length)) + (1-p**nu)*int(node.alpha[site]=="?")))
    
    def lineage_llh(self,params):
        # override the function of the base class
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
                # Auxiliary components
                v.A = [0]*self.numsites
                v.X = [-params.nu*v.edge_length + log(1-exp(-v.edge_length))]*self.numsites
                v.out_alpha = [{} for i in range(self.numsites)]
                # Main components    
                v.out0 = [-(1+params.nu)*v.edge_length]*self.numsites
                v.out1 = [log(1-exp(-v.edge_length*params.nu))]*self.numsites if params.nu>0 else [min_llh]*self.numsites
            else:
                u = v.parent
                # get the sister
                for x in u.children:
                    if x is not v:
                        w = x
                        break
                # Auxiliary components
                v.A = [None]*self.numsites
                v.X = [None]*self.numsites
                v.out_alpha = [{} for i in range(self.numsites)]
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
                        v.out1[site] = log(1-exp(-v.edge_length*params.nu)) + v.A[site] if params.nu > 0 else min_llh
                    elif w.alpha[site] == '?': # masked branch
                        v.X[site] = log_sum_exp([v.X[site],u.X[site]+w.L1[site]-params.nu*v.edge_length])
                        p = 1-exp(-v.edge_length*params.nu) # if nu=0 then p=0
                        pl = log(p) if p > 0 else min_llh
                        v.out1[site] = log_sum_exp([pl+v.A[site],pl+u.X[site]+w.L1[site],u.out1[site]])
                    else:
                        alpha0 = w.alpha[site] 
                        if alpha0 not in u.out_alpha[site]:
                            self.__out_alpha_up__(u,site,alpha0,params)
                        B = u.out_alpha[site][alpha0] + params.nu*(-v.edge_length) + w.L1[site]  
                        C = v.A[site] - params.nu*v.edge_length + log(1-exp(-v.edge_length)) + log(self.Q[site][alpha0])
                        v.out_alpha[site][alpha0] = log_sum_exp([B,C])
                        v.X[site] = log_sum_exp([v.X[site],B])
                        if params.nu > 0:
                            v.out1[site] = log(1-exp(-v.edge_length*params.nu)) + log_sum_exp([v.A[site],w.L1[site]+u.out_alpha[site][alpha0]])
                        else:
                            v.out1[site] = min_llh    

    def __out_alpha_up__(self,node,site,alpha0,params):
        # auxiliary function, shoudn't be called outside
        v = node
        path = []   
        # going up to until reaching either a z-branch or the root
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
        # going down to compute all the out_alpha along the path
        while path:
            v = path.pop()
            u = v.parent
            for x in u.children:
                if x is not v:
                    w = x
            #B = u.out_alpha[site][alpha0] + params.nu*(-u.edge_length) + w.L1[site]  
            B = u.out_alpha[site][alpha0] + params.nu*(-v.edge_length) + w.L1[site]  
            C = v.A[site] - params.nu*v.edge_length + log(1-exp(-v.edge_length)) + log(self.Q[site][alpha0])
            v.out_alpha[site][alpha0] = log_sum_exp([B,C])
            #v.out_alpha[site][alpha0] = log(exp(C) + exp(B))

    def Estep_posterior(self,params):
        # assume binary tree
        # assume az_parition, Estep_in_llh, and Estep_out_llh have been performed 
        # so that all nodes have `alpha`, `L0`, `L1`, `out0`, and `out1` attribues
        # output: add the attributes `S0-4` (refer to the paper for definitions; all S are NOT stored in log-scale)
        # and `post0` and `post1` to each node where v.post0 = log P(v=0|D) and v.post1 = log P(v=-1|D) 
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
                # compute auxiliary values: v_in1 = log P(D_v|v=-1),v_in0 = log P(D_v|v=0), v_in_alpha = log P(D_v|v=alpha0)
                v_in1 = 0 if v.alpha[site] == '?' else None
                if v.is_leaf():
                        c = self.charMtrx[v.label][site]
                        if c == 0:
                            v_in0 = log(1-params.phi)
                        else:
                            v_in0 = log(params.phi) if (c == '?' and params.phi > 0) else min_llh 
                else:    
                    v1,v2 = v.children
                    v_in0 = v1.L0[site] + v2.L0[site]                 
                # compute posterior
                v.post0[site] = v_in0 + v.out0[site] - full_llh[site]
                v.post1[site] = v_in1 + v.out1[site] - full_llh[site] if v_in1 is not None else min_llh               
                # compute S (note that all S values are NOT in log-scale)
                if v.alpha[site] == 'z': # z-branch
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
                    else: # masked branch
                        v.S2[site] = exp(u.post0[site]-v.L0[site])*(1.0-exp(-params.nu*v.edge_length))
                        v.S4[site] = (1.0-exp(u.post0[site])-exp(u.post1[site]))*(1.0-exp(-params.nu*v.edge_length))/exp(v.L1[site])
                    v.S1[site] = exp(u.post0[site]) - v.S0[site] - v.S2[site] 
                    v.S3[site] = 1.0-v.S0[site]-v.S1[site]-v.S2[site]-v.S4[site]

    def Estep(self,params):
        self.Estep_in_llh(params)
        self.Estep_out_llh(params)
        self.Estep_posterior(params)

    def __Mstep_nu0__(self,params):
    # assume without checking that params.nu is 0
    # should only be called by a function with this same assumption
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
        print("Initial nllh: " + str(-pre_llh))
        em_iter = 1
        while 1:
            self.Estep(params)
            self.__Mstep_nu0__(params)
            curr_llh = self.lineage_llh(params)
            if verbose:
                print("EM iter: " + str(em_iter) + ". Current nllh: " + str(-curr_llh))
            if curr_llh - pre_llh < conv_eps:
                break
            pre_llh = curr_llh
            em_iter += 1
        return -curr_llh    

    def Mstep(self,params,optimize_phi=True,optimize_nu=True,verbose=True,eps_nu=1e-5,eps_s=1e-6):
    # assume that Estep have been performed so that all nodes have S0-S4 attributes
    # output: optimize all parameters: branch lengths, phi, and nu
        try:
            # optimize phi
            if not optimize_phi:
                if verbose:
                    print("Fixing phi to " + str(params.phi))    
                phi_star = params.phi
            else:       
                if verbose:
                    print("Optimizing phi")
                R = []
                R_tilde = []
                for i,v in enumerate(params.tree.traverse_leaves()):
                    R.append(sum([x != '?' for x in self.charMtrx[v.label]]))
                    R_tilde.append(sum([1-exp(p) for (j,p) in enumerate(v.post1) if self.charMtrx[v.label][j] == '?'])) 
                phi_star = sum(R_tilde)/(sum(R)+sum(R_tilde))
                if abs(phi_star) < 1/(self.numsites*len(self.charMtrx)):
                    phi_star = 0
            # optimize nu and all branch lengths
            N = len(list(params.tree.traverse_preorder())) # the number of branches
            S0 = np.zeros(N)
            S1 = np.zeros(N)
            S2 = np.zeros(N)
            S3 = np.zeros(N)
            S4 = np.zeros(N)
            for i,v in enumerate(params.tree.traverse_preorder()):
                s = [sum(v.S0),sum(v.S1),sum(v.S2),sum(v.S3),sum(v.S4)]
                s = [x if abs(x) > eps_s else 0 for x in s]
                S0[i],S1[i],S2[i],S3[i],S4[i] = s
            def __optimize_brlen__(nu): # nu is a single number
                dmax = -log(1/self.numsites)*2
                dmin = -log(1-1/self.numsites)/2
                D = np.zeros(N)
                if nu <= eps_nu:
                    for i in range(N):
                        p = S0[i]/(S0[i]+S1[i])
                        if p <= exp(-dmax):
                            d = dmax
                        elif p >= exp(-dmin):
                            d = dmin
                        else:
                            d = -log(p)        
                        D[i] = d     
                    return D    
                else:
                    var_d = cp.Variable(N) # the branch length variables
                    C0 = -(nu+1)*S0.T @ var_d
                    C1 = -nu*S1.T @ var_d + S1.T @ cp.log(1-cp.exp(-var_d))
                    C2 = S2.T @ cp.log(1-cp.exp(-nu*var_d)) 
                    C3 = -nu*S3.T @ var_d
                    C4 = S4.T @ cp.log(1-cp.exp(-nu*var_d))

                    objective = cp.Maximize(C0+C1+C2+C3+C4)
                    constraints = [np.zeros(N)+dmin <= var_d, var_d <= np.zeros(N)+dmax]
                    prob = cp.Problem(objective,constraints)
                    prob.solve(verbose=False)
                    return var_d.value
            def __optimize_nu__(d): # d is a vector of all branch lengths
                var_nu = cp.Variable(1,nonneg=True) # the nu variable
                C0 = -(var_nu+1)*S0.T @ d
                C1 = -var_nu*S1.T @ d + S1.T @ cp.log(1-cp.exp(-d))
                C2 = S2.T @ cp.log(1-cp.exp(-var_nu*d)) if sum(S2) > 0 else 0
                C3 = -var_nu*S3.T @ d
                C4 = S4.T @ cp.log(1-cp.exp(-var_nu*d)) if sum(S4) > 0 else 0

                objective = cp.Maximize(C0+C1+C2+C3+C4)
                #constraints = [var_nu >= 1e-7]
                #prob = cp.Problem(objective,constraints)
                prob = cp.Problem(objective)
                prob.solve(verbose=False)
                return var_nu.value[0]

            initial_nu = params.nu
            if verbose:
                print("Optimizing branch lengths. Current phi: " + str(phi_star) + ". Current nu:" + str(initial_nu))
            d_star = __optimize_brlen__(initial_nu)
            if not optimize_nu:
                if verbose:
                    print("Fixing nu to " + str(params.nu))
                nu_star = params.nu
            else:    
                if verbose:
                    print("Optimizing nu")
                nu_star = __optimize_nu__(d_star) 
            # place the optimal value back to params
            params.phi = phi_star
            params.nu = nu_star
            for i,node in enumerate(params.tree.traverse_preorder()):
                node.edge_length = d_star[i]
            return True    
        except:
            return False        
    
    def EM_optimization(self,params,verbose=True,optimize_phi=True,optimize_nu=True):
    # assume that az_partition has been performed
    # optimize all parameters: branch lengths, phi, and nu
    # if optimize_phi is False, it is fixed to the original value in params.phi
    # the same for optimize_nu
    # caution: this function will modify params in place!
        pre_llh = self.lineage_llh(params)
        print("Initial phi: " + str(params.phi) + ". Initial nu: " + str(params.nu) + ". Initial nllh: " + str(-pre_llh))
        em_iter = 1
        while 1:
            if verbose:
                print("Starting EM iter: " + str(em_iter))
                print("Estep")
            self.Estep(params)
            if verbose:
                print("Mstep")
            if not self.Mstep(params,optimize_phi=optimize_phi,optimize_nu=optimize_nu,verbose=verbose):
                print("Fatal error: failed to optimize parameters in Mstep!")
                return None
            curr_llh = self.lineage_llh(params)
            if verbose:
                print("Finished EM iter: " + str(em_iter) + ". Current nllh: " + str(-curr_llh))
            if abs((curr_llh - pre_llh)/pre_llh) < conv_eps:
            #print(curr_llh-pre_llh,conv_eps)
            #if (curr_llh - pre_llh) < conv_eps:
                break
            pre_llh = curr_llh
            em_iter += 1
        return -curr_llh    
    
    def optimize(self,initials=20,fixed_phi=None,fixed_nu=None,verbose=True,max_trials=100,random_seeds=None):
    # override the same function of the base class
    # random_seeds can either be a single number or a list of intergers where len(random_seeds) = initials
        results = []
        all_failed = True
        all_trials = 0
        if random_seeds is None:
            rseeds = [int(random()*10000) for i in range(initials)]
        elif type(random_seeds) == int:
            print("Global random seed: " + str(random_seeds))
            seed(a=random_seeds)
            rseeds = [int(random()*10000) for i in range(initials)]
        elif type(random_seeds) == list:
            if len(random_seeds) < initials:
                print("Fatal: the number of random seeds is smaller than the number of initials!")
                return None
            elif len(random_seeds) > initials:
                print("Warning: the number of random seeds is larger than the number of initials. Ignoring the last " + str(len(random_seeds)-initials) + " seeds")
            rseeds = random_seeds[:initials]    
        else:
            print("Fatal: incorrect random_seeds type provided")        
            return None
        while all_failed and all_trials < max_trials:
            for rep in range(initials):
                randseed = rseeds[rep]
                print("Initial point " + str(rep+1) + ". Random seed: " + str(randseed))
                seed(a=randseed)
                print("Running EM with initial point " + str(rep+1))
                x0 = self.ini_all(fixed_phi=fixed_phi,fixed_nu=fixed_nu)
                self.x2params(x0,fixed_phi=fixed_phi,fixed_nu=fixed_nu)
                params = self.params
                self.az_partition(params)
                nllh = self.EM_optimization(params,verbose=verbose,optimize_phi=(fixed_phi is None),optimize_nu=(fixed_nu is None))
                if nllh is not None:
                    all_failed = False
                    print("Optimal point found for initial point " + str(rep+1))
                    print("Optimal phi: " + str(params.phi))
                    print("Optimal nu: " + str(params.nu))
                    print("Optimal tree: " + params.tree.newick())
                    print("Optimal nllh: " + str(nllh))
                    results.append((nllh,params))
                else:
                    print("Fatal: failed to optimize using initial point " + str(rep+1))    
            all_trials += initials    
        results.sort()
        best_nllh,best_params = results[0]
        self.params = best_params
        return results[0][0]

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
    #T = "((a:0.0360971597765934,b:3.339535381892265)e:0.0360971597765934,(c:0.0360971597765934,d:3.339535381892265)f:0.0360971597765934)r:0.0001;"
    #T = "((a:1,b:1)e:1,(c:1,d:1)f:1)r:1;"
    T = "(f:1,(a:1,b:1)e:1)r:1;"
    
    #S = read_sequences("../tests/seqs_m10_k" + str(k) + ".txt",filetype="fasta")
    #msa = S[6]
    #msa['c'][0] = '?'

    msa = {'f':[0],'a':['?'],'b':['?'],'c':[1],'d':[2]}

    mySolver = EM_solver(msa,Q,T,phi=1e-10,nu=0.1)
    mySolver.az_partition(mySolver.params)
    mySolver.Estep_in_llh(mySolver.params)
    print(mySolver.lineage_llh(mySolver.params))
    '''mySolver.Estep_out_llh(mySolver.params)
    mySolver.Estep_posterior(mySolver.params)
    for node in mySolver.params.tree.traverse_postorder():
        print(node.label,node.alpha,node.out0,node.out1)
        #print(node.label,node.alpha,node.S0,node.S1,node.S2,node.S3,node.S4)'''
    #print(mySolver.optimize(initials=1,verbose=True,fixed_nu=None,fixed_phi=None))
    #print("phi", mySolver.params.phi,"nu", mySolver.params.nu)
    #print(mySolver.params.tree.newick()) 

if __name__ == "__main__":
    main()        
