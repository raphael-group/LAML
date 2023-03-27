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
            # print(node.label, node.edge_length)
            p = exp(-node.edge_length)
            node.L0 = [0]*self.numsites # L0 and L1 are stored in log-scale
            node.L1 = [0]*self.numsites
            for site in range(self.numsites):   
                q = self.Q[site][node.alpha[site]] if node.alpha[site] not in ['?','z'] else 1.0
                if node.is_leaf():
                    if node.alpha[site] == "?":         
                        #masked_llh = log(1-(1-phi)*p**nu)
                        masked_llh = log(1-(1-phi)*p**nu) if self.charMtrx[node.label][site] == '?' else log(1-p**nu)
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
                    # note: l0_z, l0_alpha, and l0_masked are lists    
                    l0_z = [l0 + (nu+1)*(-node.edge_length)]
                    l0_alpha = [l1 + log(1-p)+log(q) + nu*(-node.edge_length)] if node.alpha[site] != 'z' else []
                    l0_masked = [log(1-p**nu)] if (node.alpha[site] == '?' and nu > 0) else []
                    node.L0[site] = log_sum_exp(l0_z + l0_alpha + l0_masked)
                    if node.alpha[site] == 'z':
                        node.L1[site] = min_llh
                    elif node.alpha[site] != '?' or nu == 0:
                        node.L1[site] = l1 + nu*(-node.edge_length) 
                    else:
                        node.L1[site] = log_sum_exp([l1 + nu*(-node.edge_length), log(1-p**nu)])

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
                        v.out1[site] = u.out1[site] if params.nu == 0 else log_sum_exp([pl+v.A[site],pl+u.X[site]+w.L1[site],u.out1[site]])
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
                v.post0[site] = v_in0 + v.out0[site] - full_llh[site] if v_in0 is not None else min_llh
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
                    v.S3[site] = 1.0-v.S0[site]-v.S1[site]-exp(v.post1[site])

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
        if verbose:
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

    '''def optimize(self,initials=20,fixed_phi=None,fixed_nu=None,verbose=True,max_trials=100):
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
                self.x2params(x0,fixed_phi=phi_star,fixed_nu=fixed_nu)'''

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
                #s = [x if abs(x) > eps_s else eps_s for x in s]
                s = [max(eps_s,x) for x in s]
                s = [x/sum(s)*self.numsites for x in s]
                S0[i],S1[i],S2[i],S3[i],S4[i] = s
    
            def __optimize_brlen__(nu): # nu is a single number
                #print("[__optimize_brlen__]")
                zerocount = sum([self.charMtrx[e].count(0) for e in self.charMtrx]) 
                totalcount = self.numsites * len(self.charMtrx)
                zeroprop = zerocount/totalcount
                #print("zerocount", zerocount, "totalcount", totalcount, "zeroprop", zeroprop)
                #dmax = -log(1/self.numsites)*2
                dmax = -log(zeroprop)
                #print("dmax:", dmax)
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
                    var_d = cp.Variable(N,nonneg=True) # the branch length variables
                    C0 = -(nu+1)*S0.T @ var_d
                    C1 = -nu*S1.T @ var_d + S1.T @ cp.log(1-cp.exp(-var_d))
                    C2 = S2.T @ cp.log(1-cp.exp(-nu*var_d)) if sum(S2) > 0 else 0 
                    C3 = -nu*S3.T @ var_d
                    C4 = S4.T @ cp.log(1-cp.exp(-nu*var_d)) if sum(S4) > 0 else 0

                    objective = cp.Maximize(C0+C1+C2+C3+C4)
                    constraints = [np.zeros(N)+dmin <= var_d, var_d <= np.zeros(N)+dmax]
                    prob = cp.Problem(objective,constraints)
                    prob.solve(verbose=False,solver=cp.MOSEK)
                    return var_d.value
            
            def __optimize_nu__(d): # d is a vector of all branch lengths
                var_nu = cp.Variable(1,nonneg=True) # the nu variable
                C0 = -(var_nu+1)*S0.T @ d
                C1 = -var_nu*S1.T @ d + S1.T @ cp.log(1-cp.exp(-d))
                C2 = S2.T @ cp.log(1-cp.exp(-var_nu*d)) if sum(S2) > 0 else 0
                C3 = -var_nu*S3.T @ d
                C4 = S4.T @ cp.log(1-cp.exp(-var_nu*d)) if sum(S4) > 0 else 0
                objective = cp.Maximize(C0+C1+C2+C3+C4)
                prob = cp.Problem(objective)
                prob.solve(verbose=False,solver=cp.MOSEK)
                return var_nu.value[0]

            nIters = 1
            nu_star = params.nu
            for r in range(nIters):
                if verbose:
                    print("Optimizing branch lengths. Current phi: " + str(phi_star) + ". Current nu:" + str(nu_star))
                d_star = __optimize_brlen__(nu_star)
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
                break
            pre_llh = curr_llh
            em_iter += 1
        return -curr_llh, em_iter    
   
    def optimize_one(self,randseed,fixed_phi=None,fixed_nu=None,verbose=True):
        # optimize using a specific initial point identified by the input randseed
        seed(a=randseed)
        x0 = self.ini_all(fixed_phi=fixed_phi,fixed_nu=fixed_nu)
        self.x2params(x0,fixed_phi=fixed_phi,fixed_nu=fixed_nu)
        params = self.params
        self.az_partition(params)
        nllh, em_iter = self.EM_optimization(params,verbose=verbose,optimize_phi=(fixed_phi is None),optimize_nu=(fixed_nu is None))
        print("EM finished after " + str(em_iter) + " iterations.")
        return nllh,params 
