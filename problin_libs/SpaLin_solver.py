from problin_libs.ML_solver import *

class SpaLin_solver(ML_solver):
    # at this stage, the tree topology and sig,a must be given. Only branch lengths
    # and other parameters can be optimized
    def __init__(self,charMtrx,Q,nwkTree,locations,sigma,nu=eps,phi=eps):
        super(SpaLin_solver,self).__init__(charMtrx,Q,nwkTree,nu=nu,phi=phi)
        self.given_locations = locations
        self.params.sigma = sigma
        self.inferred_locations = {}
        for x in self.given_locations:
            self.inferred_locations[x] = self.given_locations[x]
  
    def spatial_llh(self,locations):
        llh = 0
        for node in self.params.tree.traverse_preorder():
            if node.is_root() or not node.label in locations or not node.parent.label in locations:
                continue
            d = node.edge_length
            curr_sigma = self.params.sigma*sqrt(d)
            x,y = locations[node.label]
            x_par,y_par = locations[node.parent.label]
            llh -= (0.5*((x-x_par)/curr_sigma)**2 + log(curr_sigma))
            llh -= (0.5*((y-y_par)/curr_sigma)**2 + log(curr_sigma))
        return llh 

    def __llh__(self):
        return self.lineage_llh(self.params) + self.spatial_llh(self.inferred_locations)
    
    def ini_all(self,fixed_phi=None,fixed_nu=None):
        x_lin = self.ini_brlens() + [self.ini_nu(fixed_nu=fixed_nu),self.ini_phi(fixed_phi=fixed_phi)]
        x_spa = []
        for node in self.params.tree.traverse_postorder():
            if not node.label in self.given_locations:
                x_spa += [random(),random()]
        x_sigma = 22 # hard code for now        
        return x_lin + x_spa + [x_sigma]
    
    def x2params(self,x,fixed_nu=None,fixed_phi=None):
        self.x2brlen(x)
        self.x2nu(x,fixed_nu=fixed_nu)
        self.x2phi(x,fixed_phi=fixed_phi)
        i = self.num_edges + 2
        for node in self.params.tree.traverse_postorder():
            if not node.label in self.given_locations:
                self.inferred_locations[node.label] = (x[i],x[i+1])
                i += 2
        self.params.sigma = x[-1]        
               
    def bound_locations(self,lower=-np.inf,upper=np.inf):
        N = 2*len([node for node in self.params.tree.traverse_postorder() if not node.label in self.given_locations])    
        return [lower]*N,[upper]*N

    def bound_sigma(self):
        return (eps,np.inf)    

    def get_bound(self,keep_feasible=False,fixed_phi=None,fixed_nu=None):
        br_lower,br_upper = self.bound_brlen()  
        phi_lower,phi_upper = self.bound_phi(fixed_phi=fixed_phi)
        nu_lower,nu_upper = self.bound_nu(fixed_nu=fixed_nu)
        spa_lower,spa_upper = self.bound_locations()
        sigma_lower,sigma_upper = self.bound_sigma()
        bounds = optimize.Bounds(br_lower+[nu_lower,phi_lower]+spa_lower+[sigma_lower],br_upper+[nu_upper,phi_upper]+spa_upper+[sigma_upper],keep_feasible=keep_feasible)
        return bounds
