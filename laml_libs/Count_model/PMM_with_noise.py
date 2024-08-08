from .PMM_base import *

class PMMN_model(PMM_base):
    # PMMN = PMM with noise = probabilistic mixed-type missing with noise
    def __init__(self,treeList,data,prior,mu,nu,phi,eta):
        params = Param(['mu','nu','phi','eta'],[mu,nu,phi,eta],[0,0,0,0],[DEFAULT_max_mu,DEFAULT_max_nu,1,1])
        super(PMM_base,self).__init__(treeList,data,prior,params) #####*****#####

    def Psi(self,c_node,k,j,alpha,beta):
        return self._PMM_base__Psi(c_node,k,j,alpha,beta) # simply inherit from PMM_base class

    def Gamma(self,k,x,c):
        # Layer 2 transition probabilities  
        # override the PMM_base class
        ##### TODO #####
        return 1
