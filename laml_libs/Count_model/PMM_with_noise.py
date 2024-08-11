from .PMM_base import *

class PMMN_model(PMM_model):
    # PMMN = PMM with noise = probabilistic mixed-type missing with noise
    def __init__(self,treeList,data,prior,mu,nu,phi,eta):
        params = Param(['mu','nu','phi','eta'],[mu,nu,phi,eta],[0,0,0,0],[DEFAULT_max_mu,DEFAULT_max_nu,1,1])
        super(PMM_model,self).__init__(treeList,data,prior,params) #####*****#####

    def Psi(self,c_node,k,j,alpha,beta):
        #return self.Psi(c_node,k,j,alpha,beta) # simply inherit from PMM_base class
        p = super().Psi(c_node,k,j,alpha,beta)
        return p

    def Gamma(self,k,x,c):
        # Layer 2 transition probabilities  
        # override the PMM_base class
        K = self.data['alleleTable'].K
        silenced_state = tuple([-1]*K)
        x_is_silenced = (x == silenced_state)
        phi = self.params.get_value('phi')
        # overall confusion rate
        eta = self.params.get_value('eta')

        cassette_alphabet_size = len(c.keys() - set([(-1,)]))
        max_count = 0
        argmax_cassette_states = set()
        for y in c:
            if c[y] > max_count:
                max_count = c[y]
        
        for y in c:
            if c[y] == max_count:
                argmax_cassette_states.add(y)

        if max_count == 0:
            p = 1 if x_is_silenced else phi
        else:
            # if x doesn't match the argmax cassette states 
            # (of which there may be more than one due to ties)
            if x not in argmax_cassette_states:
                p = eta
            else:
                # subtract argmax state, and one for -1 silenced state already subtracted above
                # and subtract one more for 
                p = 1 - phi - (eta * (cassette_alphabet_size-len(argmax_cassette_states)))
                #p = (1 - phi - (eta * (cassette_alphabet_size-len(argmax_cassette_states)-1)))/len(argmax_cassette_states)

        #print("gamma:", p)
        return p
