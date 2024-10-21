from .PMM_base import *
import numpy as np

class PMMN_mat_model(PMM_model):
    # PMMN = PMM with noise = probabilistic mixed-type missing with noise
    def __init__(self,treeList,data,prior,**params): 
        # **params must have 'mu','nu','phi','eta_mat'

        self.num_cassettes = K
        self.site_per_cassette = J

        # not possible to observe incorrectly 
        DEFAULT_min_eta = dict()
        # uniform
        DEFAULT_max_eta = dict()

        for k in range(K):
            # get_cassette_alphabet returns a 3-way nested list (always same order)
            allele_list = self.data['alleleTable'].alphabet.get_cassette_alphabet(k)
            m_k = len(allele_list)

            tmp_k_min = np.eye(m_k)
            uniform_prob = 1/float(m_k)
            tmp_k_max = np.full((m_k, m_k), uniform_prob, dtype=float)

            DEFAULT_min_eta[k] = dict()
            DEFAULT_max_eta[k] = dict()
            
            for i, from_state in enumerate(allele_list):
                for j, to_state in enumerate(allele_list):
                    if i not in DEFAULT_min_eta[k]:
                        DEFAULT_min_eta[k][i] = dict()
                    if i not in DEFAULT_max_eta[k]:
                        DEFAULT_max_eta[k][i] = dict()
                    DEFAULT_min_eta[k][from_state][to_state] = tmp_k_min[i][j]
                    DEFAULT_max_eta[k][from_state][to_state] = tmp_k_max[i][j]

        params = Param(['mu','nu','phi','eta_mat'],[params['mu'],params['nu'],params['phi'],params['eta_mat']],[0,0,0,DEFAULT_min_eta],[DEFAULT_max_mu,DEFAULT_max_nu,1,DEFAULT_max_eta])
        super(PMM_model,self).__init__(treeList,data,prior,params) #####*****#####

    def Psi(self,c_node,k,j,alpha,beta):
        #return self.Psi(c_node,k,j,alpha,beta) # simply inherit from PMM_base class
        p = super().Psi(c_node,k,j,alpha,beta)
        return p

    def Gamma(self,k,x,c,node_label):
        # Layer 2 transition probabilities  
        
        # list of length k confusion matrices
        # get confusion matrix corresponding to cassette k
        # get probability of observing state x as state c
        eta = self.params.get_value('eta_mat')
        
        cassette_alphabet = self.data['alleleTable'].alphabet.get_cassette_alphabet(k)
        from_cassette_idx = cassette_alphabet.index(x)


        max_count = 0
        argmax_cassette_states = set()
        for y in c:
            if c[y] > max_count:
                max_count = c[y]
        for y in c:
            if c[y] == max_count:
                argmax_cassette_states.add(y)

        # TODO: to_cassette_idx should be the max
        if x in argmax_cassette_states:
            #argmax_cassette_state = argmax_cassette_states # how to handle multiple?
            to_cassette_idx = cassette_alphabet.index(x)
        else:
            # get the cassette_state with the maximum observation probability? 
            to_cassette_idx = cassette_alphabet.index()

        
        p = eta[k][from_cassette_idx][to_cassette_idx]

        return p

