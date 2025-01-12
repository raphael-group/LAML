from math import *

class Alphabet():
    """A class to represent the alphabet of a dynamic lineage tracing data. 
    Attributes:
        K: the number of cassettes
        J: the number of sites per cassette
        M: a list of K elements; element M[k] shows the size of the alphabet of cassette k
        data_struct: a three-way list, where data_struct[k][j] is a list of the alphabet (i.e. list of possible states) of site j of cassette k
    """
    def __init__(self,K,J,data_struct,silence_mechanism='convolve'):
        """Initialization of Alphabet
            K: the number of cassettes
            J: the number of sites per cassette
            data_struct: a three-way list, where data_struct[k][j] is a list of the alphabet (i.e. list of possible states) of site j of cassette k
            silence_mechanism: either 'convolve' or 'separated'
                'convolve': also is the silence mechanism introduced in the original LAML. This is site-independent silencing
                'separated': the silence mechanism unique to LAML-pro. This is cassette-independent silencing
            NOTE: we assume without checking that the data_struct has the consistent size as the input K and J
        """
        self.K = K
        self.J = J
        self.data_struct = data_struct
        
        self.M = [] # store the length of all cassette's alphabets
        for k in range(K):
            M_k = prod([len([x for x in A if x != -1]) for A in data_struct[k]]) # NOTE: this is the length of the cassette alphabet, not site alphabet
            self.M.append(M_k)
        
        self.silence_mechanism = silence_mechanism
        
        if silence_mechanism == 'separated':
            #self.J += 1
            for DS in self.data_struct:
                DS.append([0,-1]) # add one site to represent the silencing flag

    def get_site_alphabet(self,k,j): 
        """
            Get the alphabet of site j of cassette k
            return: a set of the possible site states
        """
        #print('k', k, 'j', j)
        return set([x for x in self.data_struct[k][j]])

    def get_M(self,k): 
        """ 
            Get the length of the alphabet of cassette k
            return: 
        """    
        return self.M[k]

    def get_cassette_alphabet(self,k):
        """
            Get the alphabet of cassette k
            return: a set of the possible cassette states, where a cassette state (i.e. allele) is represented by a tuple of self.J elements
        """
        a_list = self.data_struct[k]   
        def __get_alphabet(j):
            if j == 0:
                return [[x] for x in a_list[j]]
            else:
                prev = __get_alphabet(j-1)
                curr = []
                for x in prev:
                    for y in a_list[j]:        
                        curr.append(x+[y])
                return curr
        J = self.J if self.silence_mechanism == 'convolve' else self.J+1
        return set([tuple(x) for x in __get_alphabet(J-1)])
