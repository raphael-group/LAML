from math import *

class Alphabet():
    def __init__(self,K,J,data_struct):
        self.K = K
        self.J = J
        self.M = [] # store the length of all cassette's alphabets
        for k in range(K):
            M_k = prod([len(A) for A in data_struct[k]])
            self.M.append(M_k)
        self.data_struct = data_struct # data_struct is a 3-way nested list

    def get_site_alphabet(self,k,j): # get the alphabet of site j of cassette k
        return set([x for x in self.data_struct[k][j]])

    def get_M(self,k): # get the length of the alphabet of cassette k
        return self.M[k]

    def get_cassette_alphabet(self,k):
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
        return [tuple(x) for x in __get_alphabet(self.J-1)]
