import heapq

class AlleleTable:
    def __init__(self,data_struct,alphabet,max_allele_per_cassette=None):
    # K: the number of cassettes
    # J: the number of sites per cassette
    # data_struct: a mapping: cell_name -> (cassette -> (cassette_state -> count))
    # alphabet: an instance of class Alphabet
    # assume without checking: the alphabet and data_struct have the same K and J
        self.K = alphabet.K # the number of cassettes
        self.J = alphabet.J # the number of sites per cassette
        self.alphabet = alphabet
        
        if self.J == 1: # may need to correct the cassette states to be tuples, not single values
            data_struct_corrected = {}
            for cell in data_struct:
                countList = []
                for c in data_struct[cell]:
                    c_corrected = {}
                    for x in c:
                        if not isinstance(x,tuple):
                            c_corrected[tuple([x])] = c[x]
                        else:
                            c_corrected[x] = c[x]
                    countList.append(c_corrected)        
                data_struct_corrected[cell] = countList 
            self.data_struct = data_struct_corrected
        else:
            self.data_struct = data_struct    
        
        if max_allele_per_cassette is None:
            self.cassette_state_lists = []
            for k in range(self.K):
                self.cassette_state_lists.append(self.alphabet.get_cassette_alphabet(k))
        else:
            root_state = tuple([0]*self.J)
            silenced_state = tuple([-1]*self.J)
            self.cassette_state_lists = [set([root_state,silenced_state]) for _ in range(self.K)]
            for k in range(self.K):                        
                for cell in self.data_struct:
                    top_alleles = set(heapq.nlargest(max_allele_per_cassette,self.data_struct[cell][k], key=self.data_struct[cell][k].get))
                    self.cassette_state_lists[k] = self.cassette_state_lists[k].union(top_alleles)                

    def get_one_count(self,w,k,x):
    # w: a cell name/label
    # k: a cassette index
    # x: a state of the specified cassette; must be a tuple of length J
        assert len(x) == self.J
        return self.data_struct[w][k][x]
    
    def get(self,w,k): # get all counts
    # w: a cell name/label
    # k: a cassette index
        return self.data_struct[w][k]

    def is_missing(self,w,k):    
    # w: a cell name/label
    # k: a cassette index
    # check if cassette k of cell w is missing (i.e. has no UMI count)
        c = self.data_struct[w][k]
        for x in c:
            if c[x] != 0:
                return False
        return True        
