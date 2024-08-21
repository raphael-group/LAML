class AlleleTable:
    def __init__(self,K,J,data_struct,alphabet):
    # K: the number of cassettes
    # J: the number of sites per cassette
    # data_struct: a mapping: cell_name -> (cassette -> (cassette_state -> count))
    # alphabet: an instance of class Alphabet; must have the same K and J
        self.K = K
        self.J = J
        self.data_struct = data_struct
        self.alphabet = alphabet
        
    def get(self,w,k,x):
        # w: a cell name/label
        # k: a cassette index
        # x: a state of the specified cassette; must be a tuple of length J
        assert len(x) == self.J
        return self.data_struct[w][k][x]
    
    def get_all_counts(self,w,k):
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
