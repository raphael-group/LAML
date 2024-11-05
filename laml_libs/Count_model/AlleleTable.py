class AlleleTable:
    def __init__(self,data_struct,alphabet):
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
            if c[x] != 0 and c[x] is not None:
                return False
        return True        
