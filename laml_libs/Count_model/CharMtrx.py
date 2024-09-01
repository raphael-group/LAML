class CharMtrx:
    def __init__(self,data_struct,alphabet):
    # K: the number of cassettes
    # J: the number of sites per cassette
    # data_struct: a mapping: cell_name -> (cassette -> cassette_state)
    # alphabet: an instance of class Alphabet; must have the same K and J
        self.K = alphabet.K
        self.J = alphabet.J
        if self.J == 1:
            data_struct_corrected = {}
            for c in data_struct:
                seq = []
                for x in data_struct[c]:
                    if not isinstance(x,tuple):
                        seq.append(tuple([x]))
                    else:
                        seq.append(x)    
                data_struct_corrected[c] = seq                
            self.data_struct = data_struct_corrected
        else:
            self.data_struct = data_struct    
        self.alphabet = alphabet
       
    def get_cell_names(self):
        return self.data_struct.keys()   

    def get(self,w,k):
        # w: a cell name/label
        # k: a cassette index
        return self.data_struct[w][k]

    def is_missing(self,w,k,missing_symbol='?'):    
        # w: a cell name/label
        # k: a cassette index
        # check if cassette k of cell w is missing
        x = self.get(w,k)
        for x_j in x:
            if x_j != missing_symbol:
                return False
        return True        
