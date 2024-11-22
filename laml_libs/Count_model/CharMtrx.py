class CharMtrx:
    """A class to represent the character matrix
    Attributes:
        K: the number of cassettes
        J: the number of sites per cassette
        data_struct: a mapping: cell_name -> (cassette -> cassette_state). The data type is "dictionary of lists of cassette_states"
        alphabet: an instance of class Alphabet
    NOTE: We will use "cassette_state" and "allele" interchangeably throughout
    """
    def __init__(self,data_struct,alphabet):
        """Initialization of CharMtrx by the data_struct and alphabet
            data_struct: a mapping: cell_name -> (cassette -> cassette_state). The data type is "dictionary of lists of dictionaries"
            alphabet: an instance of class Alphabet
        NOTE: we assume without checking that the alphabet and data_struct have the same K and J
        """
        self.K = alphabet.K
        self.J = alphabet.J
        self.alphabet = alphabet
        self.silence_mechanism = alphabet.silence_mechanism 
        if self.J == 1: # when J=1, we may need to correct the cassette states to be tuples, not single values
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
        self.cassette_state_lists = []
        for k in range(self.K):
            self.cassette_state_lists.append(self.alphabet.get_cassette_alphabet(k))
       
    def get_cell_names(self):
        """ Get all cell names """
        return self.data_struct.keys()   

    def get(self,w,k):
        """
            Get the cassette state (i.e. allele) of cell w at cassette k
                w: a cell name/label
                k: a cassette index
            return: an allele, which is represented by a tuple of self.J elements    
        """
        return self.data_struct[w][k]

    def is_missing(self,w,k,missing_symbol='?'):    
        """Check if cassette k of cell w is missing
            w: a cell name/label
            k: a cassette index
        """
        x = self.get(w,k)
        for x_j in x:
            if x_j != missing_symbol:
                return False
        return True        
