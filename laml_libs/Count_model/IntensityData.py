class IntensityData:
    """A class to represent the intensity data (i.e. the observed data which has associated count/intensity for every allele)
    Attributes:
        K: the number of cassettes
        J: the number of sites per cassette (n the case of baseMEMOIR, consider J=1)
        data_struct: a mapping: cell_name -> (cassette -> feature_vector). The data type is "dictionary of lists of dictionaries"
        alphabet: an instance of class Alphabet
    NOTE: We will use "cassette_state" and "allele" interchangeably throughout
    """
    def __init__(self,data_struct,alphabet,max_allele_per_cassette=None):
        """Initialization of IntensityData by the data_struct, alphabet, and max_allele_per_cassette
            data_struct: a mapping: cell_name -> (cassette -> feature_vector)). The data type is "dictionary of lists of dictionaries"
            alphabet: an instance of class Alphabet
            max_allele_per_cassette: is either None or a positive integer. 
                If None: use all available alleles
                If not none: for every cell at every cassette, only use maximum `max_allele_per_cassette` top alleles 
        NOTE: we assume without checking that the alphabet and data_struct have the same K and J
        """
        self.K = alphabet.K # the number of cassettes
        self.J = alphabet.J # the number of sites per cassette
        self.alphabet = alphabet
        self.data_struct = data_struct
        
        self.cassette_state_lists = []
        for k in range(self.K):
            self.cassette_state_lists.append(self.alphabet.get_cassette_alphabet(k))

    def get_one_count(self,w,k,x):
        """
            Get the UMI count of cell w at cassette k of allele x
                w: a cell name/label
                k: a cassette index
                x: a state of the specified cassette; must be a tuple of length self.J
            return: an integer    
        """
        assert len(x) == self.J
        return self.data_struct[w][k][x]
    
    def get(self,w,k): # get all counts
        """
            Get the UMI count of cell w at cassette k of all alleles
                w: a cell name/label
                k: a cassette index
            return: a dictionary cassette_state (i.e. allele) -> count    
        """
        return self.data_struct[w][k]

    def is_missing(self,w,k):    
        """ Check if cassette k of cell w is missing (i.e. has no UMI count)
                w: a cell name/label
                k: a cassette index
            return: True/False    
        """            
        c = self.data_struct[w][k]
        for x in c:
            if c[x] != 0:
                return False
        return True        
