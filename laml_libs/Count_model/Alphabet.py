class Alphabet():
    def __init__(self,K,J,data_struct):
        self.K = K
        self.J = J
        self.data_struct = data_struct # data_struct is a 3-way nested list
    
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
