class DLT_parser: # DLT: dynamic lineage tracing
    def charMtrx_to_json(self,charMtrxFile,outputFile,<other_arguments>):
    # convert a charMtrxFile to json format
    # <other_arguments> includes all the arguments accepted 
    # by the old laml to parase a character matrix file
         # TODO @Gillian       
        # step 1: read in the charMtrxFile
        ...
        # step 2: convert to json format
        ... 
        # step 3: write the json_instance to the outputFile
        ...
        # return False if any of the steps 1,2,3 failed
        return True
    
    def parse_json(self,dataFile):
    # dataFile is a json file that contains either a character matrix a an allele table 
    # dataFile must have a field named "dataType", whose value is either "charMtrx" or "alleleTable"
        dataType = ... # TODO: parse dataFile for this field @Gillian
        if dataType == "charMtrx":
        # data_struct will be a mapping: cell_name -> (cassette -> cassette_state)
             data_struct = {}
             # TODO: @Gillian
        else: # dataType = "alleleTable"
            # data_struct will be a mapping: cell_name -> (cassette -> (cassette_state -> count))
             data_struct = {}
            # TODO: @Gillian
        return data_struct, dataType

    def parse_prior(self,priorFile):
    # parse a prior file and return Q and alphabet
    # similar to the read_Q function in sequence_lib.py will be helpful, but 
    # returns alphabet in addition to Q
        # TODO: @Gillian 
        alphabet = Alphabet(...) # TODO
        return Q,alphabet 

    def get_from_path(self,dataFile,priorFile=None):  
    # dataFile is in json format
        if priorFile is not None:
            Q,alphabet = self.parse_prior(priorFile)
        else:
            Q,alphabet = ... # TODO: fill in the uniform prior - @Gillian

        data_struct, dataType = self.parse_json(dataFile)
        if dataType == "charMtrx":
            DLT_data = CharMtrx(data_struct,alphabet)
        else:     
            DLT_data = AlleleTable(data_struct,alphabet)
        return DLT_data,Q  
