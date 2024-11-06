from laml_libs.IO_handler.sequence_lib import *
from laml_libs.Count_model.AlleleTable import AlleleTable
from laml_libs.Count_model.Alphabet import Alphabet

class DLT_parser: # DLT: dynamic lineage tracing
    def __init__(self, datafile):
        self.datafile = datafile
        self.missing_state = "?"
        self.K = None
        self.J = None
        self.alphabet = None

    def set_alphabet(self):
        ### assumes self.data_struct has been set
        alphabet_ds = get_alphabet_ds(self.data_struct, self.datatype, self.missing_state)
        self.alphabet = Alphabet(self.K, self.J, alphabet_ds, silence_mechanism='convolve')
        return alphabet_ds

    def charMtrx_to_json(self,delimiter,missing_state,outputfile): 
    # convert a charMtrxFile to json format
    # <other_arguments> includes all the arguments accepted 
    # by the old laml to parse a character matrix file
        datafile = self.datafile
        try:
            # note, still only implemented for cassette length J = 1
            # step 1: read in the charMtrxFile
            charMtrx, site_names = read_sequences(datafile, filetype="charMtrx", delimiter=delimiter)
            # step 2: convert to json format
            # assume that J = 1, K = len(site_names)
            json_data = charMtrx_to_json(charMtrx, delimiter, missing_state)
            # step 3: write the json_instance to the outputFile
            write_json(json_data,outputfile)
            # return False if any of the steps 1,2,3 failed
            self.set_json_obj(json_data)
            self.datafile = outputfile
            return True
        except Exception as e:
            raise(f"FATAL error: {str(e)}")
            return False

    def set_json_obj(self, json_data):
        self.data = json_data["cell_data"]
        self.datatype = json_data["dataType"]
        self.K = len([x["cassette_idx"] for x in self.data[0]["cassettes"]])
        if self.datatype == "charMtrx":
            self.J = len(self.data[0]["cassettes"][0]["cassette_state"])
        else:
            self.J = len(self.data[0]["cassettes"][0]["cassette_state"][0]["state"])

    def parse_json(self):
    # dataFile is a json file that contains either a character matrix a an allele table 
    # dataFile must have a field named "dataType", whose value is either "charMtrx" or "alleleTable"
        data_obj = read_json(self.datafile)
        self.set_json_obj(data_obj)

        if self.datatype == "charMtrx":
            # data_struct will be a mapping: cell_name -> (cassette -> cassette_state)
             data_struct = data_to_CharMtrx_struct(self.data)
        else: # dataType = "alleleTable"
            # data_struct will be a mapping: cell_name -> (cassette -> (cassette_state -> count))
             data_struct = data_to_AlleleTable_struct(self.data)
        self.data_struct = data_struct
        return data_struct

    def parse_prior(self,priorFile,sites_per_cassette):
    # parse a prior file and return Q and alphabet
    # sites_per_cassette will be J in the constructed alphabet
        # Note: read_priors currently does NOT fill in the unedited cassette state!
        Q = read_priors(priorFile) 
        # Q will be processed into a list of (list of dictionaries)
        alphabet_ds = get_alphabet_prior(Q,sites_per_cassette)
        if self.K:
            if len(alphabet_ds) != self.K:
                raise(f"Prior file K={len(alphabet_ds)}, which does not match K={self.K} from previously provided data.")
        else:
            self.K = len(alphabet_ds)
        alphabet = Alphabet(self.K, sites_per_cassette, alphabet_ds) 
        self.alphabet = alphabet
        return alphabet_ds,alphabet 

    def get_from_path(self,dataFile,priorFile=None):  
    # dataFile is in json format
        if priorFile is not None:
            Q,alphabet = self.parse_prior(priorFile)
        else:
            self.set_alphabet()
            Q = uniform_priors(self.alphabet, self.datatype) # TODO: fill in the uniform prior - @Gillian

        data_struct, dataType = self.parse_json(dataFile)
        if dataType == "charMtrx":
            DLT_data = CharMtrx(data_struct,alphabet)
        else:     
            DLT_data = AlleleTable(data_struct,alphabet)
        return DLT_data,Q  
