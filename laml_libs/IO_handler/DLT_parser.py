from laml_libs.Count_model.AlleleTable import AlleleTable
from laml_libs.Count_model.CharMtrx import CharMtrx
from laml_libs.Count_model.Alphabet import Alphabet
import json

recognized_missing = set(['-', '?', '-1'])

class DLT_parser: # DLT: dynamic lineage tracing

    def __init__(self, datafile=None, priorfile=None, unedited_state=0, missing_state="?", delimiter=",", outputfile=None, max_allele_per_cassette=None):
        # intended user facing
        self.datafile = datafile
        self.priorfile = priorfile 
        self.unedited_state = unedited_state
        self.missing_state = missing_state 
        self.delimiter = delimiter
        self.outputfile = outputfile
        self.max_allele_per_cassette = int(max_allele_per_cassette) if max_allele_per_cassette is not None else max_allele_per_cassette

        self.datatype = None
        self.DLT_data = None
        self.K = None
        self.J = None
        #self.alphabet = None
        self.num_cells = None
        self.priors = None
        self.alphabet = None

        if datafile is not None: 
            self.get_from_path(self.datafile, self.delimiter, self.missing_state, self.outputfile, self.priorfile, self.max_allele_per_cassette)
    
    def process_datafile(self, datafile, delimiter, missing_state, outputfile, max_allele_per_cassette):
        file_extension = datafile.strip().split(".")[-1]
        if file_extension == "csv" or file_extension == "txt":
            # overwrites the datafile to be the json file
            self.datafile_to_json(datafile, delimiter, missing_state, outputfile)
        # reads in the json datafile
        self.parse_json(self.datafile, max_allele_per_cassette)
        

    def set_alphabet(self):
        if self.data_struct is None:
            self.parse_json(self.datafile, self.max_allele_per_cassette)
            if self.data_struct is None:
                raise(f"self.data_struct has not yet been set because the provided {self.datafile} cannot be parsed.")

        ### assumes self.data_struct has been set
        #print("[set_alphabet] self.data_struct", self.data_struct)
        alphabet_ds = self.get_alphabet_ds(self.data_struct, self.datatype, self.missing_state)
        alphabet = Alphabet(K=self.K, J=self.J, data_struct=alphabet_ds, silence_mechanism='convolve')
        self.alphabet = alphabet
        return alphabet_ds, alphabet

    def datafile_to_json(self, datafile, delimiter, missing_state, outputfile): #,datafile,delimiter,missing_state,outputfile): 
        # convert a charMtrxFile to json format
        # <other_arguments> includes all the arguments accepted 
        # by the old laml to parse a character matrix file
        self.datafile = datafile
        self.delimiter = delimiter
        self.missing_state = missing_state

        if outputfile == None:
            file_path = ''.join(datafile.strip().split(".")[:-1])
            outputfile = file_path + ".json"
        self.outputfile = outputfile

        try:
            # note, still only implemented for cassette length J = 1
            # step 1: read in the charMtrxFile
            charMtrx, site_names = self.read_sequences(self.datafile, filetype="charMtrx", delimiter=self.delimiter)
            # step 2: convert to json format
            # assume that J = 1, K = len(site_names)
            json_data = self.charMtrx_to_json(charMtrx, self.delimiter, self.missing_state)
            # step 3: write the json_instance to the outputFile
            self.write_json(json_data,self.outputfile)
            #self.set_json_obj(json_data)

            # save json output as this object's datafile
            self.datafile = self.outputfile
            return True
        except Exception as e:
            raise(f"FATAL error: {str(e)}")
            return False
            # return False if any of the steps 1,2,3 failed

    def set_json_obj(self, json_data):
        self.data = json_data["cell_data"]
        self.datatype = json_data["dataType"]
        self.K = max([len(self.data[i]["cassettes"]) for i in range(len(self.data))])
        #print("[set_json_obj]", "K", self.K)
        if self.datatype == "charMtrx":
            self.J = len(self.data[0]["cassettes"][0]["cassette_state"])
            #print("[set_json_obj]", "J", self.J)
        else:
            self.J = len(self.data[0]["cassettes"][0]["cassette_state"][0]["state"])
        self.num_cells = len(self.data)

    def parse_json(self, datafile, max_allele_per_cassette=None):
        # intended user facing
        # dataFile is a json file that contains either a character matrix a an allele table 
        # dataFile must have a field named "dataType", whose value is either "charMtrx" or "alleleTable"
        self.datafile = datafile
        # will overwrite previous value for max_allele_per_cassette
        self.max_allele_per_cassette = int(max_allele_per_cassette) if max_allele_per_cassette is not None else max_allele_per_cassette

        data_obj = self.read_json(self.datafile)
        self.set_json_obj(data_obj)

        if self.datatype == "charMtrx":
            # data_struct will be a mapping: cell_name -> (cassette -> cassette_state)
             data_struct = self.data_to_CharMtrx_struct(self.data)
        else: # dataType = "alleleTable"
            # data_struct will be a mapping: cell_name -> (cassette -> (cassette_state -> count))
             data_struct = self.data_to_AlleleTable_struct(self.data)
        self.data_struct = data_struct
        return data_struct

    def parse_prior(self,priorFile,sites_per_cassette):
        # parse a prior file and return Q and alphabet
        # sites_per_cassette will be J in the constructed alphabet
        self.priorfile = priorFile
        # Note: read_priors currently does NOT fill in the unedited cassette state!
        Q = self.read_priors(priorFile) 
        # Q will be processed into a list of (list of dictionaries)
        # Q is a cassette list of site lists with dictionaries for each site
        Q = self.get_alphabet_prior(Q,sites_per_cassette,self.datatype)

        if sites_per_cassette == 1:
            alphabet_ds = [[list(set([0,-1]+list(k[0].keys())))] for k in Q]
        else:
            alphabet_ds = [[list(set([0,-1]+list(s.keys()))) for s in k] for k in Q]
        if len(alphabet_ds) != self.K:
            #print("alphabet_ds", len(alphabet_ds), "K", self.K)
            raise(f"Prior file K={len(alphabet_ds)}, which does not match K={self.K} from previously provided data.")
        #else:

        self.K = len(alphabet_ds)
        alphabet = Alphabet(self.K, sites_per_cassette, alphabet_ds) 
        self.alphabet = alphabet
        return Q #, alphabet_ds #alphabet 

    def get_from_path(self, datafile, delimiter=None, missing_state=None, outputfile=None, priorfile=None, max_allele_per_cassette=None):  
        # intended user facing
        self.datafile = datafile 
        self.priorfile = priorfile 

        self.process_datafile(datafile, delimiter, missing_state, outputfile, max_allele_per_cassette)
        #data_struct = self.parse_json(self.datafile, max_allele_per_cassette)

        alphabet_ds, charMtrx_alphabet = self.set_alphabet()
        # dataFile is in json format
        if priorfile is not None:
            # resets self.alphabet according to the priorfile
            Q = self.parse_prior(priorfile, sites_per_cassette=self.J)
        else:
            #print("getting uniform priors")
            # uses the self.alphabet according to the datafile
            Q = self.uniform_priors(self.alphabet, self.datatype, self.K, self.J, self.unedited_state, self.missing_state) 
            #print("done with uniform priors")

        self.priors = Q

        self.max_allele_per_cassette = int(max_allele_per_cassette) if max_allele_per_cassette is not None else max_allele_per_cassette
        if self.datatype == "charMtrx":
            #print("[get_from_path]", self.alphabet.M)
            DLT_data = CharMtrx(self.data_struct,self.alphabet) #,self.max_allele_per_cassette) #self.alphabet)
        else:     
            DLT_data = AlleleTable(self.data_struct,self.alphabet,self.max_allele_per_cassette)

        self.DLT_data = DLT_data

        return self.DLT_data, self.priors

    def get_alphabet_ds(self, ds, dt, missing_state):
        # Note, does not include missing states
        alphabet_ds = {}
        final_alphabet_ds = []
        #print("[get_alphabet_ds] datatype", dt)
        
        if dt == "charMtrx": 
            # ds is a dictionary of lists
            for cell_name in ds:
                cell_data = ds[cell_name]
                for cassette_idx, cassette_state in enumerate(cell_data):
                    if len(final_alphabet_ds) <= cassette_idx:
                        final_alphabet_ds.append({0, -1})

                    if cassette_state != missing_state:
                        final_alphabet_ds[cassette_idx].add(cassette_state)

                #final_alphabet_ds.append([list(cassette_alphabet_ds)])
            final_alphabet_ds = [[list(x)] for x in final_alphabet_ds]

        else:
            # ds is a dictionary of lists of dictionaries
            for cell_name in ds:
                cell_data = ds[cell_name]

                for cassette_idx, cassette_data in enumerate(cell_data): #.keys():
                    cassette_state_dict = cell_data[cassette_idx]

                    if cassette_idx not in alphabet_ds.keys():
                        alphabet_ds[cassette_idx] = {} # set of possible cassette_states

                    for cassette_state in cassette_state_dict.keys():
                        for site_idx, site_state in enumerate(cassette_state):
                            if site_idx not in alphabet_ds[cassette_idx].keys():
                                alphabet_ds[cassette_idx][site_idx] = {0,-1}
                            if site_state != missing_state:
                                alphabet_ds[cassette_idx][site_idx].add(site_state)

            # order the indices of the cassette idx
            sorted_cassette_keys = [key for key in sorted(alphabet_ds.keys())]
            for cassette_idx in sorted_cassette_keys:
                # to produce the cassette_list, sort the site indices, and sort the alphabet for each site index
                cassette_list = [sorted(list(alphabet_ds[cassette_idx][key])) for key in sorted(alphabet_ds[cassette_idx].keys())]
                final_alphabet_ds.append(cassette_list)

        #print("final_alphabet_ds", final_alphabet_ds)
        return final_alphabet_ds

    def charMtrx_to_json(self,charMtrx,delimiter,missing_char):
        result_json = {"dataType": "charMtrx", "cell_data": []}
        for cell_name in charMtrx:
            cell_data = charMtrx[cell_name]
            cell_json = {"cell_name": cell_name, "cassettes": []}
            for cassette_idx, cassette_state in enumerate(cell_data):
                if cassette_state != missing_char:
                    cassette_json = {"cassette_idx": cassette_idx,
                                     "cassette_state": (cassette_state,)
                                     }
                else:
                    cassette_json = {"cassette_idx": cassette_idx,
                                     "cassette_state": ()
                                     }
                cell_json["cassettes"].append(cassette_json)
            result_json["cell_data"].append(cell_json)
        return result_json
   
    def write_json(self,json_data,outfile):
        json_string = json.dumps(json_data, indent=4)
        with open(outfile, "w") as f:
            f.write(json_string)

    def read_json(self,datafile):
        with open(datafile, "r") as f:
            data = json.load(f)
        return data


    def data_to_CharMtrx_struct(self,data):
        # Dictionary of lists

        charmat_ds = {}
        #alphabet_ds = {}

        for cell_json in data:
            cell_name = cell_json["cell_name"]
            cassettes = cell_json["cassettes"]
            charmat_ds[cell_name] = []

            for cassette_data in cassettes:
                cassette_idx = cassette_data["cassette_idx"]
                cassette_state = cassette_data["cassette_state"]

                if cassette_state != []:
                    charmat_ds[cell_name].append(cassette_state[0])
                else:
                    charmat_ds[cell_name].append("?")

        return charmat_ds #, alphabet_ds 

    def data_to_AlleleTable_struct(self, data):
        # Dictionary of lists of dictionaries
        # missing data should still have a cassette, but cassette_states can be 
        alleletable_ds = {}
        #alphabet_ds = {}

        for cell_json in data:
            cell_name = cell_json["cell_name"]
            cassettes = cell_json["cassettes"]
            alleletable_ds[cell_name] = []

            for cassette_data in cassettes:
                cassette_idx = cassette_data["cassette_idx"]
                cassette_states = cassette_data["cassette_state"]
                cell_cassette_data = {}

                for state_data in cassette_states: 

                    if state_data != {}: 
                        state = tuple(state_data["state"])
                        observation = state_data["observation"]
                        cell_cassette_data[state] = observation

                alleletable_ds[cell_name].append(cell_cassette_data)

        return alleletable_ds #, alphabet_ds 

    def get_alphabet_prior(self, Q, J, dt):
        def convert_keys_to_tuples(d):
            new_d = {}
            for key in d:
                new_d[(key,)] = d[key]
            return new_d
        alphabet_from_priordict = []
        for idx, q_dict in enumerate(Q):
            #if dt == "charMtrx":
                #new_q_dict = convert_keys_to_tuples(q_dict)
                #new_q_dict[(0,)] = 0
            #    new_q_dict = q_dict
            #    new_q_dict[0] = 0
            #    q_dict = new_q_dict
            #else:
            q_dict = q_dict[0]
            q_dict[0] = 0
            if idx % J == 0:
                cassette_dictionary = []
            cassette_dictionary.append(q_dict)
            if idx % J == J-1:
                alphabet_from_priordict.append(cassette_dictionary)
        return alphabet_from_priordict 

    def uniform_priors(self,alphabet, dt, K, J, unedited_state, missing_state):
        def generate_q(M_i, unedited_state):
            if len(M_i) == 0:
                # add pseudo-mutated state
                m_i = 1
                q = {"1":1.0}
            else:
                m_i = len(M_i)
                q = {x:1/m_i for x in M_i}
            q[unedited_state] = 0
            return q
        Q = []
        if dt == "charMtrx":
            # list of possible site alphabets
            for k in range(K):
                M_i = set([x for x in alphabet.get_cassette_alphabet(k) if x not in [unedited_state, missing_state, -1, "-1", "?"]])
                q = generate_q(M_i, unedited_state)
                Q.append(q)
        else:
            # a list of cassettes, where each cassette is a list of site alphabets
            for k in range(K):
                cassette_q = []
                for j in range(J):
                    M_i = set([x for x in alphabet.get_site_alphabet(k,j) if x not in [unedited_state, missing_state, -1, "-1", "?"]])
                    q = generate_q(M_i, unedited_state)
                    cassette_q.append(q)
                Q.append(cassette_q)

        return Q


    
    def read_sequences(self, inFile,filetype="charMtrx",cassette_len=1,set_single_cassette_as_tuple=False,delimiter=",",masked_symbol=None, suppress_warnings=False, replace_mchar='?'):
        with open(inFile,'r') as fin:
            if filetype == "fasta":
                if not suppress_warnings: 
                    print("Warning: Reading " + str(inFile) + " as fasta file. Processing missing data in these files is not yet implemented.")
                return self.read_fasta(fin,cassette_len=cassette_len,set_single_cassette_as_tuple=set_single_cassette_as_tuple)
            elif filetype == "charMtrx":
                return self.read_charMtrx(fin,cassette_len=cassette_len,set_single_cassette_as_tuple=set_single_cassette_as_tuple,delimiter=delimiter,masked_symbol=masked_symbol,suppress_warnings=suppress_warnings,replace_mchar=replace_mchar)

    def read_fasta(self, fin,cassette_len=1,set_single_cassette_as_tuple=False):    
        if cassette_len != 1:
            ########## TODO ##########
            raise("FATAL error: read_fasta is not yet implemented for cassette_len != 1")
            return None
        # The following code only works for cassette_len = 1
        S = [] # will be a list of dictionaries
        D = {}
        for line in fin:
            if line.startswith(">"):
                name = line.strip()[1:]
            elif line.startswith("_"):
                S.append(D)
                D = {}
            else:
                seq = [int(x) for x in line.strip().split("|")]
                if set_single_cassette_as_tuple:
                    seq = [tuple([x]) for x in seq]
                D[name] = seq       
        return S

    def check_missing(self, seen_missing, x):
        # returns whether character x is a missing character
        if x in seen_missing:
            return True
        elif x in recognized_missing:
            return True
        elif x.isalpha(): # check alphanumeric 
            return True 
        else: # check positivity
            try:
                return int(x) < 0
            except:
                return False

    def read_charMtrx(self, fin,cassette_len=1,set_single_cassette_as_tuple=False,delimiter=",",masked_symbol=None,suppress_warnings=False,replace_mchar='?',convert_to_int=True,stop_key=None):
        if cassette_len != 1:
            ########## TODO ##########
            raise("FATAL error: read_charMtrx is not yet implemented for cassette_len != 1")
            return None
        
        # The following code only works for cassette_len = 1
        D = {}
        site_names = fin.readline().strip().split(delimiter)
        if site_names[0] == 'cell_name' or "cell" in site_names[0]:
            site_names = site_names[1:]

        if masked_symbol != None:
            seen_missing = set([masked_symbol])
        else: 
            seen_missing = set([])

        for line in fin:
            if stop_key is not None and line.startswith(stop_key):
                break
            line_split = line.strip().split(delimiter)
            name = line_split[0]
            # check if any are missing characters or nonnegative
            seq = []
            for x in line_split[1:]:
                if self.check_missing(seen_missing, x):
                    seen_missing.add(x)                
                    if replace_mchar is not None:
                        #seq.append(replace_mchar)
                        c = replace_mchar
                    else:
                        #seq.append(x)   
                        c = x
                else:
                    if convert_to_int:
                        #seq.append(int(x))
                        c = int(x)
                    else:    
                        #seq.append(x)
                        c = x
                if set_single_cassette_as_tuple:
                    c = tuple([c])
                seq.append(c)            
            D[name] = seq
        return D, site_names    

    def read_priors(self, pfile, msa=None, site_names=None):
        file_extension = pfile.strip().split(".")[-1]
        if file_extension == "pkl" or file_extension == "pickle": #pickled file
            infile = open(pfile, "rb")
            priors = pickle.load(infile)
            infile.close()
            Q = []
            priorkeys = sorted(priors.keys())
            mapping = dict()
            
            if site_names is not None and priorkeys != sorted([int(x[1:]) for x in site_names]):
                print("Prior keys mismatch with site names.")
                print("Prior keys:", priorkeys)
                print("Site names:", site_names)
                print("Attempting to infer mapping between site names and prior keys...")

                # check if we have the same number of keys
                if len(site_names) == len(priorkeys):
                    for i, site_name in enumerate(site_names):
                        mapping[site_name] = priorkeys[i]
                    print(mapping)
                else:
                    # inferring prior and site_name mapping
                    # we should have a mapping from every site name to a prior dictionary

                    # compute offset 
                    site_name_digits = []
                    for site_name in site_names:
                        digit_name = ''.join([x for x in site_name if x.isdigit()])
                        site_name_digits.append(int(digit_name))
               
                    # tell the difference between an offset and missing keys in the dictionary
                    all_site_names_present = True 
                    for i, site_name in enumerate(site_names):
                        digit_name = site_name_digits[i]
                        if digit_name not in priors.keys():
                            all_site_names_present = False

                    if not all_site_names_present:
                        print("Not all site names are present in the dictionary. Trying offset...")
                        offset = min(site_name_digits) - min(priorkeys)
                        print("Offset between input site names and prior keys is assumed to be", offset)

                    for i, site_name in enumerate(site_names):
                        digit_name = site_name_digits[i]
                        if not all_site_names_present:
                            prior_name = digit_name - offset
                        else:
                            prior_name = digit_name
                        mapping[site_name] = prior_name

                        if prior_name in priors.keys():
                            q = {int(x):priors[prior_name][x] for x in priors[prior_name]}
                        elif msa is not None:
                            print(f"Missing priors at site {site_name}, filling in uniform priors...")
                            # fill in uniform priors at site i
                            M_i = set(msa[x][i] for x in msa if msa[x][i] not in [0,"?"])
                            if len(M_i) == 0:
                                # add pseudo-mutated state
                                m_i = 1
                                q={"1":1.0}
                            else:
                                m_i = len(M_i)
                                q = {x:1/m_i for x in M_i}
                            q[0] = 0
                        else:
                            raise(f"Missing priors at site {site_name}, pass in MSA as input to fill in uniform priors.")
                        Q.append(q)
                    print(mapping)
                    return Q


            for i in sorted(priors.keys()):
                q = {int(x):priors[i][x] for x in priors[i]}
                q[0] = 0
                Q.append(q)

        elif file_extension == "csv":
            #k = len(site_names)
            #Q = [{0:0} for i in range(k)]
            Q = []
            Qi = {}
            seen_sites = set()
            with open(pfile, 'r') as fin:
                lines = fin.readlines()
                tokens = lines[0].split(',')
                if not tokens[1].isnumeric() and not tokens[2].isnumeric():
                    lines = lines[1:]

                # check if the first character of the character name is a string
                token = lines[0].split(',')[0]
                charname_is_str = not token.isnumeric()

                for line in lines:
                    site_idx, char_state, prob = line.strip().split(',')
                    if charname_is_str:
                        site_idx = int(site_idx[1:])
                    else:
                        site_idx = int(site_idx)
                    if site_idx not in seen_sites:
                        if len(seen_sites) > 0:
                            Qi = [Qi]
                            Q.append(Qi)
                            Qi = {}
                        seen_sites.add(site_idx)
                    char_state = int(char_state)
                    prob = float(prob)
                    #Q[len(seen_sites) - 1][char_state] = prob
                    Qi[char_state] = prob
                Qi = [Qi]
                Q.append(Qi)
        #else:
        #    Q = [{0:0} for i in range(k)]
        #    with open(pfile, 'r') as fin:
        #        for line in fin:
        #            site_idx, char_state, prob = line.strip().split()
        #            site_idx, char_state, prob = int(site_idx), int(char_state), float(prob)
        #            Q[site_idx][char_state] = prob '''
        #print("Q:", Q)
        return Q
