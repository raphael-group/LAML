#! /usr/bin/env python
from statistics import mean
import pickle
import json
import pandas as pd
import ast

recognized_missing = set(['-', '?', '-1'])

def write_sequences(char_mtrx,nsites,outFile,delimiter=","):
    with open(outFile,'w') as fout:
        # header
        fout.write("cell")
        for i in range(nsites):
            fout.write(delimiter+"r" + str(i))
            #fout.write(delimiter+"r" + str(i+1))
        fout.write("\n")
        # actual data
        for cell in char_mtrx:
            fout.write(cell)
            for x in char_mtrx[cell]:
                fout.write(delimiter+str(x))
            fout.write("\n")


def old_read_sequences(inFile,filetype="charMtrx",cassette_len=1,set_single_cassette_as_tuple=False,delimiter=",",masked_symbol=None, suppress_warnings=False, replace_mchar='?'):
    with open(inFile,'r') as fin:
        if filetype == "fasta":
            if not suppress_warnings: 
                print("Warning: Reading " + str(inFile) + " as fasta file. Processing missing data in these files is not yet implemented.")
            return read_fasta(fin,cassette_len=cassette_len,set_single_cassette_as_tuple=set_single_cassette_as_tuple)
        elif filetype == "charMtrx":
            return read_charMtrx(fin,cassette_len=cassette_len,set_single_cassette_as_tuple=set_single_cassette_as_tuple,delimiter=delimiter,masked_symbol=masked_symbol,suppress_warnings=suppress_warnings,replace_mchar=replace_mchar)

def read_sequences(inFile,filetype="charMtrx",cassette_len=1,set_single_cassette_as_tuple=False,delimiter=",",masked_symbol=None, suppress_warnings=False, replace_mchar='?'):
    with open(inFile,'r') as fin:
        if filetype == "fasta":
            if not suppress_warnings: 
                print("Warning: Reading " + str(inFile) + " as fasta file. Processing missing data in these files is not yet implemented.")
            return read_fasta(fin,cassette_len=cassette_len,set_single_cassette_as_tuple=set_single_cassette_as_tuple)
        elif filetype == "charMtrx":
            return read_charMtrx(fin,cassette_len=cassette_len,set_single_cassette_as_tuple=set_single_cassette_as_tuple,delimiter=delimiter,masked_symbol=masked_symbol,suppress_warnings=suppress_warnings,replace_mchar=replace_mchar)
        elif filetype == "alleleTab":
            return read_alleleTab(fin, missing_char=masked_symbol) 
        elif filetype == "obsFeatures":
            return read_obsFeatures(fin, missing_char=masked_symbol) 
            pass

def read_obsFeatures(fin, missing_char):
    print("Reading observed features...")
    # data_struct: a mapping: cell_name -> (cassette -> cassette_state)
    data = json.load(fin)
    # ASSUMES the target site keys start with 'site_'

    #print(type(data))
    charMtrx_data_struct = {}
    alphabet_dict = dict()
    emission_dict = dict()

    K = len([x['cassette_id'] for x in data['cell_data'][0]['cassettes']])
    J = len([key for key in data['cell_data'][0]['cassettes'][0]['count_data'][0].keys() if key.startswith('site_')])
    print("K:", K, "J:", J)

    for cell_data in data['cell_data']:
        cell_name = cell_data['cell_name']
        # save into data_struct
        charMtrx_data_struct[cell_name] = [0] * K #np.zeros(K, dtype=list)
        emission_dict[cell_name] = dict()

        for cassette_data in cell_data['cassettes']:
            cassette_id = int(cassette_data['cassette_id']) - 1
            emission_dict[cell_name][cassette_id] = dict()

            # assert(len(cassette_data['count_data']) <= 1)
            if len(cassette_data['count_data']) == 1:
                cassette_state = cassette_data['count_data'][0]['allele']
                print("cassette_id", cassette_id)
                charMtrx_data_struct[cell_name][cassette_id] = cassette_state

                # each cassette, each site has an alphabet
                if cassette_id not in alphabet_dict:
                    alphabet_dict[cassette_id] = dict()
                for key in cassette_data['count_data'][0].keys():
                    if key.startswith('site_'):
                        state = cassette_data['count_data'][0][key]
                        if key not in alphabet_dict[cassette_id]:
                            alphabet_dict[cassette_id][key] = {0, missing_char}
                        alphabet_dict[cassette_id][key].add(state)

                        
                        # each site in allele
                        site_idx = int(key[5:])
                        cassette_emission_dict = cassette_data['count_data'][0]['emissions_' + key]
                        emission_dict[cell_name][cassette_id][site_idx] = dict()
                        obs_state = state
                        for k, v in cassette_emission_dict.items():
                            true_state = int(key.split('_')[-1]) # probGen_KDE_true_state_0
                            emission_dict[cell_name][cassette_id][site_idx][obs_state] = dict()
                            emission_dict[cell_name][cassette_id][site_idx][obs_state][true_state] = v
                    
            elif len(cassette_data['count_data']) == 0:
                charMtrx_data_struct[cell_name][cassette_id] = []

    alphabet_data_struct = []
    for cassette_index in alphabet_dict:
        new_cassette_list = []
        for site_index in range(J):
            new_target_site_list = list(alphabet_dict[cassette_index][f'site_{site_index}'])
            new_cassette_list.append(new_target_site_list)
        alphabet_data_struct.append(new_cassette_list)

    if emission_dict:
        print("Detected emission priors.")

    return charMtrx_data_struct, alphabet_data_struct, emission_dict


def read_alleleTab(fin, missing_char):
    ########## TODO ##########
    data = json.load(fin)
    # ASSUMES the target site keys start with 'site_'

    data_struct = {}
    alphabet_dict = dict()
    #print(data[0])

    K = len([x['cassette_id'] for x in data['cell_data'][0]['cassettes']])
    J = len([key for key in data['cell_data'][0]['cassettes'][0]['count_data'][0].keys() if key.startswith('site_')])
    print("K:", K, "J:", J)
    
    if J == 0:
        raise("FATAL error: read_alleleTab could not find any target sites because no keys starting with 'site_' were detected.")
        return None

    for cell_data in data['cell_data']:
        cell_name = cell_data['cell_name']
        # save into data_struct
        data_struct[cell_name] = []

        for cassette in cell_data['cassettes']:
            cassette_dict = dict()
            cassette_id = cassette['cassette_id']
            if cassette_id not in alphabet_dict:
                alphabet_dict[cassette_id] = dict()

            for cassette_entry in cassette['count_data']:
                cassette_state = tuple(cassette_entry['allele'])
                umi_count = cassette_entry['umi_count']

                for key in cassette_entry.keys():
                    if key.startswith('site_'):
                        state = cassette_entry[key]
                        if key not in alphabet_dict[cassette_id]:
                            alphabet_dict[cassette_id][key] = {0, missing_char}
                        alphabet_dict[cassette_id][key].add(state)

                # print(cassette_state)
                # if missing data, report empty dictionary
                if cassette_state == tuple([missing_char for _ in range(J)]):
                    pass
                else:
                    cassette_dict[cassette_state] = umi_count

                # how do i handle missing data?
            data_struct[cell_name].append(cassette_dict)

    alphabet_data_struct = []
    for cassette_index in alphabet_dict:
        new_cassette_list = []
        for site_index in range(J):
            new_target_site_list = list(alphabet_dict[cassette_index][f'site_{site_index}'])
            new_cassette_list.append(new_target_site_list)
        alphabet_data_struct.append(new_cassette_list)

    return data_struct, alphabet_data_struct

def read_fasta(fin,cassette_len=1,set_single_cassette_as_tuple=False):    
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

def check_missing(seen_missing, x):
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

def read_charMtrx(fin,cassette_len=1,set_single_cassette_as_tuple=False,delimiter=",",masked_symbol=None,suppress_warnings=False,replace_mchar='?',convert_to_int=True,stop_key=None):
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
            if check_missing(seen_missing, x):
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
    #if len(seen_missing) > 1 and not suppress_warnings:
    #    print("Warning: Found " + str(seen_missing) + " characters and treated them as missing.")
    #elif masked_symbol == None and len(seen_missing) >= 1 and not suppress_warnings:
    #    print("Warning: Reading sequences, detected " + str(seen_missing) + " as the missing character(s). We recommend explicitly providing the missing character.")
    return D, site_names    

def read_Q(inFile,has_head=True):
    Q = {}
    with open(inFile,'r') as fin:
        #fin.readline() # skip the header
        for line in fin:
            char,state,prob = line.strip().split(',')
            if not int(char) in Q:
                Q[int(char)] = {} #{int(state):float(prob)}
            Q[int(char)][int(state)] = float(prob)
    return Q

#from treeswift import *

def extract_brlens(tfile, ofile):

    t = read_newick_tree(tfile)
    t.root.h = 1

    with open(ofile, "w+") as w:
        for nidx, node in enumerate(t.traverse_preorder()):
            if not node.is_root():
                node.h = node.parent.h + 1
                s = nidx + " " + str(node.h) + " " + str(node.edge_length)
                w.write(s)

def alphabet_size(mtx):
    alphabet_sizes = []
    for char_idx in mtx:
        keys = mtx[char_idx]
        num_unique_keys = len(set(keys))
        alphabet_sizes.append(num_unique_keys)
    return max(alphabet_sizes), min(alphabet_sizes), mean(alphabet_sizes)

    #df = pd.DataFrame.from_dict(mtx, orient='index')
    #unique_series = df.nunique()
    #return max(unique_series), min(unique_series), mean(unique_series)

# adapted from /n/fs/ragr-research/projects/scmail_experiments/Real_biodata/test_kptracer/proc_scripts
def load_pickle(f):
# returns dictionary for use with Cassiopeia
    infile = open(f, "rb")
    priors = pickle.load(infile)
    infile.close()
    Q = dict()
    for i in sorted(priors.keys()):
        # scale to sum to 1
        q = {int(x): float(priors[i][x])/sum([float(c) for c in priors[i]]) for x in priors[i]}
        #q = {int(x):float(priors[i][x])/sum([float(c) for c in priors[i]]) for x in priors[i]}
        for x in q.keys():
            # print(q[x], x, q[x] < 1.0 and q[x] >= 0.0)
            assert q[x] < 1.0 and q[x] >= 0.0
            q[0] = 0.0

            Q[i] = q
            return Q

def read_emission(efile):
    tmp_df = pd.read_csv(efile)
    emission_dict = dict()
    for _, row in tmp_df.iterrows():
        cell_name = row.iloc[0]
        emission_dict[cell_name] = [] # list of dictionaries per site

        # Extract all columns after the first one
        other_columns = row.iloc[1:].to_dict()

        for column, value in other_columns.items():
            parsed_value = ast.literal_eval(value)
            emission_dict[cell_name].append(parsed_value)
    return emission_dict


def read_priors(pfile, msa, site_names=None):
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
                    else:
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
                        Q.append(Qi)
                        Qi = {}
                    seen_sites.add(site_idx)
                char_state = int(char_state)
                prob = float(prob)
                #Q[len(seen_sites) - 1][char_state] = prob
                Qi[char_state] = prob
            Q.append(Qi)
    elif file_extension == "json":
        with open(pfile, 'r') as fin:
            Q = json.load(fin)
    #else:
    #    Q = [{0:0} for i in range(k)]
    #    with open(pfile, 'r') as fin:
    #        for line in fin:
    #            site_idx, char_state, prob = line.strip().split()
    #            site_idx, char_state, prob = int(site_idx), int(char_state), float(prob)
    #            Q[site_idx][char_state] = prob '''
    return Q
                
