#! /usr/bin/env python
from statistics import mean
import pickle

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


def read_sequences(inFile,filetype="charMtrx",delimiter=",",masked_symbol=None, suppress_warnings=False, replace_mchar='?'):
    with open(inFile,'r') as fin:
        if filetype == "fasta":
            if not suppress_warnings: 
                print("Warning: Reading " + str(inFile) + " as fasta file. Processing missing data in these files is not yet implemented.")
            return read_fasta(fin)
        elif filetype == "charMtrx":
            return read_charMtrx(fin,delimiter=delimiter,masked_symbol=masked_symbol,suppress_warnings=suppress_warnings,replace_mchar=replace_mchar)

def read_fasta(fin):    
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

def read_charMtrx(fin,delimiter=",",masked_symbol=None,suppress_warnings=False,replace_mchar='?',convert_to_int=True,stop_key=None):
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
                    seq.append(replace_mchar)
                else:
                    seq.append(x)    
            else:
                if convert_to_int:
                    seq.append(int(x))
                else:    
                    seq.append(x)
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
    #else:
    #    Q = [{0:0} for i in range(k)]
    #    with open(pfile, 'r') as fin:
    #        for line in fin:
    #            site_idx, char_state, prob = line.strip().split()
    #            site_idx, char_state, prob = int(site_idx), int(char_state), float(prob)
    #            Q[site_idx][char_state] = prob '''
    return Q
                
