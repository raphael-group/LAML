#! /usr/bin/env python
from statistics import mean
import pickle

recognized_missing = set(['-', '?', '-1'])

def write_sequences(char_mtrx,nsites,outFile,delimiter=","):
    with open(outFile,'w') as fout:
        # header
        fout.write("cell")
        for i in range(nsites):
            fout.write(delimiter+"r" + str(i+1))
        fout.write("\n")
        # actual data
        for cell in char_mtrx:
            fout.write(cell)
            for x in char_mtrx[cell]:
                fout.write(delimiter+str(x))
            fout.write("\n")


def read_sequences(inFile,filetype="charMtrx",delimiter=",",masked_symbol=None, suppress_warnings=False):
    with open(inFile,'r') as fin:
        if filetype == "fasta":
            if not suppress_warnings: 
                print("Warning: Reading " + str(inFile) + " as fasta file. Processing missing data in these files is not yet implemented.")
            return read_fasta(fin)
        elif filetype == "charMtrx":

            return read_charMtrx(fin,delimiter=delimiter,masked_symbol=masked_symbol)

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
    elif int(x) < 0: # check positivity
        return True 
    else:
        return False

def read_charMtrx(fin,delimiter=",",masked_symbol=None):
    D = {}
    site_names = fin.readline().strip().split(delimiter)[1:]
    

    if masked_symbol != '-':
        seen_missing = set([masked_symbol])
    else: 
        seen_missing = set([])

    for line in fin:
        line_split = line.strip().split(delimiter)
        name = line_split[0]
        # check if any are missing characters or nonnegative
        seq = []
        for x in line_split[1:]:
            if check_missing(seen_missing, x):
                seen_missing.add(x)
                seq.append('?')
            else:
                seq.append(int(x))
        #seq = [int(x) if x != masked_symbol else "?" for x in line_split[1:]]
        D[name] = seq
    if len(seen_missing) > 1:
        print("Warning: Found " + str(seen_missing) + " characters and treated them as missing.")
    elif masked_symbol == None:
        print("Warning: Reading sequences, detected " + str(seen_missing) + " as the missing character(s). We recommend explicitly providing the missing character.")
    return D, site_names    

def read_Q(inFile):
    Q = {}
    with open(inFile,'r') as fin:
        fin.readline() # skip the header
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

# adapted from /n/fs/ragr-research/projects/problin_experiments/Real_biodata/test_kptracer/proc_scripts
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
