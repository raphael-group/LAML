import pickle
from problin_libs.sequence_lib import read_sequences, write_sequences
import sys

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
        # q[0] = 0.0
        
        Q[i] = q
    return Q

def make_unique(fname, outfile, delimiter='\t', missing_char="?"):

    msa, site_names = read_sequences(fname, filetype="charMtrx", delimiter=delimiter)
    
    final_msa = dict()
    seen = set()

    for cellBC in msa:
        s = [str(x) for x in msa[cellBC]]
        s = ''.join(s)
        k = len(msa[cellBC])
        if s not in seen:
            # final_msa[cellBC] = msa[cellBC]
            final_msa[cellBC] = [x if x != missing_char else -1 for x in msa[cellBC] ]
            seen.add(s)
    
    write_sequences(final_msa, k, outfile, delimiter=',')



