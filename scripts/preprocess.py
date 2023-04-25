import json
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
        #for x in q.keys():
            # print(q[x], x, q[x] < 1.0 and q[x] >= 0.0)
        #    assert q[x] < 1.0 and q[x] >= 0.0
        #q[0] = 0.0
        
        Q[i] = q
    return Q

def make_unique(fname, outfile, eqfile, delimiter='\t', missing_char="?", droplenti=False):

    msa, site_names = read_sequences(fname, filetype="charMtrx", delimiter=delimiter)
    
    final_msa = dict()
    seen = dict()
    mappings = dict()

    for cellBC in msa:
        s = [str(x) for x in msa[cellBC]]
        s = ''.join(s)
        k = len(msa[cellBC])
        if s not in seen:
            # final_msa[cellBC] = msa[cellBC]
            final_msa[cellBC] = [x if x != missing_char else -1 for x in msa[cellBC] ]
            seen[s] = cellBC
            mappings[cellBC] = []
        else:
            seed_cellBC = seen[s]
            mappings[seed_cellBC].append(cellBC)

    if droplenti:
        for cellname in final_msa:
            final_msa[cellname] = final_msa[cellname][1:]
        site_names = site_names[1:]
        k -= 1
    
    write_sequences(final_msa, k, outfile, delimiter=',')
    
    with open(f"{eqfile}", "w") as f:
        json.dump(mappings, f)

    return mappings

def add_identical(tree, eqfile):
    # Destructive: Modifies tree, adding polytomies!

    l2n = tree.label_to_node(selection="all")
    with open(eqfile, "r") as f:
        mappings = json.load(f)

    for seed_label in mappings:
        # get corresponding node
        seed_node =l2n[seed_label]
        for same_label in mappings[seed_label]:
            child = Node(label=same_label)
            seed_node.add_child(child)


def norm_Q(p, outfile, cmtx):
    p = load_pickle(p)
    cmtxfile, site_names = read_sequences(cmtx, filetype="charMtrx", delimiter=",")
    Q = dict()
    for col_i, site_name in zip(p, site_names):
        Q_i = p[col_i]
        keys = list(Q_i.keys()).copy()
        for x in keys:
            if Q_i[x] <= 0 or Q_i[x] >= 1:
                Q_i.pop(x)
        s = sum([Q_i[x] for x in Q_i])
        Q_i_norm = {}
        for x in Q_i:
            if Q_i[x] > 0 and Q_i[x] < 1:
                Q_i_norm[x] = Q_i[x]/s
        Q[int(site_name[1:])] = Q_i_norm
    pickle.dump(Q, open(outfile, "wb"))

def write_pickle(p, poutfile):
    with open(poutfile, "wb") as fout:
        pickle.dump(p, fout)

