import sys, os
#import cassiopeia as cas
import pandas as pd
import numpy as np

#import cassiopeia as cas

from itertools import combinations

# fname = "/Users/gc3045/laml2_experiments/simulations/set_2/test_input/c100_noise0_rep0/c100_noise0_rep0.perturbed_n0.04_cmat.csv"
fname, outname = sys.argv[1], sys.argv[2]
cmat = (pd.read_csv(fname, index_col=0)
          .replace("?", -1))
cmat.index = cmat.index.astype(str) 
cmat = cmat.astype(np.int64)
#cmat.loc["root"] = 0
#print(cmat)
tree = cas.data.CassiopeiaTree(character_matrix=cmat)
#tree.root_sample_name = 'root'

#dmat = tree.get_dissimilarity_map()
#print(dmat)
nj_solver = cas.solver.NeighborJoiningSolver(
    cas.solver.dissimilarity.weighted_hamming_distance,
    add_root=True, #False, 
    fast=True
)

#print("tree priors:", tree.priors)

dmat = tree.compute_dissimilarity_map(cas.solver.dissimilarity.weighted_hamming_distance)
dmat = nj_solver.get_dissimilarity_map(tree)
nj_solver.solve(tree)

dmat = tree.get_dissimilarity_map()
print(dmat)
"""
missing_state_indicator = -1
def test_whd(s1, s2, weights=None):

    d = 0
    num_present = 0
    for i in range(len(s1)):
        if s1[i] == missing_state_indicator or s2[i] == missing_state_indicator:
            continue

        num_present += 1

        if s1[i] != s2[i]:
            if s1[i] == 0 or s2[i] == 0:
                if weights:
                    if s1[i] != 0:
                        d += weights[i][s1[i]]
                    else:
                        d += weights[i][s2[i]]
                else:
                    d += 1
            else:
                if weights:
                    d += weights[i][s1[i]] + weights[i][s2[i]]
                else:
                    d += 2

    if num_present == 0:
        return 0

    return d / num_present

print("run the cas.solver.dissimilarity.weighted_hamming_distance")
taxa = cmat.index.to_list()
n = len(taxa)
dmat = np.zeros((n,n), dtype=float)
dmat_test = np.zeros((n,n), dtype=float)
for i in range(n):
    si = cmat.iloc[i].to_numpy()
    for j in range(i+1, n):
        sj = cmat.iloc[j].to_numpy()

        d = cas.solver.dissimilarity.weighted_hamming_distance(si, sj)
        dmat[i, j] = dmat[j, i] = d

        test_d = test_whd(si, sj)

        print(i, si)
        print(j, sj)
        print(i, j, test_d)
        if test_d != d:
            print(d, test_d, i, j)

        dmat_test[i, j] = dmat_test[j, i] = test_d

dmat = pd.DataFrame(dmat, index=taxa, columns=taxa)
print(dmat)

dmat_test = pd.DataFrame(dmat_test, index=taxa, columns = taxa)
print(dmat_test)
"""

tstr = tree.get_newick()

outdir = os.path.dirname(outname)             
if outdir and not os.path.isdir(outdir):      
    os.makedirs(outdir, exist_ok=True)        

with open(outname, "w") as fh:
    fh.write(tstr)

print(f"Wrote NJ tree to {outname}.")
