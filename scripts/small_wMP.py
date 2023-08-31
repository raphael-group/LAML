from treeswift import *
from problin_libs.sequence_lib import * 
from small_MP import small_MP

if __name__ == "__main__":
    from sys import argv
    msa_file = argv[1]
    tree_file = argv[2]
    prior_file = argv[3]
    out_file = argv[4]

    msa, site_names = read_sequences(msa_file,filetype="charMtrx",delimiter=",",masked_symbol='?')
    with open(tree_file,'r') as fin:
        treeStr = fin.read().strip()
    
    Q = read_priors(prior_file,site_names)
    annTree,lb2seqs,wMP_cost = small_MP(treeStr,msa,Q=Q)
    with open(out_file,'w') as fout:
        fout.write(annTree + "\n")
        for lb in lb2seqs:
            fout.write(lb + "," + ",".join([str(x) for x in lb2seqs[lb]]) + "\n")
    print("wMP cost:" + str(wMP_cost))        
