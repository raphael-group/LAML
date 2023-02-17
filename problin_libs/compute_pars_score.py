from problin_libs.sequence_lib import read_sequences
import argparse
from treeswift import *



def pars_score(T, msa, mchar):
    def get_seq(n, nodedict):
        #print(n.label)
        if n.is_leaf():
            return msa[n.label]
        else:
            return nodedict[n]

    nodedict = dict()
    score = 0
    for n in T.traverse_postorder():
        if n.is_leaf():
            nodedict[n] = get_seq(n, nodedict)
        else:
            a, b = n.child_nodes()
            s1 = get_seq(a, nodedict)
            s2 = get_seq(b, nodedict)

            s = []
            # handle missing data
            for x, y in zip(s1, s2):
                if x == y: 
                    s.append(x)
                elif x == mchar or y == mchar:
                    if x == mchar:
                        s.append(y) # the one that's not mchar
                    else:
                        s.append(x)
                else:
                    score += 1
                    s.append(0)

            nodedict[n] = s
    return score

def main(args):
    tree = read_tree_newick(args.tree1)
    m, _ = read_sequences(args.msa,delimiter=",",masked_symbol="-1")
    #print(m)
    
    print(pars_score(tree, m, args.mchar))



if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-t1", "--tree1", type=str,
                        help="Input file containing tree",
                        required=True)
    parser.add_argument("-msa", type=str,
                        help="Input file containing character matrix",
                        required=True)
    parser.add_argument("-mchar", type=str,
                        help="Missing character",
                        required=True)
    main(parser.parse_args())
