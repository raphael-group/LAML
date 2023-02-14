from problin_libs.sequence_lib import read_sequences, read_Q
import argparse
from treeswift import *

def pars_score_startle():
    pass

def pars_score(T, msa, mchar, norm, priorfile):
    def get_seq(n, nodedict):
        #print(n.label)
        if n.is_leaf():
            return msa[n.label]
        else:
            return nodedict[n]

    use_weighted = False
    if priorfile != "":
        use_weighted = True
        q = read_Q(priorfile)

    nodedict = dict()
    score = 0
    num_internal = 0
    for n in T.traverse_postorder():
        if n.is_leaf():
            nodedict[n] = get_seq(n, nodedict)

        else:
            num_internal += 1
            a, b = n.child_nodes()
            s1 = get_seq(a, nodedict)
            s2 = get_seq(b, nodedict)

            if n.is_root():
                # root_state has to be 0

            s = []
            # handle missing data
            # TODO: Handle root
            for cidx, xy in enumerate(zip(s1, s2)):
                x, y = xy

                if n.is_root():
                    # root state has to be 0
                    if x != '0':
                        if use_weighted:
                            score += -np.log(q[cidx][x])
                            score += -np.log(q[cidx][y])
                        else:
                            score += 2
                else:
                    if x == y: 
                        s.append(x)
                    elif x == mchar or y == mchar:
                        if x == mchar:
                            s.append(y) # the one that's not mchar
                        else:
                            s.append(x)
                    else: # x != y
                        if use_weighted:
                            score += -np.log(q[cidx][x])
                            score += -np.log(q[cidx][y])
                        else:
                            score += 2
                        s.append(0)

            nodedict[n] = s
    if norm:
        return score / num_internal
    else:
        return score

def main(args):
    tree = read_tree_newick(args.tree1)
    m, _ = read_sequences(args.msa,delimiter="\t",masked_symbol=args.mchar)
    #print(m)
    
    print(pars_score(tree, m, args.mchar, args.norm, args.prior))



if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--delimiter", type=str,
                        help="Input character matrix delimiter.",
                        default="\t",
                        required=True)
    parser.add_argument("-t1", "--tree1", type=str,
                        help="Input file containing tree.",
                        required=True)
    parser.add_argument("-msa", type=str,
                        help="Input file containing character matrix.",
                        required=True)
    parser.add_argument("-mchar", type=str,
                        help="Missing character.",
                        required=True)
    parser.add_argument("--norm", action="store_true",
                        help="Whether to compute the normalized parsimony score.")
    parser.add_argument("--prior", type="str",
                        required=False,
                        default="",
                        help="Prior file to help compute weighted parsimony score.")
    main(parser.parse_args())
