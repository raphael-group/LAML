import treeswift
import argparse


# python ${problindir}/problin_libs/compile_results.py --tt "" -bm --pre ${odir}/results_ML_${data_name} --out_param ${paramfile}p


def main(pre, out_param, out_bl, tt, nm, bm, nd, ns, rep, datatype):
    out_dict = dict()
    if nm:
        out_dict['nm'] = parse_resultfile(pre + "_nm.txt")
    if bm:
        out_dict['bm'] = parse_resultfile(pre + "_bm.txt")
    if nd:
        out_dict['nd'] = parse_resultfile(pre + "_nd.txt")
    if ns:
        out_dict['ns'] = parse_resultfile(pre + "_nd.txt")
    
    write_outparam(out_dict, out_param, rep, datatype)

    process_outfile(out_dict, out_bl, tt, rep, datatype)

def parse_resultfile(pre):
    with open(pre, "r") as r:
        lines = r.readlines()

        out_dict = dict()
        for line in lines:

            if line[:13] == 'Optimal tree:':
                est_t = line.split("Optimal tree:")[1]
                out_dict['tree'] = est_t
            elif line[:21] == 'Optimal negative-llh:':
                nllh = line.split("Optimal negative-llh:")[1]
                out_dict['llh'] = nllh 
            elif line[:21] == 'Optimal dropout rate:':
                param_d = line.split("Optimal dropout rate:")[1]
                out_dict['dropout'] = param_d 
            elif line[:23] == 'Optimal silencing rate:':
                param_s = line.split("Optimal silencing rate:")[1]
                out_dict['silencing'] = param_s 

        return out_dict

def write_outparam(out_dict, outfile, rep, datatype):
    with open(outfile, "w+") as w:
        header = 'Rep,Datatype,Flags,LLH,Dropout,Silencing\n'
        w.write(header)
        for flag in out_dict:
            s = str(rep) + "," + str(datatype) + "," + ','.join([flag, out_dict[flag]['llh'].rstrip('\n'), out_dict[flag]['dropout'].rstrip('\n'), out_dict[flag]['silencing'].rstrip('\n')]) + '\n'
            w.write(s)

def process_outfile(out_dict, outfile, tt, rep, datatype):
    # Datatype, branchLenError
    truetree = treeswift.read_tree_newick(tt)
    with open(outfile, "a+") as w:
        w.write("Rep,Datatype,Flags,TrueBL,EstBL\n")
        for dtype in out_dict.keys():
            dtype_dict = out_dict[dtype]
            esttree = treeswift.read_tree_newick(dtype_dict['tree'])

            edge_dict = dict()
            for node, dst in truetree.root.traverse_bfs():
                edge_dict[node.label] = [node.get_edge_length()]
            for node, dst in esttree.root.traverse_bfs():
                edge_dict[node.label].append(node.get_edge_length())

            for node_label in edge_dict.keys():
                w.write(str(rep) + "," + str(datatype) + "," + dtype + "," + str(edge_dict[node_label][0]) + "," + str(edge_dict[node_label][1]) + "\n")


parser = argparse.ArgumentParser()
parser.add_argument('-p', '--pre')
parser.add_argument('-r', '--rep')
parser.add_argument('-op', '--out_param')
parser.add_argument('-obl', '--out_bl')
parser.add_argument('--tt')
parser.add_argument('--datatype')
parser.add_argument('-nm', action='store_true', default=False)
parser.add_argument('-bm', action='store_true', default=False)
parser.add_argument('-nd', action='store_true', default=False)
parser.add_argument('-ns', action='store_true', default=False)

args = parser.parse_args()

pre = args.pre
out_param = args.out_param
out_bl = args.out_bl
truetree = args.tt
rep = args.rep
datatype = args.datatype
nm, bm, nd, ns = args.nm, args.bm, args.nd, args.ns

main(pre, out_param, out_bl, truetree, nm, bm, nd, ns, rep, datatype)

