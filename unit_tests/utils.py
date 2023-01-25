import treeswift
import unittest
from problin_libs.sim_lib import *
from math import exp

def setup(tree, m, mu, k):
    Q = sim_Q(k, m)
    cmtx = simulate_seqs(tree, Q, mu)
    return cmtx

def count_missing(mtx):
    nzeros, nmissing, total = 0, 0, 0 
    for c in mtx:
        seq = mtx[c]
        nzeros += sum([1 if ch == 0 else 0 for ch in seq])
        nmissing += sum([1 if ch == -1 else 0 for ch in seq])
        total += len(seq)
    return nzeros, nmissing, total

def count_all(repdict):
    all_zero, all_missing, all_total = 0, 0, 0
    for i in repdict:
        mtx = repdict[i]
        nz, nm, total = count_missing(mtx)
        all_zero += nz
        all_missing += nm
        all_total += total
    return all_zero, all_missing, all_total

def chi_squared_tests(est_zeros, true_zeros):
    diffs = []
    for oi, ei in zip(est_zeros, true_zeros):
        diffs.append((float(oi - ei)**2)/ei)
    return sum(diffs)

def calc_expected(node_label, d, k, c, Q, allreps):
    csts = []
    for site_i in range(k):
        if c == 0:
            print("Case Not Handled: provided character is 0, and is not in Q.")
        elif c == -1:
            print("Case Not Handled: provided character is -1, and is not in Q.")
        qc = Q[site_i][c]
        site_i_seq = [allreps[rep][node_label][site_i] for rep in allreps]
        # print(site_i_seq)
        est_zeros = [1 if ch == c else 0 for ch in site_i_seq]
        true_zeros = [qc * exp(-d) * (1 - exp(-d))] * len(site_i_seq)
        cst = chi_squared_tests(est_zeros, true_zeros)
        csts.append(cst)
    return csts



