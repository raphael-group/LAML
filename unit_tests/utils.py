import scipy.stats as stats
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
        nzeros += sum([1 if (ch == 0 or ch == '0') else 0 for ch in seq])
        nmissing += sum([1 if ch == '?' else 0 for ch in seq])
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
        diffs.append(((float(oi) - ei)**2)/ei)
    return sum(diffs)

def calc_expected(node_label, d, k, c, Q, allreps, s):
    num_deg_freedom = len(allreps.keys()) - 1
    est_char = []
    for site_i in range(k):
        #if c == 0:
        #    print("Case Not Handled: provided character is 0, and is not in Q.")
        #if c == '?':
        #    print("Case Not Handled: provided character is '?', and is not in Q.")

        site_i_seq = [allreps[rep][node_label][site_i] for rep in allreps]
        # print(site_i_seq)
        est_char.append(sum([1 if ch == c else 0 for ch in site_i_seq])/len(site_i_seq))
    #print("est_char", est_char)
    if c != 0 and c != '?':
        qc = Q[site_i][c]
        exp_char = [qc * exp(-d * s) * (1 - exp(-d))] * len(site_i_seq)
    elif c == 0:
        exp_char = [ exp(-d * (1 + s) ) ] * len(site_i_seq)
        #exp_char = [1 - sum([Q[site_i][c] * exp(-d) * 1 - exp(-d)]) for c in Q[site_i].keys() - (1 - exp(-d)) ] * len(site_i_seq)
    else: # c == '?'
        exp_char = [1 - exp(-d * s)] * len(site_i_seq)
    #print("exp_char", exp_char)
    cst = chi_squared_tests(est_char, exp_char)
    #print('cst', cst)
    
    alpha = 0.05
    #print(stats.chi2.cdf(cst, num_deg_freedom))
    p_value = 1 - stats.chi2.cdf(cst, num_deg_freedom)
    #print('pval', p_value)
    if p_value <= alpha:
        # the variable does not have the expected distribution
        return False
    else:
        return True



