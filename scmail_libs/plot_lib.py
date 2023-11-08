from statistics import mean
import matplotlib.pyplot as plt
import argparse
import pandas as pd

lookup = {'nm': "model neither", 'bm': "model both", 'nd': "model heritable", 'ns': "model dropout", 'nomissing': "no missing", 'd0.1': 'dropout 0.1', 's0.01': 'silencing 0.01', 's0.32': 'silencing 0.32'}
colors = ['C0', 'C1', 'C2', 'C3']

def read_txtfile(t, cols):
    df = pd.read_csv(t)
    return df[cols]

def proc_to_sec(v):
    #for k, v in df.iteritems():
    m = v.split('m')[0]
    s = v.split('m')[1][:-1]
    m = int(m) * 60
    s = float(s)
    return m + s

def plot_time(t1, t2, o):
    fig, ax = plt.subplots()
    bp = ax.boxplot([t1, t2], notch=True, patch_artist=True, showmeans=True, showfliers=False)
    ax.set_ylim(ymin=0)
    ax.set_ylabel("Time (s)")
    ax.set_title("Compare EM vs no EM")
    ax.legend([bp["boxes"][0] for bp in bps], keys, loc='upper right')
    plt.savefig(o)

def plot_dtype(d, p, t, o):
    outfile = o + "_" + t + ".png"


    keys = []
    # fig = plt.figure(figsize=(10, 7))
    fig, ax = plt.subplots()
    bps = []

    for i, f in enumerate(d[t]):
        if f == 'nm':
           continue 
        # number of times this flag was the smallest nll
        count = sum([1 if min([p[t][f][rep]['nll'] for f in p[t].keys()]) is p[t][f][rep]['nll'] else 0 for rep in p[t][f].keys()])
        nll = round(mean([p[t][f][rep]['nll'] for rep in p[t][f].keys()]), 3)
        ed = round(mean([p[t][f][rep]['ed'] for rep in p[t][f].keys()]), 3)
        es = round(mean([p[t][f][rep]['es'] for rep in p[t][f].keys()]), 3)

        bp = ax.boxplot(d[t][f], positions=[i], notch=True, patch_artist=True, boxprops=dict(facecolor=colors[i]), showmeans=True, showfliers=False)
        #bp = ax.boxplot([mean(d[t][f][rep]) for rep in d[t][f]], positions=[i], notch=True, patch_artist=True, boxprops=dict(facecolor=colors[i]), showmeans=True)

        txt = 'count:' + str(count) + '\nblerr:' + str(round([item.get_ydata() for item in bp['means']][0][0], 3)) + '\nNLL:' + str(nll) + "\nED: " + str(ed) + "\nES: " + str(es)
        plt.text(i, 0.0, txt, bbox=dict(facecolor='pink', alpha=0.5), horizontalalignment='left', fontsize=8)
        bps.append(bp)
        keys.append(lookup[f])
    
    ax.set_ylim(ymin=0)#, ymax=0.5)
    ax.set_ylabel("Branch Length Error")
    ax.set_title("True Data (k=30): " + str(lookup[t]))
    ax.legend([bp["boxes"][0] for bp in bps], keys, loc='upper right')
    plt.savefig(outfile)

def main(f1, f2, out_pre, mu):
    datadict = dict()
    with open(f1, "r") as r:
        lines = r.readlines()[1:]
        for line in lines:
            rep, dtype, flags, tbl, ebl = line.split(',')
            if dtype not in datadict:
                datadict[dtype] = dict()
            err = abs(float(tbl) * float(mu) - float(ebl))
            if flags not in datadict[dtype]:
                datadict[dtype][flags] = []
                # datadict[dtype][flags] = dict()
            #if rep not in datadict[dtype][flags]:
            #    datadict[dtype][flags][rep] = []
            datadict[dtype][flags].append(err)
    paramdict = dict()
    with open(f2, "r") as r:
        lines = r.readlines()[1:]
        for line in lines:
            rep, dtype, flags, nll, ed, es = line.split(',')
            if dtype not in paramdict:
                paramdict[dtype] = dict()
            if flags not in paramdict[dtype]:
                paramdict[dtype][flags] = dict()
            if rep not in paramdict[dtype][flags]:
                paramdict[dtype][flags][rep] = dict()
            paramdict[dtype][flags][rep]['nll'] = float(nll)
            paramdict[dtype][flags][rep]['ed'] = float(ed)
            paramdict[dtype][flags][rep]['es'] = float(es)

    for dtype in datadict:
        plot_dtype(datadict, paramdict, dtype, out_pre)


