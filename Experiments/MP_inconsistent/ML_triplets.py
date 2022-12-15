from problin_libs.sequence_lib import read_sequences
from problin_libs.distance_based_lib import ML_pairwise_estimate,triplet_estimate

with open("ML_triplets_results.txt",'w') as fout:
    for k in [20,30,40,50,100,200,300,400,500,1000,5000]:
        S = read_sequences("seqs_m10_k" + str(k) + ".txt")
        n_total = 0
        n_correct = 0
        for D in S:
            a = D['a']
            b = D['b']
            c = D['c']
            d = D['d']

            dr_ab = ML_pairwise_estimate(a,b)[0][2]
            dr_ac = ML_pairwise_estimate(a,c)[0][2]
            dr_bc = ML_pairwise_estimate(b,c)[0][2]
            t_abc = triplet_estimate(dr_ab,dr_ac,dr_bc)
         
            dr_ad = ML_pairwise_estimate(a,d)[0][2]
            dr_cd = ML_pairwise_estimate(c,d)[0][2]
            t_acd = triplet_estimate(dr_ac,dr_ad,dr_cd)

            n_total += 1
            n_correct += (t_abc == 3 and t_acd == 1)
            print(k,n_total,n_correct)
        fout.write(str(k) + " " + str(n_correct/n_total) + "\n")
