#! /usr/bin/env python

from problin_libs.EM_solver import EM_solver
from problin_libs.ML_solver import ML_solver

Q = [{1:1}]
#T = "((a:1,b:1):1,c:1):1;"
#T = "((a:0,b:0):0,c:0):0;"

#T = "((b:9.99970114511368e-07,c:9.99970114511368e-07):9.99970114511368e-07,(d:9.999994414313543,a:9.999994414313543):2.030364475095228):9.99970114511368e-07;"
#T = "((d:9.999911853281668,a:9.999911853281668):1.1767319042305595,(b:9.99665230034214e-07,c:9.99665230034214e-07):9.99665230034214e-07):9.99665230034214e-07;"
#T = "((b:1e-06,c:1e-06):1e-06,(a:7.73987165244907,d:5.770329390739723):9.506374609533808):1e-06;"
T = "((b:1e-06,c:1e-06):0,(a:5.66798137443661e-06,d:5.736079723850543e-06):10):1e-06;"
#T = "(b:1e-06,c:1e-06,(a:5.66798137443661e-06,d:5.736079723850543e-06):10):1e-06;"

#T = "(b:0,c:0,(a:0,d:100):100):0;"

#msa = {'a':[1],'b':[1],'c':[1]}
Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
msa = {'a':[1, 1, 1, 1, 1], 'b':[0, 0, 0, 0, 0], 'c':[0, 0, 0, 0, 0], 'd':[1, 1, 1, 1, 1]}
#mySolver = ML_solver(T,{'charMtrx':msa},{'Q':Q},{'phi':0,'nu':0})
mySolver = ML_solver(T,{'charMtrx':msa},{'Q':Q},{'phi':0,'nu':0})
print(mySolver.score_tree(strategy={'ultra_constr':True,'optimize':True}))
#print(mySolver.score_tree(strategy={'ultra_constr':False,'optimize':False}))
print(mySolver.tree.newick())
#print(mySolver.score_tree(strategy={'ultra_constr':False,'optimize':False}))
