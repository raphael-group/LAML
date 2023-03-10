# problin

Notebooks can be found in `notebooks/`. Script development (files containing `main()` functions) can be found in the `scripts/` folder. The library folder `problin_libs/` will be kept clean.

<!--The notebook `likelihood_felsensteins.ipynb` provides code to calculate the likelihood of a tree using Felsenstein's pruning algorithm. Some edits have been made to the probability calculation, as well as the states which are considered, which are specific to the lineage tracing scenario. We provide three trees as an example in this notebook, where `t_C` is the maximum parsimony tree, but `t_A` and `t_B` both have better likelihood (and would be found by a maximum likelihood search over the topologies).  -->


TODOs:
[] Change all the print statements to logging statements?
[] Take startle import out 
[] Take parsimony computation out of the `run_problin.py` script.
[] Move necessary functions into library out of scripts, and make all scripts stand-alone with main() functions.
