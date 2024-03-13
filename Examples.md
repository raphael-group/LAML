# Prerequisites
To try the following examples, first do the following:
1. Download the data from [examples.zip](https://github.com/raphael-group/laml/tree/master/examples.zip)
2. Unzip the downloaded file. After unzipping, you should see a folder named ``examples``
3. Change directory to ``examples``
```
  cd examples
```
# Use Case 1: Infer time-resolved branch lengths, heritable missing rate, and dropout rates of a fixed tree topology
LAML can infer time-resolved branch lengths and the rates of the two missing data types for a fixed tree topology. If the time frame of the experiment is specified by ``--timescale``, the output tree will be scaled to the same height. Otherwise, the output tree will be scaled to the unit height 1.

For example, the following command:
```
run_laml -c examples/example1/character_matrix.csv -t examples/example1/starting.tree -p examples/example1/priors.csv -o example1 --nInitials 1 --timescale 10 -v
```
specifies the tree via ``-t`` and set ``--timescale`` to 10. Running this command will produce three output files
1. `example1_trees.nwk`: the output tree containing time-resolved branch lengths. This tree has the same topology as the starting tree specified in `-t`, but has branch lengths in time units
2. `example1_params.txt`: this file reports the dropout rate, silencing rate, the negative log-likelihood of the tree topology and parameters, and the mutation rate
3. `example1_annotations.txt`: This file has two components
   (i) the newick string of the rooted tree with internal nodes labeled and branch lengths show the infer *number of mutations*.
   (ii) imputed sequences for each node in the tree. If a site has multiple possible states, it is annotated with the probability of each possible state.


We provide sample outputs in `examples/out_example1/` for your reference. 

# Use Case 2: Infer tree topology, branch lengths, and missing data rates
LAML can simultaneously infer tree topology, branch lengths, and the missing data rates using the ``--topology_search`` option.

For example, the following command:
```
run_laml -c examples/example2/character_matrix.csv -t examples/example2/starting.tree -p examples/example2/priors.csv -o example2 --nInitials 1 --randomreps 1 --topology_search -v --timescale 10
```
enables topology search using the flag ``--topology_search``. 

Running this command will produce three output files 
1. `example2_trees.nwk`: the output tree with time-resolved branch lengths. Because topology search has been performed, this tree has a different topology from the starting tree. The new topology has higher likelihood than the starting tree.
2. `example2_params.txt`: this file reports the dropout rate, silencing rate, the negative log-likelihood of the tree topology and parameters, and the mutation rate
3. `example2_annotations.txt`: This file has two components
   (i) the newick string of the rooted tree with internal nodes labeled and branch lengths show the infer *number of mutations*.
   (ii) imputed sequences for each node in the tree. For sites with multiple possible states, that site is annotated with the probability of each possible state.

In addition, a checkpoint file `example2._ckpt.<randomnumber>.txt` is produced, which is important for running LAML on large data. Every 50 NNI iterations, this file is updated with a checkpoint containing (i) the NNI iteration number, (ii) the current best newick tree, (iii) the current best negative log-likelihood, (iv) the current best dropout rate, and (v) the current best silencing rate.

We provide sample outputs in `examples/out_example2/` for your reference. 
