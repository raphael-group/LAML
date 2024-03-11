# :camel: LAML: Lineage Analysis via Maximum Likelihood

LAML is a maximum likelihood algorithm under the Probabilistic Mixed-type Missing (PMM) model. Given a lineage tracing experiment character matrix with heterogeneous per-site alphabets and mutation probabilities, LAML will find a maximum likelihood tree topology and estimate branch lengths as well as stochastic dropout and heritable silencing missing data rates. 

For additional information about the method, you can refer to the [website](https://raphael-group.github.io/laml/).

# Installation
## Precursors 
The following precursors **are required** to install and run LAML
### Python
The software requires python >= 3.8.

### [IMPORTANT] MOSEK License
The software uses [MOSEK](https://www.mosek.com) for numerical optimization, which requires a license. Please do the following 2 steps:
1. Visit [this page](https://www.mosek.com/products/academic-licenses/) to get a **free academic license**. 
2. After obtaining the license file ``mosek.lic``, follow  [this page](https://docs.mosek.com/latest/licensing/client-setup.html) to place the license file in the correct place. 

## Install from PyPI
```
pip install laml 
```
After installation, type the following
```
run_laml -h
```
to see the commandline help of LAML.

## [Optional] Testing

Unit tests are available to ensure the success of installation. We highly recommend that the user performs the following step to test the installation.

In your terminal, type the following

```
laml_tests.py 
```
If LAML was installed properly, you would see on the screen `Running tests for LAML...` to begin, and print progress dots (one for each test passed). 
At the end, you should see the following message:
```
----------------------------------------------------------------------
Ran 80 tests in 13.486s

OK
```
## [For developers] Install from source
If you wish to install from source, do the following steps:
1. Clone the LAML github to your machine:
``git clone https://github.com/raphael-group/LAML.git``
2. Change directory to the ``LAML`` folder. Then use ``pip`` to install from source.
``pip install .``

# Usage
```
run_laml -c <character_matrix> -t <tree_topology> 
```
## Input
LAML requires the following two input files:
1. A file containing the character matrix, a [comma-separated values (CSV) file](https://en.wikipedia.org/wiki/Comma-separated_values) that has rows representing  cells and columns representing target sites. This file must have a header showing a list of site names and every subsequent line must begin with the cell name. Values of the character matrix must be either non-negative integers or '?', with 0 indicating the unmutated state, other integers indicating mutated state, and '?' as the missing data character. Refer to the paper for more details. 
2. A tree topology, given in [newick format](https://en.wikipedia.org/wiki/Newick_format). 

See an example character matrix in [examples/example1/character_matrix.csv](https://github.com/raphael-group/LAML/tree/laml/examples/example1/character_matrix.csv) and an example tree topology in [examples/example1/starting.tree](https://github.com/raphael-group/LAML/tree/laml/examples/example1/starting.tree)


## Output
There are four output files: 

1. `LAML_output_trees.nwk`: The output tree with time-resolved branch lengths
2. `LAML_output_params.txt`: This file reports the dropout rate, silencing rate, and negative log-likelihood.
3. `LAML_output_annotations.txt`: This file contains the inferred maximum likelihood sequences for all internal nodes and leaf nodes, with possible characters and associated probabilities for sites with more than one possibility.
4. `LAML_output.log`: The LAML logfile.

## Advanced I/O options
LAML has the following additional options for I/O
```
input options:
    --delimiter DELIMITER    The delimiter of the input character matrix. Can be one of {'comma','tab','whitespace'} .Default: 'comma'.
    -m MISSING_DATA, --missing_data MISSING_DATA Missing data character. Default: if not specified, assumes '?'.
    -p PRIORS, --priors PRIORS    The input prior matrix Q. Default: if not specified, use a uniform prior.
output options:
    -o OUTPUT, --output OUTPUT    Output prefix. Default: LAML_output
```

### Customize the format of the input character matrix 
The software allows some flexibility on the format of the input character matrix, using `-m` and `--delimiter` options. For example, if the character matrix is in the [tab-separated values (CSV) file](https://en.wikipedia.org/wiki/Tab-separated_values) format, it will still be accepted if `--delimiter Tab` is specified. The placeholder of the missing character can also be adjusted using `-m`. For instance, if the input file has missing entries represented by "-" instead of "?", it will still be accepted if `-m -` is specified. 

Note: LAML also accepts a character matrix that contains negative integers and/or non-alphanumeric values and treats them all as a placeholder for missing entries. However, for best practices, the user should explicitly specify their missing data character using `-m`.

### The mutation priors
While not strictly required, **mutation priors** can have a large effect on the outputs. If no mutation priors are provided, LAML uses uniform priors
by default. However, if possible we highly recommend specifying mutation prior using `-p`. We accept the following two formats for mutation priors 

**Recommended** A file containing the prior matrix, a [comma-separated values (CSV) file](https://en.wikipedia.org/wiki/Comma-separated_values), with three columns: site index, character state, and probability. The site index and character states must be integers, and the probability must be a float. We do *not* expect the unmutated state to appear in the alphabet. See an example input prior file in
[examples/example1/priors.csv](https://github.com/raphael-group/LAML/tree/laml/examples/example1/priors.csv).

**Not recommended** We also accept [Python-pickled files](https://docs.python.org/3/library/pickle.html#data-stream-format), as this is the indel prior output format for [Cassiopeia](https://cassiopeia-lineage.readthedocs.io/en/latest/notebooks/reconstruct.html). We print a warning if the keys of the pickled prior dictionary do not match the site names in your provided character matrix file. 

### The output prefix
The user can change the output prefix using `-o`. The default prefix is `LAML_output`.

## Examples
To try the following examples, first do the followings:
1. Download the data from [examples.zip](https://github.com/raphael-group/laml/tree/master/examples.zip)
2. Unzip the downloaded file. After unzipping, you should see a folder named ``examples``
3. Change directory to ``examples``
```
  cd examples
```
### Use Case 1: Infer time-resolved branch lengths and heritable missing and dropout rates on a fixed topology
LAML can infer time-resolved branch lengths and the rates of the two missing data types for a fixed tree topology. If the time frame of the experiment is specified by ``--timescale``, the output tree will be scaled to the same height. Otherwise, the output tree will be scaled to the unit height 1.

For example, the following command:
```
run_laml -c examples/example1/character_matrix.csv -t examples/example1/starting.tree -p examples/example1/priors.csv -o example1 --nInitials 1 --timescale 10
```
specifies the tree via ``-t`` and set ``--timescale`` to 10. Running this command will produce three output files
1. `example1_trees.nwk`: the output tree containing time-resolved branch lengths. This tree has the same topology as the starting tree specified in `-t`, but has branch lengths in time units
2. `example1_params.txt`: this file reports the dropout rate, silencing rate, the negative log-likelihood of the tree topology and parameters, and the mutation rate
3. `example1_annotations.txt`: This file has two components
   (i) the newick string of the rooted tree with internal nodes labeled and branch lengths show the infer *number of mutations*.
   (ii) imputed sequences for each node in the tree. If a site has multiple possible states, it is annotated with the probability of each possible state.


We provide sample outputs in `examples/out_example1/` for your reference. 

### Use Case 2: Infer tree topology, branch lengths, and missing data rates
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

## Other options
Below are some other important options available in LAML. For full documentation, please run `run_laml -h`.

### Numerical optimization
```
  -L COMPUTE_LLH, --compute_llh COMPUTE_LLH Compute log-likelihood of the input tree using the input (phi,nu). Will NOT optimize branch lengths, phi, or nu. The input tree MUST have branch lengths. This option has higher priority than --topology_search and --resolve_search.
  --noSilence         Assume there is no gene silencing, but allow missing data by dropout in sc-sequencing. Does not necessarily produce ultrametric trees, and cannot be time-scaled. This option has higher priority than --timescale or --ultrametric.
  --noDropout           Assume there is no sc-sequencing dropout, but allow missing data by gene silencing.
  --timescale TIMESCALE Timeframe of experiment. Scales ultrametric output tree branches to this timescale. Default: 1.0.
  --solver SOLVER       Specify a solver. Options are 'Scipy' or 'EM'. Default: EM
  --nInitials NINITIALS    The number of initial points. Default: 20.
```
### Topology search
```
  --topology_search     Perform topology search using NNI operations. Always returns a fully resolved (i.e. binary) tree.
  --resolve_search      Resolve polytomies by performing topology search ONLY on branches with polytomies. This option has higher priority than --topology_search.
  --keep_polytomies     Keep polytomies while performing topology search. This option only works with --topology_search.
  --parallel            Turn on parallel version of topology search.
  --randomreps RANDOMREPS    Number of replicates to run for the random strategy of topology search.
  --maxIters MAXITERS   Maximum number of iterations to run topology search.
```
### Other options
```
  -v, --verbose         Show verbose messages.
```
