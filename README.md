# :camel: LAML: Lineage Analysis via Maximum Likelihood

LAML is a maximum likelihood algorithm under the Probabilistic Mixed-type Missing (PMM) model. Given a lineage tracing experiment character matrix with heterogeneous per-site alphabets and mutation probabilities, LAML will find a maximum likelihood tree topology and estimate branch lengths as well as stochastic dropout and heritable silencing missing data rates. 

For additional information about the method refer to the [paper](https://www.biorxiv.org/content/10.1101/2024.03.05.583638v1) and the [website](https://raphael-group.github.io/LAML/).

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

## Examples
We provide two examples for the two common use cases of LAML. See [Examples.md](https://github.com/raphael-group/laml/tree/master/Examples.md) for more details.

## Advanced I/O options
LAML has the following additional options for I/O
```
input options:
    --delimiter DELIMITER    The delimiter of the input character matrix. Can be one of {'comma','tab','whitespace'} .Default: 'comma'.
    -m MISSING_DATA, --missing_data MISSING_DATA Missing data character. Default: if not specified, assumes '?'.
    -p PRIORS, --priors PRIORS    The input prior matrix Q. Default: if not specified, use a uniform prior.
output options:
    -o OUTPUT, --output OUTPUT    Output prefix. Default: LAML_output
    -v, --verbose         Show verbose messages.
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

### The output prefix and verbose option
The user can change the output prefix using `-o`. The default prefix is `LAML_output`. The software can be run in verbose mode using `-v`.

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
