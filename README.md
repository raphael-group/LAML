# sc-MAIL: single-cell Maximum-likelihood Ancestral Inference for Lineage-tracing

sc-MAIL is a maximum likelihood algorithm under the Probabilistic Mixed-type Missing (PMM) model. Given a lineage tracing experiment character matrix with heterogeneous per-site alphabets and mutation probabilities, sc-MAIL will find a maximum likelihood tree topology and estimate branch lengths as well as stochastic dropout and heritable silencing missing data rates. 

# Precursors

1. First, we ask users to set up a MOSEK license. Please refer to the MOSEK installation page [here](https://www.mosek.com/products/academic-licenses/).
2. Add the license file to your path by adding an environment variable. For instance, add the following to your `.bashrc` and load [This page](https://docs.mosek.com/latest/licensing/client-setup.html) may be useful to reference.

```
export MOSEKLM_LICENSE_FILE=<path_to_folder_containing_mosek_license>
```

# Installation

## Installing from source

For users:

1. Please clone the repository with:

```
git clone https://github.com/raphael-group/sc-mail.git
```

2. Run the setup script.
```
python setup.py install 
```
You can run it with `--prefix=<your_preferred_install_dir>` but be sure to set this prefix to your preferred PYTHONPATH.

3. (optional) Please run the unit tests with:

```
python scmail_tests.py
```
You can comment out lines 1 and 2 if you'd like the unit tests to run faster. The full test suite runs in about ~9 minutes on my Linux machine.

## Installing from pip/conda

in progress...


# Running sc-MAIL

Although there are many more options available, sc-MAIL only strictly requires three arguments, using the following command.
```
python run_scmail.py -t <topology> -c <characters> -o <output> 
```

The output consists of three files: 

1. `<prefix>_annotations.txt`: this file contains the inferred maximum likelihood sequences for all internal nodes and leaf nodes, with possible characters and associated probabilities for sites with more than one possibility.
2. `<prefix>_params.txt`: this file contains the dropout rate, silencing rate, and negative log likelihood.
3. `<prefix>_trees.nwk`: this file contains the rooted estimated tree with branch lengths.



## Examples

Note that if you get an error in the following runs (especially use case 3), please make sure you have installed the MOSEK license file.

### Use Case 1: Infer branch lengths on a topology

From the `sc-mail/` directory, please run the following code:
```
python run_scmail.py -c examples/character_matrix.csv -t examples/starting.tree -p examples/priors.csv --delimiter comma -o example1 --nInitials 1
```

This will output three files. You can compare these outputs with those in `examples/out_example1/`.


### Use Case 2: Compute the likelihood of an existing tree

From the `sc-mail/` directory, please run the following code:
```
python run_scmail.py -c examples/character_matrix.csv -t examples/starting.tree -p examples/priors.csv --delimiter comma -o example2 -L "0 4.879273344239771e-07" --solver Scipy
```

This will output three files. You can compare these outputs with those in `examples/out_example2/`.

### Use Case 3: Infer a topology

From the `sc-mail/` directory, please run the following code:
```
python run_scmail.py -c examples/character_matrix.csv -t examples/starting.tree -p examples/priors.csv --delimiter comma -o example3 --nInitials 1 --randomreps 1 --topology_search -v --ultrametric --parallel --randseeds 1984
```

This will output four files. You can compare **the likelihood** of the resulting tree with the results in `examples/out_example3/`. Note that when performing topology search, a checkpoint file will be generated (and updated) as well. Note that this will resolve all polytomies, run in parallel, and return an ultrametric tree.



## Documentation


```
-t TOPOLOGY, --topology TOPOLOGY    Binary input tree topology in newick format. Branch lengths will be ignored.
-c CHARACTERS, --characters CHARACTERS  The input character matrix. Must have header.
-p PRIORS, --priors PRIORS  The input prior matrix Q. Default: if not specified, use a uniform prior.
--delimiter DELIMITER   The delimiter of the input character matrix. Can be one of {'comma','tab','whitespace'} .Default: 'tab'.
-m MASKEDCHAR, --maskedchar MASKEDCHAR  Masked character. Default: if not specified, assumes '-'.
-o OUTPUT, --output OUTPUT  Output prefix.
--solver SOLVER       Specify a solver. Options are 'Scipy' or 'EM'. Default: EM
--topology_search     Perform topology search using NNI operations. Always return fully resolved (i.e. binary) tree.
--resolve_search      Resolve polytomies by performing topology search ONLY on branches with polytomies. This option has higher
priority than --topology_search.
--keep_polytomies     Keep polytomies while performing topology search. This option only works with --topology_search.
-L COMPUTE_LLH, --compute_llh COMPUTE_LLH   Compute likelihood of the input tree using the input (phi,nu). Will NOT optimize branch lengths, phi, or nu. The input tree MUST have branch lengths. This option has higher priority than --topoloy_search and
--resolve_search.
--ultrametric         Enforce ultrametricity to the output tree.
--noSilence           Assume there is no gene silencing, but allow missing data by dropout in sc-sequencing.
--noDropout           Assume there is no sc-sequencing dropout, but allow missing data by gene silencing.
-v, --verbose         Show verbose messagesRandom seeds. Can be a single integer number or a list of integers whose length is equal to the number of random seeds. Can be a single integer number or a list of integers whose length is equal to the number of.
--nInitials NINITIALS   The number of initial points. Default: 20.
--randseeds RANDSEEDS   Random seeds. Can be a single integer number or a list of integers whose length is equal to the number of initial points (see --nInitials).
--randomreps RANDOMREPS
Number of replicates to run for the random strategy of topology search.
--maxIters MAXITERS   Maximum number of iterations to run topology search.
--parallel            Turn on parallel version of topology search.
```

