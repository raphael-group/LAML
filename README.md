# sc-MAIL: single-cell Maximum-likelihood Ancestral Inference for Lineage-tracing

sc-MAIL is a maximum likelihood algorithm under the Probabilistic Mixed-type Missing (PMM) model. Given a lineage tracing experiment character matrix with heterogeneous per-site alphabets and mutation probabilities, sc-MAIL will find a maximum likelihood tree topology and estimate branch lengths as well as stochastic dropout and heritable silencing missing data rates. 

For additional information about the method, you can refer to the [website](https://raphael-group.github.io/sc-mail/).
# Precursors (required before installation)

We ask that users use python >= 3.6 with our code. 

## MOSEK License
1. First, we ask users to set up a MOSEK license. The preferred option is to place the license file mosek.lic in the directory mosek in the user’s home directory. Please refer to the MOSEK installation page [here](https://www.mosek.com/products/academic-licenses/).
2. (Optional): If you decide to add the license file elsewhere, please add the license file to your path by adding an environment variable. For instance, add the following line to your `.bashrc` and load (`source ~/.bashrc` for Linux/Unix and `. ~/.bashrc` for Windows). [This page](https://docs.mosek.com/latest/licensing/client-setup.html) may be useful to reference.

```
export MOSEKLM_LICENSE_FILE=<path_to_folder_containing_mosek_license>
```

# Installation
If you've fulfilled the required precursor steps, you can pick one of two ways to install sc-MAIL. We recommend installing using `pip`, but you can also install from source.

### Install using pip (recommended)

0. Please set up the MOSEK license. 

1. sc-MAIL is available on the Python Package Index (PyPI). To install, use `pip` as follows:
```
pip install scmail"
```
This will install the scmail package to your default package location. If you would like to specify a separate installation location, you can use the flag: `--prefix="<your_preferred_install_dir>`, then set this prefix in your PATH and PYTHONPATH (see below for help).

If `pip` installs the package in a directory which is not on path, `pip` will throw a warning and ask the user to consider adding this directory to PATH. This should be heeded (see below for help).

2. If you open a python interpreter as follows:

```
$ python
Python 3.9.16 (main, Sep 12 2023, 00:00:00)
[GCC 11.3.1 20221121 (Red Hat 11.3.1-4)] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import scmail_libs
>>>
```

You can now import functions from `scmail_libs`.

3. You can also run the following:

```
run_scmail
```
to see the commandline help of sc-MAIL.

## Install from source

1. Please clone the repository with:

```
git clone https://github.com/raphael-group/sc-mail.git
```

2. Run the setup script from inside the sc-mail directory. 
```
cd sc-mail
python setup.py install 
```
You can (for example, if you are running on a server and get permission denied when you try to install it in the default location), run it with `--prefix=<your_preferred_install_dir>` but be sure to set this prefix to your preferred PYTHONPATH (see below for help).

3. Set up the MOSEK license file. See above (section Precursors).

4. (optional) Please run the unit tests with:

```
$ run_scmail_tests
```
The lightweight test suite runs in about 2 minutes on my Linux machine.

After installation, run:

```
run_scmail -h
```
to see the commandline help of sc-MAIL.

## (optional) Setting prefix in environment variable 

To set the prefix `<your_preferred_install_dir>` in your PATH and PYTHONPATH, follow these steps:

(Unix/Linux users) In a Terminal:
```
vi ~/.bashrc
export PATH=$PATH:<your_preferred_install_dir>"
export PYTHONPATH=$PYTHONPATH:<your_preferred_install_dir>"
```

<!--export PATH=$PATH/pkgs/scmail-0.5-py311_0/bin:<your_preferred_install_dir>"-->
(Windows users) In a Command Prompt or PowerShell:
```
setx PATH "%PATH%;C:<your_preferred_install_dir>"
setx PYTHONPATH "%PYTHONPATH%;C:<your_preferred_install_dir>"
```

For both users, be sure to restart any applications or shells you want to use the updated PATH variable.


# Usage

If you downloaded using `pip` or `conda`, you should download the examples from github from [examples.zip](https://github.com/raphael-group/sc-mail/tree/master/examples.zip), unzip it, and run each command below from the directory containing this examples folder. If you downloaded from source, you can run `run_scmail` after the setup.

Although there are many more options available, sc-MAIL only strictly requires three arguments, using the following command:
```
$ run_scmail -t <topology> -c <characters> -o <output> 
```

The output consists of three files: 

1. `<prefix>_annotations.txt`: this file contains the inferred maximum likelihood sequences for all internal nodes and leaf nodes, with possible characters and associated probabilities for sites with more than one possibility.
2. `<prefix>_params.txt`: this file contains the dropout rate, silencing rate, and negative log likelihood.
3. `<prefix>_trees.nwk`: this file contains the rooted estimated tree with branch lengths.
4. (optional) `<prefix>._ckpt.<randomnumber>.txt`: this file contains checkpoints in the topology search process.


## Examples

Note that if you get an error in trying the following commands (especially use case 3), please make sure you have set up the MOSEK license file.

### Use Case 1: Infer branch lengths on a topology

From the `sc-mail/` directory, please run the following code:
```
$ run_scmail -c examples/character_matrix.csv -t examples/starting.tree -p examples/priors.csv --delimiter comma -o example1 --nInitials 1 --randseeds 1984
```

This will output three files (`example1_annotations.txt`, `example1_params.txt`, `example1_trees.nwk`). You can compare these outputs with those in `examples/out_example1/`. For instance, in order to compare the likelihoods, display the contents of the two files using the following (if on Linux/Unix):
```
$ cat example1_params.txt
$ cat examples/out_example1/example1_params.txt
```
or (if on Windows in Command Prompt):
```
$ type example1_params.txt examples/out_example1/example1_params.txt
```


### Use Case 2: Compute the likelihood of an existing tree

From the `sc-mail/` directory, please run the following code:
```
$ run_scmail -c examples/character_matrix.csv -t examples/starting.tree -p examples/priors.csv --delimiter comma -o example2 -L "0 4.879273344239771e-07" --solver Scipy
```

This will output three files (`example2_annotations.txt`, `example2_params.txt`, `example2_trees.nwk`). You can compare these outputs with those in `examples/out_example2/`. For instance, in order to compare the likelihoods, display the contents of the two files using the following (if on Linux/Unix):
```
$ cat example2_params.txt
$ cat examples/out_example2/example2_params.txt
```
or (if on Windows in Command Prompt):
```
$ type example2_params.txt examples/out_example2/example2_params.txt
```

### Use Case 3: Infer a topology

From the `sc-mail/` directory, please run the following code:
```
$ run_scmail -c examples/character_matrix.csv -t examples/starting.tree -p examples/priors.csv --delimiter comma -o example3 --nInitials 1 --randomreps 1 --topology_search -v --ultrametric --parallel
```

This will output four files (`example3_annotations.txt`, `example3_params.txt`, `example3_trees.nwk`, `example3._ckpt.<randnumber>.txt`). When performing topology search, a checkpoint file is also generateed. Note that this command will resolve all polytomies, run in parallel, and return an ultrametric tree.

You can compare these outputs with those in `examples/out_example3/`. For instance, in order to compare the likelihoods, display the contents of the two files using the following (if on Linux/Unix):
```
$ cat example3_params.txt
$ cat examples/out_example3/example3_params.txt
```
or (if on Windows in Command Prompt):
```
$ type example3_params.txt examples/out_example3/example3_params.txt
```


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

## Troubleshooting

If you installed using `pip` and would like to check the version of scmail you have installed:
```
pip show scmail
```

