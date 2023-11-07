# sc-MAIL: single-cell Maximum-likelihood Ancestral Inference for Lineage-tracing

sc-MAIL is a maximum likelihood algorithm under the Probabilistic Mixed-type Missing (PMM) model. Given a lineage tracing experiment character matrix with heterogeneous per-site alphabets and mutation probabilities, sc-MAIL will find a maximum likelihood tree topology and estimate branch lengths as well as stochastic dropout and heritable silencing missing data rates. 

# Dependencies
1. Please set up a MOSEK license.
2. `pip install cvxpy`


# Installation
For users:

Please clone the repository with:

```
git clone https://github.com/raphael-group/sc-mail.git
```
Please run the unit tests with:

```
python problin_tests.py
```

Although there are many more options available, sc-MAIL only strictly requires three arguments, using the following command.
```
python run_problin.py -t <topology> -c <characters> -o <output> 
```
<!-- # Reference -->

