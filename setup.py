from setuptools import setup, find_packages
import scmail_libs
from os import walk, listdir
from os.path import join, normpath, isfile

param = {
        'name': scmail_libs.PROGRAM_NAME,
        'version': scmail_libs.PROGRAM_VERSION,
        'description': scmail_libs.PROGRAM_DESCRIPTION,
        'author': scmail_libs.PROGRAM_AUTHOR,
        'license': scmail_libs.PROGRAM_LICENSE,
        'packages': find_packages(),
        'include_package_data': True,
        'scripts': ['run_scmail.py', 'scmail_tests.py'],
        'zip_safe': True,
        'install_requires': ['numpy>=1.16', 'treeswift>=1.1.37', 'scipy>=1.3.1', 'cvxpy>=1.4', 'mosek>=10.1.16'],
        'keywords': 'Phylogenetics Evolution Computational Maximum-likelihood Single-cell Lineage Tracing',
        'long_description': """sc-MAIL is a maximum likelihood algorithm under the Probabilistic Mixed-type Missing (PMM) model. Given a lineage tracing experiment character matrix with heterogeneous per-site alphabets and mutation probabilities, sc-MAIL will find a maximum likelihood tree topology and estimate branch lengths as well as stochastic dropout and heritable silencing missing data rates.""",
        'classifiers': [
            "Environment :: Console",
            "Intended Audience :: Developers",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: GNU General Public License (GPL)",
            "Natural Language :: English",
            "Operating System :: OS Independent",
            "Programming Language :: Python",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            ],
        }

setup(**param)
