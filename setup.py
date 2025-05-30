from setuptools import setup, find_packages
from setuptools.command.install import install
import sysconfig
import os
import laml_libs 
from os import walk, listdir
from os.path import join, normpath, isfile

def recursive_list_dir(path):
    listing=[]
    for x in walk(path):
        if isfile(x[0]):
            listing.append(x[0].split(path+'/')[1])
        for y in listdir(x[0]):
            z = normpath(join(x[0],y))
            if isfile(z):
                listing.append(z.split(path+'/')[1])
    return listing

# Note from GC
    # 1. If installing from python setup.py install, please `pip install cvxpy` first. 
    # 2. SCS requires numpy < 1.2 to be installed: https://github.com/bodono/scs-python/issues/32#issuecomment-802977693
    # 3. cvxpy fails to install on python 3.12, requires SCS, ECOS and OSQP: https://github.com/cvxpy/cvxpy/issues/1367
    # 4. qdldl (dependency of osqp has no wheel for python 3.12): https://github.com/cvxpy/cvxpy/issues/2269
    # 5. osqp < 0.6.2 does not depend on qdldl to install cvxpy

param = {
        'name': laml_libs.PROGRAM_NAME,
        'version': laml_libs.PROGRAM_VERSION,
        'description': laml_libs.PROGRAM_DESCRIPTION,
        'author': laml_libs.PROGRAM_AUTHOR,
        'license': laml_libs.PROGRAM_LICENSE,
        'packages': find_packages(),
        'package_data': {'':recursive_list_dir('laml_unit_tests')},
        'include_package_data': True,
        'scripts': ['run_laml.py', 'laml_tests.py', 'extra_laml_tests.py'],
        'zip_safe': True,
        'install_requires': ['scipy>=1.3.1', 'cvxpy>=1.4', 'treeswift>=1.1.39', 'Mosek>=10.1.16', 
                             'Biopython>=1.71','matplotlib>=2.2.2','setuptools',"jax>=0.4.30,<0.5",
                             "loguru>=0.7.3,<1.0", "networkx>=3.2.1,<4.0", "optax>=0.2.4,<1.0", "pandas>=2.2.3,<3.0"], 
        'keywords': 'Phylogenetics Evolution Computational Maximum-likelihood Lineage Tracing',
        'long_description': """LAML is a maximum likelihood algorithm under the Probabilistic Mixed-type Missing (PMM) model. Given a lineage tracing experiment character matrix with heterogeneous per-site alphabets and mutation probabilities, LAML finds a maximum likelihood tree topology and estimates parameters including time-resolved branch lengths, heritable silencing rate and non-heritable dropout probability.""",
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
        'entry_points': {
            'console_scripts': [
                'run_laml= run_laml:main',
             ],
        },
}

setup(**param)
