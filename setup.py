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


param = {
        'name': laml_libs.PROGRAM_NAME,
        'version': laml_libs.PROGRAM_VERSION,
        'description': laml_libs.PROGRAM_DESCRIPTION,
        'author': laml_libs.PROGRAM_AUTHOR,
        'license': laml_libs.PROGRAM_LICENSE,
        'packages': find_packages(),
        'package_data': {'':recursive_list_dir('laml_unit_tests')},
        'include_package_data': True,
        'scripts': ['run_laml.py', 'laml_tests.py'],
        'zip_safe': True,
        #'setup_requires': ['numpy'],
        # SCS requires numpy < 1.2 to be installed: https://github.com/bodono/scs-python/issues/32#issuecomment-802977693
        # cvxpy fails to install on python 3.12, requires SCS, ECOS and OSQP: https://github.com/cvxpy/cvxpy/issues/1367
        # qdldl (dependency of osqp has no wheel for python 3.12): https://github.com/cvxpy/cvxpy/issues/2269
        # osqp < 0.6.2 does not depend on qdldl to install cvxpy
        # NOTE: if installing from python setup.py install, please `pip install cvxpy` first. 
        'install_requires': ['cvxpy', 'treeswift>=1.1.39', 'Mosek', 'setuptools'], #, 'treeswift', 'Mosek', 'setuptools', 'scipy==1.11.4', 'scs==3.2.4.post1', 'clarabel==0.7.1'], #'ecos==2.0.13'], #'osqp==0.6.5'], # 'qdldl==0.1.7.post0'], #osqp>=0.6.1', 'ecos>=2', 'clarabel>=0.5.0', 'scs>=3.0', 'cvxpy'],
        #'install_requires': ['numpy>=1.16', 'treeswift>=1.1.39', 'scipy>=1.3.1', 'cvxpy>=1.4', 'Mosek>=10.1.16', 'cmake>=3.18.0', 'setuptools', 'pybind11'], #'osqp==0.6.1'],
        'keywords': 'Phylogenetics Evolution Computational Maximum-likelihood Lineage Tracing',
        'long_description': """LAML is a maximum likelihood algorithm under the Probabilistic Mixed-type Missing (PMM) model. Given a lineage tracing experiment character matrix with heterogeneous per-site alphabets and mutation probabilities, LAML will find a maximum likelihood tree topology and estimate branch lengths as well as stochastic dropout and heritable silencing missing data rates.""",
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
