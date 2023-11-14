from setuptools import setup, find_packages
from setuptools.command.install import install
import sysconfig
import os
import scmail_libs
from os import walk, listdir
from os.path import join, normpath, isfile

class PostInstallCommand(install):
    def run(self):
        install.run(self)
        self.execute(post_install, (self.install_lib,), msg="Running post install task")

def post_install(install_lib):

    # Get the installation paths
    install_paths = sysconfig.get_paths()
    
    # Use the 'scripts' path, which typically contains executable scripts
    install_path = install_paths['scripts']
    print('install_path', install_path)

    # Modify the user's shell profile file to include the install path in the PATH variable
    profile_file = os.path.expanduser('~/.bashrc')  # Change '~/.bashrc' to the appropriate profile file
    with open(profile_file, 'a') as f:
        print(f'export PATH="{install_path}:$PATH"\n')
        f.write(f'export PATH="{install_path}:$PATH"\n')

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
        'install_requires': ['numpy>=1.16', 'treeswift>=1.1.37', 'scipy>=1.3.1', 'cvxpy>=1.4', 'mosek>=10.1.16', 'cmake>=3.18.0'],
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
        'entry_points': {
            'console_scripts': [
                'run_scmail = run_scmail:main',
                'run_scmail_tests = scmail_tests:main',
             ],
        },
        #'cmdclass': {'install': PostInstallCommand},
}

setup(**param)
