#!/usr/bin/env python3

__name__ =       'ThermoPlotter'
__author__ =     'Kieran B. Spooner'
__copyright__ =  'Copyright Scanlon Materials Theory Group (2020)'
__version__ =    '1.0.0'
__maintainer__ = 'Kieran B. Spooner'
__email__ =      'kieran.spooner.14@ucl.ac.uk'
__date__ =       'February 4th 2021'

from glob import glob
import os
import setuptools
from setuptools.command.install import install
from shutil import copy

with open('README.rst', 'r') as f:
    long_description=f.read()

def load_test_suite():
    test_loader = unittest.TestLoader()
    test_suite = test_loader.discover('tests', pattern='test*.py')
    return test_suite

def install_style():
    from matplotlib import get_configdir
    style = ['tp.mplstyle', 'tp-large.mplstyle']
    styledir = os.path.join(get_configdir(), 'stylelib')
    if not os.path.exists(styledir):
        os.makedirs(styledir)

    for plotting in style:
        copy(os.path.join(os.path.dirname(__file__), plotting),
             os.path.join(styledir, plotting))

class PostInstallMoveFile(install):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        install_style()

setuptools.setup(
    name='tp',
    version=__version__,
    author=__author__,
    description='streamlined analysis of thermoelectric properties',
    long_description=long_description,
    long_description_content_type='text/x-rst',
    url='https://smtg-ucl.github.io/ThermoPlotter/',
    packages=setuptools.find_packages(),
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Visualization'],
    keywords='chemistry materials thermoelectric dft phonopy phono3py '
             'amset boltztrap tp te matplotlib',
    test_suite='setup.load_test_suite',
    install_requires=['click', 'h5py', 'matplotlib', 'numpy', 'pymatgen',
                      'pyyaml', 'scipy'],
    python_requires='>=3.6',
    cmdclass={'install': PostInstallMoveFile},
    entry_points={'console_scripts': ['tp = tp.cli.cli:tp_cli']})
