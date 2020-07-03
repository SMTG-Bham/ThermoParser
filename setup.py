import setuptools

with open('README.md', 'r') as f:
    long_description=f.read()

def load_test_suite():
    test_loader = unittest.TestLoader()
    test_suite = test_loader.discover('tests', pattern='test*.py')
    return test_suite

setuptools.setup(
    name='ThermoPlotter',
    version='0.0.3.post1',
    author='Kieran B Spooner',
    description='A simple thermoelectrics plotting tool',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://test.pypi.org/project/ThermoPlotter/',
    packages=setuptools.find_packages(),
    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Visualization'],
    keywords='chemistry materials thermoelectric dft phonopy phono3py amset',
    test_suite='setup.load_test_suite',
    install_requires=['h5py', 'matplotlib',
                      'numpy', 'pymatgen', 'scipy', 'pyyaml'],
    python_requires='>=3')
