#!/usr/bin/env python

"""
This is a setup script for Intervene: a toolfor intersection and visualization of multiple genomic region sets

This code is free software; you can redistribute it and/or modify it under the terms of the 
BSD License (see the file LICENSE.md included with the distribution).

@author: Aziz Khan
@email: aziz.khan@ncmm.uio.no
"""
import os
import re
from setuptools import setup, find_packages

with open('intervene/__init__.py') as f:
    VERSION = re.search(r"__version__ = '(.+)'", f.read()).group(1)


CLASSIFIERS = [
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: MIT License',
    'Operating System :: OS Independent',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
    'Programming Language :: Python :: 3.10',
    'Programming Language :: Python :: 3.11',
    'Programming Language :: Python :: 3.12',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Software Development :: Libraries :: Python Modules',
]

install_requires = [
    'pybedtools',
    'matplotlib',
    'pandas',
    'numpy',
    'scipy',
    'seaborn',
]

def readme(fname):
    with open(os.path.join(os.path.dirname(__file__), fname)) as f:
        return f.read()

setup(
    name="intervene",
    description="A tool for intersection and visualization of multiple gene or genomic region sets",
    version=VERSION,
    author="Aziz Khan",
    license='MIT',
    platforms='linux/unix',
    author_email="azez.khan@gmail.com",
    url="https://github.com/asntech/intervene",
    long_description=readme("README.rst"),
    package_dir={'intervene': 'intervene'},

    packages=['intervene',
        'intervene.modules',
        'intervene.modules.pairwise',
        'intervene.modules.venn',
        'intervene.modules.upset',
        'intervene.example_data',
        'intervene.example_data.dbSUPER_mm9',
        'intervene.example_data.ENCODE_hESC',
        'intervene.example_data.Gene_list'],

    scripts=['intervene/intervene',
                   ],
    package_data={'intervene': ['example_data/dbSUPER_mm9/*.bed', 'example_data/ENCODE_hESC/*.bed','example_data/Gene_list/*.txt',]},
    #package_data={'intervene': ['example_data/*']},
    include_package_data=True,
    install_requires = install_requires,
    python_requires='>=3.8',
    classifiers=CLASSIFIERS,
)
