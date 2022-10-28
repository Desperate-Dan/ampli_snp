from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import glob
import os
import pkg_resources

from ampli_snp import __version__, _program

setup(
    name='ampli_snp',
    version='0.0.1',    
    description='ampli_snp',
    url='https://github.com/Desperate-Dan/ampli_snp',
    author='Daniel Maloney',
    author_email='dmaloney@ed.ac.uk',
    license='...',
    packages=find_packages(),
    scripts=['ampli_snp/function_file.py'],
    install_requires=['biopython'
                      ],
    entry_points="""
    [console_scripts]
    {program} = ampli_snp.command:main
    """.format(program = _program)
)
