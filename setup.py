from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import glob
import os
import pkg_resources

from primalscreen import __version__, _program

setup(
    name='primalscreen',
    version='0.0.1',    
    description='primalscreen',
    url='https://github.com/Desperate-Dan/primalscreen',
    author='Daniel Maloney',
    author_email='dmaloney@ed.ac.uk',
    license='...',
    packages=find_packages(),
    scripts=['primalscreen/function_file.py'],
    install_requires=['biopython'
                      ],
    entry_points="""
    [console_scripts]
    {program} = primalscreen.command:main
    """.format(program = _program)
)
