import os
import sys
try:
    import ez_setup
    ez_setup.use_setuptools()
except ImportError:
    pass

from distutils.core import setup

setup(
    name='goby',
    version='1.9.8.3',
    packages=['goby'],
    author='Campagne Lab',
    author_email='icb@med.cornell.edu',
    scripts=[
        'GobyAlignmentStats.py',
        'GobyAlignmentToText.py',
        'GobyCompactToFasta.py',
        'GobyReadsStats.py'
        ],
    url='http://goby.campagnelab.org/',
    description='Python API for reading binary data files created with the Goby next-gen data management framework.',
    license='GNU General Public License (GPL)',
    long_description=open('README.txt').read(),
    requires=[
        'google.protobuf (>=2.3)',
        'pyjavaproperties (>=0.3)',
        ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
        ],
)
