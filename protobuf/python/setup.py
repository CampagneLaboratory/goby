from distutils.core import setup

setup(
    name='Goby',
    version='1.6',
    packages=['goby'],
    author='Fabien Campagne',
    author_email='fac2003@med.cornell.edu',
    scripts=['GobyAlignmentStats.py', 'GobyAlignmentToText.py', 'GobyCompactToFasta.py', 'GobyReadsStats.py'],
    url='http://goby.campagnelab.org/',
    description='Goby is a next-gen data management framework designed to facilitate the implementation of efficient data analysis pipelines.',
    license='GNU General Public License (GPL)',
    long_description=open('README.txt').read(),
    requires=['pyjavaproperties (>=0.3)', 'google.protobuf (>=2.3)'],
    classifiers=['Development Status :: 4 - Beta', 'Intended Audience :: Developers', 'License :: OSI Approved :: GNU General Public License (GPL)', 'Topic :: Scientific/Engineering :: Bio-Informatics'],
)
