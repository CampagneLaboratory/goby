from distutils.core import setup

setup(
    name='Goby',
    version='1.6',
    packages=['goby'],
    author='Fabien Campagne, Nyasha Chambwe, Kevin Dorff, Marko Srdanovic',
    author_email='fac2003@med.cornell.edu',
    scripts=['AlignmentStats.py', 'AlignmentToText.py', 'CompactToFasta.py', 'ReadsStats.py'],
    url='http://goby.campagnelab.org/',
    description='Goby is a next-gen data management framework designed to facilitate the implementation of efficient data analysis pipelines.',
    license='GNU General Public License (GPL)',
    long_description=open('README.txt').read(),
)
