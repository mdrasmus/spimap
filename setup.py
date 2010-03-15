#!/usr/bin/env python
# 
# setup for SPIDIR library package
#
# use the following to install:
#   python setup.py install
#

from distutils.core import setup, Extension

VERSION = '2.0'

setup(
    name='spidir',
    version=VERSION,
    description='Species informed gene tree reconstruction',
    long_description = """
            """,
    author='Matt Rasmussen',
    author_email='rasmus@mit.edu',
    url='http://compbio.mit.edu/spimap/',
    download_url='http://compbio.mit.edu/pub/spimap/spimap-%s.tar.gz' % VERSION,
    
    classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Environment :: Console',
          'Intended Audience :: Developers',
          'Intended Audience :: Education',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License (GPL)',
          'Operating System :: POSIX',
          'Programming Language :: Python',
          'Programming Language :: C++',
          'Topic :: Education',
          ],
    
    package_dir = {'': 'python'},
    packages=['spidir', 'rasmus', 'compbio'],
    py_modules=[],
    scripts=[],
    #ext_modules=[
    #    Extension(
    #        '', 
    #        [],
    #        include_dirs=[],
    #        libraries=[]
    #        )]
    )


