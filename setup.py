# pmx  Copyright Notice
# ============================
#
# The pmx source code is copyrighted, but you can freely use and
# copy it as long as you don't change or remove any of the copyright
# notices.
#
# ----------------------------------------------------------------------
# pmx is Copyright (C) 2006-2013 by Daniel Seeliger
#
#                        All Rights Reserved
#
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name of Daniel Seeliger not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
#
# DANIEL SEELIGER DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
# SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS.  IN NO EVENT SHALL DANIEL SEELIGER BE LIABLE FOR ANY
# SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
# RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
# CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
# ----------------------------------------------------------------------

from os import path
# Always prefer setuptools over distutils
from setuptools import setup
# For extension managements
from setuptools.extension import Extension
from Cython.Build import cythonize

here = path.abspath(path.dirname(__file__))

pmx = Extension('pmx/_pmx',
                libraries=['m'],
                include_dirs=['src/pmx'],
                sources=['src/pmx/Geometry.c', 'src/pmx/wrap_Geometry.c',
                         'src/pmx/init.c', 'src/pmx/Energy.c']
                )

xdrio = Extension('pmx/_xdrio',
                  libraries=['m'],
                  include_dirs=['src/xdr'],
                  sources=['src/xdr/xdrfile.c', 'src/xdr/xdrfile_trr.c',
                           'src/xdr/xdrfile_xtc.c']
                  )

## cpp_test = Extension('pmx/_cpp_test',
##                    libraries = ['m'],
##                    include_dirs = ['src'],
##                    sources = ['src/cpp_test.cpp']

##                  )
extenstions = [pmx, xdrio]

with open('README.rst', 'r') as f:
    long_description = f.read()

setup(
    name='pmx',
    version='1.1.1dev',
    description='Python Toolbox structure file editing and \
        writing simulation setup/analysis tools',
    author='Arthur Zalevsky <aozalevsky.fbb.msu.ru>, \
        Daniel Seeliger <seeliger.biosoft@gmail.de>',
    author_email='aozalevsky@fbb.msu.ru',
    url='https://github.com/dseeliger/pmx/',
    long_description=long_description,
    packages=['pmx'],
    ext_modules=cythonize(extenstions)
)
