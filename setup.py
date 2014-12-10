from setuptools import setup
import os
import sys

from distutils.core import Extension

version = '0.12.0'

extLevensthein = Extension('Levenshtein._levenshtein',
                           sources = ['Levenshtein/_levenshtein.c'],
                           )

if sys.version_info >= (3, 0):
    _open = lambda f: open(f, encoding='utf8')
else:
    _open = open

setup(name='python-Levenshtein',
      version=version,
      description="Python extension for computing string edit distances and similarities.",
      long_description=_open("README.rst").read() + "\n" +
                       _open(os.path.join("HISTORY.txt")).read(),
      # Get more strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      classifiers=[
        "License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: Implementation :: CPython"
        ],
      keywords='string Levenshtein comparison edit-distance',
      author='Antti Haapala',
      author_email='antti@haapala.name',
      url='http://github.com/ztane/python-Levenshtein',
      license='GPL',
      packages=['Levenshtein'],
      namespace_packages=[],
      include_package_data=True,
      zip_safe=False,
      ext_modules = [extLevensthein],
      install_requires=[
          'setuptools',
          # -*- Extra requirements: -*-
      ],
      entry_points="""
      """,
      )
