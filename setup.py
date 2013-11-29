from setuptools import setup, find_packages
import os

from distutils.core import Extension

version = '0.11.0'

extLevensthein = Extension('Levenshtein',
                           sources = ['Levenshtein.c'],
                           )

setup(name='python-Levenshtein',
      version=version,
      description="Python extension for computing string edit distances and similarities.",
      long_description=open("README.rst").read() + "\n" +
                       open(os.path.join("HISTORY.txt")).read(),
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
      packages=find_packages(exclude=['ez_setup']),
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
