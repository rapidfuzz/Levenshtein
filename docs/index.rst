.. RapidFuzz documentation master file, created by
   sphinx-quickstart on Fri Jan  1 19:02:29 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Levenshtein's documentation!
=======================================

A C extension module for fast computation of:

- Levenshtein (edit) distance and edit sequence manipulation
- string similarity
- approximate median strings, and generally string averaging
- string sequence and set similarity

Levenshtein has a some overlap with difflib (SequenceMatcher).  It
supports only strings, not arbitrary sequence types, but on the
other hand it's much faster.

It supports both normal and Unicode strings, but can't mix them, all
arguments to a function (method) have to be of the same type (or its
subclasses).

.. toctree::
   :maxdepth: 2
   :caption: Installation:

   installation

.. toctree::
   :maxdepth: 2
   :caption: Usage:

   levenshtein

.. toctree::
   :maxdepth: 2
   :caption: Changelog:

   changelog

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
