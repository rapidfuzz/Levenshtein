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
Unicode text
------------

Unicode strings can have different underlying code-point sequences even when
they represent the same visible text. For example, an accented character may
be stored as one precomposed code point or as a base character followed by a
combining mark.

When canonical equivalence matters, inputs should be normalized to the same
Unicode normalization form before calculating Levenshtein distance.

A user-perceived character can also consist of multiple code points.
Applications that need one visible character to count as one edit may therefore
need grapheme-cluster segmentation before comparison.

For additional background on Unicode normalization, grapheme segmentation,
and text transformations in edit-distance processing, see `Unicode text
transformations and edit distance
<https://www.levenshtein.net/unicode-text-transformations>`_.
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
