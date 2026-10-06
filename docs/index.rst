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

This package accepts both Python ``str`` and ``bytes`` values, but arguments
passed to the same function must use compatible string types.

For Python ``str`` input, the extension reads the CPython Unicode representation
using ``PyUnicode_GET_LENGTH`` and ``PyUnicode_DATA``. The resulting
Levenshtein sequence therefore follows the Unicode code-point sequence of the
Python string rather than its UTF-8 byte encoding.

For ``bytes`` input, the sequence elements are byte values instead.

This distinction matters because code-point distance is not the same as
grapheme-cluster distance. A user-perceived character can contain multiple
Unicode code points, and canonically equivalent text can also have different
underlying code-point sequences.

The extension does not automatically apply Unicode normalization or
grapheme-cluster segmentation before calculating the distance. Applications
that require those semantics should preprocess both inputs consistently before
comparison.

For additional background on how runtime string representation defines the
sequence supplied to edit-distance algorithms, see `Levenshtein implementations
and Unicode sequence units
<https://www.levenshtein.net/levenshtein-implementations>`_.

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
