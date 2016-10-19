.. contents ::

Introduction
------------

The Levenshtein Python C extension module contains functions for fast
computation of

* Levenshtein (edit) distance, and edit operations

* string similarity

* approximate median strings, and generally string averaging

* string sequence and set similarity

It supports both normal and Unicode strings.

Python 2.2 or newer is required; Python 3 is supported.

StringMatcher.py is an example SequenceMatcher-like class built on the top of
Levenshtein.  It misses some SequenceMatcher's functionality, and has some
extra OTOH.

Levenshtein.c can be used as a pure C library, too.  You only have to define
NO_PYTHON preprocessor symbol (-DNO_PYTHON) when compiling it.  The
functionality is similar to that of the Python extension.  No separate docs
are provided yet, RTFS.  But they are not interchangeable:

* C functions exported when compiling with -DNO_PYTHON (see Levenshtein.h)
  are not exported when compiling as a Python extension (and vice versa)

* Unicode character type used with -DNO_PYTHON is wchar_t, Python extension
  uses Py_UNICODE, they may be the same but don't count on it

Documentation
--------------

* `Documentation for the current version <https://rawgit.com/ztane/python-Levenshtein/master/docs/Levenshtein.html>`_

gendoc.sh generates HTML API documentation,
you probably want a selfcontained instead of includable version, so run
in ``./gendoc.sh --selfcontained``.  It needs Levenshtein already installed
and genextdoc.py.

License
-------

Levenshtein is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your option)
any later version.

See the file COPYING for the full text of GNU General Public License version 2.

History
-------

This package was long missing from the Python Package Index and available as source checkout only, but can now `be found on PyPI again <https://pypi.python.org/pypi/python-Levenshtein>`_.

We needed to restore this package for `Go Mobile for Plone <http://webandmobile.mfabrik.com>`_
and `Pywurfl <http://celljam.net/>`_ projects which depend on this.

Source code
-----------

* http://github.com/ztane/python-Levenshtein/

Authors
-------

* Maintainer: `Antti Haapala <antti@haapala.name>`

* Python 3 compatibility: Esa Määttä

* Jonatas CD: Fixed documentation generation

* Previous maintainer: `Mikko Ohtamaa <http://opensourcehacker.com>`_

* Original code: David Necas (Yeti) <yeti at physics.muni.cz>
