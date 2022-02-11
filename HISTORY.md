## Changelog

### v0.18.1
#### Fixed
- Indel.normalized_similarity was broken in RapidFuzz v2.0.0 (see #20)

### v0.18.0
#### Fixed
* Fixed memory leak in error path of setratio
* Fixed out of bound reads due to uninitialized variable in median
  * e.g. quickmedian(["test", "teste"], [0, 0]) caused out of bound reads

#### Changed
* Use a faster editops implementation provided by RapidFuzz
* Reduce code duplication
* reuse implementations from rapidfuzz-cpp
* Transition to scikit-build 

### v0.17.0
* Removed support for Python 3.5

### v0.16.1
* Add support for RapidFuzz v1.9.*

### v0.16.0
* Add support for Python 3.10

### v0.15.0
* Update SequenceMatcher interface to support the autojunk parameter

### v0.14.0
* Drop Python 2 support
* Fixed free of non heap object due caused by zero offset on a heap object
* Fixed warnings about missing type conversions
* Fix segmentation fault in subtract_edit when incorrect input types are used
* Fixed unchecked memory allocations
* Implement distance/ratio/hamming/jaro/jaro_winkler
  using rapidfuzz instead of providing a own implementation
* Implement Wrapper for inverse/editops/opcodes/matching_blocks/subtract_edit/apply_edit
  using Cython to simplify support for new Python versions

### v0.13.0
* Maintainership passed to Max Bachmann
* use faster bitparallel implementations for distance and ratio
* avoid string copies in distance, ratio and hamming
* Fix usage of deprecated Unicode APIs in distance, ratio and hamming
* Fixed incorrect window size inside Jaro and Jaro-Winkler implementation
* Fixed incorrect exception messages
* Removed unused functions and compiler specific hacks
* Split the Python and C implementations to simplify building of
  the C library
* Fixed multiple bugs which prevented the use as C library, since some functions
  only got defined when compiling for Python
* Build and deliver python wheels for the library
* Fixed incorrect allocation size in lev_editops_matching_blocks and
  lev_opcodes_matching_blocks

### v0.12.1
* Fixed handling of numerous possible wraparounds in calculating the size
  of memory allocations; incorrect handling of which could cause denial
  of service or even possible remote code execution in previous versions
  of the library.

### v0.12.0
* Fixed a bug in StringMatcher.StringMatcher.get_matching_blocks /
  extract_editops for Python 3; now allow only `str` editops on
  both Python 2 and Python 3, for simpler and working code.
* Added documentation in the source distribution and in GIT
* Fixed the package layout: renamed the .so/.dll to _levenshtein,
  and made it reside inside a package, along with the StringMatcher
  class.
* Fixed spelling errors.

### v0.11.2
* Fixed a bug in setup.py: installation would fail on Python 3 if the locale
  did not specify UTF-8 charset (Felix Yan).

* Added COPYING, StringMatcher.py, gendoc.sh and NEWS in MANIFEST.in, as they
  were missing from source distributions.

### v0.11.1
* Added Levenshtein.h to MANIFEST.in

### v0.11.0
* Python 3 support, maintainership passed to Antti Haapala

### v0.10.2
* Made python-Lehvenstein Git compatible and use setuptools for PyPi upload
* Created HISTORY.txt and made README reST compatible

### v0.10.1
* apply_edit() broken for Unicodes was fixed (thanks to Radovan Garabik)
* subtract_edit() function was added

### v0.10.0
* Hamming distance, Jaro similarity metric and Jaro-Winkler similarity
      metric were added
* ValueErrors raised on wrong argument types were fixed to TypeErrors

### v0.9.0
* a poor-but-fast generalized median method quickmedian() was added
* some auxiliary functions added to the C api (lev_set_median_index,
      lev_editops_normalize, ...)

### v0.8.2
* fixed missing `static' in the method list

### v0.8.1
* some compilation problems with non-gcc were fixed

v0.8.0
* median_improve(), a generalized median improving function, was added
* an arbitrary length limitation imposed on greedy median() result was
      removed
* out of memory should be handled more gracefully (on systems w/o memory
      overcomitting)
* the documentation now passes doctest

### v0.7.0
* fixed greedy median() for Unicode characters > U+FFFF, it's now usable
      with whatever integer type wchar_t happens to be
* added missing MANIFEST
* renamed exported C functions, all public names now have lev_, LEV_ or
      Lev prefix; defined lev_byte, lev_wchar, and otherwise santinized
      the (still unstable) C interface
* added edit-ops group of functions, with two interfaces: native, useful
      for string averaging, and difflib-like for interoperability
* added an example SequenceMatcher-like class StringMatcher

### v0.6.0
* a segfault in seqratio()/setratio() on invalid input has been fixed
      to an exception
* optimized ratio() and distance() (about 20%)
* Levenshtein.h header file was added to make it easier to actually use
      it as a C library

### v0.5.0
* a segfault in setratio() was fixed
* median() handles all empty strings situation more gracefully

### v0.4.0
* new functions seqratio() and setratio() computing similarity between
      string sequences and sets
* Levenshtein optimizations (affects all routines except median())
* all Sequence objects are accepted, not just Lists

### v0.3.0
* setmedian() finding set median was added
* median() initial overhead for Unicodes was reduced

### v0.2.0
* ratio() and distance() now accept both Strings and Unicodes
* removed uratio() and udistance()
* Levenshtein.c is now compilable as a C library (with -DNO_PYTHON)
* a median() function finding approximate weighted median of a string
      set was added

### v0.1.0
* Inital release

