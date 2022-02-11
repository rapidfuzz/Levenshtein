"""
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
"""

__author__ = "Max Bachmann"
__license__ = "GPL"
__version__ = "0.18.1"

from rapidfuzz.distance.Levenshtein import distance
from rapidfuzz.distance.Indel import normalized_similarity as ratio
from rapidfuzz.distance.Hamming import distance as hamming
from rapidfuzz.distance.Jaro import similarity as jaro
from rapidfuzz.distance.JaroWinkler import similarity as jaro_winkler

from Levenshtein._levenshtein import (
    median,
    median_improve,
    quickmedian,
    setmedian,
    seqratio,
    setratio
)

from Levenshtein.cpp_levenshtein import (
    inverse,
    editops,
    opcodes,
    matching_blocks,
    subtract_edit,
    apply_edit
)
