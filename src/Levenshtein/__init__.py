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
__version__ = "0.16.0"

from rapidfuzz import string_metric as _string_metric

from Levenshtein._levenshtein import (
    median,
    median_improve,
    quickmedian,
    setmedian,
    seqratio,
    setratio
)

from Levenshtein.c_levenshtein import (
    inverse,
    editops,
    opcodes,
    matching_blocks,
    subtract_edit,
    apply_edit
)

def distance(string1, string2):
    """
    Compute absolute Levenshtein distance of two strings.

    Parameters
    ----------
    string1 : str
        First string to compare.
    string2 : str
        Second string to compare.

    Returns
    -------
    distance : int
        distance between s1 and s2

    Examples
    --------
    it's hard to spell Levenshtein correctly
    
    >>> distance('Levenshtein', 'Lenvinsten')
    4
    >>> distance('Levenshtein', 'Levensthein')
    2
    >>> distance('Levenshtein', 'Levenshten')
    1
    >>> distance('Levenshtein', 'Levenshtein')
    0
    
    Yeah, we've managed it at last.
    """
    return _string_metric.levenshtein(string1, string2)

def ratio(string1, string2):
    """
    Compute similarity of two strings using the InDel distance.
    The InDel distance is similar to the Levenshtein distance with the following
    weights:
    - Insertion: 1
    - Deletion: 1
    - Substitution: 2
    
    Parameters
    ----------
    string1 : str
        First string to compare.
    string2 : str
        Second string to compare.

    Returns
    -------
    similarity : float
        similarity between s1 and s2 as a float between 0 and 1

    Notes
    -----
    The similarity is usually equal or somewhat higher than
    difflib.SequenceMatcher.ratio(), because it's
    based on the real minimal edit distance.
    
    Examples
    --------
    >>> ratio('Hello world!', 'Holly grail!')
    0.583333
    >>> ratio('Brian', 'Jesus')
    0.0
    
    Really?  I thought there was some similarity.
    """
    return _string_metric.normalized_levenshtein(string1, string2, weights=(1,1,2)) / 100

def hamming(string1, string2):
    """
    Compute Hamming distance of two strings.
    The Hamming distance is simply the number of differing characters.
    That means the length of the strings must be the same.

    Parameters
    ----------
    string1 : str
        First string to compare.
    string2 : str
        Second string to compare.

    Returns
    -------
    distance : int
        distance between s1 and s2

    Raises
    ------
    ValueError
        If s1 and s2 have a different length

    Examples
    --------
    >>> hamming('Hello world!', 'Holly grail!')
    7
    >>> hamming('Brian', 'Jesus')
    5
    """
    return _string_metric.hamming(string1, string2)

def jaro(string1, string2):
    """
    Compute Jaro string similarity metric of two strings.
    The Jaro string similarity metric is intended for short strings like
    personal last names.
    
    Parameters
    ----------
    string1 : str
        First string to compare.
    string2 : str
        Second string to compare.

    Returns
    -------
    similarity : float
        similarity between s1 and s2 as a float between 0 and 1.
    
    Examples
    --------
    >>> jaro('Brian', 'Jesus')
    0.0
    >>> jaro('Thorkel', 'Thorgier')  # doctest: +ELLIPSIS
    0.779761...
    >>> jaro('Dinsdale', 'D')  # doctest: +ELLIPSIS
    0.708333...
    """
    return _string_metric.jaro_similarity(string1, string2) / 100

def jaro_winkler(string1, string2, prefix_weight=0.1):
    """
    Compute Jaro string similarity metric of two strings.
    The Jaro-Winkler string similarity metric is a modification of Jaro
    metric giving more weight to common prefix, as spelling mistakes are
    more likely to occur near ends of words.
    
    The prefix weight is inverse value of common prefix length sufficient
    to consider the strings *identical*.  If no prefix weight is
    specified, 1/10 is used.

    
    Parameters
    ----------
    string1 : str
        First string to compare.
    string2 : str
        Second string to compare.
    prefix_weight : float, optional
        weight used for the common prefix of the two strings.
        Has to be between 0 and 0.25. Default is 0.1.

    Returns
    -------
    similarity : float
        similarity between s1 and s2 as a float between 0 and 1.

    Raises
    ------
    ValueError
        If prefix_weight is invalid

    Examples
    --------
    >>> jaro_winkler('Brian', 'Jesus')
    0.0
    >>> jaro_winkler('Thorkel', 'Thorgier')
    0.867857
    >>> jaro_winkler('Dinsdale', 'D')
    0.7375
    
    When using a prefix_weight of 0.25 any all strings with a common prefix of 4
    will have a similarity of 1.0

    >>> jaro_winkler('Thorkel', 'Thorgier', 0.25)
    1.0
    """
    return _string_metric.jaro_winkler_similarity(string1, string2, prefix_weight=prefix_weight) / 100
