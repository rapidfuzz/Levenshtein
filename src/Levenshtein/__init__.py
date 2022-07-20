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

__author__: str = "Max Bachmann"
__license__: str = "GPL"
__version__: str = "0.19.3"

from rapidfuzz.distance.Levenshtein import distance
from rapidfuzz.distance.Indel import normalized_similarity as ratio
from rapidfuzz.distance.Hamming import distance as hamming
from rapidfuzz.distance.Jaro import similarity as jaro
from rapidfuzz.distance.JaroWinkler import similarity as jaro_winkler
from rapidfuzz.distance.Levenshtein import (
    editops as _editops,
    opcodes as _opcodes,
)
from rapidfuzz.distance import (
    Editops as _Editops,
    Opcodes as _Opcodes,
)

from Levenshtein.levenshtein_cpp import (
    quickmedian,
    median,
    median_improve,
    setmedian,
    setratio,
    seqratio,
)


def editops(*args):
    """
    Find sequence of edit operations transforming one string to another.

    editops(source_string, destination_string)
    editops(edit_operations, source_length, destination_length)

    The result is a list of triples (operation, spos, dpos), where
    operation is one of 'equal', 'replace', 'insert', or 'delete';  spos
    and dpos are position of characters in the first (source) and the
    second (destination) strings.  These are operations on signle
    characters.  In fact the returned list doesn't contain the 'equal',
    but all the related functions accept both lists with and without
    'equal's.

    Examples
    --------
    >>> editops('spam', 'park')
    [('delete', 0, 0), ('insert', 3, 2), ('replace', 3, 3)]

    The alternate form editops(opcodes, source_string, destination_string)
    can be used for conversion from opcodes (5-tuples) to editops (you can
    pass strings or their lengths, it doesn't matter).
    """
    # convert: we were called (bops, s1, s2)
    if len(args) == 3:
        arg1, arg2, arg3 = args
        len1 = arg2 if isinstance(arg2, int) else len(arg2)
        len2 = arg3 if isinstance(arg3, int) else len(arg3)
        return _Editops(arg1, len1, len2)

    # find editops: we were called (s1, s2)
    arg1, arg2 = args
    return _editops(arg1, arg2)


def opcodes(*args):
    """
    Find sequence of edit operations transforming one string to another.

    opcodes(source_string, destination_string)
    opcodes(edit_operations, source_length, destination_length)

    The result is a list of 5-tuples with the same meaning as in
    SequenceMatcher's get_opcodes() output.  But since the algorithms
    differ, the actual sequences from Levenshtein and SequenceMatcher
    may differ too.

    Examples
    --------
    >>> for x in opcodes('spam', 'park'):
    ...     print(x)
    ...
    ('delete', 0, 1, 0, 0)
    ('equal', 1, 3, 0, 2)
    ('insert', 3, 3, 2, 3)
    ('replace', 3, 4, 3, 4)

    The alternate form opcodes(editops, source_string, destination_string)
    can be used for conversion from editops (triples) to opcodes (you can
    pass strings or their lengths, it doesn't matter).
    """
    # convert: we were called (ops, s1, s2)
    if len(args) == 3:
        arg1, arg2, arg3 = args
        len1 = arg2 if isinstance(arg2, int) else len(arg2)
        len2 = arg3 if isinstance(arg3, int) else len(arg3)
        return _Opcodes(arg1, len1, len2)

    # find editops: we were called (s1, s2)
    arg1, arg2 = args
    return _opcodes(arg1, arg2)


def inverse(edit_operations):
    """
    Invert the sense of an edit operation sequence.

    In other words, it returns a list of edit operations transforming the
    second (destination) string to the first (source).  It can be used
    with both editops and opcodes.

    Parameters
    ----------
    edit_operations : list[]
        edit operations to invert

    Returns
    -------
    edit_operations : list[]
        inverted edit operations

    Examples
    --------
    >>> inverse(editops('spam', 'park'))
    [('insert', 0, 0), ('delete', 2, 3), ('replace', 3, 3)]
    >>> editops('park', 'spam')
    [('insert', 0, 0), ('delete', 2, 3), ('replace', 3, 3)]
    """
    return edit_operations.inverse()


def subtract_edit(edit_operations, subsequence):
    """
    Subtract an edit subsequence from a sequence.
    
    subtract_edit(edit_operations, subsequence)
    
    The result is equivalent to
    editops(apply_edit(subsequence, s1, s2), s2), except that is
    constructed directly from the edit operations.  That is, if you apply
    it to the result of subsequence application, you get the same final
    string as from application complete edit_operations.  It may be not
    identical, though (in amibuous cases, like insertion of a character
    next to the same character).
    
    The subtracted subsequence must be an ordered subset of
    edit_operations.
    
    Note this function does not accept difflib-style opcodes as no one in
    his right mind wants to create subsequences from them.

    Examples
    --------
    >>> e = editops('man', 'scotsman')
    >>> e1 = e[:3]
    >>> bastard = apply_edit(e1, 'man', 'scotsman')
    >>> bastard
    'scoman'
    >>> apply_edit(subtract_edit(e, e1), bastard, 'scotsman')
    'scotsman'
    """
    return edit_operations - subsequence


def apply_edit(edit_operations, source_string, destination_string):
    """
    Apply a sequence of edit operations to a string.
    
    apply_edit(edit_operations, source_string, destination_string)
    
    In the case of editops, the sequence can be arbitrary ordered subset
    of the edit sequence transforming source_string to destination_string.
    
    Examples
    --------
    >>> e = editops('man', 'scotsman')
    >>> apply_edit(e, 'man', 'scotsman')
    'scotsman'
    >>> apply_edit(e[:3], 'man', 'scotsman')
    'scoman'
    
    The other form of edit operations, opcodes, is not very suitable for
    such a tricks, because it has to always span over complete strings,
    subsets can be created by carefully replacing blocks with 'equal'
    blocks, or by enlarging 'equal' block at the expense of other blocks
    and adjusting the other blocks accordingly.

    >>> a, b = 'spam and eggs', 'foo and bar'
    >>> e = opcodes(a, b)
    >>> apply_edit(inverse(e), b, a)
    'spam and eggs'
    >>> e[4] = ('equal', 10, 13, 8, 11)
    >>> e
    [('delete', 0, 1, 0, 0), ('replace', 1, 4, 0, 3), ('equal', 4, 9, 3, 8), ('delete', 9, 10, 8, 8), ('equal', 10, 13, 8, 11)]
    >>> apply_edit(e, a, b)
    'foo and ggs'
    """
    return edit_operations.apply(source_string, destination_string)


def matching_blocks(edit_operations, source_string, destination_string):
    """
    Find identical blocks in two strings.

    Parameters
    ----------
    edit_operations : list[]
        editops or opcodes created for the source and destination string
    source_string : str | int
        source string or the length of the source string
    destination_string : str | int
        destination string or the length of the destination string

    Returns
    -------
    matching_blocks : list[]
        List of triples with the same meaning as in SequenceMatcher's
        get_matching_blocks() output.

    Examples
    --------
    >>> a, b = 'spam', 'park'
    >>> matching_blocks(editops(a, b), a, b)
    [(1, 0, 2), (4, 4, 0)]
    >>> matching_blocks(editops(a, b), len(a), len(b))
    [(1, 0, 2), (4, 4, 0)]
    
    The last zero-length block is not an error, but it's there for
    compatibility with difflib which always emits it.
    
    One can join the matching blocks to get two identical strings:

    >>> a, b = 'dog kennels', 'mattresses'
    >>> mb = matching_blocks(editops(a,b), a, b)
    >>> ''.join([a[x[0]:x[0]+x[2]] for x in mb])
    'ees'
    >>> ''.join([b[x[1]:x[1]+x[2]] for x in mb])
    'ees'
    """
    return _MatchingBlocks(edit_operations, source_string, destination_string)
