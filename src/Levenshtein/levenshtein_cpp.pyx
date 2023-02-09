# distutils: language=c++
# cython: language_level=3
# cython: binding=True

from libc.stdint cimport uint32_t
from libcpp.vector cimport vector
from libcpp cimport bool
from libcpp.utility cimport move

cdef extern from *:
    int PyUnicode_4BYTE_KIND
    object PyUnicode_FromKindAndData(int kind, const void *buffer, Py_ssize_t size)

cdef extern from "<string>" namespace "std" nogil:
    cdef cppclass basic_string[T]:
        ctypedef size_t size_type

        basic_string() except +

        void resize(size_type) except +

        T& operator[](size_type)

        const T* data()
        size_type size()

cdef extern from "_levenshtein.hpp":
    cdef basic_string[uint32_t] lev_greedy_median(const vector[RF_String]& strings, const vector[double]& weights) except +
    cdef basic_string[uint32_t] lev_median_improve(const RF_String& string, const vector[RF_String]& strings, const vector[double]& weights) except +
    cdef basic_string[uint32_t] lev_quick_median(const vector[RF_String]& strings, const vector[double]& weights) except +
    cdef basic_string[uint32_t] lev_set_median(const vector[RF_String]& strings, const vector[double]& weights) except +

    cdef double lev_set_distance(const vector[RF_String]& strings1, const vector[RF_String]& strings2) except +
    cdef double lev_edit_seq_distance(const vector[RF_String]& strings1, const vector[RF_String]& strings2) except +

    ctypedef struct RF_String:
        pass

    cdef bool is_valid_string(object)
    cdef RF_String convert_string(object)


cdef inline RF_String conv_sequence(seq) except *:
    if is_valid_string(seq):
        return convert_string(seq)
    raise TypeError("Expected string or bytes")

cdef vector[double] extract_weightlist(wlist, size_t n) except *:
    cdef size_t i
    cdef double weight
    cdef vector[double] weights

    if wlist is None:
        weights.resize(n, 1.0)
    else:
        weights.resize(n)
        for i, w in enumerate(wlist):
            weight = w
            if w < 0:
                raise ValueError(f"weight {weight} is negative")
            weights[i] = w
    return weights

cdef vector[RF_String] extract_stringlist(strings) except *:
    cdef vector[RF_String] strlist

    for string in strings:
        strlist.push_back(move(conv_sequence(string)))

    return move(strlist)

def median(strlist, wlist = None, *):
    """
    Find an approximate generalized median string using greedy algorithm.

    You can optionally pass a weight for each string as the second
    argument. The weights are interpreted as item multiplicities,
    although any non-negative real numbers are accepted. Use them to
    improve computation speed when strings often appear multiple times
    in the sequence.

    Examples
    --------

    >>> median(['SpSm', 'mpamm', 'Spam', 'Spa', 'Sua', 'hSam'])
    'Spam'
    >>> fixme = ['Levnhtein', 'Leveshein', 'Leenshten',
    ...          'Leveshtei', 'Lenshtein', 'Lvenstein',
    ...          'Levenhtin', 'evenshtei']
    >>> median(fixme)
    'Levenshtein'

    Hm. Even a computer program can spell Levenshtein better than me.
    """
    if wlist is not None and len(strlist) != len(wlist):
        raise ValueError("strlist has a different length than wlist")

    weights = extract_weightlist(wlist, len(strlist))
    strings = extract_stringlist(strlist)
    median = lev_greedy_median(strings, weights)
    return PyUnicode_FromKindAndData(PyUnicode_4BYTE_KIND, median.data(), median.size())

def quickmedian(strlist, wlist = None, *):
    """
    Find a very approximate generalized median string, but fast.

    See median() for argument description.

    This method is somewhere between setmedian() and picking a random
    string from the set; both speedwise and quality-wise.

    Examples
    --------

    >>> fixme = ['Levnhtein', 'Leveshein', 'Leenshten',
    ...         'Leveshtei', 'Lenshtein', 'Lvenstein',
    ...         'Levenhtin', 'evenshtei']
    >>> quickmedian(fixme)
    'Levnshein'
    """
    if wlist is not None and len(strlist) != len(wlist):
        raise ValueError("strlist has a different length than wlist")

    weights = extract_weightlist(wlist, len(strlist))
    strings = extract_stringlist(strlist)
    median = lev_quick_median(strings, weights)
    return PyUnicode_FromKindAndData(PyUnicode_4BYTE_KIND, median.data(), median.size())

def median_improve(string, strlist, wlist = None, *):
    """
    Improve an approximate generalized median string by perturbations.

    The first argument is the estimated generalized median string you
    want to improve, the others are the same as in median(). It returns
    a string with total distance less or equal to that of the given string.

    Note this is much slower than median(). Also note it performs only
    one improvement step, calling median_improve() again on the result
    may improve it further, though this is unlikely to happen unless the
    given string was not very similar to the actual generalized median.

    Examples
    --------

    >>> fixme = ['Levnhtein', 'Leveshein', 'Leenshten',
    ...          'Leveshtei', 'Lenshtein', 'Lvenstein',
    ...          'Levenhtin', 'evenshtei']
    >>> median_improve('spam', fixme)
    'enhtein'
    >>> median_improve(median_improve('spam', fixme), fixme)
    'Levenshtein'

    It takes some work to change spam to Levenshtein.
    """
    if wlist is not None and len(strlist) != len(wlist):
        raise ValueError("strlist has a different length than wlist")

    weights = extract_weightlist(wlist, len(strlist))
    query = conv_sequence(string)
    strings = extract_stringlist(strlist)
    median = lev_median_improve(query, strings, weights)
    return PyUnicode_FromKindAndData(PyUnicode_4BYTE_KIND, median.data(), median.size())

def setmedian(strlist, wlist = None, *):
    """
    Find set median of a string set (passed as a sequence).

    See median() for argument description.

    The returned string is always one of the strings in the sequence.

    Examples
    --------

    >>> setmedian(['ehee', 'cceaes', 'chees', 'chreesc',
    ...            'chees', 'cheesee', 'cseese', 'chetese'])
    'chees'

    You haven't asked me about Limburger, sir.
    """

    if wlist is not None and len(strlist) != len(wlist):
        raise ValueError("strlist has a different length than wlist")

    weights = extract_weightlist(wlist, len(strlist))
    strings = extract_stringlist(strlist)
    median = lev_set_median(strings, weights)
    return PyUnicode_FromKindAndData(PyUnicode_4BYTE_KIND, median.data(), median.size())

def setratio(strlist1, strlist2, *):
    """
    Compute similarity ratio of two strings sets (passed as sequences).

    The best match between any strings in the first set and the second
    set (passed as sequences) is attempted. I.e., the order doesn't
    matter here.

    Examples
    --------

    >>> setratio(['newspaper', 'litter bin', 'tinny', 'antelope'],
    ...          ['caribou', 'sausage', 'gorn', 'woody'])
    0.281845

    No, even reordering doesn't help the tinny words to match the
    woody ones.
    """

    strings1 = extract_stringlist(strlist1)
    strings2 = extract_stringlist(strlist2)
    lensum = strings1.size() + strings2.size()

    if lensum == 0:
        return 1.0

    if strings1.empty():
        dist = <double>strings2.size()
    elif strings2.empty():
        dist = <double>strings1.size()
    else:
        dist = lev_set_distance(strings1, strings2)

    return (<double>lensum - dist) / <double>lensum

def seqratio(strlist1, strlist2, *):
    """
    Compute similarity ratio of two sequences of strings.

    This is like ratio(), but for string sequences. A kind of ratio()
    is used to to measure the cost of item change operation for the
    strings.

    Examples
    --------

    >>> seqratio(['newspaper', 'litter bin', 'tinny', 'antelope'],
    ...          ['caribou', 'sausage', 'gorn', 'woody'])
    0.215178
    """

    strings1 = extract_stringlist(strlist1)
    strings2 = extract_stringlist(strlist2)
    lensum = strings1.size() + strings2.size()

    if lensum == 0:
        return 1.0

    if strings1.empty():
        dist = <double>strings2.size()
    elif strings2.empty():
        dist = <double>strings1.size()
    else:
        dist = lev_edit_seq_distance(strings1, strings2)

    return (<double>lensum - dist) / <double>lensum
