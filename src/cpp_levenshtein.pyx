# distutils: language=c++
# cython: language_level=3
# cython: binding=True

cdef extern from "cpp_string_metric.hpp":
    void validate_string(object py_str, const char* err) except +
    object distance_impl(object, object) nogil except +
    double ratio_impl(object, object) nogil except +

    object hamming_impl(object, object) nogil except +

def distance(string1, string2):
    """
    Compute absolute Levenshtein distance of two strings.
    
    distance(string1, string2)
    
    Examples (it's hard to spell Levenshtein correctly):
    
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
    validate_string(string1, "string1 must be a String")
    validate_string(string2, "string2 must be a String")

    return distance_impl(string1, string2)

def ratio(string1, string2):
    """
    Compute similarity of two strings.
    
    ratio(string1, string2)
    
    The similarity is a number between 0 and 1, it's usually equal or
    somewhat higher than difflib.SequenceMatcher.ratio(), because it's
    based on real minimal edit distance.
    
    Examples:
    
    >>> ratio('Hello world!', 'Holly grail!')  # doctest: +ELLIPSIS
    0.583333...
    >>> ratio('Brian', 'Jesus')
    0.0
    
    Really?  I thought there was some similarity.
    """
    validate_string(string1, "string1 must be a String")
    validate_string(string2, "string2 must be a String")

    return ratio_impl(string1, string2)

def hamming(string1, string2):
    """
    Compute Hamming distance of two strings.
    
    hamming(string1, string2)
    
    The Hamming distance is simply the number of differing characters.
    That means the length of the strings must be the same.
    
    Examples:
    >>> hamming('Hello world!', 'Holly grail!')
    7
    >>> hamming('Brian', 'Jesus')
    5
    """
    validate_string(string1, "string1 must be a String")
    validate_string(string2, "string2 must be a String")

    return hamming_impl(string1, string2)
