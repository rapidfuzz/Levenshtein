# distutils: language=c++
# cython: language_level=3
# cython: binding=True

from libc.stdlib cimport free
from libc.string cimport strlen
from cpython.list cimport PyList_New, PyList_SET_ITEM
from cpython.object cimport PyObject
from cpython.ref cimport Py_INCREF
from cpython.unicode cimport (
    PyUnicode_CompareWithASCIIString, PyUnicode_AS_UNICODE
)
from cpython.bytes cimport PyBytes_AS_STRING
from cpython.sequence cimport PySequence_Check, PySequence_Length
from libc.stddef cimport wchar_t

cdef extern from "_levenshtein.hpp":
    ctypedef unsigned char lev_byte

    void* safe_malloc(size_t nmemb, size_t size)

    ctypedef enum LevEditType:
        LEV_EDIT_KEEP = 0,
        LEV_EDIT_REPLACE = 1,
        LEV_EDIT_INSERT = 2,
        LEV_EDIT_DELETE = 3,
        LEV_EDIT_LAST

    ctypedef struct LevEditOp:
        LevEditType type
        size_t spos
        size_t dpos

    ctypedef struct LevOpCode:
        LevEditType type
        size_t sbeg
        size_t send
        size_t dbeg
        size_t dend

    void lev_editops_invert(size_t n, LevEditOp *ops)
    void lev_opcodes_invert(size_t nb, LevOpCode *bops)

    int lev_editops_check_errors(size_t len1, size_t len2, size_t n, const LevEditOp *ops)
    int lev_opcodes_check_errors(size_t len1, size_t len2, size_t nb, const LevOpCode *bops)

    LevEditOp* lev_editops_find(size_t len1, const lev_byte *string1, size_t len2, const lev_byte *string2, size_t *n)
    LevEditOp* lev_u_editops_find(size_t len1, const wchar_t *string1, size_t len2, const wchar_t *string2, size_t *n)

    LevEditOp* lev_opcodes_to_editops(size_t nb, const LevOpCode *bops, size_t *n, int keepkeep)
    LevOpCode* lev_editops_to_opcodes(size_t n, const LevEditOp *ops, size_t *nb, size_t len1, size_t len2)

ctypedef struct OpcodeName:
    PyObject* pystring
    const char *cstring
    size_t len

cdef OpcodeName opcode_names[4]
opcode_names[0] = OpcodeName(<PyObject*>"equal",   "equal",   strlen("equal"))
opcode_names[1] = OpcodeName(<PyObject*>"replace", "replace", strlen("replace"))
opcode_names[2] = OpcodeName(<PyObject*>"insert",  "insert",  strlen("insert"))
opcode_names[3] = OpcodeName(<PyObject*>"delete",  "delete",  strlen("delete"))
cdef size_t N_OPCODE_NAMES = 4

cdef size_t get_length_of_anything(o):
    cdef Py_ssize_t length
    if isinstance(o, int):
        length = <Py_ssize_t>o
        if length < 0:
            return <size_t>-1
        return <size_t>length

    if PySequence_Check(o):
        return <size_t>PySequence_Length(o)
    
    return <size_t>-1

cdef LevEditType string_to_edittype(string):
    for i in range(N_OPCODE_NAMES):
        if <PyObject*>string == opcode_names[i].pystring:
           return <LevEditType>i

    if not isinstance(string, str):
        return LEV_EDIT_LAST

    for i in range(N_OPCODE_NAMES):
        if not PyUnicode_CompareWithASCIIString(string, <char*>opcode_names[i].cstring):
            return <LevEditType>i

    return LEV_EDIT_LAST


cdef LevEditOp* extract_editops(list editops) except *:
    cdef size_t n = <size_t>len(editops)
    cdef LevEditOp* ops = <LevEditOp*>safe_malloc(n, sizeof(LevEditOp))

    if not ops:
        raise MemoryError

    for i in range(n):
        editop = editops[i]

        if not isinstance(editop, tuple) or len(<tuple>editop) != 3:
            free(ops)
            return NULL
        
        _type, spos, dpos = <tuple>editop
        if not isinstance(spos, int) or not isinstance(dpos, int):
            free(ops)
            return NULL

        ops[i].spos = <size_t>spos
        ops[i].dpos = <size_t>dpos
        ops[i].type = string_to_edittype(_type)
        if ops[i].type == LEV_EDIT_LAST:
            free(ops)
            return NULL
    
    return ops


cdef LevOpCode* extract_opcodes(list opcodes) except *:
    cdef size_t nb = <size_t>len(opcodes)
    cdef LevOpCode *bops = <LevOpCode*>safe_malloc(nb, sizeof(LevOpCode))

    if not bops:
        raise MemoryError

    for i in range(nb):
        opcode = opcodes[i]

        if not isinstance(opcode, tuple) or len(<tuple>opcode) !=5:
            free(bops)
            return NULL
        
        _type, sbeg, send, dbeg, dend = <tuple>opcode
        if (not isinstance(sbeg, int) or not isinstance(send, int) or
               not isinstance(dbeg, int) or not isinstance(dend, int)):
            free(bops)
            return NULL

        bops[i].sbeg = <size_t>sbeg
        bops[i].send = <size_t>send
        bops[i].dbeg = <size_t>dbeg
        bops[i].dend = <size_t>dend
        bops[i].type = string_to_edittype(_type)
        if bops[i].type == LEV_EDIT_LAST:
            free(bops)
            return NULL
    
    return bops


cdef editops_to_tuple_list(size_t n, LevEditOp *ops):
    cdef list tuple_list = PyList_New(<Py_ssize_t>n)

    for i in range(n):
        result_item = (
            <object>opcode_names[<size_t>ops[i].type].pystring,
            ops[i].spos, ops[i].dpos)
        Py_INCREF(result_item)
        PyList_SET_ITEM(tuple_list, <Py_ssize_t>i, result_item)

    return tuple_list


cdef opcodes_to_tuple_list(size_t nb, LevOpCode *bops):
    cdef list tuple_list = PyList_New(<Py_ssize_t>nb)

    for i in range(nb):
        result_item = (
            <object>opcode_names[<size_t>bops[i].type].pystring,
            bops[i].sbeg, bops[i].send,
            bops[i].dbeg, bops[i].dend)
        Py_INCREF(result_item)
        PyList_SET_ITEM(tuple_list, <Py_ssize_t>i, result_item)

    return tuple_list


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
    cdef size_t n
    cdef LevEditOp *ops
    cdef LevOpCode *bops

    if not isinstance(edit_operations, list):
        raise TypeError("inverse expected a list of edit operations")

    n = <size_t>len(<list>edit_operations)
    if not n:
        return edit_operations

    ops = extract_editops(edit_operations)
    if ops:
        lev_editops_invert(n, ops)
        result = editops_to_tuple_list(n, ops)
        free(ops)
        return result

    bops = extract_opcodes(edit_operations)
    if bops:
       lev_opcodes_invert(n, bops)
       result = opcodes_to_tuple_list(n, bops)
       free(bops)
       return result


    raise TypeError("inverse expected a list of edit operations")


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
    cdef size_t n, len1, len2
    cdef LevEditOp *ops
    cdef LevOpCode *bops

    # convert: we were called (bops, s1, s2)
    if len(args) == 3:
        arg1, arg2, arg3 = args

        if not isinstance(arg1, list):
            raise ValueError("editops first argument must be a List of edit operations")
        
        n = <size_t>len(<list>arg1)
        if not n:
            return arg1

        len1 = get_length_of_anything(arg2)
        len2 = get_length_of_anything(arg3)
        if len1 == <size_t>-1 or len2 == <size_t>-1:
            raise ValueError("editops second and third argument must specify sizes")

        bops = extract_opcodes(arg1)
        if bops:
            if lev_opcodes_check_errors(len1, len2, n, bops):
                free(bops)
                raise ValueError("editops edit operation list is invalid")

            ops = lev_opcodes_to_editops(n, bops, &n, 0)
            free(bops)

            if not ops and n:
                raise MemoryError

            oplist = editops_to_tuple_list(n, ops)
            free(ops)
            return oplist

        ops = extract_editops(arg1)
        if ops:
            if lev_editops_check_errors(len1, len2, n, ops):
                free(ops)
                raise ValueError("editops edit operation list is invalid")
            
            free(ops)
            return arg1
        
        raise TypeError("editops first argument must be a List of edit operations")

    # find editops: we were called (s1, s2)
    arg1, arg2 = args
    if isinstance(arg1, bytes) and isinstance(arg2, bytes):
        len1 = len(<bytes>arg1)
        len2 = len(<bytes>arg2)

        ops = lev_editops_find(
            len1, <lev_byte*>PyBytes_AS_STRING(arg1),
            len2, <lev_byte*>PyBytes_AS_STRING(arg2),
            &n)

    elif isinstance(arg1, str) and isinstance(arg2, str):
        len1 = len(<str>arg1)
        len2 = len(<str>arg2)

        ops = lev_u_editops_find(
            len1, <wchar_t*>PyUnicode_AS_UNICODE(arg1),
            len2, <wchar_t*>PyUnicode_AS_UNICODE(arg2),
            &n)

    else:
        raise TypeError("editops expected two Strings or two Unicodes")

    if not ops and n:
        raise MemoryError
  
    oplist = editops_to_tuple_list(n, ops)
    free(ops)
    return oplist


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
    cdef size_t n, nb, len1, len2
    cdef LevEditOp *ops
    cdef LevOpCode *bops

    # convert: we were called (ops, s1, s2)
    if len(args) == 3:
        arg1, arg2, arg3 = args

        if not isinstance(arg1, list):
            raise ValueError("opcodes first argument must be a List of edit operations")
        
        n = <size_t>len(<list>arg1)
        len1 = get_length_of_anything(arg2)
        len2 = get_length_of_anything(arg3)
        if len1 == <size_t>-1 or len2 == <size_t>-1:
            raise ValueError("opcodes second and third argument must specify sizes")

        ops = extract_editops(arg1)
        if ops:
            if lev_editops_check_errors(len1, len2, n, ops):
                free(ops)
                raise ValueError("opcodes edit operation list is invalid")

            bops = lev_editops_to_opcodes(n, ops, &n, len1, len2)
            free(ops)

            if not bops and n:
                raise MemoryError

            oplist = opcodes_to_tuple_list(n, bops)
            free(bops)
            return oplist

        bops = extract_opcodes(arg1)
        if bops:
            if lev_opcodes_check_errors(len1, len2, n, bops):
                free(bops)
                raise ValueError("opcodes edit operation list is invalid")
            
            free(bops)
            return arg1
        
        raise TypeError("opcodes first argument must be a List of edit operations")

    # find editops: we were called (s1, s2)
    arg1, arg2 = args
    if isinstance(arg1, bytes) and isinstance(arg2, bytes):
        len1 = len(<bytes>arg1)
        len2 = len(<bytes>arg2)

        ops = lev_editops_find(
            len1, <lev_byte*>PyBytes_AS_STRING(arg1),
            len2, <lev_byte*>PyBytes_AS_STRING(arg2),
            &n)

    elif isinstance(arg1, str) and isinstance(arg2, str):
        len1 = len(<str>arg1)
        len2 = len(<str>arg2)

        ops = lev_u_editops_find(
            len1, <wchar_t*>PyUnicode_AS_UNICODE(arg1),
            len2, <wchar_t*>PyUnicode_AS_UNICODE(arg2),
            &n)

    else:
        raise TypeError("opcodes expected two Strings or two Unicodes")

    if not ops and n:
        raise MemoryError
  
    bops = lev_editops_to_opcodes(n, ops, &nb, len1, len2)
    free(ops)
    if not bops and nb:
        raise MemoryError

    oplist = opcodes_to_tuple_list(nb, bops)
    free(bops)
    return oplist
