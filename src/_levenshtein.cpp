/*
 * Levenshtein.c
 * @(#) $Id: Levenshtein.c,v 1.41 2005/01/13 20:05:36 yeti Exp $
 * Python extension computing Levenshtein distances, string similarities,
 * median strings and other goodies.
 *
 * Copyright (C) 2002-2003 David Necas (Yeti) <yeti@physics.muni.cz>.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
 **/
#define lev_wchar Py_UNICODE
#include <Python.h>
#include "_levenshtein.hpp"

#define LEV_UNUSED(x) ((void)x)

/* Me thinks the second argument of PyArg_UnpackTuple() should be const.
 * Anyway I habitually pass a constant string.
 * A cast to placate the compiler. */
#define PYARGCFIX(x) (char*)(x)

/* python interface and wrappers */
/* declarations and docstrings {{{ */
static PyObject* median_py(PyObject *self, PyObject *args);
static PyObject* median_improve_py(PyObject *self, PyObject *args);
static PyObject* quickmedian_py(PyObject *self, PyObject *args);
static PyObject* setmedian_py(PyObject *self, PyObject *args);
static PyObject* seqratio_py(PyObject *self, PyObject *args);
static PyObject* setratio_py(PyObject *self, PyObject *args);
static PyObject* apply_edit_py(PyObject *self, PyObject *args);
static PyObject* matching_blocks_py(PyObject *self, PyObject *args);
static PyObject* subtract_edit_py(PyObject *self, PyObject *args);

#define Levenshtein_DESC \
  "A C extension module for fast computation of:\n" \
  "- Levenshtein (edit) distance and edit sequence manipulation\n" \
  "- string similarity\n" \
  "- approximate median strings, and generally string averaging\n" \
  "- string sequence and set similarity\n" \
  "\n" \
  "Levenshtein has a some overlap with difflib (SequenceMatcher).  It\n" \
  "supports only strings, not arbitrary sequence types, but on the\n" \
  "other hand it's much faster.\n" \
  "\n" \
  "It supports both normal and Unicode strings, but can't mix them, all\n" \
  "arguments to a function (method) have to be of the same type (or its\n" \
  "subclasses).\n"

#define median_DESC \
  "Find an approximate generalized median string using greedy algorithm.\n" \
  "\n" \
  "median(string_sequence[, weight_sequence])\n" \
  "\n" \
  "You can optionally pass a weight for each string as the second\n" \
  "argument.  The weights are interpreted as item multiplicities,\n" \
  "although any non-negative real numbers are accepted.  Use them to\n" \
  "improve computation speed when strings often appear multiple times\n" \
  "in the sequence.\n" \
  "\n" \
  "Examples:\n" \
  "\n" \
  ">>> median(['SpSm', 'mpamm', 'Spam', 'Spa', 'Sua', 'hSam'])\n" \
  "'Spam'\n" \
  ">>> fixme = ['Levnhtein', 'Leveshein', 'Leenshten',\n" \
  "...          'Leveshtei', 'Lenshtein', 'Lvenstein',\n" \
  "...          'Levenhtin', 'evenshtei']\n" \
  ">>> median(fixme)\n" \
  "'Levenshtein'\n" \
  "\n" \
  "Hm.  Even a computer program can spell Levenshtein better than me.\n"

#define median_improve_DESC \
  "Improve an approximate generalized median string by perturbations.\n" \
  "\n" \
  "median_improve(string, string_sequence[, weight_sequence])\n" \
  "\n" \
  "The first argument is the estimated generalized median string you\n" \
  "want to improve, the others are the same as in median().  It returns\n" \
  "a string with total distance less or equal to that of the given string.\n" \
  "\n" \
  "Note this is much slower than median().  Also note it performs only\n" \
  "one improvement step, calling median_improve() again on the result\n" \
  "may improve it further, though this is unlikely to happen unless the\n" \
  "given string was not very similar to the actual generalized median.\n" \
  "\n" \
  "Examples:\n" \
  "\n" \
  ">>> fixme = ['Levnhtein', 'Leveshein', 'Leenshten',\n" \
  "...          'Leveshtei', 'Lenshtein', 'Lvenstein',\n" \
  "...          'Levenhtin', 'evenshtei']\n" \
  ">>> median_improve('spam', fixme)\n" \
  "'enhtein'\n" \
  ">>> median_improve(median_improve('spam', fixme), fixme)\n" \
  "'Levenshtein'\n" \
  "\n" \
  "It takes some work to change spam to Levenshtein.\n"

#define quickmedian_DESC \
  "Find a very approximate generalized median string, but fast.\n" \
  "\n" \
  "quickmedian(string[, weight_sequence])\n" \
  "\n" \
  "See median() for argument description.\n" \
  "\n" \
  "This method is somewhere between setmedian() and picking a random\n" \
  "string from the set; both speedwise and quality-wise.\n" \
  "\n" \
  "Examples:\n" \
  "\n" \
  ">>> fixme = ['Levnhtein', 'Leveshein', 'Leenshten',\n" \
  "...          'Leveshtei', 'Lenshtein', 'Lvenstein',\n" \
  "...          'Levenhtin', 'evenshtei']\n" \
  ">>> quickmedian(fixme)\n" \
  "'Levnshein'\n"

#define setmedian_DESC \
  "Find set median of a string set (passed as a sequence).\n" \
  "\n" \
  "setmedian(string[, weight_sequence])\n" \
  "\n" \
  "See median() for argument description.\n" \
  "\n" \
  "The returned string is always one of the strings in the sequence.\n" \
  "\n" \
  "Examples:\n" \
  "\n" \
  ">>> setmedian(['ehee', 'cceaes', 'chees', 'chreesc',\n" \
  "...            'chees', 'cheesee', 'cseese', 'chetese'])\n" \
  "'chees'\n" \
  "\n" \
  "You haven't asked me about Limburger, sir.\n"

#define seqratio_DESC \
  "Compute similarity ratio of two sequences of strings.\n" \
  "\n" \
  "seqratio(string_sequence1, string_sequence2)\n" \
  "\n" \
  "This is like ratio(), but for string sequences.  A kind of ratio()\n" \
  "is used to to measure the cost of item change operation for the\n" \
  "strings.\n" \
  "\n" \
  "Examples:\n" \
  "\n" \
  ">>> seqratio(['newspaper', 'litter bin', 'tinny', 'antelope'],\n" \
  "...          ['caribou', 'sausage', 'gorn', 'woody'])\n" \
  "0.21517857142857144\n"

#define setratio_DESC \
  "Compute similarity ratio of two strings sets (passed as sequences).\n" \
  "\n" \
  "setratio(string_sequence1, string_sequence2)\n" \
  "\n" \
  "The best match between any strings in the first set and the second\n" \
  "set (passed as sequences) is attempted.  I.e., the order doesn't\n" \
  "matter here.\n" \
  "\n" \
  "Examples:\n" \
  "\n" \
  ">>> setratio(['newspaper', 'litter bin', 'tinny', 'antelope'],\n" \
  "...          ['caribou', 'sausage', 'gorn', 'woody'])  # doctest: +ELLIPSIS\n" \
  "0.281845...\n" \
  "\n" \
  "No, even reordering doesn't help the tinny words to match the\n" \
  "woody ones.\n"

#define apply_edit_DESC \
  "Apply a sequence of edit operations to a string.\n" \
  "\n" \
  "apply_edit(edit_operations, source_string, destination_string)\n" \
  "\n" \
  "In the case of editops, the sequence can be arbitrary ordered subset\n" \
  "of the edit sequence transforming source_string to destination_string.\n" \
  "\n" \
  "Examples:\n" \
  "\n" \
  ">>> e = editops('man', 'scotsman')\n" \
  ">>> apply_edit(e, 'man', 'scotsman')\n" \
  "'scotsman'\n" \
  ">>> apply_edit(e[:3], 'man', 'scotsman')\n" \
  "'scoman'\n" \
  "\n" \
  "The other form of edit operations, opcodes, is not very suitable for\n" \
  "such a tricks, because it has to always span over complete strings,\n" \
  "subsets can be created by carefully replacing blocks with 'equal'\n" \
  "blocks, or by enlarging 'equal' block at the expense of other blocks\n" \
  "and adjusting the other blocks accordingly.\n" \
  "\n" \
  "Examples:\n" \
  ">>> a, b = 'spam and eggs', 'foo and bar'\n" \
  ">>> e = opcodes(a, b)\n" \
  ">>> apply_edit(inverse(e), b, a)\n" \
  "'spam and eggs'\n" \
  ">>> e[4] = ('equal', 10, 13, 8, 11)\n" \
  ">>> e\n" \
  "[('delete', 0, 1, 0, 0), ('replace', 1, 4, 0, 3), ('equal', 4, 9, 3, 8), ('delete', 9, 10, 8, 8), ('equal', 10, 13, 8, 11)]\n" \
  ">>> apply_edit(e, a, b)\n" \
  "'foo and ggs'\n"

#define matching_blocks_DESC \
  "Find identical blocks in two strings.\n" \
  "\n" \
  "matching_blocks(edit_operations, source_length, destination_length)\n" \
  "\n" \
  "The result is a list of triples with the same meaning as in\n" \
  "SequenceMatcher's get_matching_blocks() output.  It can be used with\n" \
  "both editops and opcodes.  The second and third arguments don't\n" \
  "have to be actually strings, their lengths are enough.\n" \
  "\n" \
  "Examples:\n" \
  "\n" \
  ">>> a, b = 'spam', 'park'\n" \
  ">>> matching_blocks(editops(a, b), a, b)\n" \
  "[(1, 0, 2), (4, 4, 0)]\n" \
  ">>> matching_blocks(editops(a, b), len(a), len(b))\n" \
  "[(1, 0, 2), (4, 4, 0)]\n" \
  "\n" \
  "The last zero-length block is not an error, but it's there for\n" \
  "compatibility with difflib which always emits it.\n" \
  "\n" \
  "One can join the matching blocks to get two identical strings:\n" \
  ">>> a, b = 'dog kennels', 'mattresses'\n" \
  ">>> mb = matching_blocks(editops(a,b), a, b)\n" \
  ">>> ''.join([a[x[0]:x[0]+x[2]] for x in mb])\n" \
  "'ees'\n" \
  ">>> ''.join([b[x[1]:x[1]+x[2]] for x in mb])\n" \
  "'ees'\n"

#define subtract_edit_DESC \
  "Subtract an edit subsequence from a sequence.\n" \
  "\n" \
  "subtract_edit(edit_operations, subsequence)\n" \
  "\n" \
  "The result is equivalent to\n" \
  "editops(apply_edit(subsequence, s1, s2), s2), except that is\n" \
  "constructed directly from the edit operations.  That is, if you apply\n" \
  "it to the result of subsequence application, you get the same final\n" \
  "string as from application complete edit_operations.  It may be not\n" \
  "identical, though (in amibuous cases, like insertion of a character\n" \
  "next to the same character).\n" \
  "\n" \
  "The subtracted subsequence must be an ordered subset of\n" \
  "edit_operations.\n" \
  "\n" \
  "Note this function does not accept difflib-style opcodes as no one in\n" \
  "his right mind wants to create subsequences from them.\n" \
  "\n" \
  "Examples:\n" \
  "\n" \
  ">>> e = editops('man', 'scotsman')\n" \
  ">>> e1 = e[:3]\n" \
  ">>> bastard = apply_edit(e1, 'man', 'scotsman')\n" \
  ">>> bastard\n" \
  "'scoman'\n" \
  ">>> apply_edit(subtract_edit(e, e1), bastard, 'scotsman')\n" \
  "'scotsman'\n" \

#define METHODS_ITEM(x) { #x, x##_py, METH_VARARGS, x##_DESC }
static PyMethodDef methods[] = {
  METHODS_ITEM(median),
  METHODS_ITEM(median_improve),
  METHODS_ITEM(quickmedian),
  METHODS_ITEM(setmedian),
  METHODS_ITEM(seqratio),
  METHODS_ITEM(setratio),
  METHODS_ITEM(apply_edit),
  METHODS_ITEM(matching_blocks),
  METHODS_ITEM(subtract_edit),
  { NULL, NULL, 0, NULL },
};

/* opcode names, these are to be initialized in the init func,
 * indexed by LevEditType values */
typedef struct {
  PyObject* pystring;
  const char *cstring;
  size_t len;
} OpcodeName;

static OpcodeName opcode_names[] = {
  { NULL, "equal", 0 },
  { NULL, "replace", 0 },
  { NULL, "insert", 0 },
  { NULL, "delete", 0 },
};
#define N_OPCODE_NAMES ((sizeof(opcode_names)/sizeof(opcode_names[0])))

typedef lev_byte *(*MedianFuncString)(size_t n,
                                      const size_t *lengths,
                                      const lev_byte *strings[],
                                      const double *weights,
                                      size_t *medlength);
typedef Py_UNICODE *(*MedianFuncUnicode)(size_t n,
                                         const size_t *lengths,
                                         const Py_UNICODE *strings[],
                                         const double *weights,
                                         size_t *medlength);
typedef struct {
  MedianFuncString s;
  MedianFuncUnicode u;
} MedianFuncs;

typedef lev_byte *(*MedianImproveFuncString)(size_t len, const lev_byte *s,
                                             size_t n,
                                             const size_t *lengths,
                                             const lev_byte *strings[],
                                             const double *weights,
                                             size_t *medlength);
typedef Py_UNICODE *(*MedianImproveFuncUnicode)(size_t len, const Py_UNICODE *s,
                                                size_t n,
                                                const size_t *lengths,
                                                const Py_UNICODE *strings[],
                                                const double *weights,
                                                size_t *medlength);
typedef struct {
  MedianImproveFuncString s;
  MedianImproveFuncUnicode u;
} MedianImproveFuncs;

typedef double (*SetSeqFuncString)(size_t n1,
                                   const size_t *lengths1,
                                   const lev_byte *strings1[],
                                   size_t n2,
                                   const size_t *lengths2,
                                   const lev_byte *strings2[]);
typedef double (*SetSeqFuncUnicode)(size_t n1,
                                    const size_t *lengths1,
                                    const Py_UNICODE *strings1[],
                                    size_t n2,
                                    const size_t *lengths2,
                                    const Py_UNICODE *strings2[]);

typedef struct {
  SetSeqFuncString s;
  SetSeqFuncUnicode u;
} SetSeqFuncs;

static int
extract_stringlist(PyObject *list,
                   const char *name,
                   size_t n,
                   size_t **sizelist,
                   void *strlist);

static double*
extract_weightlist(PyObject *wlist,
                   const char *name,
                   size_t n);

static PyObject*
median_common(PyObject *args,
              const char *name,
              MedianFuncs foo);

static PyObject*
median_improve_common(PyObject *args,
                      const char *name,
                      MedianImproveFuncs foo);

static double
setseq_common(PyObject *args,
              const char *name,
              SetSeqFuncs foo,
              size_t *lensum);

/* }}} */

/****************************************************************************
 *
 * Python interface and subroutines
 *
 ****************************************************************************/
/* {{{ */

static PyObject*
median_py(PyObject *self, PyObject *args)
{
  MedianFuncs engines = { lev_greedy_median, lev_u_greedy_median };
  LEV_UNUSED(self);
  return median_common(args, "median", engines);
}

static PyObject*
median_improve_py(PyObject *self, PyObject *args)
{
  MedianImproveFuncs engines = { lev_median_improve, lev_u_median_improve };
  LEV_UNUSED(self);
  return median_improve_common(args, "median_improve", engines);
}

static PyObject*
quickmedian_py(PyObject *self, PyObject *args)
{
  MedianFuncs engines = { lev_quick_median, lev_u_quick_median };
  LEV_UNUSED(self);
  return median_common(args, "quickmedian", engines);
}

static PyObject*
setmedian_py(PyObject *self, PyObject *args)
{
  MedianFuncs engines = { lev_set_median, lev_u_set_median };
  LEV_UNUSED(self);
  return median_common(args, "setmedian", engines);
}

static PyObject*
median_common(PyObject *args, const char *name, MedianFuncs foo)
{
  size_t n, len;
  void *strings = NULL;
  size_t *sizes = NULL;
  PyObject *strlist = NULL;
  PyObject *wlist = NULL;
  PyObject *strseq = NULL;
  double *weights;
  int stringtype;
  PyObject *result = NULL;

  if (!PyArg_UnpackTuple(args, PYARGCFIX(name), 1, 2, &strlist, &wlist))
    return NULL;

  if (!PySequence_Check(strlist)) {
    PyErr_Format(PyExc_TypeError,
                 "%s first argument must be a Sequence", name);
    return NULL;
  }
  strseq = PySequence_Fast(strlist, name);

  n = (size_t)PySequence_Fast_GET_SIZE(strseq);
  if (n == 0) {
    Py_INCREF(Py_None);
    Py_DECREF(strseq);
    return Py_None;
  }

  /* get (optional) weights, use 1 if none specified. */
  weights = extract_weightlist(wlist, name, n);
  if (!weights) {
    Py_DECREF(strseq);
    return NULL;
  }

  stringtype = extract_stringlist(strseq, name, n, &sizes, &strings);
  Py_DECREF(strseq);
  if (stringtype < 0) {
    free(weights);
    return NULL;
  }

  if (stringtype == 0) {
    lev_byte *medstr = foo.s(n, sizes, (const lev_byte**)strings, weights, &len);
    if (!medstr && len)
      result = PyErr_NoMemory();
    else {
      result = PyBytes_FromStringAndSize((const char*)medstr, (Py_ssize_t)len);
      free(medstr);
    }
  }
  else if (stringtype == 1) {
    Py_UNICODE *medstr = foo.u(n, sizes, (const Py_UNICODE**)strings, weights, &len);
    if (!medstr && len)
      result = PyErr_NoMemory();
    else {
      result = PyUnicode_FromUnicode(medstr, (Py_ssize_t)len);
      free(medstr);
    }
  }
  else
    PyErr_Format(PyExc_SystemError, "%s internal error", name);

  free(strings);
  free(weights);
  free(sizes);
  return result;
}

static PyObject*
median_improve_common(PyObject *args, const char *name, MedianImproveFuncs foo)
{
  size_t n, len;
  void *strings = NULL;
  size_t *sizes = NULL;
  PyObject *arg1 = NULL;
  PyObject *strlist = NULL;
  PyObject *wlist = NULL;
  PyObject *strseq = NULL;
  double *weights;
  int stringtype;
  PyObject *result = NULL;

  if (!PyArg_UnpackTuple(args, PYARGCFIX(name), 2, 3, &arg1, &strlist, &wlist))
    return NULL;

  if (PyObject_TypeCheck(arg1, &PyBytes_Type))
    stringtype = 0;
  else if (PyObject_TypeCheck(arg1, &PyUnicode_Type))
    stringtype = 1;
  else {
    PyErr_Format(PyExc_TypeError,
                 "%s first argument must be a String or Unicode", name);
    return NULL;
  }

  if (!PySequence_Check(strlist)) {
    PyErr_Format(PyExc_TypeError,
                 "%s second argument must be a Sequence", name);
    return NULL;
  }
  strseq = PySequence_Fast(strlist, name);

  n = (size_t)PySequence_Fast_GET_SIZE(strseq);
  if (n == 0) {
    Py_INCREF(Py_None);
    Py_DECREF(strseq);
    return Py_None;
  }

  /* get (optional) weights, use 1 if none specified. */
  weights = extract_weightlist(wlist, name, n);
  if (!weights) {
    Py_DECREF(strseq);
    return NULL;
  }

  if (extract_stringlist(strseq, name, n, &sizes, &strings) != stringtype) {
    PyErr_Format(PyExc_TypeError,
                 "%s argument types don't match", name);
    free(weights);
    return NULL;
  }

  Py_DECREF(strseq);
  if (stringtype == 0) {
    lev_byte *s = (lev_byte*)PyBytes_AS_STRING(arg1);
    size_t l = (size_t)PyBytes_GET_SIZE(arg1);
    lev_byte *medstr = foo.s(l, s, n, sizes, (const lev_byte**)strings, weights, &len);
    if (!medstr && len)
      result = PyErr_NoMemory();
    else {
      result = PyBytes_FromStringAndSize((const char*)medstr, (Py_ssize_t)len);
      free(medstr);
    }
  }
  else if (stringtype == 1) {
    Py_UNICODE *s = PyUnicode_AS_UNICODE(arg1);
    size_t l = (size_t)PyUnicode_GET_SIZE(arg1);
    Py_UNICODE *medstr = foo.u(l, s, n, sizes, (const Py_UNICODE**)strings, weights, &len);
    if (!medstr && len)
      result = PyErr_NoMemory();
    else {
      result = PyUnicode_FromUnicode(medstr, (Py_ssize_t)len);
      free(medstr);
    }
  }
  else
    PyErr_Format(PyExc_SystemError, "%s internal error", name);

  free(strings);
  free(weights);
  free(sizes);
  return result;
}

static double*
extract_weightlist(PyObject *wlist, const char *name, size_t n)
{
  size_t i;
  double *weights = NULL;
  PyObject *seq;

  if (wlist) {
    if (!PySequence_Check(wlist)) {
      PyErr_Format(PyExc_TypeError,
                  "%s second argument must be a Sequence", name);
      return NULL;
    }
    seq = PySequence_Fast(wlist, name);
    if ((size_t)PySequence_Fast_GET_SIZE(wlist) != n) {
      PyErr_Format(PyExc_ValueError,
                   "%s got %i strings but %i weights",
                   name, n, PyList_GET_SIZE(wlist));
      Py_DECREF(seq);
      return NULL;
    }
    weights = (double*)safe_malloc(n, sizeof(double));
    if (!weights)
      return (double*)PyErr_NoMemory();
    for (i = 0; i < n; i++) {
      PyObject *item = PySequence_Fast_GET_ITEM(wlist, i);
      PyObject *number = PyNumber_Float(item);

      if (!number) {
        free(weights);
        PyErr_Format(PyExc_TypeError,
                     "%s weight #%i is not a Number", name, i);
        Py_DECREF(seq);
        return NULL;
      }
      weights[i] = PyFloat_AS_DOUBLE(number);
      Py_DECREF(number);
      if (weights[i] < 0) {
        free(weights);
        PyErr_Format(PyExc_ValueError,
                     "%s weight #%i is negative", name, i);
        Py_DECREF(seq);
        return NULL;
      }
    }
    Py_DECREF(seq);
  }
  else {
    weights = (double*)safe_malloc(n, sizeof(double));
    if (!weights)
      return (double*)PyErr_NoMemory();
    for (i = 0; i < n; i++)
      weights[i] = 1.0;
  }

  return weights;
}

/* extract a list of strings or unicode strings, returns
 * 0 -- strings
 * 1 -- unicode strings
 * <0 -- failure
 */
static int
extract_stringlist(PyObject *list, const char *name,
                   size_t n, size_t **sizelist, void *strlist)
{
  size_t i;
  PyObject *first;

  /* check first object type.  when it's a string then all others must be
   * strings too; when it's a unicode string then all others must be unicode
   * strings too. */
  first = PySequence_Fast_GET_ITEM(list, 0);
  /* i don't exactly understand why the problem doesn't exhibit itself earlier
   * but a queer error message is better than a segfault :o) */
  if (first == (PyObject*)-1) {
    PyErr_Format(PyExc_TypeError,
                 "%s undecomposable Sequence???", name);
    return -1;
  }

  if (PyObject_TypeCheck(first, &PyBytes_Type)) {
    lev_byte **strings;
    size_t *sizes;

    strings = (lev_byte**)safe_malloc(n, sizeof(lev_byte*));
    if (!strings) {
      PyErr_Format(PyExc_MemoryError,
                   "%s cannot allocate memory", name);
      return -1;
    }
    sizes = (size_t*)safe_malloc(n, sizeof(size_t));
    if (!sizes) {
      free(strings);
      PyErr_Format(PyExc_MemoryError,
                   "%s cannot allocate memory", name);
      return -1;
    }

    strings[0] = (lev_byte*)PyBytes_AS_STRING(first);
    sizes[0] = (size_t)PyBytes_GET_SIZE(first);
    for (i = 1; i < n; i++) {
      PyObject *item = PySequence_Fast_GET_ITEM(list, i);

      if (!PyObject_TypeCheck(item, &PyBytes_Type)) {
        free(strings);
        free(sizes);
        PyErr_Format(PyExc_TypeError,
                     "%s item #%i is not a String", name, i);
        return -1;
      }
      strings[i] = (lev_byte*)PyBytes_AS_STRING(item);
      sizes[i] = (size_t)PyBytes_GET_SIZE(item);
    }

    *(lev_byte***)strlist = strings;
    *sizelist = sizes;
    return 0;
  }
  if (PyObject_TypeCheck(first, &PyUnicode_Type)) {
    Py_UNICODE **strings;
    size_t *sizes;

    strings = (Py_UNICODE**)safe_malloc(n, sizeof(Py_UNICODE*));
    if (!strings) {
      PyErr_NoMemory();
      return -1;
    }
    sizes = (size_t*)safe_malloc(n, sizeof(size_t));
    if (!sizes) {
      free(strings);
      PyErr_NoMemory();
      return -1;
    }

    strings[0] = PyUnicode_AS_UNICODE(first);
    sizes[0] = (size_t)PyUnicode_GET_SIZE(first);
    for (i = 1; i < n; i++) {
      PyObject *item = PySequence_Fast_GET_ITEM(list, i);

      if (!PyObject_TypeCheck(item, &PyUnicode_Type)) {
        free(strings);
        free(sizes);
        PyErr_Format(PyExc_TypeError,
                     "%s item #%i is not a Unicode", name, i);
        return -1;
      }
      strings[i] = PyUnicode_AS_UNICODE(item);
      sizes[i] = (size_t)PyUnicode_GET_SIZE(item);
    }

    *(Py_UNICODE***)strlist = strings;
    *sizelist = sizes;
    return 1;
  }

  PyErr_Format(PyExc_TypeError,
               "%s expected list of Strings or Unicodes", name);
  return -1;
}

static PyObject*
seqratio_py(PyObject *self, PyObject *args)
{
  SetSeqFuncs engines = { lev_edit_seq_distance, lev_u_edit_seq_distance };
  size_t lensum;
  double r = setseq_common(args, "seqratio", engines, &lensum);
  LEV_UNUSED(self);
  if (r < 0)
    return NULL;
  if (lensum == 0)
    return PyFloat_FromDouble(1.0);
  return PyFloat_FromDouble(((double)lensum - r) / (double)lensum);
}

static PyObject*
setratio_py(PyObject *self, PyObject *args)
{
  SetSeqFuncs engines = { lev_set_distance, lev_u_set_distance };
  size_t lensum;
  double r = setseq_common(args, "setratio", engines, &lensum);
  LEV_UNUSED(self);
  if (r < 0)
    return NULL;
  if (lensum == 0)
    return PyFloat_FromDouble(1.0);
  return PyFloat_FromDouble(((double)lensum - r) / (double)lensum);
}

static double
setseq_common(PyObject *args, const char *name, SetSeqFuncs foo,
              size_t *lensum)
{
  size_t n1, n2;
  void *strings1 = NULL;
  void *strings2 = NULL;
  size_t *sizes1 = NULL;
  size_t *sizes2 = NULL;
  PyObject *strlist1;
  PyObject *strlist2;
  PyObject *strseq1;
  PyObject *strseq2;
  int stringtype1, stringtype2;
  double r = -1.0;

  if (!PyArg_UnpackTuple(args, PYARGCFIX(name), 2, 2, &strlist1, &strlist2))
    return r;

  if (!PySequence_Check(strlist1)) {
    PyErr_Format(PyExc_TypeError,
                 "%s first argument must be a Sequence", name);
    return r;
  }
  if (!PySequence_Check(strlist2)) {
    PyErr_Format(PyExc_TypeError,
                 "%s second argument must be a Sequence", name);
    return r;
  }

  strseq1 = PySequence_Fast(strlist1, name);
  strseq2 = PySequence_Fast(strlist2, name);

  n1 = (size_t)PySequence_Fast_GET_SIZE(strseq1);
  n2 = (size_t)PySequence_Fast_GET_SIZE(strseq2);
  *lensum = n1 + n2;
  if (n1 == 0) {
    Py_DECREF(strseq1);
    Py_DECREF(strseq2);
    return (double)n2;
  }
  if (n2 == 0) {
    Py_DECREF(strseq1);
    Py_DECREF(strseq2);
    return (double)n1;
  }

  stringtype1 = extract_stringlist(strseq1, name, n1, &sizes1, &strings1);
  Py_DECREF(strseq1);
  if (stringtype1 < 0) {
    Py_DECREF(strseq2);
    return r;
  }
  stringtype2 = extract_stringlist(strseq2, name, n2, &sizes2, &strings2);
  Py_DECREF(strseq2);
  if (stringtype2 < 0) {
    free(sizes1);
    free(strings1);
    return r;
  }

  if (stringtype1 != stringtype2) {
    PyErr_Format(PyExc_TypeError,
                  "%s both sequences must consist of items of the same type",
                  name);
  }
  else if (stringtype1 == 0) {
    r = foo.s(n1, sizes1, (const lev_byte**)strings1, n2, sizes2, (const lev_byte**)strings2);
    if (r < 0.0)
      PyErr_NoMemory();
  }
  else if (stringtype1 == 1) {
    r = foo.u(n1, sizes1, (const Py_UNICODE**)strings1, n2, sizes2, (const Py_UNICODE**)strings2);
    if (r < 0.0)
      PyErr_NoMemory();
  }
  else
    PyErr_Format(PyExc_SystemError, "%s internal error", name);

  free(strings1);
  free(strings2);
  free(sizes1);
  free(sizes2);
  return r;
}

static LevEditType
string_to_edittype(PyObject *string)
{
  size_t i;

  for (i = 0; i < N_OPCODE_NAMES; i++) {
    if (string == opcode_names[i].pystring)
      return (LevEditType)i;
  }

  /* With Python >= 2.2, we shouldn't get here, except when the strings are
   * not Strings but subtypes. */
  /* For Python 3, the string is an unicode object; use CompareWithAsciiString */
  if (!PyUnicode_Check(string)) {
    return LEV_EDIT_LAST;
  }

  for (i = 0; i < N_OPCODE_NAMES; i++) {
    if (PyUnicode_CompareWithASCIIString(string, opcode_names[i].cstring) == 0) {
      return (LevEditType)i;
    }
  }

  return LEV_EDIT_LAST;
}

static LevEditOp*
extract_editops(PyObject *list)
{
  LevEditOp *ops;
  size_t i;
  LevEditType type;
  size_t n = (size_t)PyList_GET_SIZE(list);

  ops = (LevEditOp*)safe_malloc(n, sizeof(LevEditOp));
  if (!ops)
    return (LevEditOp*)PyErr_NoMemory();
  for (i = 0; i < n; i++) {
    PyObject *item;
    PyObject *tuple = PyList_GET_ITEM(list, i);

    if (!PyTuple_Check(tuple) || PyTuple_GET_SIZE(tuple) != 3) {
      free(ops);
      return NULL;
    }
    item = PyTuple_GET_ITEM(tuple, 0);
    if ((type = string_to_edittype(item)) == LEV_EDIT_LAST) {
      free(ops);
      return NULL;
    }
    ops[i].type = type;
    item = PyTuple_GET_ITEM(tuple, 1);
    if (!PyLong_Check(item)) {
      free(ops);
      return NULL;
    }
    ops[i].spos = (size_t)PyLong_AsLong(item);
    item = PyTuple_GET_ITEM(tuple, 2);
    if (!PyLong_Check(item)) {
      free(ops);
      return NULL;
    }
    ops[i].dpos = (size_t)PyLong_AsLong(item);
  }
  return ops;
}

static LevOpCode*
extract_opcodes(PyObject *list)
{
  LevOpCode *bops;
  size_t i;
  LevEditType type;
  size_t nb = (size_t)PyList_GET_SIZE(list);

  bops = (LevOpCode*)safe_malloc(nb, sizeof(LevOpCode));
  if (!bops)
    return (LevOpCode*)PyErr_NoMemory();
  for (i = 0; i < nb; i++) {
    PyObject *item;
    PyObject *tuple = PyList_GET_ITEM(list, i);

    if (!PyTuple_Check(tuple) || PyTuple_GET_SIZE(tuple) != 5) {
      free(bops);
      return NULL;
    }

    item = PyTuple_GET_ITEM(tuple, 0);
    if ((type = string_to_edittype(item)) == LEV_EDIT_LAST) {
      free(bops);
      return NULL;
    }
    bops[i].type = type;

    item = PyTuple_GET_ITEM(tuple, 1);
    if (!PyLong_Check(item)) {
      free(bops);
      return NULL;
    }
    bops[i].sbeg = (size_t)PyLong_AsLong(item);

    item = PyTuple_GET_ITEM(tuple, 2);
    if (!PyLong_Check(item)) {
      free(bops);
      return NULL;
    }
    bops[i].send = (size_t)PyLong_AsLong(item);

    item = PyTuple_GET_ITEM(tuple, 3);
    if (!PyLong_Check(item)) {
      free(bops);
      return NULL;
    }
    bops[i].dbeg = (size_t)PyLong_AsLong(item);

    item = PyTuple_GET_ITEM(tuple, 4);
    if (!PyLong_Check(item)) {
      free(bops);
      return NULL;
    }
    bops[i].dend = (size_t)PyLong_AsLong(item);
  }
  return bops;
}

static PyObject*
editops_to_tuple_list(size_t n, LevEditOp *ops)
{
  PyObject *list;
  size_t i;

  list = PyList_New((Py_ssize_t)n);
  for (i = 0; i < n; i++, ops++) {
    PyObject *tuple = PyTuple_New(3);
    PyObject *is = opcode_names[ops->type].pystring;
    Py_INCREF(is);
    PyTuple_SET_ITEM(tuple, 0, is);
    PyTuple_SET_ITEM(tuple, 1, PyLong_FromLong((long)ops->spos));
    PyTuple_SET_ITEM(tuple, 2, PyLong_FromLong((long)ops->dpos));
    PyList_SET_ITEM(list, i, tuple);
  }

  return list;
}

static PyObject*
matching_blocks_to_tuple_list(size_t len1, size_t len2,
                              size_t nmb, LevMatchingBlock *mblocks)
{
  PyObject *list, *tuple;
  size_t i;

  list = PyList_New((Py_ssize_t)(nmb + 1));
  for (i = 0; i < nmb; i++, mblocks++) {
    tuple = PyTuple_New(3);
    PyTuple_SET_ITEM(tuple, 0, PyLong_FromLong((long)mblocks->spos));
    PyTuple_SET_ITEM(tuple, 1, PyLong_FromLong((long)mblocks->dpos));
    PyTuple_SET_ITEM(tuple, 2, PyLong_FromLong((long)mblocks->len));
    PyList_SET_ITEM(list, i, tuple);
  }
  tuple = PyTuple_New(3);
  PyTuple_SET_ITEM(tuple, 0, PyLong_FromLong((long)len1));
  PyTuple_SET_ITEM(tuple, 1, PyLong_FromLong((long)len2));
  PyTuple_SET_ITEM(tuple, 2, PyLong_FromLong((long)0));
  PyList_SET_ITEM(list, nmb, tuple);

  return list;
}

static size_t
get_length_of_anything(PyObject *object)
{
  if (PyLong_Check(object)) {
    long int len = PyLong_AsLong(object);
    if (len < 0)
      len = -1;
    return (size_t)len;
  }
  if (PySequence_Check(object))
    return (size_t)PySequence_Length(object);
  return (size_t)-1;
}

static PyObject*
opcodes_to_tuple_list(size_t nb, LevOpCode *bops)
{
  PyObject *list;
  size_t i;

  list = PyList_New((Py_ssize_t)nb);
  for (i = 0; i < nb; i++, bops++) {
    PyObject *tuple = PyTuple_New(5);
    PyObject *is = opcode_names[bops->type].pystring;
    Py_INCREF(is);
    PyTuple_SET_ITEM(tuple, 0, is);
    PyTuple_SET_ITEM(tuple, 1, PyLong_FromLong((long)bops->sbeg));
    PyTuple_SET_ITEM(tuple, 2, PyLong_FromLong((long)bops->send));
    PyTuple_SET_ITEM(tuple, 3, PyLong_FromLong((long)bops->dbeg));
    PyTuple_SET_ITEM(tuple, 4, PyLong_FromLong((long)bops->dend));
    PyList_SET_ITEM(list, i, tuple);
  }

  return list;
}

static PyObject*
apply_edit_py(PyObject *self, PyObject *args)
{
  PyObject *list, *result, *arg1, *arg2;
  size_t n, len, len1, len2;
  LevEditOp *ops;
  LevOpCode *bops;
  LEV_UNUSED(self);

  if (!PyArg_UnpackTuple(args, PYARGCFIX("apply_edit"), 3, 3,
                         &list, &arg1, &arg2))
    return NULL;

  if (!PyList_Check(list)) {
    PyErr_SetString(PyExc_TypeError,
                 "apply_edit first argument must be a List of edit operations");
    return NULL;
  }
  n = (size_t)PyList_GET_SIZE(list);

  if (PyObject_TypeCheck(arg1, &PyBytes_Type)
      && PyObject_TypeCheck(arg2, &PyBytes_Type)) {
    lev_byte *string1, *string2, *s;

    if (!n) {
      Py_INCREF(arg1);
      return arg1;
    }
    len1 = (size_t)PyBytes_GET_SIZE(arg1);
    len2 = (size_t)PyBytes_GET_SIZE(arg2);
    string1 = (lev_byte*)PyBytes_AS_STRING(arg1);
    string2 = (lev_byte*)PyBytes_AS_STRING(arg2);

    if ((ops = extract_editops(list)) != NULL) {
      if (lev_editops_check_errors(len1, len2, n, ops)) {
        PyErr_SetString(PyExc_ValueError,
                     "apply_edit edit operations are invalid or inapplicable");
        free(ops);
        return NULL;
      }
      s = lev_editops_apply(len1, string1, len2, string2,
                            n, ops, &len);
      free(ops);
      if (!s && len)
        return PyErr_NoMemory();
      result = PyBytes_FromStringAndSize((const char*)s, (Py_ssize_t)len);
      free(s);
      return result;
    }
    if ((bops = extract_opcodes(list)) != NULL) {
      if (lev_opcodes_check_errors(len1, len2, n, bops)) {
        PyErr_SetString(PyExc_ValueError,
                     "apply_edit edit operations are invalid or inapplicable");
        free(bops);
        return NULL;
      }
      s = lev_opcodes_apply(len1, string1, len2, string2,
                            n, bops, &len);
      free(bops);
      if (!s && len)
        return PyErr_NoMemory();
      result = PyBytes_FromStringAndSize((const char*)s, (Py_ssize_t)len);
      free(s);
      return result;
    }

    if (!PyErr_Occurred())
      PyErr_SetString(PyExc_TypeError,
                  "apply_edit first argument must be "
                  "a list of edit operations");
    return NULL;
  }
  if (PyObject_TypeCheck(arg1, &PyUnicode_Type)
      && PyObject_TypeCheck(arg2, &PyUnicode_Type)) {
    Py_UNICODE *string1, *string2, *s;

    if (!n) {
      Py_INCREF(arg1);
      return arg1;
    }
    len1 = (size_t)PyUnicode_GET_SIZE(arg1);
    len2 = (size_t)PyUnicode_GET_SIZE(arg2);
    string1 = PyUnicode_AS_UNICODE(arg1);
    string2 = PyUnicode_AS_UNICODE(arg2);

    if ((ops = extract_editops(list)) != NULL) {
      if (lev_editops_check_errors(len1, len2, n, ops)) {
        PyErr_SetString(PyExc_ValueError,
                     "apply_edit edit operations are invalid or inapplicable");
        free(ops);
        return NULL;
      }
      s = lev_u_editops_apply(len1, string1, len2, string2,
                              n, ops, &len);
      free(ops);
      if (!s && len)
        return PyErr_NoMemory();
      result = PyUnicode_FromUnicode(s, (Py_ssize_t)len);
      free(s);
      return result;
    }
    if ((bops = extract_opcodes(list)) != NULL) {
      if (lev_opcodes_check_errors(len1, len2, n, bops)) {
        PyErr_SetString(PyExc_ValueError,
                     "apply_edit edit operations are invalid or inapplicable");
        free(bops);
        return NULL;
      }
      s = lev_u_opcodes_apply(len1, string1, len2, string2,
                              n, bops, &len);
      free(bops);
      if (!s && len)
        return PyErr_NoMemory();
      result = PyUnicode_FromUnicode(s, (Py_ssize_t)len);
      free(s);
      return result;
    }

    if (!PyErr_Occurred())
      PyErr_SetString(PyExc_TypeError,
                   "apply_edit first argument must be "
                   "a list of edit operations");
    return NULL;
  }

  PyErr_SetString(PyExc_TypeError,
               "apply_edit expected two Strings or two Unicodes");
  return NULL;
}

static PyObject*
matching_blocks_py(PyObject *self, PyObject *args)
{
  PyObject *list, *arg1, *arg2, *result;
  size_t n, nmb, len1, len2;
  LevEditOp *ops;
  LevOpCode *bops;
  LevMatchingBlock *mblocks;
  LEV_UNUSED(self);

  if (!PyArg_UnpackTuple(args, PYARGCFIX("matching_blocks"), 3, 3,
                         &list, &arg1, &arg2)
      || !PyList_Check(list))
    return NULL;

  if (!PyList_Check(list)) {
    PyErr_SetString(PyExc_TypeError,
                 "matching_blocks first argument must be "
                 "a List of edit operations");
    return NULL;
  }
  n = (size_t)PyList_GET_SIZE(list);
  len1 = get_length_of_anything(arg1);
  len2 = get_length_of_anything(arg2);
  if (len1 == (size_t)-1 || len2 == (size_t)-1) {
    PyErr_SetString(PyExc_ValueError,
                 "matching_blocks second and third argument "
                 "must specify sizes");
    return NULL;
  }

  if ((ops = extract_editops(list)) != NULL) {
    if (lev_editops_check_errors(len1, len2, n, ops)) {
      PyErr_SetString(PyExc_ValueError,
                   "matching_blocks edit operations are invalid or inapplicable");
      free(ops);
      return NULL;
    }
    mblocks = lev_editops_matching_blocks(len1, len2, n, ops, &nmb);
    free(ops);
    if (!mblocks && nmb)
      return PyErr_NoMemory();
    result = matching_blocks_to_tuple_list(len1, len2, nmb, mblocks);
    free(mblocks);
    return result;
  }
  if ((bops = extract_opcodes(list)) != NULL) {
    if (lev_opcodes_check_errors(len1, len2, n, bops)) {
      PyErr_SetString(PyExc_ValueError,
                   "matching_blocks edit operations are invalid or inapplicable");
      free(bops);
      return NULL;
    }
    mblocks = lev_opcodes_matching_blocks(len1, len2, n, bops, &nmb);
    free(bops);
    if (!mblocks && nmb)
      return PyErr_NoMemory();
    result = matching_blocks_to_tuple_list(len1, len2, nmb, mblocks);
    free(mblocks);
    return result;
  }

  if (!PyErr_Occurred())
    PyErr_SetString(PyExc_TypeError,
                "matching_blocks expected a list of edit operations");
  return NULL;
}

static PyObject*
subtract_edit_py(PyObject *self, PyObject *args)
{
  PyObject *list, *sub, *result;
  size_t n, ns, nr;
  LevEditOp *ops, *osub, *orem;
  LEV_UNUSED(self);

  if (!PyArg_UnpackTuple(args, PYARGCFIX("subtract_edit"), 2, 2, &list, &sub)
      || !PyList_Check(list))
    return NULL;

  ns = (size_t)PyList_GET_SIZE(sub);
  if (!ns) {
    Py_INCREF(list);
    return list;
  }
  n = (size_t)PyList_GET_SIZE(list);
  if (!n) {
    PyErr_SetString(PyExc_ValueError,
                 "subtract_edit subsequence is not a subsequence "
                 "or is invalid");
    return NULL;
  }

  if ((ops = extract_editops(list)) != NULL) {
      if ((osub = extract_editops(sub)) != NULL) {
          orem = lev_editops_subtract(n, ops, ns, osub, &nr);
          free(ops);
          free(osub);

          if (!orem && nr == (size_t)-1) {
              PyErr_SetString(PyExc_ValueError,
                           "subtract_edit subsequence is not a subsequence "
                           "or is invalid");
              return NULL;
          }
          result = editops_to_tuple_list(nr, orem);
          free(orem);

          return result;
      }
      free(ops);
  }

  if (!PyErr_Occurred())
    PyErr_SetString(PyExc_TypeError,
                "subtract_edit expected two lists of edit operations");
  return NULL;
}

static PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "_levenshtein",
  Levenshtein_DESC,
  -1,
  methods
};

PyMODINIT_FUNC PyInit__levenshtein(void)
{
  PyObject *module;
  size_t i;

  module = PyModule_Create(&moduledef);
  /* create intern strings for edit operation names */
  if (opcode_names[0].pystring)
    abort();
  for (i = 0; i < N_OPCODE_NAMES; i++) {
    opcode_names[i].pystring = PyUnicode_InternFromString(opcode_names[i].cstring);
    opcode_names[i].len = strlen(opcode_names[i].cstring);
  }
  return module;
}
/* }}} */
