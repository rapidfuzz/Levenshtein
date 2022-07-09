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
#include <Python.h>
#include "_levenshtein.hpp"

#include <algorithm>

/* Me thinks the second argument of PyArg_UnpackTuple() should be const.
 * Anyway I habitually pass a constant string.
 * A cast to placate the compiler. */
#define PYARGCFIX(x) (char*)(x)

/* python interface and wrappers */
/* declarations and docstrings {{{ */
static PyObject* quickmedian_py(PyObject* self, PyObject* args);

#define Levenshtein_DESC ""

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

#define METHODS_ITEM(x) { #x, x##_py, METH_VARARGS, x##_DESC }
static PyMethodDef methods[] = {
  METHODS_ITEM(quickmedian),
  { NULL, NULL, 0, NULL },
};

typedef std::basic_string<lev_byte> (*MedianFuncString)(size_t n,
                                      const size_t* lengths,
                                      const lev_byte** strings,
                                      const double* weights);
typedef std::basic_string<Py_UNICODE> (*MedianFuncUnicode)(size_t n,
                                         const size_t* lengths,
                                         const Py_UNICODE** strings,
                                         const double* weights);
struct MedianFuncs {
  MedianFuncString s;
  MedianFuncUnicode u;
};

static int
extract_stringlist(PyObject* list,
                   const char* name,
                   size_t n,
                   size_t** sizelist,
                   void* strlist);

static double*
extract_weightlist(PyObject* wlist,
                   const char* name,
                   size_t n);

static PyObject*
median_common(PyObject* args,
              const char* name,
              MedianFuncs foo);

/* }}} */

/****************************************************************************
 *
 * Python interface and subroutines
 *
 ****************************************************************************/
/* {{{ */

static PyObject* quickmedian_py(PyObject*, PyObject* args)
{
  MedianFuncs engines = { lev_quick_median, lev_u_quick_median };
  return median_common(args, "quickmedian", engines);
}

static PyObject* median_common(PyObject* args, const char* name, MedianFuncs foo)
{
  size_t n;
  void* strings = NULL;
  size_t* sizes = NULL;
  PyObject* strlist = NULL;
  PyObject* wlist = NULL;
  PyObject* strseq = NULL;
  double* weights;
  int stringtype;
  PyObject* result = NULL;

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
    try {
      std::basic_string<lev_byte> medstr = foo.s(n, sizes, (const lev_byte**)strings, weights);
      result = PyBytes_FromStringAndSize((const char*)medstr.data(), (Py_ssize_t)medstr.size());
    } catch (...)
    {
      result = PyErr_NoMemory();
    }
  }
  else if (stringtype == 1) {
    try {
      std::basic_string<Py_UNICODE> medstr = foo.u(n, sizes, (const Py_UNICODE**)strings, weights);
      result = PyUnicode_FromUnicode(medstr.data(), (Py_ssize_t)medstr.size());
    } catch (...)
    {
      result = PyErr_NoMemory();
    }
  }
  else
    PyErr_Format(PyExc_SystemError, "%s internal error", name);

  free(strings);
  free(weights);
  free(sizes);
  return result;
}

static double* extract_weightlist(PyObject* wlist, const char* name, size_t n)
{
  double* weights = NULL;
  PyObject* seq;

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
    for (size_t i = 0; i < n; i++) {
      PyObject* item = PySequence_Fast_GET_ITEM(wlist, i);
      PyObject* number = PyNumber_Float(item);

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
    
    std::fill(weights, weights + n, 1.0);
  }

  return weights;
}

/* extract a list of strings or unicode strings, returns
 * 0 -- strings
 * 1 -- unicode strings
 * <0 -- failure
 */
static int extract_stringlist(PyObject* list, const char* name,
                              size_t n, size_t** sizelist, void* strlist)
{
  /* check first object type.  when it's a string then all others must be
   * strings too; when it's a unicode string then all others must be unicode
   * strings too. */
  PyObject* first = PySequence_Fast_GET_ITEM(list, 0);
  /* i don't exactly understand why the problem doesn't exhibit itself earlier
   * but a queer error message is better than a segfault :o) */
  if (first == (PyObject*)-1) {
    PyErr_Format(PyExc_TypeError,
                 "%s undecomposable Sequence???", name);
    return -1;
  }

  if (PyObject_TypeCheck(first, &PyBytes_Type)) {
    lev_byte** strings = (lev_byte**)safe_malloc(n, sizeof(lev_byte*));
    if (!strings) {
      PyErr_Format(PyExc_MemoryError,
                   "%s cannot allocate memory", name);
      return -1;
    }
    size_t* sizes = (size_t*)safe_malloc(n, sizeof(size_t));
    if (!sizes) {
      free(strings);
      PyErr_Format(PyExc_MemoryError,
                   "%s cannot allocate memory", name);
      return -1;
    }

    strings[0] = (lev_byte*)PyBytes_AS_STRING(first);
    sizes[0] = (size_t)PyBytes_GET_SIZE(first);
    for (size_t i = 1; i < n; i++) {
      PyObject* item = PySequence_Fast_GET_ITEM(list, i);

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
    Py_UNICODE** strings = (Py_UNICODE**)safe_malloc(n, sizeof(Py_UNICODE*));
    if (!strings) {
      PyErr_NoMemory();
      return -1;
    }
    size_t* sizes = (size_t*)safe_malloc(n, sizeof(size_t));
    if (!sizes) {
      free(strings);
      PyErr_NoMemory();
      return -1;
    }

    strings[0] = PyUnicode_AS_UNICODE(first);
    sizes[0] = (size_t)PyUnicode_GET_LENGTH(first);
    for (size_t i = 1; i < n; i++) {
      PyObject* item = PySequence_Fast_GET_ITEM(list, i);

      if (!PyObject_TypeCheck(item, &PyUnicode_Type)) {
        free(strings);
        free(sizes);
        PyErr_Format(PyExc_TypeError,
                     "%s item #%i is not a Unicode", name, i);
        return -1;
      }
      strings[i] = PyUnicode_AS_UNICODE(item);
      sizes[i] = (size_t)PyUnicode_GET_LENGTH(item);
    }

    *(Py_UNICODE***)strlist = strings;
    *sizelist = sizes;
    return 1;
  }

  PyErr_Format(PyExc_TypeError,
               "%s expected list of Strings or Unicodes", name);
  return -1;
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
  return PyModule_Create(&moduledef);
}
/* }}} */
