#pragma once
#include "Python.h"
#define RAPIDFUZZ_PYTHON
#include <rapidfuzz/string_metric.hpp>
#include <exception>

#define PYTHON_VERSION(major, minor, micro) ((major << 24) | (minor << 16) | (micro << 8))

namespace string_metric = rapidfuzz::string_metric;

class PythonTypeError: public std::bad_typeid {
public:

    PythonTypeError(char const* error)
      : m_error(error) {}

    virtual char const* what() const noexcept {
        return m_error;
    }
private:
    char const* m_error;
};

#if PY_VERSION_HEX >= PYTHON_VERSION(3, 0, 0)
#define LIST_OF_CASES(...)   \
    X_ENUM(LEVENSHTEIN_UINT8,   uint8_t) \
    X_ENUM(LEVENSHTEIN_UINT16,  uint16_t) \
    X_ENUM(LEVENSHTEIN_UINT32,  uint32_t)
#else
#define LIST_OF_CASES(...)   \
    X_ENUM(LEVENSHTEIN_UINT8,   uint8_t) \
    X_ENUM(LEVENSHTEIN_UNICODE, Py_UNICODE)
#endif

enum LevenshteinType {
# define X_ENUM(kind, type) kind,
    LIST_OF_CASES()
# undef X_ENUM
};


struct proc_string {
    LevenshteinType kind;
    void* data;
    size_t length;
};

#if PY_VERSION_HEX >= PYTHON_VERSION(3, 0, 0)
static inline void validate_string(PyObject* py_str, const char* err)
{
    if (PyBytes_Check(py_str)) {
        return;
    }

    if (PyUnicode_Check(py_str)) {
        // PEP 623 deprecates legacy strings and therefor
        // deprecates e.g. PyUnicode_READY in Python 3.10
#if PY_VERSION_HEX < PYTHON_VERSION(3, 10, 0)
        if (PyUnicode_READY(py_str)) {
          // cython will use the exception set by PyUnicode_READY
          throw std::runtime_error("");
        }
#endif
        return;
    }

    throw PythonTypeError(err);
}

static inline LevenshteinType UnicodeToLevenshteinKind(unsigned int kind)
{
    switch(kind){
    case PyUnicode_1BYTE_KIND: return LEVENSHTEIN_UINT8;
    case PyUnicode_2BYTE_KIND: return LEVENSHTEIN_UINT16;
    default: return LEVENSHTEIN_UINT32;
    }
}

static inline proc_string convert_string(PyObject* py_str)
{
    if (PyBytes_Check(py_str)) {
        return {
            LEVENSHTEIN_UINT8,
            PyBytes_AS_STRING(py_str),
            static_cast<std::size_t>(PyBytes_GET_SIZE(py_str))
        };
    } else {

        return {
            UnicodeToLevenshteinKind(PyUnicode_KIND(py_str)),
            PyUnicode_DATA(py_str),
            static_cast<std::size_t>(PyUnicode_GET_LENGTH(py_str))
        };
    }
}
#else
static inline void validate_string(PyObject* py_str, const char* err)
{
    if (PyString_Check(py_str) || PyUnicode_Check(py_str)) {
        return;
    }

    throw PythonTypeError(err);
}

static inline proc_string convert_string(PyObject* py_str)
{
    if (PyString_Check(py_str)) {
        return {
            LEVENSHTEIN_UINT8,
            PyString_AS_STRING(py_str),
            static_cast<std::size_t>(PyString_GET_SIZE(py_str))
        };
    } else {
        return {
            LEVENSHTEIN_UNICODE,
            PyUnicode_AS_UNICODE(py_str),
            static_cast<std::size_t>(PyUnicode_GET_SIZE(py_str))
        };
    }
}
#endif

template<typename CharT>
rapidfuzz::basic_string_view<CharT> to_string_view(proc_string str)
{
    return rapidfuzz::basic_string_view<CharT>((CharT*)str.data, str.length);
}

/*
 * Levenshtein distance
 */

template<typename CharT>
size_t distance_impl_inner(proc_string s1, proc_string s2)
{
    switch(s2.kind){
# define X_ENUM(KIND, TYPE) case KIND: return string_metric::levenshtein(to_string_view<CharT>(s1), to_string_view<TYPE>(s2));
    LIST_OF_CASES()
# undef X_ENUM
    }
}

PyObject* distance_impl(PyObject* s1, PyObject* s2)
{
    size_t result = 0;
    proc_string c_s1 = convert_string(s1);
    proc_string c_s2 = convert_string(s2);

    switch(c_s1.kind){
# define X_ENUM(KIND, TYPE) case KIND: result = distance_impl_inner<TYPE>(c_s1, c_s2); break;
    LIST_OF_CASES()
# undef X_ENUM
    }

    if (result == (std::size_t)-1) {
        return PyLong_FromLong(-1);
    }
    return PyLong_FromSize_t(result);
}

/*
 *  Levenshtein ratio
 */

template<typename CharT>
inline double ratio_impl_inner(proc_string s1, proc_string s2)
{
    switch(s2.kind){
# define X_ENUM(KIND, TYPE) case KIND: return string_metric::normalized_levenshtein(to_string_view<CharT>(s1), to_string_view<TYPE>(s2), {1, 1, 2}) / 100.0;
    LIST_OF_CASES()
# undef X_ENUM
    }
}

double ratio_impl(PyObject* s1, PyObject* s2)
{
    proc_string c_s1 = convert_string(s1);
    proc_string c_s2 = convert_string(s2);

    switch(c_s1.kind){
# define X_ENUM(KIND, TYPE) case KIND: return ratio_impl_inner<TYPE>(c_s1, c_s2);
    LIST_OF_CASES()
# undef X_ENUM
    }
}

/*
 * Hamming
 */

template<typename CharT>
size_t hamming_impl_inner(proc_string s1, proc_string s2)
{
    switch(s2.kind){
# define X_ENUM(KIND, TYPE) case KIND: return string_metric::hamming(to_string_view<CharT>(s1), to_string_view<TYPE>(s2));
    LIST_OF_CASES()
# undef X_ENUM
    }
}

PyObject* hamming_impl(PyObject* s1, PyObject* s2)
{
    size_t result = 0;
    proc_string c_s1 = convert_string(s1);
    proc_string c_s2 = convert_string(s2);

    switch(c_s1.kind){
# define X_ENUM(KIND, TYPE) case KIND: result = hamming_impl_inner<TYPE>(c_s1, c_s2); break;
    LIST_OF_CASES()
# undef X_ENUM
    }

    if (result == (std::size_t)-1) {
        return PyLong_FromLong(-1);
    }
    return PyLong_FromSize_t(result);
}
