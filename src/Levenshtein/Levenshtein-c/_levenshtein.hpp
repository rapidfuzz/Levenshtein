/* @(#) $Id: Levenshtein.h,v 1.22 2005/01/13 20:02:56 yeti Exp $ */
#ifndef LEVENSHTEIN_H
#define LEVENSHTEIN_H

#include "Python.h"
#include <cstdint>
#include <numeric>
#include <memory>
#include <vector>
#include <string>
#include <set>
#include <rapidfuzz/distance/Indel.hpp>
#include <rapidfuzz/distance/Levenshtein.hpp>

#define LEV_EPSILON 1e-14
#define LEV_INFINITY 1e100

#define PYTHON_VERSION(major, minor, micro) ((major << 24) | (minor << 16) | (micro << 8))

enum RF_StringType {
    RF_UINT8,  /* uint8_t */
    RF_UINT16, /* uint16_t */
    RF_UINT32  /* uint32_t */
};

typedef struct _RF_String {
/* members */
    RF_StringType kind;
    void* data;
    int64_t length;
} RF_String;

#define LIST_OF_CASES()   \
    X_ENUM(RF_UINT8,  uint8_t ) \
    X_ENUM(RF_UINT16, uint16_t) \
    X_ENUM(RF_UINT32, uint32_t)

template <typename Func, typename... Args>
auto visit(const RF_String& str, Func&& f, Args&&... args)
{
    switch(str.kind) {
# define X_ENUM(kind, type) case kind: return f((type*)str.data, (type*)str.data + str.length, std::forward<Args>(args)...);
    LIST_OF_CASES()
# undef X_ENUM
    default:
        throw std::logic_error("Invalid string type");
    }
}

template <typename Func, typename... Args>
auto visitor(const RF_String& str1, const RF_String& str2, Func&& f, Args&&... args)
{
    return visit(str2,
        [&](auto first, auto last) {
            return visit(str1, std::forward<Func>(f), first, last, std::forward<Args>(args)...);
        }
    );
}

static inline int64_t levenshtein_distance_func(const RF_String& s1, const RF_String& s2)
{
    return visitor(s1, s2, [&](auto first1, auto last1, auto first2, auto last2) {
        return rapidfuzz::levenshtein_distance(first1, last1, first2, last2);
    });
}

static inline bool is_valid_string(PyObject* py_str)
{
    bool is_string = false;

    if (PyBytes_Check(py_str)) {
        is_string = true;
    }
    else if (PyUnicode_Check(py_str)) {
        // PEP 623 deprecates legacy strings and therefor
        // deprecates e.g. PyUnicode_READY in Python 3.10
#if PY_VERSION_HEX < PYTHON_VERSION(3, 10, 0)
        if (PyUnicode_READY(py_str)) {
          // cython will use the exception set by PyUnicode_READY
          throw std::runtime_error("");
        }
#endif
        is_string = true;
    }

    return is_string;
}

static inline RF_String convert_string(PyObject* py_str)
{
    if (PyBytes_Check(py_str)) {
        return {
            RF_UINT8,
            PyBytes_AS_STRING(py_str),
            static_cast<int64_t>(PyBytes_GET_SIZE(py_str))
        };
    } else {
        RF_StringType kind;
        switch(PyUnicode_KIND(py_str)) {
        case PyUnicode_1BYTE_KIND:
           kind = RF_UINT8;
           break;
        case PyUnicode_2BYTE_KIND:
           kind = RF_UINT16;
           break;
        default:
           kind = RF_UINT32;
           break;
        }

        return {
            kind,
            PyUnicode_DATA(py_str),
            static_cast<int64_t>(PyUnicode_GET_LENGTH(py_str))
        };
    }
}

/* Edit operation type
 * DON'T CHANGE! used as array indices and the bits are occasionally used
 * as flags */
enum LevEditType {
  LEV_EDIT_KEEP = 0,
  LEV_EDIT_REPLACE = 1,
  LEV_EDIT_INSERT = 2,
  LEV_EDIT_DELETE = 3,
  LEV_EDIT_LAST  /* sometimes returned when an error occurs */
};

/* Error codes returned by editop check functions */
enum LevEditOpError {
  LEV_EDIT_ERR_OK = 0,
  LEV_EDIT_ERR_TYPE,  /* nonexistent edit type */
  LEV_EDIT_ERR_OUT,  /* edit out of string bounds */
  LEV_EDIT_ERR_ORDER,  /* ops are not ordered */
  LEV_EDIT_ERR_BLOCK,  /* incosistent block boundaries (block ops) */
  LEV_EDIT_ERR_SPAN,  /* sequence is not a full transformation (block ops) */
  LEV_EDIT_ERR_LAST
};

/* Edit operation (atomic).
 * This is the `native' atomic edit operation.  It differs from the difflib
 * one's because it represents a change of one character, not a block.  And
 * we usually don't care about LEV_EDIT_KEEP, though the functions can handle
 * them.  The positions are interpreted as at the left edge of a character.
 */
typedef struct {
  LevEditType type;  /* editing operation type */
  size_t spos;  /* source block position */
  size_t dpos;  /* destination position */
} LevEditOp;

/* Edit operation (difflib-compatible).
 * This is not `native', but conversion functions exist.  These fields exactly
 * correspond to the codeops() tuples fields (and this method is also the
 * source of the silly OpCode name).  Sequences must span over complete
 * strings, subsequences are simply edit sequences with more (or larger)
 * LEV_EDIT_KEEP blocks.
 */
typedef struct {
  LevEditType type;  /* editing operation type */
  size_t sbeg, send;  /* source block begin, end */
  size_t dbeg, dend;  /* destination block begin, end */
} LevOpCode;

/* Matching block (difflib-compatible). */
typedef struct {
  size_t spos;
  size_t dpos;
  size_t len;
} LevMatchingBlock;

static void *
safe_malloc(size_t nmemb, size_t size) {
  /* extra-conservative overflow check */
  if (SIZE_MAX / size <= nmemb) {
    return NULL;
  }
  return malloc(nmemb * size);
}

/* compute the sets of symbols each string contains, and the set of symbols
 * in any of them (symset).  meanwhile, count how many different symbols
 * there are (used below for symlist). */
static inline std::vector<uint32_t> make_symlist(const std::vector<RF_String>& strings)
{
  std::vector<uint32_t> symlist;
  if (std::all_of(std::begin(strings), std::end(strings), [](const auto& x){ return x.length == 0; }))
  {
    return symlist;
  }

  std::set<uint32_t> symmap;
  for (const auto& string : strings) {
    visit(string, [&](auto first1, auto last1){
      for (; first1 != last1; ++first1) {
        symmap.insert(*first1);
      }
    });
  }
  /* create dense symbol table, so we can easily iterate over only characters
   * present in the strings */
  symlist.insert(std::end(symlist), std::begin(symmap), std::end(symmap));
  return symlist;
}

/**
 * lev_greedy_median:
 * @n: The size of @lengths, @strings, and @weights.
 * @lengths: The lengths of @strings.
 * @strings: An array of strings, that may contain NUL characters.
 * @weights: The string weights (they behave exactly as multiplicities, though
 *           any positive value is allowed, not just integers).
 * @medlength: Where the length of the median should be stored.
 *
 * Finds a generalized median string of @strings using the greedy algorithm.
 *
 * Note it's considerably more efficient to give a string with weight 2
 * than to store two identical strings in @strings (with weights 1).
 *
 * Returns: The generalized median, as a newly allocated string; its length
 *          is stored in @medlength.
 **/
static inline std::basic_string<uint32_t> lev_greedy_median(const std::vector<RF_String>& strings,
                         const std::vector<double>& weights)
{
  std::basic_string<uint32_t> result_median;

  /* find all symbols */
  std::vector<uint32_t> symlist = make_symlist(strings);
  if (symlist.empty()) {
    return result_median;
  }

  /* allocate and initialize per-string matrix rows and a common work buffer */
  std::vector<std::unique_ptr<size_t[]>> rows(strings.size());
  size_t maxlen = std::max_element(std::begin(strings), std::end(strings), [](const auto& a, const auto& b){
    return a.length < b.length;
  })->length;

  for (size_t i = 0; i < strings.size(); i++) {
    size_t leni = strings[i].length;
    rows[i] = std::make_unique<size_t[]>(leni + 1);
    std::iota(rows[i].get(), rows[i].get() + leni + 1, 0);
  }
  size_t stoplen = 2*maxlen + 1;
  auto row = std::make_unique<size_t[]>(stoplen + 1);

  /* compute final cost of string of length 0 (empty string may be also
   * a valid answer) */
  auto median = std::make_unique<uint32_t[]>(stoplen);
  /**
   * the total distance of the best median string of
   * given length.  warning!  mediandist[0] is total
   * distance for empty string, while median[] itself
   * is normally zero-based
   */
  auto mediandist = std::make_unique<double[]>(stoplen + 1);
  mediandist[0] = 0;
  for (size_t i = 0; i < strings.size(); i++) {
    mediandist[0] += strings[i].length + weights[i];
  }

  /* build up the approximate median string symbol by symbol
   * XXX: we actually exit on break below, but on the same condition */
  for (size_t len = 1; len <= stoplen; len++) {
    uint32_t symbol = 0;
    double minminsum = LEV_INFINITY;
    row[0] = len;
    /* iterate over all symbols we may want to add */
    for (size_t j = 0; j < symlist.size(); j++) {
      double totaldist = 0.0;
      double minsum = 0.0;
      symbol = symlist[j];
      /* sum Levenshtein distances from all the strings, with given weights */
      for (size_t i = 0; i < strings.size(); i++) {
        visit(strings[i], [&](auto first1, auto last1){
          size_t *p = rows[i].get();
          size_t *end = rows[i].get() + std::distance(first1, last1);
          size_t min = len;
          size_t x = len; /* == row[0] */
          /* compute how another row of Levenshtein matrix would look for median
          * string with this symbol added */
          while (p != end) {
            size_t D = *(p++) + (symbol != *(first1++));
            x++;
            if (x > D)
              x = D;
            if (x > *p + 1)
              x = *p + 1;
            if (x < min)
              min = x;
          }
          minsum += (double)min*weights[i];
          totaldist += (double)x*weights[i];
        });
      }
      /* is this symbol better than all the others? */
      if (minsum < minminsum) {
        minminsum = minsum;
        mediandist[len] = totaldist;
        median[len - 1] = symbol;
      }
    }
    /* stop the iteration if we no longer need to recompute the matrix rows
     * or when we are over maxlen and adding more characters doesn't seem
     * useful */
    if (len == stoplen
        || (len > maxlen && mediandist[len] > mediandist[len - 1])) {
      stoplen = len;
      break;
    }
    /* now the best symbol is known, so recompute all matrix rows for this
     * one */
    symbol = median[len - 1];
    for (size_t i = 0; i < strings.size(); i++) {
      visit(strings[i], [&](auto first1, auto last1){
        size_t* oldrow = rows[i].get();
        size_t leni = std::distance(first1, last1);
        /* compute a row of Levenshtein matrix */
        for (size_t k = 1; k <= leni; k++) {
          size_t c1 = oldrow[k] + 1;
          size_t c2 = row[k - 1] + 1;
          size_t c3 = oldrow[k - 1] + (symbol != first1[k - 1]);
          row[k] = c2 > c3 ? c3 : c2;
          if (row[k] > c1)
            row[k] = c1;
        }
        memcpy(oldrow, row.get(), (leni + 1)*sizeof(size_t));
      });
    }
  }

  /* find the string with minimum total distance */
  size_t bestlen = std::distance(mediandist.get(), std::min_element(mediandist.get(), mediandist.get() + stoplen));

  /* return result */
  result_median.insert(std::begin(result_median), median.get(), median.get() + bestlen);
  return result_median;
}


/*
 * Knowing the distance matrices up to some row, finish the distance
 * computations.
 *
 * string1, len1 are already shortened.
 */
static inline double finish_distance_computations(size_t len1, uint32_t* string1,
                                    const std::vector<RF_String>& strings,
                                    const std::vector<double>& weights, std::vector<std::unique_ptr<size_t[]>>& rows,
                                    std::unique_ptr<size_t[]>& row)
{
  double distsum = 0.0;  /* sum of distances */
  /* catch trivial case */
  if (len1 == 0) {
    for (size_t j = 0; j < strings.size(); j++)
      distsum += (double)rows[j][strings[j].length]*weights[j];
    return distsum;
  }

  /* iterate through the strings and sum the distances */
  for (size_t j = 0; j < strings.size(); j++) {
    visit(strings[j], [&](auto first1, auto last1){
      size_t* rowi = rows[j].get();  /* current row */
      size_t leni = (size_t)std::distance(first1, last1);  /* current length */
      size_t len = len1;  /* temporary len1 for suffix stripping */

      /* strip common suffix (prefix CAN'T be stripped) */
      while (len && leni && first1[leni-1] == string1[len-1]) {
        len--;
        leni--;
      }

      /* catch trivial cases */
      if (len == 0) {
        distsum += (double)rowi[leni]*weights[j];
        return;
      }
      /* row[0]; offset + len1 give together real len of string1 */
      size_t offset = rowi[0];
      if (leni == 0) {
        distsum += (double)(offset + len)*weights[j];
        return;
      }

      /* complete the matrix */
      memcpy(row.get(), rowi, (leni + 1)*sizeof(size_t));
      size_t* end = row.get() + leni;

      for (size_t i = 1; i <= len; i++) {
        size_t* p = row.get() + 1;
        const uint32_t char1 = string1[i - 1];
        auto char2p = first1;
        size_t D, x;

        D = x = i + offset;
        while (p <= end) {
          size_t c3 = --D + (char1 != *(char2p++));
          x++;
          if (x > c3)
            x = c3;
          D = *p;
          D++;
          if (x > D)
            x = D;
          *(p++) = x;
        }
      }
      distsum += weights[j]*(double)(*end);
    });
  }

  return distsum;
}

/**
 * lev_median_improve:
 * @len: The length of @s.
 * @s: The approximate generalized median string to be improved.
 * @n: The size of @lengths, @strings, and @weights.
 * @lengths: The lengths of @strings.
 * @strings: An array of strings, that may contain NUL characters.
 * @weights: The string weights (they behave exactly as multiplicities, though
 *           any positive value is allowed, not just integers).
 *
 * Tries to make @s a better generalized median string of @strings with
 * small perturbations.
 *
 * It never returns a string with larger SOD than @s; in the worst case, a
 * string identical to @s is returned.
 *
 * Returns: The improved generalized median
 **/
static inline std::basic_string<uint32_t> lev_median_improve(const RF_String& string,
                          const std::vector<RF_String>& strings, const std::vector<double>& weights)
{
  /* find all symbols */
  std::vector<uint32_t> symlist = make_symlist(strings);
  if (symlist.empty()) {
    return std::basic_string<uint32_t>();
  }

  /* allocate and initialize per-string matrix rows and a common work buffer */
  std::vector<std::unique_ptr<size_t[]>> rows(strings.size());
  size_t maxlen = 0;
  for (size_t i = 0; i < strings.size(); i++) {
    size_t leni = (size_t)strings[i].length;
    if (leni > maxlen)
      maxlen = leni;
    rows[i] = std::make_unique<size_t[]>(leni + 1);
    std::iota(rows[i].get(), rows[i].get() + leni + 1, 0);
  }

  size_t stoplen = 2*maxlen + 1;
  auto row = std::make_unique<size_t[]>(stoplen + 2);

  /* initialize median to given string */
  auto _median = std::make_unique<uint32_t[]>(stoplen + 1);
  uint32_t* median = _median.get() + 1; /* we need -1st element for insertions a pos 0 */
  size_t medlen = (size_t)string.length;

  visit(string, [&](auto first1, auto last1){
    std::copy(first1, last1, median);
  });

  double minminsum = finish_distance_computations(medlen, median, strings, weights, rows, row);

  /* sequentially try perturbations on all positions */
  for (size_t pos = 0; pos <= medlen; ) {
    uint32_t orig_symbol, symbol;
    LevEditType operation;
    double sum;

    symbol = median[pos];
    operation = LEV_EDIT_KEEP;
    /* IF pos < medlength: FOREACH symbol: try to replace the symbol
     * at pos, if some lower the total distance, chooste the best */
    if (pos < medlen) {
      orig_symbol = median[pos];
      for (size_t j = 0; j < symlist.size(); j++) {
        if (symlist[j] == orig_symbol)
          continue;
        median[pos] = symlist[j];
        sum = finish_distance_computations(medlen - pos, median + pos, strings, weights, rows, row);
        if (sum < minminsum) {
          minminsum = sum;
          symbol = symlist[j];
          operation = LEV_EDIT_REPLACE;
        }
      }
      median[pos] = orig_symbol;
    }
    /* FOREACH symbol: try to add it at pos, if some lower the total
     * distance, chooste the best (increase medlen)
     * We simulate insertion by replacing the character at pos-1 */
    orig_symbol = *(median + pos - 1);
    for (size_t j = 0; j < symlist.size(); j++) {
      *(median + pos - 1) = symlist[j];
      sum = finish_distance_computations(medlen - pos + 1, median + pos - 1, strings, weights, rows, row);
      if (sum < minminsum) {
        minminsum = sum;
        symbol = symlist[j];
        operation = LEV_EDIT_INSERT;
      }
    }
    *(median + pos - 1) = orig_symbol;
    /* IF pos < medlen: try to delete the symbol at pos, if it lowers
     * the total distance remember it (decrease medlen) */
    if (pos < medlen) {
      sum = finish_distance_computations(medlen - pos - 1, median + pos + 1, strings, weights, rows, row);
      if (sum < minminsum) {
        minminsum = sum;
        operation = LEV_EDIT_DELETE;
      }
    }
    /* actually perform the operation */
    switch (operation) {
    case LEV_EDIT_REPLACE:
      median[pos] = symbol;
      break;

    case LEV_EDIT_INSERT:
      memmove(median+pos+1, median+pos,
              (medlen - pos)*sizeof(uint32_t));
      median[pos] = symbol;
      medlen++;
      break;

    case LEV_EDIT_DELETE:
      memmove(median+pos, median + pos+1,
              (medlen - pos-1)*sizeof(uint32_t));
      medlen--;
      break;

    default:
      break;
    }
    assert(medlen <= stoplen);
    /* now the result is known, so recompute all matrix rows and move on */
    if (operation != LEV_EDIT_DELETE) {
      symbol = median[pos];
      row[0] = pos + 1;

      for (size_t i = 0; i < strings.size(); i++) {
        visit(strings[i], [&](auto first1, auto last1){
          size_t* oldrow = rows[i].get();
          size_t leni = (size_t)std::distance(first1, last1);
          /* compute a row of Levenshtein matrix */
          for (size_t k = 1; k <= leni; k++) {
            size_t c1 = oldrow[k] + 1;
            size_t c2 = row[k - 1] + 1;
            size_t c3 = oldrow[k - 1] + (symbol != first1[k - 1]);
            row[k] = c2 > c3 ? c3 : c2;
            if (row[k] > c1)
              row[k] = c1;
          }
          memcpy(oldrow, row.get(), (leni + 1)*sizeof(size_t));
        });
      }
      pos++;
    }
  }

  return std::basic_string<uint32_t>(median, medlen);
}

std::basic_string<uint32_t> lev_quick_median(const std::vector<RF_String>& strings, const std::vector<double>& weights);

/**
 * lev_set_median:
 * @n: The size of @lengths, @strings, and @weights.
 * @lengths: The lengths of @strings.
 * @strings: An array of strings, that may contain NUL characters.
 * @weights: The string weights (they behave exactly as multiplicities, though
 *           any positive value is allowed, not just integers).
 *
 * Finds the median string of a string set @strings.
 *
 * Returns: The set median
 **/
static inline std::basic_string<uint32_t> lev_set_median(const std::vector<RF_String>& strings,
                         const std::vector<double>& weights)
{
  size_t minidx = 0;
  double mindist = LEV_INFINITY;
  std::vector<long int> distances(strings.size()*(strings.size() - 1)/2, 0xff);

  for (size_t i = 0; i < strings.size(); i++) {
    size_t j = 0;
    double dist = 0.0;

    /* below diagonal */
    while (j < i && dist < mindist) {
      size_t dindex = (i - 1)*(i - 2)/2 + j;
      long int d;
      if (distances[dindex] >= 0)
        d = distances[dindex];
      else {
        d = levenshtein_distance_func(strings[j], strings[i]);
      }
      dist += weights[j] * (double)d;
      j++;
    }
    j++;  /* no need to compare item with itself */
    /* above diagonal */
    while (j < strings.size() && dist < mindist) {
      size_t dindex = (j - 1)*(j - 2)/2 + i;
      distances[dindex] = levenshtein_distance_func(strings[j], strings[i]);
      dist += weights[j] * (double)distances[dindex];
      j++;
    }

    if (dist < mindist) {
      mindist = dist;
      minidx = i;
    }
  }

  return visit(strings[minidx], [&](auto first1, auto last1) {
      return std::basic_string<uint32_t>(first1, last1);
  });
}


static inline bool is_equal(const RF_String& a, const RF_String& b)
{
  return visitor(a, b, [](auto first1, auto last1, auto first2, auto last2){
    return std::equal(first1, last1, first2, last2);
  });
}

/**
 * lev_edit_seq_distance:
 * @n1: The length of @lengths1 and @strings1.
 * @lengths1: The lengths of strings in @strings1.
 * @strings1: An array of strings that may contain NUL characters.
 * @n2: The length of @lengths2 and @strings2.
 * @lengths2: The lengths of strings in @strings2.
 * @strings2: An array of strings that may contain NUL characters.
 *
 * Finds the distance between string sequences @strings1 and @strings2.
 *
 * In other words, this is a double-Levenshtein algorithm.
 *
 * The cost of string replace operation is based on string similarity: it's
 * zero for identical strings and 2 for completely unsimilar strings.
 *
 * Returns: The distance of the two sequences.
 **/
static inline double lev_edit_seq_distance(const std::vector<RF_String>& strings1, const std::vector<RF_String>& strings2)
{
  if (strings1.size() > strings2.size())
  {
    return lev_edit_seq_distance(strings2, strings1);
  }

  auto strings1_start = std::begin(strings1);
  auto strings1_end = std::end(strings1);
  auto strings2_start = std::begin(strings2);
  auto strings2_end = std::end(strings2);

  /* strip common prefix */
  while (strings1_start != strings1_end
      && strings2_start != strings2_end
      && is_equal(*strings1_start, *strings2_start))
  {
    strings1_start++;
    strings2_start++;
  }

  /* strip common suffix */
  while (strings1_start != strings1_end
      && strings2_start != strings2_end
      && is_equal(*(strings1_end-1), *(strings2_end - 1)))
  {
    strings1_end--;
    strings2_end--;
  }

  /* catch trivial cases */
  if (strings1_start == strings1_end)
    return (double)std::distance(strings2_start, strings2_end);
  if (strings2_start == strings2_end)
    return (double)std::distance(strings1_start, strings1_end);

  /* initalize first row */
  size_t n1 = std::distance(strings1_start, strings1_end);
  size_t n2 = std::distance(strings2_start, strings2_end);
  auto row = std::make_unique<double[]>(n2 + 1);
  double* last = row.get() + n2;
  double* end = row.get() + n2 + 1;
  std::iota(row.get(), end, 0.0);

  /* go through the matrix and compute the costs.  yes, this is an extremely
   * obfuscated version, but also extremely memory-conservative and relatively
   * fast.  */
  for (size_t i = 0; i < n1; i++) {
    double* p = row.get() + 1;
    auto strings2_it = strings2_start;
    double D = (double)i;
    double x = (double)i + 1.0;

    visit(strings1[i], [&](auto first1, auto last1){
      using CharT1 = typename std::iterator_traits<decltype(first1)>::value_type;
      rapidfuzz::CachedIndel<CharT1> scorer(first1, last1);

      while (p != end) {
        size_t l = strings1[i].length + strings2_it->length;
        double q;
        if (l == 0)
          q = D;
        else {
          size_t d = visit(*strings2_it, [&](auto first2, auto last2){
            return scorer.distance(first2, last2);
          });
          strings2_it++;
          q = D + 2.0 / (double)l * (double)d;
        }
        x += 1.0;
        if (x > q)
          x = q;
        D = *p;
        if (x > D + 1.0)
          x = D + 1.0;
        *(p++) = x;
      }
    });
  }

  return *last;
}

std::unique_ptr<size_t[]> munkers_blackman(size_t n1, size_t n2, double *dists);

/**
 * lev_set_distance:
 * @n1: The length of @lengths1 and @strings1.
 * @lengths1: The lengths of strings in @strings1.
 * @strings1: An array of strings that may contain NUL characters.
 * @n2: The length of @lengths2 and @strings2.
 * @lengths2: The lengths of strings in @strings2.
 * @strings2: An array of strings that may contain NUL characters.
 *
 * Finds the distance between string sets @strings1 and @strings2.
 *
 * The difference from lev_edit_seq_distance() is that order doesn't matter.
 * The optimal association of @strings1 and @strings2 is found first and
 * the similarity is computed for that.
 *
 * Uses sequential Munkers-Blackman algorithm.
 *
 * Returns: The distance of the two sets.
 **/
static inline double lev_set_distance(const std::vector<RF_String>& strings1, const std::vector<RF_String>& strings2)
{
  /* catch trivial cases */
  if (strings1.empty())
    return (double)strings2.size();
  if (strings2.empty())
    return (double)strings1.size();

  /* make the number of columns (n1) smaller than the number of rows */
  if (strings1.size() > strings2.size()) {
    return lev_set_distance(strings2, strings1);
  }

  /* compute distances from each to each */
  if (SIZE_MAX / strings1.size() <= strings2.size())
  {
    throw std::bad_alloc();
  }
  auto dists = std::make_unique<double[]>(strings1.size() * strings2.size());
  double* r = dists.get();

  for (size_t i = 0; i < strings2.size(); i++) {
    for (size_t j = 0; j < strings1.size(); j++) {
      size_t l = strings2[i].length + strings1[j].length;
      if (l == 0)
        *(r++) = 0.0;
      else {
        size_t d = visitor(strings1[j], strings2[i], [&](auto first1, auto last1, auto first2, auto last2){
          return rapidfuzz::indel_distance(first1, last1, first2, last2);
        });
        *(r++) = (double)d / (double)l;
      }
    }
  }

  /* find the optimal mapping between the two sets */
  auto map = munkers_blackman(strings1.size(), strings2.size(), dists.get());

  /* sum the set distance */
  double sum = (double)(strings2.size() - strings1.size());
  for (size_t j = 0; j < strings1.size(); j++) {
    size_t l;
    size_t i = map[j];
    l = strings1[j].length + strings2[i].length;
    if (l > 0) {
      size_t d = visitor(strings1[j], strings2[i], [&](auto first1, auto last1, auto first2, auto last2){
        return rapidfuzz::indel_distance(first1, last1, first2, last2);
      });
      sum += 2.0 * (double)d / (double)l;
    }
  }

  return sum;
}

int
lev_editops_check_errors(size_t len1,
                         size_t len2,
                         size_t n,
                         const LevEditOp *ops);

int
lev_opcodes_check_errors(size_t len1,
                         size_t len2,
                         size_t nb,
                         const LevOpCode *bops);

void
lev_editops_invert(size_t n,
                   LevEditOp *ops);

void
lev_opcodes_invert(size_t nb,
                   LevOpCode *bops);

LevMatchingBlock*
lev_editops_matching_blocks(size_t len1,
                            size_t len2,
                            size_t n,
                            const LevEditOp *ops,
                            size_t *nmblocks);

LevMatchingBlock*
lev_opcodes_matching_blocks(size_t len1,
                            size_t len2,
                            size_t nb,
                            const LevOpCode *bops,
                            size_t *nmblocks);

template <typename CharT>
CharT* lev_editops_apply(size_t len1, const CharT* string1,
                         size_t,      const CharT* string2,
                         size_t n, const LevEditOp *ops, size_t *len)
{
    CharT *dst, *dpos;  /* destination string */
    const CharT *spos;  /* source string position */
    size_t i, j;

    /* this looks too complex for such a simple task, but note ops is not
    * a complete edit sequence, we have to be able to apply anything anywhere */
    dpos = dst = (CharT*)safe_malloc((n + len1), sizeof(CharT));
    if (!dst) {
        *len = (size_t)(-1);
        return NULL;
    }
    spos = string1;
    for (i = n; i; i--, ops++) {
        /* XXX: this is fine with gcc internal memcpy, but when memcpy is
        * actually a function, it may be pretty slow */
        j = ops->spos - (size_t)(spos - string1) + (ops->type == LEV_EDIT_KEEP);
        if (j) {
          memcpy(dpos, spos, j*sizeof(CharT));
          spos += j;
          dpos += j;
        }
        switch (ops->type) {
        case LEV_EDIT_DELETE:
            spos++;
            break;
        case LEV_EDIT_REPLACE:
            spos++;
            *(dpos++) = string2[ops->dpos];
            break;
        case LEV_EDIT_INSERT:
            *(dpos++) = string2[ops->dpos];
            break;
        default:
            break;
        }
    }
    j = len1 - (size_t)(spos - string1);
    if (j) {
        memcpy(dpos, spos, j*sizeof(CharT));
        spos += j;
        dpos += j;
    }

    *len = (size_t)(dpos - dst);
    /* possible realloc failure is detected with *len != 0 */
    return (CharT*)realloc(dst, *len*sizeof(CharT));
}


template <typename CharT>
CharT* lev_opcodes_apply(size_t len1, const CharT* string1,
                         size_t len2, const CharT* string2,
                         size_t nb, const LevOpCode *bops, size_t *len)
{
    CharT *dst, *dpos;  /* destination string */
    const CharT *spos;  /* source string position */
    size_t i;

    /* this looks too complex for such a simple task, but note ops is not
    * a complete edit sequence, we have to be able to apply anything anywhere */
    dpos = dst = (CharT*)safe_malloc((len1 + len2), sizeof(CharT));
    if (!dst) {
        *len = (size_t)(-1);
        return NULL;
    }
    spos = string1;
    for (i = nb; i; i--, bops++) {
        switch (bops->type) {
        case LEV_EDIT_INSERT:
        case LEV_EDIT_REPLACE:
            memcpy(dpos, string2 + bops->dbeg,
                  (bops->dend - bops->dbeg)*sizeof(CharT));
            break;
        case LEV_EDIT_KEEP:
            memcpy(dpos, string1 + bops->sbeg,
                (bops->send - bops->sbeg)*sizeof(CharT));
            break;
        default:
            break;
        }
        spos += bops->send - bops->sbeg;
        dpos += bops->dend - bops->dbeg;
    }

    *len = (size_t)(dpos - dst);
    /* possible realloc failure is detected with *len != 0 */
    return (CharT*)realloc(dst, *len*sizeof(CharT));
}

LevEditOp*
lev_editops_subtract(size_t n,
                     const LevEditOp *ops,
                     size_t ns,
                     const LevEditOp *sub,
                     size_t *nrem);

#endif /* not LEVENSHTEIN_H */
