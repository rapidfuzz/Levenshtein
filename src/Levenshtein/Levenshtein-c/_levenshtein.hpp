/* @(#) $Id: Levenshtein.h,v 1.22 2005/01/13 20:02:56 yeti Exp $ */
#ifndef LEVENSHTEIN_H
#define LEVENSHTEIN_H

#include <numeric>
#include <memory>
#include <vector>
#include <unordered_set>
#include <rapidfuzz/distance/Indel.hpp>
#include <rapidfuzz/distance/Levenshtein.hpp>

#define LEV_EPSILON 1e-14
#define LEV_INFINITY 1e100

#ifndef size_t
#  include <stdlib.h>
#endif

/* In C, this is just wchar_t and unsigned char, in Python, lev_wchar can
 * be anything.  If you really want to cheat, define wchar_t to any integer
 * type you like before including Levenshtein.h and recompile it. */
#ifndef lev_wchar
#  define lev_wchar wchar_t
#endif
typedef unsigned char lev_byte;

/* Edit opration type
 * DON'T CHANGE! used ad arrays indices and the bits are occasionally used
 * as flags */
typedef enum {
  LEV_EDIT_KEEP = 0,
  LEV_EDIT_REPLACE = 1,
  LEV_EDIT_INSERT = 2,
  LEV_EDIT_DELETE = 3,
  LEV_EDIT_LAST  /* sometimes returned when an error occurs */
} LevEditType;

/* Error codes returned by editop check functions */
typedef enum {
  LEV_EDIT_ERR_OK = 0,
  LEV_EDIT_ERR_TYPE,  /* nonexistent edit type */
  LEV_EDIT_ERR_OUT,  /* edit out of string bounds */
  LEV_EDIT_ERR_ORDER,  /* ops are not ordered */
  LEV_EDIT_ERR_BLOCK,  /* incosistent block boundaries (block ops) */
  LEV_EDIT_ERR_SPAN,  /* sequence is not a full transformation (block ops) */
  LEV_EDIT_ERR_LAST
} LevEditOpError;

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
template <typename CharT>
std::vector<CharT> make_symlist(size_t n, const size_t *lengths,
              const CharT *strings[])
{
  std::vector<CharT> symlist;
  if (std::all_of(lengths, lengths + n, [](size_t x){ return x == 0; }))
  {
    return symlist;
  }

  std::unordered_set<CharT> symmap;
  for (size_t i = 0; i < n; i++) {
    const CharT* stri = strings[i];
    for (size_t j = 0; j < lengths[i]; j++) {
      symmap.insert(stri[j]);
    }
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
template <typename CharT>
CharT* lev_greedy_median(size_t n, const size_t* lengths, const CharT** strings,
                         const double* weights, size_t* medlength)
{
  /* find all symbols */
  std::vector<CharT> symlist = make_symlist(n, lengths, strings);
  if (symlist.empty()) {
    *medlength = 0;
    return (CharT*)calloc(1, sizeof(CharT));
  }

  /* allocate and initialize per-string matrix rows and a common work buffer */
  std::vector<std::unique_ptr<size_t[]>> rows(n);
  size_t maxlen = *std::max_element(lengths, lengths + n);

  for (size_t i = 0; i < n; i++) {
    size_t leni = lengths[i];
    rows[i] = std::make_unique<size_t[]>(leni + 1);
    std::iota(rows[i].get(), rows[i].get() + leni + 1, 0);
  }
  size_t stoplen = 2*maxlen + 1;
  auto row = std::make_unique<size_t[]>(stoplen + 1);

  /* compute final cost of string of length 0 (empty string may be also
   * a valid answer) */
  auto median = std::make_unique<CharT[]>(stoplen);
  /**
   * the total distance of the best median string of
   * given length.  warning!  mediandist[0] is total
   * distance for empty string, while median[] itself
   * is normally zero-based
   */
  auto mediandist = std::make_unique<double[]>(stoplen + 1);
  mediandist[0] = std::inner_product(lengths, lengths + n, weights, 0.0);

  /* build up the approximate median string symbol by symbol
   * XXX: we actually exit on break below, but on the same condition */
  for (size_t len = 1; len <= stoplen; len++) {
    CharT symbol;
    double minminsum = LEV_INFINITY;
    row[0] = len;
    /* iterate over all symbols we may want to add */
    for (size_t j = 0; j < symlist.size(); j++) {
      double totaldist = 0.0;
      double minsum = 0.0;
      symbol = symlist[j];
      /* sum Levenshtein distances from all the strings, with given weights */
      for (size_t i = 0; i < n; i++) {
        const CharT* stri = strings[i];
        size_t *p = rows[i].get();
        size_t leni = lengths[i];
        size_t *end = rows[i].get() + leni;
        size_t min = len;
        size_t x = len; /* == row[0] */
        /* compute how another row of Levenshtein matrix would look for median
         * string with this symbol added */
        while (p < end) {
          size_t D = *(p++) + (symbol != *(stri++));
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
    for (size_t i = 0; i < n; i++) {
      const CharT* stri = strings[i];
      size_t* oldrow = rows[i].get();
      size_t leni = lengths[i];
      /* compute a row of Levenshtein matrix */
      for (size_t k = 1; k <= leni; k++) {
        size_t c1 = oldrow[k] + 1;
        size_t c2 = row[k - 1] + 1;
        size_t c3 = oldrow[k - 1] + (symbol != stri[k - 1]);
        row[k] = c2 > c3 ? c3 : c2;
        if (row[k] > c1)
          row[k] = c1;
      }
      memcpy(oldrow, row.get(), (leni + 1)*sizeof(size_t));
    }
  }

  /* find the string with minimum total distance */
  size_t bestlen = std::distance(mediandist.get(), std::min_element(mediandist.get(), mediandist.get() + stoplen));

  /* return result */
  {
    CharT* result = (CharT*)safe_malloc(bestlen, sizeof(CharT));
    if (!result) {
      return NULL;
    }
    memcpy(result, median.get(), bestlen*sizeof(CharT));
    *medlength = bestlen;
    return result;
  }
}

/*
 * Knowing the distance matrices up to some row, finish the distance
 * computations.
 *
 * string1, len1 are already shortened.
 */
template <typename CharT>
double finish_distance_computations(size_t len1, CharT* string1,
                                    size_t n, const size_t* lengths,
                                    const CharT** strings,
                                    const double *weights, std::vector<std::unique_ptr<size_t[]>>& rows,
                                    std::unique_ptr<size_t[]>& row)
{
  size_t *end;
  size_t i, j;
  size_t offset;  /* row[0]; offset + len1 give together real len of string1 */
  double distsum = 0.0;  /* sum of distances */

  /* catch trivia case */
  if (len1 == 0) {
    for (j = 0; j < n; j++)
      distsum += (double)rows[j][lengths[j]]*weights[j];
    return distsum;
  }

  /* iterate through the strings and sum the distances */
  for (j = 0; j < n; j++) {
    size_t* rowi = rows[j].get();  /* current row */
    size_t leni = lengths[j];  /* current length */
    size_t len = len1;  /* temporary len1 for suffix stripping */
    const CharT* stringi = strings[j];  /* current string */

    /* strip common suffix (prefix CAN'T be stripped) */
    while (len && leni && stringi[leni-1] == string1[len-1]) {
      len--;
      leni--;
    }

    /* catch trivial cases */
    if (len == 0) {
      distsum += (double)rowi[leni]*weights[j];
      continue;
    }
    offset = rowi[0];
    if (leni == 0) {
      distsum += (double)(offset + len)*weights[j];
      continue;
    }

    /* complete the matrix */
    memcpy(row.get(), rowi, (leni + 1)*sizeof(size_t));
    end = row.get() + leni;

    for (i = 1; i <= len; i++) {
      size_t* p = row.get() + 1;
      const CharT char1 = string1[i - 1];
      const CharT* char2p = stringi;
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
 * @medlength: Where the new length of the median should be stored.
 *
 * Tries to make @s a better generalized median string of @strings with
 * small perturbations.
 *
 * It never returns a string with larger SOD than @s; in the worst case, a
 * string identical to @s is returned.
 *
 * Returns: The improved generalized median, as a newly allocated string; its
 *          length is stored in @medlength.
 **/
template <typename CharT>
CharT* lev_median_improve(size_t len, const CharT* s, size_t n, const size_t* lengths,
                             const CharT** strings, const double *weights, size_t *medlength)
{
  /* find all symbols */
  std::vector<CharT> symlist = make_symlist(n, lengths, strings);
  if (symlist.empty()) {
    *medlength = 0;
    return (CharT*)calloc(1, sizeof(CharT));
  }

  /* allocate and initialize per-string matrix rows and a common work buffer */
  std::vector<std::unique_ptr<size_t[]>> rows(n);
  size_t maxlen = *std::max_element(lengths, lengths + n);

  for (size_t i = 0; i < n; i++) {
    size_t leni = lengths[i];
    rows[i] = std::make_unique<size_t[]>(leni + 1);
    std::iota(rows[i].get(), rows[i].get() + leni + 1, 0);
  }
  size_t stoplen = 2*maxlen + 1;
  auto row = std::make_unique<size_t[]>(stoplen + 1);

  /* initialize median to given string */
  auto _median = std::make_unique<CharT[]>(stoplen + 1);
  CharT* median = _median.get() + 1; /* we need -1st element for insertions a pos 0 */
  size_t medlen = len;
  memcpy(median, s, (medlen)*sizeof(CharT));
  double minminsum = finish_distance_computations(medlen, median,
                                           n, lengths, strings,
                                           weights, rows, row);

  /* sequentially try perturbations on all positions */
  for (size_t pos = 0; pos <= medlen; ) {
    CharT orig_symbol, symbol;
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
        sum = finish_distance_computations(medlen - pos, median + pos,
                                           n, lengths, strings,
                                           weights, rows, row);
        if (sum < minminsum) {
          minminsum = sum;
          symbol = symlist[j];
          operation = LEV_EDIT_REPLACE;
        }
      }
      median[pos] = orig_symbol;
    }
    /* FOREACH symbol: try to add it at pos, if some lower the total
     * distance, chooste the best (increase medlength)
     * We simulate insertion by replacing the character at pos-1 */
    orig_symbol = *(median + pos - 1);
    for (size_t j = 0; j < symlist.size(); j++) {
      *(median + pos - 1) = symlist[j];
      sum = finish_distance_computations(medlen - pos + 1, median + pos - 1,
                                          n, lengths, strings,
                                         weights, rows, row);
      if (sum < minminsum) {
        minminsum = sum;
        symbol = symlist[j];
        operation = LEV_EDIT_INSERT;
      }
    }
    *(median + pos - 1) = orig_symbol;
    /* IF pos < medlength: try to delete the symbol at pos, if it lowers
     * the total distance remember it (decrease medlength) */
    if (pos < medlen) {
      sum = finish_distance_computations(medlen - pos - 1, median + pos + 1,
                                         n, lengths, strings,
                                         weights, rows, row);
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
              (medlen - pos)*sizeof(CharT));
      median[pos] = symbol;
      medlen++;
      break;

    case LEV_EDIT_DELETE:
      memmove(median+pos, median + pos+1,
              (medlen - pos-1)*sizeof(CharT));
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
      for (size_t i = 0; i < n; i++) {
        const CharT* stri = strings[i];
        size_t* oldrow = rows[i].get();
        size_t leni = lengths[i];
        /* compute a row of Levenshtein matrix */
        for (size_t k = 1; k <= leni; k++) {
          size_t c1 = oldrow[k] + 1;
          size_t c2 = row[k - 1] + 1;
          size_t c3 = oldrow[k - 1] + (symbol != stri[k - 1]);
          row[k] = c2 > c3 ? c3 : c2;
          if (row[k] > c1)
            row[k] = c1;
        }
        memcpy(oldrow, row.get(), (leni + 1)*sizeof(size_t));
      }
      pos++;
    }
  }

  /* return result */
  CharT *result = (CharT*)safe_malloc(medlen, sizeof(CharT));
  if (!result) {
    return NULL;
  }
  *medlength = medlen;
  memcpy(result, median, medlen*sizeof(CharT));
  return result;
}

lev_byte*
lev_quick_median(size_t n,
                 const size_t *lengths,
                 const lev_byte *strings[],
                 const double *weights,
                 size_t *medlength);

lev_wchar*
lev_u_quick_median(size_t n,
                   const size_t *lengths,
                   const lev_wchar *strings[],
                   const double *weights,
                   size_t *medlength);

/**
 * lev_set_median_index:
 * @n: The size of @lengths, @strings, and @weights.
 * @lengths: The lengths of @strings.
 * @strings: An array of strings, that may contain NUL characters.
 * @weights: The string weights (they behave exactly as multiplicities, though
 *           any positive value is allowed, not just integers).
 *
 * Finds the median string of a string set @strings.
 *
 * Returns: An index in @strings pointing to the set median, -1 in case of
 *          failure.
 **/
template <typename CharT>
size_t lev_set_median_index(size_t n, const size_t* lengths, const CharT** strings, const double* weights)
{
  size_t minidx = 0;
  double mindist = LEV_INFINITY;
  std::vector<long int> distances(n*(n - 1)/2, 0xff);

  for (size_t i = 0; i < n; i++) {
    size_t j = 0;
    double dist = 0.0;
    const CharT* stri = strings[i];
    size_t leni = lengths[i];
    /* below diagonal */
    while (j < i && dist < mindist) {
      size_t dindex = (i - 1)*(i - 2)/2 + j;
      long int d;
      if (distances[dindex] >= 0)
        d = distances[dindex];
      else {
        d = rapidfuzz::levenshtein_distance(strings[j], strings[j] + lengths[j], stri, stri + leni);
      }
      dist += weights[j] * (double)d;
      j++;
    }
    j++;  /* no need to compare item with itself */
    /* above diagonal */
    while (j < n && dist < mindist) {
      size_t dindex = (j - 1)*(j - 2)/2 + i;
      distances[dindex] = rapidfuzz::levenshtein_distance(strings[j], strings[j] + lengths[j], stri, stri + leni);
      dist += weights[j] * (double)distances[dindex];
      j++;
    }

    if (dist < mindist) {
      mindist = dist;
      minidx = i;
    }
  }

  return minidx;
}

/**
 * lev_set_median:
 * @n: The size of @lengths, @strings, and @weights.
 * @lengths: The lengths of @strings.
 * @strings: An array of strings, that may contain NUL characters.
 * @weights: The string weights (they behave exactly as multiplicities, though
 *           any positive value is allowed, not just integers).
 * @medlength: Where the length of the median string should be stored.
 *
 * Finds the median string of a string set @strings.
 *
 * Returns: The set median as a newly allocate string, its length is stored
 *          in @medlength.  %NULL in the case of failure.
 **/
template <typename CharT>
CharT* lev_set_median(size_t n, const size_t *lengths, const CharT** strings,
                         const double* weights, size_t* medlength)
{
  size_t minidx = lev_set_median_index(n, lengths, strings, weights);
  CharT* result;

  if (minidx == (size_t)-1)
    return NULL;

  *medlength = lengths[minidx];
  if (!lengths[minidx])
    return (CharT*)calloc(1, sizeof(CharT));

  result = (CharT*)safe_malloc(lengths[minidx], sizeof(CharT));
  if (!result)
    return NULL;
  return (CharT*)memcpy(result, strings[minidx], lengths[minidx]*sizeof(CharT));
}

/**
 * lev_u_edit_seq_distance:
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
template <typename CharT>
double lev_edit_seq_distance(size_t n1, const size_t* lengths1, const CharT** strings1,
                             size_t n2, const size_t* lengths2, const CharT** strings2)
{
  /* strip common prefix */
  while (n1 > 0 && n2 > 0
         && *lengths1 == *lengths2
         && memcmp(*strings1, *strings2,
                   *lengths1*sizeof(CharT)) == 0) {
    n1--;
    n2--;
    strings1++;
    strings2++;
    lengths1++;
    lengths2++;
  }

  /* strip common suffix */
  while (n1 > 0 && n2 > 0
         && lengths1[n1-1] == lengths2[n2-1]
         && memcmp(strings1[n1-1], strings2[n2-1],
                   lengths1[n1-1]*sizeof(CharT)) == 0) {
    n1--;
    n2--;
  }

  /* catch trivial cases */
  if (n1 == 0)
    return (double)n2;
  if (n2 == 0)
    return (double)n1;

  /* make the inner cycle (i.e. strings2) the longer one */
  if (n1 > n2) {
    std::swap(n1, n2);
    std::swap(lengths1, lengths2);
    std::swap(strings1, strings2);
  }
  n1++;
  n2++;

  /* initalize first row */
  auto row = std::make_unique<double[]>(n2);
  double* last = row.get() + n2 - 1;
  double* end = row.get() + n2;
  std::iota(row.get(), end, 0.0);

  /* go through the matrix and compute the costs.  yes, this is an extremely
   * obfuscated version, but also extremely memory-conservative and relatively
   * fast.  */
  for (size_t i = 1; i < n1; i++) {
    double* p = row.get() + 1;
    const CharT* str1 = strings1[i - 1];
    const size_t len1 = lengths1[i - 1];
    const CharT** str2p = strings2;
    const size_t *len2p = lengths2;
    double D = (double)i - 1.0;
    double x = (double)i;
    while (p != end) {
      size_t l = len1 + *len2p;
      double q;
      if (l == 0)
        q = D;
      else {
        size_t d = rapidfuzz::indel_distance(str1, str1 + len1, *str2p, *str2p + *len2p);
        len2p++;
        str2p++;
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
template <typename CharT>
double lev_set_distance(size_t n1, const size_t* lengths1, const CharT** strings1,
                        size_t n2, const size_t *lengths2, const CharT** strings2)
{
  /* catch trivial cases */
  if (n1 == 0)
    return (double)n2;
  if (n2 == 0)
    return (double)n1;

  /* make the number of columns (n1) smaller than the number of rows */
  if (n1 > n2) {
    std::swap(n1, n2);
    std::swap(lengths1, lengths2);
    std::swap(strings1, strings2);
  }

  /* compute distances from each to each */
  if (SIZE_MAX / n1 <= n2)
  {
    throw std::bad_alloc();
  }
  auto dists = std::make_unique<double[]>(n1 * n2);
  double* r = dists.get();

  for (size_t i = 0; i < n2; i++) {
    size_t len2 = lengths2[i];
    const CharT *str2 = strings2[i];
    const size_t *len1p = lengths1;
    const CharT **str1p = strings1;
    for (size_t j = 0; j < n1; j++) {
      size_t l = len2 + *len1p;
      if (l == 0)
        *(r++) = 0.0;
      else {
        size_t d = rapidfuzz::indel_distance(str2, str2 + len2, *str1p, *str1p + *len1p);
        len1p++;
        str1p++;
        *(r++) = (double)d / (double)l;
      }
    }
  }

  /* find the optimal mapping between the two sets */
  auto map = munkers_blackman(n1, n2, dists.get());

  /* sum the set distance */
  double sum = (double)(n2 - n1);
  for (size_t j = 0; j < n1; j++) {
    size_t l;
    size_t i = map[j];
    l = lengths1[j] + lengths2[i];
    if (l > 0) {
      size_t d = rapidfuzz::indel_distance(strings1[j], strings1[j] + lengths1[j],
                                           strings2[i], strings2[i] + lengths2[i]);
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

lev_wchar*
lev_u_opcodes_apply(size_t len1,
                    const lev_wchar* string1,
                    size_t len2,
                    const lev_wchar* string2,
                    size_t nb,
                    const LevOpCode *bops,
                    size_t *len);

LevEditOp*
lev_editops_subtract(size_t n,
                     const LevEditOp *ops,
                     size_t ns,
                     const LevEditOp *sub,
                     size_t *nrem);

#endif /* not LEVENSHTEIN_H */
