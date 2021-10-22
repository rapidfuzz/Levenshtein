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

#include <string.h>
#include <math.h>
/* for debugging */
#include <stdio.h>
#include <stdint.h>

#include <assert.h>
#include "_levenshtein.hpp"

#include "rapidfuzz/string_metric.hpp"

#include <vector>
#include <memory>

/* local functions */
static size_t*
munkers_blackman(size_t n1,
                 size_t n2,
                 double *dists);

/****************************************************************************
 *
 * Generalized medians, the greedy algorithm, and greedy improvements
 *
 ****************************************************************************/
/* {{{ */

/* compute the sets of symbols each string contains, and the set of symbols
 * in any of them (symset).  meanwhile, count how many different symbols
 * there are (used below for symlist). */
static std::vector<lev_byte>
make_symlist(size_t n, const size_t *lengths,
             const lev_byte *strings[])
{
  std::vector<lev_byte> symlist;
  /* indexed by ALL symbols, contains 1 for symbols
     present in the strings, zero for others */
  std::vector<short int> symset(0x100);
  size_t symlistlen = 0;
  for (size_t i = 0; i < n; i++) {
    const lev_byte *stri = strings[i];
    for (size_t j = 0; j < lengths[i]; j++) {
      int c = stri[j];
      if (!symset[c]) {
        symlistlen++;
        symset[c] = 1;
      }
    }
  }
  if (!symlistlen) {
    return symlist;
  }

  /* create dense symbol table, so we can easily iterate over only characters
   * present in the strings */
  
  symlist.reserve(symlistlen);
  for (size_t j = 0; j < 0x100; j++) {
    if (symset[j])
      symlist.push_back((lev_byte)j);
  }

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
lev_byte*
lev_greedy_median(size_t n, const size_t *lengths,
                  const lev_byte *strings[],
                  const double *weights,
                  size_t *medlength)
{
  /* find all symbols */
  std::vector<lev_byte> symlist = make_symlist(n, lengths, strings);
  if (symlist.empty()) {
    *medlength = 0;
    return NULL;
  }

  /* allocate and initialize per-string matrix rows and a common work buffer */
  /* Levenshtein matrix rows for each string, we need to keep
     only one previous row to construct the current one */
  auto rows = std::make_unique<std::unique_ptr<size_t[]>[]>(n); 
  size_t maxlen = 0;
  for (size_t i = 0; i < n; i++) {
    size_t leni = lengths[i];
    if (leni > maxlen)
      maxlen = leni;
    
    rows[i] = std::make_unique<size_t[]>(leni + 1);
    // todo iota
    for (size_t j = 0; j < leni + 1; j++)
      rows[i][j] = j;
  }

  size_t stoplen = 2 * maxlen + 1;
  auto row = std::make_unique<size_t[]>(stoplen + 1);
  auto median = std::make_unique<lev_byte[]>(stoplen);
  auto mediandist = std::make_unique<double[]>(stoplen + 1);
  mediandist[0] = 0.0;
  for (size_t i = 0; i < n; i++)
    mediandist[0] += (double)lengths[i]*weights[i];

  /* build up the approximate median string symbol by symbol
   * XXX: we actually exit on break below, but on the same condition */
  for (size_t len = 1; len < stoplen + 1; len++) {
    lev_byte symbol;
    double minminsum = LEV_INFINITY;
    row[0] = len;
    /* iterate over all symbols we may want to add */
    for (size_t j = 0; j < symlist.size(); j++) {
      double totaldist = 0.0;
      double minsum = 0.0;
      symbol = symlist[j];
      /* sum Levenshtein distances from all the strings, with given weights */
      for (size_t i = 0; i < n; i++) {
        const lev_byte *stri = strings[i];
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
      const lev_byte *stri = strings[i];
      size_t *oldrow = rows[i].get();
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
  size_t bestlen = 0;
  for (size_t len = 1; len < stoplen + 1; len++) {
    if (mediandist[len] < mediandist[bestlen])
      bestlen = len;
  }


  /* return result */
  {
    lev_byte *result = (lev_byte*)safe_malloc(bestlen, sizeof(lev_byte));
    if (!result) {
      return NULL;
    }
    memcpy(result, median.get(), bestlen*sizeof(lev_byte));
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
static double
finish_distance_computations(size_t len1, lev_byte *string1,
                             size_t n, const size_t *lengths,
                             const lev_byte **strings,
                             const double *weights, const std::unique_ptr<std::unique_ptr<size_t[]>[]>& rows,
                             std::unique_ptr<size_t[]>& row)
{
  double distsum = 0.0;  /* sum of distances */

  /* catch trivia case */
  if (len1 == 0) {
    for (size_t j = 0; j < n; j++)
      distsum += (double)rows[j][lengths[j]]*weights[j];
    return distsum;
  }

  /* iterate through the strings and sum the distances */
  for (size_t j = 0; j < n; j++) {
    const size_t* rowi = rows[j].get();  /* current row */
    size_t leni = lengths[j];  /* current length */
    size_t len = len1;  /* temporary len1 for suffix stripping */
    const lev_byte *stringi = strings[j];  /* current string */

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
    size_t offset = rowi[0];  /* row[0]; offset + len1 give together real len of string1 */
    if (leni == 0) {
      distsum += (double)(offset + len)*weights[j];
      continue;
    }

    /* complete the matrix */
    memcpy(row.get(), rowi, (leni + 1)*sizeof(size_t));
    size_t *end = row.get() + leni;

    for (size_t i = 1; i <= len; i++) {
      size_t *p = row.get() + 1;
      const lev_byte char1 = string1[i - 1];
      const lev_byte *char2p = stringi;
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
lev_byte*
lev_median_improve(size_t len, const lev_byte *s,
                   size_t n, const size_t *lengths,
                   const lev_byte *strings[],
                   const double *weights,
                   size_t *medlength)
{
  size_t i;  /* usually iterates over strings (n) */
  size_t j;  /* usually iterates over characters */
  size_t pos;  /* the position in the approximate median string we are
                  trying to change */
  size_t maxlen;  /* maximum input string length */
  size_t stoplen;  /* maximum tried median string length -- this is slightly
                      higher than maxlen, because the median string may be
                      longer than any of the input strings */
  size_t medlen;  /* the current approximate median string length */
  double minminsum;  /* the current total distance sum */

  /* find all symbols */
  std::vector<lev_byte> symlist = make_symlist(n, lengths, strings);
  if (symlist.empty()) {
    *medlength = 0;
    return NULL;
  }

  /* allocate and initialize per-string matrix rows and a common work buffer */
  /* Levenshtein matrix rows for each string, we need to keep
     only one previous row to construct the current one */
  auto rows = std::make_unique<std::unique_ptr<size_t[]>[]>(n); 
  maxlen = 0;
  for (i = 0; i < n; i++) {
    size_t leni = lengths[i];
    if (leni > maxlen)
      maxlen = leni;
    rows[i] = std::make_unique<size_t[]>(leni + 1);
    // todo iota
    for (j = 0; j <= leni; j++)
      rows[i][j] = j;
  }
  stoplen = 2*maxlen + 1;
  auto row = std::make_unique<size_t[]>(stoplen + 2);
  auto median = std::make_unique<lev_byte[]>(stoplen + 1);
  lev_byte* median_ptr = median.get() + 1; /* we need -1st element for insertions a pos 0 */

  medlen = len;
  memcpy(median_ptr, s, (medlen)*sizeof(lev_byte));
  minminsum = finish_distance_computations(medlen, median_ptr,
                                           n, lengths, strings,
                                           weights, rows, row);

  /* sequentially try perturbations on all positions */
  for (pos = 0; pos <= medlen; ) {
    lev_byte orig_symbol, symbol;
    LevEditType operation;
    double sum;

    symbol = median[pos];
    operation = LEV_EDIT_KEEP;
    /* IF pos < medlength: FOREACH symbol: try to replace the symbol
     * at pos, if some lower the total distance, chooste the best */
    if (pos < medlen) {
      orig_symbol = median[pos];
      for (j = 0; j < symlist.size(); j++) {
        if (symlist[j] == orig_symbol)
          continue;
        median[pos] = symlist[j];
        sum = finish_distance_computations(medlen - pos, median_ptr + pos,
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
    orig_symbol = *(median_ptr + pos - 1);
    for (j = 0; j < symlist.size(); j++) {
      median_ptr[pos - 1] = symlist[j];
      sum = finish_distance_computations(medlen - pos + 1, median_ptr + pos - 1,
                                          n, lengths, strings,
                                         weights, rows, row);
      if (sum < minminsum) {
        minminsum = sum;
        symbol = symlist[j];
        operation = LEV_EDIT_INSERT;
      }
    }
    median_ptr[pos - 1] = orig_symbol;
    /* IF pos < medlength: try to delete the symbol at pos, if it lowers
     * the total distance remember it (decrease medlength) */
    if (pos < medlen) {
      sum = finish_distance_computations(medlen - pos - 1, median_ptr + pos + 1,
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
      memmove(median_ptr+pos+1, median_ptr+pos,
              (medlen - pos)*sizeof(lev_byte));
      median[pos] = symbol;
      medlen++;
      break;

    case LEV_EDIT_DELETE:
      memmove(median_ptr+pos, median_ptr + pos+1,
              (medlen - pos-1)*sizeof(lev_byte));
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
      for (i = 0; i < n; i++) {
        const lev_byte *stri = strings[i];
        size_t *oldrow = rows[i].get();
        size_t leni = lengths[i];
        size_t k;
        /* compute a row of Levenshtein matrix */
        for (k = 1; k <= leni; k++) {
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
  {
    lev_byte *result = (lev_byte*)safe_malloc(medlen, sizeof(lev_byte));
    if (!result) {
      return NULL;
    }
    *medlength = medlen;
    memcpy(result, median_ptr, medlen*sizeof(lev_byte));
    return result;
  }
}

/* used internally in make_usymlist */
typedef struct _HItem HItem;
struct _HItem {
  lev_wchar c;
  HItem *n;
};

/* free usmylist hash
 * this is a separate function because we need it in more than one place */
static void
free_usymlist_hash(HItem *symmap)
{
  for (size_t j = 0; j < 0x100; j++) {
    HItem *p = symmap + j;
    if (p->n == symmap || p->n == NULL)
      continue;
    p = p->n;
    while (p) {
      HItem *q = p;
      p = p->n;
      free(q);
    }
  }
  free(symmap);
}

/* compute the sets of symbols each string contains, and the set of symbols
 * in any of them (symset).  meanwhile, count how many different symbols
 * there are (used below for symlist). */
static lev_wchar*
make_usymlist(size_t n, const size_t *lengths,
              const lev_wchar *strings[], size_t *symlistlen)
{
  size_t j = 0;
  for (size_t i = 0; i < n; i++)
    j += lengths[i];

  *symlistlen = 0;
  if (j == 0)
    return NULL;

  /* find all symbols, use a kind of hash for storage */
  HItem *symmap = (HItem*)safe_malloc(0x100, sizeof(HItem));
  if (!symmap) {
    *symlistlen = (size_t)(-1);
    return NULL;
  }
  /* this is an ugly memory allocation avoiding hack: most hash elements
   * will probably contain none or one symbols only so, when p->n is equal
   * to symmap, it means there're no symbols yet, afters insterting the
   * first one, p->n becomes normally NULL and then it behaves like an
   * usual singly linked list */
  for (size_t i = 0; i < 0x100; i++)
    symmap[i].n = symmap;
  for (size_t i = 0; i < n; i++) {
    const lev_wchar *stri = strings[i];
    for (size_t j = 0; j < lengths[i]; j++) {
      lev_wchar c = stri[j];
      int key = ((int)c + ((int)c >> 7)) & 0xff;
      HItem *p = symmap + key;
      if (p->n == symmap) {
        p->c = c;
        p->n = NULL;
        (*symlistlen)++;
        continue;
      }
      while (p->c != c && p->n != NULL)
        p = p->n;
      if (p->c != c) {
        p->n = (HItem*)malloc(sizeof(HItem));
        if (!p->n) {
          free_usymlist_hash(symmap);
          *symlistlen = (size_t)(-1);
          return NULL;
        }
        p = p->n;
        p->n = NULL;
        p->c = c;
        (*symlistlen)++;
      }
    }
  }
  /* create dense symbol table, so we can easily iterate over only characters
   * present in the strings */
  lev_wchar *symlist= (lev_wchar*)safe_malloc((*symlistlen), sizeof(lev_wchar));
  if (!symlist) {
    free_usymlist_hash(symmap);
    *symlistlen = (size_t)(-1);
    return NULL;
  }
  size_t pos = 0;
  for (size_t j = 0; j < 0x100; j++) {
    HItem *p = symmap + j;
    while (p != NULL && p->n != symmap) {
      symlist[pos++] = p->c;
      p = p->n;
    }
  }

  /* free memory */
  free_usymlist_hash(symmap);

  return symlist;
}

/**
 * lev_u_greedy_median:
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
lev_wchar*
lev_u_greedy_median(size_t n, const size_t *lengths,
                    const lev_wchar *strings[],
                    const double *weights,
                    size_t *medlength)
{
  size_t i;  /* usually iterates over strings (n) */
  size_t j;  /* usually iterates over characters */
  size_t len;  /* usually iterates over the approximate median string */
  lev_wchar *symlist;  /* list of symbols present in the strings,
                              we iterate over it insead of set of all
                              existing symbols */
  size_t symlistlen;  /* length of symlistle */
  size_t maxlen;  /* maximum input string length */
  size_t stoplen;  /* maximum tried median string length -- this is slightly
                      higher than maxlen, because the median string may be
                      longer than any of the input strings */
  size_t **rows;  /* Levenshtein matrix rows for each string, we need to keep
                     only one previous row to construct the current one */
  size_t *row;  /* a scratch buffer for new Levenshtein matrix row computation,
                   shared among all strings */
  lev_wchar *median;  /* the resulting approximate median string */
  double *mediandist;  /* the total distance of the best median string of
                          given length.  warning!  mediandist[0] is total
                          distance for empty string, while median[] itself
                          is normally zero-based */
  size_t bestlen;  /* the best approximate median string length */

  /* find all symbols */
  symlist = make_usymlist(n, lengths, strings, &symlistlen);
  if (!symlist) {
    *medlength = 0;
    if (symlistlen != 0)
      return NULL;
    else
      return (lev_wchar*)calloc(1, sizeof(lev_wchar));
  }

  /* allocate and initialize per-string matrix rows and a common work buffer */
  rows = (size_t**)safe_malloc(n, sizeof(size_t*));
  if (!rows) {
    free(symlist);
    return NULL;
  }
  maxlen = 0;
  for (i = 0; i < n; i++) {
    size_t *ri;
    size_t leni = lengths[i];
    if (leni > maxlen)
      maxlen = leni;
    ri = rows[i] = (size_t*)safe_malloc((leni + 1), sizeof(size_t));
    if (!ri) {
      for (j = 0; j < i; j++)
        free(rows[j]);
      free(rows);
      free(symlist);
      return NULL;
    }
    for (j = 0; j <= leni; j++)
      ri[j] = j;
  }
  stoplen = 2*maxlen + 1;
  row = (size_t*)safe_malloc((stoplen + 1), sizeof(size_t));
  if (!row) {
    for (j = 0; j < n; j++)
      free(rows[j]);
    free(rows);
    free(symlist);
    return NULL;
  }

  /* compute final cost of string of length 0 (empty string may be also
   * a valid answer) */
  median = (lev_wchar*)safe_malloc(stoplen, sizeof(lev_wchar));
  if (!median) {
    for (j = 0; j < n; j++)
      free(rows[j]);
    free(rows);
    free(row);
    free(symlist);
    return NULL;
  }
  mediandist = (double*)safe_malloc((stoplen + 1), sizeof(double));
  if (!mediandist) {
    for (j = 0; j < n; j++)
      free(rows[j]);
    free(rows);
    free(row);
    free(symlist);
    free(median);
    return NULL;
  }
  mediandist[0] = 0.0;
  for (i = 0; i < n; i++)
    mediandist[0] += (double)lengths[i] * weights[i];

  /* build up the approximate median string symbol by symbol
   * XXX: we actually exit on break below, but on the same condition */
  for (len = 1; len <= stoplen; len++) {
    lev_wchar symbol;
    double minminsum = LEV_INFINITY;
    row[0] = len;
    /* iterate over all symbols we may want to add */
    for (j = 0; j < symlistlen; j++) {
      double totaldist = 0.0;
      double minsum = 0.0;
      symbol = symlist[j];
      /* sum Levenshtein distances from all the strings, with given weights */
      for (i = 0; i < n; i++) {
        const lev_wchar *stri = strings[i];
        size_t *p = rows[i];
        size_t leni = lengths[i];
        size_t *end = rows[i] + leni;
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
        minsum += (double)min * weights[i];
        totaldist += (double)x * weights[i];
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
    for (i = 0; i < n; i++) {
      const lev_wchar *stri = strings[i];
      size_t *oldrow = rows[i];
      size_t leni = lengths[i];
      size_t k;
      /* compute a row of Levenshtein matrix */
      for (k = 1; k <= leni; k++) {
        size_t c1 = oldrow[k] + 1;
        size_t c2 = row[k - 1] + 1;
        size_t c3 = oldrow[k - 1] + (symbol != stri[k - 1]);
        row[k] = c2 > c3 ? c3 : c2;
        if (row[k] > c1)
          row[k] = c1;
      }
      memcpy(oldrow, row, (leni + 1)*sizeof(size_t));
    }
  }

  /* find the string with minimum total distance */
  bestlen = 0;
  for (len = 1; len <= stoplen; len++) {
    if (mediandist[len] < mediandist[bestlen])
      bestlen = len;
  }

  /* clean up */
  for (i = 0; i < n; i++)
    free(rows[i]);
  free(rows);
  free(row);
  free(symlist);
  free(mediandist);

  /* return result */
  {
    lev_wchar *result = (lev_wchar*)safe_malloc(bestlen, sizeof(lev_wchar));
    if (!result) {
      free(median);
      return NULL;
    }
    memcpy(result, median, bestlen*sizeof(lev_wchar));
    free(median);
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
static double
finish_udistance_computations(size_t len1, lev_wchar *string1,
                             size_t n, const size_t *lengths,
                             const lev_wchar **strings,
                             const double *weights, size_t **rows,
                             size_t *row)
{
  size_t *end;
  size_t offset;  /* row[0]; offset + len1 give together real len of string1 */
  double distsum = 0.0;  /* sum of distances */

  /* catch trivia case */
  if (len1 == 0) {
    for (size_t j = 0; j < n; j++)
      distsum += (double)rows[j][lengths[j]] * weights[j];
    return distsum;
  }

  /* iterate through the strings and sum the distances */
  for (size_t j = 0; j < n; j++) {
    size_t *rowi = rows[j];  /* current row */
    size_t leni = lengths[j];  /* current length */
    size_t len = len1;  /* temporary len1 for suffix stripping */
    const lev_wchar *stringi = strings[j];  /* current string */

    /* strip common suffix (prefix CAN'T be stripped) */
    while (len && leni && stringi[leni-1] == string1[len-1]) {
      len--;
      leni--;
    }

    /* catch trivial cases */
    if (len == 0) {
      distsum += (double)rowi[leni] * weights[j];
      continue;
    }
    offset = rowi[0];
    if (leni == 0) {
      distsum += (double)(offset + len) * weights[j];
      continue;
    }

    /* complete the matrix */
    memcpy(row, rowi, (leni + 1) * sizeof(size_t));
    end = row + leni;

    for (size_t i = 1; i <= len; i++) {
      size_t *p = row + 1;
      const lev_wchar char1 = string1[i - 1];
      const lev_wchar *char2p = stringi;
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
    distsum += weights[j] * (double)(*end);
  }

  return distsum;
}

/**
 * lev_u_median_improve:
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
lev_wchar*
lev_u_median_improve(size_t len, const lev_wchar *s,
                     size_t n, const size_t *lengths,
                     const lev_wchar *strings[],
                     const double *weights,
                     size_t *medlength)
{
  size_t i;  /* usually iterates over strings (n) */
  size_t j;  /* usually iterates over characters */
  size_t pos;  /* the position in the approximate median string we are
                  trying to change */
  lev_wchar *symlist;  /* list of symbols present in the strings,
                              we iterate over it insead of set of all
                              existing symbols */
  size_t symlistlen;  /* length of symlist */
  size_t maxlen;  /* maximum input string length */
  size_t stoplen;  /* maximum tried median string length -- this is slightly
                      higher than maxlen, because the median string may be
                      longer than any of the input strings */
  size_t **rows;  /* Levenshtein matrix rows for each string, we need to keep
                     only one previous row to construct the current one */
  size_t *row;  /* a scratch buffer for new Levenshtein matrix row computation,
                   shared among all strings */
  lev_wchar *median;  /* the resulting approximate median string */
  size_t medlen;  /* the current approximate median string length */
  double minminsum;  /* the current total distance sum */

  /* find all symbols */
  symlist = make_usymlist(n, lengths, strings, &symlistlen);
  if (!symlist) {
    *medlength = 0;
    if (symlistlen != 0)
      return NULL;
    else
      return (lev_wchar*)calloc(1, sizeof(lev_wchar));
  }

  /* allocate and initialize per-string matrix rows and a common work buffer */
  rows = (size_t**)safe_malloc(n, sizeof(size_t*));
  if (!rows) {
    free(symlist);
    return NULL;
  }
  maxlen = 0;
  for (i = 0; i < n; i++) {
    size_t *ri;
    size_t leni = lengths[i];
    if (leni > maxlen)
      maxlen = leni;
    ri = rows[i] = (size_t*)safe_malloc((leni + 1), sizeof(size_t));
    if (!ri) {
      for (j = 0; j < i; j++)
        free(rows[j]);
      free(rows);
      free(symlist);
      return NULL;
    }
    for (j = 0; j <= leni; j++)
      ri[j] = j;
  }
  stoplen = 2*maxlen + 1;
  row = (size_t*)safe_malloc((stoplen + 2), sizeof(size_t));
  if (!row) {
    for (j = 0; j < n; j++)
      free(rows[j]);
    free(rows);
    free(symlist);
    return NULL;
  }

  /* initialize median to given string */
  median = (lev_wchar*)safe_malloc((stoplen+1), sizeof(lev_wchar));
  if (!median) {
    for (j = 0; j < n; j++)
      free(rows[j]);
    free(rows);
    free(row);
    free(symlist);
    return NULL;
  }
  median++;  /* we need -1st element for insertions a pos 0 */
  medlen = len;
  memcpy(median, s, (medlen)*sizeof(lev_wchar));
  minminsum = finish_udistance_computations(medlen, median,
                                            n, lengths, strings,
                                            weights, rows, row);

  /* sequentially try perturbations on all positions */
  for (pos = 0; pos <= medlen; ) {
    lev_wchar orig_symbol, symbol;
    LevEditType operation;
    double sum;

    symbol = median[pos];
    operation = LEV_EDIT_KEEP;
    /* IF pos < medlength: FOREACH symbol: try to replace the symbol
     * at pos, if some lower the total distance, chooste the best */
    if (pos < medlen) {
      orig_symbol = median[pos];
      for (j = 0; j < symlistlen; j++) {
        if (symlist[j] == orig_symbol)
          continue;
        median[pos] = symlist[j];
        sum = finish_udistance_computations(medlen - pos, median + pos,
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
    for (j = 0; j < symlistlen; j++) {
      *(median + pos - 1) = symlist[j];
      sum = finish_udistance_computations(medlen - pos + 1, median + pos - 1,
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
      sum = finish_udistance_computations(medlen - pos - 1, median + pos + 1,
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
              (medlen - pos)*sizeof(lev_wchar));
      median[pos] = symbol;
      medlen++;
      break;

    case LEV_EDIT_DELETE:
      memmove(median+pos, median + pos+1,
              (medlen - pos-1)*sizeof(lev_wchar));
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
      for (i = 0; i < n; i++) {
        const lev_wchar *stri = strings[i];
        size_t *oldrow = rows[i];
        size_t leni = lengths[i];
        size_t k;
        /* compute a row of Levenshtein matrix */
        for (k = 1; k <= leni; k++) {
          size_t c1 = oldrow[k] + 1;
          size_t c2 = row[k - 1] + 1;
          size_t c3 = oldrow[k - 1] + (symbol != stri[k - 1]);
          row[k] = c2 > c3 ? c3 : c2;
          if (row[k] > c1)
            row[k] = c1;
        }
        memcpy(oldrow, row, (leni + 1)*sizeof(size_t));
      }
      pos++;
    }
  }

  /* clean up */
  for (i = 0; i < n; i++)
    free(rows[i]);
  free(rows);
  free(row);
  free(symlist);

  /* return result */
  {
    lev_wchar *result = (lev_wchar*)safe_malloc(medlen, sizeof(lev_wchar));
    if (!result) {
      median--;
      free(median);
      return NULL;
    }
    *medlength = medlen;
    memcpy(result, median, medlen*sizeof(lev_wchar));
    median--;
    free(median);
    return result;
  }
}
/* }}} */

/****************************************************************************
 *
 * Quick (voting) medians
 *
 ****************************************************************************/
/* {{{ */

/* compute the sets of symbols each string contains, and the set of symbols
 * in any of them (symset).  meanwhile, count how many different symbols
 * there are (used below for symlist).
 * the symset is passed as an argument to avoid its allocation and
 * deallocation when it's used in the caller too */
static std::vector<lev_byte>
make_symlistset(size_t n, const size_t *lengths,
                const lev_byte *strings[],
                std::unique_ptr<double[]>& symset)
{
  std::vector<lev_byte> symlist;
  memset(symset.get(), 0, 0x100*sizeof(double));  /* XXX: needs IEEE doubles?! */
  size_t symlistlen = 0;
  for (size_t i = 0; i < n; i++) {
    const lev_byte *stri = strings[i];
    for (size_t j = 0; j < lengths[i]; j++) {
      int c = stri[j];
      if (!symset[c]) {
        symlistlen++;
        symset[c] = 1.0;
      }
    }
  }
  if (!symlistlen)
    return symlist;

  /* create dense symbol table, so we can easily iterate over only characters
   * present in the strings */
  symlist.reserve(symlistlen);
  for (size_t j = 0; j < 0x100; j++) {
    if (symset[j])
      symlist.push_back((lev_byte)j);
  }

  return symlist;
}

lev_byte*
lev_quick_median(size_t n,
                 const size_t *lengths,
                 const lev_byte *strings[],
                 const double *weights,
                 size_t *medlength)
{
  /* first check whether the result would be an empty string 
   * and compute resulting string length */
  double ml = 0.0;
  double wl = 0.0;
  for (size_t i = 0; i < n; i++) {
    ml += (double)lengths[i] * weights[i];
    wl += weights[i];
  }
  if (wl == 0.0)
    return (lev_byte*)calloc(1, sizeof(lev_byte));
  ml = floor(ml/wl + 0.499999);
  size_t len = (size_t)ml;
  *medlength = len;
  if (!len)
    return (lev_byte*)calloc(1, sizeof(lev_byte));

  /* find the symbol set;
   * now an empty symbol set is really a failure */
  auto symset = std::make_unique<double[]>(0x100);

  auto symlist = make_symlistset(n, lengths, strings, symset);
  if (symlist.empty()) {
    return NULL;
  }

  lev_byte *median = (lev_byte*)safe_malloc(len, sizeof(lev_byte));
  if (!median)
    return NULL;

  for (size_t j = 0; j < len; j++) {
    /* clear the symbol probabilities */
    if (symlist.size() < 32) {
      for (size_t i = 0; i < symlist.size(); i++)
        symset[symlist[i]] = 0.0;
    }
    else
      memset(symset.get(), 0, 0x100*sizeof(double));

    /* let all strings vote */
    for (size_t i = 0; i < n; i++) {
      const lev_byte *stri = strings[i];
      double weighti = weights[i];
      size_t lengthi = lengths[i];
      double start = (double)lengthi / ml * (double)j;
      double end = start + (double)lengthi / ml;
      size_t istart = (size_t)floor(start);
      size_t iend = (size_t)ceil(end);

      /* rounding errors can overflow the buffer */
      if (iend > lengthi)
        iend = lengthi;

      for (size_t k = istart+1; k < iend; k++)
        symset[stri[k]] += weighti;
      symset[stri[istart]] += weighti * ((double)(1 + istart) - start);
      symset[stri[iend-1]] -= weighti * ((double)iend - end);
    }

    /* find the elected symbol */
    size_t k = symlist[0];
    for (size_t i = 1; i < symlist.size(); i++) {
      if (symset[symlist[i]] > symset[k])
        k = symlist[i];
    }
    median[j] = (lev_byte)k;
  }

  return median;
}

/* used internally in make_usymlistset */
typedef struct _HQItem HQItem;
struct _HQItem {
  lev_wchar c;
  double s;
  HQItem *n;
};

/* free usmylistset hash
 * this is a separate function because we need it in more than one place */
static void
free_usymlistset_hash(HQItem *symmap)
{
  for (size_t j = 0; j < 0x100; j++) {
    HQItem *p = symmap + j;
    if (p->n == symmap || p->n == NULL)
      continue;
    p = p->n;
    while (p) {
      HQItem *q = p;
      p = p->n;
      free(q);
    }
  }
  free(symmap);
}

/* compute the sets of symbols each string contains, and the set of symbols
 * in any of them (symset).  meanwhile, count how many different symbols
 * there are (used below for symlist).
 * the symset is passed as an argument to avoid its allocation and
 * deallocation when it's used in the caller too */
static lev_wchar*
make_usymlistset(size_t n, const size_t *lengths,
                 const lev_wchar *strings[], size_t *symlistlen,
                 HQItem *symmap)
{
  size_t len_sum = 0;
  for (size_t i = 0; i < n; i++)
    len_sum += lengths[i];

  *symlistlen = 0;
  if (len_sum == 0)
    return NULL;

  /* this is an ugly memory allocation avoiding hack: most hash elements
   * will probably contain none or one symbols only so, when p->n is equal
   * to symmap, it means there're no symbols yet, afters insterting the
   * first one, p->n becomes normally NULL and then it behaves like an
   * usual singly linked list */
  for (size_t i = 0; i < 0x100; i++)
    symmap[i].n = symmap;
  for (size_t i = 0; i < n; i++) {
    const lev_wchar *stri = strings[i];
    for (size_t j = 0; j < lengths[i]; j++) {
      lev_wchar c = stri[j];
      int key = ((int)c + ((int)c >> 7)) & 0xff;
      HQItem *p = symmap + key;
      if (p->n == symmap) {
        p->c = c;
        p->n = NULL;
        (*symlistlen)++;
        continue;
      }
      while (p->c != c && p->n != NULL)
        p = p->n;
      if (p->c != c) {
        p->n = (HQItem*)malloc(sizeof(HQItem));
        if (!p->n) {
          *symlistlen = (size_t)(-1);
          return NULL;
        }
        p = p->n;
        p->n = NULL;
        p->c = c;
        (*symlistlen)++;
      }
    }
  }
  /* create dense symbol table, so we can easily iterate over only characters
   * present in the strings */
  lev_wchar *symlist = (lev_wchar*)safe_malloc((*symlistlen), sizeof(lev_wchar));
  if (!symlist) {
    *symlistlen = (size_t)(-1);
    return NULL;
  }

  size_t pos = 0;
  for (size_t j = 0; j < 0x100; j++) {
    HQItem *p = symmap + j;
    while (p != NULL && p->n != symmap) {
      symlist[pos++] = p->c;
      p = p->n;
    }
  }

  return symlist;
}

lev_wchar*
lev_u_quick_median(size_t n,
                   const size_t *lengths,
                   const lev_wchar *strings[],
                   const double *weights,
                   size_t *medlength)
{
  /* first check whether the result would be an empty string 
   * and compute resulting string length */
  double ml = 0.0;
  double wl = 0.0;
  for (size_t i = 0; i < n; i++) {
    ml += (double)lengths[i] * weights[i];
    wl += weights[i];
  }
  if (wl == 0.0)
    return (lev_wchar*)calloc(1, sizeof(lev_wchar));
  ml = floor(ml/wl + 0.499999);
  size_t len = (size_t)ml;
  *medlength = len;
  if (!len)
    return (lev_wchar*)calloc(1, sizeof(lev_wchar));

  lev_wchar *median = (lev_wchar*)safe_malloc(len, sizeof(lev_wchar));
  if (!median)
    return NULL;

  /* find the symbol set;
   * now an empty symbol set is really a failure */
  HQItem *symmap = (HQItem*)safe_malloc(0x100, sizeof(HQItem));
  if (!symmap) {
    free(median);
    return NULL;
  }
  size_t symlistlen;
  lev_wchar *symlist = make_usymlistset(n, lengths, strings, &symlistlen, symmap);
  if (!symlist) {
    free(median);
    free_usymlistset_hash(symmap);
    return NULL;
  }

  for (size_t j = 0; j < len; j++) {
    /* clear the symbol probabilities */
    for (size_t i = 0; i < 0x100; i++) {
      HQItem *p = symmap + i;
      if (p->n == symmap)
        continue;
      while (p) {
        p->s = 0.0;
        p = p->n;
      }
    }

    /* let all strings vote */
    for (size_t i = 0; i < n; i++) {
      const lev_wchar *stri = strings[i];
      double weighti = weights[i];
      size_t lengthi = lengths[i];
      double start = (double)lengthi / ml * (double)j;
      double end = start + (double)lengthi / ml;
      size_t istart = (size_t)floor(start);
      size_t iend = (size_t)ceil(end);

      /* rounding errors can overflow the buffer */
      if (iend > lengthi)
        iend = lengthi;

      /* the inner part, including the complete last character */
      for (size_t k = istart+1; k < iend; k++) {
        int c = stri[k];
        int key = (c + (c >> 7)) & 0xff;
        HQItem *p = symmap + key;
        while (p->c != c)
          p = p->n;
        p->s += weighti;
      }
      /* the initial fraction */
      {
        int c = stri[istart];
        int key = (c + (c >> 7)) & 0xff;
        HQItem *p = symmap + key;
        while (p->c != c)
          p = p->n;
        p->s += weighti * ((double)(1 + istart) - start);
      }
      /* subtract what we counted from the last character but doesn't
       * actually belong here.
       * this strategy works also when istart+1 == iend (i.e., everything
       * happens inside a one character) */
      {
        int c = stri[iend-1];
        int key = (c + (c >> 7)) & 0xff;
        HQItem *p = symmap + key;
        while (p->c != c)
          p = p->n;
        p->s -= weighti * ((double)iend - end);
      }
    }

    /* find the elected symbol */
    {
      HQItem *max = NULL;

      for (size_t i = 0; i < 0x100; i++) {
        HQItem *p = symmap + i;
        if (p->n == symmap)
          continue;
        while (p) {
          if (!max || p->s > max->s)
            max = p;
          p = p->n;
        }
      }
      median[j] = max->c;
    }
  }

  free_usymlistset_hash(symmap);
  free(symlist);

  return median;
}
/* }}} */

/****************************************************************************
 *
 * Set, sequence distances
 *
 ****************************************************************************/
/* {{{ */

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
double
lev_edit_seq_distance(size_t n1, const size_t *lengths1,
                      const lev_byte *strings1[],
                      size_t n2, const size_t *lengths2,
                      const lev_byte *strings2[])
{
  /* make the inner cycle (i.e. strings2) the longer one */
  if (n1 > n2) {
    return lev_edit_seq_distance(n2, lengths2, strings2, n1, lengths1, strings1);
  }

  /* strip common prefix */
  while (n1 > 0 && n2 > 0
         && *lengths1 == *lengths2
         && memcmp(*strings1, *strings2,
                   *lengths1*sizeof(lev_byte)) == 0) {
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
                   lengths1[n1-1]*sizeof(lev_byte)) == 0) {
    n1--;
    n2--;
  }

  /* catch trivial cases */
  if (n1 == 0)
    return (double)n2;
  if (n2 == 0)
    return (double)n1;

  n1++;
  n2++;

  /* initalize first row */
  auto row = std::make_unique<double[]>(n2);
  double *end = row.get() + n2 - 1;
  for (size_t i = 0; i < n2; i++)
    row[i] = (double)i;

  /* go through the matrix and compute the costs.  yes, this is an extremely
   * obfuscated version, but also extremely memory-conservative and relatively
   * fast.  */
  for (size_t i = 1; i < n1; i++) {
    double *p = row.get() + 1;
    const lev_byte *str1 = strings1[i - 1];
    const size_t len1 = lengths1[i - 1];
    const lev_byte **str2p = strings2;
    const size_t *len2p = lengths2;
    double D = (double)i - 1.0;
    double x = (double)i;
    while (p <= end) {
      size_t l = len1 + *len2p;
      double q;
      if (l == 0)
        q = D;
      else {
        size_t d = rapidfuzz::string_metric::levenshtein(
          rapidfuzz::basic_string_view<lev_byte>(str1, len1),
          rapidfuzz::basic_string_view<lev_byte>(*str2p, *len2p),
          {1, 1, 2}
        );
        str2p++;
        len2p++;
        if (d == (size_t)(-1)) {
          return -1.0;
        }
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

  return *end;
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
double
lev_u_edit_seq_distance(size_t n1, const size_t *lengths1,
                        const lev_wchar *strings1[],
                        size_t n2, const size_t *lengths2,
                        const lev_wchar *strings2[])
{
  /* make the inner cycle (i.e. strings2) the longer one */
  if (n1 > n2) {
    return lev_u_edit_seq_distance(n2, lengths2, strings2, n1, lengths1, strings1);
  }

  /* strip common prefix */
  while (n1 > 0 && n2 > 0
         && *lengths1 == *lengths2
         && memcmp(*strings1, *strings2,
                   *lengths1*sizeof(lev_wchar)) == 0) {
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
                   lengths1[n1-1]*sizeof(lev_wchar)) == 0) {
    n1--;
    n2--;
  }

  /* catch trivial cases */
  if (n1 == 0)
    return (double)n2;
  if (n2 == 0)
    return (double)n1;

  n1++;
  n2++;

  /* initalize first row */
  auto row = std::make_unique<double[]>(n2);
  double *end = row.get() + n2 - 1;
  for (size_t i = 0; i < n2; i++)
    row[i] = (double)i;

  /* go through the matrix and compute the costs.  yes, this is an extremely
   * obfuscated version, but also extremely memory-conservative and relatively
   * fast.  */
  for (size_t i = 1; i < n1; i++) {
    double *p = row.get() + 1;
    const lev_wchar *str1 = strings1[i - 1];
    const size_t len1 = lengths1[i - 1];
    const lev_wchar **str2p = strings2;
    const size_t *len2p = lengths2;
    double D = (double)i - 1.0;
    double x = (double)i;
    while (p <= end) {
      size_t l = len1 + *len2p;
      double q;
      if (l == 0)
        q = D;
      else {
        size_t d = rapidfuzz::string_metric::levenshtein(
          rapidfuzz::basic_string_view<lev_wchar>(str1, len1),
          rapidfuzz::basic_string_view<lev_wchar>(*str2p, *len2p),
          {1, 1, 2}
        );
        str2p++;
        len2p++;
        if (d == (size_t)(-1)) {
          return -1.0;
        }
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

  double q = *end;
  return q;
}

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
double
lev_set_distance(size_t n1, const size_t *lengths1,
                 const lev_byte *strings1[],
                 size_t n2, const size_t *lengths2,
                 const lev_byte *strings2[])
{
  /* make the number of columns (n1) smaller than the number of rows */
  if (n1 > n2) {
    return lev_set_distance(n2, lengths2, strings2, n1, lengths1, strings1);
  }

  /* catch trivial cases */
  if (n1 == 0)
    return (double)n2;
  if (n2 == 0)
    return (double)n1;

  /* compute distances from each to each */
  double *dists = (double*)safe_malloc_3(n1, n2, sizeof(double));
  double* r = dists;
  if (!r)
    return -1.0;
  for (size_t i = 0; i < n2; i++) {
    size_t len2 = lengths2[i];
    const lev_byte *str2 = strings2[i];
    const size_t *len1p = lengths1;
    const lev_byte **str1p = strings1;
    for (size_t j = 0; j < n1; j++) {
      size_t l = len2 + *len1p;
      if (l == 0)
        *(r++) = 0.0;
      else {
        size_t d = rapidfuzz::string_metric::levenshtein(
          rapidfuzz::basic_string_view<lev_byte>(str2, len2),
          rapidfuzz::basic_string_view<lev_byte>(*str1p, *len1p),
          {1, 1, 2}
        );
        str1p++;
        len1p++;
        if (d == (size_t)(-1)) {
          free(r);
          return -1.0;
        }
        *(r++) = (double)d / (double)l;
      }
    }
  }

  /* find the optimal mapping between the two sets */
  size_t *map = munkers_blackman(n1, n2, dists);
  if (!map)
    return -1.0;

  /* sum the set distance */
  double sum = (double)(n2 - n1);
  for (size_t j = 0; j < n1; j++) {
    size_t l;
    size_t i = map[j];
    l = lengths1[j] + lengths2[i];
    if (l > 0) {
      size_t d = rapidfuzz::string_metric::levenshtein(
        rapidfuzz::basic_string_view<lev_byte>(strings1[j], lengths1[j]),
        rapidfuzz::basic_string_view<lev_byte>(strings2[i], lengths2[i]),
        {1, 1, 2}
      );
      if (d == (size_t)(-1)) {
        free(map);
        return -1.0;
      }
      sum += 2.0 * (double)d / (double)l;
    }
  }
  free(map);

  return sum;
}

/**
 * lev_u_set_distance:
 * @n1: The length of @lengths1 and @strings1.
 * @lengths1: The lengths of strings in @strings1.
 * @strings1: An array of strings that may contain NUL characters.
 * @n2: The length of @lengths2 and @strings2.
 * @lengths2: The lengths of strings in @strings2.
 * @strings2: An array of strings that may contain NUL characters.
 *
 * Finds the distance between string sets @strings1 and @strings2.
 *
 * The difference from lev_u_edit_seq_distance() is that order doesn't matter.
 * The optimal association of @strings1 and @strings2 is found first and
 * the similarity is computed for that.
 *
 * Uses sequential Munkers-Blackman algorithm.
 *
 * Returns: The distance of the two sets.
 **/
double
lev_u_set_distance(size_t n1, const size_t *lengths1,
                   const lev_wchar *strings1[],
                   size_t n2, const size_t *lengths2,
                   const lev_wchar *strings2[])
{
  /* make the number of columns (n1) smaller than the number of rows */
  if (n1 > n2) {
    return lev_u_set_distance(n2, lengths2, strings2, n1, lengths1, strings1);
  }

  /* catch trivial cases */
  if (n1 == 0)
    return (double)n2;
  if (n2 == 0)
    return (double)n1;

  /* compute distances from each to each */
  double* dists = (double*)safe_malloc_3(n1, n2, sizeof(double));
  double* r = dists;
  if (!r)
    return -1.0;
  for (size_t i = 0; i < n2; i++) {
    size_t len2 = lengths2[i];
    const lev_wchar *str2 = strings2[i];
    const size_t *len1p = lengths1;
    const lev_wchar **str1p = strings1;
    for (size_t j = 0; j < n1; j++) {
      size_t l = len2 + *len1p;
      if (l == 0)
        *(r++) = 0.0;
      else {
        size_t d = rapidfuzz::string_metric::levenshtein(
          rapidfuzz::basic_string_view<lev_wchar>(str2, len2),
          rapidfuzz::basic_string_view<lev_wchar>(*str1p, *len1p),
          {1, 1, 2}
        );
        str1p++;
        len1p++;
        if (d == (size_t)(-1)) {
          free(r);
          return -1.0;
        }
        *(r++) = (double)d / (double)l;
      }
    }
  }

  /* find the optimal mapping between the two sets */
  size_t *map = munkers_blackman(n1, n2, dists);
  if (!map)
    return -1.0;

  /* sum the set distance */
  double sum = (double)(n2 - n1);
  for (size_t j = 0; j < n1; j++) {
    size_t l;
    size_t i = map[j];
    l = lengths1[j] + lengths2[i];
    if (l > 0) {
      size_t d = rapidfuzz::string_metric::levenshtein(
        rapidfuzz::basic_string_view<lev_wchar>(strings1[j], lengths1[j]),
        rapidfuzz::basic_string_view<lev_wchar>(strings1[i], lengths1[i]),
        {1, 1, 2}
      );
      if (d == (size_t)(-1)) {
        free(map);
        return -1.0;
      }
      sum += 2.0 * (double)d / (double)l;
    }
  }
  free(map);

  return sum;
}

/*
 * Munkers-Blackman algorithm.
 */
static size_t*
munkers_blackman(size_t n1, size_t n2, double *dists)
{
  size_t i, j;
  /* allocate memory */
  std::vector<size_t> covc(n1);
  size_t* zstarc = (size_t*)calloc(n1, sizeof(size_t));
  if (!zstarc) {
    return NULL;
  }
  std::vector<size_t> covr(n2);
  std::vector<size_t> zstarr(n2);
  std::vector<size_t> zprimer(n2);

  /* step 0 (subtract minimal distance) and step 1 (find zeroes) */
  for (j = 0; j < n1; j++) {
    size_t minidx = 0;
    double *col = dists + j;
    double min = *col;
    double *p = col + n1;
    for (i = 1; i < n2; i++) {
      if (min > *p) {
        minidx = i;
        min = *p;
      }
      p += n1;
    }
    /* subtract */
    p = col;
    for (i = 0; i < n2; i++) {
      *p -= min;
      if (*p < LEV_EPSILON)
        *p = 0.0;
      p += n1;
    }
    /* star the zero, if possible */
    if (!zstarc[j] && !zstarr[minidx]) {
      zstarc[j] = minidx + 1;
      zstarr[minidx] = j + 1;
    }
    else {
      /* otherwise try to find some other */
      p = col;
      for (i = 0; i < n2; i++) {
        if (i != minidx && *p == 0.0
            && !zstarc[j] && !zstarr[i]) {
          zstarc[j] = i + 1;
          zstarr[i] = j + 1;
          break;
        }
        p += n1;
      }
    }
  }

  /* main */
  while (1) {
    /* step 2 (cover columns containing z*) */
    {
      size_t nc = 0;
      for (j = 0; j < n1; j++) {
        if (zstarc[j]) {
          covc[j] = 1;
          nc++;
        }
      }
      if (nc == n1)
        break;
    }

    /* step 3 (find uncovered zeroes) */
    while (1) {
      step_3:
      /* search uncovered matrix entries */
      for (j = 0; j < n1; j++) {
        double *p = dists + j;
        if (covc[j])
          continue;
        for (i = 0; i < n2; i++) {
          if (!covr[i] && *p == 0.0) {
            /* when a zero is found, prime it */
            zprimer[i] = j + 1;
            if (zstarr[i]) {
              /* if there's a z* in the same row,
               * uncover the column, cover the row and redo */
              covr[i] = 1;
              covc[zstarr[i] - 1] = 0;
              goto step_3;
            }
            /* if there's no z*,
             * we are at the end of our path an can convert z'
             * to z* */
            goto step_4;
          }
          p += n1;
        }
      }

      /* step 5 (new zero manufacturer)
       * we can't get here, unless no zero is found at all */
      {
        /* find the smallest uncovered entry */
        double min = LEV_INFINITY;
        for (j = 0; j < n1; j++) {
          double *p = dists + j;
          if (covc[j])
            continue;
          for (i = 0; i < n2; i++) {
            if (!covr[i] && min > *p) {
              min = *p;
            }
            p += n1;
          }
        }
        /* add it to all covered rows */
        for (i = 0; i < n2; i++) {
          double *p = dists + i*n1;
          if (!covr[i])
            continue;
          for (j = 0; j < n1; j++)
            *(p++) += min;
        }
        /* subtract if from all uncovered columns */
        for (j = 0; j < n1; j++) {
          double *p = dists + j;
          if (covc[j])
            continue;
          for (i = 0; i < n2; i++) {
            *p -= min;
            if (*p < LEV_EPSILON)
              *p = 0.0;
            p += n1;
          }
        }
      }
    }

    /* step 4 (increment the number of z*)
     * i is the row number (we get it from step 3) */
    step_4:
    i++;
    do {
      size_t x = i;

      i--;
      j = zprimer[i] - 1;  /* move to z' in the same row */
      zstarr[i] = j + 1;   /* mark it as z* in row buffer */
      i = zstarc[j];       /* move to z* in the same column */
      zstarc[j] = x;       /* mark the z' as being new z* */
    } while (i);

    std::fill(zprimer.begin(), zprimer.end(), 0);
    std::fill(covr.begin(), covr.end(), 0);
    std::fill(covc.begin(), covc.end(), 0);
  }

  for (j = 0; j < n1; j++)
    zstarc[j]--;
  return zstarc;
}
/* }}} */

/****************************************************************************
 *
 * Editops and other difflib-like stuff.
 *
 ****************************************************************************/
/* {{{ */

/**
 * lev_editops_valid:
 * @len1: The length of an eventual @ops source string.
 * @len2: The length of an eventual @ops destination string.
 * @n: The length of @ops.
 * @ops: An array of elementary edit operations.
 *
 * Checks whether @ops is consistent and applicable as a partial edit from a
 * string of length @len1 to a string of length @len2.
 *
 * Returns: true if ops are valid
 **/
bool lev_editops_valid(size_t len1, size_t len2,
                         size_t n, const LevEditOp *ops)
{
  if (!n)
    return true;

  /* check bounds */
  const LevEditOp *o = ops;
  for (size_t i = n; i; i--, o++) {
    if (o->type >= LEV_EDIT_LAST)
      return false;
    if (o->spos > len1 || o->dpos > len2)
      return false;
    if (o->spos == len1 && o->type != LEV_EDIT_INSERT)
      return false;
    if (o->dpos == len2 && o->type != LEV_EDIT_DELETE)
      return false;
  }

  /* check ordering */
  o = ops + 1;
  for (size_t i = n - 1; i; i--, o++, ops++) {
    if (o->spos < ops->spos || o->dpos < ops->dpos)
      return false;
  }

  return true;
}

/**
 * lev_opcodes_valid:
 * @len1: The length of an eventual @ops source string.
 * @len2: The length of an eventual @ops destination string.
 * @nb: The length of @bops.
 * @bops: An array of difflib block edit operation codes.
 *
 * Checks whether @bops is consistent and applicable as an edit from a
 * string of length @len1 to a string of length @len2.
 *
 * Returns: true if bops are valid
 **/
bool lev_opcodes_valid(size_t len1, size_t len2,
                         size_t nb, const LevOpCode *bops)
{
  if (!nb)
    return true;

  /* check completenes */
  if (bops->sbeg || bops->dbeg
      || bops[nb - 1].send != len1 || bops[nb - 1].dend != len2)
    return false;

  /* check bounds and block consistency */
  const LevOpCode *b = bops;
  for (size_t i = nb; i; i--, b++) {
    if (b->send > len1 || b->dend > len2)
      return false;
    switch (b->type) {
    case LEV_EDIT_KEEP:
    case LEV_EDIT_REPLACE:
      if (b->dend - b->dbeg != b->send - b->sbeg || b->dend == b->dbeg)
        return false;
      break;

    case LEV_EDIT_INSERT:
      if (b->dend - b->dbeg == 0 || b->send - b->sbeg != 0)
        return false;
      break;

    case LEV_EDIT_DELETE:
      if (b->send - b->sbeg == 0 || b->dend - b->dbeg != 0)
        return false;
      break;

    default:
      return false;
      break;
    }
  }

  /* check ordering */
  b = bops + 1;
  for (size_t i = nb - 1; i; i--, b++, bops++) {
    if (b->sbeg != bops->send || b->dbeg != bops->dend)
      return false;
  }

  return true;
}

/**
 * lev_editops_invert:
 * @n: The length of @ops.
 * @ops: An array of elementary edit operations.
 *
 * Inverts the sense of @ops.  It is modified in place.
 *
 * In other words, @ops becomes a valid partial edit for the original source
 * and destination strings with their roles exchanged.
 **/
void
lev_editops_invert(size_t n, LevEditOp *ops)
{
  for (size_t i = n; i; i--, ops++) {
    size_t z = ops->dpos;
    ops->dpos = ops->spos;
    ops->spos = z;
    if (ops->type & 2)
      ops->type = (LevEditType)(ops->type ^ 1);
  }
}

/**
 * lev_editops_apply:
 * @len1: The length of @string1.
 * @string1: A string of length @len1, may contain NUL characters.
 * @len2: The length of @string2.
 * @string2: A string of length @len2, may contain NUL characters.
 * @n: The size of @ops.
 * @ops: An array of elementary edit operations.
 * @len: Where the size of the resulting string should be stored.
 *
 * Applies a partial edit @ops from @string1 to @string2.
 *
 * NB: @ops is not checked for applicability.
 *
 * Returns: The result of the partial edit as a newly allocated string, its
 *          length is stored in @len.
 **/
lev_byte*
lev_editops_apply(size_t len1, const lev_byte *string1,
                  size_t len2, const lev_byte *string2,
                  size_t n, const LevEditOp *ops,
                  size_t *len)
{
  LEV_UNUSED(len2);

  /* this looks too complex for such a simple task, but note ops is not
   * a complete edit sequence, we have to be able to apply anything anywhere */
  lev_byte *dst = (lev_byte*)safe_malloc((n + len1), sizeof(lev_byte));
  if (!dst) {
    *len = (size_t)(-1);
    return NULL;
  }
  lev_byte *dpos = dst;
  const lev_byte *spos = string1;
  for (size_t i = n; i; i--, ops++) {
    /* XXX: this fine with gcc internal memcpy, but when memcpy is
     * actually a function, it may be pretty slow */
    size_t j = ops->spos - (size_t)(spos - string1) + (ops->type == LEV_EDIT_KEEP);
    if (j) {
      memcpy(dpos, spos, j*sizeof(lev_byte));
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
  size_t j = len1 - (size_t)(spos - string1);
  if (j) {
    memcpy(dpos, spos, j*sizeof(lev_byte));
    spos += j;
    dpos += j;
  }

  *len = (size_t)(dpos - dst);
  /* possible realloc failure is detected with *len != 0 */
  return (lev_byte*)realloc(dst, *len*sizeof(lev_byte));
}

/**
 * lev_u_editops_apply:
 * @len1: The length of @string1.
 * @string1: A string of length @len1, may contain NUL characters.
 * @len2: The length of @string2.
 * @string2: A string of length @len2, may contain NUL characters.
 * @n: The size of @ops.
 * @ops: An array of elementary edit operations.
 * @len: Where the size of the resulting string should be stored.
 *
 * Applies a partial edit @ops from @string1 to @string2.
 *
 * NB: @ops is not checked for applicability.
 *
 * Returns: The result of the partial edit as a newly allocated string, its
 *          length is stored in @len.
 **/
lev_wchar*
lev_u_editops_apply(size_t len1, const lev_wchar *string1,
                    size_t len2, const lev_wchar *string2,
                    size_t n, const LevEditOp *ops,
                    size_t *len)
{
  LEV_UNUSED(len2);

  /* this looks too complex for such a simple task, but note ops is not
   * a complete edit sequence, we have to be able to apply anything anywhere */
  lev_wchar *dst = (lev_wchar*)safe_malloc((n + len1), sizeof(lev_wchar));
  if (!dst) {
    *len = (size_t)(-1);
    return NULL;
  }
  lev_wchar *dpos = dst;
  const lev_wchar *spos = string1;
  for (size_t i = n; i; i--, ops++) {
    /* XXX: this is fine with gcc internal memcpy, but when memcpy is
     * actually a function, it may be pretty slow */
    size_t j = ops->spos - (size_t)(spos - string1) + (ops->type == LEV_EDIT_KEEP);
    if (j) {
      memcpy(dpos, spos, j*sizeof(lev_wchar));
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
  size_t j = len1 - (size_t)(spos - string1);
  if (j) {
    memcpy(dpos, spos, j*sizeof(lev_wchar));
    spos += j;
    dpos += j;
  }

  *len = (size_t)(dpos - dst);
  /* possible realloc failure is detected with *len != 0 */
  return (lev_wchar*)realloc(dst, *len*sizeof(lev_wchar));
}

/**
 * editops_from_cost_matrix:
 * @len1: The length of @string1.
 * @string1: A string of length @len1, may contain NUL characters.
 * @o1: The offset where the matrix starts from the start of @string1.
 * @len2: The length of @string2.
 * @string2: A string of length @len2, may contain NUL characters.
 * @o2: The offset where the matrix starts from the start of @string2.
 * @matrix: The cost matrix.
 * @n: Where the number of edit operations should be stored.
 *
 * Reconstructs the optimal edit sequence from the cost matrix @matrix.
 *
 * The matrix is freed.
 *
 * Returns: The optimal edit sequence, as a newly allocated array of
 *          elementary edit operations, it length is stored in @n.
 **/
static LevEditOp*
editops_from_cost_matrix(size_t len1, const lev_byte *string1, size_t off1,
                         size_t len2, const lev_byte *string2, size_t off2,
                         size_t *matrix, size_t *n)
{
  int dir = 0;
  *n = matrix[len1*len2 - 1];
  size_t pos = *n;
  if (!*n) {
    free(matrix);
    return NULL;
  }
  LevEditOp *ops = (LevEditOp*)safe_malloc((*n), sizeof(LevEditOp));
  if (!ops) {
    free(matrix);
    *n = (size_t)(-1);
    return NULL;
  }
  size_t i = len1 - 1;
  size_t j = len2 - 1;
  size_t *p = matrix + len1*len2 - 1;
  while (i || j) {
    /* prefer contiuning in the same direction */
    if (dir < 0 && j && *p == *(p - 1) + 1) {
      pos--;
      ops[pos].type = LEV_EDIT_INSERT;
      ops[pos].spos = i + off1;
      ops[pos].dpos = --j + off2;
      p--;
      continue;
    }
    if (dir > 0 && i && *p == *(p - len2) + 1) {
      pos--;
      ops[pos].type = LEV_EDIT_DELETE;
      ops[pos].spos = --i + off1;
      ops[pos].dpos = j + off2;
      p -= len2;
      continue;
    }
    if (i && j && *p == *(p - len2 - 1)
        && string1[i - 1] == string2[j - 1]) {
      /* don't be stupid like difflib, don't store LEV_EDIT_KEEP */
      i--;
      j--;
      p -= len2 + 1;
      dir = 0;
      continue;
    }
    if (i && j && *p == *(p - len2 - 1) + 1) {
      pos--;
      ops[pos].type = LEV_EDIT_REPLACE;
      ops[pos].spos = --i + off1;
      ops[pos].dpos = --j + off2;
      p -= len2 + 1;
      dir = 0;
      continue;
    }
    /* we cant't turn directly from -1 to 1, in this case it would be better
     * to go diagonally, but check it (dir == 0) */
    if (dir == 0 && j && *p == *(p - 1) + 1) {
      pos--;
      ops[pos].type = LEV_EDIT_INSERT;
      ops[pos].spos = i + off1;
      ops[pos].dpos = --j + off2;
      p--;
      dir = -1;
      continue;
    }
    if (dir == 0 && i && *p == *(p - len2) + 1) {
      pos--;
      ops[pos].type = LEV_EDIT_DELETE;
      ops[pos].spos = --i + off1;
      ops[pos].dpos = j + off2;
      p -= len2;
      dir = 1;
      continue;
    }
    /* coredump right now, later might be too late ;-) */
    assert("lost in the cost matrix" == NULL);
  }
  free(matrix);

  return ops;
}

/**
 * lev_editops_find:
 * @len1: The length of @string1.
 * @string1: A string of length @len1, may contain NUL characters.
 * @len2: The length of @string2.
 * @string2: A string of length @len2, may contain NUL characters.
 * @n: Where the number of edit operations should be stored.
 *
 * Find an optimal edit sequence from @string1 to @string2.
 *
 * When there's more than one optimal sequence, a one is arbitrarily (though
 * deterministically) chosen.
 *
 * Returns: The optimal edit sequence, as a newly allocated array of
 *          elementary edit operations, it length is stored in @n.
 *          It is normalized, i.e., keep operations are not included.
 **/
LevEditOp*
lev_editops_find(size_t len1, const lev_byte *string1,
                 size_t len2, const lev_byte *string2,
                 size_t *n)
{
  /* strip common prefix */
  size_t len1o = 0;
  while (len1 > 0 && len2 > 0 && *string1 == *string2) {
    len1--;
    len2--;
    string1++;
    string2++;
    len1o++;
  }
  size_t len2o = len1o;

  /* strip common suffix */
  while (len1 > 0 && len2 > 0 && string1[len1-1] == string2[len2-1]) {
    len1--;
    len2--;
  }
  len1++;
  len2++;

  /* initialize first row and column */
  size_t *matrix = (size_t*)safe_malloc_3(len1, len2, sizeof(size_t));
  if (!matrix) {
    *n = (size_t)(-1);
    return NULL;
  }
  for (size_t i = 0; i < len2; i++)
    matrix[i] = i;
  for (size_t i = 1; i < len1; i++)
    matrix[len2*i] = i;

  /* find the costs and fill the matrix */
  for (size_t i = 1; i < len1; i++) {
    size_t *prev = matrix + (i - 1)*len2;
    size_t *p = matrix + i*len2;
    size_t *end = p + len2 - 1;
    const lev_byte char1 = string1[i - 1];
    const lev_byte *char2p = string2;
    size_t x = i;
    p++;
    while (p <= end) {
      size_t c3 = *(prev++) + (char1 != *(char2p++));
      x++;
      if (x > c3)
        x = c3;
      c3 = *prev + 1;
      if (x > c3)
        x = c3;
      *(p++) = x;
    }
  }

  /* find the way back */
  return editops_from_cost_matrix(len1, string1, len1o,
                                  len2, string2, len2o,
                                  matrix, n);
}

/**
 * ueditops_from_cost_matrix:
 * @len1: The length of @string1.
 * @string1: A string of length @len1, may contain NUL characters.
 * @o1: The offset where the matrix starts from the start of @string1.
 * @len2: The length of @string2.
 * @string2: A string of length @len2, may contain NUL characters.
 * @o2: The offset where the matrix starts from the start of @string2.
 * @matrix: The cost matrix.
 * @n: Where the number of edit operations should be stored.
 *
 * Reconstructs the optimal edit sequence from the cost matrix @matrix.
 *
 * The matrix is freed.
 *
 * Returns: The optimal edit sequence, as a newly allocated array of
 *          elementary edit operations, it length is stored in @n.
 **/
static LevEditOp*
ueditops_from_cost_matrix(size_t len1, const lev_wchar *string1, size_t o1,
                          size_t len2, const lev_wchar *string2, size_t o2,
                          size_t *matrix, size_t *n)
{
  int dir = 0;

  *n = matrix[len1*len2 - 1];
  size_t pos = *n;
  if (!*n) {
    free(matrix);
    return NULL;
  }
  LevEditOp *ops = (LevEditOp*)safe_malloc((*n), sizeof(LevEditOp));
  if (!ops) {
    free(matrix);
    *n = (size_t)(-1);
    return NULL;
  }
  size_t i = len1 - 1;
  size_t j = len2 - 1;
  size_t *p = matrix + len1*len2 - 1;
  while (i || j) {
    /* prefer contiuning in the same direction */
    if (dir < 0 && j && *p == *(p - 1) + 1) {
      pos--;
      ops[pos].type = LEV_EDIT_INSERT;
      ops[pos].spos = i + o1;
      ops[pos].dpos = --j + o2;
      p--;
      continue;
    }
    if (dir > 0 && i && *p == *(p - len2) + 1) {
      pos--;
      ops[pos].type = LEV_EDIT_DELETE;
      ops[pos].spos = --i + o1;
      ops[pos].dpos = j + o2;
      p -= len2;
      continue;
    }
    if (i && j && *p == *(p - len2 - 1)
        && string1[i - 1] == string2[j - 1]) {
      /* don't be stupid like difflib, don't store LEV_EDIT_KEEP */
      i--;
      j--;
      p -= len2 + 1;
      dir = 0;
      continue;
    }
    if (i && j && *p == *(p - len2 - 1) + 1) {
      pos--;
      ops[pos].type = LEV_EDIT_REPLACE;
      ops[pos].spos = --i + o1;
      ops[pos].dpos = --j + o2;
      p -= len2 + 1;
      dir = 0;
      continue;
    }
    /* we cant't turn directly from -1 to 1, in this case it would be better
     * to go diagonally, but check it (dir == 0) */
    if (dir == 0 && j && *p == *(p - 1) + 1) {
      pos--;
      ops[pos].type = LEV_EDIT_INSERT;
      ops[pos].spos = i + o1;
      ops[pos].dpos = --j + o2;
      p--;
      dir = -1;
      continue;
    }
    if (dir == 0 && i && *p == *(p - len2) + 1) {
      pos--;
      ops[pos].type = LEV_EDIT_DELETE;
      ops[pos].spos = --i + o1;
      ops[pos].dpos = j + o2;
      p -= len2;
      dir = 1;
      continue;
    }
    /* coredump right now, later might be too late ;-) */
    assert("lost in the cost matrix" == NULL);
  }
  free(matrix);

  return ops;
}

/**
 * lev_u_editops_find:
 * @len1: The length of @string1.
 * @string1: A string of length @len1, may contain NUL characters.
 * @len2: The length of @string2.
 * @string2: A string of length @len2, may contain NUL characters.
 * @n: Where the number of edit operations should be stored.
 *
 * Find an optimal edit sequence from @string1 to @string2.
 *
 * When there's more than one optimal sequence, a one is arbitrarily (though
 * deterministically) chosen.
 *
 * Returns: The optimal edit sequence, as a newly allocated array of
 *          elementary edit operations, it length is stored in @n.
 *          It is normalized, i.e., keep operations are not included.
 **/
LevEditOp*
lev_u_editops_find(size_t len1, const lev_wchar *string1,
                   size_t len2, const lev_wchar *string2,
                   size_t *n)
{
  /* strip common prefix */
  size_t len1o = 0;
  while (len1 > 0 && len2 > 0 && *string1 == *string2) {
    len1--;
    len2--;
    string1++;
    string2++;
    len1o++;
  }
  size_t len2o = len1o;

  /* strip common suffix */
  while (len1 > 0 && len2 > 0 && string1[len1-1] == string2[len2-1]) {
    len1--;
    len2--;
  }
  len1++;
  len2++;

  /* initalize first row and column */
  size_t *matrix = (size_t*)safe_malloc_3(len1, len2, sizeof(size_t));
  if (!matrix) {
    *n = (size_t)(-1);
    return NULL;
  }
  for (size_t i = 0; i < len2; i++)
    matrix[i] = i;
  for (size_t i = 1; i < len1; i++)
    matrix[len2*i] = i;

  /* find the costs and fill the matrix */
  for (size_t i = 1; i < len1; i++) {
    size_t *prev = matrix + (i - 1)*len2;
    size_t *p = matrix + i*len2;
    size_t *end = p + len2 - 1;
    const lev_wchar char1 = string1[i - 1];
    const lev_wchar *char2p = string2;
    size_t x = i;
    p++;
    while (p <= end) {
      size_t c3 = *(prev++) + (char1 != *(char2p++));
      x++;
      if (x > c3)
        x = c3;
      c3 = *prev + 1;
      if (x > c3)
        x = c3;
      *(p++) = x;
    }
  }

  /* find the way back */
  return ueditops_from_cost_matrix(len1, string1, len1o,
                                   len2, string2, len2o,
                                   matrix, n);
}

/**
 * lev_opcodes_to_editops:
 * @nb: The length of @bops.
 * @bops: An array of difflib block edit operation codes.
 * @n: Where the number of edit operations should be stored.
 * @keepkeep: If nonzero, keep operations will be included.  Otherwise the
 *            result will be normalized, i.e. without any keep operations.
 *
 * Converts difflib block operation codes to elementary edit operations.
 *
 * Returns: The converted edit operatiosn, as a newly allocated array; its
 *          size is stored in @n.
 **/
LevEditOp*
lev_opcodes_to_editops(size_t nb, const LevOpCode *bops,
                       size_t *n, int keepkeep)
{
  /* compute the number of atomic operations */
  *n = 0;
  if (!nb)
    return NULL;

  const LevOpCode *b = bops;
  if (keepkeep) {
    for (size_t i = nb; i; i--, b++) {
      size_t sd = b->send - b->sbeg;
      size_t dd = b->dend - b->dbeg;
      *n += (sd > dd ? sd : dd);
    }
  }
  else {
    for (size_t i = nb; i; i--, b++) {
      size_t sd = b->send - b->sbeg;
      size_t dd = b->dend - b->dbeg;
      *n += (b->type != LEV_EDIT_KEEP ? (sd > dd ? sd : dd) : 0);
    }
  }

  /* convert */
  LevEditOp *ops = (LevEditOp*)safe_malloc((*n), sizeof(LevEditOp));
  LevEditOp *o = ops;
  if (!ops) {
    *n = (size_t)(-1);
    return NULL;
  }
  b = bops;
  for (size_t i = nb; i; i--, b++) {
    switch (b->type) {
    case LEV_EDIT_KEEP:
      if (keepkeep) {
        for (size_t j = 0; j < b->send - b->sbeg; j++, o++) {
          o->type = LEV_EDIT_KEEP;
          o->spos = b->sbeg + j;
          o->dpos = b->dbeg + j;
        }
      }
      break;

    case LEV_EDIT_REPLACE:
      for (size_t j = 0; j < b->send - b->sbeg; j++, o++) {
        o->type = LEV_EDIT_REPLACE;
        o->spos = b->sbeg + j;
        o->dpos = b->dbeg + j;
      }
      break;

    case LEV_EDIT_DELETE:
      for (size_t j = 0; j < b->send - b->sbeg; j++, o++) {
        o->type = LEV_EDIT_DELETE;
        o->spos = b->sbeg + j;
        o->dpos = b->dbeg;
      }
      break;

    case LEV_EDIT_INSERT:
      for (size_t j = 0; j < b->dend - b->dbeg; j++, o++) {
        o->type = LEV_EDIT_INSERT;
        o->spos = b->sbeg;
        o->dpos = b->dbeg + j;
      }
      break;

      default:
      break;
    }
  }
  assert((size_t)(o - ops) == *n);

  return ops;
}

/**
 * lev_editops_to_opcodes:
 * @n: The size of @ops.
 * @ops: An array of elementary edit operations.
 * @nb: Where the number of difflib block operation codes should be stored.
 * @len1: The length of the source string.
 * @len2: The length of the destination string.
 *
 * Converts elementary edit operations to difflib block operation codes.
 *
 * Note the string lengths are necessary since difflib doesn't allow omitting
 * keep operations.
 *
 * Returns: The converted block operation codes, as a newly allocated array;
 *          its length is stored in @nb.
 **/
LevOpCode*
lev_editops_to_opcodes(size_t n, const LevEditOp *ops, size_t *nb,
                       size_t len1, size_t len2)
{
  /* compute the number of blocks */
  size_t nbl = 0;
  const LevEditOp *o = ops;
  size_t spos = 0;
  size_t dpos = 0;
  LevEditType type = LEV_EDIT_KEEP;
  for (size_t i = n; i; ) {
    /* simply pretend there are no keep blocks */
    while (o->type == LEV_EDIT_KEEP && --i)
      o++;
    if (!i)
      break;
    if (spos < o->spos || dpos < o->dpos) {
      nbl++;
      spos = o->spos;
      dpos = o->dpos;
    }
    nbl++;
    type = o->type;
    switch (type) {
    case LEV_EDIT_REPLACE:
      do {
        spos++;
        dpos++;
        i--;
        o++;
      } while (i && o->type == type && spos == o->spos && dpos == o->dpos);
      break;

    case LEV_EDIT_DELETE:
      do {
        spos++;
        i--;
        o++;
      } while (i && o->type == type && spos == o->spos && dpos == o->dpos);
      break;

    case LEV_EDIT_INSERT:
      do {
        dpos++;
        i--;
        o++;
      } while (i && o->type == type && spos == o->spos && dpos == o->dpos);
      break;

    default:
      break;
    }
  }
  if (spos < len1 || dpos < len2)
    nbl++;

  /* convert */
  LevOpCode *bops = (LevOpCode*)safe_malloc(nbl, sizeof(LevOpCode));
  LevOpCode *b = bops;
  if (!bops) {
    *nb = (size_t)(-1);
    return NULL;
  }
  o = ops;
  spos = dpos = 0;
  type = LEV_EDIT_KEEP;
  for (size_t i = n; i; ) {
    /* simply pretend there are no keep blocks */
    while (o->type == LEV_EDIT_KEEP && --i)
      o++;
    if (!i)
      break;
    b->sbeg = spos;
    b->dbeg = dpos;
    if (spos < o->spos || dpos < o->dpos) {
      b->type = LEV_EDIT_KEEP;
      spos = b->send = o->spos;
      dpos = b->dend = o->dpos;
      b++;
      b->sbeg = spos;
      b->dbeg = dpos;
    }
    type = o->type;
    switch (type) {
      case LEV_EDIT_REPLACE:
      do {
        spos++;
        dpos++;
        i--;
        o++;
      } while (i && o->type == type && spos == o->spos && dpos == o->dpos);
      break;

      case LEV_EDIT_DELETE:
      do {
        spos++;
        i--;
        o++;
      } while (i && o->type == type && spos == o->spos && dpos == o->dpos);
      break;

      case LEV_EDIT_INSERT:
      do {
        dpos++;
        i--;
        o++;
      } while (i && o->type == type && spos == o->spos && dpos == o->dpos);
      break;

      default:
      break;
    }
    b->type = type;
    b->send = spos;
    b->dend = dpos;
    b++;
  }
  if (spos < len1 || dpos < len2) {
    assert(len1 - spos == len2 - dpos);
    b->type = LEV_EDIT_KEEP;
    b->sbeg = spos;
    b->dbeg = dpos;
    b->send = len1;
    b->dend = len2;
    b++;
  }
  assert((size_t)(b - bops) == nbl);

  *nb = nbl;
  return bops;
}

/**
 * lev_opcodes_apply:
 * @len1: The length of the source string.
 * @string1: A string of length @len1, may contain NUL characters.
 * @len2: The length of the destination string.
 * @string2: A string of length @len2, may contain NUL characters.
 * @nb: The length of @bops.
 * @bops: An array of difflib block edit operation codes.
 * @len: Where the size of the resulting string should be stored.
 *
 * Applies a sequence of difflib block operations to a string.
 *
 * NB: @bops is not checked for applicability.
 *
 * Returns: The result of the edit as a newly allocated string, its length
 *          is stored in @len.
 **/
lev_byte*
lev_opcodes_apply(size_t len1, const lev_byte *string1,
                  size_t len2, const lev_byte *string2,
                  size_t nb, const LevOpCode *bops,
                  size_t *len)
{
  /* this looks too complex for such a simple task, but note ops is not
   * a complete edit sequence, we have to be able to apply anything anywhere */
  lev_byte *dst = (lev_byte*)safe_malloc((len1 + len2), sizeof(lev_byte));
  lev_byte *dpos = dst;
  if (!dst) {
    *len = (size_t)(-1);
    return NULL;
  }
  const lev_byte *spos = string1;
  for (size_t i = nb; i; i--, bops++) {
    switch (bops->type) {
    case LEV_EDIT_INSERT:
    case LEV_EDIT_REPLACE:
      memcpy(dpos, string2 + bops->dbeg,
             (bops->dend - bops->dbeg)*sizeof(lev_byte));
      break;

    case LEV_EDIT_KEEP:
      memcpy(dpos, string1 + bops->sbeg,
             (bops->send - bops->sbeg)*sizeof(lev_byte));
      break;

    default:
      break;
    }
    spos += bops->send - bops->sbeg;
    dpos += bops->dend - bops->dbeg;
  }

  *len = (size_t)(dpos - dst);
  /* possible realloc failure is detected with *len != 0 */
  return (lev_byte*)realloc(dst, *len*sizeof(lev_byte));
}

/**
 * lev_u_opcodes_apply:
 * @len1: The length of the source string.
 * @string1: A string of length @len1, may contain NUL characters.
 * @len2: The length of the destination string.
 * @string2: A string of length @len2, may contain NUL characters.
 * @nb: The length of @bops.
 * @bops: An array of difflib block edit operation codes.
 * @len: Where the size of the resulting string should be stored.
 *
 * Applies a sequence of difflib block operations to a string.
 *
 * NB: @bops is not checked for applicability.
 *
 * Returns: The result of the edit as a newly allocated string, its length
 *          is stored in @len.
 **/
lev_wchar*
lev_u_opcodes_apply(size_t len1, const lev_wchar *string1,
                    size_t len2, const lev_wchar *string2,
                    size_t nb, const LevOpCode *bops,
                    size_t *len)
{
  /* this looks too complex for such a simple task, but note ops is not
   * a complete edit sequence, we have to be able to apply anything anywhere */
  lev_wchar *dst = (lev_wchar*)safe_malloc((len1 + len2), sizeof(lev_wchar));
  lev_wchar *dpos = dst;
  if (!dst) {
    *len = (size_t)(-1);
    return NULL;
  }
  const lev_wchar *spos = string1;
  for (size_t i = nb; i; i--, bops++) {
    switch (bops->type) {
    case LEV_EDIT_INSERT:
    case LEV_EDIT_REPLACE:
      memcpy(dpos, string2 + bops->dbeg,
             (bops->dend - bops->dbeg)*sizeof(lev_wchar));
      break;

    case LEV_EDIT_KEEP:
      memcpy(dpos, string1 + bops->sbeg,
             (bops->send - bops->sbeg)*sizeof(lev_wchar));
      break;

    default:
      break;
    }
    spos += bops->send - bops->sbeg;
    dpos += bops->dend - bops->dbeg;
  }

  *len = (size_t)(dpos - dst);
  /* possible realloc failure is detected with *len != 0 */
  return (lev_wchar*)realloc(dst, *len*sizeof(lev_wchar));
}

/**
 * lev_opcodes_invert:
 * @nb: The length of @ops.
 * @bops: An array of difflib block edit operation codes.
 *
 * Inverts the sense of @ops.  It is modified in place.
 *
 * In other words, @ops becomes a partial edit for the original source
 * and destination strings with their roles exchanged.
 **/
void
lev_opcodes_invert(size_t nb, LevOpCode *bops)
{
  for (size_t i = nb; i; i--, bops++) {
    size_t z = bops->dbeg;
    bops->dbeg = bops->sbeg;
    bops->sbeg = z;
    z = bops->dend;
    bops->dend = bops->send;
    bops->send = z;
    if (bops->type & 2)
      bops->type = (LevEditType)(bops->type ^ 1);
  }
}

/**
 * lev_editops_matching_blocks:
 * @len1: The length of the source string.
 * @len2: The length of the destination string.
 * @n: The size of @ops.
 * @ops: An array of elementary edit operations.
 * @nmblocks: Where the number of matching block should be stored.
 *
 * Computes the matching block corresponding to an optimal edit @ops.
 *
 * Returns: The matching blocks as a newly allocated array, it length is
 *          stored in @nmblocks.
 **/
LevMatchingBlock*
lev_editops_matching_blocks(size_t len1,
                            size_t len2,
                            size_t n,
                            const LevEditOp *ops,
                            size_t *nmblocks)
{
  /* compute the number of matching blocks */
  size_t nmb = 0;
  const LevEditOp *o = ops;
  size_t spos = 0;
  size_t dpos = 0;
  LevEditType type = LEV_EDIT_KEEP;
  for (size_t i = n; i; ) {
    /* simply pretend there are no keep blocks */
    while (o->type == LEV_EDIT_KEEP && --i)
      o++;
    if (!i)
      break;
    if (spos < o->spos || dpos < o->dpos) {
      nmb++;
      spos = o->spos;
      dpos = o->dpos;
    }
    type = o->type;
    switch (type) {
      case LEV_EDIT_REPLACE:
      do {
        spos++;
        dpos++;
        i--;
        o++;
      } while (i && o->type == type && spos == o->spos && dpos == o->dpos);
      break;

      case LEV_EDIT_DELETE:
      do {
        spos++;
        i--;
        o++;
      } while (i && o->type == type && spos == o->spos && dpos == o->dpos);
      break;

      case LEV_EDIT_INSERT:
      do {
        dpos++;
        i--;
        o++;
      } while (i && o->type == type && spos == o->spos && dpos == o->dpos);
      break;

      default:
      break;
    }
  }
  if (spos < len1 || dpos < len2)
    nmb++;

  /* fill the info */
  LevMatchingBlock *mblocks = (LevMatchingBlock*)safe_malloc(nmb, sizeof(LevMatchingBlock));
  LevMatchingBlock *mb = mblocks;
  if (!mblocks) {
    *nmblocks = (size_t)(-1);
    return NULL;
  }
  o = ops;
  spos = dpos = 0;
  type = LEV_EDIT_KEEP;
  for (size_t i = n; i; ) {
    /* simply pretend there are no keep blocks */
    while (o->type == LEV_EDIT_KEEP && --i)
      o++;
    if (!i)
      break;
    if (spos < o->spos || dpos < o->dpos) {
      mb->spos = spos;
      mb->dpos = dpos;
      mb->len = o->spos - spos;
      spos = o->spos;
      dpos = o->dpos;
      mb++;
    }
    type = o->type;
    switch (type) {
      case LEV_EDIT_REPLACE:
      do {
        spos++;
        dpos++;
        i--;
        o++;
      } while (i && o->type == type && spos == o->spos && dpos == o->dpos);
      break;

      case LEV_EDIT_DELETE:
      do {
        spos++;
        i--;
        o++;
      } while (i && o->type == type && spos == o->spos && dpos == o->dpos);
      break;

      case LEV_EDIT_INSERT:
      do {
        dpos++;
        i--;
        o++;
      } while (i && o->type == type && spos == o->spos && dpos == o->dpos);
      break;

      default:
      break;
    }
  }
  if (spos < len1 || dpos < len2) {
    assert(len1 - spos == len2 - dpos);
    mb->spos = spos;
    mb->dpos = dpos;
    mb->len = len1 - spos;
    mb++;
  }
  assert((size_t)(mb - mblocks) == nmb);

  *nmblocks = nmb;
  return mblocks;
}

/**
 * lev_opcodes_matching_blocks:
 * @len1: The length of the source string.
 * @len2: The length of the destination string.
 * @nb: The size of @bops.
 * @bops: An array of difflib block edit operation codes.
 * @nmblocks: Where the number of matching block should be stored.
 *
 * Computes the matching block corresponding to an optimal edit @bops.
 *
 * Returns: The matching blocks as a newly allocated array, it length is
 *          stored in @nmblocks.
 **/
LevMatchingBlock*
lev_opcodes_matching_blocks(size_t len1,
                            size_t len2,
                            size_t nb,
                            const LevOpCode *bops,
                            size_t *nmblocks)
{
  LEV_UNUSED(len2);

  /* compute the number of matching blocks */
  size_t nmb = 0;
  const LevOpCode *b = bops;
  for (size_t i = nb; i; i--, b++) {
    if (b->type == LEV_EDIT_KEEP) {
      nmb++;
      /* adjacent KEEP blocks -- we never produce it, but... */
      while (i && b->type == LEV_EDIT_KEEP) {
        i--;
        b++;
      }
      if (!i)
        break;
    }
  }

  /* convert */
  LevMatchingBlock *mblocks = (LevMatchingBlock*)safe_malloc(nmb, sizeof(LevMatchingBlock));
  LevMatchingBlock *mb = mblocks;
  if (!mblocks) {
    *nmblocks = (size_t)(-1);
    return NULL;
  }
  b = bops;
  for (size_t i = nb; i; i--, b++) {
    if (b->type == LEV_EDIT_KEEP) {
      mb->spos = b->sbeg;
      mb->dpos = b->dbeg;
      /* adjacent KEEP blocks -- we never produce it, but... */
      while (i && b->type == LEV_EDIT_KEEP) {
        i--;
        b++;
      }
      if (!i) {
        mb->len = len1 - mb->spos;
        mb++;
        break;
      }
      mb->len = b->sbeg - mb->spos;
      mb++;
    }
  }
  assert((size_t)(mb - mblocks) == nmb);

  *nmblocks = nmb;
  return mblocks;
}

/**
 * lev_editops_normalize:
 * @n: The size of @ops.
 * @ops: An array of elementary edit operations.
 * @nnorm: Where to store then length of the normalized array.
 *
 * Normalizes a list of edit operations to contain no keep operations.
 *
 * Note it returns %NULL for an empty array.
 *
 * Returns: A newly allocated array of normalized edit operations, its length
 *          is stored to @nnorm.
 **/
LevEditOp*
lev_editops_normalize(size_t n,
                      const LevEditOp *ops,
                      size_t *nnorm)
{
  if (!n || !ops) {
    *nnorm = 0;
    return NULL;
  }

  size_t nx = 0;
  const LevEditOp *o = ops;
  for (size_t i = n; i; i--, o++)
    nx += (o->type == LEV_EDIT_KEEP);

  *nnorm = n - nx;
  if (!*nnorm)
    return NULL;

  LevEditOp *opsnorm = (LevEditOp*)safe_malloc((n - nx), sizeof(LevEditOp));
  LevEditOp *on = opsnorm;
  o = ops;
  for (size_t i = n; i; i--, o++) {
    if (o->type == LEV_EDIT_KEEP)
      continue;
    memcpy(on++, o, sizeof(LevEditOp));
  }

  return opsnorm;
}

/**
 * lev_editops_subtract:
 * @n: The size of @ops.
 * @ops: An array of elementary edit operations.
 * @ns: The size of @sub.
 * @sub: A subsequence (ordered subset) of @ops
 * @nrem: Where to store then length of the remainder array.
 *
 * Subtracts a subsequence of elementary edit operations from a sequence.
 *
 * The remainder is a sequence that, applied to result of application of @sub,
 * gives the same final result as application of @ops to original string.
 *
 * Returns: A newly allocated array of normalized edit operations, its length
 *          is stored to @nrem.  It is always normalized, i.e, without any
 *          keep operations.  On failure, %NULL is returned.
 **/
LevEditOp*
lev_editops_subtract(size_t n,
                     const LevEditOp *ops,
                     size_t ns,
                     const LevEditOp *sub,
                     size_t *nrem)
{
    static const int shifts[] = { 0, 0, 1, -1 };

    /* compute remainder size */
    *nrem = (size_t)-1;
    size_t nr = 0;
    size_t nn = 0;
    for (size_t i = 0; i < n; i++) {
        if (ops[i].type != LEV_EDIT_KEEP)
            nr++;
    }
    for (size_t i = 0; i < ns; i++) {
        if (sub[i].type != LEV_EDIT_KEEP)
            nn++;
    }
    if (nn > nr)
        return NULL;
    nr -= nn;

    /* subtract */
    /* we could simply return NULL when nr == 0, but then it would be possible
     * to subtract *any* sequence of the right length to get an empty sequence
     * -- clrealy incorrectly; so we have to scan the list to check */
    LevEditOp *rem = nr ? (LevEditOp*)safe_malloc(nr, sizeof(LevEditOp)) : NULL;
    size_t j = nn = 0;
    int shift = 0;
    for (size_t i = 0; i < ns; i++) {
        while ((ops[j].spos != sub[i].spos
                || ops[j].dpos != sub[i].dpos
                || ops[j].type != sub[i].type)
               && j < n) {
            if (ops[j].type != LEV_EDIT_KEEP) {
                rem[nn] = ops[j];
                rem[nn].spos = (size_t)((int)rem[nn].spos + shift);
                nn++;
            }
            j++;
        }
        if (j == n) {
            free(rem);
            return NULL;
        }

        shift += shifts[sub[i].type];
        j++;
    }

    while (j < n) {
        if (ops[j].type != LEV_EDIT_KEEP) {
            rem[nn] = ops[j];
            rem[nn].spos = (size_t)((int)rem[nn].spos + shift);
            nn++;
        }
        j++;
    }
    assert(nn == nr);

    *nrem = nr;
    return rem;
}
/* }}} */
