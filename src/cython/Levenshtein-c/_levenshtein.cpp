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
#include <rapidfuzz/distance/Indel.hpp>
#include <rapidfuzz/distance/Levenshtein.hpp>
#include <numeric>

#define LEV_EPSILON 1e-14
#define LEV_INFINITY 1e100

/* local functions */
static size_t* munkers_blackman(size_t n1, size_t n2, double *dists);

/****************************************************************************
 *
 * Basic stuff, Levenshtein distance
 *
 ****************************************************************************/
/* {{{ */

/**
 * @brief Wrapper for Levenshtein distance to handle exceptions
 */
template <typename CharT1, typename CharT2>
static size_t lev_levenshtein_distance(size_t len1, const CharT1* string1,
                                size_t len2, const CharT2* string2)
{
  try {
    return rapidfuzz::levenshtein_distance(string1, string1 + len1, string2, string2 + len2);
  } catch (...) {
    return (size_t)-1;
  }
}

/**
 * @brief Wrapper for Indel distance to handle exceptions
 */
template <typename CharT1, typename CharT2>
static size_t lev_indel_distance(size_t len1, const CharT1* string1,
                                size_t len2, const CharT2* string2)
{
  try {
    return rapidfuzz::indel_distance(string1, string1 + len1, string2, string2 + len2);
  } catch (...) {
    return (size_t)-1;
  }
}


/* }}} */

/****************************************************************************
 *
 * Generalized medians, the greedy algorithm, and greedy improvements
 *
 ****************************************************************************/
/* {{{ */

/* compute the sets of symbols each string contains, and the set of symbols
 * in any of them (symset).  meanwhile, count how many different symbols
 * there are (used below for symlist). */
static lev_byte*
make_symlist(size_t n, const size_t *lengths,
             const lev_byte *strings[], size_t *symlistlen)
{
  short int *symset;  /* indexed by ALL symbols, contains 1 for symbols
                         present in the strings, zero for others */
  size_t i, j;
  lev_byte *symlist;

  symset = (short int*)calloc(0x100, sizeof(short int));
  if (!symset) {
    *symlistlen = (size_t)(-1);
    return NULL;
  }
  *symlistlen = 0;
  for (i = 0; i < n; i++) {
    const lev_byte *stri = strings[i];
    for (j = 0; j < lengths[i]; j++) {
      int c = stri[j];
      if (!symset[c]) {
        (*symlistlen)++;
        symset[c] = 1;
      }
    }
  }
  if (!*symlistlen) {
    free(symset);
    return NULL;
  }

  /* create dense symbol table, so we can easily iterate over only characters
   * present in the strings */
  {
    size_t pos = 0;
    symlist = (lev_byte*)safe_malloc((*symlistlen), sizeof(lev_byte));
    if (!symlist) {
      *symlistlen = (size_t)(-1);
      free(symset);
      return NULL;
    }
    for (j = 0; j < 0x100; j++) {
      if (symset[j])
        symlist[pos++] = (lev_byte)j;
    }
  }
  free(symset);

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
lev_byte* lev_greedy_median(size_t n, const size_t *lengths,
                            const lev_byte *strings[],
                            const double *weights,
                            size_t *medlength)
{
  size_t i;  /* usually iterates over strings (n) */
  size_t j;  /* usually iterates over characters */
  size_t len;  /* usually iterates over the approximate median string */
  size_t symlistlen;  /* length of symlist */
  size_t maxlen;  /* maximum input string length */
  size_t stoplen;  /* maximum tried median string length -- this is slightly
                      higher than maxlen, because the median string may be
                      longer than any of the input strings */
  size_t **rows;  /* Levenshtein matrix rows for each string, we need to keep
                     only one previous row to construct the current one */
  size_t *row;  /* a scratch buffer for new Levenshtein matrix row computation,
                   shared among all strings */
  lev_byte *median;  /* the resulting approximate median string */
  double *mediandist;  /* the total distance of the best median string of
                          given length.  warning!  mediandist[0] is total
                          distance for empty string, while median[] itself
                          is normally zero-based */

  /* find all symbols */
  /** list of symbols present in the strings,
   * we iterate over it insead of set of all
   * existing symbols
   */
  lev_byte* symlist = make_symlist(n, lengths, strings, &symlistlen);
  if (!symlist) {
    *medlength = 0;
    if (symlistlen != 0)
      return NULL;
    else
      return (lev_byte*)calloc(1, sizeof(lev_byte));
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
  median = (lev_byte*)safe_malloc(stoplen, sizeof(lev_byte));
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
    mediandist[0] += (double)lengths[i]*weights[i];

  /* build up the approximate median string symbol by symbol
   * XXX: we actually exit on break below, but on the same condition */
  for (len = 1; len <= stoplen; len++) {
    lev_byte symbol;
    double minminsum = LEV_INFINITY;
    row[0] = len;
    /* iterate over all symbols we may want to add */
    for (j = 0; j < symlistlen; j++) {
      double totaldist = 0.0;
      double minsum = 0.0;
      symbol = symlist[j];
      /* sum Levenshtein distances from all the strings, with given weights */
      for (i = 0; i < n; i++) {
        const lev_byte *stri = strings[i];
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
    for (i = 0; i < n; i++) {
      const lev_byte *stri = strings[i];
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
  size_t bestlen = std::distance(mediandist, std::max_element(mediandist, mediandist + stoplen));

  /* clean up */
  for (i = 0; i < n; i++)
    free(rows[i]);
  free(rows);
  free(row);
  free(symlist);
  free(mediandist);

  /* return result */
  {
    lev_byte *result = (lev_byte*)safe_malloc(bestlen, sizeof(lev_byte));
    if (!result) {
      free(median);
      return NULL;
    }
    memcpy(result, median, bestlen*sizeof(lev_byte));
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
finish_distance_computations(size_t len1, lev_byte *string1,
                             size_t n, const size_t *lengths,
                             const lev_byte **strings,
                             const double *weights, size_t **rows,
                             size_t *row)
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
    size_t *rowi = rows[j];  /* current row */
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
    offset = rowi[0];
    if (leni == 0) {
      distsum += (double)(offset + len)*weights[j];
      continue;
    }

    /* complete the matrix */
    memcpy(row, rowi, (leni + 1)*sizeof(size_t));
    end = row + leni;

    for (i = 1; i <= len; i++) {
      size_t *p = row + 1;
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
  lev_byte *symlist;  /* list of symbols present in the strings,
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
  lev_byte *median;  /* the resulting approximate median string */
  size_t medlen;  /* the current approximate median string length */
  double minminsum;  /* the current total distance sum */

  /* find all symbols */
  symlist = make_symlist(n, lengths, strings, &symlistlen);
  if (!symlist) {
    *medlength = 0;
    if (symlistlen != 0)
      return NULL;
    else
      return (lev_byte*)calloc(1, sizeof(lev_byte));
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
  median = (lev_byte*)safe_malloc((stoplen+1), sizeof(lev_byte));
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
  memcpy(median, s, (medlen)*sizeof(lev_byte));
  minminsum = finish_distance_computations(medlen, median,
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
      for (j = 0; j < symlistlen; j++) {
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
    for (j = 0; j < symlistlen; j++) {
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
              (medlen - pos)*sizeof(lev_byte));
      median[pos] = symbol;
      medlen++;
      break;

      case LEV_EDIT_DELETE:
      memmove(median+pos, median + pos+1,
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
    lev_byte *result = (lev_byte*)safe_malloc(medlen, sizeof(lev_byte));
    if (!result) {
      median--;
      free(median);
      return NULL;
    }
    *medlength = medlen;
    memcpy(result, median, medlen*sizeof(lev_byte));
    median--;
    free(median);
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
  lev_wchar *symlist;

  *symlistlen = 0;
  if (std::all_of(lengths, lengths + n, [](size_t x){ return x == 0; }))
  {
    return NULL;
  }

  /* find all symbols, use a kind of hash for storage */
  HItem* symmap = (HItem*)safe_malloc(0x100, sizeof(HItem));
  if (!symmap) {
    *symlistlen = (size_t)(-1);
    return NULL;
  }
  /* this is an ugly memory allocation avoiding hack: most hash elements
   * will probably contain none or one symbols only so, when p->n is equal
   * to symmap, it means there're no symbols yet, afters inserting the
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
  {
    size_t pos = 0;
    symlist = (lev_wchar*)safe_malloc((*symlistlen), sizeof(lev_wchar));
    if (!symlist) {
      free_usymlist_hash(symmap);
      *symlistlen = (size_t)(-1);
      return NULL;
    }
    for (size_t j = 0; j < 0x100; j++) {
      HItem *p = symmap + j;
      while (p != NULL && p->n != symmap) {
        symlist[pos++] = p->c;
        p = p->n;
      }
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

  mediandist[0] = std::inner_product(lengths, lengths + n, weights, 0.0);

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
  size_t bestlen = std::distance(mediandist, std::max_element(mediandist, mediandist + stoplen));

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
  size_t i, j;
  size_t offset;  /* row[0]; offset + len1 give together real len of string1 */
  double distsum = 0.0;  /* sum of distances */

  /* catch trivia case */
  if (len1 == 0) {
    for (j = 0; j < n; j++)
      distsum += (double)rows[j][lengths[j]] * weights[j];
    return distsum;
  }

  /* iterate through the strings and sum the distances */
  for (j = 0; j < n; j++) {
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

    for (i = 1; i <= len; i++) {
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
    lev_wchar orig_symbol;
    double sum;

    lev_wchar symbol = median[pos];
    LevEditType operation = LEV_EDIT_KEEP;
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
static lev_byte*
make_symlistset(size_t n, const size_t *lengths,
                const lev_byte *strings[], size_t *symlistlen,
                double *symset)
{
  size_t i, j;
  lev_byte *symlist;

  if (!symset) {
    *symlistlen = (size_t)(-1);
    return NULL;
  }
  memset(symset, 0, 0x100*sizeof(double));
  *symlistlen = 0;
  for (i = 0; i < n; i++) {
    const lev_byte *stri = strings[i];
    for (j = 0; j < lengths[i]; j++) {
      int c = stri[j];
      if (!symset[c]) {
        (*symlistlen)++;
        symset[c] = 1.0;
      }
    }
  }
  if (!*symlistlen)
    return NULL;

  /* create dense symbol table, so we can easily iterate over only characters
   * present in the strings */
  {
    size_t pos = 0;
    symlist = (lev_byte*)safe_malloc((*symlistlen), sizeof(lev_byte));
    if (!symlist) {
      *symlistlen = (size_t)(-1);
      return NULL;
    }
    for (j = 0; j < 0x100; j++) {
      if (symset[j])
        symlist[pos++] = (lev_byte)j;
    }
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
  size_t symlistlen, len, i, j, k;
  lev_byte *symlist;
  lev_byte *median;  /* the resulting string */
  double *symset;

  /* first check whether the result would be an empty string 
   * and compute resulting string length */
  double ml = std::inner_product(weights, weights + n, lengths, 0.0);
  double wl = std::accumulate(   weights, weights + n, 0.0);

  if (wl == 0.0)
    return (lev_byte*)calloc(1, sizeof(lev_byte));
  ml = floor(ml/wl + 0.499999);
  *medlength = len = (size_t)ml;
  if (!len)
    return (lev_byte*)calloc(1, sizeof(lev_byte));
  median = (lev_byte*)safe_malloc(len, sizeof(lev_byte));
  if (!median)
    return NULL;

  /* find the symbol set;
   * now an empty symbol set is really a failure */
  symset = (double*)calloc(0x100, sizeof(double));
  if (!symset) {
    free(median);
    return NULL;
  }
  symlist = make_symlistset(n, lengths, strings, &symlistlen, symset);
  if (!symlist) {
    free(median);
    free(symset);
    return NULL;
  }

  for (j = 0; j < len; j++) {
    /* clear the symbol probabilities */
    if (symlistlen < 32) {
      for (i = 0; i < symlistlen; i++)
        symset[symlist[i]] = 0.0;
    }
    else
      memset(symset, 0, 0x100*sizeof(double));

    /* let all strings vote */
    for (i = 0; i < n; i++) {
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

      for (k = istart+1; k < iend; k++)
        symset[stri[k]] += weighti;
      symset[stri[istart]] += weighti * ((double)(1 + istart) - start);
      symset[stri[iend-1]] -= weighti * ((double)iend - end);
    }

    /* find the elected symbol */
    k = symlist[0];
    for (i = 1; i < symlistlen; i++) {
      if (symset[symlist[i]] > symset[k])
        k = symlist[i];
    }
    median[j] = (lev_byte)k;
  }

  free(symset);
  free(symlist);

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
  lev_wchar *symlist;
  *symlistlen = 0;
  if (std::all_of(lengths, lengths + n, [](size_t x){ return x == 0; }))
    return NULL;

  /* this is an ugly memory allocation avoiding hack: most hash elements
   * will probably contain none or one symbols only so, when p->n is equal
   * to symmap, it means there're no symbols yet, afters inserting the
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
  {
    size_t pos = 0;
    symlist = (lev_wchar*)safe_malloc((*symlistlen), sizeof(lev_wchar));
    if (!symlist) {
      *symlistlen = (size_t)(-1);
      return NULL;
    }
    for (size_t j = 0; j < 0x100; j++) {
      HQItem *p = symmap + j;
      while (p != NULL && p->n != symmap) {
        symlist[pos++] = p->c;
        p = p->n;
      }
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
  size_t symlistlen, len, i, j, k;
  lev_wchar *symlist;
  lev_wchar *median;  /* the resulting string */
  HQItem *symmap;

  /* first check whether the result would be an empty string 
   * and compute resulting string length */
  double ml = std::inner_product(weights, weights + n, lengths, 0.0);
  double wl = std::accumulate(   weights, weights + n, 0.0);

  if (wl == 0.0)
    return (lev_wchar*)calloc(1, sizeof(lev_wchar));
  ml = floor(ml/wl + 0.499999);
  *medlength = len = (size_t)ml;
  if (!len)
    return (lev_wchar*)calloc(1, sizeof(lev_wchar));
  median = (lev_wchar*)safe_malloc(len, sizeof(lev_wchar));
  if (!median)
    return NULL;

  /* find the symbol set;
   * now an empty symbol set is really a failure */
  symmap = (HQItem*)safe_malloc(0x100, sizeof(HQItem));
  if (!symmap) {
    free(median);
    return NULL;
  }
  symlist = make_usymlistset(n, lengths, strings, &symlistlen, symmap);
  if (!symlist) {
    free(median);
    free_usymlistset_hash(symmap);
    return NULL;
  }

  for (j = 0; j < len; j++) {
    /* clear the symbol probabilities */
    for (i = 0; i < 0x100; i++) {
      HQItem *p = symmap + i;
      if (p->n == symmap)
        continue;
      while (p) {
        p->s = 0.0;
        p = p->n;
      }
    }

    /* let all strings vote */
    for (i = 0; i < n; i++) {
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
      for (k = istart+1; k < iend; k++) {
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

      for (i = 0; i < 0x100; i++) {
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
 * Set medians
 *
 ****************************************************************************/
/* {{{ */

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
size_t
lev_set_median_index(size_t n, const size_t *lengths,
                     const lev_byte *strings[],
                     const double *weights)
{
  size_t minidx = 0;
  double mindist = LEV_INFINITY;
  size_t i;
  long int *distances;

  distances = (long int*)safe_malloc((n*(n - 1)/2), sizeof(long int));
  if (!distances)
    return (size_t)-1;

  memset(distances, 0xff, (n*(n - 1)/2)*sizeof(long int)); /* XXX */
  for (i = 0; i < n; i++) {
    size_t j = 0;
    double dist = 0.0;
    const lev_byte *stri = strings[i];
    size_t leni = lengths[i];
    /* below diagonal */
    while (j < i && dist < mindist) {
      size_t dindex = (i - 1)*(i - 2)/2 + j;
      long int d;
      if (distances[dindex] >= 0)
        d = distances[dindex];
      else {
        d = (long int)lev_levenshtein_distance(lengths[j], strings[j], leni, stri);
        if (d < 0) {
          free(distances);
          return (size_t)-1;
        }
      }
      dist += weights[j] * (double)d;
      j++;
    }
    j++;  /* no need to compare item with itself */
    /* above diagonal */
    while (j < n && dist < mindist) {
      size_t dindex = (j - 1)*(j - 2)/2 + i;
      distances[dindex] = (long int)lev_levenshtein_distance(lengths[j], strings[j], leni, stri);
      if (distances[dindex] < 0) {
        free(distances);
        return (size_t)-1;
      }
      dist += weights[j] * (double)distances[dindex];
      j++;
    }

    if (dist < mindist) {
      mindist = dist;
      minidx = i;
    }
  }

  free(distances);
  return minidx;
}

/**
 * lev_u_set_median_index:
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
size_t
lev_u_set_median_index(size_t n, const size_t *lengths,
                       const lev_wchar *strings[],
                       const double *weights)
{
  size_t minidx = 0;
  double mindist = LEV_INFINITY;
  size_t i;
  long int *distances;

  distances = (long int*)safe_malloc((n*(n - 1)/2), sizeof(long int));
  if (!distances)
    return (size_t)-1;

  memset(distances, 0xff, (n*(n - 1)/2)*sizeof(long int)); /* XXX */
  for (i = 0; i < n; i++) {
    size_t j = 0;
    double dist = 0.0;
    const lev_wchar *stri = strings[i];
    size_t leni = lengths[i];
    /* below diagonal */
    while (j < i && dist < mindist) {
      size_t dindex = (i - 1)*(i - 2)/2 + j;
      long int d;
      if (distances[dindex] >= 0)
        d = distances[dindex];
      else {
        d = (long int)lev_levenshtein_distance(lengths[j], strings[j], leni, stri);
        if (d < 0) {
          free(distances);
          return (size_t)-1;
        }
      }
      dist += weights[j] * (double)d;
      j++;
    }
    j++;  /* no need to compare item with itself */
    /* above diagonal */
    while (j < n && dist < mindist) {
      size_t dindex = (j - 1)*(j - 2)/2 + i;
      distances[dindex] = (long int)lev_levenshtein_distance(lengths[j], strings[j], leni, stri);
      if (distances[dindex] < 0) {
        free(distances);
        return (size_t)-1;
      }
      dist += weights[j] * (double)distances[dindex];
      j++;
    }

    if (dist < mindist) {
      mindist = dist;
      minidx = i;
    }
  }

  free(distances);
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
lev_byte*
lev_set_median(size_t n, const size_t *lengths,
               const lev_byte *strings[],
               const double *weights,
               size_t *medlength)
{
  size_t minidx = lev_set_median_index(n, lengths, strings, weights);
  lev_byte *result;

  if (minidx == (size_t)-1)
    return NULL;

  *medlength = lengths[minidx];
  if (!lengths[minidx])
    return (lev_byte*)calloc(1, sizeof(lev_byte));

  result = (lev_byte*)safe_malloc(lengths[minidx], sizeof(lev_byte));
  if (!result)
    return NULL;
  return (lev_byte*)memcpy(result, strings[minidx], lengths[minidx]*sizeof(lev_byte));
}

/**
 * lev_u_set_median:
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
lev_wchar*
lev_u_set_median(size_t n, const size_t *lengths,
                 const lev_wchar *strings[],
                 const double *weights,
                 size_t *medlength)
{
  size_t minidx = lev_u_set_median_index(n, lengths, strings, weights);
  lev_wchar *result;

  if (minidx == (size_t)-1)
    return NULL;

  *medlength = lengths[minidx];
  if (!lengths[minidx])
    return (lev_wchar*)calloc(1, sizeof(lev_wchar));

  result = (lev_wchar*)safe_malloc(lengths[minidx], sizeof(lev_wchar));
  if (!result)
    return NULL;
  return (lev_wchar*)memcpy(result, strings[minidx], lengths[minidx]*sizeof(lev_wchar));
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
  size_t i;
  double* row;  /* we only need to keep one row of costs */
  double* last;
  double* end;

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

  /* make the inner cycle (i.e. strings2) the longer one */
  if (n1 > n2) {
    std::swap(n1, n2);
    std::swap(lengths1, lengths2);
    std::swap(strings1, strings2);
  }
  n1++;
  n2++;

  /* initalize first row */
  row = (double*)safe_malloc(n2, sizeof(double));
  if (!row)
    return -1.0;
  last = row + n2 - 1;
  end = row + n2;
  std::iota(row, end, 0.0);

  /* go through the matrix and compute the costs.  yes, this is an extremely
   * obfuscated version, but also extremely memory-conservative and relatively
   * fast.  */
  for (i = 1; i < n1; i++) {
    double *p = row + 1;
    const lev_byte *str1 = strings1[i - 1];
    const size_t len1 = lengths1[i - 1];
    const lev_byte **str2p = strings2;
    const size_t *len2p = lengths2;
    double D = (double)i - 1.0;
    double x = (double)i;
    while (p != end) {
      size_t l = len1 + *len2p;
      double q;
      if (l == 0)
        q = D;
      else {
        size_t d = lev_indel_distance(len1, str1, *(len2p++), *(str2p++));
        if (d == (size_t)(-1)) {
          free(row);
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

  {
    double q = *last;
    free(row);
    return q;
  }
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
  size_t i;
  double* row;  /* we only need to keep one row of costs */
  double* last;
  double* end;

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

  /* make the inner cycle (i.e. strings2) the longer one */
  if (n1 > n2) {
    std::swap(n1, n2);
    std::swap(lengths1, lengths2);
    std::swap(strings1, strings2);
  }
  n1++;
  n2++;

  /* initalize first row */
  row = (double*)safe_malloc(n2, sizeof(double));
  if (!row)
    return -1.0;
  last = row + n2 - 1;
  end = row + n2;
  std::iota(row, end, 0.0);

  /* go through the matrix and compute the costs.  yes, this is an extremely
   * obfuscated version, but also extremely memory-conservative and relatively
   * fast.  */
  for (i = 1; i < n1; i++) {
    double *p = row + 1;
    const lev_wchar *str1 = strings1[i - 1];
    const size_t len1 = lengths1[i - 1];
    const lev_wchar **str2p = strings2;
    const size_t *len2p = lengths2;
    double D = (double)i - 1.0;
    double x = (double)i;
    while (p != end) {
      size_t l = len1 + *len2p;
      double q;
      if (l == 0)
        q = D;
      else {
        size_t d = lev_indel_distance(len1, str1, *(len2p++), *(str2p++));
        if (d == (size_t)(-1)) {
          free(row);
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

  {
    double q = *last;
    free(row);
    return q;
  }
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
  double *dists;  /* the (modified) distance matrix, indexed [row*n1 + col] */
  double *r;
  size_t i, j;
  size_t *map;
  double sum;

  /* catch trivial cases */
  if (n1 == 0)
    return (double)n2;
  if (n2 == 0)
    return (double)n1;

  /* make the number of columns (n1) smaller than the number of rows */
  if (n1 > n2) {
    size_t nx = n1;
    const size_t *lx = lengths1;
    const lev_byte **sx = strings1;
    n1 = n2;
    n2 = nx;
    lengths1 = lengths2;
    lengths2 = lx;
    strings1 = strings2;
    strings2 = sx;
  }

  /* compute distances from each to each */
  r = dists = (double*)safe_malloc_3(n1, n2, sizeof(double));
  if (!r)
    return -1.0;
  for (i = 0; i < n2; i++) {
    size_t len2 = lengths2[i];
    const lev_byte *str2 = strings2[i];
    const size_t *len1p = lengths1;
    const lev_byte **str1p = strings1;
    for (j = 0; j < n1; j++) {
      size_t l = len2 + *len1p;
      if (l == 0)
        *(r++) = 0.0;
      else {
        size_t d = lev_indel_distance(len2, str2, *(len1p++), *(str1p)++);
        if (d == (size_t)(-1)) {
          free(r);
          return -1.0;
        }
        *(r++) = (double)d / (double)l;
      }
    }
  }

  /* find the optimal mapping between the two sets */
  map = munkers_blackman(n1, n2, dists);
  if (!map)
    return -1.0;

  /* sum the set distance */
  sum = (double)(n2 - n1);
  for (j = 0; j < n1; j++) {
    size_t l;
    i = map[j];
    l = lengths1[j] + lengths2[i];
    if (l > 0) {
      size_t d = lev_indel_distance(lengths1[j], strings1[j], lengths2[i], strings2[i]);
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
  double *dists;  /* the (modified) distance matrix, indexed [row*n1 + col] */
  double *r;
  size_t i, j;
  size_t *map;
  double sum;

  /* catch trivial cases */
  if (n1 == 0)
    return (double)n2;
  if (n2 == 0)
    return (double)n1;

  /* make the number of columns (n1) smaller than the number of rows */
  if (n1 > n2) {
    size_t nx = n1;
    const size_t *lx = lengths1;
    const lev_wchar **sx = strings1;
    n1 = n2;
    n2 = nx;
    lengths1 = lengths2;
    lengths2 = lx;
    strings1 = strings2;
    strings2 = sx;
  }

  /* compute distances from each to each */
  r = dists = (double*)safe_malloc_3(n1, n2, sizeof(double));
  if (!r)
    return -1.0;
  for (i = 0; i < n2; i++) {
    size_t len2 = lengths2[i];
    const lev_wchar *str2 = strings2[i];
    const size_t *len1p = lengths1;
    const lev_wchar **str1p = strings1;
    for (j = 0; j < n1; j++) {
      size_t l = len2 + *len1p;
      if (l == 0)
        *(r++) = 0.0;
      else {
        size_t d = lev_indel_distance(len2, str2, *(len1p++), *(str1p)++);
        if (d == (size_t)(-1)) {
          free(r);
          return -1.0;
        }
        *(r++) = (double)d / (double)l;
      }
    }
  }

  /* find the optimal mapping between the two sets */
  map = munkers_blackman(n1, n2, dists);
  if (!map)
    return -1.0;

  /* sum the set distance */
  sum = (double)(n2 - n1);
  for (j = 0; j < n1; j++) {
    size_t l;
    i = map[j];
    l = lengths1[j] + lengths2[i];
    if (l > 0) {
      size_t d = lev_indel_distance(lengths1[j], strings1[j], lengths2[i], strings2[i]);
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
  size_t *covc, *covr;  /* 1 if column/row is covered */
  /* these contain 1-based indices, so we can use zero as `none'
   * zstarr: column of a z* in given row
   * zstarc: row of a z* in given column
   * zprimer: column of a z' in given row */
  size_t *zstarr, *zstarc, *zprimer;

  /* allocate memory */
  covc = (size_t*)calloc(n1, sizeof(size_t));
  if (!covc)
    return NULL;
  zstarc = (size_t*)calloc(n1, sizeof(size_t));
  if (!zstarc) {
    free(covc);
    return NULL;
  }
  covr = (size_t*)calloc(n2, sizeof(size_t));
  if (!covr) {
    free(zstarc);
    free(covc);
    return NULL;
  }
  zstarr = (size_t*)calloc(n2, sizeof(size_t));
  if (!zstarr) {
    free(covr);
    free(zstarc);
    free(covc);
    return NULL;
  }
  zprimer = (size_t*)calloc(n2, sizeof(size_t));
  if (!zprimer) {
    free(zstarr);
    free(covr);
    free(zstarc);
    free(covc);
    return NULL;
  }

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
    memset(zprimer, 0, n2*sizeof(size_t));
    memset(covr, 0, n2*sizeof(size_t));
    memset(covc, 0, n1*sizeof(size_t));
  }

  free(dists);
  free(covc);
  free(covr);
  /*free(zstarc);  this is the result */
  free(zstarr);
  free(zprimer);

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
 * lev_editops_check_errors:
 * @len1: The length of an eventual @ops source string.
 * @len2: The length of an eventual @ops destination string.
 * @n: The length of @ops.
 * @ops: An array of elementary edit operations.
 *
 * Checks whether @ops is consistent and applicable as a partial edit from a
 * string of length @len1 to a string of length @len2.
 *
 * Returns: Zero if @ops seems OK, a nonzero error code otherwise.
 **/
int
lev_editops_check_errors(size_t len1, size_t len2,
                         size_t n, const LevEditOp *ops)
{
  if (!n)
    return LEV_EDIT_ERR_OK;

  /* check bounds */
  const LevEditOp* o = ops;
  for (size_t i = n; i; i--, o++) {
    if (o->type >= LEV_EDIT_LAST)
      return LEV_EDIT_ERR_TYPE;
    if (o->spos > len1 || o->dpos > len2)
      return LEV_EDIT_ERR_OUT;
    if (o->spos == len1 && o->type != LEV_EDIT_INSERT)
      return LEV_EDIT_ERR_OUT;
    if (o->dpos == len2 && o->type != LEV_EDIT_DELETE)
      return LEV_EDIT_ERR_OUT;
  }

  /* check ordering */
  o = ops + 1;
  for (size_t i = n - 1; i; i--, o++, ops++) {
    if (o->spos < ops->spos || o->dpos < ops->dpos)
      return LEV_EDIT_ERR_ORDER;
  }

  return LEV_EDIT_ERR_OK;
}

/**
 * lev_opcodes_check_errors:
 * @len1: The length of an eventual @ops source string.
 * @len2: The length of an eventual @ops destination string.
 * @nb: The length of @bops.
 * @bops: An array of difflib block edit operation codes.
 *
 * Checks whether @bops is consistent and applicable as an edit from a
 * string of length @len1 to a string of length @len2.
 *
 * Returns: Zero if @bops seems OK, a nonzero error code otherwise.
 **/
int lev_opcodes_check_errors(size_t len1, size_t len2,
                             size_t nb, const LevOpCode *bops)
{
  if (!nb)
    return 1;

  /* check completenes */
  if (bops->sbeg || bops->dbeg
      || bops[nb - 1].send != len1 || bops[nb - 1].dend != len2)
    return LEV_EDIT_ERR_SPAN;

  /* check bounds and block consistency */
  const LevOpCode* b = bops;
  for (size_t i = nb; i; i--, b++) {
    if (b->send > len1 || b->dend > len2)
      return LEV_EDIT_ERR_OUT;
    switch (b->type) {
      case LEV_EDIT_KEEP:
      case LEV_EDIT_REPLACE:
      if (b->dend - b->dbeg != b->send - b->sbeg || b->dend == b->dbeg)
        return LEV_EDIT_ERR_BLOCK;
      break;

      case LEV_EDIT_INSERT:
      if (b->dend - b->dbeg == 0 || b->send - b->sbeg != 0)
        return LEV_EDIT_ERR_BLOCK;
      break;

      case LEV_EDIT_DELETE:
      if (b->send - b->sbeg == 0 || b->dend - b->dbeg != 0)
        return LEV_EDIT_ERR_BLOCK;
      break;

      default:
      return LEV_EDIT_ERR_TYPE;
      break;
    }
  }

  /* check ordering */
  b = bops + 1;
  for (size_t i = nb - 1; i; i--, b++, bops++) {
    if (b->sbeg != bops->send || b->dbeg != bops->dend)
      return LEV_EDIT_ERR_ORDER;
  }

  return LEV_EDIT_ERR_OK;
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
 * lev_opcodes_invert:
 * @nb: The length of @ops.
 * @bops: An array of difflib block edit operation codes.
 *
 * Inverts the sense of @ops.  It is modified in place.
 *
 * In other words, @ops becomes a partial edit for the original source
 * and destination strings with their roles exchanged.
 **/
void lev_opcodes_invert(size_t nb, LevOpCode *bops)
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
  size_t nmb, i, spos, dpos;
  LevEditType type;
  const LevEditOp *o;
  LevMatchingBlock *mblocks, *mb;

  /* compute the number of matching blocks */
  nmb = 0;
  o = ops;
  spos = dpos = 0;
  type = LEV_EDIT_KEEP;
  for (i = n; i; ) {
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
  mb = mblocks = (LevMatchingBlock*)safe_malloc(nmb, sizeof(LevMatchingBlock));
  if (!mblocks) {
    *nmblocks = (size_t)(-1);
    return NULL;
  }
  o = ops;
  spos = dpos = 0;
  type = LEV_EDIT_KEEP;
  for (i = n; i; ) {
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
LevMatchingBlock* lev_opcodes_matching_blocks(size_t len1, size_t,
                            size_t nb,
                            const LevOpCode *bops,
                            size_t *nmblocks)
{
  size_t nmb, i;
  const LevOpCode *b;
  LevMatchingBlock *mblocks, *mb;

  /* compute the number of matching blocks */
  nmb = 0;
  b = bops;
  for (i = nb; i; i--, b++) {
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
  mb = mblocks = (LevMatchingBlock*)safe_malloc(nmb, sizeof(LevMatchingBlock));
  if (!mblocks) {
    *nmblocks = (size_t)(-1);
    return NULL;
  }
  b = bops;
  for (i = nb; i; i--, b++) {
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
                     const LevEditOp* ops,
                     size_t ns,
                     const LevEditOp* sub,
                     size_t *nrem)
{
    static const int shifts[] = { 0, 0, 1, -1 };
    LevEditOp *rem;
    size_t i, j;
    int shift;

    /* compute remainder size */
    *nrem = (size_t)-1;

    size_t nr = std::accumulate(ops, ops + n, 0, [](size_t a, const LevEditOp& op) {
      return a + (op.type != LEV_EDIT_KEEP);
    });
  
    size_t nn = std::accumulate(sub, sub + ns, 0, [](size_t a, const LevEditOp& op) {
      return a + (op.type != LEV_EDIT_KEEP);
    });

    if (nn > nr)
        return NULL;
    nr -= nn;

    /* subtract */
    /* we could simply return NULL when nr == 0, but then it would be possible
     * to subtract *any* sequence of the right length to get an empty sequence
     * -- clrealy incorrectly; so we have to scan the list to check */
    rem = nr ? (LevEditOp*)safe_malloc(nr, sizeof(LevEditOp)) : NULL;
    j = nn = 0;
    shift = 0;
    for (i = 0; i < ns; i++) {
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
