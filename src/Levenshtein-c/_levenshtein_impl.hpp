#pragma once

#include <string.h>
#include <math.h>
/* for debugging */
#include <stdio.h>
#include <stdint.h>

#include <assert.h>
#include "_levenshtein.hpp"

#include <vector>
#include <memory>

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
size_t lev_set_median_index(size_t n, const size_t* lengths,
                     const CharT** strings,
                     const double* weights)
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
        d = (long int)rapidfuzz::string_metric::levenshtein(
          rapidfuzz::basic_string_view<CharT>(strings[j], lengths[j]),
          rapidfuzz::basic_string_view<CharT>(stri, leni)
        );
      }
      dist += weights[j] * (double)d;
      j++;
    }
    j++;  /* no need to compare item with itself */
    /* above diagonal */
    while (j < n && dist < mindist) {
      size_t dindex = (j - 1)*(j - 2)/2 + i;
      distances[dindex] = (long int)rapidfuzz::string_metric::levenshtein(
        rapidfuzz::basic_string_view<CharT>(strings[j], lengths[j]),
        rapidfuzz::basic_string_view<CharT>(stri, leni)
      );
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
CharT* lev_set_median(size_t n, const size_t* lengths,
               const CharT** strings,
               const double* weights,
               size_t* medlength)
{
  size_t minidx = lev_set_median_index(n, lengths, strings, weights);

  if (minidx == (size_t)-1)
    return NULL;

  *medlength = lengths[minidx];
  if (!lengths[minidx])
    return (CharT*)calloc(1, sizeof(CharT));

  CharT* result = (CharT*)safe_malloc(lengths[minidx], sizeof(CharT));
  if (!result)
    return NULL;
  return (CharT*)memcpy(result, strings[minidx], lengths[minidx]*sizeof(CharT));
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
template <typename CharT>
CharT* lev_opcodes_apply(size_t len1, const CharT* string1,
                  size_t len2, const CharT* string2,
                  size_t nb, const LevOpCode* bops,
                  size_t* len)
{
  /* this looks too complex for such a simple task, but note ops is not
   * a complete edit sequence, we have to be able to apply anything anywhere */
  CharT* dst = (CharT*)safe_malloc((len1 + len2), sizeof(CharT));
  CharT* dpos = dst;
  if (!dst) {
    *len = (size_t)(-1);
    return NULL;
  }
  const CharT* spos = string1;
  for (size_t i = nb; i; i--, bops++) {
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
template <typename CharT>
CharT* lev_editops_apply(size_t len1, const CharT *string1,
                  size_t len2, const CharT *string2,
                  size_t n, const LevEditOp *ops,
                  size_t *len)
{
  LEV_UNUSED(len2);

  /* this looks too complex for such a simple task, but note ops is not
   * a complete edit sequence, we have to be able to apply anything anywhere */
  CharT *dst = (CharT*)safe_malloc((n + len1), sizeof(CharT));
  if (!dst) {
    *len = (size_t)(-1);
    return NULL;
  }
  CharT *dpos = dst;
  const CharT *spos = string1;
  for (size_t i = n; i; i--, ops++) {
    /* XXX: this fine with gcc internal memcpy, but when memcpy is
     * actually a function, it may be pretty slow */
    size_t j = ops->spos - (size_t)(spos - string1) + (ops->type == LEV_EDIT_KEEP);
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
  size_t j = len1 - (size_t)(spos - string1);
  if (j) {
    memcpy(dpos, spos, j*sizeof(CharT));
    spos += j;
    dpos += j;
  }

  *len = (size_t)(dpos - dst);
  /* possible realloc failure is detected with *len != 0 */
  return (CharT*)realloc(dst, *len*sizeof(CharT));
}