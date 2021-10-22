#pragma once

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
