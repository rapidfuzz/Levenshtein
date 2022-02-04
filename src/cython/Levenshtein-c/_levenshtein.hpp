/* @(#) $Id: Levenshtein.h,v 1.22 2005/01/13 20:02:56 yeti Exp $ */
#ifndef LEVENSHTEIN_H
#define LEVENSHTEIN_H

#include <numeric>
#include <memory>
#include <vector>
#include <rapidfuzz/distance/Indel.hpp>

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

lev_byte*
lev_greedy_median(size_t n,
                  const size_t *lengths,
                  const lev_byte *strings[],
                  const double *weights,
                  size_t *medlength);

lev_wchar*
lev_u_greedy_median(size_t n,
                    const size_t *lengths,
                    const lev_wchar *strings[],
                    const double *weights,
                    size_t *medlength);

lev_byte*
lev_median_improve(size_t len, const lev_byte *s,
                   size_t n, const size_t *lengths,
                   const lev_byte *strings[],
                   const double *weights,
                   size_t *medlength);

lev_wchar*
lev_u_median_improve(size_t len, const lev_wchar *s,
                     size_t n, const size_t *lengths,
                     const lev_wchar *strings[],
                     const double *weights,
                     size_t *medlength);

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

lev_byte*
lev_set_median(size_t n,
               const size_t *lengths,
               const lev_byte *strings[],
               const double *weights,
               size_t *medlength);

size_t
lev_set_median_index(size_t n, const size_t *lengths,
                     const lev_byte *strings[],
                     const double *weights);

lev_wchar*
lev_u_set_median(size_t n,
                 const size_t *lengths,
                 const lev_wchar *strings[],
                 const double *weights,
                 size_t *medlength);

size_t
lev_u_set_median_index(size_t n, const size_t *lengths,
                       const lev_wchar *strings[],
                       const double *weights);

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
