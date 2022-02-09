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
 * Set, sequence distances
 *
 ****************************************************************************/
/* {{{ */

/*
 * Munkers-Blackman algorithm.
 */
std::unique_ptr<size_t[]>
munkers_blackman(size_t n1, size_t n2, double *dists)
{
  size_t i, j;

  /* allocate memory */
  /* 1 if column is covered */
  auto covc = std::make_unique<size_t[]>(n1);
  std::fill(covc.get(), covc.get() + n1, 0);

  /* row of a z* in given column (1-base indices, so we can use zero as `none')*/
  auto zstarc = std::make_unique<size_t[]>(n1);
  std::fill(zstarc.get(), zstarc.get() + n1, 0);

  /* 1 if row is covered */
  auto covr = std::make_unique<size_t[]>(n2);
  std::fill(covr.get(), covr.get() + n2, 0);

  /* column of a z* in given row (1-base indices, so we can use zero as `none')*/
  auto zstarr = std::make_unique<size_t[]>(n2);
  std::fill(zstarr.get(), zstarr.get() + n2, 0);

  /* column of a z' in given row (1-base indices, so we can use zero as `none')*/
  auto zprimer = std::make_unique<size_t[]>(n2);
  std::fill(zprimer.get(), zprimer.get() + n2, 0);

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

    std::fill(zprimer.get(), zprimer.get() + n2, 0);
    std::fill(covr.get(), covr.get() + n2, 0);
    std::fill(covc.get(), covc.get() + n1, 0);
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
int lev_editops_check_errors(size_t len1, size_t len2,
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
void lev_editops_invert(size_t n, LevEditOp *ops)
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
LevMatchingBlock* lev_editops_matching_blocks(size_t len1, size_t len2,
                            size_t n, const LevEditOp *ops, size_t *nmblocks)
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
