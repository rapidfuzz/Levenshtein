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

#include <limits>
#include <math.h>
#include <string.h>
/* for debugging */
#include <stdint.h>
#include <stdio.h>

#include "_levenshtein.hpp"

#define LEV_EPSILON 1e-14

/****************************************************************************
 *
 * Quick (voting) medians
 *
 ****************************************************************************/
/* {{{ */
struct HQItem {
    uint32_t c;
    double s;
    HQItem* n;
};

/* compute the sets of symbols each string contains, and the set of symbols
 * in any of them (symset).  meanwhile, count how many different symbols
 * there are (used below for symlist).
 * the symset is passed as an argument to avoid its allocation and
 * deallocation when it's used in the caller too */
class SymMap {
    std::unique_ptr<HQItem[]> symmap;

public:
    SymMap(const std::vector<RF_String>& strings)
    {
        symmap = std::make_unique<HQItem[]>(0x100);

        /* this is an ugly memory allocation avoiding hack: most hash elements
         * will probably contain none or one symbols only so, when p->n is equal
         * to symmap, it means there're no symbols yet, afters inserting the
         * first one, p->n becomes normally NULL and then it behaves like an
         * usual singly linked list */
        for (size_t i = 0; i < 0x100; i++)
            symmap[i].n = &symmap[0];
        for (size_t i = 0; i < strings.size(); i++) {
            visit(strings[i], [&](auto s1) {
                for (auto c : s1) {
                    int key = ((int)c + ((int)c >> 7)) & 0xff;
                    HQItem* p = symmap.get() + key;
                    if (p->n == symmap.get()) {
                        p->c = c;
                        p->n = NULL;
                        continue;
                    }
                    while (p->c != c && p->n != NULL)
                        p = p->n;
                    if (p->c != c) {
                        p->n = new HQItem;
                        p = p->n;
                        p->n = NULL;
                        p->c = c;
                    }
                }
            });
        }
    }

    ~SymMap()
    {
        for (size_t j = 0; j < 0x100; j++) {
            HQItem* p = &symmap[j];
            if (p->n == symmap.get() || p->n == NULL) continue;
            p = p->n;
            while (p) {
                HQItem* q = p;
                p = p->n;
                delete q;
            }
        }
    }

    void clear()
    {
        for (size_t i = 0; i < 0x100; i++) {
            HQItem* p = &symmap[i];
            if (p->n == symmap.get()) continue;
            while (p) {
                p->s = 0.0;
                p = p->n;
            }
        }
    }

    HQItem* get()
    {
        return symmap.get();
    }
};

std::basic_string<uint32_t> lev_quick_median(const std::vector<RF_String>& strings,
                                             const std::vector<double>& weights)
{
    std::basic_string<uint32_t> median; /* the resulting string */

    /* first check whether the result would be an empty string
     * and compute resulting string length */
    double ml = 0;
    double wl = 0;
    for (size_t i = 0; i < strings.size(); i++) {
        ml += weights[i] * strings[i].length;
        wl += weights[i];
    }

    if (wl == 0.0) return median;
    ml = floor(ml / wl + 0.499999);
    median.resize((size_t)ml);
    if (median.empty()) return median;

    /* find the symbol set;
     * now an empty symbol set is really a failure */
    SymMap symmap(strings);

    for (size_t j = 0; j < median.size(); j++) {
        /* clear the symbol probabilities */
        symmap.clear();

        /* let all strings vote */
        for (size_t i = 0; i < strings.size(); i++) {
            visit(strings[i], [&](auto s1) {
                double weighti = weights[i];
                size_t lengthi = (size_t)s1.size();
                double start = (double)lengthi / ml * (double)j;
                double end = start + (double)lengthi / ml;
                size_t istart = (size_t)floor(start);
                size_t iend = (size_t)ceil(end);

                /* rounding errors can overflow the buffer */
                if (iend > lengthi) iend = lengthi;

                /* the inner part, including the complete last character */
                for (size_t k = istart + 1; k < iend; k++) {
                    uint32_t c = static_cast<uint32_t>(s1[k]);
                    int key = (c + (c >> 7)) & 0xff;
                    HQItem* p = symmap.get() + key;
                    while (p->c != c)
                        p = p->n;
                    p->s += weighti;
                }
                /* the initial fraction */
                {
                    uint32_t c = static_cast<uint32_t>(s1[istart]);
                    int key = (c + (c >> 7)) & 0xff;
                    HQItem* p = symmap.get() + key;
                    while (p->c != c)
                        p = p->n;
                    p->s += weighti * ((double)(1 + istart) - start);
                }
                /* subtract what we counted from the last character but doesn't
                 * actually belong here.
                 * this strategy works also when istart+1 == iend (i.e., everything
                 * happens inside a one character) */
                {
                    uint32_t c = static_cast<uint32_t>(s1[iend - 1]);
                    int key = (c + (c >> 7)) & 0xff;
                    HQItem* p = symmap.get() + key;
                    while (p->c != c)
                        p = p->n;
                    p->s -= weighti * ((double)iend - end);
                }
            });
        }

        /* find the elected symbol */
        {
            HQItem* max = NULL;

            for (size_t i = 0; i < 0x100; i++) {
                HQItem* p = symmap.get() + i;
                if (p->n == symmap.get()) continue;
                while (p) {
                    if (!max || p->s > max->s) max = p;
                    p = p->n;
                }
            }
            median[j] = max->c;
        }
    }

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
 * Munkres-Blackman algorithm.
 */
std::vector<size_t> munkres_blackman(size_t n1, size_t n2, double* dists)
{
    size_t row = 0;

    /* allocate memory */
    /* 1 if column is covered */
    std::vector<size_t> covc(n1, 0);

    /* row of a z* in given column (1-base indices, so we can use zero as `none')*/
    std::vector<size_t> zstarc(n1, 0);

    /* 1 if row is covered */
    std::vector<size_t> covr(n2, 0);

    /* column of a z* in given row (1-base indices, so we can use zero as `none')*/
    std::vector<size_t> zstarr(n2, 0);

    /* column of a z' in given row (1-base indices, so we can use zero as `none')*/
    std::vector<size_t> zprimer(n2, 0);

    /* step 0 (subtract minimal distance) and step 1 (find zeroes) => [2] */
    auto step1 = [&]() {
        for (size_t j = 0; j < n1; j++) {
            size_t minidx = 0;
            double* col = dists + j;
            double min = *col;
            double* p = col + n1;
            for (size_t i = 1; i < n2; i++) {
                if (min > *p) {
                    minidx = i;
                    min = *p;
                }
                p += n1;
            }
            /* subtract */
            p = col;
            for (size_t i = 0; i < n2; i++) {
                *p -= min;
                if (*p < LEV_EPSILON) *p = 0.0;
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
                for (size_t i = 0; i < n2; i++) {
                    if (i != minidx && *p == 0.0 && !zstarc[j] && !zstarr[i]) {
                        zstarc[j] = i + 1;
                        zstarr[i] = j + 1;
                        break;
                    }
                    p += n1;
                }
            }
        }

        return 2;
    };

    /* step 2 (cover columns containing z*) => [0, 3] */
    auto step2 = [&]() {
        size_t nc = 0;
        for (size_t j = 0; j < n1; j++)
            if (zstarc[j]) {
                covc[j] = 1;
                nc++;
            }

        return (nc == n1) ? 0 : 3;
    };

    /* step 3 (find uncovered zeroes) => [3, 4, 5] */
    auto step3 = [&]() {
        /* search uncovered matrix entries */
        for (size_t j = 0; j < n1; j++) {
            double* p = dists + j;
            if (covc[j]) continue;

            for (size_t i = 0; i < n2; i++) {
                if (!covr[i] && *p == 0.0) {
                    /* when a zero is found, prime it */
                    zprimer[i] = j + 1;
                    if (zstarr[i]) {
                        /* if there's a z* in the same row,
                         * uncover the column, cover the row and redo */
                        covr[i] = 1;
                        covc[zstarr[i] - 1] = 0;
                        return 3;
                    }
                    /* if there's no z*,
                     * we are at the end of our path an can convert z'
                     * to z* */
                    row = i;
                    return 4;
                }
                p += n1;
            }
        }

        return 5;
    };

    /* step 4 (increment the number of z*)
     * i is the row number (we get it from step 3) => [2] */
    auto step4 = [&]() {
        row++;
        do {
            size_t x = row;

            row--;
            size_t j = zprimer[row] - 1; /* move to z' in the same row */
            zstarr[row] = j + 1;         /* mark it as z* in row buffer */
            row = zstarc[j];             /* move to z* in the same column */
            zstarc[j] = x;               /* mark the z' as being new z* */
        } while (row);

        std::fill(std::begin(zprimer), std::end(zprimer), 0);
        std::fill(std::begin(covr), std::end(covr), 0);
        std::fill(std::begin(covc), std::end(covc), 0);

        return 2;
    };

    /* step 5 (new zero manufacturer)
     * we can't get here, unless no zero is found at all => [3] */
    auto step5 = [&]() {
        /* find the smallest uncovered entry */
        double min = std::numeric_limits<double>::max();
        for (size_t j = 0; j < n1; j++) {
            double* p = dists + j;
            if (covc[j]) continue;
            for (size_t i = 0; i < n2; i++) {
                if (!covr[i] && min > *p) {
                    min = *p;
                }
                p += n1;
            }
        }
        /* add it to all covered rows */
        for (size_t i = 0; i < n2; i++) {
            double* p = dists + i * n1;
            if (!covr[i]) continue;
            for (size_t j = 0; j < n1; j++)
                *(p++) += min;
        }
        /* subtract if from all uncovered columns */
        for (size_t j = 0; j < n1; j++) {
            double* p = dists + j;
            if (covc[j]) continue;
            for (size_t i = 0; i < n2; i++) {
                *p -= min;
                if (*p < LEV_EPSILON) *p = 0.0;
                p += n1;
            }
        }

        return 3;
    };

    /* main */
    int next_step = 1;
    while (next_step) {
        switch (next_step) {
        case 1: next_step = step1(); break;
        case 2: next_step = step2(); break;
        case 3: next_step = step3(); break;
        case 4: next_step = step4(); break;
        case 5: next_step = step5(); break;
        default: next_step = 0; break;
        }
    }

    for (size_t j = 0; j < n1; j++)
        zstarc[j]--;
    return zstarc;
}

/* }}} */
