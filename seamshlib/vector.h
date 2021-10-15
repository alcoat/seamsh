/*
 * MigFlow - Copyright (C) <2010-2020>
 * <Universite catholique de Louvain (UCL), Belgium
 *  Universite de Montpellier, France>
 * 	
 * List of the contributors to the development of MigFlow: see AUTHORS file.
 * Description and complete License: see LICENSE file.
 * 	
 * This program (MigFlow) is free software: 
 * you can redistribute it and/or modify it under the terms of the GNU Lesser General 
 * Public License as published by the Free Software Foundation, either version
 * 3 of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program (see COPYING and COPYING.LESSER files).  If not, 
 * see <http://www.gnu.org/licenses/>.
 */

#ifndef _VECTOR_H_
#define _VECTOR_H_
#include <string.h>
#include <stdlib.h>

static size_t _vectorSize(void *m) {
  return m == NULL ? 0 : (*((size_t*)m - 1));
}

static void _vectorFree(void **m) {
  if (*m != NULL)
    free(((size_t *)*m) - 2);
  *m = NULL;
}
static void *_vectorPush(void **m, size_t s) {
  if (s == 0) return *m;
  if (*m == NULL) {
    size_t *n = (size_t*)malloc(s * 2 + 2 * sizeof(size_t));
    n[0] = 2 * s;
    n[1] = s;
    *m = n + 2;
    return *m;
  }
  size_t *n = (*(size_t**)m) - 2;
  n[1] += s;
  if (n[0] < n[1]) {
    n[0] *= 2;
    n = (size_t*)realloc(n, n[0] + 2 * sizeof(size_t));
    *m= n + 2;
  }
  return ((char*) *m) + n[1] - s;
}

static void *_vectorInsert(void **m, size_t p, size_t s) {
  _vectorPush(m, s);
  memmove(((char*)*m) + p + s, ((char*)*m) + p, _vectorSize(*m) - s - p);
  return ((char*)*m) + p;
}

static void *_vectorDup(void *m) {
  if (m == NULL)
    return NULL;
  size_t *n = ((size_t*)m - 2);
  size_t N = n[1];
  size_t *a = (size_t*) malloc(sizeof(size_t) * 2 + N);
  memcpy(a, n, sizeof(size_t) * 2 + N);
  a[0] = a[1] = N;
  return a + 2;
}
static void vectorClear(void *m) {
  if (m != NULL)
    *(((size_t*) m) - 1) = 0;
}
static void _vectorPop(void *m, size_t s) {
  if (m != NULL) {
    *((size_t*)m - 1)-= s;
  }
}
static void _vectorRemoveFlag(void *m, const int *flag, int size) {
  size_t r = 0;
  for (size_t i = 0; i < _vectorSize(m)/size; ++i) {
    if(flag[i]) {
      if (r != 0) {
        memcpy(((char*)m)+size*(i-r),((char*)m)+size*i,size);
      }
    }
    else {
      r += 1;
    }
  }
  _vectorPop(m,size*r);
}
void quicksort (void *const pbase, size_t total_elems, size_t size,
            int (*cmp)(const void*, const void*, void*), void *arg);

#define vector_remove_flag(v,f,repeat) _vectorRemoveFlag((void*)v,f,repeat*sizeof(*v))
#define vector_size(v) (_vectorSize((void*)v)/sizeof(*v))
#define vector_push(v) ((__typeof__(*v))_vectorPush((void**)v, sizeof(**v)))
#define vector_push_n(v, x) ((__typeof__(*v))_vectorPush((void**)v, sizeof(**v) * (x)))
#define vector_insert(v, p) ((__typeof__(*v))_vectorInsert((void**)v, p * (sizeof(**v)), sizeof(**v)))
#define vector_pop(v) _vectorPop((void*)v, sizeof(*v))
#define vector_pop_n(v, x) _vectorPop((void*)v, sizeof(*v) * (x)))
#define vector_free(m) _vectorFree((void**)&m)
#define vector_dup(m) ((__typeof__(m))_vectorDup(m))
#endif
