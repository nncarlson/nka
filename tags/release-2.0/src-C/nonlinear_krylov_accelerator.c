/*
 *  NONLINEAR_KRYLOV_ACCELERATOR
 *
 *  Neil N. Carlson <neil.n.carlson@gmail.com>
 *
 *  This module implements the nonlinear Krylov accelerator (NKA) introduced
 *  in [1] for fixed point or Picard iterations.  Placed in the iteration loop,
 *  this black-box accelerator listens to the sequence of solution updates and
 *  replaces them with accelerated updates.  More generally, NKA can accelerate
 *  typical quasi-Newton iterations, which can usually be viewed as a fixed
 *  point iteration for a preconditioned function.
 *
 *  This code is a straightforward translation of the original Fortran 95
 *  implementation into C.
 *
 *  [1] N.N.Carlson and K.Miller, "Design and application of a gradient-
 *      weighted moving finite element code I: in one dimension", SIAM J.
 *      Sci. Comput;, 19 (1998), pp. 728-765.  See section 9.
 *
 *******************************************************************************
 *
 *  Copyright (c) 2009, 2013  Neil N. Carlson
 *
 *  Permission is hereby granted, free of charge, to any person obtaining a
 *  copy of this software and associated documentation files (the "Software"),
 *  to deal in the Software without restriction, including without limitation
 *  the rights to use, copy, modify, merge, publish, distribute, sublicense,
 *  and/or sell copies of the Software, and to permit persons to whom the
 *  Software is furnished to do so, subject to the following conditions:
 *
 *  The above copyright notice and this permission notice shall be included
 *  in all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 *  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 *  DEALINGS IN THE SOFTWARE.
 *
 *******************************************************************************
 *
 *  PROGRAMING INTERFACE
 *
 *  #include "nonlinear_krylov_accelerator.h"
 *
 *  NKA state = nka_init (int vlen, int mvec, double vtol,
 *                          double (*dp)(int len, double *x, double *y));
 *
 *      The function nka_init() creates a new instance of the accelerator
 *      and returns an opaque handle state (of type NKA -- a structure pointer)
 *      to the accelerator.  The acceleration subspace will contain up to mvec
 *      vectors (old vectors will be dropped if necessary) of length vlen.
 *      The argument vtol specifies the vector drop tolerance: a vector is
 *      dropped when the sine of the angle between the vector and the subspace
 *      spanned by the preceding vectors is less than this value.  Acceptable
 *      values are in the interval (0,1) and a reasonable default is 0.01 or
 *      smaller.
 *
 *      The function pointer dp should either be NULL or point to a function
 *      that returns the dot product of two vectors x and y of length len.
 *      If NULL, an internal dot product function will be used; otherwise
 *      the specified function will be used to compute vector dot products.
 *      In a parallel implementation of an interative nonlinear solver where
 *      the vector components are distributed across processors, a global dot
 *      product function that performs the require communication internally
 *      should be specified, for example.
 *
 *  void nka_delete (NKA state);
 *
 *      The function nka_delete() frees all the memory associated with the
 *      accelerator state.
 *
 *  void nka_accel_update (NKA state, double *f);
 *
 *      The function nka_accel_update() takes the function value f, which would
 *      be the update vector in a fixed point iteration, and overwrites it with
 *      the accelerated update computed from the acceleration subspace stored
 *      in the accelerator state.  This acceleration subspace is updated prior
 *      to computing the update using f and the previous function value and
 *      update that were cached on the preceding call to nka_accel_update(),
 *      if any.  The input function value f and the returned update are then
 *      cached in the accelerator for use on the next call.
 *
 *  void nka_relax (NKA state);
 *
 *      The function nka_relax deletes the pending vectors that were cached by
 *      the preceding call to nka_relax(), if any.  This modifies the behavior
 *      of the next call to nka_relax() in that the acceleration subspace will
 *      not be updated prior to computing the accelerated update.  This could
 *      be used, for example, to carry over the subspace from one nonlinear
 *      solve to another.  (Whether this is an effective strategy is an open
 *      question.)  nka_relax() expects that the passed function value is
 *      connected to the preceding update (if it exists), but this is not
 *      normally true for the first call in a subsequent nonlinear solve, and
 *      would result in the subspace being updated with bogus information.
 *      A call to nka_relax() at the end of a nonlinear solve prevents this
 *      from occuring.
 *
 *  void nka_restart (NKA state);
 *
 *      The function nka_restart() flushes the acceleration subspace from the
 *      accelerator state, restoring it to its newly created condition; the
 *      next call to nka_accel_update() begins the process of accumulating a
 *      new subspace.  Typical usage is to call nka_restart() at the start of
 *      each nonlinear solve in a sequence of solves.  This allows a single
 *      instance of the accelerator to be reused and eliminates the overhead
 *      of repeated memory allocation and deallocation that would otherwise
 *      occur.
 *
 *  int nka_num_vec (NKA state);
 *
 *      The function nka_num_vec() returns the number of vectors in the
 *      acceleration subspace.
 * 
 *  int nka_max_vec (NKA state);
 * 
 *      The function nka_max_vec() returns the max number of vectors in the
 *      acceleration subspace.
 * 
 *  int nka_vec_len (NKA state);
 * 
 *      The function nka_vec_len returns the length of the vectors.
 * 
 *  double nka_vec_tol (NKA state);
 * 
 *      The function nka_vec_tol() returns the vector drop tolerance.
 *
 *  USAGE
 *
 *  The following simple example shows the usage of this acceleration
 *  procedure.  For more details, see the associated documentation.
 *  Consider an quasi-Newton iteration for solving the nonlinear system
 *  F(x) = 0.  Suppose PC(y) is some preconditioning procedure; for example
 *  the application of some approximation of the inverse of the Jacobian of
 *  F(x) to the vector y.  The original quasi-Newton iteration (equivalent
 *  to the fixed point iteration for PC(F(x)) = 0) would look something like
 *
 *    x = 0
 *    do <until converged>
 *      dx = PC(F(x))
 *      x = x - dx
 *    end do
 *
 *  The accelerated iteration would look something like
 *
 *    state = nka_init(size(v), 5, 0.01, NULL)
 *    x = 0
 *    do <until converged>
 *      dx = PC(F(x))
 *      nka_accel_update(state, dx)
 *      x = x - dx
 *    end do
 *    nka_delete(state)
 *
 *  The init and delete can of course be moved outside any nonlinear
 *  solution procedure containing this iteration, and a single state instance
 *  used for repeated calls to the procedure.  This avoids the repeated
 *  allocations and deallocations of arrays associated with the state variable.
 *  In this case, one should either include a call to nka_restart() before the
 *  loop so that each iterative solve starts with clean slate, or include a
 *  call to nka_relax() after the loop so that first call to nka_accel_update()
 *  in the next iterative solve doesn't update the acceleration subspace with
 *  bogus information.
 *
 */

#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "nonlinear_krylov_accelerator.h"

#define TRUE 1
#define FALSE 0
#define EOL -1  /* end-of-list marker */

struct nka_state {
  int subspace;       /* boolean: a nonempty subspace */
  int pending;        /* contains pending vectors -- boolean */
  int vlen;           /* vector length */
  int mvec;           /* maximum number of subspace vectors */
  double vtol;        /* vector drop tolerance */
  /* Subspace storage */
  double **v;   /* correction vectors */
  double **w;   /* function difference vectors */
  double **h;   /* matrix of w vector inner products */
  /* Linked-list organization of the vector storage. */
  int first;  /* index of first subspace vector */
  int last;   /* index of last subspace vector */
  int free;   /* index of the initial vector in free storage linked list */
  int *next;  /* next index link field */
  int *prev;  /* previous index link field in doubly-linked subspace v */
  /* Dot product function pointer */
  double (*dp)(int, double *, double *);
};


double dot_product (int len, double *a, double *b)
{
  int j;
  double s;
  assert(len >= 0);
  for (s = 0.0, j = 0; j < len; j++)
    s += a[j] * b[j];
  return s;
}


NKA nka_init (int vlen, int mvec, double vtol,
    double (*dp)(int, double *, double *))
{
  int j, n;
  NKA state;

  assert(mvec > 0);
  assert(vlen >= 0);
  assert(vtol > 0.0);

  state = (NKA)malloc(sizeof(*state));

  state->vlen = vlen;
  state->mvec  = mvec;
  state->vtol  = vtol;

  if (dp == NULL) {
    state->dp = &dot_product;
  } else {
    state->dp = dp;
  }

  n = mvec + 1;
  state->v = (double**)malloc(n*sizeof(double *));
  state->v[0] = (double*)malloc(n * vlen * sizeof(double));
  for (j = 1; j < n; j++) {
    state->v[j] = state->v[j-1] + vlen;
  }

  state->w = (double**)malloc(n*sizeof(double *));
  state->w[0] = (double*)malloc(n * vlen * sizeof(double));
  for (j = 1; j < n; j++) {
    state->w[j] = state->w[j-1] + vlen;
  }

  state->h = (double**)malloc(n * sizeof(double *));
  state->h[0] = (double*)malloc(n * n * sizeof(double));
  for (j = 1; j < n; j++) {
    state->h[j] = state->h[j-1] + n;
  }

  state->next = (int*)malloc(n * sizeof(int));
  state->prev = (int*)malloc(n * sizeof(int));

  nka_restart(state);

  return state;
}


void nka_delete (NKA state)
{
  if (state->v) {
    free(state->v[0]);
    free(state->v);
  }

  if (state->w) {
    free(state->w[0]);
    free(state->w);
  }

  if (state->h) {
    free(state->h[0]);
    free(state->h);
  }

  free(state->next);
  free(state->prev);

  free(state);
}


void nka_accel_update (NKA state, double *f)
{
  int i, j, k, nvec, newv;
  double s, hkk, hkj, cj;
  double *v, *w, *hk, *hj, *c;

  /*
   *  UPDATE THE ACCELERATION SUBSPACE
   */

  if (state->pending) {

    /* next function difference w_1 */
    w = state->w[state->first];
    for (j = 0; j < state->vlen; j++)
      w[j] -= f[j];
    s = sqrt(state->dp(state->vlen, w, w));

  /* If the function difference is 0, we can't update the subspace with
     this data; so we toss it out and continue.  In this case it is likely
     that the outer iterative solution procedure has gone badly awry
     (unless the function value is itself 0), and we merely want to do
     something reasonable here and hope that situation is detected on the
     outside. */
    if (s == 0.0) nka_relax(state);

  }

  if (state->pending) {

    /* Normalize w_1 and apply same factor to v_1. */
    v = state->v[state->first];
    for (j = 0; j < state->vlen; j++) {
      v[j] /= s;
      w[j] /= s;
    }

    /* Update H. */
    for (k = state->next[state->first]; k != EOL; k = state->next[k])
      state->h[state->first][k] = state->dp(state->vlen, w, state->w[k]);

    /*
     *  CHOLESKI FACTORIZATION OF H = W^t W
     *  original matrix kept in the upper triangle (implicit unit diagonal)
     *  lower triangle holds the factorization
     */

    /* Trivial initial factorization stage. */
    nvec = 1;
    state->h[state->first][state->first] = 1.0;

    for (k = state->next[state->first]; k != EOL; k = state->next[k]) {

      /* Maintain at most MVEC vectors. */
      if (++nvec > state->mvec) {
        /* Drop the last vector and update the free storage list. */
        assert(state->last == k);
        state->next[state->last] = state->free;
        state->free = k;
        state->last = state->prev[k];
        state->next[state->last] = EOL;
        break;
      }

      /* Single stage of Choleski factorization. */
      hk = state->h[k];   /* row k of H */
      hkk = 1.0;
      for (j = state->first; j != k; j = state->next[j]) {
        hj = state->h[j];   /* row j of H */
        hkj = hj[k];
        for (i = state->first; i != j; i = state->next[i])
          hkj -= hk[i] * hj[i];
        hkj /= hj[j];
        hk[j] = hkj;
        hkk -= hkj*hkj;
      }

      if (hkk > pow(state->vtol,2)) {
        hk[k] = sqrt(hkk);
      } else  {
        /* The current w nearly lies in the span of the previous vectors: */
        /* Drop this vector, */
        assert(state->prev[k] != EOL);
        state->next[state->prev[k]] = state->next[k];
        if (state->next[k] == EOL)
          state->last = state->prev[k];
        else
          state->prev[state->next[k]] = state->prev[k];
        /* update the free storage list, */
        state->next[k] = state->free;
        state->free = k;
        /* back-up and move on to the next vector. */
        k = state->prev[k];
        nvec--;
      }
    }

    assert(state->first != EOL);
    state->subspace = TRUE; /* the acceleration subspace isn't empty */

  }

  /*
   *  ACCELERATED CORRECTION
   */

  /* Locate storage for the new vectors. */
  assert(state->free != EOL);
  newv = state->free;
  state->free = state->next[state->free];

  /* Save the original f for the next call. */
  for (j = 0; j < state->vlen; j++)
    state->w[newv][j] = f[j];

  if (state->subspace) {
    c = (double*)malloc((state->mvec + 1) * sizeof(*c));
    assert(c != NULL);
    /* Project f onto the span of the w vectors: */
    /* forward substitution */
    for (j = state->first; j != EOL; j = state->next[j]) {
      cj = state->dp(state->vlen, f, state->w[j]);
      for (i = state->first; i != j; i = state->next[i])
        cj -= state->h[j][i] * c[i];
      c[j] = cj / state->h[j][j];
    }
    /* backward substitution */
    for (j = state->last; j != EOL; j = state->prev[j]) {
      cj = c[j];
      for (i = state->last; i != j; i = state->prev[i])
        cj -= state->h[i][j] * c[i];
      c[j] = cj / state->h[j][j];
    }
    /* The accelerated correction */
    for (k = state->first; k != EOL; k = state->next[k]) {
      w = state->w[k];
      v = state->v[k];
      for (j = 0; j < state->vlen; j++)
        f[j] += c[k] * (v[j] - w[j]);
    }
    free(c);
  }

  /* Save the accelerated correction for the next call. */
  for (j = 0; j < state->vlen; j++)
    state->v[newv][j] = f[j];

  /* Prepend the new vectors to the list. */
  state->prev[newv] = EOL;
  state->next[newv] = state->first;
  if (state->first == EOL) {
    state->last = newv;
  } else {
    state->prev[state->first] = newv;
  }
  state->first = newv;

  /* The original f and accelerated correction are cached for the next call. */
  state->pending = TRUE;
}


void nka_restart (NKA state)
{
  int k;

  /* No vectors are stored. */
  state->first    = EOL;
  state->last     = EOL;
  state->subspace = FALSE;
  state->pending  = FALSE;

  /* Initialize the free storage linked list. */
  state->free = 0;
  for (k = 0; k < state->mvec; k++) {
    state->next[k] = k + 1;
  }
  state->next[state->mvec] = EOL;
}


void nka_relax (NKA state)
{
  int newv;

  if (state->pending) {
    /* Drop the initial slot where the pending vectors are stored. */
    assert(state->first >= 0);
    newv = state->first;
    state->first = state->next[state->first];
    if (state->first == EOL) {
      state->last = EOL;
    } else {
      state->prev[state->first] = EOL;
    }
    /* Update the free storage list. */
    state->next[newv] = state->free;
    state->free = newv;
    state->pending = FALSE;
  }
}


int nka_num_vec (NKA state)
{
  int k, n;
  n = 0;
  k = state->first;
  while (k != EOL) {
    n = n + 1;
    k = state->next[k];
  }
  if (state->pending) n = n - 1;
  return n;
}


int nka_max_vec (NKA state) { return state->mvec; }

int nka_vec_len (NKA state) { return state->vlen; }

double nka_vec_tol (NKA state) { return state->vtol; }
