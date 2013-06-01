/*
 *  NKA_C_TEST
 *
 *  A test program for nonlinear_krylov_accelerator.c.
 *
 *  Neil N. Carlson <neil.n.carlson@gmail.com>
 *
 *  This program exercises the nonlinear Krylov acceleration procedure
 *  implemented in nonlinear_krylov_accelerator.c by applying it to the
 *  solution of the nonlinear elliptic problem
 *
 *      -Div[(a + u)Grad[u]] = q,  u:[0,1]^2 -> R
 *
 *  with zero boundary values.  We use a mimetic discretization over a
 *  completely regular rectangular mesh.  The scalar field u is discretized
 *  in the cell-based space.  The discrete equations are of the form r(u) = 0,
 *  where r(u) = A(u) u - q.  A(u) is a symmetric positive definite matrix
 *  that depends on the value of the scalar field u.  To solve r(u) = 0 we
 *  use a fixed point iteration applied to a preconditioned r(u).  The
 *  preconditioner is one or more passes of SSOR to approximately solve
 *  A(u)^{-1} r(u) (here u is fixed).  We use a=0.02 and a uniform unit
 *  source field q.
 *
 *  We first solve the equation with an accelerated iteration, and then,
 *  for comparison, solve again without acceleration.  In each case, the
 *  l2 norm of r(u), the reduction factor in the norm, and the net convergence
 *  rate are printed for each iteration.
 *
 *  This code is a straightforward translation of the original Fortran test
 *  program into C, and the resulting code is fairly inscrutible.  Use the
 *  Fortran code as a reference to understand the details of this code.
 *
 *******************************************************************************
 *
 *  Copyright (c) 2009  Neil N. Carlson
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
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "nonlinear_krylov_accelerator.h"

typedef struct {
  int nx, ny;
  double a, hx, hy;
  double *ax, *ay, *ac, *q;
  NKA nka;
} SYSTEM;

void solve (SYSTEM *, int, double, int, double *);
void residual (SYSTEM *, double *, double *);
void update_system (SYSTEM *, double *);
void ssor_pc (SYSTEM *, int, double, double *);
double l2norm (double *, int);


int main (void)
{
  int N=50;
  double A=0.02;
  int i, nsweep=2;
  double omega=1.4;
  SYSTEM sys;
  double *u;

  /* Define the SYSTEM structure */
  sys.a  = A;
  sys.nx = N;
  sys.ny = N;
  sys.hx = 1.0/sys.nx;
  sys.hy = 1.0/sys.ny;
  sys.ax = malloc(sys.ny*(sys.nx+1)*sizeof(double));
  sys.ay = malloc(sys.nx*(sys.ny+1)*sizeof(double));
  sys.ac = malloc(sys.nx*sys.ny*sizeof(double));
  sys.q  = malloc(sys.nx*sys.ny*sizeof(double));
  sys.nka = nka_init(sys.nx*sys.ny, 5, 0.01, NULL);
  /* uniform unit source */
  for (i = 0; i < sys.nx*sys.ny; i++) sys.q[i] = sys.hx*sys.hy;

  /* Solve with acceleration and then again without */
  u = malloc(sys.nx*sys.ny*sizeof(double));
  solve(&sys, nsweep, omega, 1, u);
  solve(&sys, nsweep, omega, 0, u);
  nka_delete(sys.nka);

  return 0;
}


void solve (SYSTEM *sys, int nsweep, double omega, int use_nka, double *sol)
{
  int i, j, k, n, nx, ny, itr;
  double *upad, *u, *r, *t;
  double rnorm0, rnorm, red, rate;
  int MAXITR=999;
  double TOL=1.0e-6;

  if (use_nka) {
    printf("\nACCELERATED SOLVE\n\n");
  } else {
    printf("\nUNACCELERATED SOLVE\n\n");
  }
  printf("Iter Residual Norm    Reduction    Rate\n");

  nx = sys->nx;
  ny = sys->ny;

  n = (nx+2)*(ny+2);
  upad = malloc(n*sizeof(double));
  for (i = 0; i < n; i++) upad[i] = 0.0;

  r = malloc(nx*ny*sizeof(double));

  residual(sys, upad, r);
  rnorm0 = l2norm(r, nx*ny);
  printf("%3d:%14.6E\n", 0, rnorm0);

  for (itr = 1; itr <= MAXITR; itr++) {
    ssor_pc(sys, nsweep, omega, r);
    if (use_nka) nka_accel_update(sys->nka, r);

    /* solution iterate update: u_{i+1} := u_{i} - v */
    t = r;
    u = &upad[nx+2];
    for (k = 0; k < ny; k++) {
      u++;
      for (j = 0; j < nx; j++) {
        *u -= *t;
        u++;
        t++;
      }
      u++;
    }

    residual(sys, upad, r);
    rnorm = l2norm(r, nx*ny);
    red = rnorm / rnorm0;
    rate = pow(red, 1.0/itr);
    printf("%3d:%14.6E%13.3E%8.3f\n",itr, rnorm, red, rate);
    if (rnorm < TOL*rnorm0) break;
  }

  /* copy result to return array */
  u = &upad[nx+2];
  for (k = 0; k < ny; k++) {
    u++;
    for (j = 0; j < nx; j++) {
      *sol = *u;
      u++;
      sol++;
    }
    u++;
  }
}


void residual (SYSTEM *sys, double *upad, double *r)
{
  int j, k, nx, ny;
  double *u, *ac, *ax, *ay, *q;

  update_system (sys, upad);

  nx = sys->nx;
  ny = sys->ny;

  ac = sys->ac;
  ax = sys->ax;
  ay = sys->ay;
  q  = sys->q;
  u  = &upad[nx+2];
  for (k = 0; k < ny; k++) {
    u++;
    for (j = 0; j < nx; j++) {
      r[0] = ac[0]*u[0] - ax[0]*u[-1] - ax[1]*u[1] - ay[0]*u[-nx-2] - ay[nx]*u[nx+2] - q[0];
      r++;
      u++;
      q++;
      ac++;
      ax++;
      ay++;
    }
    u++;
    ax++;
  }

}


void update_system (SYSTEM *sys, double *upad)
{
  int i, j, k, nx, ny;
  double rx, ry, t;
  double *u, *ax, *ay, *ac;

  nx = sys->nx;
  ny = sys->ny;

  rx = sys->hx / sys->hy;
  ry = sys->hy / sys->hx;

  for (i = 0; i < ny*(nx+1); i++) {
    sys->ax[i] = 0.0;
  }
  for (i = 0; i < nx*(ny+1); i++) {
    sys->ay[i] = 0.0;
  }

  ax = sys->ax;
  ay = sys->ay;
  u  = &upad[nx+2];
  for (k = 0; k < ny; k++){
    u++;
    for (j = 0; j < nx; j++) {
      t = 1.0 / (sys->a + u[0]);
      ax[0] += rx*t;
      ax[1] += rx*t;
      ay[0] += ry*t;
      ay[nx] += ry*t;
      ax++;
      ay++;
      u++;
    }
    u++;
    ax++;
  }

  for (i = 0; i < ny*(nx+1); i++) {
    sys->ax[i] = 2.0 / sys->ax[i];
  }

  for (i = 0; i < nx*(ny+1); i++) {
    sys->ay[i] = 2.0 / sys->ay[i];
  }

  ac = sys->ac;
  ax = sys->ax;
  ay = sys->ay;
  for (k = 0; k < ny; k++) {
    for (j = 0; j < nx; j++) {
      ac[0] = ax[0] + ax[1] + ay[0] + ay[nx];
      ac++;
      ax++;
      ay++;
    }
    ax++;
  }
}


void ssor_pc (SYSTEM *sys, int nsweep, double omega, double *r)
{
  int i, j, k, n, nx, ny, nbuf;
  double *ac, *ax, *ay, *zpad, *z;

  nx = sys->nx;
  ny = sys->ny;

  n = (nx+2)*(ny+2);
  zpad = malloc(n*sizeof(double));
  for (i = 0; i < n; i++) zpad[i] = 0.0;

  ac = sys->ac;
  ax = sys->ax;
  ay = sys->ay;
  z  = &zpad[nx+2];

  for (i = 0; i < nsweep; i++) {
    /* Forward sweep */
    for (k = 0; k < ny; k++) {
      z++;
      for (j = 0; j < nx; j++) {
        z[0] = (1.0 - omega) * z[0]
             + omega * (r[0] + ax[0] * z[-1] + ax[1] * z[1]
                             + ay[0] * z[-nx-2] + ay[nx] * z[nx+2]) / ac[0];
        z++;
        r++;
        ac++;
        ax++;
        ay++;
      }
      z++;
      ax++;
    }
    /* Backward sweep */
    for (k = 0; k < ny; k++) {
      z--;
      ax--;
      for (j = 0; j < nx; j++) {
        z--;
        r--;
        ac--;
        ax--;
        ay--;
        z[0] = (1.0 - omega) * z[0]
             + omega * (r[0] + ax[0] * z[-1] + ax[1] * z[1]
                             + ay[0] * z[-nx-2] + ay[nx] * z[nx+2]) / ac[0];
      }
      z--;
    }
  }

  /* copy result to return array */
  for (k = 0; k < ny; k++) {
    z++;
    for (j = 0; j < nx; j++) {
      *r = *z;
      r++;
      z++;
    }
    z++;
  }
  free(zpad);
}


double l2norm (double *x, int len)
{
  int i;
  double a;
  a = 0.0;
  for (i = 0; i < len; i++) {
    a += x[i]*x[i];
  }
  return sqrt(a);
}




