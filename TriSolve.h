// TriSolve.h - Solve tridiagonal equation systems using LU-factorization.

#ifndef TRISOLVE_H
#define TRISOLVE_H

class TriSolve
{
public:

  double *l, *d, *u;
  int N;

  TriSolve() : l(0), d(0), u(0) {}

  TriSolve(double* l_, double* d_, double* u_, int N_) { init(l_,d_,u_,N_); }

  void init(double* l_, double* d_, double* u_, int N_) {
    l = l_; d = d_; u = u_;
    N = N_;
  }

  void lu(double* l_, double* d_, double* u_, int N_) {
    l = l_; d = d_; u = u_;
    N = N_;
    lu();
  }

  void lu() {
    d[0] = 1.0/d[0];
    for (int n = 1; n < N; n++) {
      l[n-1] *= d[n-1];
      d[n] -= l[n-1]*u[n-1];
      d[n] = 1.0/d[n];
    }
      // d[n] -= l[n-1]*d[n-1]*u[n-1];
      // d[n] = 1.0/d[n];

    // d now contains 1/diagonal elements of U in LU factorization.
  }

  void solve(double* x) { // Right hand side b of the equation system, overwritten by solution x.
    for (int n = 1; n < N; n++)
      x[n] -= l[n-1]*x[n-1];

    x[N-1] *= d[N-1];
    for (int n = N-2; n >= 0; n--)
      x[n] = d[n]*(x[n] - u[n]*x[n+1]);
  }

  // Solve Ax=b directly.
  // On entry x should contain the right hand side b,
  // on exit x is overwritten with the solution x.
  // l, d, and u are preserved.
 void direct_solve(double* x) {
    double D[N];
    double dn = D[0] = 1.0/d[0];
    for (int n = 0; n < N-1; n++) {
      x[n+1] -= l[n]*dn*x[n];
      D[n+1] = dn = 1.0/(d[n+1] - l[n]*dn*u[n]);
    }
    x[N-1] *= dn;
    for (int n = N-2; n >= 0; n--)
      x[n] = D[n]*(x[n] - u[n]*x[n+1]);
  }

  // Solve Ax=b directly without first LU factorization:
  // On entry x should contain the right hand side b,
  // on exit x is overwritten with the solution x.
  // l and u are preserved, d is destroyed.
  void direct_solve_destructive(double* x) { // destroys diagonal d
    d[0] = 1.0/d[0];
    for (int n = 0; n < N-1; n++) {
      d[n+1] -= l[n]*d[n]*u[n];
      x[n+1] -= l[n]*d[n]*x[n];
      d[n+1] = 1.0/d[n+1];
    }
    x[N-1] *= d[N-1];
    for (int n = N-2; n >= 0; n--)
      x[n] = d[n]*(x[n] - u[n]*x[n+1]);
  }

};

// Solve Ax=b directly.
// On entry x should contain the right hand side b,
// on exit x is overwritten with the solution x.
// l, d, and u are preserved.
void tri_solve(double* l, double* d, double* u, int N, double* x) {
  double D[N];
  double dn = D[0] = 1.0/d[0];
  for (int n = 0; n < N-1; n++) {
    x[n+1] -= l[n]*dn*x[n];
    D[n+1] = dn = 1.0/(d[n+1] - l[n]*dn*u[n]);
  }
  x[N-1] *= dn;
  for (int n = N-2; n >= 0; n--)
    x[n] = D[n]*(x[n] - u[n]*x[n+1]);
}

// Solve Ax=b directly.
// On entry x should contain the right hand side b,
// on exit x is overwritten with the solution x.
// l and u are preserved, d is destroyed.
void tri_solve_destructive(double* l, double* d, double* u, int N, double* x) {
  double dn = d[0] = 1.0/d[0];
  for (int n = 0; n < N-1; n++) {
    x[n+1] -= l[n]*dn*x[n];
    d[n+1] = dn = 1.0/(d[n+1] - l[n]*dn*u[n]);
  }
  x[N-1] *= dn;
  for (int n = N-2; n >= 0; n--)
    x[n] = d[n]*(x[n] - u[n]*x[n+1]);
}

#endif
