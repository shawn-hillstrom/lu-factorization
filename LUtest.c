/* ACKNOWLEDGEMENT
 * ---------------
 * Written by Paul Bonamy
 * Computer Science Department, WSUV
 *
 * Written for the purpose of testing the correctness of LUfact.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "LUfact.h"

double A5_[5][5] = {
  {-1, -2, -3,  4,  5},
  { 2,  3, -4, -5, -6},
  {-3,  4, -5,  6, -7},
  {-4,  5,  6, -7,  8},
  { 5, -6,  7, -8,  9}
};

double A10_[10][10] = {
  { 8,  0,  8,-10,-10, -8,  9,  2,-10, -5},
  { 3,  5, -7, -5, -3,  2, -6,  1, -6,  4},
  { 1,  3,  9,  2, -1, -9,  0, 10, -9, 10},
  { 3,  3, -9, -5,  1, -2,  0, -2,  0, -6},
  { 2,  5, 10, -8, -4, -7,  4, -8, -4,  7},
  { 7,  2, -8,  1,  6, -3,  2, -2, 10,  5},
  {-2,  3, -7,  9, -4,  7,-10,  3,  5, -3},
  { 1, -1,  8,  0, -8, -6,  3, -5, -6, -2},
  {-7, -5,  0,  6,  6, -7,  8, -6, -2,  8},
  {-7,  8, -6,  6, -8, -1,  8, -8,  9,  7}
};

int main() {
  const double *A5[5] = {
    A5_[0], A5_[1], A5_[2], A5_[3], A5_[4],
  };
  const double b5[] = {1, 2, 3, 4, 5};

  //
  // Solution from wolfram-alpha 
  // http://goo.gl/Nwe6f5
  //
  printf("solving 5x5 system...\n");
  LUfact *LU = LUfactor(5, A5);
  double x5[5];
  LUsolve(LU, b5, x5);
  double x5_soln[5] = {
    263.0/12, 107.0/6, 61.0/20, 139.0/15, 92.0/15
  };
  for (int i = 0; i < 5; i++) {
    const double err = x5[i] - x5_soln[i];
    printf("x[i] = %11.7f (%11.7f, error=%0.7e)\n", 
	   x5[i], x5_soln[i], err);
  }
  LUdestroy(LU);

  //
  // compute inverse of A10.
  // Solve A*X = I where I is 10x10 identity matrix.
  //
  printf("computing 10x10 inverse...\n");
  const double *A10[10] = {
    A10_[0], A10_[1], A10_[2], A10_[3], A10_[4], 
    A10_[5], A10_[6], A10_[7], A10_[8], A10_[9],
  };
  LU = LUfactor(10, A10);
  double X10_[10][10]; 
  double *X10[10] = {  // will hold transpose of inverse
    X10_[0], X10_[1], X10_[2], X10_[3], X10_[4], 
    X10_[5], X10_[6], X10_[7], X10_[8], X10_[9],
  };
  for (int i = 0; i < 10; i++) {
    double b[10];
    for (int j = 0; j < 10; j++)
      b[j] = (i == j) ? 1.0 : 0.0;
    LUsolve(LU, b, X10[i]);
  }
  printf("checking solution (A*A^-1 = I)...\n");
  double I[10][10];
  double maxError = 0.0;
  for (int i = 0; i < 10; i++)
    for (int j = 0; j < 10; j++) {
      double sum = 0.0;
      for (int k = 0; k < 10; k++)
	sum += A10[i][k]*X10[j][k];
      I[i][j] = sum;
      const double error  = fabs(sum - ((i == j) ? 1.0 : 0.0));
      if (error > maxError)
	maxError = error;
    }
  printf("max error = %0.12e\n", maxError);
  /* XXX
  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 10; j++) {
      printf("%7.3f ", I[i][j]);
    }
    printf("\n");
  }
  */
  LUdestroy(LU);

  //
  // Big NxN system
  //
  const int N = 600;
  printf("solving %dx%d system...\n", N, N);
  double **S = (double **) malloc(N*sizeof(double*));
  for (int i = 0; i < N; i++) {
    S[i] = (double *) malloc(N*sizeof(double));
    const double f = i + 1;
    const double s = sin(f);
    for (int j = 0; j < N; j++)
      S[i][j] = sin(s*exp(j));
  }
  LU = LUfactor(N, (const double **) S);
  double *B = (double *) malloc(N*sizeof(double));
  for (int i = 0; i < N; i++)
    B[i] = i;
  double *X = (double *) malloc(N*sizeof(double));
  LUsolve(LU, B, X);
  printf("checking solution...\n");
  maxError = 0.0;
  for (int i = 0; i < N; i++) {
    double sum = 0.0;
    for (int j = 0; j < N; j++)
      sum += S[i][j]*X[j];
    const double error  = fabs(B[i] - sum);
    if (error > maxError)
      maxError = error;
  }
  LUdestroy(LU);
  printf("max error = %0.12e\n", maxError);
  
  return 0;
}
