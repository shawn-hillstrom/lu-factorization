/* Shawn Hillstrom -- CS 330 -- Program 2 Library Implementation 
 * -------------------------------------------------------------
 * Includes all necessary implementation for the LU factorization library.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "LUfact.h"

double **createMatrix(int N) {
	double **M = (double **) malloc(N*sizeof(double*));
	for (int i = 0; i < N; i++)
		M[i] = (double*) malloc(N*sizeof(double));
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			M[i][j] = (i == j) ? 1.0 : 0.0;
	return M;
}

void destroyMatrix(int N, double **M) {
	for (int i = 0; i < N; i++)
		free(M[i]);
	free(M);
}

LUfact *LUfactor(int N, const double **A) {
	
/* Initialization */
	LUfact *LU = (LUfact*) malloc(sizeof(LUfact));
  	LU->N = N;
  	LU->LU = createMatrix(N);
  	LU->mutate = (short *) malloc(N*sizeof(short));

  	double **A_ = LU->LU; // Local pointer to LU matrix.
  
/* Clone matrix A into LU matrix */
  	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			A_[i][j] = A[i][j];

/* Form the mutate array */
  	for (int i = 0; i < N; i++)
    	LU->mutate[i] = (short) i;
  
/* Begin LU factorization */
  	for (int i = 0; i < N; i++) {
	 	  
	/* Find our pivot row */
		double pivotVal = 0.0;
		int pivotIndex = i;
		
	  	for (int k = i; k < N; k++) {
		  	if (fabs(A_[k][i]) > pivotVal) {
			  	pivotVal = fabs(A_[k][i]);
				pivotIndex = k;
		  	}
	  	}
		
		if (pivotVal <= 0) { // Bad system, terminate immediately.
			LUdestroy(LU);
			return NULL;
		}
		
	/* Pivot */
		if (pivotIndex != i) {
		/* Pivot the mutate array */
			short tempShort = LU->mutate[i];
			LU->mutate[i] = LU->mutate[pivotIndex];
			LU->mutate[pivotIndex] = tempShort;
		/* Pivot the LU matrix */
			double *tempPtr = A_[i];
			A_[i] = A_[pivotIndex];
			A_[pivotIndex] = tempPtr;
		}
		
	/* Fill in L and U in matrix A_ */
		for (int j = i + 1; j < N; j++) {
			A_[j][i] /= A_[i][i];
			for (int k = i + 1; k < N; k++) {
				A_[j][k] -= A_[j][i] * A_[i][k];
			}
		}
	  
  }

  return LU; //  Return dat boi.
  
}

void LUdestroy(LUfact *fact) {
	free(fact->mutate);
	destroyMatrix(fact->N, fact->LU);
	free(fact);
}

void LUsolve(LUfact *fact, const double *b, double *x) {
	
/* Initializaton */
	double **A_ = fact->LU;
	
/* Forward substitution */
	for (int i = 0; i < fact->N; i++) {
		x[i] = b[fact->mutate[i]];
		for (int k = 0; k < i; k++) {
			x[i] -= A_[i][k] * x[k];
		}
	}
	
/* Back substitution */
	for (int i = fact->N - 1; i >= 0; i--) {
		for (int k = i + 1; k < fact->N; k++) {
			x[i] -= A_[i][k] * x[k];
		}
		x[i] = x[i]/A_[i][i];
	}
	
}
