/* Shawn Hillstrom -- CS 330 -- Program 2 Header File
 * --------------------------------------------------
 * Includes all necessary definitions, structures, and prototypes to implement the LU
 *	factorization library
 */

#ifndef LUFACT_H
#define LUFACT_H

typedef struct {
	int N;         /* Size of input NxN matrix A */
	double **LU;   /* NxN matrix holding combined L and U matrix */
	short *mutate; /* Row permutations of A */
} LUfact;

/*
 * Given NxN matrix A (stored as an array of N row ptrs),
 * returns LU factorization information.
 * If A is singluar, NULL is returned.
 */
LUfact *LUfactor(int N, const double **A);

/*
 * Deallocate decomposition info created by LUdecompose.
 */
void LUdestroy(LUfact *fact);

/*
 * Given LU decomposition info for A, 
 * solves linear system Ax = b for x.
 * x and b vectors are assumed to be length N.
 */
void LUsolve(LUfact *fact, const double *b, double *x);

#endif /* LUFACT_H */
