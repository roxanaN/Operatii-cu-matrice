/*
 * Tema 2 ASC
 * 2020 Spring
 */
#include "utils.h"
#include "cblas.h"
#include <string.h>

double* my_solver(int N, double *A, double *B) {
	int i, j, row, index, NMAX = N * N;
	double *BxAt = (double *)calloc(NMAX, sizeof(double));
	double *A2 = (double *)calloc(NMAX, sizeof(double));

	// Copiez matricea B in BxAt,
	// pentru a nu o deteriora in timpul operatiei B x At
	memcpy(BxAt, B, NMAX * sizeof(double));
	// Fac inmultirea B * At si retin rezultatul in matricea BxAt
	cblas_dtrmm(CblasRowMajor, CblasRight, CblasUpper, CblasTrans, CblasNonUnit,
				N, N, 1, A, N, BxAt, N);

	// Copiez matricea A in A2 (A^2),
	// pentru a nu pierde valorile din A, in timpul ridicarii la patrat
	memcpy(A2, A, NMAX * sizeof(double));
	// Ridic matricea A la patrat
	cblas_dtrmm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit,
				N, N, 1, A, N, A2, N);
	
	// Calculez A^2 * B si retin rezultatul in B
	cblas_dtrmm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit,
				N, N, 1, A2, N, B, N);

	// Eliberez memoria alocata pentru matricea A2
	free(A2); 
	
	// B * At + A^2 * B = BxAt + B = BxAt
	for (i = 0; i < N; ++i) {
		row = i * N;
		for (j = 0; j < N; ++j) {
			index = row + j;
			BxAt[index] += B[index];
		}
	}

	// Returnez rezultatul retinut in BxAt
	return BxAt;
}
