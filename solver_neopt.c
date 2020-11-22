/*
 * Tema 2 ASC
 * 2020 Spring
 */
#include "utils.h"

double* my_solver(int N, double *A, double* B) {
	int i, j, k;
	int NMAX = N * N;
	double *C = (double *)calloc(NMAX, sizeof(double));
	double *At = (double *)calloc(NMAX, sizeof(double));
	double *A2 = (double *)calloc(NMAX, sizeof(double));
	double *BxAt = (double *)calloc(NMAX, sizeof(double));
	double *A2xB = (double *)calloc(NMAX, sizeof(double));

	// Obtin matricele At si A^2
    for (i = 0; i < N; ++i) {
		// j >= i (superior triunghiulara)
        for (j = i; j < N; ++j) {
			// Calculez A transpus
			At[j * N + i] = A[i * N + j];

			// Calculez A^2
			// k >= i (superior triunghiulara)
			for (k = i; k < N; ++k) {
				A2[i * N + j] += A[i * N + k] * A[k * N + j];
			}
        }
    }

	// C = B * At + A2 * B
	for (i = 0; i < N; ++i) {
		for (j = 0; j < N; ++j) {
			for(k = 0; k < N; ++k) {
				// A matrice superior triunghiulara => At inferior triunghiulara
				// Calculez B * At
				if (k >= j)
					BxAt[i * N + j] += B[i * N + k] * At[k * N + j];

				// Tin cont de faptul ca A^2 este superior triunghiulara
				// Calculez A^2 * B
				if (i <= k)
					A2xB[i * N + j] += A2[i * N + k] * B[k * N + j];
			}

			// Retin adunarea celor 2 inmultiri in C
			C[i * N + j] = BxAt[i * N + j] + A2xB[i * N + j];
		}
	}

	free(At);
	free(A2);
	free(BxAt);
	free(A2xB);

	return C;
}
