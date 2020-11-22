/*
 * Tema 2 ASC
 * 2020 Spring
 */
#include "utils.h"

double* my_solver(int N, double *A, double* B) {
	register int i, j, k, row, NMAX = N * N;
	register double sum = 0.0;
	register double *pA, *pi, *pj, *pAt, *pA2;
	double *At = (double *)calloc(NMAX, sizeof(double));
	double *A2 = (double *)calloc(NMAX, sizeof(double));
	double *C = (double *)calloc(NMAX, sizeof(double));

	// Obtin matricele At si A^2
    for (i = 0; i < N; ++i) {
		// pointer pentru linie, in matricea A
		row = i * N;
		pA = &A[row];
        for (j = 0; j < N; ++j) {
			// Tin cont de faptul ca A este superior triunghiulara
            if (i <= j) {
				// initializez A transpus
                At[j * N + i] = *pA;

				// ponteri pentru linia i (pi) si coloana j (pj)
				pi = &A[row];
				pj = &A[j];

				// Calculez A^2, dupa detectarea constantei A2[row + j],
				// pe care o accesez la fiecare N operatii
				sum = 0.0;
				for (k = 0; k < N; ++k) {
					// tin cont de faptul ca A este superior triunghiulara
					if (i <= k)
						sum += *pi * *pj;

					// Urmatoarea linie
					++pi;
					// urmatoarea coloana
					pj += N;
				}

				// Accesez matricea doar pentru a pune suma finala
				A2[row + j] = sum;
            }

			// Trec la urmatoarea linie
			++pA;
        }
    }

	// C = B * At + A2 * B
	for (i = 0; i < N; ++i) {
		row = i * N;
		for (j = 0; j < N; ++j) {
			// pointer pentru linia in matricea B
			pi = &B[row];
			// pointer pentru linie in matricea A
			pA2 = &A2[row];
			// pointer pentru coloana in matricea B
			pj = &B[j];
			// pointer pentru coloana in matricea At
			pAt = &At[j];

			// sum = BxAt[i][j] + A2xB[i][j]
			sum = 0.0;
			for(k = 0; k < N; ++k) {
				// BxAt cu At inferior triunghiulara
				if (k >= j)
					sum += *pi * *pAt;

				// A2xB cu A2 superior triunghiulara
				if (i <= k)
					sum += *pA2 * *pj;

				// Urmatoarea linie in matricea B
				++pi;
				// Urmatoarea linie in matricea A2
				++pA2;
				// Urmatoarea coloana in matricea B
				pj += N;
				// Urmatoarea coloana in matricea At
				pAt += N;
			}

			// C[i][j] = sum
			C[row + j] = sum;
		}
	}

	free(At);
	free(A2);

	return C;
}