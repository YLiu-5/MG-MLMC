#pragma once
#include <random>
void mlmc_test(void(*mlmc_l)(int, int, double *, std::mt19937&), int M, int N, int L,
	int N0, float *Eps, int Lmin, int Lmax, FILE *fp, std::mt19937&);