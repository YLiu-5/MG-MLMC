#pragma once

/*
P = mlmc(Lmin,Lmax,N0,eps, mlmc_l, alpha,beta,gamma, Nl)

multilevel Monte Carlo control routine

Lmin  = minimum level of refinement       >= 2
Lmax  = maximum level of refinement       >= Lmin
N0    = initial number of samples         > 0
eps   = desired accuracy (rms error)      > 0

alpha -> weak error is  O(2^{-alpha*l})
beta  -> variance is    O(2^{-beta*l})
gamma -> sample cost is O(2^{gamma*l})

if alpha, beta, gamma are not positive then they will be estimated

mlmc_l(l,N,sums)   low-level function
l       = level
N       = number of paths
sums[0] = sum(cost)
sums[1] = sum(Y)
sums[2] = sum(Y.^2)
where Y are iid samples with expected value:
E[P_0]           on level 0
E[P_l - P_{l-1}] on level l>0

P     = value
Nl    = number of samples at each level
NlCl  = total cost of samples at each level

*/
#include <random>
void regression(int, float *, float *, float &a, float &b);

float mlmc(int Lmin, int Lmax, int N0, float eps,
	void(*mlmc_l)(int, int, double *, std::mt19937&),
	float alpha_0, float beta_0, float gamma_0, int *Nl, float *Cl, std::mt19937& _rng);
