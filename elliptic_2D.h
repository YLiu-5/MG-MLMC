#pragma once
#include "CS_2D.h"
#include <iostream>

#include <random>
#include <chrono>

void elliptic_2D(int l, int N, double *sums, std::mt19937& _rng);

double* construct_cov_mat(int nx, int ny);

void print_matrix_2D(double *array_2d, int nx, int ny);
void print_matrix_1D(double *array_1d, int n);

double* lowertri_times_vector(double *lowertri, double *vec, int N); //Lower triangular matrix times a vector, but the lower triangular matrix as well as the unused upper part are stored in a 1-d vector

double* upscale(double* k_fine, int nx, int ny);


