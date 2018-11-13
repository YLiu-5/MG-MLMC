#include "elliptic_2D.h"
#include "mkl.h"

//#include "mathimf.h"
//#define K5
//#define DEBUG_Covariance
void elliptic_2D(int l, int N, double * sums, std::mt19937& _rng)  //level l, number of samples N
{

	int   M, nf, nc;
	float T, r, sig, B, hf, hc, X0, Xf, Xc, Af, Ac, Mf, Mc, Bf, Bc,
		Xf0, Xc0, Xc1, vf, vc, dWc, ddW, Pf, Pc, dP;

	std::uniform_real_distribution<> uniform_dis(0, 1);
	std::normal_distribution<> normal_dis(0, 1);
	//std::mt19937 rng(1234);

	nf = 1 << (l+2);  //in 2D model, the DoF is 2^(l+2) * 2^(l+2)
	nc = nf / 2;      //

	for (int k = 0; k<7; k++) sums[k] = 0.0;
#ifdef DEBUG_Covariance
double* ksum = new double[nf*nf]();
double* xisum = new double[nf*nf]();
#endif
	for (int np = 0; np < N; np++) 
	{
		// Generate random number for nf, then upscale to nc
		// Construct covariance matrix
		double* cov_mat = construct_cov_mat(nf,nf);
		//print_matrix_2D(cov_mat,nf*nf,nf*nf);
		
		int	info = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L', nf*nf, cov_mat, nf*nf);

		//print_matrix_2D(cov_mat, nf*nf, nf*nf);
		//print_matrix_1D(cov_mat, nf*nf*nf*nf);

		//Generate random normal Gaussian vector
		double* xi = new double[nf*nf]();
		for (int i = 0; i < nf*nf; i++)
		{
			xi[i] = 0.01*normal_dis(_rng);
		}

		double* kf = new double[nf*nf]();
		kf = lowertri_times_vector(cov_mat, xi, nf*nf);

		//print_matrix_1D(kf, nf*nf);
		for (int i = 0; i < nf*nf; i++)
		{
			kf[i] = exp(kf[i]);
		}
		//print_matrix_1D(kf, nf*nf);
#ifdef  K5
		for (int i = 0; i < nf*nf; i++)
		{
			kf[i] = 5;
		}
#endif //  K5

#ifdef DEBUG_Covariance
		for (int i = 0; i < nf*nf; i++)
		{
			xisum[i] += xi[i];
			ksum[i] += kf[i];
		}
			
#endif // DEBUG_Covariance

		
		//print_matrix_1D(kf, nf*nf);
		//print_matrix_1D(*cov_mat, nf*nf);
		//std::cout << uniform_dis(_rng) << "\n";
		delete cov_mat;
		delete xi;
		if (l == 0)
		{
			//Compute Pf only
			std::auto_ptr<CS_2D> Modelf(new CS_2D(l + 1, 4, 4, kf));
			Pf = Modelf->compute(l,np);
		}
		else
		{
			//Compute Pf and Pc
			std::auto_ptr<CS_2D> Modelf(new CS_2D(l + 1, 4, 4, kf));
			Pf = Modelf->compute(l, np);
			double *kc = upscale(kf, nf, nf);
			std::auto_ptr<CS_2D> Modelc(new CS_2D(l, 4, 4, kc));
			Pc = Modelc->compute(l, np);
		}

		if (l == 0) dP = Pf;
		else dP = Pf - Pc;

		sums[0] += nf;     // add number of timesteps as cost
		sums[1] += dP;
		sums[2] += dP*dP;
		sums[3] += dP*dP*dP;
		sums[4] += dP*dP*dP*dP;
		sums[5] += Pf;
		sums[6] += Pf*Pf;
		delete kf;
	}
#ifdef DEBUG_Covariance
	for (int i = 0; i < nf*nf; i++)
	{
		ksum[i] /= N;
		xisum[i] /= N;
	}
	print_matrix_2D(ksum, nf, nf);
	print_matrix_2D(xisum, nf, nf);
#endif // DEBUG_Covariance
}




double * construct_cov_mat(int nx, int ny)
{
	// 1*1 physical field, with nx*ny grid
	double dx = 1.0 / nx;
	double dy = 1.0 / ny;

	double* array2D = 0;
	array2D = new double[(nx*ny)*(nx*ny)];
	
	for (int i = 0; i < nx*ny; i++)
	{
		for (int j = 0; j < nx*ny; j++)
		{
			double hi = (i/nx) * dy;
			double hj = (j/nx) * dy;
			double wi = (i%nx) * dx;
			double wj = (j%nx) * dx; //height,width of i and j

			double height = abs(hi-hj);//abs(i-j)/nx * dy;
			double width = abs(wi-wj);//abs(i-j)%nx * dx;
			array2D[i*(nx*ny) + j] = exp(-pow(height,2) - pow(width,2)); //exp(-height^2 - width^2)
		}
	}
	return array2D;
}

void print_matrix_2D(double *array_2d, int nx, int ny)
{
	std::cout << "print 2d array" << std::endl;
	for (int i = 0; i < (nx); i++)
	{
		for (int j = 0; j < (ny); j++)
		{
			std::cout << array_2d[i*ny+j] << std::setw(15);
		}
		std::cout << "\n";
	}
}

void print_matrix_1D(double *array_1d, int n)
{
	std::cout << "print 1d array" << std::endl;
	for (int i = 0; i < n; i++)
		std::cout << array_1d[i] << std::setw(10);
}

double * lowertri_times_vector(double * lowertri, double * vec, int N)
{
	double* result = new double[N]();
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			result[i] += lowertri[i*N + j] * vec[j];
		}
	}
	return result;
}

double * upscale(double * k_fine, int nx, int ny)
{
	double* k_coarse = new double[nx / 2 * ny / 2]();
	for (int i = 0; i < nx / 2; i++)
	{
		for (int j = 0; j < ny / 2; j++)
		{
			k_coarse[i * ny / 2 + j] = 1 * k_fine[(2 * i + 1)*ny + (2 * j + 1)] * k_fine[(2 * i + 1)*ny + (2 * j)] * k_fine[(2 * i)*ny + (2 * j + 1)] * k_fine[(2 * i)*ny + (2 * j)];
			k_coarse[i * ny / 2 + j] = pow(k_coarse[i * ny / 2 + j], 1.0 / 4.0); // 1 / 4
		}
	}
	return k_coarse;
}
