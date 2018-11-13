#include "CS_2D.h"
#include "mlmc.h"
#include "mlmc_test.h"
#include "elliptic_2D.h"
int main()
{
	int M = 2;     // refinement cost factor
	int N0 = 2000;   // initial samples on each level
	int Lmin = 2;   // minimum refinement level
	int Lmax = 10;  // maximum refinement level

	int   N=N0, L=3; //mlmc_test parameters
	float Eps[]={ 0.005, 0.01, 0.02, 0.05, 0.1, 0.0 };
	char  filename[32];
	FILE *fp;

	sprintf(filename, "elliptic_2D.txt");
	//sprintf(filename, "elliptic_%d D.txt",d);
	fp = fopen(filename, "w");
	std::mt19937 rng(1234);
	mlmc_test(elliptic_2D, M, N, L, N0, Eps, Lmin, Lmax, fp, rng);

	fclose(fp);


	/*CS_2D test(22,22,3);

	test.GRID();
	test.BOUND(1);
	test.SOLVE();
	test.BOUND(1);
	test.OUTPUT(1);*/

	return 0;
}