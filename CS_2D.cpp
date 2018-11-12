#include "CS_2D.h"





CS_2D::CS_2D(int _NGIT, int _NICV, int _NJCV, double* _Gamma)
{

	NGIT = _NGIT;
	NICV = _NICV;
	NJCV = _NJCV;
	NX = NICV * int(pow(2, NGIT - 1)) + 2;
	NY = NJCV * int(pow(2, NGIT - 1)) + 2;
	// Allocate Memory (NGIT+1) * (NX+1)

	X = new double*[NGIT+1];
	XCV = new double*[NGIT + 1];
	XU = new double*[NGIT + 1];
	XDIF = new double*[NGIT + 1];
	for (int i = 0; i < NGIT + 1; i++)
	{
		X[i] = new double[NX + 1];
		XCV[i] = new double[NX + 1];
		XU[i] = new double[NX + 1];
		XDIF[i] = new double[NX + 1];
	}

	Y = new double*[NGIT + 1];
	YCV = new double*[NGIT + 1];
	YU = new double*[NGIT + 1];
	YDIF = new double*[NGIT + 1];
	for (int i = 0; i < NGIT + 1; i++)
	{
		Y[i] = new double[NY + 1];
		YCV[i] = new double[NY + 1];
		YU[i] = new double[NY + 1];
		YDIF[i] = new double[NY + 1];
	}
		
	//Allocate 3D array (NGIT+1) * (NX+1) * (NY+1)
	F = new double**[NGIT + 1];
	F0 = new double**[NGIT + 1];
	FF = new double**[NGIT + 1];
	AP = new double**[NGIT + 1];
	AE = new double**[NGIT + 1];
	AW = new double**[NGIT + 1];
	AN = new double**[NGIT + 1];
	AS = new double**[NGIT + 1];
	CON = new double**[NGIT + 1];
	GAM = new double**[NGIT + 1];
	ERROR = new double**[NGIT + 1];
	RES = new double**[NGIT + 1];
	for (int i = 0; i < NGIT+1; i++)
	{
		F[i] = new double*[NX + 1];
		F0[i] = new double*[NX + 1];
		FF[i] = new double*[NX + 1];
		AP[i] = new double*[NX + 1];
		AE[i] = new double*[NX + 1];
		AW[i] = new double*[NX + 1];
		AN[i] = new double*[NX + 1];
		AS[i] = new double*[NX + 1];
		CON[i] = new double*[NX + 1];
		GAM[i] = new double*[NX + 1];
		ERROR[i] = new double*[NX + 1];
		RES[i] = new double*[NX + 1];
		for (int j = 0; j < NX + 1; j++)
		{
			F[i][j] = new double[NY + 1]();
			F0[i][j] = new double[NY + 1];
			FF[i][j] = new double[NY + 1];
			AP[i][j] = new double[NY + 1];
			AE[i][j] = new double[NY + 1];
			AW[i][j] = new double[NY + 1];
			AN[i][j] = new double[NY + 1];
			AS[i][j] = new double[NY + 1];
			CON[i][j] = new double[NY + 1];
			GAM[i][j] = new double[NY + 1];
			ERROR[i][j] = new double[NY + 1];
			RES[i][j] = new double[NY + 1];
		}
	}
	for (int i = 0; i < NX - 2; i++)
	{
		for (int j = 0; j < NY - 2; j++)
		{
			GAM[1][i + 2][j + 2] = _Gamma[i*(NY - 2) + j];
			if (j == 0) GAM[1][i + 2][j + 1] = GAM[1][i + 2][j + 2];
			if (j == NY - 1) GAM[1][i + 2][j + 1] = GAM[1][i + 2][j];
		}
	}
	
	for (int j = 0; j < NY - 2; j++)
	{
		GAM[1][1][j + 2] = GAM[1][2][j + 2];
		GAM[1][NX][j + 2] = GAM[1][NX-1][j + 2];
	}

	// Construct Gamma for other levels (NGIT)
	// Upscaling
	for (int N = 2; N <= NGIT; N++)
	{
		int ix_l = NICV * pow(2, NGIT - N);
		int iy_l = NJCV * pow(2, NGIT - N);
		for (int i = 0; i < ix_l; i++)
		{
			for (int j = 0; j < iy_l; j++)
			{
				GAM[N][i+2][j+2] = 1 * GAM[N-1][2*i +1 +2][2*j+1 +2] * GAM[N - 1][2 * i +  2][2 * j + 1 + 2] * GAM[N - 1][2 * i + 1 + 2][2 * j + 2] * GAM[N - 1][2 * i + 2][2 * j + 2];
				GAM[N][i + 2][j + 2] = pow(GAM[N][i + 2][j + 2], 1.0 / 4.0); //1.0/4.0	

				if (j == 0) GAM[1][i + 2][j + 1] = GAM[1][i + 2][j + 2];
				if (j == iy_l - 1) GAM[1][i + 2][j + 3] = GAM[1][i + 2][j + 2];
			}
		}
		for (int j = 0; j < NY - 2; j++)
		{
			GAM[1][1][j + 2] = GAM[1][2][j + 2];
			GAM[1][ix_l+2][j + 2] = GAM[1][ix_l + 1][j + 2];
		}
	}
}

CS_2D::~CS_2D()
{
	
}

void CS_2D::SETIND(int _KGR1)
{
	L1 = NICV*int(pow(2, (NGIT - _KGR1))) + 2;
	L2 = L1 - 1;
	L3 = L2 - 1;
	M1 = NJCV*int(pow(2, (NGIT - _KGR1))) + 2;
	M2 = M1 - 1;
	M3 = M2 - 1;
	KGR = _KGR1;
}

void CS_2D::GRID()
{
	double TEMP;
	int NN;
	
	for (NN = 1; NN < NGIT+1; NN++)
	{
		TEMP = 0.0;
		//std::cout << "NN = " << NN;  //Construct grid
		SETIND(NN);
		XCV[KGR][1] = 0.0;
		for (I = 2; I < L2 + 1; I++)
		{
			XCV[KGR][I] = 1.0;
			TEMP = TEMP + XCV[KGR][I];
		}
		XU[KGR][2] = 0.0;
		for (I = 2;  I<L2+1;  I++)
		{
			XCV[KGR][I] = XCV[KGR][I] * XL / TEMP;
			XU[KGR][I + 1] = XU[KGR][I] + XCV[KGR][I];
		}
		XCV[KGR][L1] = 0.0;
		XU[KGR][L1] = XL;
		X[KGR][1] = 0.0;
		for (I = 2; I < L1 + 1; I++)
		{
			XDIF[KGR][I] = (XCV[KGR][I - 1] + XCV[KGR][I]) / 2.0;
			X[KGR][I] = X[KGR][I - 1] + XDIF[KGR][I];
		}
		TEMP = 0.0;
		YCV[KGR][1] = 0.0;
		for (J = 2; J < M2 + 1; J++)
		{
			YCV[KGR][J] = 1.0;
			TEMP = TEMP + YCV[KGR][J];
		}
		YU[KGR][2] = 0.0;
		for (J = 2; J < M2 + 1; J++)
		{
			YCV[KGR][J] = YCV[KGR][J] * YL / TEMP;
			YU[KGR][J + 1] = YU[KGR][J] + YCV[KGR][J];
		}
		YCV[KGR][M1] = 0.0;
		YU[KGR][M1] = YL;
		Y[KGR][1] = 0;
		for (J = 2; J < M1 + 1; J++)
		{
			YDIF[KGR][J] = (YCV[KGR][J - 1] + YCV[KGR][J]) / 2.0;
			Y[KGR][J] = Y[KGR][J - 1] + YDIF[KGR][J];
		}

	}
}

void CS_2D::OUTPUT(int _KGR1, int _l, int _np)
{
	SETIND(_KGR1);
	std::ofstream output;
	char filename[50];
	sprintf(filename, "res_level_%d_n_%d",_KGR1, _np);
	std::string sfilename = filename;
	output.open(sfilename+".dat");
	output << "TITLE = \"FAI\" \n";
	output << "VARIABLE = \"X\", \"Y\", \"FAI\" \n";
	output << "ZONE T = \"T\" " << "I=" << L1 <<std::setw(5) << "J=" << M1 <<std::setw(5)<< "  C=BLACK\n";
	for (I = 1; I < L1+1; I++)
	{
		for (J = 1; J < M1+1; J++)
		{
			output << "\n" << X[_KGR1][I] <<std::setw(15) << Y[_KGR1][J] << std::setw(15) << F[_KGR1][I][J];
		}
	}
}

void CS_2D::OUTPUT_Gamma(int _KGR1, int _l, int _np)
{
	SETIND(_KGR1);
	std::ofstream output;
	char filename[50];
	sprintf(filename, "gam_level_%d_n_%d", _KGR1, _np);
	std::string sfilename = filename;
	output.open(sfilename + ".dat");
	output << "TITLE = \"FAI\" \n";
	output << "VARIABLE = \"X\", \"Y\", \"FAI\" \n";
	output << "ZONE T = \"T\" " << "I=" << L1 << std::setw(5) << "J=" << M1 << std::setw(5) << "  C=BLACK\n";
	for (I = 1; I < L1 + 1; I++)
	{
		for (J = 1; J < M1 + 1; J++)
		{
			output << "\n" << X[_KGR1][I] << std::setw(15) << Y[_KGR1][J] << std::setw(15) << GAM[_KGR1][I][J];
		}
	}
}

void CS_2D::BOUND(int _KGR1)
{
	SETIND(_KGR1);
	if (KGR == 1)
	{
		for (J = 1; J < M1 + 1 ; J++)
		{
			F[KGR][1][J] = 100.0;
			F[KGR][L1][J] = 0.0;
		}
		for (I = 1; I < L1 + 1; I++)
		{
			F[KGR][I][1] = 10.0;
			F[KGR][I][M1] = 50.0;
		}
	}


	if (KGR != 1)
	{
		for (J = 1; J < M1 + 1; J++)
		{
			F[KGR][1][J] = 0.0;
			F[KGR][L1][J] = 0.0;
		}
		for (I = 1; I < L1 + 1; I++)
		{
			F[KGR][I][1] = 0.0;
			F[KGR][I][M1] = 0.0;
		}
	}
}

void CS_2D::GAMSOR(int _KGR1)
{
	SETIND(_KGR1);
	// Construct Gamma
	/*for (I = 1; I < L1 + 1; I++)
		for (J = 1; J < M1 + 1; J++)
			GAM[KGR][I][J] = 5.0;*/
	//Upscaling



	if (KGR == 1)
	{
		for (I = 2; I < L2 + 1; I++)
			for (J = 2; J < M2 + 1; J++)
				CON[KGR][I][J] = 20000.0*XCV[KGR][I] * YCV[KGR][J];

		for (J = 2; J < M2+1; J++)
		{
			I = 2;
			CON[KGR][I][J] = CON[KGR][I][J] + GAM[KGR][I - 1][J] / XDIF[KGR][I] * YCV[KGR][J] * F[KGR][I - 1][J];
			I = L2;
			CON[KGR][I][J] = CON[KGR][I][J] + GAM[KGR][I + 1][J] / XDIF[KGR][I+1] * YCV[KGR][J] * F[KGR][I + 1][J];
		}

		for (I = 2; I < L2+1; I++)
		{
			J = 2;
			CON[KGR][I][J] = CON[KGR][I][J] + GAM[KGR][I][J - 1] / YDIF[KGR][J] * XCV[KGR][I] * F[KGR][I][J - 1];
			J = M2;
			CON[KGR][I][J] = CON[KGR][I][J] + GAM[KGR][I][J + 1] / YDIF[KGR][J+1] * XCV[KGR][I] * F[KGR][I][J + 1];
		}
	}
}

void CS_2D::COEF(int _KGR1)
{
	SETIND(_KGR1);
	GAMSOR(KGR);
	BOUND(KGR);
	for (I = 2; I < L2+1; I++)
	{
		for (J = 2; J < M2 + 1; J++)
		{
			AW[KGR][I][J] = GAM[KGR][I - 1][J] / XDIF[KGR][I] * YCV[KGR][J];
			AE[KGR][I][J] = GAM[KGR][I + 1][J] / XDIF[KGR][I + 1] * YCV[KGR][J];
			AS[KGR][I][J] = GAM[KGR][I][J - 1] / YDIF[KGR][J] * XCV[KGR][I];
			AN[KGR][I][J] = GAM[KGR][I][J + 1] / YDIF[KGR][J + 1] * XCV[KGR][I];
			AP[KGR][I][J] = AW[KGR][I][J] + AE[KGR][I][J] + AS[KGR][I][J] + AN[KGR][I][J];
		}
	}

	if (KGR == 1)
	{
		for (J = 2; J < M2 + 1; J++)
		{
			I = 2;
			AW[KGR][I][J] = 0.0;
			I = L2;
			AE[KGR][I][J] = 0.0;
		}
		for (I = 2; I < L2 + 1; I++)
		{
			J = 2;
			AS[KGR][I][J] = 0.0;
			J = M2;
			AN[KGR][I][J] = 0.0;
		}
	}
} 

void CS_2D::GAUSSSEIDEL1(int _KGR1)
{	
	SETIND(_KGR1);
	if (KGR == 1) RESMAX = 0.0;
	for (int N = 1; N < NITER + 1; N++)
	{
		for (I = 2; I < L2 + 1; I++)
		{
			for (J = 2; J < M2 + 1; J++)
			{
				F[KGR][I][J] = SORG*((AW[KGR][I][J] * F[KGR][I - 1][J] + AE[KGR][I][J] * F[KGR][I + 1][J] + AS[KGR][I][J] * F[KGR][I][J - 1] + AN[KGR][I][J] * F[KGR][I][J + 1] + CON[KGR][I][J]) / AP[KGR][I][J]) + (1.0 - SORG) * F[KGR][I][J];
			}
		}
	}

	

	for (I = 2; I < L2 + 1; I++)
	{
		for (J = 2; J < M2 + 1; J++)
		{
			RES[KGR][I][J] = AE[KGR][I][J] * F[KGR][I + 1][J] + AW[KGR][I][J] * F[KGR][I - 1][J] + AN[KGR][I][J] * F[KGR][I][J + 1] + AS[KGR][I][J] * F[KGR][I][J - 1] + CON[KGR][I][J] - AP[KGR][I][J] * F[KGR][I][J];
			if (KGR == 1) RESMAX = RESMAX + pow(RES[KGR][I][J], 2);
		}
	}
	RESMAX = sqrt(RESMAX / float(L2 - 1) / float(M2 - 1));
}

double CS_2D::compute(int _l, int _np)
{
	GRID();
	BOUND(1);
	SOLVE();
	BOUND(1);
	//OUTPUT(1,_l,_np);
	//OUTPUT_Gamma(1, _l, _np);
	return (F[1][NX/2][NY/2]+ F[1][NX / 2 + 1][NY / 2]+ F[1][NX / 2][NY / 2 + 1]+ F[1][NX / 2 +1 ][NY / 2 + 1])/4.0;
}

void CS_2D::SOLVE()
{
	// TODO : measure CPU time

	for (int NVC = 1; NVC < 10000; NVC++)
	{
		NITER = 5;
		COEF(1);
		GAUSSSEIDEL1(1);
		for (int N = 2; N < NGIT + 1; N++)
		{
			NITER = 10;
			COEF(N);
			RCN(N);
			GAUSSSEIDEL1(N);
		}
		for (int N = NGIT; N >= 2; N--)
		{
			NITER = 2;
			PLG(N);
			GAUSSSEIDEL1(N - 1);
		}
		//std::cout << "\n" << NVC << " RESMAX= " << RESMAX;
		if (RESMAX < 1.0e-10) break;
		memset(CON, 0, sizeof(CON));
	}
}

void CS_2D::RCN(int _KGR1) //restriction
{
	SETIND(_KGR1);
	for (I = 2; I < L2 + 1; I++)
	{
		for (J = 2; J < M2 + 1; J++)
		{
			RES[KGR][I][J] = SORR * (RES[KGR - 1][2 * I - 1][2 * J - 2] + RES[KGR - 1][2 * I - 2][2 * J - 2] + RES[KGR - 1][2 * I - 1][2 * J - 1] + RES[KGR - 1][2 * I - 2][2 * J - 1]);
			if (KGR != 1) CON[KGR][I][J] = RES[KGR][I][J];
		}
	}
}

void CS_2D::PLG(int _KGR1)
{
	SETIND(_KGR1);
	for (J = 2; J < M2 + 1; J++)
	{
		for (I = 2; I < L2 + 1; I++)
		{
			F[KGR - 1][2 * I - 2][2 * J - 2] = F[KGR - 1][2 * I - 2][2 * J - 2] + SORP*F[KGR][I][J];
			F[KGR - 1][2 * I - 1][2 * J - 2] = F[KGR - 1][2 * I - 1][2 * J - 2] + SORP*F[KGR][I][J];
			F[KGR - 1][2 * I - 2][2 * J - 1] = F[KGR - 1][2 * I - 2][2 * J - 1] + SORP*F[KGR][I][J];
			F[KGR - 1][2 * I - 1][2 * J - 1] = F[KGR - 1][2 * I - 1][2 * J - 1] + SORP*F[KGR][I][J];
		}
	}
}