#pragma once
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
class CS_2D
{
	int NX, NY, NGIT;
	double XL=1.0, YL=1.0, RESMAX;
	double** X, **XCV, **XU, **XDIF; 
	double** Y, **YCV, **YU, **YDIF;
	double*** F, ***F0, ***FF;
	double*** AP, ***AE, ***AW, ***AN, ***AS, ***CON;
	double*** GAM, ***ERROR, ***RES;
	double T_INTEVAL, T_END, T_START, SORG=1.0, SORR=1.0, SORP=1.0;
	int I, J, L1, L2, L3, M1, M2, M3;
	int KGR, NICV, NJCV, NITER, NMETHOD;
	std::string filename;

	

public:
	CS_2D(int _NGIT, int _NICV, int _NJCV, double* _Gamma);
	~CS_2D();
	void SETIND(int _KGR1);
	void GRID();
	void OUTPUT(int _KGR1, int _l, int _np);
	void OUTPUT_Gamma(int _KGR1, int _l, int _np);

	void BOUND(int _KGR1);
	void SOLVE();
	void COEF(int _KGR1);
	void GAMSOR(int _KGR1);
	void RCN(int _KGR1);
	void PLG(int _KGR1);
	void GAUSSSEIDEL1(int _KGR1);
	double compute(int _l, int _np);
};

