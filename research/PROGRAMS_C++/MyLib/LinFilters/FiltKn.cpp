//---------------------------------------------------------------------------


#pragma hdrstop

#include "FiltKn.h"

#include <stdio.h>
#include <string.h>
#include <vcl.h>
#include <math.h>
#include "LinFilt.h"
#include "Measure.h"
#include "Comp.h"
#include "MatrixProccess.h"
#include "Equations.h"



//---------------------------------------------------------------------------
TFiltKn::TFiltKn():TLinFilt()
{
}



TFiltKn::TFiltKn(const double ValT,   double *parrEstX, double *parrMtrxK):TLinFilt(2,  1,  ValT,  parrEstX, parrMtrxK)
{
}



void TFiltKn::createMtrxCorrU(const TMeasure Meas,const double Val_h,const double Val_Tau2, double *pMtrxCorrU)
   {
	 pMtrxCorrU[0] = Val_h * Val_h * Val_h /3. ;
	 pMtrxCorrU[1] = Val_h * Val_h  /2. ;
	 pMtrxCorrU[2] = pMtrxCorrU[1];
	 pMtrxCorrU[3] = Val_h  ;
	 for (int i =0; i < 4; i++)
	 {
	   pMtrxCorrU[i] *= Val_Tau2;
	 }
   }
 void TFiltKn::createMtrxC(const TMeasure Meas,const double Val_h, double *pMtrxC)
	{
	 pMtrxC[0] = 1.;
	 pMtrxC[1] = 0.;
	}


 void TFiltKn::createMtrxA(const TMeasure Measure, const double Val_h, double *pMtrxA)
{
	pMtrxA[0] = 1.;
	pMtrxA[1] = Val_h;
	pMtrxA[2] = 0.;
	pMtrxA[3] = 1.;
}


// ЭТО НЕ ФИЛЬТР КНЯЗЕВА!!!!!
// аналитический расчет фильтра Калмана в установиви режиме
// Input:
// Val_Tau2 - дисперсия шума в правой части уравн движения (объекта)
// ValDispKsi- дисперсия  шума измерителя
// Output:
// arrK - коррелл матрица
// arrP - коефф усиления
// pcmparrRoots[0], pcmparrRoots[1]- корни характерист уравнения
// arrT[0], arrT[1]- память фильтра по положению и по скорости
// если система неустойчивая, то  arrT[1] = arrT[0] =  -1
void TFiltKn::StabSolutionKalm(const double Val_h,const double Val_Tau2,const double ValDispKsi
  , double *arrK,double * arrP, TComp *pcmparrRoots
	 , double *arrT)
{

	const double Sigm_w = sqrt(Val_Tau2);
	const double Sig_bmo = sqrt(ValDispKsi);
	double c =  Sigm_w * Val_h * Val_h / Sig_bmo ;
	double x = SolvCharactEqKlm(c);
	arrK[0] =  Sig_bmo * Sig_bmo * (x * x -1 )/ x / x;
	arrK[1] = arrK[2] =  Sig_bmo * Sigm_w * Val_h / x ;
	arrK[3] = Sigm_w * Sig_bmo * ( x * x -1 ) /x - Sigm_w * Sigm_w * Val_h * Val_h /2 ;
	arrP [0] =  (x * x -1 )/ x / x;
	arrP [1] =  Sigm_w  * Val_h/ Sig_bmo / x ;
	// расчет корней характ уравнения
	double a1 =1;
	double b1 = 2 - arrP [0] - Val_h * arrP [1] ;
	double c1 = 1- arrP [0] ;

	SolvEq2( a1, b1, c1,pcmparrRoots[0], pcmparrRoots[1] ) ;
	double vmod1 = pcmparrRoots[0].modul();
	double vmod2 = pcmparrRoots[1].modul();
	if ((vmod1 > 1) || (vmod2 > 1))
	{
	  arrT[0] = -1 ;
	  arrT[1] = -1 ;
	  return ;
	}
	// расчет памяти фильтра
	double arrL[4]={0},arrV[4] = {0},arrLamb[4]={0} ;
	arrL[0] =  1- arrP [0] ;
	arrL[1] =  Val_h * (1- arrP [0] );
	arrL[2] =  - arrP [1] ;
	arrL[3] =  1 - Val_h *  arrP [1] ;
	CalcProperVectors2(arrL,arrV , arrLamb);
	arrT[0] = -0. ;
	arrT[1] = -0. ;
	double vmod1_n =  1;
	double vmod2_n =  1;

		double arr_q[2] = {0} ;
		MtrxTranspMultMatrx(arrV,2, 2, arrP,1, arr_q) ;  // arr_q = Ф(0)T * P
		double val1 = 0.05 * (fabs(arrV[0]  * arr_q[0] ) + fabs(arrV[1] *  arr_q[1]) ) ;
		double val2 = 0.05 * (fabs(arrV[2]  * arr_q[0] ) + fabs(arrV[3] *  arr_q[1]) ) ;
   // расчет памяти по положению
	while (1)
	{
		arrT[0] += Val_h;

		 vmod1_n *=vmod1;
		 vmod2_n *=vmod2;

		if ((vmod1_n * fabs(arrV[0] * arr_q[0]) + vmod2_n * fabs(arrV[1] * arr_q[1])  ) <=  val1) break ;

	}
	 // расчет памяти по скорости
	  vmod1_n =  1;
	  vmod2_n =  1;
	while (1)
	{
		arrT[1] += Val_h;

		 vmod1_n *=vmod1;
		 vmod2_n *=vmod2;
		if ((vmod1_n * fabs(arrV[2] * arr_q[0]) + vmod2_n * fabs(arrV[3] * arr_q[1])  ) <=  val2) break ;

	}

}

double  TFiltKn::SolvCharactEqKlm(const double c)
{
	if (c >= 4.)
	{
	   ShowMessage(L"C >= 4") ;
	   return - 1001;
	}
	return 1 + c /4. + sqrt(c /2. + c * c /16.) ;
}


// аналитический расчет фильтра Калмана сопровождения в установиви режиме
// Input:

// pparrWeight - двойной указатель на массив весовой йункции
// lenparrWeight - начальная длина массива весов функции
// Output:
// arrK - коррелл матрица
// arrP - коефф усиления
// rt1, rt2 - корни характерист уравнения
// valT1, valT2 - память фильтра по положению и по скорости
// если система неустойчивая, то  valT1 = valT2 =  -1

void TFiltKn::fncAnalyticSolution(const double Val_h,const double Val_Tau2,const double ValDispKsi
  , double *arrK,double * arrP, TComp *pcmparrRoots
	 , double *arrT)
{
	const double Sigm_w = sqrt(Val_Tau2);
	const double Sig_bmo = sqrt(ValDispKsi);
	double mu =  Sigm_w* Val_h * sqrt(Val_h) / Sig_bmo;
	double x = SolvCharactEqKlm_Real(mu);
	arrK[0] =  Sig_bmo* Sig_bmo* (x * x -1 )/ x / x;
	arrK[1] =  Sig_bmo* Sigm_w* sqrt(Val_h )/ x ;
	arrK[2] =arrK[1] ;
	arrK[3] = Sigm_w* Sig_bmo* ( x * x -1 ) /x / sqrt(Val_h ) - Sigm_w* Sigm_w* Val_h  /2 ;
	arrP [0] =  (x * x -1 )/ x / x;
	arrP [1] =  Sigm_w * sqrt( Val_h ) / Sig_bmo/ x ;
	// расчет корней характ уравнения
	double a1 =1;
	double b1 = 2 - arrP [0] - Val_h * arrP [1] ;
	double c1 = 1- arrP [0] ;
	SolvEq2( a1, b1, c1,pcmparrRoots[0],pcmparrRoots[1] ) ;
	double vmod1 = pcmparrRoots[0].modul();
	double vmod2 = pcmparrRoots[1].modul();
	if ((vmod1 > 1) || (vmod2 > 1))
	{
	  arrT[0] = -1 ;
	  arrT[1] = -1 ;
	  return ;
	}
	// расчет памяти фильтра
	double arrL[4]={0},arrLamb[4]={0}, arrV[4] ={0} ;
	arrL[0] =  1- arrP [0] ;
	arrL[1] =  Val_h * (1- arrP [0] );
	arrL[2] =  - arrP [1] ;
	arrL[3] =  1 - Val_h *  arrP [1] ;
	CalcProperVectors2(arrL,arrV , arrLamb);
	arrT[0] = 0;
	arrT[1] = 0;
		double vmod1_n =  1;
		double vmod2_n =  1;

		double arr_q[2] = {0} ;
		MtrxTranspMultMatrx(arrV,2, 2, arrP,1, arr_q) ;  // arr_q = Ф(0)T * P
		double val1 = 0.05 * (fabs(arrV[0]  * arr_q[0] ) + fabs(arrV[1] *  arr_q[1]) ) ;
		double val2 = 0.05 * (fabs(arrV[2]  * arr_q[0] ) + fabs(arrV[3] *  arr_q[1]) ) ;
   // расчет памяти по положению
	while (1)
	{
		arrT[0] += Val_h;

		 vmod1_n *=vmod1;
		 vmod2_n *=vmod2;

		if ((vmod1_n * fabs(arrV[0] * arr_q[0]) + vmod2_n * fabs(arrV[1] * arr_q[1])  ) <=  val1) break ;

	}
	 // расчет памяти по скорости
	  vmod1_n =  1;
	  vmod2_n =  1;
	while (1)
	{
		arrT[1] += Val_h;

		 vmod1_n *=vmod1;
		 vmod2_n *=vmod2;
		if ((vmod1_n * fabs(arrV[2] * arr_q[0]) + vmod2_n * fabs(arrV[3] * arr_q[1])  ) <=  val2) break ;

	}


}
//решение нелинейного уравнения для реального фильтра калмана
double  TFiltKn::SolvCharactEqKlm_Real(const double mu)
{
	double t0 = mu/2. + sqrt(mu * mu /12. + 4.) ;
	return t0/2. + sqrt(t0*t0/4. -1.) ;
}



#pragma package(smart_init)
