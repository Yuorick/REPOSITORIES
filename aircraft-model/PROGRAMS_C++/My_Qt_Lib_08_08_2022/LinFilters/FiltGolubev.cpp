//---------------------------------------------------------------------------


#pragma hdrstop

#include "FiltGolubev.h"
#include <string.h>
#include <math.h>
#include <vcl.h>
#include "Comp.h"
#include "Equations.h"
#include "MatrixProccess.h"

//---------------------------------------------------------------------------
 TFiltGolubev::~TFiltGolubev()
{
	if(mparrEstX) delete mparrEstX ;
	mparrEstX = NULL ;
	if(mparrMtrxK) delete mparrMtrxK ;
	mparrMtrxK = NULL ;
}
TFiltGolubev::TFiltGolubev()
{
	// ����������� ������� ��������� �������
	mDimX = 2;

	// ����������� ������� ����������
	mDimY = 1;
	// ������������ ����������
//  1.����� �������� ������ ������� ���������
	 mTf = 0.;
//  2.������ ������� ���������   �� ������ mTf
	mparrEstX = new double [mDimX];

 // 3.�������������� ������� ������ ����������
	mparrMtrxK = new double [mDimX * mDimX];


}


TFiltGolubev::TFiltGolubev(const double ValT,   double *parrEstX, double *parrMtrxK)
{
	// ����������� ������� ��������� �������
  	mDimX = 2;

	// ����������� ������� ����������
	mDimY = 1;
	mTf = ValT;


  //	2.������ ������� ���������   �� ������ mTf
	mparrEstX = new double [mDimX];
	memcpy(mparrEstX, parrEstX, mDimX * sizeof(double));

 // 3.�������������� ������� ������ ����������
	mparrMtrxK = new double [mDimX * mDimX];
	memcpy(mparrMtrxK, parrMtrxK, mDimX *mDimX * sizeof(double));

}

void TFiltGolubev::createMtrxCorrU(const TMeasure Meas,const double Val_h,const double Val_Tau2, double *pMtrxCorrU)
   {
	 double sq = sqrt(mparrMtrxK[3] - mparrMtrxK[1]* mparrMtrxK[1]/mparrMtrxK[0]);
	 double sigm_w = sqrt(Val_Tau2);
	 double vala =  sqrt(Val_Tau2) * sq ;
	 pMtrxCorrU[0] = Val_Tau2 * Val_h * Val_h * Val_h * Val_h  /4. + sigm_w * Val_h * Val_h * Val_h *sq;
	 pMtrxCorrU[1] = Val_h * Val_h* Val_h* Val_Tau2 /2. + 3.* sigm_w  * Val_h* Val_h *sq /2. ;
	 pMtrxCorrU[2] = pMtrxCorrU[1];
	 pMtrxCorrU[3] = 2. *sigm_w * sq  * Val_h  + Val_Tau2 * Val_h * Val_h;

   }
 void TFiltGolubev::createMtrxC(const TMeasure Measure,const double Val_h, double *pMtrxC)
	{
	 double valAlf = fncMinDouble(1., sqrt(Measure.mparrMMO_K[0]/ mparrMtrxK[0]));
	 pMtrxC[0] = 1. -valAlf;
	 pMtrxC[1] = 0.;
	}


 void TFiltGolubev::createMtrxA(const TMeasure Measure, const double Val_h, double *pMtrxA)
{
	pMtrxA[0] = 1.;
	pMtrxA[1] = Val_h;
	pMtrxA[2] = 0.;
	pMtrxA[3] = 1.;
}

double  TFiltGolubev::fncMinDouble(double a, double b)
{
	if (a < b)
	{
	  return a;
	}
	return b;
}

int TFiltGolubev::fncStep( TMeasure Measure, const double Val_Tau2)
{
	const double Val_h = Measure.mTYZv - mTf;
	if (Val_h <= 0.)
	{
	return -1;
	}
	double sigm_mmo = sqrt(Measure.mparrMMO_K[0]);
	double sigm_bmo = sqrt(Measure.mparrBMO_K[0]);
	double sigm_w = sqrt( Val_Tau2);
	double arrP[2] ={0.};
		// �������������
	double arrKExtr[4] ;
	double sq = sqrt(mparrMtrxK [3] - mparrMtrxK [1] * mparrMtrxK [1] /mparrMtrxK [0]) ;
	arrKExtr[0] =  mparrMtrxK [0] + 2. * Val_h * mparrMtrxK[1] + Val_h * Val_h * mparrMtrxK [3]
	   +  sigm_w * sigm_w * Val_h* Val_h * Val_h* Val_h/4. + sigm_w * Val_h * Val_h * Val_h * sq;
	arrKExtr[1] =  mparrMtrxK [1] + Val_h *mparrMtrxK [3] + 3. * sigm_w * Val_h * Val_h * sq/2.
		  + sigm_w * sigm_w * Val_h * Val_h * Val_h /2.;
	arrKExtr[2] = arrKExtr[1] ;
	arrKExtr[3] =  mparrMtrxK[3] + 2. * sigm_w * Val_h * sq + sigm_w * sigm_w * Val_h * Val_h ;
	// ����� ���
	double alf =  sigm_mmo/ sqrt( arrKExtr[0]) ;
	if (alf > 1) alf = 1 ;
	double valF = ( 1.- alf)/( ( 1.- alf)* ( 1.- alf) * arrKExtr[0] + sigm_bmo * sigm_bmo) ;

	arrP[0] = arrKExtr[0]* valF ;
	arrP[1] =  valF * arrKExtr[1] ;
	if (arrP[0] > 1.)
	{
	  arrP[0] = 1. ;
	  arrP[1] = arrKExtr[1] / arrKExtr[0];
	  alf = ( sigm_mmo * sigm_mmo + sigm_bmo * sigm_bmo) / arrKExtr[0] ;
	}

	double valP = 1. - (1.- alf) * arrP[0] ;
	mparrMtrxK[ 0] = valP * arrKExtr[0] ;
	mparrMtrxK[ 1] = valP * arrKExtr[1] ;
	mparrMtrxK[ 2] = mparrMtrxK[ 1] ;
	mparrMtrxK[ 3] = arrKExtr[3] -  ( 1. - alf) * arrP[1] * arrKExtr[1] ;
	mTf =  Measure.mTYZv;
}


void TFiltGolubev::fncAnalyticSolution(const double Val_h,const double Val_Tau2,const double ValDispKsi
  , double *arrK,double * arrP, TComp *pcmparrRoots
	 , double *arrT)
{

}
#pragma package(smart_init)
