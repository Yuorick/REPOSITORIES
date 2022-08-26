//---------------------------------------------------------------------------


#pragma hdrstop

#include "FiltAlfa.h"
#include <string.h>
#include <math.h>
#include <vcl.h>
#include "Comp.h"
#include "Equations.h"
#include "MatrixProccess.h"

//---------------------------------------------------------------------------
 TFiltAlfa::~TFiltAlfa()
{
	if(mparrEstX) delete mparrEstX ;
	mparrEstX = NULL ;
	if(mparrMtrxK) delete mparrMtrxK ;
	mparrMtrxK = NULL ;
}
TFiltAlfa::TFiltAlfa()
{
	// размерность вектора состояния объекта
	mDimX = 2;

	// размерность вектора наблюдений
	mDimY = 1;
	// Динамическая информация
//  1.Время привязки оценок вектора состояния
	 mTf = 0.;
//  2.Оценка вектора состояния   на момент mTf
	mparrEstX = new double [mDimX];

 // 3.корреляционная матрица ошибок оценивания
	mparrMtrxK = new double [mDimX * mDimX];


}


TFiltAlfa::TFiltAlfa(const double ValT,   double *parrEstX, double *parrMtrxK)
{
	// размерность вектора состояния объекта
  	mDimX = 2;

	// размерность вектора наблюдений
	mDimY = 1;
	mTf = ValT;


  //	2.Оценка вектора состояния   на момент mTf
	mparrEstX = new double [mDimX];
	memcpy(mparrEstX, parrEstX, mDimX * sizeof(double));

 // 3.корреляционная матрица ошибок оценивания
	mparrMtrxK = new double [mDimX * mDimX];
	memcpy(mparrMtrxK, parrMtrxK, mDimX *mDimX * sizeof(double));

}

void TFiltAlfa::createMtrxCorrU(const TMeasure Meas,const double Val_h,const double Val_Tau2, double *pMtrxCorrU)
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
 void TFiltAlfa::createMtrxC(const TMeasure Measure,const double Val_h, double *pMtrxC)
	{

	 pMtrxC[0] = 1. ;
	 pMtrxC[1] = 0.;
	}


 void TFiltAlfa::createMtrxA(const TMeasure Measure, const double Val_h, double *pMtrxA)
{
	pMtrxA[0] = 1.;
	pMtrxA[1] = Val_h;
	pMtrxA[2] = 0.;
	pMtrxA[3] = 1.;
}



int TFiltAlfa::fncStep( TMeasure Measure, const double Val_Tau2)
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
		// Экстраполяция
	 double arrKExtr[4] = {0} , arrKt1[4] = {0},arrF[4] = {0} ;

	 double arrKtemp1[4], arrA[4] = {0.};;

	 createMtrxCorrU(Measure, Val_h, Val_Tau2, arrF) ;
	 createMtrxA( Measure,  Val_h, arrA) ;

  /*arrF[0] = msigm_w* msigm_w* mh * mh * mh /3 ;
  arrF [ 1] =  msigm_w* msigm_w* mh * mh  /2 ;
  arrF [ 2] =arrF [ 1] ;
  arrF [ 3] = msigm_w* msigm_w* mh  ;*/

	MtrxMultMatrx(arrA ,2, 2, mparrMtrxK,2, arrKt1) ;
	MtrxMultMatrxTransp(arrKt1,2, 2,arrA ,2, arrKtemp1) ;
	MtrxSumMatrx(arrKtemp1, arrF,2,2, arrKExtr) ;
	// Выбор ММО
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

void TFiltAlfa::fncAnalyticSolution(const double Val_h,const double Val_Tau2,const double ValDispKsi
  , double *arrK,double * arrP, TComp *pcmparrRoots
	 , double *arrT)
{

}
#pragma package(smart_init)
