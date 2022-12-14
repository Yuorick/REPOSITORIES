


#pragma hdrstop
#include <vcl.h>
#include "Traject.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "MatrixProccess.h"
#include <time.h>
#include "Gauss.h"


// ??? ??????? ???? ?????
#define DEBUG



//---------------------------------------------------------------------------
TTraject::TTraject()
{
/*
	mTCur = 0.;
	mTBegin = 0;


	mSigW = 33. ;

	memset(marrVectSostGSK, 0, LEN_TARG_VS*sizeof(double)) ;
	memset(marrVectDeviations, 0, LEN_TARG_VS*sizeof(double)) ;
	memset(marrVectSostGSK_Begin, 0, LEN_TARG_VS*sizeof(double)) ;
   */


}

//---------------------------------------------------------------------------


// ??????????? ???????????
 TTraject ::TTraject (const TTraject &R)
 {
	mTCur = R.mTCur ;

	mTBegin = R.mTBegin;
	memcpy(marrSigW,R.marrSigW, 3 * sizeof(double)) ;
	memcpy(marrVectSostGSK,R.marrVectSostGSK, 9 * sizeof(double)) ;
	memcpy(marrVectDeviations,R.marrVectDeviations, 9 * sizeof(double)) ;
	memcpy(marrVectSostGSK_Begin,R.marrVectSostGSK_Begin, 9 * sizeof(double)) ;
 //	memcpy(marrVS_OSKPSb, R.marrVS_OSKPSb, 6 * sizeof(double)) ;
 //	mHomingHead3D = R.mHomingHead3D ;

 }
 // ???????? ????????????
 TTraject &TTraject::operator=(const TTraject  &R)
 {

	mTCur = R.mTCur ;

	mTBegin = R.mTBegin;

	memcpy(marrSigW,R.marrSigW, 3 * sizeof(double)) ;
	memcpy(marrVectSostGSK,R.marrVectSostGSK, 9 * sizeof(double)) ;
	memcpy(marrVectDeviations,R.marrVectDeviations, 9 * sizeof(double)) ;
	memcpy(marrVectSostGSK_Begin,R.marrVectSostGSK_Begin, 9 * sizeof(double)) ;
  //	memcpy(marrVS_OSKPSb, R.marrVS_OSKPSb, 6 * sizeof(double)) ;
//	mHomingHead3D = R.mHomingHead3D ;


	return *this ;
 }

  // ????? ???????????
 TTraject::TTraject (const double TCur, const double 	SigW
				,double *arrVectDeviations,double *arrVectSostGSK
			,double *arrVectSostGSK_Begin )
 {


	mTCur = TCur;
	mTBegin  = TCur;

	marrSigW[0] =  SigW ;
	marrSigW[1] =  SigW ;
	marrSigW[2] =  SigW ;

	memcpy(marrVectSostGSK,arrVectSostGSK, 9 * sizeof(double)) ;
	memcpy(marrVectDeviations,arrVectDeviations, 9 * sizeof(double)) ;
	memcpy(marrVectSostGSK_Begin,arrVectSostGSK_Begin, 9 * sizeof(double)) ;

 }

  // ????? ???????????
 TTraject::TTraject (const double TBegin,  const double SigW ,double *arrVectSostGSK_Begin )
 {


	mTCur = TBegin;
	mTBegin  = TBegin;
	marrSigW[0] =  SigW ;
	marrSigW[1] =  SigW ;
	marrSigW[2] =  SigW ;


	memcpy(marrVectSostGSK,arrVectSostGSK_Begin, 9 * sizeof(double)) ;
	memcpy(marrVectSostGSK_Begin,arrVectSostGSK_Begin, 9 * sizeof(double)) ;
	memset(marrVectDeviations,0, 9 * sizeof(double)) ;


 }

  // ????? ???????????
 TTraject::TTraject (const double TCur,  const double 	SigW,	TInitTargData InitData )
 {


	mTCur = TCur;
	mTBegin  = TCur;
	marrSigW[0] =  SigW ;
	marrSigW[1] =  SigW ;
	marrSigW[2] =  SigW ;
//	memcpy(marrParams,arrParams, 25 * sizeof(double)) ;

if (fabs(InitData.mH /InitData.mR) > 0.9999999)
{
	ShowMessage(L"ERROR");
	int iii = 0;
	InitData.mH  =0.;

}


	 double valeps = asin(InitData.mH /InitData.mR) ;
   double  arrVectSostGSK[9] = {0} ;
   arrVectSostGSK[0] = InitData.mR * cos(valeps) * sin(InitData.mBearing  ) ;
   arrVectSostGSK[1] = InitData.mR * cos(valeps) * cos(InitData.mBearing ) ;
   arrVectSostGSK[2] = InitData.mH ;
   arrVectSostGSK[3] = InitData.mV * sin(InitData.mTargZenitAng) * sin(InitData.mTargCourse) ;
   arrVectSostGSK[4] = InitData.mV * sin(InitData.mTargZenitAng) * cos(InitData.mTargCourse) ;
   arrVectSostGSK[5] = InitData.mV * cos(InitData.mTargZenitAng) ;

	 double arrVectDeviations [9] ={0.};

	memcpy(marrVectSostGSK,arrVectSostGSK, 9 * sizeof(double)) ;
	memcpy(marrVectDeviations,arrVectDeviations, 9 * sizeof(double)) ;
	memcpy(marrVectSostGSK_Begin,arrVectSostGSK, 9 * sizeof(double)) ;

 }

 /////


  // ????? ???????????
 TTraject::TTraject (const double TCur,  double 	*arrSigW,	TInitTargData InitData )
 {


	mTCur = TCur;
	mTBegin  = TCur;
	marrSigW[0] =  arrSigW[0] ;
	marrSigW[1] =  arrSigW[1] ;
	marrSigW[2] =  arrSigW[2] ;
//	memcpy(marrParams,arrParams, 25 * sizeof(double)) ;

if (fabs(InitData.mH /InitData.mR) > 0.99999999)
{
	ShowMessage(L"ERROR");
	int iii = 0;
	InitData.mH  =0.;

}

	 double valeps = asin(InitData.mH /InitData.mR) ;
   double  arrVectSostGSK[9] = {0} ;
   arrVectSostGSK[0] = InitData.mR * cos(valeps) * sin(InitData.mBearing  ) ;
   arrVectSostGSK[1] = InitData.mR * cos(valeps) * cos(InitData.mBearing ) ;
   arrVectSostGSK[2] = InitData.mH ;
	 arrVectSostGSK[3] = InitData.mV * sin(InitData.mTargZenitAng) * sin(InitData.mTargCourse) ;
   arrVectSostGSK[4] = InitData.mV * sin(InitData.mTargZenitAng) * cos(InitData.mTargCourse) ;
   arrVectSostGSK[5] = InitData.mV * cos(InitData.mTargZenitAng) ;

	 double arrVectDeviations [9] ={0.};

	memcpy(marrVectSostGSK,arrVectSostGSK, 9 * sizeof(double)) ;
	memcpy(marrVectDeviations,arrVectDeviations, 9 * sizeof(double)) ;
	memcpy(marrVectSostGSK_Begin,arrVectSostGSK, 9 * sizeof(double)) ;

 }

 /////


// ?????? ???????????? ????????
void TTraject::recalcTrajPoint(const double tNext)
{
	double h = tNext - mTCur;
  if ( h < 0)
  {
	 ShowMessage(L"Error.TTraject::recalcTrajPoint_0. h<0") ;
	 return ;
  }
  double arrTemp[9] = {0} ;
	for (int i = 0; i < 3; i++) arrTemp[6 + i] = getGauss(0,marrSigW[i]) ;
  for (int i = 0; i < 3; i++)
  {
	 arrTemp[3 + i] = marrVectSostGSK  [3 + i]+ h *arrTemp[6 + i];
	 arrTemp[    i] = marrVectSostGSK  [ i] + marrVectSostGSK  [3 + i]* h
				+ arrTemp[6 + i]*h*h/2. ;;
  }
	memcpy( marrVectSostGSK, arrTemp, 9 * sizeof(double)) ;
  mTCur = tNext;


}

// ???????? ????????????? ??????? ????????? ?? ?????  VAlTExtr ???????
 void TTraject::extrapolateTargVS(const double VAlTExtr, double *arrTargExtrapVS)
 {
	 memcpy(arrTargExtrapVS, marrVectSostGSK, 9 * sizeof(double));
	 double arrT[3] = {0.};
	 MatrxMultScalar(&marrVectSostGSK[3], 1, 3, VAlTExtr, arrT);
	 MtrxSumMatrx(arrT, marrVectSostGSK,3, 1, arrTargExtrapVS) ;
 }


 // ???????? ????????????? ?????????? ??????? ????????? ?? ?????  VAlTExtr ???????
 void TTraject::extrapolateTarg_BeginVS(const double VAlTExtr, double *arrTargExtrapVS)
 {
	 memcpy(arrTargExtrapVS, marrVectSostGSK_Begin, 9 * sizeof(double));
	 double arrT[3] = {0.};
	 MatrxMultScalar(&marrVectSostGSK_Begin[3], 1, 3, VAlTExtr, arrT);
	 MtrxSumMatrx(arrT, marrVectSostGSK_Begin,3, 1, arrTargExtrapVS) ;
 }


 /*

// ????????? ?????????? ????? ? ??????????? ?? ??????? [0;1]  ??????????????
double TTraject::getRand01( )
{

 // return ((double)rand())/ ((double)RAND_MAX);
  return (Rand_())/ ((double)RAND_MAX);
}
//
double TTraject::Rand_()
{
  #ifndef DEBUG
  static int rand_init = 0;
  if( 0  == rand_init )
  {
	  rand_init = time( NULL );
	  srand((unsigned) time(NULL));
  }

  #endif

 // return (double)rand();
 // double breturn =  (double)rand(); // ??? ???????
  return (double)rand();                  // ??? ???????

}
 */
#pragma package(smart_init)
