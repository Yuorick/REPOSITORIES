//---------------------------------------------------------------------------


#pragma hdrstop

#include "SimpleSumDgrm.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "YrWriteShapeFile.h"
#include "URPolyLine.h"
#include "UrPointXY.h"
#include "DiagrSinX.h"
#include "MatrixProccess.h"
#include "Comp.h"

_fastcall  TSimpleSumDgrm ::TSimpleSumDgrm()
{

 // ?????????? ????? ???????? ?????????
  mKSum = 1.;
 //???? ???????????? ????????? ?????????
  mScnSum = -10. * M_PI/ 3000.;
 // ??????? ?????????? ????????? ?????????
   mTension = 127.;

}

// ??????????? ???????????
__fastcall  TSimpleSumDgrm::TSimpleSumDgrm (const TSimpleSumDgrm &R)
 {

	// ?????????? ????? ???????? ?????????
	mKSum = R.mKSum;
	//???? ???????????? ????????? ?????????
	mScnSum = R.mScnSum;
	// ??????? ?????????? ????????? ?????????
	mTension = R.mTension;
 }
// ???????? ????????????
  TSimpleSumDgrm TSimpleSumDgrm::operator=(TSimpleSumDgrm  R)
 {

	// ?????????? ????? ???????? ?????????
	mKSum = R.mKSum;
	//???? ???????????? ????????? ?????????
	mScnSum = R.mScnSum;
	// ??????? ?????????? ????????? ?????????
	mTension = R.mTension;
	 return *this ;
 }

 // ????? ??????
 __fastcall TSimpleSumDgrm::TSimpleSumDgrm(const double KSum , const double ScnSum, const double Tension )
 {

	 mKSum = KSum ;
	 mScnSum =ScnSum ;
	 mTension = Tension ;
 }

//----------------------------------------------------------
void  _fastcall TSimpleSumDgrm::calcPartial_vectG_and_mtrxH(TComp cmpSZv, double alfTrg,double  alfAnp,double valK11, double valK12
	 ,double valK21,double valK22,double* arr_gradK11, double* arr_gradK12, double* arr_gradK21
	 ,double* arr_gradK22, double*arr_HessK11, double* arr_HessK12, double* arr_HessK21
	 ,double* arr_HessK22,double*   arr_FGreek, double*  arr_dFGreek )
{
  double valF1 = mKSum *fncDiagrSinx_div_x(mTension *( alfTrg -mScnSum )) ;
  double val_dF1 = mKSum * mTension *fncDerivDiagrSinx_div_x(mTension *( alfTrg -mScnSum)) ;
  double val_d2F1 = mKSum *mTension *mTension *fncDeriv2DiagrSinx_div_x(mTension *( alfTrg -mScnSum)) ;

  double valF2 = mKSum *fncDiagrSinx_div_x(mTension *( alfAnp - mScnSum )) ;
  double val_dF2 = mKSum * mTension *fncDerivDiagrSinx_div_x(mTension *(alfAnp - mScnSum )) ;
  double val_d2F2 = mKSum * mTension *mTension *fncDeriv2DiagrSinx_div_x(mTension *(alfAnp - mScnSum )) ;

  double valA1 = valK11 * valF1 + valK21 * valF2 -  cmpSZv.m_Re ;
  double valA2 = valK12* valF1 + valK22 * valF2 -  cmpSZv.m_Im ;

  // ?????????? ????????? G
  double arrB1[2] ={0.}, arrB2[2]={0.}, arrT0[2]= {0.}, arrT1[2]= {0.};
  arrB1[0] = valF1 * arr_gradK11[0] +  valF2 * arr_gradK21[0] + valK11 *val_dF1;
  arrB1[1] = valF1 * arr_gradK11[1] +  valF2 * arr_gradK21[1] + valK21 *val_dF2;

  arrB2[0] = valF1 * arr_gradK12[0] +  valF2 * arr_gradK22[0] + valK12 *val_dF1;
  arrB2[1] = valF1 * arr_gradK12[1] +  valF2 * arr_gradK22[1] + valK22 *val_dF2;

   MatrxMultScalar(arrB1, 1, 2, valA1,arrT0);
   MatrxMultScalar(arrB2, 1, 2, valA2,arrT1);
   MtrxSumMatrx(arrT0, arrT1,1, 2, arr_FGreek) ;
   ///

   // ?????????? ????????
   double arr_dB1_po_dx[4] = {0.}, arr_dB2_po_dx[4] = {0.};
	// arr_dB1_po_dx:
   arr_dB1_po_dx[0] = arr_HessK11[0] * valF1 + val_dF1* arr_gradK11[0] + valF2 * arr_HessK21[0]
	 + arr_gradK11[0]* val_dF1 + valK11 * val_d2F1;

   arr_dB1_po_dx[1] = valF1 * arr_HessK11[1] + valF2 * arr_HessK21[1] + val_dF2 * arr_gradK21[0]
	 + arr_gradK11[1] * val_dF1 ;

   arr_dB1_po_dx[2] =  valF1 * arr_HessK11[1] + val_dF1 * arr_gradK11[1] + valF2 * arr_HessK21[1]
	 + arr_gradK21[0] * val_dF2;

   arr_dB1_po_dx[3] = valF1 * arr_HessK11[3] + valF2 * arr_HessK21[3] + val_dF2 *arr_gradK21[1]
	+ val_dF2 *arr_gradK21[1] +  valK21 * val_d2F2;
	///

	// arr_dB2_po_dx:
	arr_dB2_po_dx[0] =  valF1 * arr_HessK12[0] + val_dF1 * arr_gradK12[0] +  valF2 * arr_HessK22[0]
	  +arr_gradK12[0] * val_dF1 + valK12 * val_d2F1;

	arr_dB2_po_dx[1] = valF1 * arr_HessK12[1] + arr_HessK22[1] * valF2 +  val_dF2* arr_gradK22[0]
	  + arr_gradK12[1] * val_dF1;

	arr_dB2_po_dx[2] = arr_HessK12[1] * valF1 + val_dF1 * arr_gradK12[1] + valF2 * arr_HessK22[1]
	  + arr_gradK22[0] * val_dF2;

	arr_dB2_po_dx[3] = valF1 * arr_HessK12[3] + valF2 * arr_HessK22[3] + val_dF2 * arr_gradK22[1]
	  + val_dF2 * arr_gradK22[1] + valK22 * val_d2F2;
	  ///

	  double arrT2[4] = {0.}, arrT3[4] = {0.}, arrT4[4] = {0.}, arrT5[4] = {0.};
	  arrT2[0] = arrB1[0]* arrB1[0];
	  arrT2[1] = arrB1[1]* arrB1[0];
	  arrT2[2] = arrB1[0]* arrB1[1];
	  arrT2[3] = arrB1[1]* arrB1[1];

	  MatrxMultScalar(arr_dB1_po_dx, 2, 2, valA1,arrT3);
	  MtrxSumMatrx(arrT3, arrT2,2,2,  arr_dFGreek ) ;

	  ///

	  arrT2[0] = arrB2[0]* arrB2[0];
	  arrT2[1] = arrB2[1]* arrB2[0];
	  arrT2[2] = arrB2[0]* arrB2[1];
	  arrT2[3] = arrB2[1]* arrB2[1];

	  MatrxMultScalar(arr_dB2_po_dx, 2, 2, valA2,arrT3);
	  MtrxSumMatrx(arrT3, arrT2,2,2,  arrT4 ) ;

	  MtrxSumMatrx(arrT4, arr_dFGreek,2,2,  arrT5 ) ;
	  memcpy(arr_dFGreek,  arrT5, 4 * sizeof(double));


}

#pragma package(smart_init)
