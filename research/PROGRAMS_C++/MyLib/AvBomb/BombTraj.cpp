//---------------------------------------------------------------------------


#pragma hdrstop
#include <string.h>
#include <math.h>
#include "Atmosphere.h"
#include "MatrixProccess.h"
#include "YrWriteShapeFile.h"
#include "URPolyLine.h"
#include "URPolygon.h"
#include "FinalSituation.h"
#include "URPointXY.h"
#include "ShipTraj.h"

#include "BombTraj.h"
#include "CntlFuncPar.h"
const int QUANT_COLS = 10 , QUANT_POINTS_MAX = 3000;

  // ����������� ��
  // ������� ���� ��� ������������� �� ��� OX ����������� �� ������ ������� ������� !!!

//---------------------------------------------------------------------------
TBombTraj::TBombTraj()
{
	mTet0 = 40. * M_PI / 180;  // ��� ���� ��������  40 ����
   //	miControlType  = 0; // ����������

	mAltit = 6000.; // ��� ������
	mX0 = 0. ;
	mTStart = 0.;
	mTCur = 0.;  // ������ ��������� �������� �������
	mStepInt = 0.0005; // ���� �������������� �� ������� ��� ��������� ��������
	memset(marrCfW, 0, 5 * sizeof(long double)) ;

	mBomb =TBomb();  // ���

	// �������  ������  � ����������� ��:
//  marrStrSK_VS [0]- x
//  marrStrSK_VS [1]-  y
//  marrStrSK_VS [2]-  z
//  marrStrSK_VS [3]-  ���� ���
//  marrStrSK_VS [4]-  ������� ��������� �����
//  marrStrSK_VS [5]-  ��������� ��������� ��
//  marrStrSK_VS [6]-  ������� �������� V
//  marrStrSK_VS [7]-  ���� �����
	fncFillNachalnieUsloviaVS();


}

//---------------------------------------------------------------------------


// ����������� �����������
 TBombTraj ::TBombTraj (const TBombTraj &R)
 {
	mTet0 = R.mTet0;  // ��� ���� ��������  17 ����
	mV0 = R.mV0 ;
	mStepInt = R.mStepInt ;
	mTStart = R.mTStart ;
	mAltit = R.mAltit; // ��� ������
	mX0 = R.mX0;
	mTCur = R.mTCur ;  // ������ ��������� �������� �������
	mBomb = R.mBomb ;  // ������
	memcpy(marrStrSK_VS, R.marrStrSK_VS, 5 * sizeof(long double));
	memcpy(marrCfW, R.marrCfW, 5 * sizeof(long double));


 }

 // �������� ������������
 TBombTraj TBombTraj::operator=(TBombTraj  R)
 {
	mTet0 = R.mTet0;  // ��� ���� ��������  17 ����
	mV0 = R.mV0 ;
	mStepInt = R.mStepInt ;
	mTStart = R.mTStart ;
	mAltit = R.mAltit; // ��� ������
	mX0 = R.mX0;  // ���. ��������� �� ��� X
	mTCur = R.mTCur ;  // ������ ��������� �������� �������
	mBomb = R.mBomb ;  // ������
	memcpy(marrStrSK_VS, R.marrStrSK_VS, 5 * sizeof(long double));
	memcpy(marrCfW, R.marrCfW, 5 * sizeof(long double));
	return *this ;
 }

	// ����� �����������
	TBombTraj::TBombTraj (	long double Tet0  // ��� ���� ��������
	,long double V0  // ��� ��������
	,long double Altit  // ��� ������
	,long double X0  // ���. ��������� �� ��� X
	,long double StepInt  // ���� �������������� �� ������� ��� ��������� ��������
	 ,long double TStart
	)
	{
	   mBomb  = TBomb() ;
	   mTStart = TStart ;
	   mTCur = TStart ;
	   mTet0 = Tet0 ;
	   mV0 = V0 ;
	   mAltit =  Altit ;
	   mX0 =  X0 ;

	   mStepInt = StepInt ;
	   memset(marrCfW, 0, 5 * sizeof(long double)) ;
	   fncFillNachalnieUsloviaVS();

	}
 // ����� ����������� 2
 TBombTraj::TBombTraj (
	long double Tet0  // ��� ���� ��������
	,long double V0  // ��� ��������
	,long double Altit  // ��� ������
	 ,long double X0  // ���. ��������� �� ��� X
	,long double StepInt  // ���� �������������� �� ������� ��� ��������� ��������
	 ,long double TStart
	,TBomb Bomb
	)
	{
	   mBomb  = Bomb ;
	   mTStart = TStart ;
	   mTCur = TStart ;

	   mX0 =  X0 ;
	   mTet0 = Tet0 ;
	   mV0 = V0 ;
	   mAltit =  Altit ;

	   mStepInt = StepInt ;
	   memset(marrCfW, 0, 5 * sizeof(long double)) ;
	   fncFillNachalnieUsloviaVS();


	}

	 // ����� ����������� 3
 TBombTraj::TBombTraj (
	long double Tet0  // ��� ���� ��������
	,long double V0  // ��� ��������
	,long double Altit  // ��� ������
	 ,long double X0  // ���. ��������� �� ��� X
	,long double StepInt  // ���� �������������� �� ������� ��� ��������� ��������
	 ,long double TStart
	,TBomb Bomb
	,long double *arrCfW
	)
	{
	   mBomb  = Bomb ;
	   mTStart = TStart ;
	   mTCur = TStart ;

	   mX0 =  X0 ;
	   mTet0 = Tet0 ;
	   mV0 = V0 ;
	   mAltit =  Altit ;

	   mStepInt = StepInt ;
	   memcpy(marrCfW, arrCfW, 5 * sizeof(long double)) ;
	   fncFillNachalnieUsloviaVS();


	}
///
 void TBombTraj::fncFillNachalnieUsloviaVS()
{
  marrStrSK_VS[0] = mX0;
  marrStrSK_VS[1] = mAltit ;
  marrStrSK_VS[2] = fncCalcOtnositP(marrStrSK_VS[1]) ;//
  marrStrSK_VS[3] = mV0 ;
  marrStrSK_VS[4] = mTet0 ;

}

// ���������� ������ ����� �������� ��� ���������
void TBombTraj::fncCalc_F(const long double valU, long double *arrF)
{

  // ������� ������
  // const long double valSm = M_PI * mBomb.mSurfFuzSq * mBomb.mDm / 4. ;
  // ���������� ���������� ����������� �����������  � �� ���������
  long double valTay = 0.,valDerivTay = 0. ;
  fncCalcNormTemperature(marrStrSK_VS[1],valTay, valDerivTay)  ;
  ///

  // ���������� ����� ���� � ������� ��� ������� �����������

  long double valMach = 0. ;
  fncCalcMach(valTay,valDerivTay, valMach);
  ///

  // ���������� ������ q � ������� �� ���������
  long double val_q = 0. , arrGrad_q [5] = {0.} ;;

  fncCalc_q_and_Grad_q(val_q, arrGrad_q) ;
  ///


   long double valCy = mBomb.fncCy( valMach);
 // long double valW =  valCy * val_q * mBomb.mMidSq   /mBomb.mMass  *  valU ;
   ///////////////////////
   arrF[0] =  marrStrSK_VS [3] * cosl(marrStrSK_VS [4]) ;
   arrF[1] =  marrStrSK_VS [3] * sinl(marrStrSK_VS [4]) ;
   arrF[2] = -marrStrSK_VS [2] * marrStrSK_VS [3] * sinl(marrStrSK_VS [4]) / valTay*(G_ZEMLI /ATM_R_UNIVER + valDerivTay);//HT[0][2]);
   arrF[3] = -G_ZEMLI * sinl(marrStrSK_VS [4])   -  mBomb.fncCx( valMach) * val_q * mBomb.mMidSq  /mBomb.mMass   ;
   arrF[4] = -G_ZEMLI * cosl(marrStrSK_VS [4])/ marrStrSK_VS [3]
   - mBomb.fncCy( valMach) * val_q * mBomb.mMidSq   /mBomb.mMass /  marrStrSK_VS [3] *   valU ;


}

// arrF[5] -  ������ ����� �������� ��� ���������
// arr_dF_po_dx[25] -������� ������� ����������� ������� ������ ������ ����� �� ������� ���������
// arr_dF_po_dU[5] - ������ ������� ����������� ������� ������ ����� �� ���������� ���������� U (= AlfaP)
void TBombTraj::fncCalc_F_and_dF_po_dx_and_dF_po_dU(const long double valU, long double *arrF,long double *arr_dF_po_dx,long double *arr_dF_po_dU)
{
  // ���������� ���������� ����������� �����������  � �� ���������
  long double valTay = 0.,valDerivTay = 0. ;
  fncCalcNormTemperature(marrStrSK_VS[1],valTay, valDerivTay)  ;
  ///

  // ���������� ����� ���� � ������� ��� ������� �����������
  long double valMach = 0.,arrGradMach[5] ={0.} ;
  fncCalcMach_and_GradMach(valTay, valDerivTay ,valMach, arrGradMach);
  ///

  // ���������� ������ q � ������� �� ���������
  long double val_q = 0., arrGrad_q[5] = {0.} ;
  fncCalc_q_and_Grad_q( val_q, arrGrad_q);
  ///

  // ���������� ������ Cx � ������� �� ���������
  long double val_Cx = 0., arrGrad_Cx[5] = {0.} ;
  fncCalc_Cx_and_Grad_Cx(valMach,arrGradMach, val_Cx, arrGrad_Cx);
  ///

  // ���������� ������ Cy � ������� �� ���������
  long double val_Cy = 0., arrGrad_Cy[5] = {0.} ;
  fncCalc_Cy_and_Grad_Cy(valMach,arrGradMach, val_Cy, arrGrad_Cy);
  ///



	 ///

   arrF[0] =  marrStrSK_VS [3] * cosl(marrStrSK_VS [4]) ;
   arrF[1] =  marrStrSK_VS [3] * sinl(marrStrSK_VS [4]) ;
   arrF[2] = - marrStrSK_VS [2] * marrStrSK_VS [3] * sinl(marrStrSK_VS [4]) / valTay*(G_ZEMLI /ATM_R_UNIVER + valDerivTay);//HT[0][2]);
   arrF[3] = -G_ZEMLI * sinl(marrStrSK_VS [4])   -  val_Cx * val_q * mBomb.mMidSq  /mBomb.mMass   ;
   arrF[4] = -G_ZEMLI * cosl(marrStrSK_VS [4])/ marrStrSK_VS [3]
   - val_Cy * val_q * mBomb.mMidSq   /mBomb.mMass /  marrStrSK_VS [3] *  valU  ;

   // ���������� ������� ������� ������������
   memset(arr_dF_po_dx , 0, 25 * sizeof(long double ));
   memset(arr_dF_po_dU , 0,  5 * sizeof(long double ));
   fncGradF0( arr_dF_po_dx);
   fncGradF1( &arr_dF_po_dx[5]);
   fncGradF2( arrF[2], valTay,  valDerivTay , &arr_dF_po_dx[10]);
   fncGradF3(  val_q , arrGrad_q, val_Cx,  arrGrad_Cx, &arr_dF_po_dx[15]);
   fncGradF4(  val_q , arrGrad_q, val_Cy,  arrGrad_Cy, valU, &arr_dF_po_dx[20]);
   fncCalc_dF_po_dU(val_Cy,  val_q,arr_dF_po_dU) ;
}

//���������� ������� ������� ����������� ������� ������ ����� �� ���������� ���������� U (= AlfaP)
void TBombTraj::fncCalc_dF_po_dU(long double val_Cy, long double val_q, long double *arr_dF_po_dU)
{
  memset( arr_dF_po_dU, 0, 5 * sizeof( long double)) ;
  arr_dF_po_dU[4] = - val_Cy * val_q * mBomb.mMidSq   /mBomb.mMass /  marrStrSK_VS [3]  ;

}
// ���������� ��������� ������� F0 �� ������� ����������
void TBombTraj::fncGradF0(long double *arr_dF0_po_dx)
{
  memset( arr_dF0_po_dx, 0, 5 * sizeof( long double)) ;
  arr_dF0_po_dx[3] =  cosl(marrStrSK_VS [4]) ;
  arr_dF0_po_dx[4] = -marrStrSK_VS [3] * sinl(marrStrSK_VS [4] );
}
// ���������� ��������� ������� F1 �� ������� ����������
void TBombTraj::fncGradF1(long double *arr_dF1_po_dx)
{
  memset( arr_dF1_po_dx, 0, 5 * sizeof( long double)) ;
  arr_dF1_po_dx[3] =  sinl(marrStrSK_VS [4]) ;
  arr_dF1_po_dx[4] =  marrStrSK_VS [3] * cosl(marrStrSK_VS [4] );
}

// ���������� ��������� ������� F2 �� ������� ����������
void TBombTraj::fncGradF2(long double valF2,long double valTay, long double valDerivTay , long double *arr_dF2_po_dx)
{
  memset( arr_dF2_po_dx, 0, 5 * sizeof( long double)) ;
  arr_dF2_po_dx[1] = - valF2 * valDerivTay / valTay;
  arr_dF2_po_dx[2] =   valF2  / marrStrSK_VS [2];
  arr_dF2_po_dx[3] =   valF2  / marrStrSK_VS [3];
  arr_dF2_po_dx[4] =   valF2  / tanl(marrStrSK_VS [4]);
}

// ���������� ��������� ������� F3 �� ������� ����������
void TBombTraj::fncGradF3(long double val_q , long double *arrGrad_q
	, long double val_Cx, long double* arrGrad_Cx, long double *arr_dF3_po_dx)
{
  memset( arr_dF3_po_dx, 0, 5 * sizeof( long double)) ;
  arr_dF3_po_dx[1] = - mBomb.mMidSq / mBomb.mMass * arrGrad_Cx[1] * val_q  ;
  arr_dF3_po_dx[2] = - mBomb.mMidSq / mBomb.mMass * val_Cx * arrGrad_q[2] ;
  arr_dF3_po_dx[3] = - mBomb.mMidSq / mBomb.mMass *( arrGrad_Cx[3] * val_q  + val_Cx * arrGrad_q[3]);
  arr_dF3_po_dx[4] = -G_ZEMLI * cosl(marrStrSK_VS [4] );
}

// ���������� ��������� ������� F4 �� ������� ����������
void TBombTraj::fncGradF4(long double val_q , long double *arrGrad_q
	, long double val_Cy, long double* arrGrad_Cy, long double valU, long double *arr_dF4_po_dx)
{
  memset( arr_dF4_po_dx, 0, 5 * sizeof( long double)) ;
  arr_dF4_po_dx[1] = - mBomb.mMidSq / mBomb.mMass / marrStrSK_VS [3] *valU  * arrGrad_Cy[1] * val_q  ;
  arr_dF4_po_dx[2] = - mBomb.mMidSq / mBomb.mMass / marrStrSK_VS [3] *  val_Cy * arrGrad_q[2] *valU;
  arr_dF4_po_dx[3] = G_ZEMLI * cosl(marrStrSK_VS [4] )/ marrStrSK_VS [3]/ marrStrSK_VS [3]
		- mBomb.mMidSq / mBomb.mMass   / marrStrSK_VS [3]/ marrStrSK_VS [3] *valU
		* (( arrGrad_Cy[3] * val_q  + val_Cy * arrGrad_q[3]) * marrStrSK_VS [3] - val_Cy * val_q);
  arr_dF4_po_dx[4] = G_ZEMLI * sinl(marrStrSK_VS [4] )/ marrStrSK_VS [3];
}

 // ����������� ����������� ������
// INPUT:
// valMach - ����� ����
// OUTPUT :
// val_q - ���������� �����

void TBombTraj::fncCalc_q(long double valMach,long double  &val_q)
{


   val_q = ATM_RoN0  * marrStrSK_VS [2]* marrStrSK_VS [3]* marrStrSK_VS [3] / 2. ;


}



 // ���������� ����� ����  �� x
// INPUT:
//valTay - ���������� �����������  �����������
// valDerivTay - ����������� �� ������ ���������� �����������  �����������
// OUTPUT
// valMach - ����� ����

void TBombTraj::fncCalcMach(long double valTay, long double valDerivTay ,long double &valMach)
{

   valMach = sqrtl(ATM_TAYN0/ valTay)/ ATM_AN0 * marrStrSK_VS [3] ;


}

// ���������� ����� ���� � ������� ������ ����������� ����� ���� �� x
// INPUT:
//valTay - ���������� �����������  �����������
// valDerivTay - ����������� �� ������ ���������� �����������  �����������
// OUTPUT
// valMach - ����� ����
// arrGradMach - ������ ��������� ����� ���� �� ������� ����������  ����������� 5
void TBombTraj::fncCalcMach_and_GradMach(long double valTay, long double valDerivTay
	 ,long double &valMach, long double *arrGradMach)
{
   memset( arrGradMach, 0, 5 * sizeof(long double));
   arrGradMach [3] = sqrtl(ATM_TAYN0/ valTay)/ ATM_AN0 ;
   valMach = marrStrSK_VS[3] * arrGradMach [3] ;
   arrGradMach [1] = - 0.5 * valMach * valDerivTay/ valTay ;

}
// ����������� ����������� ������ � ������� ������� ����������� �������� ������ �� x
// INPUT:
// valMach - ����� ����
// arrGradMach[5] - ������ ��������� ����� ���� �� ������� ����������
// OUTPUT :
// val_q - ���������� �����
// arrGrad_q[5]  -  ������ ��������� c���������� ������  �� ������� ����������
void TBombTraj::fncCalc_q_and_Grad_q(long double  &val_q, long double *arrGrad_q)
{
  memset( arrGrad_q, 0, 5 * sizeof(long double));
  arrGrad_q [2] = ATM_RoN0 * marrStrSK_VS [3]* marrStrSK_VS [3]  / 2. ;  // �� �������� ����
  arrGrad_q [3] =  ATM_RoN0 * marrStrSK_VS [2]* marrStrSK_VS [3] ;  // �� ��������
  val_q = ATM_RoN0  * marrStrSK_VS [2]* marrStrSK_VS [3]* marrStrSK_VS [3] / 2. ;

}

// ����������� ������������ �������� �������������  �
// ��� ������� ������� �����������  �� x
// INPUT:
// valMach - ����� ����
// arrGradMach[5] - ������ ��������� ����� ���� �� ������� ����������
// OUTPUT :
// val_�x - ����������� �������� �������������
// arrGrad_Cx[5]  -  ������ ��������� ������������ �������� �������������  �� ������� ����������
void TBombTraj::fncCalc_Cx_and_Grad_Cx(long double valMach,long double  *arrGradMach
		,long double  &val_Cx, long double *arrGrad_Cx)
{
  val_Cx = mBomb.fncCx(valMach) ;
  memset( arrGrad_Cx, 0,  5 * sizeof(long double));
  long double dCx_po_dMach = mBomb.fnc_dCx_po_dMach(valMach) ;
  arrGrad_Cx[1] =  dCx_po_dMach *  arrGradMach[1] ;
  arrGrad_Cx[3] =  dCx_po_dMach *  arrGradMach[3] ;
}

// ����������� ������������ ��������� ���� �
// ��� ������� ������� �����������  �� x
// INPUT:
// valMach - ����� ����
// arrGradMach[5] - ������ ��������� ����� ���� �� ������� ����������
// OUTPUT :
// val_�y - ����������� ��������� ����
// arrGrad_Cy[5]  -  ������ ������������ ��������� ���� �� ������� ����������
void TBombTraj::fncCalc_Cy_and_Grad_Cy(long double valMach,long double  *arrGradMach
		,long double  &val_Cy, long double *arrGrad_Cy)
{
  val_Cy = mBomb.fncCy(valMach) ;
  memset( arrGrad_Cy, 0,  5 * sizeof(long double));
  long double dCy_po_dMach = mBomb.fnc_dCy_po_dMach(valMach) ;
  arrGrad_Cy[1] =  dCy_po_dMach *  arrGradMach[1] ;
  arrGrad_Cy[3] =  dCy_po_dMach *  arrGradMach[3] ;
}
// ��� ����� ������ �� ���� valStepInt
 void TBombTraj::fncEilerStep(const long double valU, const long double valStepInt)
 {
   long double arrF[5] ={0.} ;
   // ��� �������������� �������� �������
   fncCalc_F(valU, arrF) ;
   for (int i = 0; i < 5; i++)
   {
	arrF[i] += TShipTraj::getGauss(0., arrF[i] * marrCfW[i] ) ;
   }
   long double arrT0[25] = {0.},arrT1[25] = {0.};
   MatrxMultScalar(arrF, 5, 1, valStepInt,arrT0); // f * dt
   MtrxSumMatrx(arrT0, marrStrSK_VS,5, 1, arrT1) ;
   memcpy(marrStrSK_VS, arrT1, 5 * sizeof(long double)) ;   ///

  mTCur += valStepInt ;
 }




#pragma package(smart_init)
