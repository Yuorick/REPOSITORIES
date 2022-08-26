//---------------------------------------------------------------------------


#pragma hdrstop
#include <vcl.h>
#include "Filt.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "ShipTarg.h"
#include "MatrixProccess.h"
#include "BombTraj.h"

//---------------------------------------------------------------------------
TFilt::TFilt()
{
	 mTf = 0;
// 2. Оценка вектора состояния угла визирования на момент mTf
	memset(marrEstAngVis, 0, 2 * sizeof(long double )) ;


//  3.корреляционная матрица ошибок оценивания
	memset(marrK_AngVis, 0, 4 * sizeof(long double )) ;


	marrK_AngVis[ 0 ] = 10000. ;
	marrK_AngVis[ 1] = 4000.;
	marrK_AngVis[ 2] = 4000.;
	marrK_AngVis[ 3] = 2500. ;

// 4. Оценка вектора состояния угла визирования на момент mTf
	memset(marrEstR, 0, 2 * sizeof(long double )) ;


//  5. корреляционная матрица ошибок оценивания в ГСК
	memset(marrK_R, 0, 4 * sizeof(long double )) ;


	marrK_R[ 0 ] = 10000. ;
	marrK_R[ 1] = 4000.;
	marrK_R[ 2] = 4000.;
	marrK_R[ 3] = 2500. ;

 // 6. Оценка вектора состояния цели  в ГСК по оси X на момент mTf
	memset(marrEstTargX, 0, 2 * sizeof(long double )) ;

  // 7.корреляционная матрица ошибок оценивания вектора состояния цели  в ГСК по оси X на момент mTf

	marrK_TargX[ 0 ] = 10000. ;
	marrK_TargX[ 1] = 4000.;
	marrK_TargX[ 2] = 4000.;
	marrK_TargX[ 3] = 2500. ;
	// 8.
	mbInit = true;

}




// конструктор копирования
 TFilt ::TFilt (const TFilt &R)
 {
	mTf = R.mTf;
	memcpy(marrEstAngVis, R.marrEstAngVis, 2 * sizeof(long double )) ;
	memcpy(marrK_AngVis, R.marrK_AngVis, 4 * sizeof(long double )) ;
	memcpy(marrEstR, R.marrEstR, 2 * sizeof(long double )) ;
	memcpy(marrK_R, R.marrK_R, 4 * sizeof(long double )) ;
	memcpy(marrEstTargX, R.marrEstTargX, 2 * sizeof(long double )) ;
	memcpy(marrK_TargX, R.marrK_TargX, 4 * sizeof(long double )) ;
	mbInit = R.mbInit ;
 }

 // оператор присваивания
 TFilt TFilt::operator=(TFilt  R)
 {
	mTf = R.mTf;
	memcpy(marrEstAngVis, R.marrEstAngVis, 2 * sizeof(long double )) ;
	memcpy(marrK_AngVis, R.marrK_AngVis, 4 * sizeof(long double )) ;
	memcpy(marrEstR, R.marrEstR, 2 * sizeof(long double )) ;
	memcpy(marrK_R, R.marrK_R, 4 * sizeof(long double )) ;
	memcpy(marrEstTargX, R.marrEstTargX, 2 * sizeof(long double )) ;
	memcpy(marrK_TargX, R.marrK_TargX, 4 * sizeof(long double )) ;
	mbInit = R.mbInit ;
	return *this ;
 }


   // парам конструктор 1
   TFilt::TFilt(TShipTarg ShipTarg, TBombTraj BombTraj)
{

	 mTf = ShipTarg.mT;
	 // вектор состояния корабля-цели
	 long double arrVS_Ship_GSK[4] = {0.};
	 arrVS_Ship_GSK[0 ] =  ShipTarg.mTraject.marrVectSostGSK[0] ;
	 arrVS_Ship_GSK[1 ] =  ShipTarg.mTraject.marrVectSostGSK[2] ;
	 arrVS_Ship_GSK[2] =   ShipTarg.mTraject.marrVectSostGSK[3] ;
	 arrVS_Ship_GSK[3 ] =  ShipTarg.mTraject.marrVectSostGSK[5] ;
	 ///

	// вектор состояния УАБ
   long double arrVS_Bomb_GSK[4] = {0.};
	 arrVS_Bomb_GSK[0 ] =  BombTraj.marrStrSK_VS[0] ;
	 arrVS_Bomb_GSK[1 ] =  BombTraj.marrStrSK_VS[1] ;
	 arrVS_Bomb_GSK[2] =   BombTraj.marrStrSK_VS [3] * cosl(BombTraj.marrStrSK_VS [4]) ;
	 arrVS_Bomb_GSK[3 ] =  BombTraj.marrStrSK_VS [3] * sinl(BombTraj.marrStrSK_VS [4]) ;
	 ///

	 // коррел матрица вектора сотояния ошибок оценивания угла визирования
	marrK_AngVis[ 0 ] = 10000. ;
	marrK_AngVis[ 1] = 4000.;
	marrK_AngVis[ 2] = 4000.;
	marrK_AngVis[ 3] = 2500. ;
	///

	// РАсчет оценок вектора сотояния угла визирования
   long double arrDelS[4] = {0.};
   MtrxMinusMatrx( arrVS_Ship_GSK, arrVS_Bomb_GSK,4, 1, arrDelS);// вектор разности
   long double temp0 = OuterProduct_2(arrDelS , &(arrVS_Ship_GSK[2])  )/NormVect2(arrDelS)/NormVect2(&(arrVS_Ship_GSK[2]) );
  if (fabs(temp0)> 1.)
  {
	ShowMessage(L"ERROR0");
	int iii = 0;
  }
 // marrEstAngVis[0]=  asinl(OuterProduct_2(arrDelS , &(arrVS_Ship_GSK[2])  )/NormVect2(arrDelS)/NormVect2(&(arrVS_Ship_GSK[2]) ));
  marrEstAngVis[0]=  asinl(temp0);

  long double  arrDelV [2] ={0.};
  MtrxMinusMatrx(&(arrVS_Ship_GSK[2]), &(arrVS_Bomb_GSK[2]),2, 1, arrDelV);

  long double temp1 = OuterProduct_2(arrDelS , arrDelV)/NormVect2(arrDelS)/NormVect2(arrDelS);
  if (fabs(temp1)> 1.)
  {
   //	ShowMessage(L"ERROR1");
   //	int iii = 0;
   temp1 = temp1/ fabsl(temp1);
  }
  marrEstAngVis[1]= asinl(temp1);

	 // коррел матрица вектора сотояния ошибок оценивания дальности
	marrK_R[ 0 ] = 10000. ;
	marrK_R[ 1] = 4000.;
	marrK_R[ 2] = 4000.;
	marrK_R[ 3] = 2500. ;

  // Расчет векора оценок по дальности

  marrEstR[0] =  NormVect2(arrDelS);
  marrEstR[1] =ScalProduct(arrDelS , &arrDelS[2], 2) / marrEstR[0] ;

  // оценка вектора состояния цели по оси X
	marrK_TargX[ 0 ] = 10000. ;
	marrK_TargX[ 1] = 4000.;
	marrK_TargX[ 2] = 4000.;
	marrK_TargX[ 3] = 2500. ;

	marrEstTargX[0] = arrVS_Ship_GSK[0 ] ;
	marrEstTargX[1] = arrVS_Ship_GSK[2] ;

	// 4.
  mbInit = true;
}



// шаг пересчета фильтра Голубева в чистом виде  , в соответсвии с формулами Борисенко
void 	TFilt ::OneStepGolubev(long double valFiVisZv,long double valRZv,long double valXZv,long double valTZv,long double sigmAngVis_w
	   ,long double sigmAngVis_bmo 	,long double sigmAngVis_mmo, long double sigmR_w,long double sigmR_bmo ,long double sigmR_mmo
	   , long double sigmX_w,long double sigmX_bmo ,long double sigmX_mmo)
{
	fncFiltStep(valFiVisZv, valTZv,sigmAngVis_w,sigmAngVis_bmo 	,sigmAngVis_mmo,marrEstAngVis, marrK_AngVis );
	fncFiltStep(valRZv, valTZv,sigmR_w,sigmR_bmo 	,sigmR_mmo,marrEstR, marrK_R );
	fncFiltStep(valXZv, valTZv,sigmX_w,sigmX_bmo 	,sigmX_mmo,marrEstTargX, marrK_TargX );


 

	mTf =  valTZv ;
}
// фильтр 2-го порядка
void 	TFilt ::fncFiltStep(long double valYZv,long double valTZv,long double sigm_w,long double sigm_bmo
,long double sigm_mmo,long double *arrXEst, long double *arrK )
{

	 long double h = valTZv - mTf;
	// Экстраполяция
	long double arrKExtr[4] ;
	long double sq = sqrt(arrK  [3] - arrK  [1] * arrK  [1] /arrK  [0]) ;
	arrKExtr[0] =  arrK  [0] + 2. * h * arrK [1] + h * h * arrK [3]
	 +  sigm_w * sigm_w * h* h * h * h/4. + sigm_w * h * h * h * sq;

	arrKExtr[1] =  arrK  [1] + h *arrK  [3] + 3. * sigm_w * h * h * sq/2. + sigm_w * sigm_w * h * h * h /2;
	arrKExtr[2] = arrKExtr[1] ;
	arrKExtr[3] =  arrK  [3] + 2. * sigm_w * h * sq + sigm_w * sigm_w * h * h ;
	// Выбор ММО
	long double alf =  sigm_mmo/ sqrtl( arrKExtr[0]) ;
	if (alf > 1.) alf = 1. ;
	long double valF = ( 1.- alf)/( ( 1.- alf)* ( 1.- alf) * arrKExtr[0] + sigm_bmo * sigm_bmo) ;

	long double arrP[2] = {0.};
	arrP[0] = arrKExtr[0]* valF ;
	arrP[1] =  valF * arrKExtr[1] ;
	if (arrP[0] > 1.)
	{
	  arrP[0] = 1. ;
	  arrP[1] = arrKExtr[1] / arrKExtr[0];
	  alf = ( sigm_mmo * sigm_mmo + sigm_bmo * sigm_bmo) / arrKExtr[0] ;
	}
	long double valP = 1. - (1.- alf) * arrP[0] ;
	arrK  [ 0] = valP * arrKExtr[0] ;
	arrK  [ 1] = valP * arrKExtr[1] ;
	arrK  [ 2] =  arrK  [ 1] ;
	arrK  [ 3] = arrKExtr[3] -  ( 1 - alf) * arrP[1] * arrKExtr[1] ;

	long double arrXExtrap[2] ={0.};
	arrXExtrap[0] = arrXEst[0] +  h * arrXEst[1];
	arrXExtrap[1] = arrXEst[1] ;

	long double del = valYZv - arrXExtrap[0];
	arrXEst[0] = arrXExtrap[0] + arrP[0]* del;
	arrXEst[1] = arrXExtrap[1] + arrP[1]* del;
	int iiii = 0;
}



void 	TFilt::OneStepFltrKlm_Real(long double valYZv,long double valTZv,long double sigm_w,long double sigm_bmo
,long double sigm_mmo,long double *arrXEst, long double *arrK )
{
  long double h = valTZv - mTf;


  long double arrKextr[4] = {0} , arrKt1[4] = {0},arrF[4] = {0} ;


  long double arrKtemp[4] ,arrKtemp1[4];
  memcpy(arrKtemp ,arrK, 4 * sizeof(long double)) ;
   long double arrT0[4] = {1.,0,0,1.}, arrA[4] ={0.}, arrB[2] = {0.},arrC[2] = {0.} ;

	arrT0[1] = h ;
	memcpy(arrA, arrT0, 4 * sizeof(long double)) ;    // матрица перехода за один такт;
	 arrB [0] = h * h/2. ;    // а атрицм помехи в канале объекта
	 arrB [1] = h ;
	 arrC[0] = 1. ;
	 arrC[1] = 0;

  arrF[0] = sigm_w* sigm_w* h * h * h /3 ;
  arrF [ 1] =  sigm_w* sigm_w* h * h  /2 ;
  arrF [ 2] = arrF [ 1] ;
  arrF [ 3] = sigm_w* sigm_w* h  ;
	long double arrP[2] = {0.};
	MtrxMultMatrx(arrA ,2, 2, arrK,2, arrKt1) ;
	MtrxMultMatrxTransp(arrKt1,2, 2,arrA ,2, arrKtemp1) ;
	MtrxSumMatrx(arrKtemp1, arrF,2,2, arrKextr) ;
	arrP [0] = arrKextr[ 0 ]/( sigm_bmo * sigm_bmo + arrKextr[ 0 ]) ;
	arrP [1] = arrKextr[ 1 ]/( sigm_bmo * sigm_bmo + arrKextr[ 0 ]) ;

	arrK [0] =  arrKextr[ 0 ] - arrKextr[ 0 ] * arrKextr[ 0 ] / ( sigm_bmo * sigm_bmo + arrKextr[ 0 ]) ;
	arrK [1] =  arrKextr[ 1 ] - arrKextr[ 0 ] * arrKextr[ 1 ] / ( sigm_bmo * sigm_bmo + arrKextr[ 0 ]) ;
	arrK [2] = arrK [1] ;
	arrK [3] =  arrKextr[ 3 ] - arrKextr[ 1 ] * arrKextr[ 1 ] / ( sigm_bmo * sigm_bmo + arrKextr[ 0 ]) ;
	long double arrXExtrap[2] ={0.};
	arrXExtrap[0] = arrXEst[0] +  h * arrXEst[1];
	arrXExtrap[1] = arrXEst[1] ;

	long double del = valYZv - arrXExtrap[0];
	arrXEst[0] = arrXExtrap[0] + arrP[0]* del;
	arrXEst[1] = arrXExtrap[1] + arrP[1]* del;
}



#pragma package(smart_init)






