//---------------------------------------------------------------------------


#pragma hdrstop

#include "CalcCorMatrx.h"
#include <string.h>
#include <math.h>
#include "MatrixProccess.h"


//---------------------------------------------------------------------------
void calcMatrP1(const double a, double*matrP)
{
  memset(matrP,0,9*sizeof(double));
  matrP[0] = 1;
  matrP[4] = cos(a);
  matrP[5] = sin(a);
  matrP[7] = -sin(a);
  matrP[8] = cos(a);

}
void calcMatrP2(const double a, double*matrP)
{
  memset(matrP,0,9*sizeof(double));
  matrP[4] = 1;
  matrP[0] = cos(a);
  matrP[2] = sin(a);
  matrP[6] = -sin(a);
  matrP[8] = cos(a);

}
void calcMatrP3(const double a, double*matrP)
{
  memset(matrP,0,9*sizeof(double));
  matrP[8] = 1;
  matrP[0] = cos(a);
  matrP[1] = sin(a);
  matrP[3] = -sin(a);
  matrP[4] = cos(a);
}

//---------------------------------------------------------------------------
void calcMatr_dP1_po_da(const double a, double*matrP)
{
  memset(matrP,0,9*sizeof(double));

  matrP[4] = -sin(a);
  matrP[5] = cos(a);
  matrP[7] = -cos(a);
  matrP[8] = -sin(a);

}
//---------------------------------------------------------------------------
void calcMatr_dP2_po_da(const double a, double*matrP)
{
  memset(matrP,0,9*sizeof(double));

  matrP[0] = -sin(a);
  matrP[2] = cos(a);
  matrP[6] = -cos(a);
  matrP[8] = -sin(a);

}
//---------------------------------------------------------------------------
void calcMatr_dP3_po_da(const double a, double*matrP)
{
  memset(matrP,0,9*sizeof(double));

  matrP[0] = -sin(a);
  matrP[1] = cos(a);
  matrP[3] = -cos(a);
  matrP[4] = -sin(a);

}



void calcMatr_ASK_v_KGSK(double *arrMu,double *matrPereh_ASK_V_KGSK)
{
	/*const double valQ   = arrMu[0] ;
	const double valPsi = arrMu[1] ;
	const double valTet = arrMu[2] ;
	const double valBet = arrMu[3] ;
	const double valEps = arrMu[4] ;
   double arr1[9]={0},arr0[9]={0},arr2[9]={0}
	 ,matrP3Q[9]={0},matrP1Psi[9]={0},matrP2Tet[9]={0},matrP3Bet[9]={0},matrP1Eps[9]={0};

	 calcMatrP3(valQ, matrP3Q) ;
	 calcMatrP1(valPsi, matrP1Psi) ;
	 calcMatrP2(valTet, matrP2Tet) ;
	 calcMatrP3(valBet, matrP3Bet) ;
	 calcMatrP1(-valEps, matrP1Eps) ;
	 MtrxMultMatrx( matrP3Q,3, 3, matrP1Psi,3, arr0);
	 MtrxMultMatrx( arr0,3, 3, matrP2Tet,3, arr1);
	 MtrxMultMatrx( arr1,3, 3, matrP3Bet,3, arr2);
	 MtrxMultMatrx( arr2,3, 3, matrP1Eps,3, matrPereh_ASK_V_KGSK); */
	 double  matrPereh_ASK_V_PSK[9] = {0.}, matrPereh_PSK_V_KGSK[9] ={0.} ;
	 calcMatr_ASK_v_PSK( &arrMu[3],matrPereh_ASK_V_PSK) ;
	 calcMatr_PSK_v_KGSK(arrMu,matrPereh_PSK_V_KGSK);
	 MtrxMultMatrx(matrPereh_PSK_V_KGSK,3, 3, matrPereh_ASK_V_PSK,3, matrPereh_ASK_V_KGSK ) ;


}

void calcMatr_ASK_v_PSK(double *arrMu,double *matrPereh_ASK_V_PSK)
{

	const double valBet = arrMu[0] ;
	const double valEps = arrMu[1] ;
	double matrP3Bet[9]={0},matrP1Eps[9]={0};
	calcMatrP3(valBet, matrP3Bet) ;
	calcMatrP1(-valEps, matrP1Eps) ;
	MtrxMultMatrx( matrP3Bet,3, 3, matrP1Eps,3, matrPereh_ASK_V_PSK);

}

void calcMatr_PSK_v_KGSK(double *arrMu,double *matrPereh_PSK_V_KGSK)
{
	const double valQ   = arrMu[0] ;
	const double valPsi = arrMu[1] ;
	const double valTet = arrMu[2] ;

	double arr0[9]={0},matrP3Q[9]={0},matrP1Psi[9]={0},matrP2Tet[9]={0};

	calcMatrP3(valQ, matrP3Q) ;
	calcMatrP1(valPsi, matrP1Psi) ;
	calcMatrP2(valTet, matrP2Tet) ;
	MtrxMultMatrx( matrP3Q,3, 3, matrP1Psi,3, arr0);
	MtrxMultMatrx( arr0,3, 3, matrP2Tet,3, matrPereh_PSK_V_KGSK);

}


void calc_dF_po_dQ_sq(const double valBet,const double valEps,double *matrRez)
{
	double arr[3]={0},arr0[3]={0};
	arr[0] = cos(valBet) * cos ( valEps) ;
	arr[1] = -sin(valBet) * cos ( valEps) ;
	memcpy(arr0,arr,3*sizeof(double));
	MtrxMultMatrxTransp(arr,3, 1, arr0,3, matrRez);
}
void calc_dF_po_dPsi_sq(const double valBet,const double valEps,double *matrRez)
{
	double arr[3]={0},arr0[3]={0};
	arr[1] = sin(valEps) ;
	arr[2] = -cos(valBet) * cos ( valEps) ;
	memcpy(arr0,arr,3*sizeof(double));
	MtrxMultMatrxTransp(arr,3, 1, arr0,3, matrRez);

}
void calc_dF_po_dTet_sq(const double valBet,const double valEps,double *matrRez)
{
	double arr[3]={0},arr0[3]={0};
	arr[0] = sin(valEps) ;
	arr[2] = -sin(valBet) * cos ( valEps) ;
	memcpy(arr0,arr,3*sizeof(double));
	MtrxMultMatrxTransp(arr,3, 1, arr0,3, matrRez);

}
void calc_dF_po_dBet_sq(const double valBet,const double valEps,double *matrRez)
{
	double arr[3]={0},arr0[3]={0};
	arr[0] = cos(valBet) * cos ( valEps) ;
	arr[1] = -sin(valBet) * cos ( valEps) ;
	memcpy(arr0,arr,3*sizeof(double));
	MtrxMultMatrxTransp(arr,3, 1, arr0,3, matrRez);

}
void calc_dF_po_dEps_sq(const double valBet,const double valEps,double *matrRez)
{
	double arr[3]={0},arr0[3]={0};
	arr[0] =  -sin(valBet) *sin(valEps) ;
	arr[1] = -cos(valBet) * sin ( valEps) ;
	arr[2] =  cos ( valEps) ;
	memcpy(arr0,arr,3*sizeof(double));
	MtrxMultMatrxTransp(arr,3, 1, arr0,3, matrRez);

}
void createExtendMtrx(double *arrInp, double * arrOut)
{
  memset(arrOut, 0, 36 * sizeof(double));
 for (int i =0; i < 3; i++)
 for (int j =0; j < 3; j++)
 {
  arrOut[ i * 6 + j] =  arrInp [3 * i + j];
  arrOut[ (3 + i) * 6 + 3 + j] =  arrInp [3 * i + j];
 }
}
void calcMatrJ1(const double valBet,const double valEps,double *arrJ1)
{
 memset(arrJ1, 0, 9 * sizeof(double)) ;
 arrJ1[1] = cos(valEps);
 arrJ1[2] = -sin(valEps);
 arrJ1[3] = -cos(valEps);
 arrJ1[6] = sin(valEps);

}
void calcMatrJ2(const double valBet,const double valEps,double *arrJ2)
{
 memset(arrJ2, 0, 9 * sizeof(double)) ;
 arrJ2[1] = -sin(valBet)*sin(valEps);
 arrJ2[2] = -sin(valBet)*cos(valEps);
 arrJ2[3] = sin(valBet)*sin(valEps);
 arrJ2[5] = cos(valBet);
 arrJ2[6] = sin(valBet)*cos(valEps);
 arrJ2[7] = -cos(valBet);


}
void calcMatrJ3(const double valBet,const double valEps,double *arrJ3)
{
 memset(arrJ3, 0, 9 * sizeof(double)) ;
 arrJ3[1] = cos(valBet)*sin(valEps);
 arrJ3[2] = cos(valBet)*cos(valEps);
 arrJ3[3] = -cos(valBet)*sin(valEps);
 arrJ3[5] = sin(valBet);
 arrJ3[6] = -cos(valBet)*cos(valEps);
 arrJ3[7] = -sin(valBet);

}
void calcMatrJ4(const double valBet,const double valEps,double *arrJ4)
{
 calcMatrJ1( valBet, valEps,arrJ4) ;

}
void calcMatrJ5(const double valBet,const double valEps,double *arrJ5)
{
 memset(arrJ5, 0, 9 * sizeof(double)) ;

 arrJ5 [ 5] = -1 ;
 arrJ5 [ 7] = 1 ;

}
void calcMatrJ1W(const double valBet,const double valEps,double *arrJ1W)
{
   memset(arrJ1W, 0, 36 * sizeof(double)) ;
   double arr[9]={0} ;
   calcMatrJ1( valBet, valEps,arr);
   createExtendMtrx(arr, arrJ1W);
}
void calcMatrJ2W(const double valBet,const double valEps,double *arrJ2W)
{
   memset(arrJ2W, 0, 36 * sizeof(double)) ;
   double arr[9]={0} ;
   calcMatrJ2( valBet, valEps,arr);
   createExtendMtrx(arr, arrJ2W);
}
void calcMatrJ3W(const double valBet,const double valEps,double *arrJ3W)
{
   memset(arrJ3W, 0, 36 * sizeof(double)) ;
   double arr[9]={0} ;
   calcMatrJ3( valBet, valEps,arr);
   createExtendMtrx(arr, arrJ3W);
}
void calcMatrJ4W(const double valBet,const double valEps,double *arrJ4W)
{
   memset(arrJ4W, 0, 36 * sizeof(double)) ;
   double arr[9]={0} ;
   calcMatrJ4( valBet, valEps,arr);
   createExtendMtrx(arr, arrJ4W);
}
void calcMatrJ5W(const double valBet,const double valEps,double *arrJ5W)
{
   memset(arrJ5W, 0, 36 * sizeof(double)) ;
   double arr[9]={0} ;
   calcMatrJ5( valBet, valEps,arr);
   createExtendMtrx(arr, arrJ5W);
}

//--------------------------------------------------------------------------------
// Вычисление корреляционнной матрицы ошибок , возникающих при переводе вектора из ПСК в КГСК.
// ОШибки возникают в следствие того, что палубные углы известны неточно, а с ошибками измерений.
// Задан вектор arrVS[3], измерения палубных углов и дисперсии ошибок измерения этих углов.
// Требуется определить коррел матрицу ошибок преобразованного вектора arrVS
// INPUT:
// VAlQ,  VAlPsi,  VAlTet  - измерения курсового угда, угла килевой и бортовой качек
// VAlDispQ, VAlDispPsi	 , VAlDispTet  - дисперсии ошибок измерений этих углов
// arrVS[3] - вектор в ПСК
// OUTPUT:
// arrVS_KGSK[3] - вектор в КГСК
// arrMtxCorr[9] - корреляционная матрица ошибок
void calcCorMatrx_PSK_KGSK(const double VAlQ, const double VAlPsi, const double VAlTet
	 , const double VAlDispQ,const double VAlDispPsi	 , const double VAlDispTet
	 , double  *arrVS_PSK,double  *arrVS_KGSK, double *arrMtxCorr)
{
	double arrP1[9] = {0.},arrP2[9] = {0.},arrP3[9] = {0.},
		arr_dP1[9] = {0.},arr_dP2[9] = {0.},arr_dP3[9] = {0.};
	calcMatrP1(VAlPsi, arrP1);
	calcMatr_dP1_po_da(VAlPsi, arr_dP1);

	calcMatrP2(VAlTet, arrP2);
	calcMatr_dP2_po_da(VAlTet, arr_dP2);

	calcMatrP3(VAlQ, arrP3);
	calcMatr_dP3_po_da(VAlQ, arr_dP3);

	double arrJ1[9] = {0.}, arrJ2[9] = {0.}, arrJ3[9] = {0.};
	MtrxMultMatrx_MultMatrx(arr_dP3, arrP1,arrP2,3, arrJ1)  ;
	MtrxMultMatrx_MultMatrx(arrP3, arr_dP1,arrP2,3, arrJ2)  ;
	MtrxMultMatrx_MultMatrx(arrP3, arrP1,arr_dP2,3, arrJ3)  ;
	///
 memset(arrMtxCorr, 0, 9 * sizeof(double));
 double arrQ[9] = {0.}, arrPsi[9] = {0.}, arrTet[9] = {0.}, arrTemp0[9] = {0.}, arrVS_PSK_Copy[3] ={0.}, arrS_mult_ST[9] = {0.};
 memcpy(arrVS_PSK_Copy, arrVS_PSK, 3 * sizeof(double));
 MtrxMultMatrxTransp(arrVS_PSK,3, 1, arrVS_PSK_Copy,3, arrS_mult_ST) ;
 MtrxMultMatrx_MultMatrxTransp(arrJ1,arrS_mult_ST, arrJ1, 3, arrQ);
 MtrxMultMatrx_MultMatrxTransp(arrJ2,arrS_mult_ST, arrJ2, 3, arrPsi);
 MtrxMultMatrx_MultMatrxTransp(arrJ3,arrS_mult_ST, arrJ3, 3, arrTet);

 MatrxMultScalar(arrQ, 3, 3,VAlDispQ,arrQ) ;
 MatrxMultScalar(arrPsi, 3, 3,VAlDispPsi,arrPsi) ;
 MatrxMultScalar(arrTet, 3, 3,VAlDispTet,arrTet) ;

 MtrxSumMatrx(arrQ, arrPsi,3, 3, arrTemp0) ;
 MtrxSumMatrx(arrTemp0, arrTet,3, 3, arrMtxCorr) ;

 ///

MtrxMultMatrx_MultMatrx(arrP3, arrP1,arrP2,3, arrJ1)  ;
MtrxMultMatrx(arrJ1,3, 3, arrVS_PSK, 1, arrVS_KGSK) ;
}

// матрица перехода parrMtrxPer из исходной прямоугольной сиситемы координат
// в скоростную сиситему координат (CCК), задаваемую вектором arrV
// ось X  ССК направлена по вектору V, ось Y CCК праллельна плоскости OXY
// исходной прямоугольной сиситемы координат, а ось Z дополняет до правой тройки
void calcMatrxPer_from_DecartPrSK_To_SSK(double *arrV, double * parrMtrxPer)
{
	double val_v0 = NormVect2(arrV);
	memset(parrMtrxPer, 0, 9 * sizeof(double));
	if (val_v0 < 0.0000000001)
  {
	parrMtrxPer[2] = 1.;
	parrMtrxPer[3] = 1.;
	parrMtrxPer[7] = 1.;
	return;
  }
  double arr_v[3] = {0.};
  memcpy(arr_v, arrV, 3 *  sizeof(double));
	NormalizeVect3(arr_v) ;
	double val_v = NormVect2(arr_v);
  parrMtrxPer[0] = arr_v [0];
  parrMtrxPer[1] = arr_v [1];
	parrMtrxPer[2] = arr_v[2];
	parrMtrxPer[3] = -arr_v [1] / val_v;
	parrMtrxPer[4] = arr_v [0] / val_v;
  parrMtrxPer[5] = 0.;
	parrMtrxPer[6] = -arr_v [2] * arr_v [0] / val_v;
	parrMtrxPer[7] = -arr_v [2] * arr_v [1] / val_v;
  parrMtrxPer[8] =  val_v ;
  return;
}


// матрица перехода parrMtrxPer из скоростной  сиситемы координат (CCК), задаваемой вектором arrV
// в исходную  прямоугольную  сиситему координат
// ось X  ССК направлена по вектору V, ось Y CCК праллельна плоскости OXY
// исходной прямоугольной сиситемы координат, а ось Z дополняет до правой тройки
// parrMtrxPer - это матрица состоящая из столбцов вектров единичных ортов
// CCК  в исходной прямоуг сиситемы координат
void calcMatrxPer_from_SSK_To_DecartPrSK(double *arrV, double * parrMtrxPer)
{
  double arrT[9] = {0.};
  calcMatrxPer_from_DecartPrSK_To_SSK(arrV, arrT) ;
  MatrTransp(arrT, 3, 3, parrMtrxPer);
  return;
}

#pragma package(smart_init)
