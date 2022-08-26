//---------------------------------------------------------------------------


#pragma hdrstop

#include "Veer41.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "YrWriteShapeFile.h"
#include "URPolyLine.h"
#include "UrPointXY.h"
#include "DiagrSinX.h"
#include "MatrixProccess.h"
#include "Comp.h"
#include "Diagr.h"
#include "Far.h"

_fastcall  TVeer41 ::TVeer41()
{
  marrDgr[0] = TDiagr();
  marrDgr[1] = TDiagr();

   mTang = 60.;;
 // угол сканирования разностной диаграммы
  marrScnDif[0] = -5. *  M_PI/ 3000.;
  marrScnDif[1] = 5. * M_PI/ 3000.;


}
// Конструктор копирования
__fastcall  TVeer41::TVeer41 (const TVeer41 &R)
 {
	marrDgr[0] = R.marrDgr[0] ;
	marrDgr[1] = R.marrDgr[1] ;
	mTang =  R.mTang;
	marrScnDif[0] = R.marrScnDif[0] ;
	marrScnDif[1] = R.marrScnDif[1] ;

 }
// оператор присваивания
  TVeer41 TVeer41::operator=(TVeer41  R)
 {
	marrDgr[0] = R.marrDgr[0] ;
	marrDgr[1] = R.marrDgr[1] ;
	mTang =  R.mTang;
	marrScnDif[0] = R.marrScnDif[0] ;
	marrScnDif[1] = R.marrScnDif[1] ;
   return *this ;
 }

 // парам констр
 __fastcall TVeer41::TVeer41(const double Tang,  double *arrScnDif, const double ang )
 {
	TFar Far;
	int iTypeOfDiagram = 0;
	double parrPar[2] = {0.005, 0.};

	parrPar[0] = 0.;
	parrPar[1] =  ang;
	marrDgr[0] = TDiagr( Far, iTypeOfDiagram, parrPar) ;
	parrPar[1] = - ang;
	marrDgr[1] = TDiagr( Far, iTypeOfDiagram, parrPar) ;

	mTang = Tang ;
	marrScnDif[0] = arrScnDif[0] ;
	marrScnDif[1] = arrScnDif[1] ;
}

//--------------------------------------------------------------------------------------
// Нахождение оптимальных углов цели и антипода
// метод работает с измерениями по 3-м диаграммам
// 2 разностные и одна комплексная( первая по счету)
// INPUT:
// pSDifZv[2]   - массив измерений по разностным диаграммам (комплексные числа)
// cmpSZv   - измерение по первой диаграмме (комплексное число)
// ZTarg -  коэффициент отражения цели нач значение
// ZAnt   - коэффициент отражения  антипода  нач значение
// *alfTrg  - УМ цели   нач значение
// *alfAnp  - УМ антипода  нач значение
// OUTPUT:
// ZTarg -  коэффициент отражения цели
// ZAnt   - коэффициент отражения  антипода
// *alfTrg  - УМ цели
// *alfAnp  - УМ антипода
// возвращает:
// -2 - если матрица якоби вырождена
// -3 - Если метод не сошелся

int _fastcall TVeer41::solv3(TComp *pSDifZv,TComp *pcmpSZv , TComp *pZTarg, TComp *pZAnt, double *palfTrg, double *palfAnp )
{

  int breturn = -3;
  int i = 0;
  double arr_FGreek[2] ={0.},arr_dFGreek[4] ={0.},arr_dFGreekInv[4] ={0.} ;
  double arrX[2] ={0.}, arrXT[2] ={0.}; // вектор с решениями
  arrX[0] = *palfTrg ;
  arrX[1] = *palfAnp;
  double del = -2.;

  for (i = 0; i < 1000; i++)
  {

	double arr_F00[2] ={0.},arr_dF00[4] ={0.}, arrDelX[2] ={0.}, arrDelX1[2] ={0.} ;

	calc_vectG_and_mtrxH_solv3 (pSDifZv,pcmpSZv
		,  arrX[0], arrX[1], arr_FGreek,  arr_dFGreek ) ;
	  //	double parrRez[4] ={0.}, parrRez1[4] ={0.};
	   //	MatrxMultScalar(arr_dFGreek, 2, 2, 0.000001,parrRez);
		// проверка положит определ для отлажки!!!
	   //	double valPos = arr_dFGreek[0] *arr_dFGreek[3] - arr_dFGreek[1] * arr_dFGreek[2];
	///

   //	if(!InverseMtrx2(parrRez, parrRez1))
	if(!InverseMtrx2(arr_dFGreek, arr_dFGreekInv))
	{
	 *palfTrg = arrX[0]   ;
	 *palfAnp = arrX[1]  ;
	 return -2;
	}
   //	MatrxMultScalar(parrRez1, 2, 2, 1000000.,arr_dFGreekInv);
	MtrxMultMatrx(arr_dFGreekInv ,2, 2, arr_FGreek,1, arrDelX) ;
	del = NormVect(arrDelX, 2);
	double arrT[2] ={0.};
	MatrxMultScalar(arrDelX, 2, 1, 0.1,arrDelX1);
	memcpy( arrDelX, arrDelX1, 2 * sizeof(double));
	MtrxMinusMatrx(arrX, arrDelX,1, 2, arrXT) ;
	memcpy( arrX, arrXT, 2 * sizeof(double));
	if (del< 0.0000001)
	{

	  *palfTrg = arrX[0]   ;
	  *palfAnp = arrX[1]  ;
	   breturn= 0;
	  break;
	}
 }
 double x1,x2,x3,x4,arr_gradK11[2],  arr_gradK12[2],  arr_gradK21[2]
	 , arr_gradK22[2], arr_HessK11[4],  arr_HessK12[4],  arr_HessK21[4], arr_HessK22[4];
 calc_CoeffOtrag(pSDifZv, *palfTrg, *palfAnp, &x1,  &x2
	 , &x3, &x4, arr_gradK11,  arr_gradK12,  arr_gradK21
	 , arr_gradK22, arr_HessK11,  arr_HessK12,  arr_HessK21, arr_HessK22);
  *pZTarg =TComp(x1,x2);
  *pZAnt =TComp(x3,x4);
  return breturn;

}


int _fastcall TVeer41::calc_vectG_and_mtrxH_solv3 (TComp *pSDifZv, TComp *pcmpSZv
		,  double alfTrg, double alfAnp, double* arr_FGreek, double* arr_dFGreek )
{
	double arr_Part_F[2] ={0.},arr_Part_dF[4] ={0.}, arrT[2]={0.}, arrT1[4]={0.};
	memset(arr_FGreek, 0, 2 * sizeof(double));
	memset(arr_dFGreek , 0, 4 * sizeof(double));

	double valK11 = 0., valK12 =0., valK21 =0., valK22 =0., arr_gradK11[2] ={0.}, arr_gradK12[2] ={0.}, arr_gradK21[2] ={0.},
	  arr_gradK22[2] ={0.}, arr_HessK11[4] ={0.}, arr_HessK12[4] ={0.}, arr_HessK21[4] ={0.},  arr_HessK22 [4] ={0.};

	calc_CoeffOtrag(pSDifZv, alfTrg, alfAnp, &valK11, &valK12, &valK21, &valK22, arr_gradK11
	, arr_gradK12, arr_gradK21, arr_gradK22, arr_HessK11, arr_HessK12
	, arr_HessK21, arr_HessK22);

	marrDgr[0].calcPartial_vectA_and_mtrx_dA_po_dx_solv3 (pcmpSZv[0],alfTrg, alfAnp, valK11, valK12, valK21, valK22, arr_gradK11
	, arr_gradK12, arr_gradK21, arr_gradK22, arr_HessK11, arr_HessK12
	, arr_HessK21, arr_HessK22,  arr_Part_F,  arr_Part_dF ) ;

	MtrxSumMatrx(arr_FGreek, arr_Part_F,1, 2, arrT) ;

	memcpy(arr_FGreek,arrT, 2 * sizeof(double));

	MtrxSumMatrx(arr_dFGreek, arr_Part_dF,2, 2, arrT1) ;

	memcpy(arr_dFGreek,arrT1, 4 * sizeof(double));


	return 0;
}

void _fastcall TVeer41::calc_CoeffOtrag(TComp *pSDifZv,double alfTrg,double  alfAnp,double* valK11, double* valK12
	 ,double* valK21,double* valK22,double* arr_gradK11, double* arr_gradK12, double* arr_gradK21
	 ,double* arr_gradK22, double*arr_HessK11, double* arr_HessK12, double* arr_HessK21,double* arr_HessK22)
{
   // преобразование массива измерений
   TComp qZv[2], cmpTang(mTang, 0.);
   TComp cmpi(0.,1.);
   TComp cmpZero(0.,0.);
  // TComp cmpC0( marrDgr[0].mKSum, 0.);
 //  TComp cmpC1( marrDgr[1].mKSum, 0.);
   qZv[0] = cmpZero - (cmpi * pSDifZv[0] /cmpTang) ;
   qZv[1] = cmpZero - (cmpi * pSDifZv[1] /cmpTang) ;
   ///

  // вычисление вспомогательных функций fi1, fi2, grad_f1, grad_f2, hess_fi1, hess_fi2
  TComp cmpX1(alfTrg, 0.);
  TComp cmpX2(alfAnp, 0.);
  TComp cmpFi1, cmpFi2, cmpArrGradFi1[2], cmpArrGradFi2[2],cmpArrHessFi1[4], cmpArrHessFi2[4];
  calcFncs_Fi1_and_Fi2( cmpX1, cmpX2, cmpFi1, cmpFi2);
  calcVectrsGrad_Fi1_and_Fi2( cmpX1, cmpX2, cmpArrGradFi1, cmpArrGradFi2);
  calcMtrxHess_Fi1_and_Fi2( cmpX1, cmpX2, cmpArrHessFi1, cmpArrHessFi2);
  ///
  // вычисление z1, z2, grad_z1, grad_z2, hess_z1, hess_z2
  TComp cmpAlf1(marrScnDif[0], 0.), cmpAlf2(marrScnDif[1], 0.);
  TComp cA1 = ( qZv[0] - qZv[1]) / (cmpAlf2 - cmpAlf1) ;
  TComp cA2 = ( cmpAlf1 * qZv[1] -  cmpAlf2 *  qZv[0]) / (cmpAlf2 - cmpAlf1) ;
  TComp cZ1 = cA1 * cmpFi1 +  cA2 * cmpFi2;
  TComp cZ2 = cA1 - cZ1;
  TComp  grad_z1[2], grad_z2[2], hess_z1[4], hess_z2[4], carrT0[2], carrT1[2]
	, carrT2[4], carrT3[4];
  MatrxMultScalar(cmpArrGradFi1, 2, 1, cA1,carrT0);
  MatrxMultScalar(cmpArrGradFi2, 2, 1, cA2,carrT1);
  MtrxSumMatrx(carrT0, carrT1, 2, 1, grad_z1) ;
  TComp cmpMinus1(-1.,0.);
  MatrxMultScalar(grad_z1, 2, 1, cmpMinus1,grad_z2);

  MatrxMultScalar(cmpArrHessFi1, 4, 1, cA1,carrT2);
  MatrxMultScalar(cmpArrHessFi2, 4, 1, cA2,carrT3);
  MtrxSumMatrx(carrT2, carrT3, 4, 1, hess_z1) ;
  MatrxMultScalar( hess_z1, 4, 1, cmpMinus1, hess_z2);
  ///

  *valK11 =  cZ1.m_Re;
  *valK12 =  cZ1.m_Im;
  *valK21 =  cZ2.m_Re;
  *valK22 =  cZ2.m_Im;

  MatrRe(grad_z1, 1, 2, arr_gradK11);
  MatrIm(grad_z1, 1, 2, arr_gradK12);
  MatrRe(grad_z2, 1, 2, arr_gradK21);
  MatrIm(grad_z2, 1, 2, arr_gradK22);

   MatrRe(hess_z1, 1, 4, arr_HessK11);
   MatrIm(hess_z1, 1, 4, arr_HessK12);
   MatrRe(hess_z2, 1, 4, arr_HessK21);
   MatrIm(hess_z2, 1, 4, arr_HessK22);


}

void _fastcall TVeer41::calcFncs_Fi1_and_Fi2( TComp  cmpX1, TComp cmpX2,TComp& cmpFi1,TComp & cmpFi2)
{
  TComp  temp =  cmpX2 - cmpX1;
  cmpFi1 =   cmpX2 / ( cmpX2 - cmpX1);
  TComp  cmp_1(1., 0.);
  cmpFi2 =  cmp_1 / ( cmpX2 - cmpX1);
}
void _fastcall TVeer41::calcVectrsGrad_Fi1_and_Fi2( TComp   cmpX1, TComp  cmpX2,TComp* cmpArrGradFi1, TComp* cmpArrGradFi2)
{
  TComp carrT0[2];
  TComp cmpNOLL(0.,0.);
  TComp  cmp_1(1., 0.);
  carrT0[0] =  cmpX2;
  carrT0[1] = cmpNOLL - cmpX1;
  TComp temp = cmp_1 / ((cmpX2 - cmpX1) * (cmpX2 - cmpX1));
  MatrxMultScalar(carrT0, 1, 2, temp,cmpArrGradFi1);
  carrT0[0] = cmp_1;
  carrT0[1] = cmpNOLL -cmp_1;
  MatrxMultScalar(carrT0, 1, 2, temp,cmpArrGradFi2);


}
void _fastcall TVeer41::calcMtrxHess_Fi1_and_Fi2( TComp  cmpX1, TComp  cmpX2
	   ,TComp* cmpArrHessFi1,TComp* cmpArrHessFi2)
{
  TComp carrT[4];
  TComp cmpNOLL(0.,0.);
  TComp  cmp_1(1., 0.);
  TComp  cmp_2(2., 0.);
  carrT[0] =  cmp_2 * cmpX2 ;
  carrT[1] = cmpNOLL - (cmpX1 + cmpX2);
  carrT[2] = carrT[1];
  carrT[3] = cmp_2 *cmpX1 ;
  TComp temp = cmp_1 / ((cmpX2 - cmpX1) * (cmpX2 - cmpX1)* (cmpX2 - cmpX1));
  MatrxMultScalar(carrT, 1, 4, temp,cmpArrHessFi1);

  carrT[0] =  cmp_2;
  carrT[1] = cmpNOLL - cmp_2;
  carrT[2] = carrT[1];
  carrT[3] =  cmp_2;
  MatrxMultScalar(carrT, 1, 4, temp,cmpArrHessFi2);

}

#pragma package(smart_init)
