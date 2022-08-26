//---------------------------------------------------------------------------


#pragma hdrstop
#include <math.h>
#include <string.h>
#include "Comp.h"
#include "QuadrEqMeth.h"
#include "Equations.h"
#include "MatrixProccess.h"
#include "Faceta.h"


 //INPUT:
 // pcmpSZv - массив измерений дианрамм, комплексный
 // OUTPUT:
 // pZTarg - комплексный коэфф  отражени€ цели
 // pZAnt -  комплексный коэфф  отражени€ антипода
 // * z1- цель
 // *z2 -  антипода
 // arrMtrxCorr_fi - коррел€ционна€ матирица
 int   solvQuadrEqMeth(TComp *pcmpSZv, TComp *pZTarg, TComp *pZAnt
	   , TComp * z1, TComp *z2 , double *arrMtrxCorr_fi )
{
  TComp cmpa(1., 0.);
  TComp cmpb = (pcmpSZv[1] * pcmpSZv[2] - pcmpSZv[0] * pcmpSZv[3] )/ ( pcmpSZv[0] * pcmpSZv[2] - pcmpSZv[1] * pcmpSZv[1] ) ;
  TComp cmpc = ( pcmpSZv[1] * pcmpSZv[3] - pcmpSZv[2] * pcmpSZv[2] )/ ( pcmpSZv[0] * pcmpSZv[2] - pcmpSZv[1] * pcmpSZv[1] ) ;
  SolvEq2( cmpa, cmpb,  cmpc, *z1, *z2);

  // вычисление коррел€ц матрицы ошибок аргументов z1, z2

  calcMtrxCorr_z1_z2(pcmpSZv, *z1, *z2,  arrMtrxCorr_fi )  ;


  return 0;
}


// ѕ≈–≈√–”∆≈ЌЌјя
//  pZTarg  и  pZAnt не нужны
//INPUT:
 // pcmpSZv - массив измерений дианрамм, комплексный
 // OUTPUT:
 // * z1- цель
 // *z2 -  антипода
 // arrMtrxCorr_fi - коррел€ционна€ матирица
 int   solvQuadrEqMeth(TComp *pcmpSZv, TComp * z1, TComp *z2 , double *arrMtrxCorr_fi )
{
  TComp cmpa(1., 0.);
  TComp cmpb = (pcmpSZv[1] * pcmpSZv[2] - pcmpSZv[0] * pcmpSZv[3] )/ ( pcmpSZv[0] * pcmpSZv[2] - pcmpSZv[1] * pcmpSZv[1] ) ;
  TComp cmpc = ( pcmpSZv[1] * pcmpSZv[3] - pcmpSZv[2] * pcmpSZv[2] )/ ( pcmpSZv[0] * pcmpSZv[2] - pcmpSZv[1] * pcmpSZv[1] ) ;
  SolvEq2( cmpa, cmpb,  cmpc, *z1, *z2);

  // вычисление коррел€ц матрицы ошибок аргументов z1, z2

  calcMtrxCorr_z1_z2(pcmpSZv, *z1, *z2,  arrMtrxCorr_fi )  ;


  return 0;
}

//-----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------
//----------- ѕ–ќ¬≈– ј ƒ–”√»ћ —ѕќ—ќЅќћ ------------------------------------------------------------------------
//-----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------
// вычисление коррел€ц матрицы ошибок z1, z2
void calcMtrxCorr_z1_z2(TComp *pcmpSZv, TComp z1, TComp z2,  double *arrMtrxCorr )
{
  TComp cmparr_k1[4], cmparr_k2[4];
  calcArr_k(pcmpSZv,  z1, cmparr_k1);
  calcArr_k(pcmpSZv,  z2, cmparr_k2);
  arrMtrxCorr[0] = calcPhaseDisp( pcmpSZv,  z1, cmparr_k1) ;
  arrMtrxCorr[3] = calcPhaseDisp( pcmpSZv,  z2, cmparr_k2);
  arrMtrxCorr[1] = calcSmeshMoment(pcmpSZv,   z1, cmparr_k1, z2, cmparr_k2);
  arrMtrxCorr[2] = arrMtrxCorr[1];
}

// дисперси€ фазового угла
double calcPhaseDisp( TComp *pcmpSZv, TComp z, TComp *cmparr_k)
{

	TComp cmpFig2(2., 0.);
	TComp cmpA = (pcmpSZv[0] * pcmpSZv[2] - pcmpSZv[1]* pcmpSZv[1]);
	TComp cmpB = (pcmpSZv[1] * pcmpSZv[2] - pcmpSZv[0]* pcmpSZv[3]);
	TComp cmpC = (pcmpSZv[1] * pcmpSZv[3] - pcmpSZv[2]* pcmpSZv[2]);
	TComp cmpFDeriv_Po_z = cmpFig2 * z * cmpA  + cmpB;
	double sum = 0.;
	for (int i = 0; i < 4; i++)
	{
	  sum +=  cmparr_k[i].modul() * cmparr_k[i].modul();
	  TComp cmpT0 =  cmparr_k[i]/ cmpFDeriv_Po_z ; ///////
	  int ii=0;

	}

	return sum / (cmpFDeriv_Po_z.modul() * cmpFDeriv_Po_z.modul() * z.modul()* z.modul());
}

TComp calcFDeriv_po_z( TComp *pcmpSZv, TComp z)
{
	TComp cmpFig2(2., 0.);
	TComp cmpA = (pcmpSZv[0] * pcmpSZv[2] - pcmpSZv[1]* pcmpSZv[1]);
	TComp cmpB = (pcmpSZv[1] * pcmpSZv[2] - pcmpSZv[0]* pcmpSZv[3]);
	TComp cmpC = (pcmpSZv[1] * pcmpSZv[3] - pcmpSZv[2]* pcmpSZv[2]);
	return ( cmpFig2 * z * cmpA  + cmpB);
}

// ¬ычисление массива коэффициентов k дл€ корн€ квадратного уравнени€
//INPUT:
//pcmpSZv[4] - массив измерений
//z - корень
// OUTPUT:
// cmparr_k[4] - массив коэффициентов
void   calcArr_k(TComp *pcmpSZv, TComp z, TComp *cmparr_k)
{
 TComp cmpA = (pcmpSZv[0] * pcmpSZv[2] - pcmpSZv[1]* pcmpSZv[1]);
 TComp cmpB = (pcmpSZv[1] * pcmpSZv[2] - pcmpSZv[0]* pcmpSZv[3]);
 TComp cmpFig1(1.,0.);
 TComp cmpFig2(2.,0.);
 TComp cmpMinFig2(-2.,0.);

 cmparr_k[0] =   pcmpSZv[2] * z* z - (pcmpSZv[3] * z);

 cmparr_k[1] = cmpMinFig2 * pcmpSZv[1] * z*z + pcmpSZv[2] * z +  pcmpSZv[3] ;

 cmparr_k[2] = pcmpSZv[0] * z* z + pcmpSZv[1] * z   - (cmpFig2 * pcmpSZv[2]);

 cmparr_k[3] =  pcmpSZv[1]  -  pcmpSZv[0] * z;

}

// сиешанный момоент ошибок фазовых углов
double calcSmeshMoment( TComp *pcmpSZv,  TComp z1,TComp* cmparr_k1,TComp z2, TComp* cmparr_k2)

{
	double arr_k1[4] = {0.},arr_k2[4] = {0.}, arrsum[4] ={0.}, arrT0[4] = {0.}, arrT1[4] = {0.};
	for (int i = 0; i < 4; i++)
	{
	 doArr(cmparr_k1[i], arr_k1);
	 doArr(cmparr_k2[i], arr_k2);
	 MtrxMultMatrxTransp(arr_k1,2, 2, arr_k2, 2, arrT0) ;
	 MtrxSumMatrx(arrT0, arrsum,2, 2, arrT1) ;
	 memcpy( arrsum, arrT1, 4 * sizeof(double));
	}

	double arrz1[2] ={0.},arrz2[2] ={0.}, arrFDeriv1[4] = {0.}, arrFDeriv2[4] = {0.}, arrFDerivT1[4] = {0.}, arrFDerivT2[4] = {0.};
	TComp cmpDeriv1 = calcFDeriv_po_z( pcmpSZv, z1);
	TComp cmpDeriv2 = calcFDeriv_po_z( pcmpSZv, z2);
	TComp cmpMult = z1 * z2 * cmpDeriv1 * cmpDeriv2;
	doArr(cmpDeriv1, arrFDerivT1);
	MatrxDivideScalar(arrFDerivT1, 2, 2, cmpMult.modul()  ,arrFDeriv1);
	doArr(cmpDeriv2, arrFDerivT2);
	MatrxDivideScalar(arrFDerivT2, 2, 2, cmpMult.modul()  ,arrFDeriv2);

	arrz1[0] = -z1.m_Im;
	arrz1[1] =  z1.m_Re;
	arrz2[0] = -z2.m_Im;
	arrz2[1] =  z2.m_Re;
	double arrT2[2] ={0.}, arrT3[2] ={0.}, arrT4[2] ={0.};
	MtrxMultMatrxTransp(arrz1, 1, 2, arrFDeriv1,2, arrT2) ;
	MtrxMultMatrx(arrT2, 1, 2, arrsum,2, arrT3) ;
	MtrxMultMatrx(arrT3, 1, 2, arrFDeriv2,2, arrT4) ;
	double rez =0.;
	MtrxMultMatrx(arrT4, 1, 2, arrz2,2, &rez) ;


	return rez;
}

void doArr(TComp cmpInp, double *arrRez)
{
  arrRez[0] =  cmpInp.m_Re;
  arrRez[1] = -cmpInp.m_Im ;
  arrRez[2] =  cmpInp.m_Im ;
  arrRez[3] =  cmpInp.m_Re;
}


//INPUT:
 // pcmpSZv - массив измерений дианрамм, комплексный
 // arrDisp[4] -  массив диспенрсий шума в строковых диаграммах
 //        arrDisp[i] - ƒиспероси€ шума в строке с номером i
 // valAMDist - рассто€ние между јћ
 // OUTPUT:
 // * valEstAngTarg - оценка угла цели
 // *valEstAngAntp -  оценка угла антипода
 // *cmpKTarg - оценка коеффиц отражени€ цели
 // *cmpKAntp - оценка коеффиц отражени€ артипода
 // arrMtrxCorr - коррел€ционна€ матирица ошибок оценивани€ углов
 int   fncEstimateMsd(double valAMDist, TFaceta Faceta, TComp *pcmpSZv, double * arrDisp, double *valEstAngTarg, double *valEstAngAntp
	  , TComp *cmpKTarg , TComp *cmpKAntp , double *arrMtrxCorr )
{
  TComp z1(0.,0.), z2(0.,0.);
  double arrMtrxCorr_fi[4]= {0.};
  int irez =   solvQuadrEqMeth(pcmpSZv, &z1, &z2 , arrMtrxCorr_fi );

	double ph1 = z1.phase();
	double ph2 = z2.phase();

   //	if ((ph1 > M_PI/ 4. )&& (fabs(ph2)< M_PI/ 4.  ) )
	if (ph1 > M_PI/ 4.* 1.2 )
	{
	ph1 -= M_PI * 2.;
	}

   //	if ((ph2 > M_PI/ 4.  )&& (fabs(ph1)< M_PI/ 4. ) )
	if (ph2 > M_PI/ 4.* 1.2 )
	{
	ph2 -= M_PI * 2.;
	}
	if( ph1< ph2)
	{
	TComp cmpt = z2;
	z2 = z1;
	z1 = cmpt;
	double temp = ph2 ;
	ph2 = ph1;
	ph1 = temp;

	temp = arrMtrxCorr_fi[3] ;
	arrMtrxCorr_fi[3] = arrMtrxCorr_fi[0] ;
	arrMtrxCorr_fi[0] = temp;

	}
	// оценки углов:
	const double VAL_WAVE_CONST = 2. * M_PI * valAMDist / Faceta.mLambda;

	*valEstAngTarg = asin( ph1/VAL_WAVE_CONST);
	*valEstAngAntp  = asin(ph2/VAL_WAVE_CONST);

	TComp cmpY1 = ( pcmpSZv[ 3] - pcmpSZv[2]* z2)/z1/z1/(z1-z2);
	TComp cmpY2 = (pcmpSZv[2] * z1 - pcmpSZv[3])/z2/z2/(z1-z2);

	TComp cmpFTarg ( Faceta.fncFFaceta (*valEstAngTarg), 0.);
	TComp cmpFAntp ( Faceta.fncFFaceta (*valEstAngAntp ), 0.);
	*cmpKTarg= cmpY1 / cmpFTarg;
	*cmpKAntp = cmpY2 / cmpFAntp ;

	//
	// пересчет корррел матрицы
	const double constdTarg = sqrt(1. - ph1 * ph1/ VAL_WAVE_CONST/VAL_WAVE_CONST );
	const double constdAntp = sqrt(1. - ph2 * ph2/ VAL_WAVE_CONST/VAL_WAVE_CONST );
	double valSumDispNoise = 0;
	for (int i = 0; i < 4; i++)
	{
	  valSumDispNoise +=  arrDisp[i] ;
	}
	arrMtrxCorr[0] = 1./ VAL_WAVE_CONST/ VAL_WAVE_CONST / constdTarg/ constdTarg * arrMtrxCorr_fi[0] *valSumDispNoise/2. ;
	arrMtrxCorr[3] = 1./ VAL_WAVE_CONST/ VAL_WAVE_CONST / constdAntp/ constdAntp * arrMtrxCorr_fi[3] *valSumDispNoise /2.;
	arrMtrxCorr[1] = 1./ VAL_WAVE_CONST/ VAL_WAVE_CONST / constdTarg/ constdAntp * arrMtrxCorr_fi[1] * valSumDispNoise  /2.;
	arrMtrxCorr[2] = arrMtrxCorr[1];
	double coefCor = arrMtrxCorr[1]/ sqrt(arrMtrxCorr[0] * arrMtrxCorr[3]);



  return 0;
}
#pragma package(smart_init)
