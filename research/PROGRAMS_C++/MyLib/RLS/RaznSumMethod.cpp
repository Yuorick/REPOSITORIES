//---------------------------------------------------------------------------


#pragma hdrstop
#include "Comp.h"
#include "RaznSumMethod.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "URPolyLine.h"
#include "UrPointXY.h"
#include "Gauss.h"

 TRaznSumMethod::~TRaznSumMethod()
{

 if(mpiarrNum)
 {
	 free (mpiarrNum);
	 mpiarrNum = NULL;
 }
}

TRaznSumMethod::TRaznSumMethod()
{
 // количество модулей (фасет)
  m_NRows = 4;
 // расстояние между модулями
   m_D = 32.8;
 // длина волны
   mLambda = 3.15;;
 // массив фасет
  mFaceta = TFaceta();

  mpiarrNum = (int*)malloc(m_NRows * sizeof(int));
  mpiarrNum [0] = 6;
  mpiarrNum [1] = 8;
  mpiarrNum [2] = 8;
  mpiarrNum [3] = 6;

  mparrDisp =  (double*)malloc(m_NRows * sizeof(double));
  mparrDisp[0] = 400./ 28.;
  mparrDisp[1] = 400./ 28.;
  mparrDisp[2] = 400./ 28.;
  mparrDisp[3] = 400./ 28.;

}
// Конструктор копирования
TRaznSumMethod::TRaznSumMethod (const TRaznSumMethod &R2)
 {

// количество модулей (фасет)
	m_NRows = R2.m_NRows;
	// расстояние между модулями
	m_D = R2.m_D;
	// длина волны
	mLambda = R2.mLambda;

	mFaceta = R2.mFaceta;
	if (mpiarrNum)
	{
   //	free (mpiarrNum);
	mpiarrNum = NULL;
	}
	if (R2.mpiarrNum)
	{
	mpiarrNum = (int*)malloc(m_NRows * sizeof(int));
	memcpy(mpiarrNum, R2.mpiarrNum,m_NRows * sizeof(int));
	}

	if (mparrDisp)
	{
   //	free (mparrDisp);
	mparrDisp = NULL;
	}
	if (R2.mparrDisp)
	{
	mparrDisp = (double*)malloc(m_NRows * sizeof(double));
	memcpy(mparrDisp, R2.mparrDisp,m_NRows * sizeof(double));
	}
 }



 // оператор присваивания
 TRaznSumMethod TRaznSumMethod::operator=(TRaznSumMethod  R2)
{
	// количество модулей (фасет)
	m_NRows = R2.m_NRows;
	// расстояние между модулями
	m_D = R2.m_D;
	// длина волны
	mLambda = R2.mLambda;

	mFaceta = R2.mFaceta;
	if (mpiarrNum)
	{
	free (mpiarrNum);
	mpiarrNum = NULL;
	}
	if (R2.mpiarrNum)
	{
	mpiarrNum = (int*)malloc(m_NRows * sizeof(int));
	memcpy(mpiarrNum, R2.mpiarrNum,m_NRows * sizeof(int));
	}
	if (mparrDisp)
	{
	free (mparrDisp);
	mparrDisp = NULL;
	}
	if (R2.mparrDisp)
	{
	mparrDisp = (double*)malloc(m_NRows * sizeof(double));
	memcpy(mparrDisp, R2.mparrDisp,m_NRows * sizeof(double));
	}


  return *this ;
}

// парам констр
  TRaznSumMethod::TRaznSumMethod(const int N,const double D,const double Lambda
   ,TFaceta Faceta, int *piarrNum, double *parrDisp)
 {
	 m_NRows = N ;
	 m_D= D;
	 mLambda = Lambda ;
	 mFaceta = Faceta;
	 if (piarrNum)
	 {
		mpiarrNum =  (int *)malloc( N * sizeof(int));
		memcpy(mpiarrNum, piarrNum, N * sizeof(int));
	 }

	 if (parrDisp)
	 {
		mparrDisp =  (double *)malloc( N * sizeof(double));
		memcpy(mparrDisp, parrDisp, N * sizeof(double));
	 }

 }

// пеленгационная функция СРМ
// INPUT:
// pcmpSZv - массив замеров по строкам, комплексный
//  OUTPUT -  массив  частных производных пеленгац ф-ции по строковым диаграммам
// pcarr_dPel_po_dSi
 double TRaznSumMethod::fncPeleng(TComp *pcmpSZv, TComp *pcarr_dPel_po_dSi)
 {
 // TComp cmpUp =  (pcmpSZv[0] + pcmpSZv[1]);
 // TComp cmpDown =  (pcmpSZv[2] + pcmpSZv[3]);
  TComp cmpUp(0.,0.);// =  (pcmpSZv[0] + pcmpSZv[1]);
  TComp cmpDown (0.,0.);//=  (pcmpSZv[2] + pcmpSZv[3]);
  for (int i = 0; i < m_NRows/2; i++)
  {
	cmpUp += pcmpSZv[i];
	cmpDown += pcmpSZv[m_NRows/2 + i];
  }
  TComp cmpRazn = cmpUp -   cmpDown;
  TComp cmpSum  = cmpUp +   cmpDown;
  TComp cmpPel =  cmpRazn/cmpSum;

   if(pcarr_dPel_po_dSi)
   {
	TComp cmp2(2.,0.);
	TComp cmpMin2(-2.,0.);
	for (int i = 0; i < m_NRows/2; i++)
	{
	pcarr_dPel_po_dSi[           i] = cmp2 *cmpDown / cmpSum/cmpSum ;
	pcarr_dPel_po_dSi[m_NRows/2 +i] = cmpMin2 *cmpUp / cmpSum/cmpSum ;
	}
  }

  return cmpPel.m_Im;
 }

double TRaznSumMethod::calcDispPelengFnc(TComp *pcmpSZv)
{
  TComp *pcarr_dPel_po_dSi = (TComp*)malloc(m_NRows * sizeof(TComp));
  double valPel = fncPeleng(pcmpSZv, pcarr_dPel_po_dSi);
  double rez = 0.;
  for (int i = 0; i < m_NRows; i++)
  {
	rez +=  pcarr_dPel_po_dSi[i].modul() * pcarr_dPel_po_dSi[i].modul() * mpiarrNum[i] * mparrDisp[i] ;
  }

  free(pcarr_dPel_po_dSi);
  return rez;

}
// формирование вектора диаграмм без шумов для зхаданного угла цели
// INPUT:
// ValAng  - угол цели, рад
// cmpKTarg - коэфф отражения цели 
// OUTPUT:
// pcmpS - вектор измерений диаграмм
void  TRaznSumMethod::imitateVectorIdealDiagrams(const double ValTetta, TComp cmpKTarg, TComp *pcmpS)
{
  int iNumAM = 0; // общее к-во АМ 
  for (int i = 0; i < m_NRows; i++)
  {
   iNumAM += mpiarrNum[i];   
  }
  ///
  TComp cmpModulAmp ( mFaceta.fncFFacetaApprox(ValTetta), 0.) ;
  double valMu = m_D * sin (ValTetta)/ mLambda  * 2. * M_PI ;
   for (int i = 0; i < m_NRows; i++)
   {
	 TComp cmpRotate = exp_(valMu * ((double)i));
	 TComp cmpCoef(((double)mpiarrNum[i])/ ((double)iNumAM), 0.);
	 pcmpS[i] = cmpModulAmp * cmpCoef * cmpRotate *cmpKTarg;
   }
  
}

// формирование вектора измерений диаграмм с шумами для заданного угла цели
// INPUT:
// ValAng  - угол цели, рад
// OUTPUT:
// pcmpSZv - вектор измерений диаграмм
void  TRaznSumMethod::imitateVectorMeasures(const double ValTetta, TComp cmpKTarg, TComp *pcmpS, TComp *pcmpSZv)
{
  
  imitateVectorIdealDiagrams( ValTetta, cmpKTarg, pcmpS) ;
   for (int i = 0; i < m_NRows; i++)
   { 	
	 pcmpSZv[i].m_Re = pcmpS[i].m_Re +   getGauss(0., sqrt(mparrDisp[i] *((double)mpiarrNum[i]/ 2.) ));
	 pcmpSZv[i].m_Im = pcmpS[i].m_Im +   getGauss(0., sqrt(mparrDisp[i] *((double)mpiarrNum[i]/ 2.) ));

   }

}

void TRaznSumMethod::createIdealPelengFncGraph(wchar_t *wchFileName, const double ValRang )
{
  const double VAL_STEP = 1. * M_PI/ 3000./10.;
  const int I_NUM = ValRang / VAL_STEP * 2 +1;;
  // заготовка полилиний для графиков
  TURPolyLine plnBearing(1,   I_NUM ) ;

  TComp *pcmpS =  (TComp*) malloc(m_NRows * sizeof(TComp ));
  TComp cmpKTarg(800.,0.);
  for (int i= 0; i < I_NUM ; i++)
  {

	double valAngMin = -((double)((I_NUM-1)/2) ) * VAL_STEP;
	double valTet = valAngMin + ((double)i)* VAL_STEP;
	imitateVectorIdealDiagrams(valTet,cmpKTarg,  pcmpS);
   //	plnBearing.Points[i].Y = fncPeleng(pcmpS, NULL)*10.;
   //	plnBearing.Points[i].X =  valTet *3000./M_PI ;
	plnBearing.Points[i].X = fncPeleng(pcmpS, NULL)*10.;
	plnBearing.Points[i].Y =  valTet *3000./M_PI ;

  }

  plnBearing.WriteSetSHPFiles(wchFileName,&plnBearing, 1) ;

  free (pcmpS);

}


void TRaznSumMethod::createSKZGraph(wchar_t *wchFileName, const double ValRang )
{
  TRaznSumMethod rsmImproved =  correctFirstAndLastRows();
  const double VAL_STEP = 1. * M_PI/ 3000./10.;
  const int I_NUM = ValRang / VAL_STEP * 2 +1;;
  // заготовка полилиний для графиков
  TURPolyLine plnSKZ(1,   I_NUM ) ;

  TComp *pcmpS =  (TComp*) malloc(m_NRows * sizeof(TComp ));
  TComp cmpKTarg(800.,0.);
  for (int i= 0; i < I_NUM ; i++)
  {

	double valAngMin = -((double)((I_NUM-1)/2) ) * VAL_STEP;
	double valTet = valAngMin + ((double)i)* VAL_STEP;
	imitateVectorIdealDiagrams(valTet,cmpKTarg,  pcmpS);
	double valDisp = 0.;
	double valTetZv = fncEstimateRsmTetta(pcmpS,  &valDisp ) ;
	double valTetZv1 = rsmImproved.fncEstimateRsmTetta(pcmpS,  &valDisp ) ;

	plnSKZ.Points[i].Y = sqrt(valDisp +(valTetZv - valTet) * (valTetZv - valTet))* 3000./ M_PI *10.;
	plnSKZ.Points[i].X =  valTet *3000./M_PI ;

  }

  plnSKZ.WriteSetSHPFiles(wchFileName,&plnSKZ, 1) ;

  free (pcmpS);

}


TRaznSumMethod TRaznSumMethod::correctFirstAndLastRows()
{								  
  TRaznSumMethod rsmRez = *this;   
  rsmRez.mparrDisp[0] = rsmRez.mparrDisp[0] * ((double)rsmRez.mpiarrNum[1] ) /((double)rsmRez.mpiarrNum[0] )
	   * ((double)rsmRez.mpiarrNum[1] ) /((double)rsmRez.mpiarrNum[0] ) ;
  rsmRez.mparrDisp[rsmRez.m_NRows -1] = rsmRez.mparrDisp[0];

  rsmRez.mpiarrNum[0] =  rsmRez.mpiarrNum[1];
  rsmRez.mpiarrNum[m_NRows -1] =  rsmRez.mpiarrNum[1];
  
  return rsmRez;
}

// вычисление оценки угла РСМ и расчет дисперсии ошибки для пеленг 
// функции tg(y/2)
//  INPUT:
 // pcmpSZv - массив замеров по строкам, комплексный
 // OUTPUT:
 // pDisp - дисперсия ошибки оценивания
double TRaznSumMethod::fncEstimateRsmTetta(TComp *pcmpSZv,  double *pDisp )
{
	double valPelFnc = fncPeleng(pcmpSZv, NULL);
	double valDispPelFnc = calcDispPelengFnc(pcmpSZv);
	double valCoef =  mLambda/( M_PI *m_D * ((double)m_NRows)/ 2.);
   	double valTet = -atan(valPelFnc)* valCoef;
	double valDeriv = 1./ (1. + valPelFnc * valPelFnc );
	*pDisp =  valDispPelFnc * valDeriv * valDeriv * valCoef* valCoef;
	return valTet;
}

// вычисление оценки угла места РСМ и расчет дисперсии по аналитической формуле для
// скомпенсированной диаграммы по верхней и нижней строке
// функции tg(y/2)
//  INPUT:
 // pcmpSZv - массив замеров по строкам, комплексный
 // OUTPUT:
 // pDisp - дисперсия ошибки оценивания
 //
 //Attention!! 4 строковые диаграммы   изделие 5П10 только!!!!
double TRaznSumMethod::fncEstimateRsmUM_with_ApproxDisp(TComp *pcmpSZv,  double *pDisp )
{ 	
	double valNoiseSKZ = sqrt(calcTotalDisp());
	double valPelFnc = 0.;
	double valDispPelFnc = 0.;
	fncBearingRSM( valNoiseSKZ, pcmpSZv, &valPelFnc, &valDispPelFnc) ;
	double valCoef =  mLambda/( 2. * M_PI *m_D);
	double valTet = - atan(valPelFnc)* valCoef;
	
	double valDeriv = 1./ (1. + valPelFnc * valPelFnc );
	*pDisp =  valDispPelFnc * valDeriv * valDeriv * valCoef* valCoef;

	return valTet;
}

// общее число АМ
int TRaznSumMethod::calcQuantAM()
{
  int iNumAM = 0; // общее к-во АМ
  for (int i = 0; i < m_NRows; i++)
  {
   iNumAM += mpiarrNum[i];   
  }
 return iNumAM ;
}

// суммарная дисперсия
double TRaznSumMethod::calcTotalDisp()
{
  double valTotalDisp = 0; // 
  for (int i = 0; i < m_NRows; i++)
  {
   valTotalDisp  += ((double)mpiarrNum[i]) * mparrDisp[i];   
  }
 return valTotalDisp  ;
}

double TRaznSumMethod::clcSystErr(const double ValTetta)
{
	TComp *pcmpS =  (TComp*) malloc(m_NRows * sizeof(TComp ));
	TComp cmpKTarg (800.,0.);
	imitateVectorIdealDiagrams( ValTetta,cmpKTarg,  pcmpS);
	double pDisp = 0.;
	double valTettaEst = fncEstimateRsmTetta(pcmpS,  &pDisp )  ;
	free(pcmpS);
	return (ValTetta-  valTettaEst);
	return 0.;
}
 // пеленгационная функция разностно суммарного метода и ее дисперсия
 //  INPUT:
 // valNoiseSKZ - скз шума
 // pcmpSZv - массив замеров по строкам, комплексный
 // OUTPUT:
// pvalPelFunc   - величина пеленгационной функции
// pvalDispPelFunc  - дисперсия пеленгационной функции
void fncBearingRSM(double valNoiseSKZ, TComp *pcmpSZv, double *pvalPelFunc, double *pvalDispPelFunc)
{
  TComp cmpUp =  (pcmpSZv[0] + pcmpSZv[1]);
  TComp cmpDown =  (pcmpSZv[2] + pcmpSZv[3]);
  TComp cmpRazn = cmpUp -   cmpDown;
  TComp cmpSum  = cmpUp +   cmpDown;
  TComp cmpPel =  cmpRazn/cmpSum;
  *pvalPelFunc = cmpPel.m_Im;
  *pvalDispPelFunc =  valNoiseSKZ * valNoiseSKZ * 4./3.*(cmpUp.modul()*cmpUp.modul()
	 + cmpDown.modul() * cmpDown.modul())/ cmpSum.modul()/ cmpSum.modul()/ cmpSum.modul()/ cmpSum.modul() ;

}


// вычисление функции хи квадрат для изделия 5П10
//по углу места
// INPUT:
// valNoiseSKZ - шум в суммарной диаграмме
// pcmpSZv - замеры строковых диаграмм, скомпенсированные по модулю
//valMuZv - величина Mu = D* sin(tet)*2* M_PI / lambda
//
double  fncXi2_5P10_UM(double valNoiseSKZ, TComp *pcmpSZv, double valMuZv)
{
	TComp cmparrQZv[4];
	for (int i=0; i < 4; i++)
	{
	  TComp cmpRotate = exp_(-((double)i) * valMuZv);
	  cmparrQZv[i] = pcmpSZv[i] * cmpRotate;
	}
	TComp cmp0_75(3./4.,0);
	TComp cmpB =  cmparrQZv[0] +  cmparrQZv[3];
	cmpB *= cmp0_75;
	cmpB  += cmparrQZv[1];
	cmpB  += cmparrQZv[2];

	double temp0 = 21. *(cmparrQZv[0].modul()*cmparrQZv[0].modul() + cmparrQZv[3].modul()*cmparrQZv[3].modul())/8. ;
	double temp1 = 7. *(cmparrQZv[1].modul()*cmparrQZv[1].modul() + cmparrQZv[2].modul()*cmparrQZv[2].modul())/2. ;

	double temp2 = temp0 + temp1 -   cmpB.modul()*cmpB.modul();
	double temp3 = temp2/valNoiseSKZ/valNoiseSKZ ;

	// альтернативный рассчет
	double sig0Sq = valNoiseSKZ*valNoiseSKZ * 8./ 21.;
	double sig1Sq = valNoiseSKZ*valNoiseSKZ * 2./ 7.;
	double sig0Sqmin1 = 1./ sig0Sq;
	TComp cmpsig0Sqmin1(sig0Sqmin1,0.);
	double sig1Sqmin1 = 1./ sig1Sq;
	TComp cmpsig1Sqmin1(sig1Sqmin1,0.);
	double temp4 = sig0Sqmin1 *(cmparrQZv[0].modul()*cmparrQZv[0].modul() + cmparrQZv[3].modul()*cmparrQZv[3].modul()) ;
	double temp5 = sig1Sqmin1 *(cmparrQZv[1].modul()*cmparrQZv[1].modul() + cmparrQZv[2].modul()*cmparrQZv[2].modul()) ;

	TComp cmpsigSqSum ( 1./(sig0Sqmin1 +  sig1Sqmin1),0.);
	TComp cmpT1 = cmparrQZv[0] +  cmparrQZv[3];
	TComp cmpT2 = cmpT1 *cmpsig0Sqmin1;

	TComp cmpT3 =  cmparrQZv[1] +  cmparrQZv[2];
	TComp cmpT4 =  cmpT3 *cmpsig1Sqmin1;

	TComp cmpT5 = cmpT2 + cmpT4;

	double temp6 =   cmpT5.modul()*cmpT5.modul()/ 2./(sig0Sqmin1 +sig1Sqmin1);

	double temp7 =  temp4 + temp5 - temp6;



	return temp3 ;
}

#pragma package(smart_init)
