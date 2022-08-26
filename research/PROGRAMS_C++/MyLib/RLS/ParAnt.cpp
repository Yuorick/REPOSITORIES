//---------------------------------------------------------------------------

#pragma hdrstop
#include <float.h>
#include <vcl.h>
#include <dir.h>
#include "ParAnt.h"
#include <math.h>
#include "Comp.h"
#include <stdio.h>
#include <stdlib.h>
#include "YrWriteShapeFile.h"
#include "URPolyLine.h"
#include "UrPointXY.h"
#include "Equations.h"
#include "MatrixProccess.h"
#include "Gauss.h"
#include "Diagr.h"
#include "URPolygon.h"
#include "DiagrSinX.h"
//#include "TargRootsTab.h"
#include "Line.h"
#include "ArcEllipse.h"
#include "Hyperbola.h"
extern const double TET0707;
extern const double VAL_C ;
extern const double VAlCritAngSdvig =  2.928;
extern const double VAL_C = 299792458.;

TParAnt::TParAnt()
{
 // к-во парциальных диагшрамм в веере
  mQuantDiagr = 7;
 // апертура антенны
  mAppert = 120.;
 // длина волны
   mLambda = 3. ;
 // угол сдвига парциальных диаграмм  радианы
  mAngSdvig = 1.5026 *M_PI/ 180.;
 // дисперсиz шума в парциальных диаграммах
  mNoiseDisp = 400;

 // разброс коэффиц усиления k в парциальной диаграмме
  mAmplFactSig = 0.01;
}
// Конструктор копирования
TParAnt::TParAnt (const TParAnt &R2)
 {
// к-во парциальных диагшрамм в веере
  mQuantDiagr = R2.mQuantDiagr;
 // апертура антенны
  mAppert = R2.mAppert;
 // длина волны
   mLambda = R2. mLambda ;
 // угол сдвига парциальных диаграмм  радианы
  mAngSdvig = R2.mAngSdvig;
 // дисперсиz шума в парциальных диаграммах
  mNoiseDisp = R2.mNoiseDisp;

 // разброс коэффиц усиления k в парциальной диаграмме
  mAmplFactSig = R2. mAmplFactSig;
 }

  // оператор присваивания
  TParAnt TParAnt::operator=(TParAnt  R2)
{
// к-во парциальных диагшрамм в веере
  mQuantDiagr = R2.mQuantDiagr;
 // апертура антенны
  mAppert = R2.mAppert;
 // длина волны
   mLambda = R2. mLambda ;
 // угол сдвига парциальных диаграмм  радианы
  mAngSdvig = R2.mAngSdvig;
 // дисперсиz шума в парциальных диаграммах
  mNoiseDisp = R2.mNoiseDisp;

 // разброс коэффиц усиления k в парциальной диаграмме
  mAmplFactSig = R2.mAmplFactSig;

  return *this ;
}



// парам констр  1
 __fastcall TParAnt::TParAnt(const int quantDiagr,const double Appert,const double Lambda
   ,const double AngSdvig, const double NoiseDisp, const double AmplFactSig)
 {
	 mAppert = Appert ;
	 mLambda = Lambda ;
	 mNoiseDisp = NoiseDisp;
	 mAngSdvig = AngSdvig;
	 mAmplFactSig = AmplFactSig;
	 mQuantDiagr = quantDiagr;
 }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ширина парциальной диаграммы в обобщенном угле  по уровню 0,7
double TParAnt::findPartDiagrWidth07_GeneralizedAng()
{
	return  2.* TET0707;
}



// ширина парциальной диаграммы в обычном угле по уровню 0,7, рад
double TParAnt::findPartDiagrWidth07_Rad()
{
	return  2.* asin(TET0707 * mLambda/ mAppert/ M_PI);
}

// построение графиков 2-х парциальных диаграмм и суммарной в обобщенных угловых координатах и реальных
void   TParAnt::createGraphsPartial_and_Sum_Diagrams(wchar_t *wchFoldName1)
{
	const double  valGenAngSdvig =  transformAngToGeneralizedAng ( mAngSdvig );

	wchar_t wchFoldName[300] ={0};
	wcscpy(  wchFoldName,  wchFoldName1);
	wcscat(wchFoldName, L"\\");
	double step =0.0001;

	const double VAlDiap = 8. * M_PI;  // диапахзон графиков по обобщ координате плюс минус
	const int nBuffRows = 2. *VAlDiap/ step ;
	const int nBuffCols = 8;
	double  *parrBuff = new double [nBuffRows  * nBuffCols] ;
	memset(parrBuff, 0, nBuffRows * nBuffCols * sizeof(double));

	const int lenName = 30;
	wchar_t *wcharrFileNames  = new wchar_t [30 * nBuffCols];
	memset(wcharrFileNames, 0, 30 * nBuffCols * sizeof(wchar_t));

	wcscpy( wcharrFileNames, L"GeneralizedTetta");
	wcscpy( &wcharrFileNames[30], L"Tetta");
	wcscpy( &wcharrFileNames[2 * 30], L"PartialF");

	wcscpy( &wcharrFileNames[3 * 30], L"PartialF1");
	wcscpy( &wcharrFileNames[4 * 30], L"PartialF2");
	wcscpy( &wcharrFileNames[5 * 30], L"SumF");
	wcscpy( &wcharrFileNames[6 * 30], L"SumFSdvigMin");
	wcscpy( &wcharrFileNames[7 * 30], L"SumFSdvigPlus");

  double tet = -VAlDiap;
  for (int i=0 ; i < nBuffRows; i++)
  {
   tet = -VAlDiap + ((double)i) * step;
   parrBuff[ i * nBuffCols] = tet ; // обобщ угол
   parrBuff[ i * nBuffCols + 1] = asin ( tet * mLambda/ mAppert/ M_PI);
   parrBuff[ i * nBuffCols + 2] =  fncDiagrSinx_div_x( tet);
   parrBuff[ i * nBuffCols + 3] =  fncDiagrSinx_div_x( tet - valGenAngSdvig/2.);
   parrBuff[ i * nBuffCols + 4] =  fncDiagrSinx_div_x( tet + valGenAngSdvig/2.);
   parrBuff[ i * nBuffCols + 5] =  parrBuff[ i * nBuffCols + 3] + parrBuff[ i * nBuffCols + 4];
   parrBuff[ i * nBuffCols + 6] =  fncDiagrSinx_div_x( tet - valGenAngSdvig/2.)
							+ fncDiagrSinx_div_x( tet - valGenAngSdvig/2.*3.);
   parrBuff[ i * nBuffCols + 7] =  fncDiagrSinx_div_x( tet + valGenAngSdvig/2.)
							+ fncDiagrSinx_div_x( tet + valGenAngSdvig/2.*3.);



  }

 // double scalex = 100.;
  double *pscaley = new double [nBuffCols] ;
  pscaley[0] = 100.;
  pscaley[1] = 100.;
  pscaley[2] = 10.;
  pscaley[3] = 10.;
  pscaley[4] = 10.;
  pscaley[5] = 10.;
  pscaley[6] = 10.;
  pscaley[7] = 10.;
 // pscaley[6] = 1;
 // pscaley[7] = 1;
  for (int i=2; i < nBuffCols; i++)
  {

  TYrWriteShapeFile::WriteOneReport(                 wchFoldName  // путь к папке
								  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  , nBuffCols // - к-во переменных о корорых накоплена информация в буфере
								  , nBuffRows //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,lenName // максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i  // номер переменной по оси Y
								  ,1. //  масштаб по оси X
								  ,pscaley[i]// масштаб по оси Y
								   ) ;
  }

   for (int i=2; i < nBuffCols; i++)
  {

  TYrWriteShapeFile::WriteOneReport(                 wchFoldName  // путь к папке
								  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  , nBuffCols // - к-во переменных о корорых накоплена информация в буфере
								  , nBuffRows //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,lenName // максимальная длина имени переменной
								  ,1  // номер переменной по оси X
								  ,i  // номер переменной по оси Y
								  ,100. //  масштаб по оси X
								  ,pscaley[i]// масштаб по оси Y
								   ) ;
  }

	delete parrBuff;
	delete pscaley;
	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  wchFoldName);
	wcscat(wchAxesFileName, L"AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-150., 150.
	 ,-pscaley[0] * 1.1,pscaley[0] * 1.1, 50000.*step) ;



  delete wcharrFileNames;
}

// пересчет угла в рад в обобщенный угол
double TParAnt::transformAngToGeneralizedAng (const double  valAng  )
{
	return  mAppert * M_PI/ mLambda  * sin (valAng);
}

// пересчет обобщенного угла в обычный в рад
double TParAnt::transformGeneralizedAngToAng (const double  GeneralizedAng  )
{
	return GeneralizedAng  * mLambda/ mAppert/ M_PI;//asin ( GeneralizedAng  * mLambda/ mAppert/ M_PI);
}

// нахождение обощенного угла при котором амплитуда суммарной
// диаграммы 2 парциальных диаграмм равна sqrt(2)/2
double TParAnt::findTet07_For_SumDiagr()
{

 double tet0 = 1.8;//M_PI/2;

 int i =0;
 double valTemp = fncSumDiagr(   0.);
 double a = sqrt(2.)/2. * valTemp;
 for ( i =0; i < 10; i++)
 {
  double del = -(fncSumDiagr(tet0) - a)/ fncDerivSumDiagr(tet0);
  tet0 += del;
  if (fabs(del) < 0.0000001) break;

 }
 return tet0;
}

// нахождене первого нуля суммарной диаграммы  обобщ угол
double TParAnt::findFirstZero_For_SumDiagr()
{

	double valMu = 4.8; // первое пнриближение
	for (int i = 0; i < 100; i++)
	{
	  double valMu1 = valMu - fncSumDiagr( valMu)/fncDerivSumDiagr( valMu) ;
	  if (fabs(valMu1 - valMu)< 0.001)
	  {
		return  valMu1;
	  }
	  valMu = valMu1;
	}
	return -10000.0;
}

// нахождене второго нуля суммарной диаграммы  обобщ угол
double TParAnt::findSecondZero_For_SumDiagr()
{

	double valMu = 8.1; // первое пнриближение
	for (int i = 0; i < 100; i++)
	{
	  double valMu1 = valMu - fncSumDiagr( valMu)/fncDerivSumDiagr( valMu) ;
	  if (fabs(valMu1 - valMu)< 0.001)
	  {
		return  valMu1;
	  }
	  valMu = valMu1;
	}
	return -10000.0;
}


// функция суммарной диаграммы двух последовательных парциальных диаграмм
// INPUT:
// tetGeneralized  - аргумент, выраженный в обобщенном угле, для которого расвсчитывается
//            диаграмма. Угол исчисляется относительно центра суммарной диаграммы
//
double TParAnt::fncSumDiagr( double tetGeneralized)
{
 const double  valGenAngSdvig =  transformAngToGeneralizedAng ( mAngSdvig );
 return fncDiagrSinx_div_x( tetGeneralized - valGenAngSdvig/2.) + fncDiagrSinx_div_x( tetGeneralized + valGenAngSdvig/2.);

}

// функция производной суммарной диаграммы двух последовательных парциальных диаграмм
// INPUT:
// tetGeneralized  - аргумент, выраженный в обобщенном угле, для которого расвсчитывается
//            диаграмма. Угол исчисляется относительно центра суммарной диаграммы
double TParAnt::fncDerivSumDiagr( double tetGeneralized)
{
 const double  valGenAngSdvig =  transformAngToGeneralizedAng ( mAngSdvig );
 return fncDerivDiagrSinx_div_x( tetGeneralized - valGenAngSdvig/2.) + fncDerivDiagrSinx_div_x( tetGeneralized + valGenAngSdvig/2.);
}


// ширина суммарной диаграммы в обобщенном угле  по уровню 0,7
double TParAnt::findSumDiagrWidth07_GeneralizedAng()
{

  return 2. *findTet07_For_SumDiagr() ;

}

// ширина  суммарной  диаграммы в обычном угле по уровню 0,7, рад
double TParAnt::findSumDiagrWidth07_Rad()
{
	double valTet07 = findTet07_For_SumDiagr() ;
	return  2.* asin(valTet07 * mLambda/ mAppert/ M_PI);
}


// нахождение уровня пересечения двух суммарных сдвинутых диаграмм
double TParAnt::findCrossLevel_For_SumDiagr()
{
  const double  valGenAngSdvig =  transformAngToGeneralizedAng ( mAngSdvig );
  return  fncSumDiagr(valGenAngSdvig/2.)/ fncSumDiagr(0.);
}

// имитация замеров  первичных диаграмм
// INPUT:
// valalfUMTrg - угол цели в рад
// valalfUMAntp - угол антип в рад
// cmpKAntp
// OUTPUT:
// cmparrPartS []
// cmparrPartSZv []
void TParAnt::ImitateMeasureArrayPartialDiagrams( double valalfUMTrg, TComp cmpKTarg, double valalfUMAntp,TComp cmpKAntp
		, TComp *cmparrPartS, TComp *cmparrPartSZv)
{

   for (int i =0; i < mQuantDiagr; i++)
   {
	  // вычисление угла целей относителбно оси диаграммы с номером i в обобщенных координатах
	  double valDiagrAng = calcRadAngSdvigaPartDiagr(i) ; // угол оси диагр с номером i
	 //  ( -(((double)mQuantDiagr) -1.)/2. + i)*  mAngSdvig;
	  double valTrg =  valalfUMTrg - valDiagrAng; // угол цели в диаграмме с номером i, рад.
	  double valTrgGen = transformAngToGeneralizedAng (  valTrg ) ; // угол цели в диаграмме с номером i, обощ угол
	  double valAntp =  valalfUMAntp - valDiagrAng; // угол цели в диаграмме с номером i, рад.
	  double valAntpGen = transformAngToGeneralizedAng ( valAntp ) ; // угол цели в диаграмме с номером i, обощ угол
	  ///
	  cmparrPartS[i] = (cmpKTarg * TComp(fncDiagrSinx_div_x(valTrgGen), 0.)) + ( cmpKAntp * TComp(fncDiagrSinx_div_x(valAntpGen), 0.) );

	  // диаграмма с направления цели
	  double valDiagrTarg = fncDiagrSinx_div_x(valTrgGen);

	  // диаграмма с навления антипода
	  double valDiagrAntp = fncDiagrSinx_div_x(valAntpGen);

	  // дисперсия действительной части сигнала
	   double valDispRe = ( (valDiagrTarg * cmpKTarg.m_Re * valDiagrTarg * cmpKTarg.m_Re)
						  + (valDiagrAntp * cmpKAntp.m_Re * valDiagrAntp * cmpKAntp.m_Re) )
						   * mAmplFactSig * mAmplFactSig + mNoiseDisp;
	   // дисперсия мнипмой части сигнала
	   double valDispIm = ( (valDiagrTarg * cmpKTarg.m_Im * valDiagrTarg * cmpKTarg.m_Im)
						  + (valDiagrAntp * cmpKAntp.m_Im * valDiagrAntp * cmpKAntp.m_Im) )
						   * mAmplFactSig * mAmplFactSig + mNoiseDisp;


	  cmparrPartSZv[i].m_Re = cmparrPartS[i].m_Re  +  getGauss(0., sqrt(valDispRe));
	  cmparrPartSZv[i].m_Im = cmparrPartS[i].m_Im +  getGauss(0., sqrt(valDispIm));
   }
}

// имитация замеров  суммарных  диаграмм
void TParAnt::ImitateMeasureArraySumDiagrams(  TComp *cmparrPartS, TComp *cmparrPartSZv,  TComp *cmparrSumS, TComp *cmparrSumSZv)
{
  for (int i =0; i < (mQuantDiagr - 1); i++)
  {
	cmparrSumS[i] =  cmparrPartS[i] + cmparrPartS[i + 1];
	cmparrSumSZv[i] =  cmparrPartSZv[i] + cmparrPartSZv[i + 1];
  }
}

// оценивание обобщенных углов цели и антипода по измерениям 3 суммарных диаграмм
//INPUT:
// wchFoldName1 - путь к папаке с графиками
// cmparrS[3] -  массив замеров
//  OUTPUT:
//  *pvalGenTargEps -   угол цели, обобщ коорд
// *pvalGenAntpEps -    угол антипода, обобщ коорд
int TParAnt::EstGenAngsThreeSumDiagr(TComp *cmparrS
  ,  double  *pvalGenTargEps,  double  *pvalGenAntpEps, double *pval_b0, double *pval_b1)
{
  // вычисленние вектора b
   double arr_b[2] = {0.}  ;
   bool brez = calcVect_b(cmparrS,arr_b) ;
   ///

   *pval_b0 = arr_b[0];
   *pval_b1 = arr_b[1];
	double arrRoots[2] ={0.};
   int iNUmRoots = findRootsFgr_For_3SumDiagr(   arr_b,  arrRoots);
   *pvalGenTargEps =  arrRoots[1];
   *pvalGenAntpEps =  arrRoots[0];
   return iNUmRoots ;
}

// нахождение корней функции FGr длшя метода 3 суммарных диаграмм
// возвращает к-во корней
int TParAnt::findRootsFgr_For_3SumDiagr(double *arr_b,  double *arrRoots)
 {
   const double  valGenAngSdvig =  transformAngToGeneralizedAng (mAngSdvig );
   int iNUmRoots = 0;

   double valMu0 =  findFirstZero_For_SumDiagr() + valGenAngSdvig - 0.1;
   double valStep = 0.25;
	const int lenBuff = 2. *valMu0/ valStep ;
   double *parrBuff  = new double [lenBuff];

  double valmu = -valMu0;

  for (int i=0 ; i < lenBuff; i++)
  {
   valmu = -valMu0 + ((double)i) * valStep;
   parrBuff[ i ] = fncFGr_SumDiagr(arr_b, valmu);

  }

  for (int i = 0; i < (lenBuff-1); i++)
  {
   if (parrBuff[ i ] * parrBuff[ i + 1] <= 0.)
   {
	 if (iNUmRoots >1)
	 {
	   break;
	 }
	 double valX0 = -valMu0 + ((double)i) * valStep;
	 double valX1 = -valMu0 + ((double) (i +1)) * valStep;
	 arrRoots [iNUmRoots ] = findRootMethChord_For_Fgr_For_3SumDiagr(arr_b,  valX0, valX1) ;
	 iNUmRoots++;


   }
  }
	delete parrBuff;
  return iNUmRoots;

 }

 double TParAnt::fncFGr_SumDiagr(double *arr_b, double valmu)
{
  const double  valGenAngSdvig =  transformAngToGeneralizedAng (mAngSdvig );
  return (fncSumDiagr( valmu) - arr_b[0] * fncSumDiagr( valmu +valGenAngSdvig)
		 - arr_b[1] * fncSumDiagr( valmu - valGenAngSdvig));
}


// вычисление вектора b
// INPUT:
// cmparrS[3] - массив измерений тройки диаграмм
// OUTPUT:
// arr_b[2] -  коэффициенты уравнения
bool TParAnt::calcVect_b(TComp *cmparrS, double *arr_b)
{
 // вычисленние вектора b

/* double s11 = cmparrS[0].m_Re;
 double s12 = cmparrS[0].m_Im;
 double s21 = cmparrS[1].m_Re;
 double s22 = cmparrS[1].m_Im;
 double s31 = cmparrS[2].m_Re;
 double s32 = cmparrS[2].m_Im; */

 double temp = cmparrS[0].m_Re * cmparrS[2].m_Im - cmparrS[2].m_Re * cmparrS[0].m_Im;
 if (fabs(temp)< 0.00000000001)
 {
   return false;
 }

 arr_b [0] = (cmparrS[2].m_Im * cmparrS[1].m_Re - cmparrS[2].m_Re* cmparrS[1].m_Im)/ temp;
 arr_b [1] = (-cmparrS[0].m_Im * cmparrS[1].m_Re + cmparrS[0].m_Re * cmparrS[1].m_Im)/ temp;


   return true ;
}

double 	  TParAnt::findRootMethChord_For_Fgr_For_3SumDiagr(double *arr_b, double  valX0, double valX1)
{
  double valMu0 = valX0 - 0.00001;
  double valMu1 = valX1 + 0.00001;
  double valMuRez =  valMu0 ;
  double eps = 0.001;
  for (int i = 0; i < 100; i++)
  {
	double f0 = fncFGr_SumDiagr(arr_b,  valMu0 ) ;
	double f1 = fncFGr_SumDiagr(arr_b,   valMu1 ) ;
	double valMuRez0 = valMu0 - f0 * ( valMu1 - valMu0)/ (f1 - f0);
	if (fabs(valMuRez0 - valMuRez) < eps)
	{
	 return valMuRez0;
	}
	double f2 = fncFGr_SumDiagr(arr_b,   valMuRez0 ) ;
	if( f2* f0 < 0.)
	{
	   valMu1 = valMuRez0;

	}
	else
	{
      valMu0 = valMuRez0;
    }
     valMuRez = valMuRez0;
  }
	return -100000.;
}




// оценивание обобщенных углов цели и антипода по измерениям 3 суммарных диаграмм
//INPUT:
// cmparrS[6] -  массив замеров
// NumRayTriple-  номер тройки рабочих лучей (нумерация лучей идет снизу начиная с нуля)
//                номер тройки это порядковый номер нижнего луча
//  OUTPUT:
//  *valEstAngTarg -   угол цели, рад
// *valEstAngAntp -    угол антипода, рад
// cmpKTarg - коэф отражения цели
// cmpKAntp - коэф отражения антипода
// arrMtrxCorr  - коррел матрица ошибок измерения угла цели и антипода
int TParAnt::estimateMethThreeSumDiagr( TComp *cmparrS
, int iNumRayTriple , double *valEstAngTarg, double *valEstAngAntp, TComp *cmpKTarg , TComp *cmpKAntp , double *arrMtrxCorr
  , double *pval_b0, double *pval_b1 )
{
	if (iNumRayTriple > (mQuantDiagr -3))
	{
	  return -2;
	}
	const double  valGenAngSdvig =  transformAngToGeneralizedAng (mAngSdvig );
	// нахождение обобщенных углов цели и антипода
	double  valGenTargEps = -100.,  valGenAntpEps = -100.0;  // обобщ углы цели и антипода
	if (EstGenAngsThreeSumDiagr(  &cmparrS[iNumRayTriple], &valGenTargEps
	,  &valGenAntpEps, pval_b0, pval_b1) !=2)
	{
	   return -1;
	}
	//
	// вычисление коэффиц отражения
	TComp cmparrA[4],cmparrAInv[4], cmparrK[2], cmparrSTemp[2];

	cmparrA[0] = TComp(fncSumDiagr( valGenTargEps + valGenAngSdvig), 0.);
	cmparrA[1] = TComp(fncSumDiagr( valGenAntpEps + valGenAngSdvig), 0.);
	cmparrA[2] = TComp(fncSumDiagr( valGenTargEps - valGenAngSdvig), 0.);
	cmparrA[3] = TComp(fncSumDiagr( valGenAntpEps - valGenAngSdvig), 0.);

	bool brez =   InverseMtrx2( cmparrA, cmparrAInv);
	if (!brez)
	{
	 return -2;
	}

	cmparrSTemp[0] = cmparrS [iNumRayTriple ];
	cmparrSTemp[1] = cmparrS [iNumRayTriple + 2 ];
	MtrxMultMatrx(cmparrAInv,2, 2, cmparrSTemp,1, cmparrK);
	*cmpKTarg = cmparrK[0] ;
	*cmpKAntp = cmparrK[1] ;

	///

	// вычисление углоы цели и антипода в рад
		  // относиьтельный угол цетрального луча тройки суммарных диаграмм в радианах

			double valRadCentreRayPos = calcRadAngSdvigaSumDiagr(iNumRayTriple +1);
			//  ( -(((double)mQuantDiagr) - 1.)/2. + 1.5 + ((double)iNumRayTriple)) * mAngSdvig;

		  // вычисление угла цели в радианах в АСК
			double valRadTargEps =  transformGeneralizedAngToAng ( valGenTargEps ) ;// в радианы
			*valEstAngTarg = valRadTargEps + valRadCentreRayPos;
		  // вычисление угла антипода в радианах в АПСК
			double valRadAntpАEps =  transformGeneralizedAngToAng ( valGenAntpEps ) ;// в радианы
			*valEstAngAntp = valRadAntpАEps + valRadCentreRayPos;
	///

	// вычисление коррел матрицы ошибок определения обобщенных углов цели и антипода
	double arrMtrxCorrMu[4] = {0.};
	calcMtrxCorrGenAngs_Meth3SumDiagr(iNumRayTriple,valGenTargEps,  valGenAntpEps
	,*cmpKTarg,*cmpKAntp, arrMtrxCorrMu) ;
	///

	// вычисление корреляц матрицы ошибок определения углов цели и антипода в радианаах
	double arrQ[4] = {0.}, arrTemp[4] = {0.};
	double coef = mLambda / mAppert / M_PI;
	arrQ[0] =  coef / sqrt( 1.- coef * valGenTargEps * coef * valGenTargEps);
	arrQ[3] =  coef / sqrt( 1.- coef * valGenAntpEps * coef * valGenAntpEps);
	MtrxMultMatrx(arrQ,2, 2, arrMtrxCorrMu,2, arrTemp)  ;
	MtrxMultMatrxTransp(arrTemp,2, 2, arrQ,2, arrMtrxCorr) ;
	return 1;

}

bool TParAnt::calcMtrxCorrGenAngs_Meth3SumDiagr(const int NumAnsamble, double valGenTargEps,  double valGenAntpEps
  ,TComp cmpKTarg ,TComp cmpKAntp , double *arrMtrxCorrGenAngs )
{
  const double  valGenAngSdvig =  transformAngToGeneralizedAng (mAngSdvig );
  // формирование матирицы частных производных dG/dx
  double arrdG_po_dX[36] = {0.};
  arrdG_po_dX[0] =  cmpKTarg.m_Re * fncDerivSumDiagr(  valGenAngSdvig + valGenTargEps);
  arrdG_po_dX[1] =  cmpKAntp.m_Re * fncDerivSumDiagr(  valGenAngSdvig + valGenAntpEps);
  arrdG_po_dX[2] =  fncSumDiagr(  valGenAngSdvig + valGenTargEps);
  arrdG_po_dX[4] =  fncSumDiagr(  valGenAngSdvig + valGenAntpEps);

  arrdG_po_dX[6] =  cmpKTarg.m_Im * fncDerivSumDiagr( valGenAngSdvig + valGenTargEps);
  arrdG_po_dX[7] =  cmpKAntp.m_Im * fncDerivSumDiagr( valGenAngSdvig + valGenAntpEps);
  arrdG_po_dX[9] =  fncSumDiagr( valGenAngSdvig + valGenTargEps);
  arrdG_po_dX[11] =  fncSumDiagr( valGenAngSdvig + valGenAntpEps);

  arrdG_po_dX[12] =  cmpKTarg.m_Re * fncDerivSumDiagr(  valGenTargEps);
  arrdG_po_dX[13] =  cmpKAntp.m_Re * fncDerivSumDiagr(  valGenAntpEps);
  arrdG_po_dX[14] =  fncSumDiagr(  valGenTargEps);
  arrdG_po_dX[16] =  fncSumDiagr(  valGenAntpEps);

  arrdG_po_dX[18] =  cmpKTarg.m_Im * fncDerivSumDiagr(  valGenTargEps);
  arrdG_po_dX[19] =  cmpKAntp.m_Im * fncDerivSumDiagr(  valGenAntpEps);
  arrdG_po_dX[21] =  fncSumDiagr(  valGenTargEps);
  arrdG_po_dX[23] =  fncSumDiagr(  valGenAntpEps);

  arrdG_po_dX[24] =  cmpKTarg.m_Re * fncDerivSumDiagr( -valGenAngSdvig + valGenTargEps);
  arrdG_po_dX[25] =  cmpKAntp.m_Re * fncDerivSumDiagr( -valGenAngSdvig + valGenAntpEps);
  arrdG_po_dX[26] =  fncSumDiagr( -valGenAngSdvig + valGenTargEps);
  arrdG_po_dX[28] =  fncSumDiagr( -valGenAngSdvig + valGenAntpEps);

  arrdG_po_dX[30] =  cmpKTarg.m_Im * fncDerivSumDiagr( -valGenAngSdvig + valGenTargEps);
  arrdG_po_dX[31] =  cmpKAntp.m_Im * fncDerivSumDiagr( -valGenAngSdvig + valGenAntpEps);
  arrdG_po_dX[33] =  fncSumDiagr( -valGenAngSdvig + valGenTargEps);
  arrdG_po_dX[35] =  fncSumDiagr( -valGenAngSdvig + valGenAntpEps);
  ///

  // формирование корреляц матрицы измерений
 /* double arrCorrPartDiagr[64] = {0.} // коррел матрица измерений 4 парциальных диаграмм
  , arrLinTransf[48] ={0.}, arrMtrxCorrSumDiagr[36] ={0.}, arrTEmp[48] = {0.};
  int iNumAnsamble = 1;
  int lenAnsamble = 4;
  createMtrxCorrForMeasuresForAnsambleOfPartDiagrs( valGenTargEps + valGenAngSdvig/2.,   valGenAntpEps + valGenAngSdvig/2.
   , cmpKTarg , cmpKAntp ,  iNumAnsamble,  lenAnsamble  , arrCorrPartDiagr);


  for (int i =0; i < 6; i++)
  {
	arrLinTransf [ 8 * i + i] = 1.;
	arrLinTransf [ 8 * i + i + 2] = 1.;
  }

  MtrxMultMatrx(arrLinTransf,6, 8, arrCorrPartDiagr,8, arrTEmp) ;
  MtrxMultMatrxTransp(arrTEmp,6, 8, arrLinTransf,6, arrMtrxCorrSumDiagr) ;  */



   double arrMtrxCorrSumDiagr[36] ={0.};
   int lenAnsamble1 = 3;
   createMtrxCorrForMeasuresForAnsambleOfSumDiagrs( valGenTargEps,   valGenAntpEps
  , cmpKTarg , cmpKAntp,  NumAnsamble, lenAnsamble1
  , arrMtrxCorrSumDiagr ) ;
  ///

  // вычисление матрицы dX/dS
  double arrMtrx_dX_po_dS[36] = {0.};
  bool brez = InverseMtrx6(arrdG_po_dX,arrMtrx_dX_po_dS);
  if (!brez)
  {
	return false;
  }
  ///

  // вычисление коррел матрицы вектора дельтах
  double arrMtrxCorrX[36] ={0.}, arrtemp0[36] = {0.};
  MtrxMultMatrx(arrMtrx_dX_po_dS,6, 6, arrMtrxCorrSumDiagr, 6, arrtemp0) ;
  MtrxMultMatrxTransp(arrtemp0,6, 6, arrMtrx_dX_po_dS,6, arrMtrxCorrX) ;
  ///

  // формирование корреляц матрицы ошибок обобщенных углов
  arrMtrxCorrGenAngs[0]  = arrMtrxCorrX[0] ;
  arrMtrxCorrGenAngs[1]  = arrMtrxCorrX[1] ;
  arrMtrxCorrGenAngs[2]  = arrMtrxCorrX[6] ;
  arrMtrxCorrGenAngs[3]  = arrMtrxCorrX[7] ;
  return true;

}

//формирование корреляционной матрицы ошибок первичных измерений
// ансамбля первичных диаграмм
// вектор первичных ошибок имеет вид ...S[i].m_Re S , S[i].m_Im...
// INPUT:
// valGenTargEps  - обобщ угол цели
// valGenAntpEps - обобщ угол антипода
// valGenAngSdvig - обобщ угол сдвига диаграмм
// cmpKTarg, cmpKAntp - сигналы(комлексные) цели и антипода
//  iNumAnsamble   - намер первой диаграммы в ансамбле
//          (все диаграммы перенумерованы в порядке возрастания угла сдвига начиная с нулевого номера)
// lenAnsamble  - количество диаграмм в ансамбле
// OTPUT:
// arrMtrxCorr[2 *lenAnsamble *2 *lenAnsamble] - корреляц матрица
void TParAnt::createMtrxCorrForMeasuresForAnsambleOfPartDiagrs(double valGenTargEps,  double valGenAntpEps
  ,TComp cmpKTarg ,TComp cmpKAntp, int iNumAnsamble, int lenAnsamble
  , double *arrMtrxCorr )
 {
	const double  valGenAngSdvig =  transformAngToGeneralizedAng (mAngSdvig );
	memset(arrMtrxCorr, 0, 4 * lenAnsamble * lenAnsamble * sizeof(double));

  // вычисление коррел матрицы, обусловленной разбросом коэффиц усиления диаграмм
	for (int i = 0; i < lenAnsamble; i++)  // цикл по диаграммам в ансамбле
	{
	   // вычисление угла целей относителбно оси диаграммы с номером i в обобщенных координатах
	  double valDiagrGenAng = calGenAngSdvigaPartDiagr(iNumAnsamble + i); // угол оси диагр с номером i
	 // ((double)(-mQuantDiagr/2 + iNumAnsamble + i))*  valGenAngSdvig;
	  double valTargGen =  valGenTargEps - valDiagrGenAng; // угол цели в диаграмме с номером i, рад.
	  double valAntpGen =  valGenAntpEps - valDiagrGenAng; // угол цели в диаграмме с номером i, рад.
	  ///
	  // диаграмма с направления цели
	  double valDiagrTarg = fncDiagrSinx_div_x(valTargGen);

	  // диаграмма с навления антипода
	  double valDiagrAntp = fncDiagrSinx_div_x(valAntpGen);

	  // дисперсия действительной части сигнала
	   double valDispRe = ( (valDiagrTarg * cmpKTarg.m_Re * valDiagrTarg * cmpKTarg.m_Re)
						  + (valDiagrAntp * cmpKAntp.m_Re * valDiagrAntp * cmpKAntp.m_Re) )
						   * mAmplFactSig * mAmplFactSig;
	   // дисперсия мнипмой части сигнала
	   double valDispIm = ( (valDiagrTarg * cmpKTarg.m_Im * valDiagrTarg * cmpKTarg.m_Im)
						  + (valDiagrAntp * cmpKAntp.m_Im * valDiagrAntp * cmpKAntp.m_Im) )
						   * mAmplFactSig * mAmplFactSig;
	 ///
	  // заполнение диагональных элементов матрицы arrMtrxCorr в соответствующих строках
	  // действительная и комплек5сная части сигнала диаграммы с номером i занимают строки 2*i и 2*i+1 соответственно
	  arrMtrxCorr[2*lenAnsamble * (2 * i) + 2 *i] =  valDispRe  + mNoiseDisp;
	  arrMtrxCorr[2*lenAnsamble * (2 * i + 1) + 2 *i + 1] =  valDispIm  + mNoiseDisp;

	}
	return;
 }

//формирование корреляционной матрицы ошибок первичных измерений
// ансамбля уммарных диаграмм
// вектор первичных ошибок имеет вид ...S[i].m_Re S , S[i].m_Im...
// INPUT:
// valGenTargEps  - обобщ угол цели
// valGenAntpEps - обобщ угол антипода

// cmpKTarg, cmpKAntp - сигналы(комлексные) цели и антипода
//  iNumAnsamble   - намер первой диаграммы в ансамбле
//          (все диаграммы перенумерованы в порядке возрастания угла сдвига начиная с нулевого номера)
// lenAnsamble  - количество диаграмм в ансамбле
// OTPUT:
// arrMtrxCorr[2 *lenAnsamble *2 *lenAnsamble] - корреляц матрица
void TParAnt::createMtrxCorrForMeasuresForAnsambleOfSumDiagrs( double valGenTargEps,  double valGenAntpEps
  ,TComp cmpKTarg ,TComp cmpKAntp, int iNumAnsamble, int lenAnsamble
  , double *arrMtrxCorrSumDiagr )
{


  double *arrCorrPartDiagr = new double [4 * (lenAnsamble + 1 )*(lenAnsamble + 1 )]  ;
  double *arrLinTransf = new double [2 * (lenAnsamble + 1) * 2 * lenAnsamble];
  double *arrTEmp = new double[2 * (lenAnsamble + 1) * 2 * lenAnsamble];
  memset(arrLinTransf, 0, 2 * (lenAnsamble + 1) * 2 * lenAnsamble * sizeof(double));
  createMtrxCorrForMeasuresForAnsambleOfPartDiagrs( valGenTargEps ,   valGenAntpEps
   , cmpKTarg , cmpKAntp ,  iNumAnsamble,  lenAnsamble + 1  , arrCorrPartDiagr);


  for (int i =0; i < 2 * lenAnsamble; i++)
  {
	arrLinTransf [ 2 * (lenAnsamble + 1) * i + i] = 1.;
	arrLinTransf [ 2 * (lenAnsamble + 1) * i + i + 2] = 1.;
  }

  MtrxMultMatrx(arrLinTransf,2 * lenAnsamble, 2 * (lenAnsamble + 1), arrCorrPartDiagr,2 * (lenAnsamble + 1), arrTEmp) ;
  MtrxMultMatrxTransp(arrTEmp,2 * lenAnsamble,2 * (lenAnsamble + 1), arrLinTransf,2 * lenAnsamble, arrMtrxCorrSumDiagr) ;
  delete arrCorrPartDiagr;
  delete arrLinTransf;
  delete arrTEmp;
}
// INPUT
//  valScaleY  - растяжение пол оси Y
//
//
//
//
//


void TParAnt::createPtesentGraphs_for_SumDiagr( wchar_t *wchrPresntSumDiagrams, int iNumRayTriple
		,double  alfUMTrg,TComp cmpKTarg,double  alfUMAntp, TComp cmpKAntp
		,double  valEstRadAngTarg,double  valEstRadAngAntp
		  ,TComp cmpKTargEst , TComp cmpKAntpEst, double val_b0, double val_b1, double valSMeshenieVert
		  , double valScaleY )
{
   // построение 6 суммарных ьдиаграмм в зависимости от угла в рад
   wchar_t Fold[400] ={0.};
   wcscpy( Fold, wchrPresntSumDiagrams);
   wcscat(Fold, L"\\6_SumDiagr");
	_wmkdir(Fold);

   for (int i =0; i < (mQuantDiagr -1); i++)
   {
	wchar_t FileName[400] ={0.};
	wcscpy( FileName, Fold);
	wcscat(FileName, L"\\SumDiagrNo_");
	wchar_t string[10] = {0.};
	_itow(i,string, 10);
	wcscat(FileName, string);
	wcscat(FileName, L".shp");

	// угол смещения суммарной диаграммы с номером i относительно нормали к антенне в рад
	double valSmeshenieGoriz =  calcRadAngSdvigaSumDiagr(i);
	// ( - ((double)mQuantDiagr -1.) /2. + 1.5 + ((double) i)) * mAngSdvig;
	double valSMeshenieVert = 0.;
	createGraphSumDiagr_from_rad(FileName, valSmeshenieGoriz, valSMeshenieVert, valScaleY);
   }
   ///

   // построение 3 рабочих суммарных диаграмм в зависимости от угла в рад
   wcscpy( Fold, wchrPresntSumDiagrams);
   wcscat(Fold, L"\\3_WorkingDiagr");
	_wmkdir(Fold);
	for (int i =0; i < 3; i++)
   {
	wchar_t FileName[400] ={0.};
	wcscpy( FileName, Fold);
	wcscat(FileName, L"\\WorkingSumDiagrNo_");
	wchar_t string[10] = {0.};
	_itow(i,string, 10);
	wcscat(FileName, string);
	wcscat(FileName, L".shp");
	double valSmeshenieGoriz = calcRadAngSdvigaSumDiagr(iNumRayTriple + i);
	//( - ((double)mQuantDiagr -1.) /2. + 0.5 + ((double)(iNumRayTriple + i))) * mAngSdvig;

	createGraphSumDiagr_from_rad(FileName, valSmeshenieGoriz, 0., valScaleY ) ;
   }
   ///

   // построение графика функции FGreece в зависимости от угла рад
   wchar_t Fold1[400] ={0.};
   wcscpy( Fold1, wchrPresntSumDiagrams);
   wcscat(Fold1, L"\\GraphFGreece_from_rad_For_3SumDiagr");
	_wmkdir(Fold1);
	wcscat(Fold1, L"\\");

   // относиьтельный угол цетрального луча тройки в радианах
   double valRadCentreRayPos =  calcRadAngSdvigaSumDiagr(iNumRayTriple + 1);
   // ( -((double)mQuantDiagr -1.)/2. + 1.5 + ((double)iNumRayTriple)) * mAngSdvig;

   createGraphFGreece_from_rad_For_3SumDiagr(Fold1, valRadCentreRayPos , valSMeshenieVert, val_b0, val_b1)  ;
   ///

   // построение графиков вертикальных линий антипода и цели
	 wchar_t File[400] ={0.};
	// истинная цель
	TURPointXY  pnt1(alfUMTrg * 100., 0.);
	TURPointXY  pnt2(alfUMTrg * 100., cmpKTarg.modul()/ 800. * 5.);
	TURPolyLine plnTargTrue(   pnt1,   pnt2) ;
	wcscpy( File, Fold);
	wcscat(File, L"\\plnTargTrue.shp");
	plnTargTrue.WriteSetSHPFiles(File ,&plnTargTrue, 1) ;

   ///
   // истинный антипод
	TURPointXY  pnt3(alfUMAntp * 100., 0.);
	TURPointXY  pnt4(alfUMAntp * 100., cmpKAntp.modul()/ 800. * 5.);
	TURPolyLine plnAntpTrue(   pnt3,   pnt4) ;
	wcscpy( File, Fold);
	wcscat(File, L"\\plnAntpTrue.shp");
	plnAntpTrue.WriteSetSHPFiles(File ,&plnAntpTrue, 1) ;
	///

	// оценка цели
	TURPointXY  pnt5(valEstRadAngTarg * 100., 0.);
	TURPointXY  pnt6(valEstRadAngTarg * 100., cmpKTargEst.modul()/ 800. * 5.);
	TURPolyLine plnTargEst(   pnt5,   pnt6) ;
	wcscpy( File, Fold);
	wcscat(File, L"\\plnTargEst.shp");
	plnTargEst.WriteSetSHPFiles(File ,&plnTargEst, 1) ;
	///

	// оценка антипода
   TURPointXY  pnt7(valEstRadAngAntp * 100., 0.);
	TURPointXY  pnt8(valEstRadAngAntp * 100., cmpKAntpEst.modul()/ 800. * 5.);
	TURPolyLine plnAntpEst(   pnt7,   pnt8) ;
	wcscpy( File, Fold);
	wcscat(File, L"\\plnAntpEst.shp");
	plnAntpEst.WriteSetSHPFiles(File ,&plnAntpEst, 1) ;
	///

	// вертик линия оценки цели
	TURPointXY  pnt9(valEstRadAngTarg * 100., cmpKTargEst.modul()/ 800. * 5.);
	TURPointXY  pnt10(valEstRadAngTarg * 100., valSMeshenieVert);
	TURPolyLine plnTargEst1(   pnt9,   pnt10) ;
	wcscpy( File, Fold);
	wcscat(File, L"\\plnTargEst1.shp");
	plnTargEst1.WriteSetSHPFiles(File ,&plnTargEst1, 1) ;
	///

	// вертик линия оценки антипода
	TURPointXY  pnt11(valEstRadAngAntp * 100., cmpKAntpEst.modul()/ 800. * 5.);
	TURPointXY  pnt12(valEstRadAngAntp * 100., valSMeshenieVert);
	TURPolyLine plnAntpEst1(   pnt11,   pnt12) ;
	wcscpy( File, Fold);
	wcscat(File, L"\\plnAntpEst1.shp");
	plnAntpEst1.WriteSetSHPFiles(File ,&plnAntpEst1, 1) ;
	///
   //-------------------------------------
	// изображение углов в виде полигонов
	// полигон оценки цели
	TURPointXY  pntTopRight0(valEstRadAngTarg * 100. + 0.025, cmpKTargEst.modul()/ 800. * 5.);
	TURPointXY  pntBottomLeft0(valEstRadAngTarg * 100. - 0.025, 0.);
	TURPolygon plgTargEst(pntTopRight0, pntBottomLeft0);
	wcscpy( File, Fold);
	wcscat(File, L"\\plgTargEst.shp");
	plgTargEst.WriteSetSHPFiles(File ,&plgTargEst, 1) ;
	///

	// полигон оценки антипода
	TURPointXY  pntTopRight1(valEstRadAngAntp * 100. + 0.025, cmpKAntpEst.modul()/ 800. * 5.);
	TURPointXY  pntBottomLeft1(valEstRadAngAntp * 100. - 0.025, 0.);
	TURPolygon plgAntpEst(pntTopRight1, pntBottomLeft1);
	wcscpy( File, Fold);
	wcscat(File, L"\\plgAntpEst.shp");
	plgAntpEst.WriteSetSHPFiles(File ,&plgAntpEst, 1) ;
	///

	// полигон истинная цель
	TURPointXY  pntTopRight2(alfUMTrg * 100.+ 0.025, cmpKTarg.modul()/ 800. * 5.);
	TURPointXY  pntBottomLeft2(alfUMTrg * 100.- 0.025, 0.);
	TURPolygon plgTargTrue(pntTopRight2, pntBottomLeft2);
	wcscpy( File, Fold);
	wcscat(File, L"\\plgTargTrue.shp");
	plgTargTrue.WriteSetSHPFiles(File ,&plgTargTrue, 1) ;

   ///
   //  полигон  истинный антипод
	TURPointXY  pntTopRight3(alfUMAntp * 100.+ 0.025, cmpKAntp.modul()/ 800. * 5.);
	TURPointXY  pntBottomLeft3(alfUMAntp * 100.- 0.025, 0.);
	TURPolygon plgAntpTrue(pntTopRight3, pntBottomLeft3);
	wcscpy( File, Fold);
	wcscat(File, L"\\plgAntpTrue.shp");
	plgAntpTrue.WriteSetSHPFiles(File ,&plgAntpTrue, 1) ;

	///
	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  Fold);
	wcscat(wchAxesFileName, L"\\AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-13, 13.5
	 ,-10.,10., 0.25) ;

}

// valAngSdvig - сдвиг диаграмм в рад
void TParAnt::createGraphSumDiagr_from_rad(wchar_t *FileName,double  valSmeshenieGoriz, double  valSMeshenieVert,  double valScaleY )
{

	double step =0.0001;
	double  valGenAngSdvig = transformAngToGeneralizedAng (mAngSdvig  ) ;
	double valZero =  findFirstZero_For_SumDiagr();
	const double VAlDiap = 1.5 *valZero;  // диапахзон графиков по обобщ координате плюс минус
	const int nBuffRows = 2. *VAlDiap/ step ;
	 TURPolyLine pln( 1, nBuffRows) ;



  double tet = -VAlDiap;
  for (int i=0 ; i < nBuffRows; i++)
  {
   tet = -VAlDiap + ((double)i) * step;
   pln.Points[i].X = (transformGeneralizedAngToAng ( tet  ) + valSmeshenieGoriz)* 100. ;
   pln.Points[i].Y =fncSumDiagr( tet) * valScaleY + valSMeshenieVert;

  }

 pln.WriteSetSHPFiles(FileName, &pln, 1);
}

void TParAnt::createGraphFGreece_from_rad_For_3SumDiagr(wchar_t *Fold,double  valSmeshenieGoriz
	, double  valSMeshenieVert, double val_b0, double val_b1 )
{
	double arr_b[2] = {0.};
	arr_b[0] =  val_b0 ;
	arr_b[1]  = val_b1;
   // построение графика функции FGreece
   double valGenAngSdvig = transformAngToGeneralizedAng (mAngSdvig ) ;
   double valMu0 =  findFirstZero_For_SumDiagr() + valGenAngSdvig - 0.1;
   double valStep = 0.00001;
	const int nBuffRows = 2. *valMu0/ valStep ;
	const int nBuffCols = 2;
	double  *parrBuff = new double [nBuffRows  * nBuffCols] ;
	memset(parrBuff, 0, nBuffRows * nBuffCols * sizeof(double));

	const int lenName = 30;
	wchar_t *wcharrFileNames  = new wchar_t [30 * nBuffCols];
	memset(wcharrFileNames, 0, 30 * nBuffCols * sizeof(wchar_t));

	wcscpy( wcharrFileNames, L"Mu");
	wcscpy( &wcharrFileNames[30], L"FGr_SumDiagr");


  double valmu = -valMu0;
  for (int i=0 ; i < nBuffRows; i++)
  {
   valmu = -valMu0 + ((double)i) * valStep;
   parrBuff[ i * nBuffCols] = (transformGeneralizedAngToAng ( valmu ) + valSmeshenieGoriz) * 100.;
   parrBuff[ i * nBuffCols + 1] = fncFGr_SumDiagr(arr_b,  valmu) + valSMeshenieVert;

  }

 // double scalex = 100.;
  double *pscaley = new double [nBuffCols] ;
  pscaley[0] = 100.;
  pscaley[1] = 1.;
 // wchar_t wchFoldName[] = L"E:\\PROJECTS_C++\\НЛЦ-5П10-03\\Новая папка\\GraphFGr\\";
  for (int i=1 ; i < nBuffCols; i++)
  {

  TYrWriteShapeFile::WriteOneReport(Fold  // путь к папке
								  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  , nBuffCols // - к-во переменных о корорых накоплена информация в буфере
								  , nBuffRows //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,lenName // максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i  // номер переменной по оси Y
								  ,1. //  масштаб по оси X
								  ,pscaley[i]// масштаб по оси Y
								   ) ;
  }



	TURPointXY pntSdvig(valSmeshenieGoriz, valSMeshenieVert);
	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  Fold);
	wcscat(wchAxesFileName, L"\\AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-13., 13.5
	 ,-13.,13.5, -parrBuff[ 0] * 0.05, pntSdvig) ;
	 delete  parrBuff;
	 delete wcharrFileNames;
	 delete pscaley;
}


// построение графика функции F греческое для метода 3 суммарных диаграмм в зависимости от обобщенного угла
void TParAnt::createGraphFGreece_For_3SumDiagr_from_GenAng(wchar_t *Fold,double  valSmeshenieGoriz
	, double  valSMeshenieVert, double val_b0, double val_b1 )
{
	double arr_b[2] = {0.};
	arr_b[0] =  val_b0 ;
	arr_b[1]  = val_b1;
   // построение графика функции FGreece
   double valGenAngSdvig = transformAngToGeneralizedAng (mAngSdvig ) ;
   double valMu0 =  findFirstZero_For_SumDiagr() + valGenAngSdvig - 0.1;
   double valStep = 0.00001;
	const int nBuffRows = 2. *valMu0/ valStep ;
	const int nBuffCols = 2;
	double  *parrBuff = new double [nBuffRows  * nBuffCols] ;
	memset(parrBuff, 0, nBuffRows * nBuffCols * sizeof(double));

	const int lenName = 30;
	wchar_t *wcharrFileNames  = new wchar_t [30 * nBuffCols];
	memset(wcharrFileNames, 0, 30 * nBuffCols * sizeof(wchar_t));

	wcscpy( wcharrFileNames, L"Mu");
	wcscpy( &wcharrFileNames[30], L"FGr_SumDiagr");


  double valmu = -valMu0;
  for (int i=0 ; i < nBuffRows; i++)
  {
   valmu = -valMu0 + ((double)i) * valStep;
   parrBuff[ i * nBuffCols] =  valmu  + valSmeshenieGoriz;
   parrBuff[ i * nBuffCols + 1] = fncFGr_SumDiagr(arr_b,  valmu) + valSMeshenieVert;

  }

 // double scalex = 100.;
  double *pscaley = new double [nBuffCols] ;
  pscaley[0] = 1.;
  pscaley[1] = 1.;
 // wchar_t wchFoldName[] = L"E:\\PROJECTS_C++\\НЛЦ-5П10-03\\Новая папка\\GraphFGr\\";
  for (int i=1 ; i < nBuffCols; i++)
  {

  TYrWriteShapeFile::WriteOneReport(Fold  // путь к папке
								  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  , nBuffCols // - к-во переменных о корорых накоплена информация в буфере
								  , nBuffRows //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,lenName // максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i  // номер переменной по оси Y
								  ,1. //  масштаб по оси X
								  ,pscaley[i]// масштаб по оси Y
								   ) ;
  }



	TURPointXY pntSdvig(valSmeshenieGoriz, valSMeshenieVert);
	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  Fold);
	wcscat(wchAxesFileName, L"\\AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-13., 13.5
	 ,-13.,13.5, -parrBuff[ 0] * 0.05, pntSdvig) ;
	 delete  parrBuff;
	 delete wcharrFileNames;
	 delete pscaley;
}


// построение графико СКЗ оши бок оценивания в звисимости от углового расстояния
// INPUT:
// Fold
// mLambda - раб частота
// mAppert - апертура антиенны
// valAngSdvig - угол сдвига между диаграммами в рад
// valSmeshenieGoriz  - смещение графика по горизонтали
// valSMeshenieVert   - смещение графика по вертикали
// cmpKTarg  - комплексный	 сигнал цели
// cmpKAntp  - комплексный	 сигнал антипода
// valNoiseSkz - скз шума в парциальной диаграмме
// mAmplFactSig - скз шума коэффиц усиления диаграммы
// INPUT:
// графики зависимости ско ошибок оценивания цели и антипода в зависимости
// от углового расстояния между целью и антиподом
// цельнаходится в центре центрального луча
// антипод приближается к цели
void TParAnt::createGraphsSKZ_from_AngDiffer_For_3SumDiagr(wchar_t *Fold,double  valSmeshenieGoriz
	, double  valSMeshenieVert, TComp cmpKTarg, TComp cmpKAntp )
{

   double valGenAngSdvig = transformAngToGeneralizedAng (mAngSdvig ) ;

   double valStep = 0.001;
	const int nBuffRows = 3. * findFirstZero_For_SumDiagr()/ valStep ;
	const int nBuffCols = 4;
	double  *parrBuff = new double [nBuffRows  * nBuffCols] ;
	memset(parrBuff, 0, nBuffRows * nBuffCols * sizeof(double));

	const int lenName = 30;
	wchar_t *wcharrFileNames  = new wchar_t [30 * nBuffCols];
	memset(wcharrFileNames, 0, 30 * nBuffCols * sizeof(wchar_t));

	wcscpy( wcharrFileNames, L"DeltaAng");
	wcscpy( &wcharrFileNames[30], L"SigmaTarg");
	wcscpy( &wcharrFileNames[60], L"SigmaAntp");
	wcscpy( &wcharrFileNames[90], L"CoefCorrel");

  double valRadAngTarg = mAngSdvig /2.;
  double valGenAngTarg =  transformAngToGeneralizedAng (valRadAngTarg ) ;
  double valGenAngAntp = 0.;
  double   arrMtrxCorr[4] ={0.}, arrMtrxCorrMu[4] = {0.};

  for (int i=0 ; i < nBuffRows; i++)
  {
   valGenAngAntp = valGenAngTarg - ((double)i + 10.) * valStep;


   	// вычисление коррел матрицы ошибок определения обобщенных углов цели и антипода
	calcMtrxCorrGenAngs_Meth3SumDiagr(2, valGenAngTarg ,  valGenAngAntp
	  , cmpKTarg, cmpKAntp, arrMtrxCorrMu) ;
	///

	// вычисление корреляц матрицы ошибок определения углов цели и антипода в радианаах
	double arrQ[4] = {0.}, arrTemp[4] = {0.};
	double coef = mLambda / mAppert / M_PI;
	arrQ[0] =  coef / sqrt( 1.- coef * valGenAngTarg * coef * valGenAngTarg);
	arrQ[3] =  coef / sqrt( 1.- coef *  valGenAngAntp * coef *  valGenAngAntp);
	MtrxMultMatrx(arrQ,2, 2, arrMtrxCorrMu,2, arrTemp)  ;
	MtrxMultMatrxTransp(arrTemp,2, 2, arrQ,2, arrMtrxCorr) ;


   parrBuff[ i * nBuffCols] = valRadAngTarg - transformGeneralizedAngToAng (valGenAngAntp);

   parrBuff[ i * nBuffCols + 1] =  sqrt(arrMtrxCorr[0])* 1000. ;
   parrBuff[ i * nBuffCols + 2] = sqrt(arrMtrxCorr[3])* 1000. ;
   if ((arrMtrxCorr[0] >0.) && (arrMtrxCorr[3] >0.))
   {
   parrBuff[ i * nBuffCols + 3] = arrMtrxCorr[1] / sqrt(arrMtrxCorr[0])/sqrt(arrMtrxCorr[3]) ;

   }
   else
   {
	   parrBuff[ i * nBuffCols + 3] = 0.;
   }

  }


  double *pscaley = new double [nBuffCols] ;
  pscaley[0] = 1000.;
  pscaley[1] = 1.;
  pscaley[2] = 1.;
  pscaley[3] = 1.;

  for (int i=1 ; i < nBuffCols; i++)
  {

  TYrWriteShapeFile::WriteOneReport(Fold  // путь к папке
								  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  , nBuffCols // - к-во переменных о корорых накоплена информация в буфере
								  , nBuffRows //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,lenName // максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i  // номер переменной по оси Y
								  ,100. //  масштаб по оси X
								  ,pscaley[i]// масштаб по оси Y
								   ) ;
  }



	TURPointXY pntSdvig(valSmeshenieGoriz, valSMeshenieVert);
	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  Fold);
	wcscat(wchAxesFileName, L"AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-13., 13.5
	 ,-13.,13.5, -parrBuff[ 0] * 0.05, pntSdvig) ;

	 TURPointXY pntSdvig1(valSmeshenieGoriz, 0.);
	wchar_t wchAxesFileName1[300] ={0};
	wcscpy(  wchAxesFileName1,  Fold);
	wcscat(wchAxesFileName1, L"AxesArr1.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName1,-13., 13.5
	 ,-13.,13.5, -parrBuff[ 0] * 0.05, pntSdvig1) ;

	 // построение графиков 3-х суммарных диаграмм снизу

	wcscpy(  wchAxesFileName,  Fold);
	wcscat(wchAxesFileName, L"DigarSum1.shp");
	 createGraphSumDiagr_from_rad(wchAxesFileName, valSmeshenieGoriz,   valSMeshenieVert,  4. );

		wcscpy(  wchAxesFileName,  Fold);
	wcscat(wchAxesFileName, L"DigarSum0.shp");
	 createGraphSumDiagr_from_rad(wchAxesFileName, valSmeshenieGoriz -mAngSdvig,   valSMeshenieVert,  4. );

		wcscpy(  wchAxesFileName,  Fold);
	wcscat(wchAxesFileName, L"DigarSum2.shp");
	 createGraphSumDiagr_from_rad(wchAxesFileName,valSmeshenieGoriz+ mAngSdvig,   valSMeshenieVert,  4. );

	// полигон изображающий истинная цель
	TURPointXY  pntTopRight2(0.05, cmpKTarg.modul()/ 800. * 5.);
	TURPointXY  pntBottomLeft2(- 0.05, 0.);
	TURPolygon plgTargTrue(pntTopRight2, pntBottomLeft2);
	wcscpy( wchAxesFileName, Fold);
	wcscat(wchAxesFileName, L"\\plgTargTrue.shp");
	plgTargTrue.WriteSetSHPFiles(wchAxesFileName ,&plgTargTrue, 1) ;
	 delete  parrBuff;
	 delete wcharrFileNames;
	 delete pscaley;

}

// графики для сравнения суммарной диаграммы м ее аппроксимации полинрмом 4 степени
void TParAnt::createGraphs_for_Compare_Diagrams(wchar_t *Fold  )
{
	double valGenAngSdvig = transformAngToGeneralizedAng (mAngSdvig ) ;
	double valStep =0.01;

	double valZero =  findFirstZero_For_SumDiagr();
	const double VAlDiap = valZero + valGenAngSdvig;  // диапахзон графиков по обобщ координате плюс минус
	const int nBuffRows = 2. *VAlDiap/ valStep ;
	const int nBuffCols = 3;
	double  *parrBuff = new double [nBuffRows  * nBuffCols] ;
	memset(parrBuff, 0, nBuffRows * nBuffCols * sizeof(double));

	const int lenName = 30;
	wchar_t *wcharrFileNames  = new wchar_t [30 * nBuffCols];
	memset(wcharrFileNames, 0, 30 * nBuffCols * sizeof(wchar_t));

	wcscpy( wcharrFileNames, L"GenAng");
	wcscpy( &wcharrFileNames[30], L"FSum");
	wcscpy( &wcharrFileNames[60], L"Polinom4");

	double valF0 = fncSumDiagr( 0.) ;
	double root1 = findFirstZero_For_SumDiagr();
	double root2 = findSecondZero_For_SumDiagr();
   for (int i=0 ; i < nBuffRows; i++)
  {
   double valGenAngCur = VAlDiap - ((double)i ) * valStep;
   parrBuff[ i * nBuffCols] =  valGenAngCur;
   parrBuff[ i * nBuffCols +1] =fncSumDiagr( valGenAngCur)/valF0 ;
   parrBuff[ i * nBuffCols +2] = ( valGenAngCur * valGenAngCur -  root1 * root1)
		* (valGenAngCur * valGenAngCur -  root2 * root2)/ (root1 * root1 * root2 * root2);

  }


  double *pscaley = new double [nBuffCols] ;
  pscaley[0] = 1000.;
  pscaley[1] = 1.;
  pscaley[2] = 1.;

  for (int i=1 ; i < nBuffCols; i++)
  {

  TYrWriteShapeFile::WriteOneReport(Fold  // путь к папке
								  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  , nBuffCols // - к-во переменных о корорых накоплена информация в буфере
								  , nBuffRows //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,lenName // максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i  // номер переменной по оси Y
								  ,1. //  масштаб по оси X
								  ,pscaley[i]// масштаб по оси Y
								   ) ;
  }




	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  Fold);
	wcscat(wchAxesFileName, L"AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-VAlDiap, VAlDiap + 0.5
	 ,-2.,2., 0.25) ;
	 delete pscaley ;
	 delete  parrBuff;
	 delete wcharrFileNames;

}

// графики для сравнения парциальной  диаграммы sinx/x и ее аппроксимации полинрмом 4,6 и 8 степени
void TParAnt::createGraphs_for_Compare_PartDiagrams(wchar_t *Fold )
{
	double valStep =0.01;
	const double VAlDiap = 4. * M_PI;  // диапахзон графиков по обобщ координате плюс минус
	const int nBuffRows = 2. *VAlDiap/ valStep ;
	const int nBuffCols = 7;
	double  *parrBuff = new double [nBuffRows  * nBuffCols] ;
	memset(parrBuff, 0, nBuffRows * nBuffCols * sizeof(double));

	const int lenName = 30;
	wchar_t *wcharrFileNames  = new wchar_t [30 * nBuffCols];
	memset(wcharrFileNames, 0, 30 * nBuffCols * sizeof(wchar_t));

	wcscpy( wcharrFileNames, L"GenAng");
	wcscpy( &wcharrFileNames[30], L"sinc");
	wcscpy( &wcharrFileNames[60], L"Polinom4");
	wcscpy( &wcharrFileNames[90], L"Polinom6");
	wcscpy( &wcharrFileNames[120], L"Polinom8");
	wcscpy( &wcharrFileNames[150], L"Polinom10");
	wcscpy( &wcharrFileNames[180], L"Polinom12");


   for (int i=0 ; i < nBuffRows; i++)
  {
   double valGenAngCur = VAlDiap - ((double)i ) * valStep;
   parrBuff[ i * nBuffCols] =  valGenAngCur;
   parrBuff[ i * nBuffCols +1] = fncDiagrSinx_div_x(valGenAngCur) ;
   parrBuff[ i * nBuffCols +2] = ( valGenAngCur * valGenAngCur -  M_PI * M_PI)
		* ( valGenAngCur * valGenAngCur -  4. *M_PI * M_PI)/(M_PI * M_PI * 4. *M_PI * M_PI);
   parrBuff[ i * nBuffCols +3] = -( valGenAngCur * valGenAngCur -  M_PI * M_PI)
		* ( valGenAngCur * valGenAngCur -  4. *M_PI * M_PI)* ( valGenAngCur * valGenAngCur -  9. *M_PI * M_PI)
		/(M_PI * M_PI * 4. *M_PI * M_PI * 9. *M_PI * M_PI);
   parrBuff[ i * nBuffCols +4] = ( valGenAngCur * valGenAngCur -  M_PI * M_PI)
		* ( valGenAngCur * valGenAngCur -  4. *M_PI * M_PI)* ( valGenAngCur * valGenAngCur -  9. *M_PI * M_PI)
		* ( valGenAngCur * valGenAngCur -  16. *M_PI * M_PI)/(M_PI * M_PI * 4. *M_PI * M_PI * 9. *M_PI * M_PI *16. *M_PI * M_PI);
   parrBuff[ i * nBuffCols +5] = -parrBuff[ i * nBuffCols +4]
		* ( valGenAngCur * valGenAngCur -  25. *M_PI * M_PI)/(M_PI * M_PI * 25. );
   parrBuff[ i * nBuffCols +6] = -parrBuff[ i * nBuffCols +5]
		* ( valGenAngCur * valGenAngCur -  36. *M_PI * M_PI)/(M_PI * M_PI * 36. );


  }

 // double scalex = 100.;
  double *pscaley = new double [nBuffCols] ;
  pscaley[0] = 1000.;
  pscaley[1] = 1.;
  pscaley[2] = 1.;
  pscaley[3] = 1.;
  pscaley[4] = 1.;
  pscaley[5] = 1.;
  pscaley[6] = 1.;
  for (int i=1 ; i < nBuffCols; i++)
  {

  TYrWriteShapeFile::WriteOneReport(Fold  // путь к папке
								  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  , nBuffCols // - к-во переменных о корорых накоплена информация в буфере
								  , nBuffRows //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,lenName // максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i  // номер переменной по оси Y
								  ,1. //  масштаб по оси X
								  ,pscaley[i]// масштаб по оси Y
								   ) ;
  }




	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  Fold);
	wcscat(wchAxesFileName, L"AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-VAlDiap, VAlDiap + 0.5
	 ,-2.,2., 0.25) ;
	 delete pscaley ;
	 delete  parrBuff;
	 delete wcharrFileNames;

}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//////  ВЫЧИСЛИТЕЛЬНЫЕ ФУНКЦИИ ДЛЯ МЕТОДА 3 ПАРЦИАЛЬНЫХ ДИАГРАММ  /////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// оценивание обобщенных углов цели и антипода по измерениям 3 ПАРЦИАЛЬНЫХ диаграмм
//INPUT:
// cmparrS[6] -  массив замеров
// NumRayTriple-  номер тройки рабочих лучей (нумерация лучей идет снизу начиная с нуля)
//                номер тройки это порядковый номер нижнего луча
//  OUTPUT:
//  *valEstAngTarg -   угол цели, рад
// *valEstAngAntp -    угол антипода, рад
// cmpKTarg - коэф отражения цели
// cmpKAntp - коэф отражения антипода
// arrMtrxCorr  - коррел матрица ошибок измерения угла цели и антипода
int TParAnt::estimateMethThreePartDiagr( TComp *cmparrS
  , int iNumRayTriple , double *valEstAngTarg, double *valEstAngAntp, TComp *cmpKTarg , TComp *cmpKAntp , double *arrMtrxCorr
  , double *pval_b0, double *pval_b1 )
{
	double valGenAngSdvig = transformAngToGeneralizedAng (mAngSdvig ) ;
	// нахождение обобщенных углов цели и антипода
	double  valGenTargEps = -100.,  valGenAntpEps = -100.0;  // обобщ углы цели и антипода
	if (EstGenAngsThreePartDiagr(  &cmparrS[iNumRayTriple]
	,  &valGenTargEps,  &valGenAntpEps, pval_b0, pval_b1) !=2)
	{
	   return -1;
	}

	// вычисление коэффиц отражения
	TComp cmparrA[4],cmparrAInv[4], cmparrK[2], cmparrSTemp[2];

	cmparrA[0] = TComp(fncDiagrSinx_div_x(valGenTargEps + valGenAngSdvig), 0.);
	cmparrA[1] = TComp(fncDiagrSinx_div_x (valGenAntpEps + valGenAngSdvig ), 0.);
	cmparrA[2] = TComp(fncDiagrSinx_div_x (valGenTargEps - valGenAngSdvig ), 0.);
	cmparrA[3] = TComp(fncDiagrSinx_div_x (valGenAntpEps - valGenAngSdvig ), 0.);

	bool brez =   InverseMtrx2( cmparrA, cmparrAInv);
	if (!brez)
	{
	 return -2;
	}

	cmparrSTemp[0] = cmparrS [iNumRayTriple ];
	cmparrSTemp[1] = cmparrS [iNumRayTriple + 2 ];
	MtrxMultMatrx(cmparrAInv,2, 2, cmparrSTemp,1, cmparrK);
	*cmpKTarg = cmparrK[0] ;
	*cmpKAntp = cmparrK[1] ;

	///

	// вычисление углоы цели и антипода в рад
		  // относиьтельный угол цетрального луча тройки в радианах
			double valRadCentreRayPos = calcRadAngSdvigaPartDiagr(iNumRayTriple + 1);
			// ( - ((double)mQuantDiagr -1.) /2. + 1. + ((double)iNumRayTriple))				*mAngSdvig ;
		  // вычисление угла цели в радианах в АСК
			const double ValRadTargEps =  transformGeneralizedAngToAng ( valGenTargEps ) ;// в радианы
			*valEstAngTarg = ValRadTargEps + valRadCentreRayPos;
		  // вычисление угла антипода в радианах в АПСК
			const double ValRadAntpАEps =  transformGeneralizedAngToAng (valGenAntpEps ) ;// в радианы
			*valEstAngAntp = ValRadAntpАEps + valRadCentreRayPos;
	///

	// вычисление коррел матрицы ошибок определения обобщенных углов цели и антипода
	double arrMtrxCorrMu[4] = {0.};
	calcMtrxCorrGenAngs_Meth3PartDiagr(iNumRayTriple ,valGenTargEps,  valGenAntpEps
	  ,*cmpKTarg,*cmpKAntp, arrMtrxCorrMu) ;


				///

	// вычисление корреляц матрицы ошибок определения углов цели и антипода в радианаах
	double arrQ[4] = {0.}, arrTemp[4] = {0.};
	double coef = mLambda / mAppert / M_PI;
	arrQ[0] =  coef / sqrt( 1.- coef * valGenTargEps * coef * valGenTargEps);
	arrQ[3] =  coef / sqrt( 1.- coef * valGenAntpEps * coef * valGenAntpEps);
	MtrxMultMatrx(arrQ,2, 2, arrMtrxCorrMu,2, arrTemp)  ;
	MtrxMultMatrxTransp(arrTemp,2, 2, arrQ,2, arrMtrxCorr) ;
	///
	return 0;
}


// оценивание обобщенных углов цели и антипода по измерениям 3 парциальных диаграмм
//INPUT:
// wchFoldName1 - путь к папаке с графиками
// cmparrS[3] -  массив замеров
//  OUTPUT:
//  *pvalGenTargEps -   угол цели, обобщ коорд
// *pvalGenAntpEps -    угол антипода, обобщ коорд
int TParAnt::EstGenAngsThreePartDiagr(TComp *cmparrS ,  double  *pvalGenTargEps
   ,  double  *pvalGenAntpEps, double *pval_b0, double *pval_b1)
{

  // вычисленние вектора b
   double arr_b[2] = {0.}  ;
   bool brez = calcVect_b(cmparrS,arr_b) ;
   ///
   *pval_b0 = arr_b[0];
   *pval_b1 = arr_b[1];


   double arrRoots[2] ={0.};
   int iNUmRoots = findRootsFgr_For_3PartDiagr(arr_b,  arrRoots);
   *pvalGenTargEps =  arrRoots[1];
   *pvalGenAntpEps =  arrRoots[0];
   return iNUmRoots ;
}


// нахождение корней функции FGr длшя метода 3 суммарных диаграмм
// возвращает к-во корней
 int TParAnt::findRootsFgr_For_3PartDiagr(double *arr_b,  double *arrRoots)
 {
	double valGenAngSdvig = transformAngToGeneralizedAng (mAngSdvig ) ;
   int iNUmRoots = 0;

   double valMu0 =  M_PI + valGenAngSdvig - 0.3;
   double valStep = 0.25;
	const int lenBuff = 2. *valMu0/ valStep ;
   double *parrBuff  = new double [lenBuff];

  double valmu = -valMu0;

  for (int i=0 ; i < lenBuff; i++)
  {
   valmu = -valMu0 + ((double)i) * valStep;
   parrBuff[ i ] = fncFGr_PartDiagr(arr_b, valmu);

  }

  for (int i = 0; i < (lenBuff-1); i++)
  {
   if (parrBuff[ i ] * parrBuff[ i + 1] <= 0.)
   {
	 if (iNUmRoots >1)
	 {
	   break;
	 }
	 double valX0 = -valMu0 + ((double)i) * valStep;
	 double valX1 = -valMu0 + ((double) (i +1)) * valStep;
	 arrRoots [iNUmRoots ] = findRootMethChord_For_Fgr_For_3PartDiagr(  arr_b,  valX0, valX1) ;
	 iNUmRoots++;


   }
  }

/*  if ((fabs (arr_b[0])<=700.)&& (fabs (arr_b[1])<=700.) )
	 {
	 int numRow =  (arr_b[0] + 700.)/2.;
	 int numCol =  (arr_b[1] + 700.)/2.;
	 double temp0 =  ArrTargRoots [numRow * 700 + numCol] + (ArrTargRoots [numRow * 700 + numCol + 1] - ArrTargRoots [numRow * 700 + numCol])
				   /2. * (arr_b[1] + 350. -  ((double)numCol) )
				   + (ArrTargRoots [(numRow+1) * 700 + numCol ] - ArrTargRoots [numRow * 700 + numCol])
				   /2. * (arr_b[0] + 350. -  ((double)numRow) );
	 int iii=0;
	 }  */
  delete parrBuff;
  return iNUmRoots;

 }


double TParAnt::fncFGr_PartDiagr(double *arr_b, double valmu)
{
  double valGenAngSdvig = transformAngToGeneralizedAng (mAngSdvig ) ;
  return (fncDiagrSinx_div_x(valmu) - arr_b[0] * fncDiagrSinx_div_x(valmu +valGenAngSdvig)
		 - arr_b[1] * fncDiagrSinx_div_x( valmu - valGenAngSdvig));
}


double TParAnt::findRootMethChord_For_Fgr_For_3PartDiagr(double *arr_b, double  valX0, double valX1)
{
  double valMu0 = valX0 - 0.00001;
  double valMu1 = valX1 + 0.00001;
  double valMuRez =  valMu0 ;
  double eps = 0.001;
  for (int i = 0; i < 100; i++)
  {
	double f0 = fncFGr_PartDiagr(arr_b, valMu0 ) ;
	double f1 = fncFGr_PartDiagr(arr_b, valMu1 ) ;
	double valMuRez0 = valMu0 - f0 * ( valMu1 - valMu0)/ (f1 - f0);
	if (fabs(valMuRez0 - valMuRez) < eps)
	{
	 return valMuRez0;
	}
	double f2 = fncFGr_PartDiagr(arr_b, valMuRez0 ) ;
	if( f2* f0 < 0.)
	{
	   valMu1 = valMuRez0;

	}
	else
	{
      valMu0 = valMuRez0;
    }
	 valMuRez = valMuRez0;
  }
	return -100000.;
}


// вычисление коррел матрицы ошибок оценивания обобщенных углов цели и антипода
// INPUT:
// valGenTargEps, valGenAntpEps - оценки обобщенных углов цели и антипода относительно центральной диаграммы тройки
// valGenAngSdvig  -обобщ уголл сдвига диаграмм
// cmpKTarg, cmpKAntp - оценки  комплексных сигналов цели и антипода
// valNoiseSkz - скз внутр шума антеныы
// mAmplFactSig  - скз шума коэффиц усиления (в долях)
// OUTPUT:
// arrMtrxCorrGenAngs[4] - коррел матрица ошибок оценивания углов цели и антоипода
bool TParAnt::calcMtrxCorrGenAngs_Meth3PartDiagr(const int NumAnsamble, double valGenTargEps,  double valGenAntpEps
  ,TComp cmpKTarg ,TComp cmpKAntp, double *arrMtrxCorrGenAngs )
{
	 double valGenAngSdvig = transformAngToGeneralizedAng (mAngSdvig ) ;

  // формирование матирицы частных производных dG/dx
  double arrdG_po_dX[36] = {0.};
  arrdG_po_dX[0] =  cmpKTarg.m_Re * fncDerivDiagrSinx_div_x(  valGenAngSdvig + valGenTargEps);
  arrdG_po_dX[1] =  cmpKAntp.m_Re * fncDerivDiagrSinx_div_x(  valGenAngSdvig + valGenAntpEps);
  arrdG_po_dX[2] =  fncDiagrSinx_div_x( valGenAngSdvig + valGenTargEps);
  arrdG_po_dX[4] =  fncDiagrSinx_div_x( valGenAngSdvig + valGenAntpEps);

  arrdG_po_dX[6] =  cmpKTarg.m_Im * fncDerivDiagrSinx_div_x(  valGenAngSdvig + valGenTargEps);
  arrdG_po_dX[7] =  cmpKAntp.m_Im * fncDerivDiagrSinx_div_x(  valGenAngSdvig + valGenAntpEps);
  arrdG_po_dX[9] =  fncDiagrSinx_div_x( valGenAngSdvig + valGenTargEps);
  arrdG_po_dX[11] =  fncDiagrSinx_div_x( valGenAngSdvig + valGenAntpEps);

  arrdG_po_dX[12] =  cmpKTarg.m_Re * fncDerivDiagrSinx_div_x(   valGenTargEps);
  arrdG_po_dX[13] =  cmpKAntp.m_Re * fncDerivDiagrSinx_div_x(  valGenAntpEps);
  arrdG_po_dX[14] =  fncDiagrSinx_div_x(  valGenTargEps);
  arrdG_po_dX[16] =  fncDiagrSinx_div_x(  valGenAntpEps);

  arrdG_po_dX[18] =  cmpKTarg.m_Im * fncDerivDiagrSinx_div_x(   valGenTargEps);
  arrdG_po_dX[19] =  cmpKAntp.m_Im * fncDerivDiagrSinx_div_x(   valGenAntpEps);
  arrdG_po_dX[21] =  fncDiagrSinx_div_x(  valGenTargEps);
  arrdG_po_dX[23] =  fncDiagrSinx_div_x(  valGenAntpEps);

  arrdG_po_dX[24] =  cmpKTarg.m_Re * fncDerivDiagrSinx_div_x(  -valGenAngSdvig + valGenTargEps);
  arrdG_po_dX[25] =  cmpKAntp.m_Re * fncDerivDiagrSinx_div_x(  -valGenAngSdvig + valGenAntpEps);
  arrdG_po_dX[26] =  fncDiagrSinx_div_x( -valGenAngSdvig + valGenTargEps);
  arrdG_po_dX[28] =  fncDiagrSinx_div_x( -valGenAngSdvig + valGenAntpEps);

  arrdG_po_dX[30] =  cmpKTarg.m_Im * fncDerivDiagrSinx_div_x(  -valGenAngSdvig + valGenTargEps);
  arrdG_po_dX[31] =  cmpKAntp.m_Im * fncDerivDiagrSinx_div_x( -valGenAngSdvig + valGenAntpEps);
  arrdG_po_dX[33] =  fncDiagrSinx_div_x( -valGenAngSdvig + valGenTargEps);
  arrdG_po_dX[35] =  fncDiagrSinx_div_x( -valGenAngSdvig + valGenAntpEps);
  ///

  // формирование корреляц матрицы измерений
  double arrCorrPartDiagr[36] = {0.};

  createMtrxCorrForMeasuresForAnsambleOfPartDiagrs( valGenTargEps,  valGenAntpEps
  , cmpKTarg , cmpKAntp ,NumAnsamble, 3
  , arrCorrPartDiagr );


  ///

  // вычисление матрицы dX/dS
  double arrMtrx_dX_po_dS[36] = {0.};
  bool brez = InverseMtrx6(arrdG_po_dX,arrMtrx_dX_po_dS);
  if (!brez)
  {
	return false;
  }
  ///

  // вычисление коррел матрицы вектора дельтах
  double arrMtrxCorrX[36] ={0.}, arrtemp0[36] = {0.};
  MtrxMultMatrx(arrMtrx_dX_po_dS,6, 6, arrCorrPartDiagr, 6, arrtemp0) ;
  MtrxMultMatrxTransp(arrtemp0,6, 6, arrMtrx_dX_po_dS,6, arrMtrxCorrX) ;
  ///

  // формирование корреляц матрицы ошибок обобщенных углов
  arrMtrxCorrGenAngs[0]  = arrMtrxCorrX[0] ;
  arrMtrxCorrGenAngs[1]  = arrMtrxCorrX[1] ;
  arrMtrxCorrGenAngs[2]  = arrMtrxCorrX[6] ;
  arrMtrxCorrGenAngs[3]  = arrMtrxCorrX[7] ;
  return true;

}

// построение графика функции F греческое для метода 3 парциальных диаграмм в зависимости от обобщенного угла
void TParAnt::createGraphFGreece_For_3PartDiagr_from_GenAng(wchar_t *Fold, double  valSmeshenieGoriz
	, double  valSMeshenieVert, double val_b0, double val_b1 )
{

    double valGenAngSdvig = transformAngToGeneralizedAng (mAngSdvig ) ;
	double arr_b[2] = {0.};
	arr_b[0] =  val_b0 ;
	arr_b[1]  = val_b1;
   // построение графика функции FGreece

   double valMu0 =  M_PI + valGenAngSdvig - 0.1;
   double valStep = 0.00001;
	const int nBuffRows = 2. *valMu0/ valStep ;
	const int nBuffCols = 2;
	double  *parrBuff = new double [nBuffRows  * nBuffCols] ;
	memset(parrBuff, 0, nBuffRows * nBuffCols * sizeof(double));

	const int lenName = 30;
	wchar_t *wcharrFileNames  = new wchar_t [30 * nBuffCols];
	memset(wcharrFileNames, 0, 30 * nBuffCols * sizeof(wchar_t));

	wcscpy( wcharrFileNames, L"Mu");
	wcscpy( &wcharrFileNames[30], L"FGr_PartDiagr");


  double valmu = -valMu0;
  for (int i=0 ; i < nBuffRows; i++)
  {
   valmu = -valMu0 + ((double)i) * valStep;
   parrBuff[ i * nBuffCols] =  valmu  + valSmeshenieGoriz;
   parrBuff[ i * nBuffCols + 1] = fncFGr_PartDiagr(arr_b,  valmu) + valSMeshenieVert;

  }

 // double scalex = 100.;
  double *pscaley = new double [nBuffCols] ;
  pscaley[0] = 1.;
  pscaley[1] = 1.;
 // wchar_t wchFoldName[] = L"E:\\PROJECTS_C++\\НЛЦ-5П10-03\\Новая папка\\GraphFGr\\";
  for (int i=1 ; i < nBuffCols; i++)
  {

  TYrWriteShapeFile::WriteOneReport(Fold  // путь к папке
								  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  , nBuffCols // - к-во переменных о корорых накоплена информация в буфере
								  , nBuffRows //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,lenName // максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i  // номер переменной по оси Y
								  ,1. //  масштаб по оси X
								  ,pscaley[i]// масштаб по оси Y
								   ) ;
  }



	TURPointXY pntSdvig(valSmeshenieGoriz, valSMeshenieVert);
	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  Fold);
	wcscat(wchAxesFileName, L"\\AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-13., 13.5
	 ,-13.,13.5, -parrBuff[ 0] * 0.05, pntSdvig) ;
	 delete pscaley ;
	 delete  parrBuff;
	 delete wcharrFileNames;

}


// построение основной картинки презентации
// INPUT
//  valScaleY  - растяжение пол оси Y
//
//
//
//
//


void TParAnt::createPtesentGraphs_for_PartDiagr( wchar_t *wchrPresntPartDiagrams, int iNumPartRayTriple
		,double  alfUMTrg,TComp cmpKTarg,double  alfUMAntp, TComp cmpKAntp
		,double  valEstRadAngTarg,double  valEstRadAngAntp
		  ,TComp cmpKTargEst , TComp cmpKAntpEst, double val_b0, double val_b1, double valSMeshenieVert
		  , double valScaleY )
{
   // построение 6 суммарных ьдиаграмм в зависимости от угла в рад
   wchar_t Fold[400] ={0.};
   wcscpy( Fold, wchrPresntPartDiagrams);
   wcscat(Fold, L"\\7_PartDiagr");
	_wmkdir(Fold);

   for (int i =0; i < mQuantDiagr; i++)
   {
	wchar_t FileName[400] ={0.};
	wcscpy( FileName, Fold);
	wcscat(FileName, L"\\PartDiagrNo_");
	wchar_t string[10] = {0.};
	_itow(i,string, 10);
	wcscat(FileName, string);
	wcscat(FileName, L".shp");
	double valSmeshenieGoriz = calcRadAngSdvigaPartDiagr(i);
	//( -((double)mQuantDiagr -1.)/2. + ((double) i)) * mAngSdvig;
	double valSMeshenieVert = 0.;
	createGraphPartDiagr_from_rad(FileName, valSmeshenieGoriz, valSMeshenieVert, valScaleY);
   }
   ///

   // построение 3 рабочих суммарных диаграмм в зависимости от угла в рад
   wcscpy( Fold, wchrPresntPartDiagrams);
   wcscat(Fold, L"\\3_WorkingDiagr");
	_wmkdir(Fold);
	for (int i =0; i < 3; i++)
   {
	wchar_t FileName[400] ={0.};
	wcscpy( FileName, Fold);
	wcscat(FileName, L"\\WorkingSumDiagrNo_");
	wchar_t string[10] = {0.};
	_itow(i,string, 10);
	wcscat(FileName, string);
	wcscat(FileName, L".shp");
	double valSmeshenieGoriz = calcRadAngSdvigaPartDiagr(iNumPartRayTriple + i);
	//( -((double)mQuantDiagr -1.) /2. + ((double)(iNumPartRayTriple + i))) * mAngSdvig;

	createGraphPartDiagr_from_rad(FileName, valSmeshenieGoriz, 0., valScaleY ) ;
   }
   ///

   // построение графика функции FGreece в зависимости от угла рад
   wchar_t Fold1[400] ={0.};
   wcscpy( Fold1, wchrPresntPartDiagrams);
   wcscat(Fold1, L"\\GraphFGreece_from_rad_For_PartDiagr");
	_wmkdir(Fold1);
	wcscat(Fold1, L"\\");

   // относиьтельный угол цетрального луча тройки в радианах
   double valRadCentreRayPos = calcRadAngSdvigaPartDiagr(iNumPartRayTriple + 1);
   // ( -((double)mQuantDiagr -1. )/2. + 1. + ((double)iNumPartRayTriple)) * mAngSdvig;

   createGraphFGreece_from_rad_For_3PartDiagr(Fold1 ,valRadCentreRayPos , valSMeshenieVert, val_b0, val_b1)  ;
   ///

   // построение графиков вертикальных линий антипода и цели
	 wchar_t File[400] ={0.};
	// истинная цель
	TURPointXY  pnt1(alfUMTrg * 100., 0.);
	TURPointXY  pnt2(alfUMTrg * 100., cmpKTarg.modul()/ 800. * 5.);
	TURPolyLine plnTargTrue(   pnt1,   pnt2) ;
	wcscpy( File, Fold);
	wcscat(File, L"\\plnPartDiagrTargTrue.shp");
	plnTargTrue.WriteSetSHPFiles(File ,&plnTargTrue, 1) ;

   ///
   // истинный антипод
	TURPointXY  pnt3(alfUMAntp * 100., 0.);
	TURPointXY  pnt4(alfUMAntp * 100., cmpKAntp.modul()/ 800. * 5.);
	TURPolyLine plnAntpTrue(   pnt3,   pnt4) ;
	wcscpy( File, Fold);
	wcscat(File, L"\\plnPartDiagrAntpTrue.shp");
	plnAntpTrue.WriteSetSHPFiles(File ,&plnAntpTrue, 1) ;
	///

	// оценка цели
	TURPointXY  pnt5(valEstRadAngTarg * 100., 0.);
	TURPointXY  pnt6(valEstRadAngTarg * 100., cmpKTargEst.modul()/ 800. * 5.);
	TURPolyLine plnTargEst(   pnt5,   pnt6) ;
	wcscpy( File, Fold);
	wcscat(File, L"\\plnPartDiagrTargEst.shp");
	plnTargEst.WriteSetSHPFiles(File ,&plnTargEst, 1) ;
	///

	// оценка антипода
   TURPointXY  pnt7(valEstRadAngAntp * 100., 0.);
	TURPointXY  pnt8(valEstRadAngAntp * 100., cmpKAntpEst.modul()/ 800. * 5.);
	TURPolyLine plnAntpEst(   pnt7,   pnt8) ;
	wcscpy( File, Fold);
	wcscat(File, L"\\plnPartDiagrAntpEst.shp");
	plnAntpEst.WriteSetSHPFiles(File ,&plnAntpEst, 1) ;
	///

	// вертик линия оценки цели
	TURPointXY  pnt9(valEstRadAngTarg * 100., cmpKTargEst.modul()/ 800. * 5.);
	TURPointXY  pnt10(valEstRadAngTarg * 100., valSMeshenieVert);
	TURPolyLine plnTargEst1(   pnt9,   pnt10) ;
	wcscpy( File, Fold);
	wcscat(File, L"\\plnPartDiagrTargEst1.shp");
	plnTargEst1.WriteSetSHPFiles(File ,&plnTargEst1, 1) ;
	///

	// вертик линия оценки антипода
	TURPointXY  pnt11(valEstRadAngAntp * 100., cmpKAntpEst.modul()/ 800. * 5.);
	TURPointXY  pnt12(valEstRadAngAntp * 100., valSMeshenieVert);
	TURPolyLine plnAntpEst1(   pnt11,   pnt12) ;
	wcscpy( File, Fold);
	wcscat(File, L"\\plnPartDiagrAntpEst1.shp");
	plnAntpEst1.WriteSetSHPFiles(File ,&plnAntpEst1, 1) ;
	///
   //-------------------------------------
	// изображение углов в виде полигонов
	// полигон оценки цели
	TURPointXY  pntTopRight0(valEstRadAngTarg * 100. + 0.025, cmpKTargEst.modul()/ 800. * 5.);
	TURPointXY  pntBottomLeft0(valEstRadAngTarg * 100. - 0.025, 0.);
	TURPolygon plgTargEst(pntTopRight0, pntBottomLeft0);
	wcscpy( File, Fold);
	wcscat(File, L"\\plgPartDiagrTargEst.shp");
	plgTargEst.WriteSetSHPFiles(File ,&plgTargEst, 1) ;
	///

	// полигон оценки антипода
	TURPointXY  pntTopRight1(valEstRadAngAntp * 100. + 0.025, cmpKAntpEst.modul()/ 800. * 5.);
	TURPointXY  pntBottomLeft1(valEstRadAngAntp * 100. - 0.025, 0.);
	TURPolygon plgAntpEst(pntTopRight1, pntBottomLeft1);
	wcscpy( File, Fold);
	wcscat(File, L"\\plgPartDiagrAntpEst.shp");
	plgAntpEst.WriteSetSHPFiles(File ,&plgAntpEst, 1) ;
	///

	// полигон истинная цель
	TURPointXY  pntTopRight2(alfUMTrg * 100.+ 0.025, cmpKTarg.modul()/ 800. * 5.);
	TURPointXY  pntBottomLeft2(alfUMTrg * 100.- 0.025, 0.);
	TURPolygon plgTargTrue(pntTopRight2, pntBottomLeft2);
	wcscpy( File, Fold);
	wcscat(File, L"\\plgPartDiagrTargTrue.shp");
	plgTargTrue.WriteSetSHPFiles(File ,&plgTargTrue, 1) ;

   ///
   //  полигон  истинный антипод
	TURPointXY  pntTopRight3(alfUMAntp * 100.+ 0.025, cmpKAntp.modul()/ 800. * 5.);
	TURPointXY  pntBottomLeft3(alfUMAntp * 100.- 0.025, 0.);
	TURPolygon plgAntpTrue(pntTopRight3, pntBottomLeft3);
	wcscpy( File, Fold);
	wcscat(File, L"\\plgPartDiagrAntpTrue.shp");
	plgAntpTrue.WriteSetSHPFiles(File ,&plgAntpTrue, 1) ;

	///
	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  Fold);
	wcscat(wchAxesFileName, L"\\AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-13, 13.5
	 ,-10.,10., 0.25) ;
}

 // построние графика парциальной диаграммы смещенной на  valSmeshenieGoriz и valSMeshenieVert
void TParAnt::createGraphPartDiagr_from_rad(wchar_t *FileName,double  valSmeshenieGoriz, double  valSMeshenieVert,  double valScaleY )
{

	double step =0.0001;
	double  valGenAngSdvig = transformAngToGeneralizedAng (mAngSdvig  ) ;
	double valZero =  M_PI;
	const double VAlDiap = 3. *valZero;  // диапахзон графиков по обобщ координате плюс минус
	const int nBuffRows = 2. *VAlDiap/ step ;
	 TURPolyLine pln( 1, nBuffRows) ;



  double tet = -VAlDiap;
  for (int i=0 ; i < nBuffRows; i++)
  {
   tet = -VAlDiap + ((double)i) * step;
   pln.Points[i].X = (transformGeneralizedAngToAng (tet) + valSmeshenieGoriz)* 100. ;
   pln.Points[i].Y =fncDiagrSinx_div_x(tet) * valScaleY + valSMeshenieVert;

  }

 pln.WriteSetSHPFiles(FileName, &pln, 1);
}

// построние графика парциальной диаграммы смещенной на  valSmeshenieGoriz и valSMeshenieVert
void TParAnt::createGraphPartDiagr_from_Gen(wchar_t *FileName,double  valSmeshenieGoriz, double  valSMeshenieVert,  double valScaleY )
{

	double step =0.0001;
	double  valGenAngSdvig = transformAngToGeneralizedAng (mAngSdvig  ) ;
	double valZero =  M_PI;
	const double VAlDiap = 3. *valZero;  // диапахзон графиков по обобщ координате плюс минус
	const int nBuffRows = 2. *VAlDiap/ step ;
	 TURPolyLine pln( 1, nBuffRows) ;



  double tet = -VAlDiap;
  for (int i=0 ; i < nBuffRows; i++)
  {
   tet = -VAlDiap + ((double)i) * step;
   pln.Points[i].X = (tet + valSmeshenieGoriz) ;
   pln.Points[i].Y =fncDiagrSinx_div_x(tet) * valScaleY + valSMeshenieVert;

  }

 pln.WriteSetSHPFiles(FileName, &pln, 1);
}

// построение графика функции F греческое для метода 3 парциальных диаграмм в зависимости от угла в рад
void TParAnt::createGraphFGreece_from_rad_For_3PartDiagr(wchar_t *Fold,double  valSmeshenieGoriz
	, double  valSMeshenieVert, double val_b0, double val_b1 )
{
	double arr_b[2] = {0.};
	arr_b[0] =  val_b0 ;
	arr_b[1]  = val_b1;
   // построение графика функции FGreece
   double valGenAngSdvig = transformAngToGeneralizedAng (mAngSdvig ) ;
   double valMu0 =  M_PI + valGenAngSdvig - 0.1;
   double valStep = 0.00001;
	const int nBuffRows = 2. *valMu0/ valStep ;
	const int nBuffCols = 2;
	double  *parrBuff = new double [nBuffRows  * nBuffCols] ;
	memset(parrBuff, 0, nBuffRows * nBuffCols * sizeof(double));

	const int lenName = 30;
	wchar_t *wcharrFileNames  = new wchar_t [30 * nBuffCols];
	memset(wcharrFileNames, 0, 30 * nBuffCols * sizeof(wchar_t));

	wcscpy( wcharrFileNames, L"Mu");
	wcscpy( &wcharrFileNames[30], L"FGr_PartDiagr");


  double valmu = -valMu0;
  for (int i=0 ; i < nBuffRows; i++)
  {
   valmu = -valMu0 + ((double)i) * valStep;
   parrBuff[ i * nBuffCols] = (transformGeneralizedAngToAng (valmu ) + valSmeshenieGoriz) * 100.;
   parrBuff[ i * nBuffCols + 1] = fncFGr_PartDiagr(arr_b, valmu) + valSMeshenieVert;

  }

 // double scalex = 100.;
  double *pscaley = new double [nBuffCols] ;
  pscaley[0] = 100.;
  pscaley[1] = 1.;

  for (int i=1 ; i < nBuffCols; i++)
  {

  TYrWriteShapeFile::WriteOneReport(Fold  // путь к папке
								  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  , nBuffCols // - к-во переменных о корорых накоплена информация в буфере
								  , nBuffRows //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,lenName // максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i  // номер переменной по оси Y
								  ,1. //  масштаб по оси X
								  ,pscaley[i]// масштаб по оси Y
								   ) ;
  }



	TURPointXY pntSdvig(valSmeshenieGoriz, valSMeshenieVert);
	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  Fold);
	wcscat(wchAxesFileName, L"\\AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-13., 13.5
	 ,-13.,13.5, -parrBuff[ 0] * 0.05, pntSdvig) ;
	 delete pscaley ;
	 delete  parrBuff;
	 delete wcharrFileNames;

}

// построение графико СКЗ оши бок оценивания в звисимости от углового расстояния
// INPUT:
// Fold
// mLambda - раб частота
// mAppert - апертура антиенны
// valAngSdvig - угол сдвига между диаграммами в рад
// valSmeshenieGoriz  - смещение графика по горизонтали
// valSMeshenieVert   - смещение графика по вертикали
// cmpKTarg  - комплексный	 сигнал цели
// cmpKAntp  - комплексный	 сигнал антипода
// valNoiseSkz - шум в парциальной диаграмме
// INPUT:
// графики зависимости ско ошибок оценивания цели и антипода в зависимости
// от углового расстояния между целью и антиподом
// цельнаходится в центре центрального луча
// антипод приближается к цели
void TParAnt::createGraphsSKZ_from_AngDiffer_For_3PartDiagr(wchar_t *Fold ,double  valSmeshenieGoriz
	, double  valSMeshenieVert, TComp cmpKTarg, TComp cmpKAntp )
{

   double valGenAngSdvig = transformAngToGeneralizedAng (mAngSdvig ) ;
  // double valMu0 =  3. * findFirstZero_For_SumDiagr_5P10_03(valGenAngSdvig) + valGenAngSdvig ;
   double valStep = 0.001;
	const int nBuffRows = 3. * M_PI/ valStep ;
	const int nBuffCols = 4;
	double  *parrBuff = new double [nBuffRows  * nBuffCols] ;
	memset(parrBuff, 0, nBuffRows * nBuffCols * sizeof(double));

	const int lenName = 30;
	wchar_t *wcharrFileNames  = new wchar_t [30 * nBuffCols];
	memset(wcharrFileNames, 0, 30 * nBuffCols * sizeof(wchar_t));

	wcscpy( wcharrFileNames, L"DeltaAng");
	wcscpy( &wcharrFileNames[30], L"SigmaTarg_3PartDiagr");
	wcscpy( &wcharrFileNames[60], L"SigmaAntp_3PartDiagr");
	wcscpy( &wcharrFileNames[90], L"CoefCorrel_3PartDiagr");

  double valRadAngTarg = mAngSdvig /2.;
  double valGenAngTarg =  transformAngToGeneralizedAng ( valRadAngTarg ) ;
  double valGenAngAntp = 0.;
  double   arrMtrxCorr[4] ={0.}, arrMtrxCorrMu[4] = {0.};

  for (int i=0 ; i < nBuffRows; i++)
  {
   valGenAngAntp = valGenAngTarg - ((double)i + 10.) * valStep;


   	// вычисление коррел матрицы ошибок определения обобщенных углов цели и антипода
	calcMtrxCorrGenAngs_Meth3PartDiagr(2, valGenAngTarg ,  valGenAngAntp
	  , cmpKTarg, cmpKAntp, arrMtrxCorrMu) ;
	///

	// вычисление корреляц матрицы ошибок определения углов цели и антипода в радианаах
	double arrQ[4] = {0.}, arrTemp[4] = {0.};
	double coef = mLambda / mAppert / M_PI;
	arrQ[0] =  coef / sqrt( 1.- coef * valGenAngTarg * coef * valGenAngTarg);
	arrQ[3] =  coef / sqrt( 1.- coef *  valGenAngAntp * coef *  valGenAngAntp);
	MtrxMultMatrx(arrQ,2, 2, arrMtrxCorrMu,2, arrTemp)  ;
	MtrxMultMatrxTransp(arrTemp,2, 2, arrQ,2, arrMtrxCorr) ;


   parrBuff[ i * nBuffCols] = valRadAngTarg - transformGeneralizedAngToAng (valGenAngAntp);

   parrBuff[ i * nBuffCols + 1] =  sqrt(arrMtrxCorr[0])* 1000. ;
   parrBuff[ i * nBuffCols + 2] = sqrt(arrMtrxCorr[3])* 1000. ;
   if ((arrMtrxCorr[0] >0.) && (arrMtrxCorr[3] >0.))
   {
   parrBuff[ i * nBuffCols + 3] = arrMtrxCorr[1] / sqrt(arrMtrxCorr[0])/sqrt(arrMtrxCorr[3]) ;

   }
   else
   {
	   parrBuff[ i * nBuffCols + 3] = 0.;
   }

  }

 // double scalex = 100.;
  double *pscaley = new double [nBuffCols] ;
  pscaley[0] = 1000.;
  pscaley[1] = 1.;
  pscaley[2] = 1.;
  pscaley[3] = 1.;
  for (int i=1 ; i < nBuffCols; i++)
  {

  TYrWriteShapeFile::WriteOneReport(Fold  // путь к папке
								  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  , nBuffCols // - к-во переменных о корорых накоплена информация в буфере
								  , nBuffRows //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,lenName // максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i  // номер переменной по оси Y
								  ,100. //  масштаб по оси X
								  ,pscaley[i]// масштаб по оси Y
								   ) ;
  }



	TURPointXY pntSdvig(valSmeshenieGoriz, valSMeshenieVert);
	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  Fold);
	wcscat(wchAxesFileName, L"AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-13., 13.5
	 ,-13.,13.5, -parrBuff[ 0] * 0.05, pntSdvig) ;

	 TURPointXY pntSdvig1(valSmeshenieGoriz, 0.);
	wchar_t wchAxesFileName1[300] ={0};
	wcscpy(  wchAxesFileName1,  Fold);
	wcscat(wchAxesFileName1, L"AxesArr1.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName1,-13., 13.5
	 ,-13.,13.5, -parrBuff[ 0] * 0.05, pntSdvig1) ;

	 // построение графиков 3-х суммарных диаграмм снизу

	wcscpy(  wchAxesFileName,  Fold);
	wcscat(wchAxesFileName, L"DigarPart1.shp");
	 createGraphPartDiagr_from_rad(wchAxesFileName,  valSmeshenieGoriz,   valSMeshenieVert,  4. );

		wcscpy(  wchAxesFileName,  Fold);
	wcscat(wchAxesFileName, L"DigarPart0.shp");
	 createGraphPartDiagr_from_rad(wchAxesFileName, valSmeshenieGoriz -mAngSdvig,   valSMeshenieVert,  4. );

		wcscpy(  wchAxesFileName,  Fold);
	wcscat(wchAxesFileName, L"DigarPart2.shp");
	 createGraphPartDiagr_from_rad(wchAxesFileName, valSmeshenieGoriz+ mAngSdvig,   valSMeshenieVert,  4. );

	// полигон изображающий истинная цель
	TURPointXY  pntTopRight2(0.05, cmpKTarg.modul()/ 800. * 5.);
	TURPointXY  pntBottomLeft2(- 0.05, 0.);
	TURPolygon plgTargTrue(pntTopRight2, pntBottomLeft2);
	wcscpy( wchAxesFileName, Fold);
	wcscat(wchAxesFileName, L"\\plgTargTrue.shp");
	plgTargTrue.WriteSetSHPFiles(wchAxesFileName ,&plgTargTrue, 1) ;
	delete pscaley ;
	delete  parrBuff;
	delete wcharrFileNames;

}

// вычисление угла ценнтральной оси парциальной диаграммы в радианах с номером numPartDiagr
// номер диаграммы исчисляется от 0 до (mQuantDiagr -1)
double  TParAnt::calcRadAngSdvigaPartDiagr(const int numPartDiagr)
{
	return (- (((double) mQuantDiagr)-1.)/2. + (double)numPartDiagr) *   mAngSdvig;
}

// вычисление угла ценнтральной оси парциальной диаграммы в обобщенном угле с номером numPartDiagr
// номер диаграммы исчисляется от 0 до (mQuantDiagr -1)
double  TParAnt::calGenAngSdvigaPartDiagr(const int numPartDiagr)
{
	double valRad = calcRadAngSdvigaPartDiagr( numPartDiagr);
	return transformAngToGeneralizedAng (valRad  ) ;;
}

// вычисление угла ценнтральной оси суцммарной диаграммы в радианах с номером numSumDiagr
// номер диаграммы исчисляется от 0 до (mQuantDiagr -2)
double  TParAnt::calcRadAngSdvigaSumDiagr(const int numSumDiagr)
{
    // угол центральной оси нулевой парциальной диаграммы
	double valPart0_Rad = calcRadAngSdvigaPartDiagr( 0);
	///
	// угол центральной оси нулевой суммарной  диаграммы
	double valSum0_Rad =valPart0_Rad + mAngSdvig/ 2.;
	///
	return valSum0_Rad + ((double)numSumDiagr) *   mAngSdvig;
}

// вычисление угла ценнтральной оси суммарной диаграммы в обобщенном угле с номером numSumDiagr
// номер диаграммы исчисляется от 0 до (mQuantDiagr -2)
double  TParAnt::calGenAngSdvigaSumDiagr(const int numSumDiagr)
{
	double valRad = calcRadAngSdvigaSumDiagr( numSumDiagr);
	return transformAngToGeneralizedAng (valRad  ) ;
}

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
////// ФУНКЦИИ АМПЛИТУДНОГО РАЗНОСТНО-СУММАРНОГО МЕТОДА ///////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
void TParAnt::createGraphs_For_BearingFuncs( wchar_t *wchFoldName1)
{
double valGenAngSdvig = transformAngToGeneralizedAng (mAngSdvig ) ;
 // графики для парциальных диаграмм
	wchar_t wchFoldName[300] ={0};
	wcscpy(  wchFoldName,  wchFoldName1);
	wcscat(wchFoldName, L"\\");
	double step =0.0001;

	const double VAlDiap = M_PI;  // диапазон графиков по обобщ координате плюс минус
								   // для метода парциальных диаграмм
	const int nBuffRows = 2. *VAlDiap/ step ;
	const int nBuffCols = 5;
	double  *parrBuff = new double [nBuffRows  * nBuffCols] ;
	memset(parrBuff, 0, nBuffRows * nBuffCols * sizeof(double));

	const int lenName = 30;
	wchar_t *wcharrFileNames  = new wchar_t [30 * nBuffCols];
	memset(wcharrFileNames, 0, 30 * nBuffCols * sizeof(wchar_t));

	wcscpy( wcharrFileNames, L"GeneralizedTetta");

	wcscpy( &wcharrFileNames[ 30], L"fncBear_Part_ARSM");

	wcscpy( &wcharrFileNames[2 * 30], L"fncBear_Part_Division");
	wcscpy( &wcharrFileNames[3 * 30], L"fncFirstApprox_Bear_Part_ARSM");
	wcscpy( &wcharrFileNames[4* 30], L"fncThirdApprox_Bear_Part_ARSM");

  double valCoefTalor1 =  calcFirstCoefTalorBearFunc_For_PartDiagr();
  double valCoefTalor3 = calcThirdCoefTalorBearFunc_For_PartDiagr() ;
  double tet = -VAlDiap;
  for (int i=0 ; i < nBuffRows; i++)
  {
   tet = -VAlDiap + ((double)i) * step;
   double temp0 = fncDiagrSinx_div_x( tet - valGenAngSdvig/2.);
   double temp1 = fncDiagrSinx_div_x( tet + valGenAngSdvig/2.);

   parrBuff[ i * nBuffCols] = tet ; // обобщ угол
   parrBuff[ i * nBuffCols + 1] = (temp0 - temp1)/ (temp0 + temp1);
   parrBuff[ i * nBuffCols + 2] =  temp0/ temp1 -1.;
   parrBuff[ i * nBuffCols + 3] =  -tet * valCoefTalor1;
   parrBuff[ i * nBuffCols + 4] =  -tet * valCoefTalor1 + valCoefTalor3* tet* tet*tet;



  }

 // double scalex = 100.;
  double *pscaley = new double [nBuffCols] ;
  pscaley[0] = 1.;
  pscaley[1] = 1.;
  pscaley[2] = 1.;
  pscaley[3] = 1.;
  pscaley[4] = 1.;
//  pscaley[5] = 10.;
//  pscaley[6] = 10.;
//  pscaley[7] = 10.;
 // pscaley[6] = 1;
 // pscaley[7] = 1;
  for (int i=1; i < nBuffCols; i++)
  {

  TYrWriteShapeFile::WriteOneReport(                 wchFoldName  // путь к папке
								  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  , nBuffCols // - к-во переменных о корорых накоплена информация в буфере
								  , nBuffRows //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,lenName // максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i  // номер переменной по оси Y
								  ,1. //  масштаб по оси X
								  ,1.// масштаб по оси Y
								   ) ;
  }

  TURPointXY  arrPoints[2];
  arrPoints[0] =  TURPointXY ( -VAlDiap, 0.);
  arrPoints[1] =  TURPointXY (  VAlDiap, 0.);
  TURPolyLine plnDiap0(arrPoints,2) ;
  wchar_t wchFileName[300] ={0};
  wcscpy(  wchFileName,  wchFoldName);
  wcscat(wchFileName, L"Diap0.shp");
  TURPolyLine::WriteSetSHPFiles(wchFileName, &plnDiap0, 1) ;
  delete parrBuff;
  delete pscaley;
   ///

  // графики для суммарных диаграмм
	const double VAlDiap1 =  findFirstZero_For_SumDiagr();  // диапазон графиков по обобщ координате плюс минус
								   // для метода сум  диаграмм
	const int nBuffRows1 = 2. *VAlDiap1/ step ;
	const int nBuffCols1 = 3;
	double  *parrBuff1 = new double [nBuffRows1  * nBuffCols1] ;
	memset(parrBuff1, 0, nBuffRows1 * nBuffCols1 * sizeof(double));


	wchar_t *wcharrFileNames1  = new wchar_t [30 * nBuffCols1];
	memset(wcharrFileNames1, 0, 30 * nBuffCols1 * sizeof(wchar_t));

	wcscpy( wcharrFileNames1, L"GeneralizedTetta");

	wcscpy( &wcharrFileNames1[ 30], L"fncBear_Sum_ARSM");
	wcscpy( &wcharrFileNames1[ 2 *30], L"fncFirstApprox_For_Bear_Sum_ARSM");


  double valSumDEriv1 =  calcFirstCoefTalorBearFunc_For_SumDiagr();
  tet = -VAlDiap1;
  for (int i=0 ; i < nBuffRows1; i++)
  {
   tet = -VAlDiap1 + ((double)i) * step;
   double temp0 = fncSumDiagr( tet - valGenAngSdvig/2.) ;
   double temp1 = fncSumDiagr( tet + valGenAngSdvig/2.) ;

   parrBuff1[ i * nBuffCols1] = tet ; // обобщ угол
   parrBuff1[ i * nBuffCols1 + 1] = (temp0 - temp1)/ (temp0 + temp1);
   parrBuff1[ i * nBuffCols1 + 2] = -tet * valSumDEriv1;
  }


  double *pscaley1 = new double [nBuffCols1] ;
  pscaley1[0] = 1.;
  pscaley1[1] = 1.;
  pscaley[2] = 1.;
//  pscaley[3] = 10.;
//  pscaley[4] = 10.;
//  pscaley[5] = 10.;
//  pscaley[6] = 10.;
//  pscaley[7] = 10.;
 // pscaley[6] = 1;
 // pscaley[7] = 1;
  for (int i=1; i < nBuffCols1; i++)
  {

  TYrWriteShapeFile::WriteOneReport(                 wchFoldName1  // путь к папке
								  , parrBuff1 // массив с информацией - матрица nBuffRows x nBuffCols
								  , nBuffCols1 // - к-во переменных о корорых накоплена информация в буфере
								  , nBuffRows1 //  - к-во точек
								  ,wcharrFileNames1 //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,lenName // максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i  // номер переменной по оси Y
								  ,1. //  масштаб по оси X
								  ,1.// масштаб по оси Y
								   ) ;
  }

  TURPointXY  arrPoints1[2];
  arrPoints1[0] =  TURPointXY ( -VAlDiap1, 0.);
  arrPoints1[1] =  TURPointXY (  VAlDiap1, 0.);
  TURPolyLine plnDiap1(arrPoints1,2) ;

  wcscpy(  wchFileName,  wchFoldName);
  wcscat(wchFileName, L"Diap1.shp");
  TURPolyLine::WriteSetSHPFiles(wchFileName, &plnDiap1, 1) ;
  delete parrBuff1;
  delete pscaley1;

  ///

	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  wchFoldName);
	wcscat(wchAxesFileName, L"AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-(VAlDiap1 + 0.05), VAlDiap1 + 0.5
	 ,-3. ,3.5, 0.3) ;

  delete wcharrFileNames;
  delete wcharrFileNames1;
}

// нахождение тангенса  гла наклона пеленг функции для суммарных диаграмм
// в обобщенных координатах (первой производной в нуле)
double  TParAnt::calcFirstCoefTalorBearFunc_For_SumDiagr()
{
	double valGenAngSdvig = transformAngToGeneralizedAng (mAngSdvig ) ;
	return fncDerivSumDiagr(valGenAngSdvig/2.) / fncSumDiagr(valGenAngSdvig/2.);
}

// нахождение тангенса угла наклона пеленг функции для парц  диаграмм
// в обобщенных координатах (первой производной в нуле)
double  TParAnt::calcFirstCoefTalorBearFunc_For_PartDiagr()
{
	double valGenAngSdvig = transformAngToGeneralizedAng (mAngSdvig ) ;
	return  fncDerivDiagrSinx_div_x(valGenAngSdvig/2.)/fncDiagrSinx_div_x(valGenAngSdvig/2.) ;
}

// нахождение второй производной  пеленг функции для парц  диаграмм
// в обобщенных координатах в точке  tet
double  TParAnt::calcSecondDerivBearFunc_For_PartDiagr(const double tet)
{
	double bet = transformAngToGeneralizedAng (mAngSdvig )/2.;
	double temp0 =  fncDeriv2DiagrSinx_div_x( tet -bet );
	double temp1 =  fncDeriv2DiagrSinx_div_x( tet +bet );

	double temp2 = fncDerivDiagrSinx_div_x( tet -bet );
	double temp3 = fncDerivDiagrSinx_div_x( tet +bet );

	double temp4 = fncDiagrSinx_div_x( tet -bet );
	double temp5 = fncDiagrSinx_div_x( tet +bet );

	double tempSum =  temp4 +  temp5;

	double temp6 =  (temp0* temp5   - temp1 * temp4) * tempSum;

	double temp7 =   -2.*(temp2 * temp5 - temp3 * temp4) * ( temp2 + temp3);

	double temp8 = temp6 + temp7;

	return 2. * temp8 /(tempSum * tempSum * tempSum);
}

double  TParAnt::calcThirdCoefTalorBearFunc_For_PartDiagr()
{
	double temp = calcSecondDerivBearFunc_For_PartDiagr(0.01)/0.01/6.;
	return  temp ;
}

 //оценка угла цели АРСМ в АСК
 //INPUT:
 // cmparrPartSZv[ mQuantDiagr] - комплексный массив измерений веера диаграмм
 // NumPartRayPare  - порядковый1 номер пары рабочих диаграмм (нумерация начинается с нуля)
 // OUTPUT:
 // *pvalEstPartDiaRadAngTarg - оценка угла цели, рад
 // *pcmpPartDiaKTarg  - оценка сигнала цели
 // *pvalPartDiaDisp - дисперсия ошибки оценки угла цели, рад*рад
 // *pvalPsi2 - нормированная мнимая часть моноимпульсного отношения

int TParAnt::estimateARSMethPartDiagr( TComp *cmparrPartSZv
		  ,const int NumPartRayPare,  double *pvalEstPartDiaRadAngTarg,  TComp *pcmpPartDiaKTarg
		   , double *pvalPartDiaDisp, double *pvalPsi2)
{
		double valGenAngSdvig = transformAngToGeneralizedAng (mAngSdvig ) ;
	// вычисление коэффиц тейлора пеленгационной функции
	double val_a0 =  calcThirdCoefTalorBearFunc_For_PartDiagr(); // при тетта в кубе
	double val_a2 =  -calcFirstCoefTalorBearFunc_For_PartDiagr(); // при тетта
	///

	// вычисление правой части пеленгационного уравнения
	const TComp CmpTemp = (cmparrPartSZv [NumPartRayPare + 1] - cmparrPartSZv [NumPartRayPare] )
					/(cmparrPartSZv [NumPartRayPare + 1] + cmparrPartSZv [NumPartRayPare] );
	double val_a3 = CmpTemp.m_Re;
	///

	double valGenAng = val_a3/val_a2;
	if(abs(valGenAng)>= M_PI)
	{
		return -1;
	}
	// оценка угла цели в обобщ координатах
	double valF = (val_a0 * valGenAng * valGenAng * valGenAng + val_a2 * valGenAng  -val_a3);
	double valFDeriv = (3. * val_a0 * valGenAng * valGenAng  + val_a2 );
	valGenAng = valGenAng - valF / valFDeriv;
	   ///
		// оценка угла цели в радианах
   *pvalEstPartDiaRadAngTarg = transformGeneralizedAngToAng (valGenAng +  calGenAngSdvigaPartDiagr(NumPartRayPare) + valGenAngSdvig/2.) ;
	   ///

	   // оценка сигнала цели
	   double  temp = fncDiagrSinx_div_x(valGenAng + valGenAngSdvig/2.);
	   if (valGenAng < 0.)
	   {
		 *pcmpPartDiaKTarg =   cmparrPartSZv [NumPartRayPare]/ TComp(fncDiagrSinx_div_x(valGenAng + valGenAngSdvig/2.), 0.);
	   }
	   else
	   {
		 *pcmpPartDiaKTarg =   cmparrPartSZv [NumPartRayPare +1]/ TComp(fncDiagrSinx_div_x(valGenAng - valGenAngSdvig/2.), 0.);

	   }
	   ///

	   // формирование корреляц матрицы arrMtrxCorr[4] ошибок измерения рабочих диаграмм
	   double arrMtrxCorr[16] = {0.};
	   double valGenAntpEps = 0.;
	   TComp cmpKAntp(0.,0.);
	   int lenAnsamble =  2;
	   createMtrxCorrForMeasuresForAnsambleOfPartDiagrs(valGenAng,   valGenAntpEps
	   ,*pcmpPartDiaKTarg ,cmpKAntp, NumPartRayPare,lenAnsamble, arrMtrxCorr )  ;
		///

		// вычисление дисперсии ошибки   определения обобщенного угла


		double valS1 = (*pcmpPartDiaKTarg).modul() * fncDiagrSinx_div_x(valGenAng + valGenAngSdvig/2.) ;
		double valS2 = (*pcmpPartDiaKTarg).modul() * fncDiagrSinx_div_x(valGenAng - valGenAngSdvig/2.) ;
		double valSum2= (valS1 + valS2) * (valS1 + valS2);
		double arra[4] = {0.}, arrT[4] ={0.};
		double valDispPsi = 0.;
		arra[0] =  2. *valS2/valSum2;
		arra[2] = -2. * valS1/valSum2;
		MtrxMultMatrx(arra,1, 4, arrMtrxCorr,4, arrT) ;
		MtrxMultMatrx(arrT,1, 4, arra,1, &valDispPsi);

		double valDispGen =  valDispPsi / (valFDeriv*valFDeriv);
		///

		double valSKZ = transformGeneralizedAngToAng (sqrt(valDispGen) ) ;
		 *pvalPartDiaDisp =  valSKZ * valSKZ;
		 *pvalPsi2 = fabs(CmpTemp.m_Im)/sqrt(valDispPsi);
	return 0;
}




//оценка угла цели АРСМ в АСК
 //INPUT:
 // cmparrSumSZv[ mQuantDiagr] - комплексный массив измерений веера диаграмм
 // NumSumRayPare  - порядковый1 номер пары рабочих диаграмм (нумерация начинается с нуля)
 // OUTPUT:
 // *pvalEstSumDiaRadAngTarg - оценка угла цели, рад
 // *pcmpSumDiaKTarg  - оценка сигнала цели
 // *pvalSumDiaDisp - дисперсия ошибки оценки угла цели, рад*рад
 // *pvalPsi2 - нормированная мнимая часть моноимпульсного отношения
int TParAnt::estimateARSMethSumDiagr(TComp *cmparrSumSZv ,const int NumSumRayPare
   ,double *pvalEstSumDiaRadAngTarg, TComp *pcmpSumDiaKTarg ,double *pvalSumDiaDisp, double *pPsi2)
{
	double valGenAngSdvig = transformAngToGeneralizedAng (mAngSdvig ) ;
	// вычисление коэффиц тейлора пеленгационной функции

	double val_a2 =  -calcFirstCoefTalorBearFunc_For_SumDiagr(); //
	///

	// вычисление правой части пеленгационного уравнения
	const TComp CmpTemp = (cmparrSumSZv [NumSumRayPare + 1] - cmparrSumSZv [NumSumRayPare] )
					/(cmparrSumSZv [NumSumRayPare + 1] + cmparrSumSZv [NumSumRayPare] );
	double val_a3 = CmpTemp.m_Re;
	///
	 // оценка угла цели в обобщ координатах
	double valGenAng = val_a3/val_a2;
	if(abs(valGenAng)>= M_PI)
	{
		return -1;
	}

  //	double valF = (val_a0 * valGenAng * valGenAng * valGenAng + val_a2 * valGenAng  -val_a3);
	double valFDeriv = val_a2;
   //	valGenAng = valGenAng - valF / valFDeriv;
	   ///
		// оценка угла цели в радианах
   *pvalEstSumDiaRadAngTarg = transformGeneralizedAngToAng (valGenAng +  calGenAngSdvigaSumDiagr(NumSumRayPare) + valGenAngSdvig/2.) ;
	   ///

	   // оценка сигнала цели

	   if (valGenAng < 0.)
	   {
		 *pcmpSumDiaKTarg =   cmparrSumSZv [NumSumRayPare]/ TComp(fncSumDiagr(valGenAng + valGenAngSdvig/2.), 0.);
	   }
	   else
	   {
		 *pcmpSumDiaKTarg =   cmparrSumSZv [NumSumRayPare +1]/ TComp(fncSumDiagr(valGenAng - valGenAngSdvig/2.), 0.);

	   }

	   ///

	   // формирование корреляц матрицы arrMtrxCorr[4] ошибок измерения рабочих диаграмм
	   double arrMtrxCorr[16] = {0.};
	   double valGenAntpEps = 0.;
	   TComp cmpKAntp(0.,0.);
	   int lenAnsamble =  2;
	   createMtrxCorrForMeasuresForAnsambleOfSumDiagrs(valGenAng,   valGenAntpEps
	   ,*pcmpSumDiaKTarg ,cmpKAntp, NumSumRayPare,lenAnsamble, arrMtrxCorr )  ;
		///

		// вычисление дисперсии ошибки   определения обобщенного угла
		double valS1 = (*pcmpSumDiaKTarg).modul() * fncSumDiagr(valGenAng + valGenAngSdvig/2.) ;
		double valS2 = (*pcmpSumDiaKTarg).modul() * fncSumDiagr(valGenAng - valGenAngSdvig/2.) ;
	  //	double valS1 = cmparrSumSZv [NumSumRayPare].modul() ;
	 //	double valS2 = cmparrSumSZv [NumSumRayPare +1 ].modul() ;
		double valSum2= (valS1 + valS2) * (valS1 + valS2);
        double arra[4] = {0.}, arrT[4] ={0.};
		double valDispPsi = 0.;
		arra[0] =  2. *valS2/valSum2;
		arra[2] = -2. * valS1/valSum2;
		MtrxMultMatrx(arra,1, 4, arrMtrxCorr,4, arrT) ;
		MtrxMultMatrx(arrT,1, 4, arra,1, &valDispPsi);


	   //	const double ValDispPsi = 4. * (valS1 * valS1 * arrMtrxCorr[ 2 * 4 + 2] + valS2 * valS2 * arrMtrxCorr[0])
	//	   / (valSum2*valSum2);
		double valDispGen =  valDispPsi / (valFDeriv*valFDeriv);
		///

		double valSKZ = transformGeneralizedAngToAng (sqrt(valDispGen) ) ;
		 *pvalSumDiaDisp =  valSKZ * valSKZ;
	 *pPsi2 = fabs(CmpTemp.m_Im)/sqrt(valDispPsi);
	return 0;

}

//  создание графиков для результатов имит моделирования АРСМ
void TParAnt::createGraphsRezImitARSM(wchar_t *wchFoldName1, const double AlfRadUMTrg
	   , const TComp CmpKTarg, const int NumSumRayPare, const int NumPartRayPare, const int QuantIspit)
{
		double alfUMAntp = 0.;
		// подготовка временных рабочих массивов   и переменнных
		const TComp CmpKAntp(0.,0.);
	   TComp *cmparrPartS = new TComp[ mQuantDiagr];
		memset(cmparrPartS, 0, mQuantDiagr * sizeof(TComp));

		TComp *cmparrPartSZv = new TComp[ mQuantDiagr];
		memset(cmparrPartSZv, 0, mQuantDiagr * sizeof(TComp));

		TComp * cmparrSumS = new TComp[ mQuantDiagr-1];
		memset( cmparrSumS, 0, (mQuantDiagr -1) * sizeof(TComp));

		TComp * cmparrSumSZv = new TComp[ mQuantDiagr-1];
		memset( cmparrSumSZv, 0, (mQuantDiagr -1) * sizeof(TComp));

		double valEstPartDiaRadAngTarg = 0., valPartDiaDisp = 0., valPartPsi2 = 0.;
		TComp cmpPartDiaKTarg (0.,0.);

		double valEstSumDiaRadAngTarg = 0., valSumDiaDisp = 0., valSumPsi2 = 0.;
		TComp cmpSumDiaKTarg (0.,0.);

		double valAverPart =0.,valAverSquarePart =0., valAverSum =0.,valAverSquareSum =0. ;
		///

		// подготовка директории дл графиков
	wchar_t wchFoldName[300] ={0};
	wcscpy(  wchFoldName,  wchFoldName1);
	wcscat(wchFoldName, L"\\");
	///
	// подготовка буфера информации
	const int nBuffRows = QuantIspit ;
	const int nBuffCols = 13;
	double  *parrBuff = new double [nBuffRows  * nBuffCols] ;
	memset(parrBuff, 0, nBuffRows * nBuffCols * sizeof(double));

	const int lenName = 30;
	wchar_t *wcharrFileNames  = new wchar_t [30 * nBuffCols];
	memset(wcharrFileNames, 0, 30 * nBuffCols * sizeof(wchar_t));
	// номер испытания
	wcscpy( wcharrFileNames, L"i");
	// оценка угла методом парц диаграмм
	wcscpy( &wcharrFileNames[ 30], L"EstAng_PartMeth");
	// оценка СКЗ ошибки измерения угла
	wcscpy( &wcharrFileNames[2 * 30], L"CurSKZ_PartMeth");
	// нормированная оценка мнимой части моноимпульсного отношения  метода парц диаграмм
	wcscpy( &wcharrFileNames[3 * 30], L"Part_PsiIm");
	//  ошибка ищзмерения угла методом парц диаграмм
	wcscpy( &wcharrFileNames[4 * 30], L"Part_Err");
	//  статистическое среднее метода парц диаграмм
	wcscpy( &wcharrFileNames[5 * 30], L"Part_StatMean");
	// статистическое СКЗ метода парц диаграмм
	wcscpy( &wcharrFileNames[6 * 30], L"Part_StatSKZ");

	// оценка угла методом сум диаграмм
	wcscpy( &wcharrFileNames[7 * 30], L"EstAng_SumMeth");
	// оценка СКЗ ошибки измерения угла
	wcscpy( &wcharrFileNames[8 * 30], L"CurSKZ_SumMeth");
	// нормированная оценка мнимой части моноимпульсного отношения  метода сум диаграмм
	wcscpy( &wcharrFileNames[9 * 30], L"Sum_PsiIm");
	//  ошибка ищзмерения угла методом сум диаграмм
	wcscpy( &wcharrFileNames[10 * 30], L"Sum_Err");
	//  статистическое среднее метода сум диаграмм
	wcscpy( &wcharrFileNames[11 * 30], L"Sum_StatMean");
	// статистическое СКЗ метода сум диаграмм
	wcscpy( &wcharrFileNames[12 * 30], L"Sum_StatSKZ");

	// Испытания и накопление буфера
   for (int i =0; i < nBuffRows; i++)
   {
		ImitateMeasureArrayPartialDiagrams(AlfRadUMTrg,CmpKTarg
		, alfUMAntp, CmpKAntp, cmparrPartS,cmparrPartSZv);

		ImitateMeasureArraySumDiagrams( cmparrPartS,cmparrPartSZv,  cmparrSumS, cmparrSumSZv);

		estimateARSMethPartDiagr( cmparrPartSZv, NumPartRayPare,  &valEstPartDiaRadAngTarg,  &cmpPartDiaKTarg
		   , &valPartDiaDisp, &valPartPsi2) ;

		estimateARSMethSumDiagr( cmparrSumSZv, NumSumRayPare,  &valEstSumDiaRadAngTarg,  &cmpSumDiaKTarg
		   , &valSumDiaDisp, &valSumPsi2) ;

		   // заполнение строки буфера
		   parrBuff[ i * nBuffCols] = (double)i;
			// оценка угла методом парц диаграмм
		   parrBuff[ i * nBuffCols + 1 ] = valEstPartDiaRadAngTarg;
		   // оценка СКЗ ошибки измерения угла
		   parrBuff[ i * nBuffCols + 2 ] = sqrt(valPartDiaDisp);
			// нормированная оценка мнимой части моноимпульсного отношения  метода парц диаграмм
		   parrBuff[ i * nBuffCols + 3 ] =  valPartPsi2;
		   //  ошибка ищзмерения угла методом парц диаграмм
		   parrBuff[ i * nBuffCols + 4 ] = valEstPartDiaRadAngTarg -   AlfRadUMTrg;
		   //  стат среднее и среднее квадрата

		   CalcStatParams_(valEstPartDiaRadAngTarg, i+1
			, valAverPart,valAverSquarePart);
		   parrBuff[ i * nBuffCols + 5 ] = valAverPart;  //стат среднее
		   parrBuff[ i * nBuffCols + 6 ] = fncSKO_(valAverPart ,valAverSquarePart); //стат SKO


			// оценка угла методом сум диаграмм
		   parrBuff[ i * nBuffCols + 7 ] = valEstSumDiaRadAngTarg;
		   // оценка СКЗ ошибки измерения угла
		   parrBuff[ i * nBuffCols + 8 ] = sqrt(valSumDiaDisp);
			// нормированная оценка мнимой части моноимпульсного отношения  метода сум диаграмм
		   parrBuff[ i * nBuffCols + 9 ] =  valSumPsi2;
		   //  ошибка ищзмерения угла методом сум диаграмм
		   parrBuff[ i * nBuffCols + 10] = valEstSumDiaRadAngTarg -   AlfRadUMTrg;
		   //  стат среднее и среднее квадрата

		   CalcStatParams_(valEstSumDiaRadAngTarg, i+1
			, valAverSum,valAverSquareSum);
		   parrBuff[ i * nBuffCols + 11] = valAverSum;  //стат среднее
		   parrBuff[ i * nBuffCols + 12] = fncSKO_(valAverSum ,valAverSquareSum); //стат SKO

   }

   ///

   // вывод информации накопленной в буфере в графики

  double *pscaley = new double [nBuffCols] ;
  pscaley[0] = 1.;

  pscaley[1] = 1000.;
  pscaley[2] = 10000.;
  pscaley[3] = 1.;
  pscaley[4] = 10000.;
  pscaley[5] = 1000.;
  pscaley[6] = 10000.;

  pscaley[7] = 1000.;
  pscaley[8] = 10000.;
  pscaley[9] = 1.;
  pscaley[10] = 10000.;
  pscaley[11] = 1000.;
  pscaley[12] = 10000.;
  for (int i=1; i < nBuffCols; i++)
  {

  TYrWriteShapeFile::WriteOneReport(                 wchFoldName // путь к папке
								  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  , nBuffCols // - к-во переменных о корорых накоплена информация в буфере
								  , nBuffRows //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,lenName // максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i  // номер переменной по оси Y
								  ,1. //  масштаб по оси X
								  ,pscaley[i]// масштаб по оси Y
								   ) ;
  }


   delete cmparrPartS;
   delete wcharrFileNames;
   delete parrBuff;
   delete cmparrPartSZv;
   delete cmparrSumS;
   delete cmparrSumSZv;


	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  wchFoldName);
	wcscat(wchAxesFileName, L"AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-5., (((double)nBuffRows) + 5.) * pscaley[0]
	 ,-20. ,20., 5. * pscaley[0]) ;

  TURPointXY  arrPoints1[2];
  arrPoints1[0] =  TURPointXY ( 0., AlfRadUMTrg * pscaley[1]);
  arrPoints1[1] =  TURPointXY ((((double)nBuffRows) + 1.) * pscaley[0], AlfRadUMTrg * pscaley[1]);
  TURPolyLine plnAngTrue(arrPoints1,2) ;

 wcscpy(  wchAxesFileName,  wchFoldName);
 wcscat(wchAxesFileName, L"AngTrue.shp");
 plnAngTrue.WriteSetSHPFiles(wchAxesFileName, &plnAngTrue, 1);
		delete pscaley;
}


// построение графиков СКЗ ошибки измерения угла АРСМ для суммарных и парциальных диграмм
// при различных соотношениях сигнал/шум

void TParAnt::createGraphsSKZ_ARSM(wchar_t *wchFoldName1, TComp cmpKTarg)
{

   	// подготовка директории дл графиков
	wchar_t wchFoldName[300] ={0};
	wcscpy(  wchFoldName,  wchFoldName1);
	wcscat(wchFoldName, L"\\");
  // парциальные диаграммы
  TParAnt ParAnt2Ray(2, mAppert, mLambda
   , mAngSdvig, mNoiseDisp, mAmplFactSig);
  double valGenAngSdvig = transformAngToGeneralizedAng (mAngSdvig ) ;
  ///

  //СУмм диаграммы

  TParAnt ParAnt3Ray(3, mAppert, mLambda
   , mAngSdvig, mNoiseDisp, mAmplFactSig);



   double valStep = 0.0001;
	const int nBuffRows = valGenAngSdvig/ valStep ;
	const int nBuffCols = 4;
	double  *parrBuff = new double [nBuffRows  * nBuffCols] ;
	memset(parrBuff, 0, nBuffRows * nBuffCols * sizeof(double));

	const int lenName = 30;
	wchar_t *wcharrFileNames  = new wchar_t [30 * nBuffCols];
	memset(wcharrFileNames, 0, 30 * nBuffCols * sizeof(wchar_t));

	wcscpy( wcharrFileNames, L"AngRad");
	wcscpy( &wcharrFileNames[30], L"SigRadPartDiagr");
	wcscpy( &wcharrFileNames[60], L"SigRadSumDiagr");

	wcscpy( &wcharrFileNames[90], L"AngGen");





  double valGenAngAntp = 0.;

  TComp cmpKAntp(0.,0.), cmparrPartSZv[2], cmparrPartS[2], cmparrPart3S[3], cmparrPart3SZv[3];
  double valEstPartDiaRadAngTarg =0.;
  TComp cmpPartDiaKTarg(0.,0.);
  double valPsi2 = 0.;
  double valPartDiaDisp = 0.;
  for (int i=0 ; i < nBuffRows; i++)
  {
   double valGenAngTarg = -valGenAngSdvig/ 2. + ((double)i ) * valStep;// -valGenAngSdvig/ 2. + M_PI;//
   double valalfUMTrg = ParAnt2Ray.transformGeneralizedAngToAng (valGenAngTarg ) ;


   ParAnt2Ray.ImitateMeasureArrayPartialDiagrams(  valalfUMTrg,  cmpKTarg, valGenAngAntp, cmpKAntp
		, cmparrPartS, cmparrPartSZv)  ;
   const int NumPartRayPare = 0;
   ParAnt2Ray.estimateARSMethPartDiagr( cmparrPartS, NumPartRayPare, &valEstPartDiaRadAngTarg
		  ,  &cmpPartDiaKTarg , &valPartDiaDisp, &valPsi2);




   parrBuff[ i * nBuffCols] = valalfUMTrg ;

   parrBuff[ i * nBuffCols + 1] =  sqrt(valPartDiaDisp
		 + (valEstPartDiaRadAngTarg - valalfUMTrg)*(valEstPartDiaRadAngTarg - valalfUMTrg))  ;


	ParAnt3Ray.ImitateMeasureArrayPartialDiagrams(  valalfUMTrg,  cmpKTarg, valGenAngAntp, cmpKAntp
		, cmparrPart3S, cmparrPart3SZv)  ;
	ParAnt3Ray.ImitateMeasureArraySumDiagrams(  cmparrPart3S, cmparrPart3SZv, cmparrPartS, cmparrPartSZv);
	ParAnt3Ray.estimateARSMethSumDiagr( cmparrPartS, NumPartRayPare, &valEstPartDiaRadAngTarg
		  ,  &cmpPartDiaKTarg , &valPartDiaDisp, &valPsi2);

   parrBuff[ i * nBuffCols + 2] =  sqrt(valPartDiaDisp
           + (valEstPartDiaRadAngTarg - valalfUMTrg)*(valEstPartDiaRadAngTarg - valalfUMTrg))  ;
   parrBuff[ i * nBuffCols + 3] =  valGenAngTarg   ;

  }

 // double scalex = 100.;
  double *pscaley = new double [nBuffCols] ;
  pscaley[0] = 100.;
  pscaley[1] = 1000.;
  pscaley[2] = 1000.;

  for (int i=1 ; i < nBuffCols-1; i++)
  {

  TYrWriteShapeFile::WriteOneReport(wchFoldName  // путь к папке
								  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  , nBuffCols // - к-во переменных о корорых накоплена информация в буфере
								  , nBuffRows //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,lenName // максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i  // номер переменной по оси Y
								  ,pscaley[0] //  масштаб по оси X
								  ,pscaley[i]// масштаб по оси Y
								   ) ;
  }

  for (int i=1 ; i < nBuffCols-1; i++)
  {

  TYrWriteShapeFile::WriteOneReport(wchFoldName  // путь к папке
								  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  , nBuffCols // - к-во переменных о корорых накоплена информация в буфере
								  , nBuffRows //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,lenName // максимальная длина имени переменной
								  ,3  // номер переменной по оси X
								  ,i  // номер переменной по оси Y
								  ,1. //  масштаб по оси X
								  ,pscaley[i]// масштаб по оси Y
								   ) ;
  }




	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  wchFoldName);
	wcscat(wchAxesFileName, L"AxesArr1.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-mAngSdvig/2. *  pscaley[0], mAngSdvig/2.*  pscaley[0] + 0.25
	 ,0.,1.3, 0.15) ;



	 double valSvigVert = -2.;
	 double valScaleY = 1.;
	wcscpy(  wchAxesFileName,  wchFoldName);
	wcscat(wchAxesFileName, L"PartDiagr0.shp");
	createGraphPartDiagr_from_Gen(wchAxesFileName, -valGenAngSdvig/2., valSvigVert,   valScaleY  )  ;

	wcscpy(  wchAxesFileName,  wchFoldName);
	wcscat(wchAxesFileName, L"PartDiagr1.shp");
	createGraphPartDiagr_from_Gen(wchAxesFileName, valGenAngSdvig/2., valSvigVert,   valScaleY  )  ;

	TURPointXY point(0.,valSvigVert);
    wcscpy(  wchAxesFileName,  wchFoldName);
	wcscat(wchAxesFileName, L"AxesArr2.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-valGenAngSdvig/2., valGenAngSdvig/2.
	 ,-1.5,1.5, 0.15, point) ;
	 delete pscaley;
	 delete  parrBuff;
	 delete wcharrFileNames;



}

// ПОСТРОЕНИЕ ГРАФИКОВ ТОЧЕК И ЛИНИЙ ЦЕЛИ И АНТИПОДА В ПЛОСКОСТИ В И ОБРАЗЕ ПЛОСКОСТИ В
// INPUT
// wchrTableFold
// valTargAngRad, valAntpAngRad
// pntSdvig, arrMtxPer_From_b_To_Image_b
// OUTPUT:
// PointB.shp  - точка В
//  PointImageB.shp - образ точки В
// pntEvelopeAntp.shp  -  точка касания линии антипода к огибающей
// pntEvelopeImageAntp.shp -  образ точка касания линии антипода к огибающей
// pntEvelopeTarg.shp  -  точка касания линии цели к огибающей
// pntEvelopeImageTarg.shp -образ точка касания линии цели к огибающей
// lineTarg.shp    - линия цели в плоскости В
// lineTargImage.shp - образ линии цели
//lineAntp.shp    - линия антипода  в плоскости В
// lineAntpImage.shp - образ линии антипода
void  TParAnt::createGraphsForTabulation1(wchar_t *wchrTableFold,  double valTargAngRad
,  double valAntpAngRad, TURPointXY  pntSdvig, double *arrMtxPer_From_b_To_Image_b)
{
	const double VAlGenAngSdvig= transformAngToGeneralizedAng (mAngSdvig ) ;
	TComp *cmparrPartS = new TComp[ mQuantDiagr];
	memset(cmparrPartS, 0, mQuantDiagr * sizeof(TComp));

	TComp *cmparrPartSZv = new TComp[ mQuantDiagr];
	memset(cmparrPartSZv, 0, mQuantDiagr * sizeof(TComp));
	double alfGenTrg = transformAngToGeneralizedAng (valTargAngRad  ) ;
	double alfGenAntp = transformAngToGeneralizedAng (valAntpAngRad  ) ;

	TComp cmpKTarg(800.,0.);
	TComp cmpKAntp(0.,800.);
	ImitateMeasureArrayPartialDiagrams(valTargAngRad,cmpKTarg
		, valAntpAngRad,cmpKAntp, cmparrPartS,cmparrPartSZv);
	 // вычисленние вектора b
   double arr_b0[2] = {0.}  ;
   bool brez = calcVect_b(&cmparrPartS[2],arr_b0) ;

   TURPointXY pntB(arr_b0[0], arr_b0[1]);

   TURPointXY pntImageBt =  pntB.SdvigTransform(pntSdvig )   ;
   TURPointXY pntImageB =   pntImageBt.fncLinTrasform(arrMtxPer_From_b_To_Image_b);

   TLine lineTarg;
	double val_b0= -1000.,  val_b1 =0.;
	bool br1 = calc_b2(  alfGenTrg , val_b0, &val_b1)  ;
	lineTarg.Points[0].X = val_b0;
	lineTarg.Points[0].Y = val_b1;

	lineTarg.Points[1].X = 1000.;
	br1 = calc_b2( alfGenTrg , lineTarg.Points[1].X, &lineTarg.Points[1].Y)  ;

	 TURPolyLine lineTargImaget = lineTarg.SdvigTransform(pntSdvig )   ;
	 TURPolyLine lineTargImage = lineTargImaget.fncLinTransform(arrMtxPer_From_b_To_Image_b);

	TLine lineAntp;
	lineAntp.Points[0].X = -1000.;
	 br1 = calc_b2(  alfGenAntp , lineAntp.Points[0].X, &lineAntp.Points[0].Y)  ;


	lineAntp.Points[1].X =   1000.;// arr_b0[0];//
	br1 = calc_b2( alfGenAntp , lineAntp.Points[1].X, &lineAntp.Points[1].Y)  ;

	TURPointXY pntEvelopeAntp;
	bool br2= calcEnvelopePoint(  alfGenAntp, pntEvelopeAntp)  ;

	TURPointXY pntEvelopeTarg;
	 br2= calcEnvelopePoint(  alfGenTrg, pntEvelopeTarg)  ;

	 wchar_t wchrPointB [400] = {0};
wcscpy( wchrPointB,wchrTableFold);
 wcscat(wchrPointB, L"\\PointB.shp");
 TURPointXY::WriteSetSHPFiles(wchrPointB, &pntB,  1) ;

  wchar_t wchrPointImageB [400] = {0};
wcscpy(wchrPointImageB,wchrTableFold);
 wcscat(wchrPointImageB, L"\\PointImageB.shp");
 TURPointXY::WriteSetSHPFiles(wchrPointImageB, &pntImageB,  1) ;


	 wchar_t wchrpntEvelopeAntp [400] = {0};
wcscpy( wchrpntEvelopeAntp,wchrTableFold);
 wcscat(wchrpntEvelopeAntp, L"\\pntEvelopeAntp.shp");
 TURPointXY::WriteSetSHPFiles(wchrpntEvelopeAntp, &pntEvelopeAntp,  1) ;


  TURPointXY pntEvelopeImageAntp =  pntEvelopeAntp.LinTransform( pntSdvig,arrMtxPer_From_b_To_Image_b );
	 wchar_t wchrpntEvelopeImageAntp [400] = {0};
wcscpy( wchrpntEvelopeImageAntp,wchrTableFold);
 wcscat(wchrpntEvelopeImageAntp, L"\\pntEvelopeImageAntp.shp");
 TURPointXY::WriteSetSHPFiles(wchrpntEvelopeImageAntp, &pntEvelopeImageAntp,  1) ;

	 wchar_t wchrpntEvelopeTarg [400] = {0};
wcscpy( wchrpntEvelopeTarg,wchrTableFold);
 wcscat(wchrpntEvelopeTarg, L"\\pntEvelopeTarg.shp");
 TURPointXY::WriteSetSHPFiles(wchrpntEvelopeTarg, &pntEvelopeTarg,  1) ;

   TURPointXY pntEvelopeImageTarg =  pntEvelopeTarg.LinTransform( pntSdvig,arrMtxPer_From_b_To_Image_b );
	 wchar_t wchrpntEvelopeImageTarg [400] = {0};
wcscpy( wchrpntEvelopeImageTarg,wchrTableFold);
 wcscat(wchrpntEvelopeImageTarg, L"\\pntEvelopeImageTarg.shp");
 TURPointXY::WriteSetSHPFiles(wchrpntEvelopeImageTarg, &pntEvelopeImageTarg,  1) ;



	 wchar_t wchrlineTarg [400] = {0};
wcscpy( wchrlineTarg,wchrTableFold);
 wcscat(wchrlineTarg, L"\\lineTarg.shp");
 TURPolyLine::WriteSetSHPFiles(wchrlineTarg, &lineTarg,  1) ;

	 wchar_t wchrlineTargImage [400] = {0};
wcscpy( wchrlineTargImage,wchrTableFold);
 wcscat(wchrlineTargImage, L"\\lineTargImage.shp");
 TURPolyLine::WriteSetSHPFiles(wchrlineTargImage, &lineTargImage,  1) ;



  wchar_t wchrlineAntp [400] = {0};
wcscpy( wchrlineAntp,wchrTableFold);
 wcscat(wchrlineAntp, L"\\lineAntp.shp");
 TURPolyLine::WriteSetSHPFiles(wchrlineAntp, &lineAntp,  1) ;

	TURPolyLine lineAntpImaget = lineAntp.SdvigTransform(pntSdvig )   ;
	 TURPolyLine lineAntpImage = lineAntpImaget.fncLinTransform(arrMtxPer_From_b_To_Image_b);
	 wchar_t wchrlineAntpImage [400] = {0};
wcscpy( wchrlineAntpImage,wchrTableFold);
 wcscat(wchrlineAntpImage, L"\\lineAntpImage.shp");
 TURPolyLine::WriteSetSHPFiles(wchrlineAntpImage, &lineAntpImage,  1) ;


 TURPointXY pointCentre(0.,0.);
  TURPolygon plgCircle = TURPolygon::fncCreateCircle( pointCentre
	 ,1., 1000) ;
  wchar_t wchrplgCircle [400] = {0};
wcscpy( wchrplgCircle,wchrTableFold);
 wcscat(wchrplgCircle, L"\\plgCircle.shp");
 TURPolygon::WriteSetSHPFiles(wchrplgCircle, &plgCircle,  1) ;

   delete [] cmparrPartS;
   delete [] cmparrPartSZv;
 }


// ПОСТРОЕНИЕ ГРАФИКОВ ТОЧЕК И ЛИНИЙ ЦЕЛИ И АНТИПОДА В ПЛОСКОСТИ В И ОБРАЗЕ ПЛОСКОСТИ В
// INPUT
// wchrTableFold
// valTargAngRad, valAntpAngRad
// pntSdvig, arrMtxPer_From_b_To_Image_b
// OUTPUT:
// PointB.shp  - точка В
//  PointImageB.shp - образ точки В
// pntEvelopeAntp.shp  -  точка касания линии антипода к огибающей
// pntEvelopeImageAntp.shp -  образ точка касания линии антипода к огибающей
// pntEvelopeTarg.shp  -  точка касания линии цели к огибающей
// pntEvelopeImageTarg.shp -образ точка касания линии цели к огибающей
// lineTarg.shp    - линия цели в плоскости В
// lineTargImage.shp - образ линии цели
//lineAntp.shp    - линия антипода  в плоскости В
// lineAntpImage.shp - образ линии антипода
void  TParAnt::createGraphsForTabulation2(wchar_t *wchrTableFold, TComp *cmparrPartS
	 , TURPointXY  pntSdvig, double *arrMtxPer_From_b_To_Image_b)
{
	const double VAlGenAngSdvig= transformAngToGeneralizedAng (mAngSdvig ) ;




	 // вычисленние вектора b
   double arr_b0[2] = {0.}  ;
   bool brez = calcVect_b(cmparrPartS,arr_b0) ;

   TURPointXY pntB(arr_b0[0], arr_b0[1]);

   TURPointXY pntImageBt =  pntB.SdvigTransform(pntSdvig )   ;
   TURPointXY pntIdealImageB =   pntImageBt.fncLinTrasform(arrMtxPer_From_b_To_Image_b);

	 wchar_t wchrPointB [400] = {0};
wcscpy( wchrPointB,wchrTableFold);
 wcscat(wchrPointB, L"\\PointMeasureB.shp");
 TURPointXY::WriteSetSHPFiles(wchrPointB, &pntB,  1) ;

  wchar_t wchrPointImageB [400] = {0};
wcscpy(wchrPointImageB,wchrTableFold);
 wcscat(wchrPointImageB, L"\\PointIdealImageMeasureB.shp");
 TURPointXY::WriteSetSHPFiles(wchrPointImageB, &pntIdealImageB,  1) ;



   // нахожэдение порядкового нгомера частоты, нумерация с 0
   double valFreq= VAL_C * 100./ mLambda/ 1000000000.;
	int  numFreq = -1;
   for (int i = 0; i < 8; i++)
   {
	 if (fabs(valFreq - constArrMtrxTransf[i * 5]) < 0.01)
	 {
	   numFreq = i;
	   break;
	 }
   }
   if (numFreq == -1)
   {
	return ;
   }
   ///
  // линейное преобразование вектора b
	// сдвиг
  double arrbT[2] = {0.}, arrTransf_b[2] = {0.};
  arrbT[0] = arr_b0[0] + constArrPointsSdvig[ 3 * numFreq +1];
  arrbT[1] = arr_b0[1] + constArrPointsSdvig[ 3 * numFreq +2];
	//матрицей
	double arrMtrxTransf[4] ;
	memcpy(arrMtrxTransf, &constArrMtrxTransf[5 * numFreq +1], 4 * sizeof(double));
  MtrxMultMatrx(arrMtrxTransf, 2, 2, arrbT, 1, arrTransf_b) ;

   TURPointXY pntApprtoximatedImageB(arrTransf_b[0], arrTransf_b[1]);
	wchar_t wchrApprtoximatedImageB [400] = {0};
wcscpy(wchrApprtoximatedImageB,wchrTableFold);
 wcscat(wchrApprtoximatedImageB, L"\\PointApproximatedImageMeasureB.shp");
 TURPointXY::WriteSetSHPFiles(wchrApprtoximatedImageB, &pntApprtoximatedImageB,  1) ;
 ///

 //
	double arrRoots[2] ={0.};
	int iNUmRoots = 0;
 if (mLambda < 3.14)
 { // эллиптический случай
	// нахождение углов точек касания прямых, проходящих через точку arrTransf_b
	//   единичной  окружности
	double valFi0 = atan2(arrTransf_b[1], arrTransf_b[0]);
	double valD = sqrt(arrTransf_b[0] * arrTransf_b[0] + arrTransf_b[1] * arrTransf_b[1]);
	if (valD < 1.)
	{
	  return;
	}
	double gam  =  acos(1. / valD);
	double arrt[2]= {0.};
	arrt[0] = valFi0  + gam;
	arrt[1] = valFi0  - gam;
	fncCorrectArg(arrt);
	fncCorrectArg(&arrt[1]);

	 TURPointXY pntCasat0(cos(arrt[0]), sin(arrt[0]));
	wchar_t wchrCasat0 [400] = {0};
wcscpy( wchrCasat0,wchrTableFold);
 wcscat( wchrCasat0, L"\\PointCasat0.shp");
 TURPointXY::WriteSetSHPFiles( wchrCasat0, &pntCasat0,  1) ;

  TURPointXY pntCasat1(cos(arrt[1]), sin(arrt[1]));
	wchar_t wchrCasat1 [400] = {0};
wcscpy( wchrCasat1,wchrTableFold);
 wcscat( wchrCasat1, L"\\PointCasat1.shp");
 TURPointXY::WriteSetSHPFiles( wchrCasat1, &pntCasat1,  1) ;
 }
}



 void TParAnt::calcApproxPolinomCoeffArray_ForHyperbolicCase (wchar_t * wchrTableFold, const int NUmValuePoints
  , const int NPolinom,double *arrCoeff
  , double *pvalLeftTabDiapMin,double *pvalLeftTabDiapMax, double *pvalLeftRazreshenieMinRad, double *pvalLeftRazreshenieMaxRad
  , double *pvalRightTabDiapMin,double *pvalRightTabDiapMax, double *pvalRightRazreshenieMinRad, double *pvalRightRazreshenieMaxRad
  , double *pvalGenMaxError)
{
  const double VAlGenAngSdvig= transformAngToGeneralizedAng (mAngSdvig ) ;
  TURPointXY pntSdvig;
  double arrMtxPer_From_b_To_Image_b [4] = {0.};
  bool brez = calcTransformationParams( &pntSdvig
   ,arrMtxPer_From_b_To_Image_b);
   // вычисление угла при котором происходит переход касательной с одгной ветви
  //  гиперболы (левой) на другую (правую)
/*  double arrb0[2] = {0.},arrb1[2] = {1000.,1001.}, arrt0[2] = {0.}, arrInvert [4] = {0.};
  InverseMtrx2(arrMtxPer_From_b_To_Image_b, arrInvert);
  MtrxMultMatrx(arrInvert,2, 2, arrb1,1, arrt0) ;
  arrb0[0] =  arrb1[0] -  pntSdvig.X;
  arrb0[1] =  arrb1[1] -  pntSdvig.Y;
  double arrRoots[2]= {0.};
  int  iRez =  findRootsFgr_For_3PartDiagr(arrb0,  arrRoots);
  const double VAlMuCrit = fabs(arrRoots[0]);
  */
  // другой способ
  double valGenOverSwitch = -1.;
  bool brez0 =calcAngOverSwitchForHyperbBranch (NULL,  &valGenOverSwitch) ;
  const double VAlMuCrit = valGenOverSwitch;


  ///
   ///////////////////////////////////////////////////////////////////////////////////////////////
	 ///////////////////////////////////////////////////////////////////////////////////////////////
	 ////// ВЫЧИСЛЕНИЯ ДЛЯ ЗНАЧЕНИЙ МЮ МЕНЬШИХ КРИТИЧЕСКОГО (ЛЕВАЯ ВЕТВЬ ГИПЕРБОЛЫ)/////////////////////////////////////////////////////////////////////////////////////////
	 ///////////////////////////////////////////////////////////////////////////////////////////////
	 ///////////////////////////////////////////////////////////////////////////////////////////////
	 ///////////////////////////////////////////////////////////////////////////////////////////////

  // аппроксимация на отрезке [0; VAlMuCrit]  полиномом степени  NPolinom

  double *arrValueArgs = new double[NUmValuePoints];
  const double VAlDiapMax = VAlMuCrit-0.01; //M_PI/3.-0.25 +  VAlGenAngSdvig ;
  const double VAlDiapMin = 0.;
  double *arrFunc = new double[NUmValuePoints];
  int iNumValuePointsReal=0;
 // CreateValueArraysHyperbolicType(VAlGenAngSdvig ,mLambda,  VAlDiap
  // ,arrValueArgs,arrFunc, NUmValuePoints , &iNumValuePointsReal);
  CreateValueArraysHyperbolicType(  VAlDiapMin,  VAlDiapMax
   ,arrValueArgs,arrFunc,  NUmValuePoints,  &iNumValuePointsReal);
	TURPolyLine plnValue( 1, iNumValuePointsReal) ;
   for (int i = 0; i < iNumValuePointsReal; i++)
   {
	 plnValue.Points[i].X = arrValueArgs[i];
	 plnValue.Points[i].Y =  arrFunc[i];
   }

	if(wchrTableFold != NULL)
	 {
	 // создание shape  файла линии значений
	wchar_t wchrValuePoliLine[400] = {0};
	wcscpy( wchrValuePoliLine, wchrTableFold);
	wcscat(wchrValuePoliLine, L"\\ValuesLine.shp");
	TURPolyLine::WriteSetSHPFiles(wchrValuePoliLine,&plnValue,  1) ;
	}
	///
		*pvalLeftTabDiapMin = 0.;
		*pvalLeftTabDiapMax = arrValueArgs[iNumValuePointsReal -1];
		*pvalLeftRazreshenieMinRad = 0.;
		*pvalLeftRazreshenieMaxRad =  transformGeneralizedAngToAng(VAlMuCrit );

	 const double STepTabTemp =  0.1;


	const int NUmPointsRegularTab = fabs(arrValueArgs[iNumValuePointsReal -1] - arrValueArgs[0]) / STepTabTemp +1;// 2.5/ STepTabTemp +1;
	const double STepTab = (arrValueArgs[iNumValuePointsReal -1] - arrValueArgs[0])/ ((double)NUmPointsRegularTab -1.);

	double * arrRegularTabFunc = new double[NUmPointsRegularTab];



	CreateRegTab(arrValueArgs, arrFunc ,iNumValuePointsReal
		 ,arrRegularTabFunc , NUmPointsRegularTab, STepTab );

	TURPolyLine plnRegTargTab( 1, NUmPointsRegularTab) ;
	double *arrArg1 = new double[NUmPointsRegularTab];
   for (int i = 0; i < NUmPointsRegularTab; i++)
   {
	 arrArg1[i]  = arrValueArgs [0] + ((double)i)* STepTab;
	 plnRegTargTab.Points[i].X = arrValueArgs [0] +((double)i)* STepTab;

	 plnRegTargTab.Points[i].Y =  arrRegularTabFunc[i];
   }

   if(wchrTableFold != NULL)
	 {
   // полилилиния регулярной таблицы
	wchar_t wchrRegTargTablePoliLine[400] = {0};
	wcscpy( wchrRegTargTablePoliLine, wchrTableFold);
	wcscat(wchrRegTargTablePoliLine, L"\\RegTargTablePoliLine.shp");
	TURPolyLine::WriteSetSHPFiles(wchrRegTargTablePoliLine,&plnRegTargTab,  1) ;
	}
	 ///

	///

	  // полилиния полинома
	double valMaxDelta = 0.;
	CalcPolinomCoefForOddFunc(2 *NPolinom +1, arrArg1, arrRegularTabFunc , arrCoeff, NUmPointsRegularTab, &valMaxDelta);

     if(wchrTableFold != NULL)
	 {
	wchar_t wchrPolinomPoliLine1[400] = {0};
	wcscpy( wchrPolinomPoliLine1,wchrTableFold);
	wcscat(wchrPolinomPoliLine1, L"\\Polinom1.shp");
	TURPolyLine::createOddPolinomGraph(wchrPolinomPoliLine1, 2 *NPolinom +1, arrCoeff, arrArg1[NUmPointsRegularTab-1] ) ;
	///
	}
// вычисление максимальной ошибки
	TURPointXY *pPntArrDel =  new TURPointXY[plnValue.NumPoints];
	int iLenArrDel = 0;

	 for (int i =0; i < plnValue.NumPoints; i++)
	 {

	   pPntArrDel[iLenArrDel].Y = plnValue.Points[i].Y
			-fncOddPolinom(arrCoeff, 2 *NPolinom +1, plnValue.Points[i].X);
	   pPntArrDel[iLenArrDel].X = plnValue.Points[i].X;
	   iLenArrDel++;
	 }
	 TURPolyLine plnErr1(pPntArrDel,iLenArrDel);
	 plnErr1.calcBoundBox();
	 *pvalGenMaxError = fabs(plnErr1.Box[1]);
	 if (fabs(plnErr1.Box[3]) > fabs(plnErr1.Box[1]) )
	 {
	   *pvalGenMaxError = fabs(plnErr1.Box[3]);
	 }

	 if(wchrTableFold != NULL)
	 {

	// постронеие графика ошибки аппроксимации полиномом
	 wchar_t wchrDeltaApprox[400] = {0};
	wcscpy( wchrDeltaApprox, wchrTableFold);
	wcscat(wchrDeltaApprox, L"\\DelMilRad.shp");
	TURPolyLine::WriteSetSHPFiles(wchrDeltaApprox,&plnErr1,  1) ;
   ///////////////////////////////////////////////////////////////////////////////////////////////
	 }
	delete arrArg1;
	delete [] pPntArrDel ;
	delete arrRegularTabFunc;
	 ///////////////////////////////////////////////////////////////////////////////////////////////
	 ///////////////////////////////////////////////////////////////////////////////////////////////
	 ////// ВЫЧИСЛЕНИЯ ДЛЯ ЗНАЧЕНИЙ МЮ БОЛЬШИХ КРИТИЧЕСКОГО (ПРАВАЯ ВЕТВЬ ГИПЕРБОЛЫ)/////////////////////////////////////////////////////////////////////////////////////////
	 ///////////////////////////////////////////////////////////////////////////////////////////////
	 ///////////////////////////////////////////////////////////////////////////////////////////////
	 ///////////////////////////////////////////////////////////////////////////////////////////////
	  //


  double *arrValueArgs1 = new double[NUmValuePoints];
  const double VAlDiapMax1 = M_PI/3.-0.3 +  VAlGenAngSdvig ;
  const double VAlDiapMin1 = VAlMuCrit +0.01;


		*pvalRightRazreshenieMinRad = transformGeneralizedAngToAng(VAlMuCrit );
		*pvalRightRazreshenieMaxRad = transformGeneralizedAngToAng(VAlDiapMax1 );

 // double *arrFunc = new double[NUmValuePoints];
  int iNumValuePointsReal1=0;

  CreateValueArraysHyperbolicType(  VAlDiapMin1,  VAlDiapMax1
   ,arrValueArgs,arrFunc,  NUmValuePoints,  &iNumValuePointsReal1);
	TURPolyLine plnValue11( 1, iNumValuePointsReal1) ;
   for (int i = 0; i < iNumValuePointsReal1; i++)
   {
	 plnValue11.Points[i].X = arrValueArgs[i];
	 plnValue11.Points[i].Y =  arrFunc[i];
   }


	///
		double valTabDiap1 = arrValueArgs[iNumValuePointsReal1 -1];

	 const double STepTabTemp1 =  0.15;
	const int NUmPointsRegularTab1 = fabs(arrValueArgs[iNumValuePointsReal1 -1] - arrValueArgs[0]) / STepTabTemp1 +1;// 2.5/ STepTabTemp +1;
	const double STepTab1 = (arrValueArgs[iNumValuePointsReal1 -1] - arrValueArgs[0] )/ ((double)NUmPointsRegularTab1 -1.);

	double * arrRegularTabFunc1 = new double[NUmPointsRegularTab1];



	CreateRegTab(arrValueArgs, arrFunc ,iNumValuePointsReal1
		 ,arrRegularTabFunc1 , NUmPointsRegularTab1, STepTab1 );

		*pvalRightTabDiapMin = arrValueArgs [0];
		*pvalRightTabDiapMax = arrValueArgs[iNumValuePointsReal1 -1];

 // полилилиния регулярной таблицы
   TURPolyLine plnRegTargTab10( 1, NUmPointsRegularTab1) ;

   for (int i = 0; i < NUmPointsRegularTab1; i++)
   {

	 plnRegTargTab10.Points[i].X = arrValueArgs [0] +((double)i)* STepTab1;


	 plnRegTargTab10.Points[i].Y =  arrRegularTabFunc1[i] ;


   }




	TURPolyLine plnRegTargTab1( 1, NUmPointsRegularTab1) ;
	double *arrArg11 = new double[NUmPointsRegularTab1];
   for (int i = 0; i < NUmPointsRegularTab1; i++)
   {
	 arrArg11[i]  = arrValueArgs [0] + ((double)i)* STepTab1;
	 plnRegTargTab1.Points[i].X = arrValueArgs [0] +((double)i)* STepTab1;


   //	 plnRegTargTab1.Points[i].Y =  arrRegularTabFunc1[i] *plnRegTargTab1.Points[i].X;
	 arrRegularTabFunc1[i] = arrRegularTabFunc1[i] *plnRegTargTab1.Points[i].X*plnRegTargTab1.Points[i].X;
	 plnRegTargTab1.Points[i].Y =  arrRegularTabFunc1[i];
   }



	  // полилиния полинома
	double valMaxDelta1 = 0.;

	//CalcPolinomCoefForOddFuncMinusZ(NPolinom, arrArg11, arrRegularTabFunc1 , &arrCoeff[NPolinom/2 +1], NUmPointsRegularTab1, &valMaxDelta1);
  //	CalcPolinomCoefForOddFunc(NPolinom, arrArg11, arrRegularTabFunc1 , arrCoeff, NUmPointsRegularTab, &valMaxDelta);

 // int nPolinom = 5;
	  CalcPolinomCoefApprox( NPolinom, arrArg11, arrRegularTabFunc1 ,NUmPointsRegularTab1
   ,&arrCoeff[NPolinom +1],  &valMaxDelta1);
// вычисление максимальной ошибки
	TURPointXY *pPntArrDel1 =  new TURPointXY[plnValue11.NumPoints];
	int iLenArrDel1 = 0;

	 for (int i =0; i < plnValue11.NumPoints; i++)
	 {

	   pPntArrDel1[iLenArrDel1].Y = plnValue11.Points[i].Y
			-fncPolinom(&arrCoeff[NPolinom +1], NPolinom +1, plnValue11.Points[i].X)/ plnValue11.Points[i].X/ plnValue11.Points[i].X;
	   pPntArrDel1[iLenArrDel1].X = plnValue11.Points[i].X;
	   iLenArrDel1++;
	 }
	 TURPolyLine plnValue111(pPntArrDel1,iLenArrDel1);
	 plnValue111.calcBoundBox();
	 double valGenMaxError1 = fabs(plnValue111.Box[1]);
	 if (fabs(plnValue111.Box[3]) > fabs(plnValue111.Box[1]) )
	 {
	   valGenMaxError1 = fabs(plnValue111.Box[3]);
	 }
	 ///

	 if ( valGenMaxError1 >(*pvalGenMaxError)  )
	 {
		(*pvalGenMaxError) = valGenMaxError1  ;
	 }

	 double valXTemp0 = plnValue11.Points[plnValue11.NumPoints-1].X ;
	  double valTemp1= fncPolinom(&arrCoeff[NPolinom +1], NPolinom +1, valXTemp0)/valXTemp0/ valXTemp0;

	//  *pvalRazreshenieRad = fabs(transformGeneralizedAngToAng (valTemp1)) ;
	 // ВЫВОД ГРАФИКИ
	 if(wchrTableFold != NULL)
	 {
   	 // создание shape  файла линии значений
	wchar_t wchrValuePoliLine1[400] = {0};
	wcscpy( wchrValuePoliLine1, wchrTableFold);
	wcscat(wchrValuePoliLine1, L"\\ValuesLine1.shp");
	TURPolyLine::WriteSetSHPFiles(wchrValuePoliLine1,&plnValue11,  1) ;

	///

	wchar_t wchrRegTargTablePoliLine11[400] = {0};
	wcscpy( wchrRegTargTablePoliLine11, wchrTableFold);
	wcscat(wchrRegTargTablePoliLine11, L"\\RegTargTablePoliLine10.shp");
	TURPolyLine::WriteSetSHPFiles(wchrRegTargTablePoliLine11,&plnRegTargTab10,  1) ;
///
    	// полилилиния регулярной таблицы
	wchar_t wchrRegTargTablePoliLine12[400] = {0};
	wcscpy( wchrRegTargTablePoliLine12, wchrTableFold);
	wcscat(wchrRegTargTablePoliLine12, L"\\RegTargTablePoliLine12.shp");
	TURPolyLine::WriteSetSHPFiles(wchrRegTargTablePoliLine12,&plnRegTargTab1,  1) ;

	///

	wchar_t wchrPolinomPoliLine11[400] = {0};
	wcscpy( wchrPolinomPoliLine11,wchrTableFold);
	wcscat(wchrPolinomPoliLine11, L"\\Polinom11.shp");
	TURPolyLine::createPolinomGraph(wchrPolinomPoliLine11, NPolinom
	  , &arrCoeff[NPolinom +1], arrArg11[NUmPointsRegularTab1-1], arrArg11[0] ) ;
	///
	// постронеие графика ошибки аппроксимации полиномом
	 wchar_t wchrDeltaApprox1[400] = {0};
	wcscpy( wchrDeltaApprox1, wchrTableFold);
	wcscat(wchrDeltaApprox1, L"\\DelMilRad1.shp");
	TURPolyLine::WriteSetSHPFiles(wchrDeltaApprox1,&plnValue111,  1) ;
	 }
	///




	delete arrValueArgs ;
	delete arrFunc ;
	delete arrArg11;
	delete [] pPntArrDel1 ;
	delete arrRegularTabFunc1;


}

bool TParAnt::calcAngOverSwitchForHyperbBranch (wchar_t * wchrTableFold, double *pvalGenOverSwitch)
{
	  const double VAlGenAngSdvig= transformAngToGeneralizedAng (mAngSdvig ) ;
	 bool brez = false;


	 double valStep = 0.0005;

	 double valBegin = 0.01;
	 double valEnd = 50.;
	 int iNc = (valEnd - valBegin)/ valStep;
	 double valGenAng = valBegin;
	 double valF = fncFTemp1(valGenAng, VAlGenAngSdvig);
	 valF = (valF < 0.)?-1.:1.;
	for (int i =0 ; i < iNc; i++)
	{
	 double valGenAng1 = valBegin + ((double)i) * valStep;
	 double valF1 = fncFTemp1(valGenAng1, VAlGenAngSdvig);
	 if (valF1 *valF <= 0.)
	 {
	   *pvalGenOverSwitch = valGenAng1 ;
	   brez = true;
	   break;
	 }

	 }

   return brez ;
}

double fncFTemp1(double valGenAng,double  VAlGenAngSdvig)
{
	return fncDiagrSinx_div_x(valGenAng + VAlGenAngSdvig) * fncDerivDiagrSinx_div_x(valGenAng - VAlGenAngSdvig)
	 - fncDiagrSinx_div_x(valGenAng - VAlGenAngSdvig) * fncDerivDiagrSinx_div_x(valGenAng + VAlGenAngSdvig) ;
}

bool  TParAnt::CreateRegTab(double *arrArg, double *arrFunc, const int LEnArr
  ,double *arrRegularTabFunc , const int LEnRegTab, const double STepTab)
{

 for (int i = 0; i < LEnRegTab; i++)
 {
	double arg = arrArg[0] + ((double)i) * STepTab;
	int numCur = 0;
	for (int j = 0; j < LEnArr; j++)
	{
	 // if ((arg >= arrArg[j])&& (arg <= arrArg[j +1]) )
	  if ((arg - arrArg[j])* (arg -arrArg[j +1])<= 0. )
	  {
	   numCur = j;
	   break;
	  }
	}
	TURPointXY  pnt1(arrArg[numCur],arrFunc[numCur]);
	TURPointXY  pnt2(arrArg[numCur +1],arrFunc[numCur +1]);
	TLine LineTemp(  pnt1,   pnt2);
	if (!LineTemp.calcY(arg, &arrRegularTabFunc[i]) )
	{
		break;
		return false;
    }
 }
 return true;
}

// создание массива точек для гиперболического типа
// INPUT:
//VAlGenAngSdvig - угол сдвига обобщенный
// VAlLambda  - длина волны
// VAlDiap  - диапазон изменения угла (цели) в обобщенных координатах
// NUmPoints0  - к-во точек
// OUTPUT:
// arrArgs[ NUmPoints0] - массив значений параметра точки касания к гиперболе
// arrFunc [ NUmPoints0]- массив значений исходных углов целей на регулярной сетке
//                       от -VAlDiap до VAlDiap
void TParAnt::CreateValueArraysHyperbolicType( const double VAlDiapMin, const double VAlDiapMax
   ,double *arrArgs,double *arrFunc, const int NUmPoints0,  int *piNumPointsReal )
{
   const double VAlGenAngSdvig= transformAngToGeneralizedAng (mAngSdvig ) ;
  TURPointXY pntSdvig ;
  double arrMtxPer_From_b_To_Image_b[4] ={0.};
  calcTransformationParams( &pntSdvig, arrMtxPer_From_b_To_Image_b);
  const double STep =  (VAlDiapMax - VAlDiapMin)/ ((double)NUmPoints0 -1.);
  TURPointXY pntEnvel;
  *piNumPointsReal = 0;
  for (int i = 0; i < NUmPoints0; i++)
  {
	double valMu = VAlDiapMin + ((double)i) *  STep;
	calcEnvelopePoint(valMu, pntEnvel);

	TURPointXY pntPoject =  pntEnvel.LinTransform( pntSdvig,arrMtxPer_From_b_To_Image_b );

   //	if ((pntPoject.X* pntPoject.X -1.)< 0.0000001)
   //	{
   //	  continue;
   //	}
	double tamp = 0. ;
	if (( pntPoject.X* pntPoject.X -1.) > 0.)
	{
	  tamp =  sqrt( pntPoject.X* pntPoject.X -1.) ;
	}
	pntPoject.Y = pntPoject.Y / fabs(pntPoject.Y) * tamp;
	double temp = pntPoject.Y/pntPoject.X;
   //	if (fabs(temp) > 0.99)
  //	{
  //	  continue;
  //	}
	arrArgs[(*piNumPointsReal)] = log((1. + temp)/ (1. - temp))/ 2.;//temp;//!
	arrFunc[(*piNumPointsReal)] = valMu;

	(*piNumPointsReal)++;
  }
}


bool  TParAnt::calcTransformationParams_old( TURPointXY *pntSdvig, double *arrMtxPer_From_b_To_Image_b)
{
	const double VAlGenAngSdvig= transformAngToGeneralizedAng (mAngSdvig ) ;


	double valBSingular = 0., valMuSingular =0.;
	if( VAlGenAngSdvig > VAlCritAngSdvig)
	{
	findEnvelopeSigularPoint (  &valMuSingular, &valBSingular) ;
	}
	else
	{
	 valBSingular = -0.5;
    }
	 double valBRegular = 1./ fncDiagrSinx_div_x(VAlGenAngSdvig) /2.;
	TURPointXY pntVertex1(valBRegular, valBRegular);
	TURPointXY pntVertex2(valBSingular, valBSingular);
	*pntSdvig = TURPointXY::ParamPoint(pntVertex1, pntVertex2,0.5);
	(*pntSdvig).X = -(*pntSdvig).X;
	(*pntSdvig).Y = -(*pntSdvig).Y;
	double alf = M_PI/4. + M_PI;
	double arrMtxRotate[4] = {0.};
	arrMtxRotate [0] = cos(alf);
	arrMtxRotate [1] = sin (alf);
	arrMtxRotate [2] = -sin (alf);
	arrMtxRotate [3] = cos(alf);

	double val_k = fncCurvation( ) ;
	double vala = fabs(( valBRegular -  valBSingular)) / sqrt(2.);
	double valb = sqrt(vala /fabs(val_k ));
	if (VAlGenAngSdvig > VAlCritAngSdvig)
	{
		vala -= 0.1675;
		valb += 0.13 ;
	}
	else
	{
	  valb += -0.544;
	}

	double  arrMtxCompress[4] ={0.};
	arrMtxCompress[0] = 1./vala;
	arrMtxCompress[3] = 1./ valb;
	MtrxMultMatrx(arrMtxCompress,2 ,2 , arrMtxRotate,2, arrMtxPer_From_b_To_Image_b)  ;

}


bool  TParAnt::calcTransformationParams( TURPointXY *pntSdvig, double *arrMtxPer_From_b_To_Image_b)
{
	const double VAlGenAngSdvig= transformAngToGeneralizedAng (mAngSdvig ) ;
	if( VAlGenAngSdvig > VAlCritAngSdvig)
	{
		calcTransformationParamsEllipticCase( pntSdvig, arrMtxPer_From_b_To_Image_b) ;
	}
	else
	{
	  calcTransformationParamsHyperbolicCase( pntSdvig, arrMtxPer_From_b_To_Image_b) ;
	}
	return true;

  /*	double valBSingular = 0., valMuSingular =0.;
	if( VAlGenAngSdvig > VAlCritAngSdvig)
	{
	findEnvelopeSigularPoint (  &valMuSingular, &valBSingular) ;
	}
	else
	{
	 valBSingular = -0.5;
    }
	 double valBRegular = 1./ fncDiagrSinx_div_x(VAlGenAngSdvig) /2.;
	TURPointXY pntVertex1(valBRegular, valBRegular);
	TURPointXY pntVertex2(valBSingular, valBSingular);
	*pntSdvig = TURPointXY::ParamPoint(pntVertex1, pntVertex2,0.5);
	(*pntSdvig).X = -(*pntSdvig).X;
	(*pntSdvig).Y = -(*pntSdvig).Y;
	double alf = M_PI/4. + M_PI;
	double arrMtxRotate[4] = {0.};
	arrMtxRotate [0] = cos(alf);
	arrMtxRotate [1] = sin (alf);
	arrMtxRotate [2] = -sin (alf);
	arrMtxRotate [3] = cos(alf);

	double val_k = fncCurvation( ) ;
	double vala = fabs(( valBRegular -  valBSingular)) / sqrt(2.);
	double valb = sqrt(vala /fabs(val_k ));
	if (VAlGenAngSdvig > VAlCritAngSdvig)
	{
		vala -= 0.1675;
		valb += 0.13 ;
	}
	else
	{
	  valb += -0.544;
	}

	double  arrMtxCompress[4] ={0.};
	arrMtxCompress[0] = 1./vala;
	arrMtxCompress[3] = 1./ valb;
	MtrxMultMatrx(arrMtxCompress,2 ,2 , arrMtxRotate,2, arrMtxPer_From_b_To_Image_b)  ;
  */
}
bool  TParAnt::calcTransformationParamsEllipticCase( TURPointXY *pntSdvig, double *arrMtxPer_From_b_To_Image_b)
{
   const double VAlGenAngSdvig= transformAngToGeneralizedAng (mAngSdvig ) ;


	double valBSingular = 0., valMuSingular =0.;
	 findEnvelopeSigularPoint (  &valMuSingular, &valBSingular) ;
	 double valBRegular = 1./ fncDiagrSinx_div_x(VAlGenAngSdvig) /2.;
	TURPointXY pntVertex1(valBRegular, valBRegular);
	TURPointXY pntVertex2(valBSingular, valBSingular);
	*pntSdvig = TURPointXY::ParamPoint(pntVertex1, pntVertex2,0.5);
	(*pntSdvig).X = -(*pntSdvig).X;
	(*pntSdvig).Y = -(*pntSdvig).Y;
	double alf = M_PI/4. + M_PI;
	double arrMtxRotate[4] = {0.};
	arrMtxRotate [0] = cos(alf);
	arrMtxRotate [1] = sin (alf);
	arrMtxRotate [2] = -sin (alf);
	arrMtxRotate [3] = cos(alf);

	double val_k = fncCurvation( ) ;
	double vala = fabs(( valBRegular -  valBSingular)) / sqrt(2.);
	double valb = sqrt(vala /fabs(val_k ));

	vala -= 0.1675;
	valb += 0.13 ;
	double  arrMtxCompress[4] ={0.};
	arrMtxCompress[0] = 1./vala;
	arrMtxCompress[3] = 1./ valb;
	MtrxMultMatrx(arrMtxCompress,2 ,2 , arrMtxRotate,2, arrMtxPer_From_b_To_Image_b)  ;
	return true;
}
bool  TParAnt::calcTransformationParamsHyperbolicCase( TURPointXY *pntSdvig, double *arrMtxPer_From_b_To_Image_b)
{
	const double VAlGenAngSdvig= transformAngToGeneralizedAng (mAngSdvig ) ;


	double valBSingular = 0., valMuSingular =0.;
	valBSingular = -0.5;

	 double valBRegular = 1./ fncDiagrSinx_div_x(VAlGenAngSdvig) /2.;
	TURPointXY pntVertex1(valBRegular, valBRegular);
	TURPointXY pntVertex2(valBSingular, valBSingular);
	*pntSdvig = TURPointXY::ParamPoint(pntVertex1, pntVertex2,0.5);
	(*pntSdvig).X = -(*pntSdvig).X;
	(*pntSdvig).Y = -(*pntSdvig).Y;
	double alf = M_PI/4. + M_PI;
	double arrMtxRotate[4] = {0.};
	arrMtxRotate [0] = cos(alf);
	arrMtxRotate [1] = sin (alf);
	arrMtxRotate [2] = -sin (alf);
	arrMtxRotate [3] = cos(alf);

  //	double val_k = fncCurvation( ) ;
	double vala = fabs(( valBRegular -  valBSingular)) / sqrt(2.);
  //	double valb = sqrt(vala /fabs(val_k ));

	double  arrMtxCompress1[4] ={0.};
	arrMtxCompress1[0] = 1./vala;
	arrMtxCompress1[3] = 1.;
	double arrMtxPerTemp[4] ={0.};
	MtrxMultMatrx(arrMtxCompress1,2 ,2 , arrMtxRotate,2, arrMtxPerTemp)  ;


   const double VAlDiap = M_PI/3.-0.25 +  VAlGenAngSdvig ;
   const int NUmPoints0 = 10000;
   int numPointsActual = -1;
   TURPointXY *arrPoints = new TURPointXY[NUmPoints0];
   //TURPointXY *arrCentredPoints = new TURPointXY[NUmPoints0];


 calcEvelopePointsArray( VAlDiap,  NUmPoints0, arrPoints, &numPointsActual );
 for (int i =0; i < numPointsActual; i++)
 {
 arrPoints[i] = arrPoints[i].LinTransform( *pntSdvig ,arrMtxPerTemp );
 }
 double maxx = 0.;
 int iarg = -1;
 for (int i =0; i < numPointsActual; i++)
 {
	if( fabs(arrPoints[i].Y) >  maxx)
	{
	 maxx = fabs(arrPoints[i].X);
	 iarg = i;

	}
 }

	arrMtxCompress1[3] = fabs(arrPoints[iarg].X/ arrPoints[iarg].Y);
	MtrxMultMatrx(arrMtxCompress1,2 ,2 , arrMtxRotate,2, arrMtxPer_From_b_To_Image_b)  ;





	delete []arrPoints ;

  }
/*
void TParAnt::calcApproxPolinomCoeffArray_ForEllipticCase (wchar_t * wchrTableFold, const int NUmValuePoints
  , const int NPolinom,double *arrCoeff , double *pvalTabDiap, double *pvalGenMaxError, double *pvalRazreshenie)
{
   const double VAlGenAngSdvig= transformAngToGeneralizedAng (mAngSdvig ) ;
  double *arrValueArgs = new double[NUmValuePoints];
  const double VAlDiap = M_PI/3.-0.9 +  VAlGenAngSdvig ;
  double *arrFunc = new double[NUmValuePoints];
  CreateValueArraysEllipsType__(VAlDiap ,arrValueArgs,arrFunc, NUmValuePoints );

	///

	// флиппование
   fncFlipDblArr(arrFunc , NUmValuePoints);
   fncFlipDblArr(arrValueArgs , NUmValuePoints);
   ///

   TURPolyLine plnValue( 1, NUmValuePoints) ;
   for (int i = 0; i < NUmValuePoints; i++)
   {
	 plnValue.Points[i].X = arrValueArgs[i];
	 plnValue.Points[i].Y =  arrFunc[i];
   }


	// создание таблицы с регулярным шагом
	// аргумент меняется от 0 до ARgTabMax с шагом STepTab

		//*pvalTabDiap = 0.975 *fabs(arrValueArgs[ NUmValuePoints - 1]) ;
		*pvalTabDiap =fabs(arrValueArgs[ NUmValuePoints - 1]) ;

   //	 const double STepTabTemp =  0.03;
	  const double STepTabTemp =  0.1;
	const int NUmPointsRegularTab = (*pvalTabDiap) / STepTabTemp +1;// 2.5/ STepTabTemp +1;
	const double STepTab =  (*pvalTabDiap)/ ((double)NUmPointsRegularTab -1.);

	double * arrRegularTabFunc = new double[NUmPointsRegularTab];



	CreateRegularTab(&arrValueArgs[(NUmValuePoints - 1) /2 ], &arrFunc[(NUmValuePoints - 1) /2 ]
		, (NUmValuePoints - 1) /2 +1 ,arrRegularTabFunc , NUmPointsRegularTab, STepTab );

	TURPolyLine plnRegTargTab( 1, NUmPointsRegularTab) ;
	double *arrArg1 = new double[NUmPointsRegularTab];
   for (int i = 0; i < NUmPointsRegularTab; i++)
   {
	 arrArg1[i]  = ((double)i)* STepTab;
	 plnRegTargTab.Points[i].X = ((double)i)* STepTab;

	 plnRegTargTab.Points[i].Y =  arrRegularTabFunc[i];
   }

	///


	  // полилиния полинома
	double valMaxDelta1 = 0.;
	CalcPolinomCoefForOddFunc(NPolinom, arrArg1, arrRegularTabFunc , arrCoeff, NUmPointsRegularTab, &valMaxDelta1);



	//double kTang0 = (plnTargTab.Points[NUmPointsTab - 1].Y - plnTargTab.Points[NUmPointsTab - 2 ].Y )
	 //	 / (plnTargTab.Points[NUmPointsTab- 1].X - plnTargTab.Points[NUmPointsTab- 2 ].X );

	// вычисление максимальной ошибки
	TURPointXY *pPntArrDel =  new TURPointXY[plnValue.NumPoints];
	int iLenArrDel = 0;

	 for (int i =0; i < plnValue.NumPoints; i++)
	 {
	   if (plnValue.Points[i].X > (*pvalTabDiap))
	   {
		break;
	   }
	   if (plnValue.Points[i].X < -(*pvalTabDiap))
	   {
		continue;
	   }
	 //  pPntArrDel[iLenArrDel].Y = plnValue.Points[i].Y
	  //		-fncOddPolinom(arrCoeff, NPolinom / 2 +1, plnValue.Points[i].X);
	   pPntArrDel[iLenArrDel].Y = tan(tan(plnValue.Points[i].Y) )
			-tan(tan(fncOddPolinom(arrCoeff, NPolinom / 2 +1, plnValue.Points[i].X)));
	   pPntArrDel[iLenArrDel].X = plnValue.Points[i].X;
	   iLenArrDel++;
	 }
	 TURPolyLine plnValue1(pPntArrDel,iLenArrDel);
	 plnValue1.calcBoundBox();
	 *pvalGenMaxError = fabs(plnValue1.Box[1]);
	 if (fabs(plnValue1.Box[3]) > fabs(plnValue1.Box[1]) )
	 {
	   *pvalGenMaxError = fabs(plnValue1.Box[3]);
	 }

   //	double valTemp0= tan(tan(arrRegularTabFunc[NUmPointsRegularTab-1]));
	double valTemp0= tan(tan(arrRegularTabFunc[NUmPointsRegularTab-1]));
	*pvalRazreshenie = transformGeneralizedAngToAng (valTemp0  ) ;
	  //
	 // ВЫВОД ГРАФИКИ
	 if(wchrTableFold != NULL)
	 {
	 // создание shape  файла линии значений
	wchar_t wchrTargTablePoliLine[400] = {0};
	wcscpy( wchrTargTablePoliLine, wchrTableFold);
	wcscat(wchrTargTablePoliLine, L"\\ValuesLine.shp");
	TURPolyLine::WriteSetSHPFiles(wchrTargTablePoliLine,&plnValue,  1) ;
	///

	// полилилиния регулярной таблицы
	wchar_t wchrRegTargTablePoliLine[400] = {0};
	wcscpy( wchrRegTargTablePoliLine, wchrTableFold);
	wcscat(wchrRegTargTablePoliLine, L"\\RegTargTablePoliLine.shp");
	TURPolyLine::WriteSetSHPFiles(wchrRegTargTablePoliLine,&plnRegTargTab,  1) ;
	 ///


	wchar_t wchrPolinomPoliLine1[400] = {0};
	wcscpy( wchrPolinomPoliLine1,wchrTableFold);
	wcscat(wchrPolinomPoliLine1, L"\\Polinom1.shp");
	TURPolyLine::createOddPolinomGraph(wchrPolinomPoliLine1, NPolinom, arrCoeff, arrArg1[NUmPointsRegularTab-1] ) ;
	///

	// постронеие графика ошибки аппроксимации полиномом
	 wchar_t wchrDeltaApprox[400] = {0};
	wcscpy( wchrDeltaApprox, wchrTableFold);
	wcscat(wchrDeltaApprox, L"\\DelMilRad.shp");
	TURPolyLine::WriteSetSHPFiles(wchrDeltaApprox,&plnValue1,  1) ;
	 }
	///
	delete arrValueArgs ;
	delete arrFunc ;
	delete arrArg1;
	delete [] pPntArrDel ;
	delete arrRegularTabFunc;
}
*/

void TParAnt::calcApproxPolinomCoeffArray_ForEllipticCase (wchar_t * wchrTableFold, const int NUmValuePoints
  , const int NPolinom,double *arrCoeff , double *pvalTabDiap, double *pvalGenMaxError
	, double *pvalRazreshenie, double *valSlop)
{
   const double VAlGenAngSdvig= transformAngToGeneralizedAng (mAngSdvig ) ;
  double *arrValueArgs = new double[NUmValuePoints];
  const double VAlDiap = M_PI/3.-0.55 +  VAlGenAngSdvig ;
  double *arrFunc0 = new double[NUmValuePoints];
  CreateValueArraysEllipsType(VAlDiap ,arrValueArgs,arrFunc0, NUmValuePoints );

	///

	// флиппование
   fncFlipDblArr(arrFunc0 , NUmValuePoints);
   fncFlipDblArr(arrValueArgs , NUmValuePoints);
   ///

   TURPolyLine plnValue0( 1, NUmValuePoints) ;
   for (int i = 0; i < NUmValuePoints; i++)
   {
	 plnValue0.Points[i].X = arrValueArgs[i];
	 plnValue0.Points[i].Y =  arrFunc0[i];
   }

   TLine  lineTemp(plnValue0.Points[NUmValuePoints/2], plnValue0.Points[NUmValuePoints/2 +1]) ;

   lineTemp.calcTang(valSlop);

   TURPolyLine plnValue =  plnValue0;
   double *arrFunc = new double[NUmValuePoints];
   for (int i = 0; i < NUmValuePoints; i++)
   {

	 plnValue.Points[i].Y = atan( plnValue.Points[i].Y / (*valSlop) );
	 arrFunc [i] = atan(arrFunc0[i] / (*valSlop)) ;
   }

	// создание таблицы с регулярным шагом
	// аргумент меняется от 0 до ARgTabMax с шагом STepTab

		//*pvalTabDiap = 0.975 *fabs(arrValueArgs[ NUmValuePoints - 1]) ;
		*pvalTabDiap =fabs(arrValueArgs[ NUmValuePoints - 1]) ;

   //	 const double STepTabTemp =  0.03;
	  const double STepTabTemp =  0.1;
	const int NUmPointsRegularTab = (*pvalTabDiap) / STepTabTemp +1;// 2.5/ STepTabTemp +1;
	const double STepTab =  (*pvalTabDiap)/ ((double)NUmPointsRegularTab -1.);

	double * arrRegularTabFunc = new double[NUmPointsRegularTab];



	CreateRegularTab(&arrValueArgs[(NUmValuePoints - 1) /2 ], &arrFunc[(NUmValuePoints - 1) /2 ]
		, (NUmValuePoints - 1) /2 +1 ,arrRegularTabFunc , NUmPointsRegularTab, STepTab );

	TURPolyLine plnRegTargTab( 1, NUmPointsRegularTab) ;
	double *arrArg1 = new double[NUmPointsRegularTab];
   for (int i = 0; i < NUmPointsRegularTab; i++)
   {
	 arrArg1[i]  = ((double)i)* STepTab;
	 plnRegTargTab.Points[i].X = ((double)i)* STepTab;

	 plnRegTargTab.Points[i].Y =  arrRegularTabFunc[i];
   }

	///


	  // полилиния полинома
	double valMaxDelta1 = 0.;
   //	CalcPolinomCoefForOddFunc(NPolinom, arrArg1, arrRegularTabFunc , arrCoeff, NUmPointsRegularTab, &valMaxDelta1);
	 //int nPolin = 6;
	 CalcPolinomCoefApprox(NPolinom, arrArg1, arrRegularTabFunc ,NUmPointsRegularTab,arrCoeff, &valMaxDelta1);


	//double kTang0 = (plnTargTab.Points[NUmPointsTab - 1].Y - plnTargTab.Points[NUmPointsTab - 2 ].Y )
	 //	 / (plnTargTab.Points[NUmPointsTab- 1].X - plnTargTab.Points[NUmPointsTab- 2 ].X );

	// вычисление максимальной ошибки
	TURPointXY *pPntArrDel =  new TURPointXY[plnValue.NumPoints];
	int iLenArrDel = 0;

	 for (int i =0; i < plnValue.NumPoints; i++)
	 {
	   if (plnValue.Points[i].X > (*pvalTabDiap))
	   {
		break;
	   }
	   if (plnValue.Points[i].X < -(*pvalTabDiap))
	   {
		continue;
	   }

	  double valArgTemp = plnValue.Points[i].X;
	  pPntArrDel[iLenArrDel].X = valArgTemp;
	  if (valArgTemp >0)
	  {
	   pPntArrDel[iLenArrDel].Y = (tan(plnValue.Points[i].Y) ) *(*valSlop)

			-(tan(fncPolinom(arrCoeff, NPolinom +1, plnValue.Points[i].X)))*(*valSlop);

	  }
	  else
	  {
		pPntArrDel[iLenArrDel].Y = (tan(plnValue.Points[i].Y) ) *(*valSlop)

			+(tan(fncPolinom(arrCoeff, NPolinom +1, -valArgTemp)))*(*valSlop);
	  }



	   iLenArrDel++;
	 }
	 TURPolyLine plnValue1(pPntArrDel,iLenArrDel);
	 plnValue1.calcBoundBox();
	 *pvalGenMaxError = fabs(plnValue1.Box[1]);
	 if (fabs(plnValue1.Box[3]) > fabs(plnValue1.Box[1]) )
	 {
	   *pvalGenMaxError = fabs(plnValue1.Box[3]);
	 }

   //	double valTemp0= tan(tan(arrRegularTabFunc[NUmPointsRegularTab-1]));
	double valTemp0= (tan(arrRegularTabFunc[NUmPointsRegularTab-1])) * (*valSlop);
	*pvalRazreshenie = transformGeneralizedAngToAng (valTemp0  ) ;
	  //
	 // ВЫВОД ГРАФИКИ
	 if(wchrTableFold != NULL)
	 {
	 // создание shape  файла линии значений
	wchar_t wchrValue0[400] = {0};
	wcscpy( wchrValue0, wchrTableFold);
	wcscat(wchrValue0, L"\\Value0.shp");
	TURPolyLine::WriteSetSHPFiles(wchrValue0,&plnValue0,  1) ;
	///

		 // создание shape  файла линии значений
	wchar_t wchrTargTablePoliLine[400] = {0};
	wcscpy( wchrTargTablePoliLine, wchrTableFold);
	wcscat(wchrTargTablePoliLine, L"\\ValuesLine.shp");
	TURPolyLine::WriteSetSHPFiles(wchrTargTablePoliLine,&plnValue,  1) ;
	///

	// полилилиния регулярной таблицы
	wchar_t wchrRegTargTablePoliLine[400] = {0};
	wcscpy( wchrRegTargTablePoliLine, wchrTableFold);
	wcscat(wchrRegTargTablePoliLine, L"\\RegTargTablePoliLine.shp");
	TURPolyLine::WriteSetSHPFiles(wchrRegTargTablePoliLine,&plnRegTargTab,  1) ;
	 ///


	wchar_t wchrPolinomPoliLine1[400] = {0};
	wcscpy( wchrPolinomPoliLine1,wchrTableFold);
	wcscat(wchrPolinomPoliLine1, L"\\Polinom1.shp");
	TURPolyLine::createPolinomGraph(wchrPolinomPoliLine1, NPolinom, arrCoeff
	    , arrArg1[NUmPointsRegularTab-1], arrArg1[0] ) ;
	///
	// постронеие графика ошибки аппроксимации полиномом
	 wchar_t wchrDeltaApprox[400] = {0};
	wcscpy( wchrDeltaApprox, wchrTableFold);
	wcscat(wchrDeltaApprox, L"\\DelMilRad.shp");
	TURPolyLine::WriteSetSHPFiles(wchrDeltaApprox,&plnValue1,  1) ;
	 }
	///
	delete arrValueArgs ;
	delete arrFunc ;
	delete arrFunc0 ;
	delete arrArg1;
	delete [] pPntArrDel ;
	delete arrRegularTabFunc;
}


void TParAnt::CreateValueArraysEllipsType( const double VAlDiap,double *arrArgs,double *arrFunc, const int NUmPoints0 )
{
	const double VAlGenAngSdvig= transformAngToGeneralizedAng (mAngSdvig ) ;
  TURPointXY pntSdvig ;
  double arrMtxPer_From_b_To_Image_b[4] ={0.};
  calcTransformationParams(  &pntSdvig, arrMtxPer_From_b_To_Image_b);
  const double STep = 2. * VAlDiap/ ((double)NUmPoints0 -1.);
  TURPointXY pntEnvel;
  for (int i = 0; i < NUmPoints0; i++)
  {
  	double valMu = -VAlDiap + ((double)i) *  STep;
	calcEnvelopePoint(valMu, pntEnvel);
	arrFunc[i] = valMu;
	TURPointXY pntImageEnvel =  pntEnvel.LinTransform( pntSdvig,arrMtxPer_From_b_To_Image_b );
	arrArgs[i] = atan2(pntImageEnvel.Y, pntImageEnvel.X);
  }
}

void TParAnt::CreateValueArraysEllipsType__( const double VAlDiap,double *arrArgs,double *arrFunc, const int NUmPoints0 )
{
	const double VAlGenAngSdvig= transformAngToGeneralizedAng (mAngSdvig ) ;
  TURPointXY pntSdvig ;
  double arrMtxPer_From_b_To_Image_b[4] ={0.};
  calcTransformationParams(  &pntSdvig, arrMtxPer_From_b_To_Image_b);
  const double STep = 2. * VAlDiap/ ((double)NUmPoints0 -1.);
  TURPointXY pntEnvel;
  for (int i = 0; i < NUmPoints0; i++)
  {
  	double valMu = -VAlDiap + ((double)i) *  STep;
	calcEnvelopePoint(valMu, pntEnvel);
	arrFunc[i] = atan(atan(valMu));
	TURPointXY pntImageEnvel =  pntEnvel.LinTransform( pntSdvig,arrMtxPer_From_b_To_Image_b );
	arrArgs[i] = (atan2(pntImageEnvel.Y, pntImageEnvel.X));
  }
}


// задан массив значений некоторой функции
// arrArg [LEnArr] - массив знаячений аргументов
// arrFunc   - массив соответствующих знаячений функции
//требуется создать массив значений функции на регулярной сетке аргументов с шагом  STepTab
// OUTPUT:
// arrRegularTabFunc -массив значений функции на сетке с шагом  STepTab
//  рассчитывается линейной интерполяцией
// ВНИМАНИЕ!!!
//  перед обращением к функции для массива arrRegularTabFunc должна быть зарезервирована память
//  LEnRegTab * sizeof(double)
// кроме того, (LEnRegTab -1) *  STepTab должно быть равно диапазону изменения аргумента, то есть
// равно  fabs(arrArg[0] - arrArg[LEnArr-1])
bool  TParAnt::CreateRegularTab(double *arrArg, double *arrFunc, const int LEnArr
  ,double *arrRegularTabFunc , const int LEnRegTab, const double STepTab)
{

 for (int i = 0; i < LEnRegTab; i++)
 {
	double arg = ((double)i) * STepTab;
	int numCur = 0;
	for (int j = 0; j < LEnArr; j++)
	{
	  if ((arg >= arrArg[j])&& (arg <= arrArg[j +1]) )
	  {
	   numCur = j;
	   break;
	  }
	}
	TURPointXY  pnt1(arrArg[numCur],arrFunc[numCur]);
	TURPointXY  pnt2(arrArg[numCur +1],arrFunc[numCur +1]);
	TLine LineTemp(  pnt1,   pnt2);
	if (!LineTemp.calcY(arg, &arrRegularTabFunc[i]) )
	{
		break;
		return false;
    }
 }
 return true;
}

// запись 2 мерного массива типа double в текстовый файл подряд
 // разделитель между числами - ,
 // знак десятичного разделителя - . (десятичная точка)
 int TParAnt::WriteTXTReportForPolinomCoefArray(const wchar_t*FileName,	double *arrOutMtrxTransf, double *arrOutVect00
 , double *arrOutPolinomCoefEllipticCase, const int numHyperbolic ,const int numElliptic
  , double *arrOutBoundsHyperbCaseLeft,double *arrOutPolinomCoefHyperbCaseLeft
  ,double *arrOutBoundsHyperbCaseRight,double *arrOutPolinomCoefHyperbCaseRight,const int NPolinom)

 {

 char ch[300] ;
	FILE *fw ;


	if ((fw = _wfopen(FileName,L"w"))== NULL)
	{
	 String St =  FileName ;
	 ShowMessage(L"TYrWrite::WriteMassiveInFltFile\nERROR ! Not possible to open " +St) ;
	 return 1 ;
	}

		   fprintf(fw,"ТАБЛИЦА МАТРИЦ ПЕРЕХОДА ПО ЧАСТОТАМ\n\n");
		   fprintf(fw, "ЧАСТОТА       A11       A12       A21       A22 \n");
		   fncPrintTab(fw, arrOutMtrxTransf,  numElliptic + numHyperbolic, 5 );

			fprintf(fw,"\n\nТАБЛИЦА ТОЧЕК СДВИГА НАЧАЛА КООРДИНАТ ПО ЧАСТОТАМ\n\n");
		   fprintf(fw, "ЧАСТОТА       X0         Y0\n");
		   fncPrintTab(fw, arrOutVect00,  numElliptic + numHyperbolic, 3 );


		   char str[2000] ={0};
		   fprintf(fw,"\n\n       ЭЛЛИПТИЧЕСКИЙ СЛУЧАЙ\n\n");
		   fprintf(fw,"  ТАБЛИЦА КОЭФФИЦИЕНТОВ ПОЛИНОМОВ СТЕПЕНИ %i \n\n", NPolinom);
		   fprintf(fw, "ЧАСТОТА   Диап.арг.  Разреш,мрад   valKSlop       1 ");
			for (int ii = 0; ii < NPolinom; ii++)
			{
			  fprintf(fw, "       X^%i", ii+1);
			}
			fprintf(fw, " \n");
			fncPrintTab(fw, arrOutPolinomCoefEllipticCase,  numElliptic, NPolinom+1 + 4 );

			fprintf(fw,"\n\n   ГИПЕРБОЛИЧЕСКИЙ СЛУЧАЙ\n\n");
			fprintf(fw,"\n\nТАБЛИЦА ОБЛАСТЕЙ ОПРЕДЕЛЕНИЯ ДЛЯ ЛЕВОЙ ВЕТВИ\n");
			fprintf(fw,"\n\nЛевая ветвь гиперболы. Формат таблицы:\n");
			fprintf(fw,"первый столбец - частота\n");
			fprintf(fw,"второй столбец - tmin -нижняя грантца области определения переменной полинома P(t)\n");
			fprintf(fw,"третий столбец - tmax -верхняя грантца области определения переменной полинома P(t)\n");
			fprintf(fw,"четвертый столбец - TettaMin -нижняя граница области значений полинома P(t)\n");
			fprintf(fw,"пятый столбец - TettaMax - верхняя граница области значений полинома P(t)\n\n");
			fprintf(fw, "ЧАСТОТА     tmin      tmax      TettaMin  TettaMax \n");
			fncPrintTab(fw, arrOutBoundsHyperbCaseLeft,  numHyperbolic, 5 );

			 fprintf(fw,"\n\nТАБЛИЦА КОЭФФИЦИЕНТОВ АППРОКСИМИРУЮЩЕГО ПОЛИНОМА НЕЧЕТНОЙ СТЕПЕНИ \n");
			 fprintf(fw,"ДЛЯ ЛЕВОЙ ВЕТВИ ГИПЕРБОЛЫ. АППРОКСИМАЦИЯ НА ИНТЕРВАЛЕ [tmin;tmax]\n\n");
			 fprintf(fw, "ЧАСТОТА ");
			 for (int ii = 0; ii <= NPolinom; ii++)
			{
			  fprintf(fw, "       X^%i", ii*2 +1);
			}
			fprintf(fw, " \n");
			fncPrintTab(fw, arrOutPolinomCoefHyperbCaseLeft,  numHyperbolic, NPolinom+2 );

			fprintf(fw,"\n\nТАБЛИЦА ОБЛАСТЕЙ ОПРЕДЕЛЕНИЯ ДЛЯ ПРАВОЙ ВЕТВИ\n");
			fprintf(fw,"\n\nПравая ветвь гиперболы. Формат таблицы:\n");
			fprintf(fw,"первый столбец - частота\n");
			fprintf(fw,"второй столбец - tmin -нижняя грантца области определения переменной полинома P(t)\n");
			fprintf(fw,"третий столбец - tmax -верхняя грантца области определения переменной полинома P(t)\n");
			fprintf(fw,"четвертый столбец - TettaMin -нижняя граница области значений полинома P(t)\n");
			fprintf(fw,"пятый столбец - TettaMax - верхняя граница области значений полинома P(t)\n\n");
			fprintf(fw, "ЧАСТОТА     tmin      tmax      TettaMin  TettaMax \n");
			fncPrintTab(fw, arrOutBoundsHyperbCaseRight,  numHyperbolic, 5 );

			 fprintf(fw,"\n\nТАБЛИЦА КОЭФФИЦИЕНТОВ АППРОКСИМИРУЮЩЕГО ОТРЕЗКА РЯДА ЛОРАНА \n");
			 fprintf(fw,"ДЛЯ ПРАВОЙ ВЕТВИ ГИПЕРБОЛЫ. АППРОКСИМАЦИЯ НА ИНТЕРВАЛЕ [tmin;tmax]\n\n");
			 fprintf(fw, "ЧАСТОТА     1/(X^2)       1/X       1");
			 for (int ii = 0; ii < NPolinom-2; ii++)
			{
			  fprintf(fw, "        X^%i", ii);
			}
			fprintf(fw, " \n");
			fncPrintTab(fw, arrOutPolinomCoefHyperbCaseRight,  numHyperbolic, NPolinom+2 );



			fprintf(fw, "АЛГОРИТМ ВЫЧИСЛЕНИЙ\n");

			fprintf(fw, "1. ВЫЧИСЛЕНИЕ ВЕКТОРА b.\n ");
			fprintf(fw, "INPUT:\n");
			fprintf(fw, "cmparrS[3] - массив измерений тройки диаграммn");
			fprintf(fw, "OUTPUT:\n");
			fprintf(fw, "arr_b[2] - вектор b \n");
			fprintf(fw, "bool TParAnt::calcVect_b(TComp *cmparrS, double *arr_b)\n");
			fprintf(fw, "{ double arrS[4] = {0.};\n   arrS[0] = cmparrS[0].m_Re;\n arrS[1] = cmparrS[2].m_Re;\n  arrS[2] = cmparrS[0].m_Im;\n arrS[3] = cmparrS[2].m_Im;\n double arrs[2] = {0.};\n ");
			fprintf(fw, "{ arrs[0] = cmparrS[1].m_Re;\n arrs[1] = cmparrS[1].m_Im;\n");
			fprintf(fw,"double  arrSInv[4] = {0.};\n  bool brez =   InverseMtrx2(arrS, arrSInv);\n MtrxMultMatrx(arrSInv,2, 2, arrs,1, arr_b) ;  \n");
			fprintf(fw," return brez;\n}\n");




			fprintf(fw, "INPUT:\n");
			fprintf(fw, "INPUT:\n");
			fprintf(fw, "INPUT:\n");
			fprintf(fw, "INPUT:\n");
			fprintf(fw, "INPUT:\n");
			fprintf(fw, "INPUT:\n");
			fprintf(fw, "INPUT:\n");
			fprintf(fw, "INPUT:\n");
			fprintf(fw, "1. ВЫЧИСЛЕНИЕ ВЕКТОРА b");
			fprintf(fw, "1. ВЫЧИСЛЕНИЕ ВЕКТОРА b");
			fprintf(fw, "1. ВЫЧИСЛЕНИЕ ВЕКТОРА b");
			fprintf(fw, "1. ВЫЧИСЛЕНИЕ ВЕКТОРА b");
			fprintf(fw, "1. ВЫЧИСЛЕНИЕ ВЕКТОРА b");
			fprintf(fw, "1. ВЫЧИСЛЕНИЕ ВЕКТОРА b");

			 /*

// вычисление вектора b
// INPUT:
// cmparrS[3] - массив измерений тройки диаграмм
// OUTPUT:
// arr_b[2] -  коэффициенты уравнения
bool TParAnt::calcVect_b(TComp *cmparrS, double *arr_b)
{
 // вычисленние вектора b
   double arrS[4] = {0.};
   arrS[0] = cmparrS[0].m_Re;
   arrS[1] = cmparrS[2].m_Re;
   arrS[2] = cmparrS[0].m_Im;
   arrS[3] = cmparrS[2].m_Im;
   double arrs[2] = {0.};
   arrs[0] = cmparrS[1].m_Re;
   arrs[1] = cmparrS[1].m_Im;

   double  arrSInv[4] = {0.};
   bool brez =   InverseMtrx2(arrS, arrSInv);
   MtrxMultMatrx(arrSInv,2, 2, arrs,1, arr_b) ;
   return brez;
}


			 */
			fclose(fw);

 }


void fncPrintTab(FILE *fw, double *parr,  int nrows, int ncols )
{

char ch[100] = {0};
 for( int i =0; i< nrows; i++)
{
			 for (int j = 0; j < ncols; j++)
			 {


				if ((i + 1) * (j +1) != nrows *ncols)
				{
				  if (j == (ncols -1))
				  {
				   sprintf(ch,"%-6.10f,\n",parr[ i * ncols + j]) ;
				  }
				  else
				  {
				  sprintf(ch,"%-6.10f, ",parr[ i * ncols + j]) ;
				  }
				}
			   else
			   {
				sprintf(ch,"%-6.10f\n",parr[ i * ncols + j]) ;
			   }

				fprintf(fw,"%s",ch);

			 }
 }
}
// построение огибающей фигуры гиперболического типа и ее
// преобразование  в гиперболу xy=1
// INPUT:
// VAlGenAngSdvig - обобщенный угол сдвига менжду диаграммаи
// VAlDiap - диапазон сетки, от  -VAlDiap до VAlDiap
// NUmPoints0 - к-во точек
// OUTPUT:
// *ppntSdvig - точка сдвига начала координат (-центр огибающей)
// arrMtxPer_From_b_To_Image_b - матрица линейного преобразования, переводящего
// огибающую фигуру в  в гиперболу xy=1
//
void TParAnt::transformArrayEnvelopePoints(const double VAlDiap
  , const int NUmPoints0, TURPointXY *arrPoints0, int *numPointsActual, TURPointXY *pntSdvig, double *arrMtxPer_From_b_To_Image_b)
{
	const double VAlGenAngSdvig= transformAngToGeneralizedAng (mAngSdvig ) ;
  TURPointXY *arrPoints = new TURPointXY[ NUmPoints0];
  *numPointsActual = 0;
  calcEvelopePointsArray( VAlDiap,  NUmPoints0, arrPoints, numPointsActual );

 calcTransformationParams( pntSdvig , arrMtxPer_From_b_To_Image_b);
	for (int i = 0; i < *numPointsActual; i++)
	{
	  TURPointXY pntTemp = arrPoints[i].SdvigTransform((*pntSdvig) );
	  arrPoints0 [i] = pntTemp.fncLinTrasform(arrMtxPer_From_b_To_Image_b);
	}

	delete []arrPoints;
   return ;

}
// функции для огибающей гиперболического типа
// вычисление массива точек огибающей на сетке углов mu
// INPUT:
// VAlGenAngSdvig - обобщенный угол сдвига менжду диаграммаи
// VAlDiap - диапазон сетки, от  -VAlDiap до VAlDiap
// NUmPoints0 - длинна зарезервированной памаяти массива точек arrPoints
// OUTPUT:
// arrPoints [NUmPoints0 ] - массив точек
// *numPointsActual - реальное количество тточек
void TParAnt::calcEvelopePointsArray( const double VAlDiap, const int NUmPoints0, TURPointXY *arrPoints, int *numPointsActual )
{
	const double VAlGenAngSdvig= transformAngToGeneralizedAng (mAngSdvig ) ;
   *numPointsActual = 0;
  const double STep = 2. * VAlDiap/ ((double)NUmPoints0);
  for (int i = 0; i < NUmPoints0; i++)
  {
	double valMu = -VAlDiap + ((double)i) *  STep;

	if (calcEnvelopePoint(valMu, arrPoints[*numPointsActual]))
	{
	  (*numPointsActual)++;
	}
  }

  return ;
}


// нахождение точки (b1,b2) огибающей семейства кривых
// F(mu)- b1 F(mu+ alf)- b2 F(mu- alf) =0
// в точке mu
// для метода 3 парциальных дщиаграмм
// INPUT:
//  VAlMu - обобще6нный угол
// VAlGenAngSdvig - обобщенный угол сдвига менжду диаграммаи
//  OUTPUT:
// pntB - точb1ка
bool TParAnt::calcEnvelopePoint(const double VAlMu, TURPointXY &pntB)
{
	const double VAlGenAngSdvig= transformAngToGeneralizedAng (mAngSdvig ) ;
 // формирование матрицы A и векирора правой части b
 double arrA[4] = {0.}, arrb[2] = {0.}, arrX[2] ={0.};
  arrA[0] = fncDiagrSinx_div_x(VAlMu + VAlGenAngSdvig) ;
  arrA[1] = fncDiagrSinx_div_x(VAlMu - VAlGenAngSdvig) ;
  arrA[2] = fncDerivDiagrSinx_div_x(VAlMu + VAlGenAngSdvig) ;
  arrA[3] = fncDerivDiagrSinx_div_x(VAlMu - VAlGenAngSdvig) ;

  arrb[0] =  fncDiagrSinx_div_x(VAlMu);
  arrb[1] =  fncDerivDiagrSinx_div_x(VAlMu);
  ///

   ///
  if( !SolvLinEq2(arrA, arrb,arrX))
  return false;
  pntB.X = arrX[0];
  pntB.Y = arrX[1];

  return true;
}

// вычисление массива точек огибающей на сетке углов mu
// INPUT:
// VAlGenAngSdvig - обобщенный угол сдвига менжду диаграммаи
// VAlDiap - диапазон сетки, от  -VAlDiap до VAlDiap
// NUmPoints0 - длинна зарезервированной памаяти массива точек arrPoints
// OUTPUT:
// arrPoints [NUmPoints0 ] - массив точек
// numPointsActual - реальное количество тточек
TURPolygon TParAnt::calcEvelopePolyg( const double VAlDiap, const int NUmPoints0 )
{
	const double VAlGenAngSdvig= transformAngToGeneralizedAng (mAngSdvig ) ;
  TURPointXY *arrPoints = new TURPointXY[ NUmPoints0];
  int numPointsActual = 0;

  calcEvelopePointsArray( VAlDiap ,  NUmPoints0, arrPoints, &numPointsActual );


   TURPolygon plgEnvel(  numPointsActual,arrPoints);

  delete []arrPoints;
  return plgEnvel;

}

// построение огибающей фигуры и ее
// преобразование  в единичную окружнсть с центром (0,0)
// INPUT:
// VAlGenAngSdvig - обобщенный угол сдвига менжду диаграммаи
// VAlDiap - диапазон сетки, от  -VAlDiap до VAlDiap
// NUmPoints0 - к-во точек
// OUTPUT:
// plgEnvelopeCentred  - отцентрованный  в начало координат огибающий пполигон
// *ppntSdvig - точка сдвига начала координат (-центр огибающей)
// arrMtxPer_From_b_To_Image_b - матрица линейного преобразования, переводящего
// огибающую фигуру(усеченный эллипс) в усеченную окружность

TURPolygon TParAnt::transformEnvelopePolygon( const double VAlDiap, const int NUmPoints0
  , TURPolygon &plgEnvelopeRotated, TURPointXY *ppntSdvig, double *arrMtxPer_From_b_To_Image_b)
{
	const double VAlGenAngSdvig= transformAngToGeneralizedAng (mAngSdvig ) ;
  TURPolygon plgEnvelope = calcEvelopePolyg(   VAlDiap, NUmPoints0 ) ;

   TURPointXY pntCentre =  calcCentrePoint(plgEnvelope);
   *ppntSdvig = TURPointXY (-pntCentre.X, -pntCentre.Y);
   TURPolygon plgEnvelopeCentred =   plgEnvelope.SdvigTransform(*ppntSdvig );
   int iNumSemiAxe = findMaxSemiAxe(plgEnvelopeCentred);
   TURPointXY pntMax = plgEnvelopeCentred.Points[ iNumSemiAxe ] ;

   // мееньшая полуось
   TURPointXY pntSegmTemp0,pntMin, pnt0(0.,0.);
   double arrVect[2] ={0.};
   arrVect[0]=  pntMax.Y;
   arrVect[1]=  -pntMax.X;
   plgEnvelopeCentred.fncFindSecantSegm(  pnt0, arrVect,pntSegmTemp0,pntMin);    ///

  double  arrMtxRotate[4] = {1., 0.
   ,0.,1. };
   double arrMtxCompress[4] = {1., 0.
   ,0.,1.
   };
	// маьтрица поворота
	 double valMaxSemiAxe =  pntMax.Norm() ;
	 double valMinSemiAxe =  pntMin.Norm() ;

	double alf = M_PI/4. + M_PI;

	arrMtxRotate [0] = cos(alf);
	arrMtxRotate [1] = sin (alf);
	arrMtxRotate [2] = -sin (alf);
	arrMtxRotate [3] = cos(alf);
	plgEnvelopeRotated =  plgEnvelopeCentred.fncLinTransform(arrMtxRotate )  ;
	arrMtxCompress[0] = 1./valMaxSemiAxe;
	arrMtxCompress[3] = 1./ valMinSemiAxe;
	MtrxMultMatrx(arrMtxCompress,2 ,2 , arrMtxRotate,2, arrMtxPer_From_b_To_Image_b)  ;
	///

   TURPolygon plgTransformed = plgEnvelopeCentred.fncLinTransform(arrMtxPer_From_b_To_Image_b )  ;

   return plgTransformed;
}

bool TParAnt::calc_b2( const double VAlMu, const double VAl_b1, double *pval_b2)
{
	const double VAlGenAngSdvig= transformAngToGeneralizedAng (mAngSdvig ) ;
  double valA0 = fncDiagrSinx_div_x(VAlMu + VAlGenAngSdvig) ;
  double valA1 = fncDiagrSinx_div_x(VAlMu - VAlGenAngSdvig) ;
  double valA2 = fncDiagrSinx_div_x(VAlMu ) ;
	if (fabs(valA1) < DBL_MIN)
	{
	return false;
	}
	*pval_b2 = ( valA2 -  VAl_b1 * valA0) /valA1;
	return true;
}


double TParAnt::fncCurvation()
{
	const double VAlGenAngSdvig= transformAngToGeneralizedAng (mAngSdvig ) ;
	double valF     = fncDiagrSinx_div_x(VAlGenAngSdvig) ;
	double valFder1 = fncDerivDiagrSinx_div_x(VAlGenAngSdvig) ;
	double valFder2 = fncDeriv2DiagrSinx_div_x(VAlGenAngSdvig) ;
	double valFder3 = fncDeriv3DiagrSinx_div_x(VAlGenAngSdvig) ;
	double val_k = -sqrt(2.) * valFder1 * (valF *(valFder3 + 1./3.* valFder1 ) +valFder1*valFder2  )
	/(valFder2 + 1./3.* valF)	/(valFder2 + 1./3.* valF);
 return  val_k ;
}

bool TParAnt::findEnvelopeSigularPoint (double *pvalMuSingular, double *pvalBSingular)
{
	const double VAlGenAngSdvig= transformAngToGeneralizedAng (mAngSdvig ) ;
	 double valMu0 = M_PI/3.;
	double f0 = fncFi(valMu0 ) ;
	 double valMu1 = valMu0;
	 double f1 =0.;
	 double valMuRez =  valMu0 ;
  double eps = 0.0001;
  for (int i = 0; i < 100000; i++)
  {
	valMu1 = valMu0 +0.001;
	 f1 = fncFi(valMu1 ) ;
	 if(f1 * f0 <= 0.)
	 {
	  *pvalMuSingular =  solvSingularEquationMethChord( valMu0, valMu1) ;
	   *pvalBSingular = fncDiagrSinx_div_x(*pvalMuSingular)
	  /(fncDiagrSinx_div_x(*pvalMuSingular +VAlGenAngSdvig) + fncDiagrSinx_div_x(*pvalMuSingular -VAlGenAngSdvig));

	   return true;
	 }
	 valMu0 = valMu1;
   }
  return false;
}

double TParAnt::fncFi(const double VAlMu0 )
{
	const double VAlGenAngSdvig= transformAngToGeneralizedAng (mAngSdvig ) ;

	double temp0 = fncDerivDiagrSinx_div_x(VAlMu0)
				  *(fncDiagrSinx_div_x(VAlMu0 + VAlGenAngSdvig) + fncDiagrSinx_div_x(VAlMu0 - VAlGenAngSdvig));
	double temp1 = fncDiagrSinx_div_x(VAlMu0)
				  *(fncDerivDiagrSinx_div_x(VAlMu0 + VAlGenAngSdvig) + fncDerivDiagrSinx_div_x(VAlMu0 - VAlGenAngSdvig));

	return temp0 - temp1;

}



double 	  TParAnt::solvSingularEquationMethChord( double  valX0, double valX1)
{
  const double VAlGenAngSdvig= transformAngToGeneralizedAng (mAngSdvig ) ;
  double valMu0 = valX0 ;
  double valMu1 = valX1 ;
  double valMuRez =  valMu0 ;
  double eps = 0.0001;
  for (int i = 0; i < 100; i++)
  {
	double f0 = fncFi( valMu0 ) ;
	double f1 = fncFi(valMu1 ) ;
	double valMuRez0 = valMu0 - f0 * ( valMu1 - valMu0)/ (f1 - f0);
	if (fabs(valMuRez0 - valMuRez) < eps)
	{
	 return valMuRez0;
	}
	double f2 = fncFi(  valMuRez0 ) ;
	if( f2* f0 < 0.)
	{
	   valMu1 = valMuRez0;

	}
	else
	{
      valMu0 = valMuRez0;
    }
     valMuRez = valMuRez0;
  }
	return -100000.;
}

void TParAnt::ShowTheoreticRotatedEnvelopeEllipseType(wchar_t *FileName )
{
  const double VAlGenAngSdvig= transformAngToGeneralizedAng (mAngSdvig ) ;

	double val_k = fncCurvation( ) ;
	double valMuSingular = 0., valBSingular =0.;
	bool brez11 = findEnvelopeSigularPoint ( &valMuSingular, &valBSingular) ;
	double valBRegular = 1./ fncDiagrSinx_div_x(VAlGenAngSdvig) /2.;
	double vala = fabs(( valBRegular -  valBSingular)) / sqrt(2.)-0.1675;
 double valb = sqrt(vala /fabs(val_k ))+0.13 ;

  TArcEllipse  ellips( vala, valb);
  ellips.ShowFullGraph(FileName) ;

}

 void TParAnt::ShowTheoreticRotatedEnvelopeHyperbolaType(wchar_t *FileName )
{
  const double VAlGenAngSdvig= transformAngToGeneralizedAng (mAngSdvig ) ;

	double val_k = fncCurvation() ;
	double valBSingular = -0.5;

	double valBRegular = 1./ fncDiagrSinx_div_x(VAlGenAngSdvig) /2.;
	double vala = fabs(( valBRegular -  valBSingular)) / sqrt(2.);
 double valb = sqrt(vala /fabs(val_k ));// -0.13 ;

  THyperbola  hyperb( vala, valb);
  hyperb.ShowMe(FileName, 10000, 700.) ;

}


//  нахождение    центра огибающей фигуры

TURPointXY  calcCentrePoint( TURPolygon plgEnvelope)
{
	 TURPointXY pntCent0(-5.5, -5.5), pntCur;

	 double arrVectVert[2] = {0.,  1.
	                         };
	 double  arrVectHor[2] = {1.,  0.
							  };
	 for (int i =0; i < 1000; i++)
	 {
		TURPointXY pntSegmTemp0,pntSegmTemp1;

		plgEnvelope.fncFindSecantSegm(  pntCent0, arrVectVert,pntSegmTemp0,pntSegmTemp1);
		TURPointXY pntCent1 =  TURPointXY::ParamPoint(pntSegmTemp0,pntSegmTemp1,0.5);
		plgEnvelope.fncFindSecantSegm(  pntCent1, arrVectHor
		,pntSegmTemp0,pntSegmTemp1);
		TURPointXY pntCent2 =  TURPointXY::ParamPoint(pntSegmTemp0,pntSegmTemp1,0.5);
		double distTEmp = TURPointXY::dist(pntCent0,pntCent2);
		pntCent0 = pntCent2 ;

		if (distTEmp < 0.001)
		{
		break;
		}

	 }

   return pntCent0;
}


int findMaxSemiAxe( TURPolygon plgEnvelopeCentred)
{
  double valMax = 0.;
  int ireturn = -1;
  for (int i =0 ; i < plgEnvelopeCentred.NumPoints; i++)
  {
	double valTemp =  plgEnvelopeCentred.Points[i].Norm();
	if (valTemp > valMax )
	{
	   valMax = valTemp;
	   ireturn = i;
	}
  }
  return ireturn;
}
// рекуррентный пересчет статистических средненго и дисперсии случайной величины
// valX - значение случ величины , полученное на испытании с номером N
// N - номиер текущего испытания (начинается с 2)
// valAver - среднее на испытании  N-1
// valAverSquare - среднее квадрата сл величины на испытании  N-1
//
// OUTPUT:
// valAver - среднее на испытании  N
// valAverSquare - среднее квадрата сл величины на испытании  N
//
void CalcStatParams_(const double valX, const int N
	, double &valAver,double  &valAverSquare)
{
   if (N == 1)
   {
	valAver = valX;
	valAverSquare = valX * valX;
	return;
   }
   valAver = valAver * (N - 1.)/ N + valX / N ;
   valAverSquare = valAverSquare * (N - 1.)/ N +  valX * valX / N ;

}

// вычисление СКО через мат ожидание и матожидание квадрата сл. вел
// CKO = sqrt( M(X*X) - ( MX ) * ( MX ))
double fncSKO_(
		 const double valXAver  //INPUT:  мат. ожидание сл. вел X
		,const double valX2Aver  //INPUT:  мат. ожидание сл. вел X*X
		)
{
	double temp =  valX2Aver - valXAver * valXAver;
	double tenp0 = (temp >= 0.)? temp:0.;
	return sqrt( tenp0) ;
}


void fncFlipDblArr(double *arrFunc , const int NUmPointsTab)
{
	for (int i = 0; i  < NUmPointsTab /2; i++)
	{
	 double temp = arrFunc [i];
	  arrFunc [i] =  arrFunc [NUmPointsTab - 1- i];
	  arrFunc [NUmPointsTab - 1- i] = temp;
	}
}

 void CalcPolinomCoefForOddFunc(const int NPolinom, double *arrArg,double *arrF
   ,double * arrCoeff, int lenArr, double *pMaxDelta)
 {
   const int NRows= (NPolinom-1) /2 +1;
   long double *arrMtrxX = new long double [NRows * lenArr ];
   long double *arrMtrxX1 = new long double [NRows * lenArr ];
   long double *arrLongF  = new long double [NRows * lenArr ];
   long double *arrLongCoeff  = new long double [NRows  ];
   memcpy(arrMtrxX, arrArg, lenArr * sizeof(double));
   for (int i =0; i < lenArr; i++)
   {
	 arrMtrxX[i] = (long double)  arrArg[i];
	 arrLongF[i] =  (long double) arrF[i];
   }
	for(int i=1;i < NRows;i++)
	{
		for (int j = 0; j < lenArr; j++)
		{
		  arrMtrxX[ i * lenArr + j] = arrMtrxX[ (i -1) * lenArr + j] * arrMtrxX[j]* arrMtrxX[j];
		}
	}
	long double coef =1.;
	long double *parrRez = new long double [NRows * lenArr ];
	MatrxDivideScalar(arrMtrxX, NRows * lenArr, 1, coef, parrRez);
	memcpy( arrMtrxX, parrRez ,  NRows * lenArr * sizeof(long double));

	memcpy(arrMtrxX1, arrMtrxX, NRows * lenArr* sizeof(long double));
	long double *arrMtrxMult = new long double [NRows * NRows];
	MtrxMultMatrxTransp(arrMtrxX,NRows, lenArr, arrMtrxX1,NRows, arrMtrxMult)  ;
	long double *arrXF = new long double [NRows];
	MtrxMultMatrx(arrMtrxX,NRows, lenArr, arrLongF,1, arrXF);

	long double *arrMtrxMultInv = new long double [NRows * NRows];
	InverseMtrx( arrMtrxMult, NRows, arrMtrxMultInv ) ;


	MtrxMultMatrx(arrMtrxMultInv,NRows, NRows, arrXF,1, parrRez);
	MatrxDivideScalar(parrRez,  NRows , 1, coef, arrLongCoeff);
	for (int i =0; i < NRows; i++)
	{
	  arrCoeff[i] = (double)arrLongCoeff[i];
	}
   delete arrMtrxX;
   delete arrMtrxMult;
   delete arrXF;
   delete arrMtrxMultInv;
   delete arrMtrxX1;
   delete parrRez;
   delete arrLongF;
   delete arrLongCoeff;
	*pMaxDelta = 0.;
	for (int j = 0; j < lenArr; j++)
	{
	   double valDelta = fabs( arrF[j] - fncOddPolinom(arrCoeff,  NRows,arrArg[j]));
	   if ((*pMaxDelta) < valDelta)
	   {
		 (*pMaxDelta) = valDelta ;
	   }
	}

 }
 void CalcPolinomCoefApprox(const int NPolinom, double *arrArg,double *arrF,const int  lenArr
   ,double * arrCoeff,  double *pMaxDelta)
 {
   const int NRows= NPolinom +1;

   long double *arrMtrxX = new long double [NRows * lenArr ];
   long double *arrMtrxX1 = new long double [NRows * lenArr ];
   long double *arrLongF  = new long double [NRows * lenArr ];
   long double *arrLongCoeff  = new long double [NRows  ];
   memcpy(arrMtrxX, arrArg, lenArr * sizeof(double));
   for (int i =0; i < lenArr; i++)
   {
	 arrMtrxX[i] = 1.;
	 arrLongF[i] =  (long double) arrF[i];
   }
	for(int i=1;i < NRows;i++)
	{
		for (int j = 0; j < lenArr; j++)
		{
		  arrMtrxX[ i * lenArr + j] = arrMtrxX[ (i -1) * lenArr + j] * ((long double) arrArg[j]);
		}
	}
	long double coef =1.;
	long double *parrRez = new long double [NRows * lenArr ];
	MatrxDivideScalar(arrMtrxX, NRows * lenArr, 1, coef, parrRez);
	memcpy( arrMtrxX, parrRez ,  NRows * lenArr * sizeof(long double));

	memcpy(arrMtrxX1, arrMtrxX, NRows * lenArr* sizeof(long double));
	long double *arrMtrxMult = new long double [NRows * NRows];
	MtrxMultMatrxTransp(arrMtrxX,NRows, lenArr, arrMtrxX1,NRows, arrMtrxMult)  ;
	long double *arrXF = new long double [NRows];
	MtrxMultMatrx(arrMtrxX,NRows, lenArr, arrLongF,1, arrXF);

 	long double *arrMtrxMultInv = new long double [NRows * NRows];
 //	InverseMtrx( arrMtrxMult, NRows, arrMtrxMultInv ) ;
 //	MtrxMultMatrx(arrMtrxMultInv,NRows, NRows, arrXF,1, parrRez);

   /////////////////////////////////
   /////////////////////////////////
  //	long double *parrRez1 = new long double [NRows  ];
 //	long double *parrRez2 = new long double [NRows  ];
//	long double *parrRez3 = new long double [NRows  ];
   bool brez = GaussMeth(arrMtrxMult,NRows,  arrXF, parrRez);
 //  MtrxMultMatrx(arrMtrxMult,NRows, NRows, parrRez1,1, parrRez2) ;
  //	  MtrxMultMatrx(arrMtrxMult,NRows, NRows, parrRez,1, parrRez3);
//   delete parrRez1;delete parrRez2;
 //  delete parrRez3;


   /////////////////////////////////
   /////////////////////////////////



	MatrxDivideScalar(parrRez,  NRows , 1, coef, arrLongCoeff);
	for (int i =0; i < NRows; i++)
	{
	  arrCoeff[i] = (double)arrLongCoeff[i];
	}
   delete arrMtrxX;
   delete arrMtrxMult;
   delete arrXF;
   delete arrMtrxMultInv;
   delete arrMtrxX1;
   delete parrRez;
   delete arrLongF;
   delete arrLongCoeff;
	*pMaxDelta = 0.;
	for (int j = 0; j < lenArr; j++)
	{
	   double valDelta = fabs( arrF[j] - fncPolinom(arrCoeff,  NRows,arrArg[j]));
	   if ((*pMaxDelta) < valDelta)
	   {
		 (*pMaxDelta) = valDelta ;
	   }
	}

 }

 double fncOddPolinom(double *arrCoeff, int lenarr, double valArg)
 {

   double valTemp =  valArg;
   double valSum = 0.;
   for (int i = 0; i < lenarr; i++)
   {
	 valSum +=  valTemp *  arrCoeff[i];
	 valTemp = valTemp * valArg * valArg;
   }
   return valSum;
 }

 double fncPolinom(double *arrCoeff, int lenarrCoeff, double valArg)
 {

   double valTemp =  1.;
   double valSum = 0.;
   for (int i = 0; i < lenarrCoeff; i++)
   {
	 valSum +=  valTemp *  arrCoeff[i];
	 valTemp = valTemp * valArg ;
   }
   return valSum;
 }


// оценивание обобщенных углов цели и антипода по измерениям 3 ПАРЦИАЛЬНЫХ диаграмм
//INPUT:
// cmparrS[6] -  массив замеров
// NumRayTriple-  номер тройки рабочих лучей (нумерация лучей идет снизу начиная с нуля)
//                номер тройки это порядковый номер нижнего луча
//  OUTPUT:
//  *valEstAngTarg -   угол цели, рад
// *valEstAngAntp -    угол антипода, рад
// cmpKTarg - коэф отражения цели
// cmpKAntp - коэф отражения антипода
// arrMtrxCorr  - коррел матрица ошибок измерения угла цели и антипода
int TParAnt::tabulatedSolution_5P10_03( TComp *cmparrS
  , int iNumRayTriple , double *valEstAngTarg, double *valEstAngAntp, TComp *cmpKTarg , TComp *cmpKAntp , double *arrMtrxCorr
  , double *pval_b0, double *pval_b1 )
{
	double valGenAngSdvig = transformAngToGeneralizedAng (mAngSdvig ) ;
	// нахождение обобщенных углов цели и антипода
	double  valGenTargEps = -100.,  valGenAntpEps = -100.0;  // обобщ углы цели и антипода
	int irez = tabulatedEstimationGenAngs_5P10_03(  &cmparrS[iNumRayTriple]
	,  &valGenTargEps,  &valGenAntpEps, pval_b0, pval_b1);
	if (irez < 1)
	{
	   return -1;
	}
	//
	// вычисление коэффиц отражения
	TComp cmparrA[4],cmparrAInv[4], cmparrK[2], cmparrSTemp[2];

	cmparrA[0] = TComp(fncDiagrSinx_div_x(valGenTargEps + valGenAngSdvig), 0.);
	cmparrA[1] = TComp(fncDiagrSinx_div_x (valGenAntpEps + valGenAngSdvig ), 0.);
	cmparrA[2] = TComp(fncDiagrSinx_div_x (valGenTargEps - valGenAngSdvig ), 0.);
	cmparrA[3] = TComp(fncDiagrSinx_div_x (valGenAntpEps - valGenAngSdvig ), 0.);

	bool brez =   InverseMtrx2( cmparrA, cmparrAInv);
	if (!brez)
	{
	 return -2;
	}

	cmparrSTemp[0] = cmparrS [iNumRayTriple ];
	cmparrSTemp[1] = cmparrS [iNumRayTriple + 2 ];
	MtrxMultMatrx(cmparrAInv,2, 2, cmparrSTemp,1, cmparrK);
	*cmpKTarg = cmparrK[0] ;
	*cmpKAntp = cmparrK[1] ;

	///

	// вычисление углоы цели и антипода в рад
		  // относиьтельный угол цетрального луча тройки в радианах
			double valRadCentreRayPos = calcRadAngSdvigaPartDiagr(iNumRayTriple + 1);
			// ( - ((double)mQuantDiagr -1.) /2. + 1. + ((double)iNumRayTriple))				*mAngSdvig ;
		  // вычисление угла цели в радианах в АСК
			const double ValRadTargEps =  transformGeneralizedAngToAng ( valGenTargEps ) ;// в радианы
			*valEstAngTarg = ValRadTargEps + valRadCentreRayPos;
		  // вычисление угла антипода в радианах в АПСК
			const double ValRadAntpАEps =  transformGeneralizedAngToAng (valGenAntpEps ) ;// в радианы
			*valEstAngAntp = ValRadAntpАEps + valRadCentreRayPos;
	///

	// вычисление коррел матрицы ошибок определения обобщенных углов цели и антипода
	double arrMtrxCorrMu[4] = {0.}, arrMtrxCorrMu1[4] = {0.};
	calcMtrxCorrGenAngs_Meth3PartDiagr(iNumRayTriple ,valGenTargEps,  valGenAntpEps
	  ,*cmpKTarg,*cmpKAntp, arrMtrxCorrMu1) ;
	calcMtrxCorrGenAngs_AlternativeMeth3PartDiagr(cmparrS,   valGenTargEps,   valGenAntpEps
  ,*cmpKTarg ,*cmpKAntp ,iNumRayTriple, arrMtrxCorrMu);
	///

	// вычисление корреляц матрицы ошибок определения углов цели и антипода в радианаах
	double arrQ[4] = {0.}, arrTemp[4] = {0.};
	double coef = mLambda / mAppert / M_PI;
	arrQ[0] =  coef / sqrt( 1.- coef * valGenTargEps * coef * valGenTargEps);
	arrQ[3] =  coef / sqrt( 1.- coef * valGenAntpEps * coef * valGenAntpEps);
	MtrxMultMatrx(arrQ,2, 2, arrMtrxCorrMu,2, arrTemp)  ;
	MtrxMultMatrxTransp(arrTemp,2, 2, arrQ,2, arrMtrxCorr) ;
	///
	return irez;
}


// оценивание обобщенных углов цели и антипода по измерениям 3 парциальных диаграмм
//INPUT:
// wchFoldName1 - путь к папаке с графиками
// cmparrS[3] -  массив замеров
//  OUTPUT:
//  *pvalGenTargEps -   угол цели, обобщ коорд
// *pvalGenAntpEps -    угол антипода, обобщ коорд
//код возврата
//-1 -  частота не соответствует заданным
//
//
//
int TParAnt::tabulatedEstimationGenAngs_5P10_03(TComp *cmparrS ,  double  *pvalGenTargEps
   ,  double  *pvalGenAntpEps, double *pval_b0, double *pval_b1)
{
   double valGenAngSdvig = transformAngToGeneralizedAng (mAngSdvig ) ;
  // вычисленние вектора b
   double arr_b[2] = {0.}  ;
   bool brez = calcVect_b(cmparrS,arr_b) ;
   ///
   *pval_b0 = arr_b[0];
   *pval_b1 = arr_b[1];
   ///

   // нахожэдение порядкового нгомера частоты, нумерация с 0
   double valFreq= VAL_C * 100./ mLambda/ 1000000000.;
	int  numFreq = -1;
   for (int i = 0; i < ( QUantCasesEllipticType + QUantCasesHyperbolicType); i++)
   {
	 if (fabs(valFreq - constArrMtrxTransf[i * 5]) < 0.01)
	 {
	   numFreq = i;
	   break;
	 }
   }
   if (numFreq == -1)
   {
	return -1;
   }
   ///
  // линейное преобразование вектора b
	// сдвиг
  double arrbT[2] = {0.}, arrTransf_b[2] = {0.};
  arrbT[0] = arr_b[0] + constArrPointsSdvig[ 3 * numFreq +1];
  arrbT[1] = arr_b[1] + constArrPointsSdvig[ 3 * numFreq +2];
	//матрицей
	double arrMtrxTransf[4] ;
	memcpy(arrMtrxTransf, &constArrMtrxTransf[5 * numFreq +1], 4 * sizeof(double));
  MtrxMultMatrx(arrMtrxTransf, 2, 2, arrbT, 1, arrTransf_b) ;
 ///

 //
	double arrGenAngs[2] ={0.};
	int ireturn = 0;
 if ( valGenAngSdvig > VAlCritAngSdvig)
 { // эллиптический случай
	// нахождение углов точек касания прямых, проходящих через точку arrTransf_b
	//   единичной  окружности
	double valFi0 = atan2(arrTransf_b[1], arrTransf_b[0]);
	double valD = sqrt(arrTransf_b[0] * arrTransf_b[0] + arrTransf_b[1] * arrTransf_b[1]);
	if (valD < 1.)
	{
	  return 0;
	}
	double gam  =  acos(1. / valD);
	double arrt[2]= {0.};
	arrt[0] = valFi0  + gam;
	arrt[1] = valFi0  - gam;
	fncCorrectArg(arrt);
	fncCorrectArg(&arrt[1]);

	///

	//


	ireturn = 0;
	// нахождение номера строки в таблтице    ЭЛЛИПТИЧЕСКИЙ СЛУЧАЙ//  ТАБЛИЦА КОЭФФИЦИЕНТОВ ПОЛИНОМОВ СТЕПЕНИ 6
	int iNumRow1 = findNumRow(constArrCoeffPolinEllipticCase, QUantCasesEllipticType, 4 + LEnArrCoefPolinEllipticCase,  valFreq) ;
	for (int i =0; i < 2; i++)
	{
	  if (fabs(arrt[ireturn])> constArrCoeffPolinEllipticCase[iNumRow1 * (4 + LEnArrCoefPolinEllipticCase) +1])
	  {
		 continue;
	  }
		 double temp = fabs( arrt[i]);
		 double valPolinom = fncPolinom(&constArrCoeffPolinEllipticCase[iNumRow1 * (4 + LEnArrCoefPolinEllipticCase) + 4], LEnArrCoefPolinEllipticCase, temp);
		 if (fabs(valPolinom) > M_PI/2.)
		 {
		   int iii =0;
		 }
		 arrGenAngs [ireturn] = -tan(valPolinom)*constArrCoeffPolinEllipticCase[iNumRow1 * (4 + LEnArrCoefPolinEllipticCase) + 3];

	 if (arrt[i] >0.)
	 {
	   arrGenAngs[ireturn] = -arrGenAngs [ireturn];
	 }
	 ireturn++;

	}



 }
 else
 {
   // ГИперболический случай
   THyperbola hyperb( 1., 1.);
   TURPointXY  pntArrOut[2], pntB(arrTransf_b[0], arrTransf_b[1]);
   int irez = hyperb.findTangencyPoints(pntB,  pntArrOut[0],  pntArrOut[1]);
   if (irez == 0)
   {
	 return 0;
   }
	ireturn = 0;
	for (int i = 0; i < 2; i++)
	{
	  bool brez = findGenAng_UsingApproximation_ForHyperbolicCase(pntArrOut[i], &arrGenAngs[ireturn]);
	  if(brez)
	  {
		ireturn++ ;
	  }
	}


 }



   if (ireturn ==2)
   {
	 *pvalGenTargEps =  arrGenAngs[0];
	*pvalGenAntpEps =  arrGenAngs[1];
	if ((*pvalGenTargEps) < (*pvalGenAntpEps) )
	 {
	 double temp = *pvalGenAntpEps;
	 (*pvalGenAntpEps) = (*pvalGenTargEps);
	 (*pvalGenTargEps) = temp;
	 }
   }

   if (ireturn == 1)
   {
        *pvalGenTargEps =  arrGenAngs[0];
   }

   return ireturn ;
}


///


bool TParAnt::findGenAng_UsingApproximation_ForHyperbolicCase(TURPointXY pntInp, double *pGenAng)
{
	double valFreq= VAL_C * 100./ mLambda/ 1000000000.;
	double tang = pntInp.Y/pntInp.X;
	double valt = log((1.+ tang)/(1. - tang))/ 2.;
	if (pntInp.X <=0.)
	{// левая ветвь гиперболы

	 int iNumRow0 = findNumRow(constArrCoeffPolinHypCaseLeft, QUantCasesHyperbolicType, LEnArrCoefPolinHypCaseLeft+1,  valFreq) ;
	 double *arrCoeff = new double  [LEnArrCoefPolinHypCaseLeft];
	 memcpy(arrCoeff, &constArrCoeffPolinHypCaseLeft[(LEnArrCoefPolinHypCaseLeft+1)*iNumRow0 +1], LEnArrCoefPolinHypCaseLeft * sizeof(double));
	 *pGenAng = fncOddPolinom(arrCoeff, LEnArrCoefPolinHypCaseLeft, valt);
	 delete arrCoeff ;
	 return true;
	}
	else
	{ // правая ветвь гиперболы
	  int iNumRow0 = findNumRow(constArrLoranCoeffPolinHypCaseRight, QUantCasesHyperbolicType, LEnArrLoranCoeff + 1,  valFreq) ;
	   double *arrCoeff = new double  [LEnArrLoranCoeff];
	   memcpy(arrCoeff, &constArrLoranCoeffPolinHypCaseRight[(LEnArrLoranCoeff + 1)*iNumRow0 +1], LEnArrLoranCoeff * sizeof(double));

	   int iNumRow1= findNumRow(constArrAreaCoeffHypCaseRight, QUantCasesHyperbolicType,5,  valFreq) ;
	   double bound = constArrAreaCoeffHypCaseRight[iNumRow1 * 5 + 2];
	   if (valt < bound)
	   {
		*pGenAng = fncPolinom(arrCoeff, LEnArrLoranCoeff, valt) / valt/ valt;
		 delete arrCoeff ;
		 return true;
	   }
	   else
	   {
		   if (valt > -bound)
		   {
			 *pGenAng = -fncPolinom(arrCoeff, LEnArrLoranCoeff, -valt) / valt/ valt;
			 delete arrCoeff ;
			 return true;
		   }
		   else
		   {
			 delete arrCoeff ;
			 return false;
           }
       }
	 // fncPolinom(&arrCoeff[NPolinom +1], NPolinom +1, plnValue11.Points[i].X)/ plnValue11.Points[i].X/ plnValue11.Points[i].X;
	  delete arrCoeff ;
	}
}

///
int  findNumRow(const double *arr, int  numRows, int  numCols, double valFreq)
{
for (int i = 0; i < numRows; i++)
   {
	 if (fabs(valFreq - arr[i * numCols]) < 0.01)
	 {
	   return i;
	   break;
	 }
   }

	return -1;

}

void fncCorrectArg(double *t)
{
	 if ((*t) > M_PI)
	{
	  (*t)  -= 2. *M_PI;
	}

	if ((*t)  < -M_PI)
	{
	  (*t)  += 2. *M_PI;
	}

}

double fncPoinomApproximatedEstimationEllipticCase(double valArgTemp ,double valSlop
	 ,double *arrCoeff, int lenArrCoeff)
{
  if (valArgTemp >0)
	  {


			return tan(fncPolinom(arrCoeff, lenArrCoeff, valArgTemp))*valSlop;

	  }
	  else
	  {
			return -tan(fncPolinom(arrCoeff, lenArrCoeff, valArgTemp))*valSlop;
	  }
}



// вычисление коррел матрицы ошибок оценивания обобщенных углов цели и антипода
// альтернативным способом через измерения диаграмм
// INPUT:
// valGenTargEps, valGenAntpEps - оценки обобщенных углов цели и антипода относительно центральной диаграммы тройки
// OUTPUT:
// arrMtrxCorrGenAngs[4] - коррел матрица ошибок оценивания углов цели и антоипода
bool TParAnt::calcMtrxCorrGenAngs_AlternativeMeth3PartDiagr(TComp *cmparrS,  double valGenTargEps,  double valGenAntpEps
  ,TComp cmpKTarg ,TComp cmpKAntp ,const int NumAnsamble, double *arrMtrxCorrGenAngs )
{
	 double valGenAngSdvig = transformAngToGeneralizedAng (mAngSdvig ) ;

	 // проверки
	 double arr_b[2] ={0.};
	 calcVect_b(&cmparrS[NumAnsamble], arr_b);
	 double val_b0 = arr_b[0];
	 double val_b1 = arr_b[1];

  // формирование матирицы частных производных dB/dS
  // первая строка это градиент функции b0 по S
  // вторая строка это градиент функции b1 по S
  double arrdB_po_dS[12] = {0.};
  calc_dB0_po_dS(&cmparrS[NumAnsamble], valGenAngSdvig,  val_b0,  val_b1, arrdB_po_dS);
  calc_dB1_po_dS(&cmparrS[NumAnsamble], valGenAngSdvig,  val_b0,  val_b1, &arrdB_po_dS[6]);
  ///

  double arrGenAngs[2] = {0.};
  arrGenAngs[0] =  valGenTargEps;
  arrGenAngs[1] =  valGenAntpEps;
  ///

  // формирование матрицы частных производных углов по измерениям
  double arr_dMu_po_dS[12] = {0.};
  for (int i =0; i < 2; i++)
  {

	 double coeff = fncDerivDiagrSinx_div_x(arrGenAngs[i])
				   -  val_b0 *  fncDerivDiagrSinx_div_x(arrGenAngs[i] + valGenAngSdvig)
				   -  val_b1 *  fncDerivDiagrSinx_div_x(arrGenAngs[i] - valGenAngSdvig) ;
	 double arrF[2] = {0.};
	 arrF[0] = -fncDiagrSinx_div_x(arrGenAngs[i] + valGenAngSdvig) / coeff;
	 arrF[1] = -fncDiagrSinx_div_x(arrGenAngs[i] - valGenAngSdvig) / coeff;
	// MtrxTranspMultMatrx(arrdB_po_dS,2, 6, arrF, 1, &arr_dMu_po_dS[ i * 6]) ;
	 MtrxMultMatrx(arrF,1, 2, arrdB_po_dS, 6, &arr_dMu_po_dS[ i * 6]) ;
  }

  ///

  // формирование корреляц матрицы измерений
  double arrCorrPartDiagr[36] = {0.};

  createMtrxCorrForMeasuresForAnsambleOfPartDiagrs( valGenTargEps,  valGenAntpEps
  , cmpKTarg , cmpKAntp ,NumAnsamble, 3
  , arrCorrPartDiagr );
  ///

  // вычисление коррел матрицы вектора дельтах
  double  arrtemp0[12] = {0.};
  MtrxMultMatrx(arr_dMu_po_dS,2, 6, arrCorrPartDiagr, 6, arrtemp0) ;
  MtrxMultMatrxTransp(arrtemp0,2, 6, arr_dMu_po_dS,2, arrMtrxCorrGenAngs) ;
  ///
  return true;
}

void calc_dB0_po_dS(TComp *cmparrS, double valGenAngSdvig, double val_b0, double val_b1, double *arrdB_po_dS)
{
/* double s11 = cmparrS[0].m_Re;
 double s12 = cmparrS[0].m_Im;
 double s21 = cmparrS[1].m_Re;
 double s22 = cmparrS[1].m_Im;
 double s31 = cmparrS[2].m_Re;
 double s32 = cmparrS[2].m_Im;*/

 double temp = cmparrS[0].m_Re * cmparrS[2].m_Im - cmparrS[2].m_Re * cmparrS[0].m_Im;
 arrdB_po_dS[0] = -val_b0 * cmparrS[2].m_Im/ temp;
 arrdB_po_dS[1] =  val_b0 * cmparrS[2].m_Re/ temp;
 arrdB_po_dS[2] =  cmparrS[2].m_Im/ temp;
 arrdB_po_dS[3] =  -cmparrS[2].m_Re/ temp;
 arrdB_po_dS[4] =  -cmparrS[2].m_Im * val_b1/ temp;
 arrdB_po_dS[5] =  cmparrS[2].m_Re * val_b1/ temp;

}


void calc_dB1_po_dS(TComp *cmparrS, double valGenAngSdvig, double val_b0, double val_b1, double *arrdB_po_dS)
{
/* double s11 = cmparrS[0].m_Re;
 double s12 = cmparrS[0].m_Im;
 double s21 = cmparrS[1].m_Re;
 double s22 = cmparrS[1].m_Im;
 double s31 = cmparrS[2].m_Re;
 double s32 = cmparrS[2].m_Im; */

 double temp = cmparrS[0].m_Re * cmparrS[2].m_Im - cmparrS[2].m_Re * cmparrS[0].m_Im;
 arrdB_po_dS[0] = val_b0 * cmparrS[0].m_Im/ temp;
 arrdB_po_dS[1] =  -val_b0 * cmparrS[0].m_Re/ temp;
 arrdB_po_dS[2] = - cmparrS[0].m_Im/ temp;
 arrdB_po_dS[3] =  cmparrS[0].m_Re/ temp;
 arrdB_po_dS[4] =  cmparrS[0].m_Im * val_b1/ temp;
 arrdB_po_dS[5] =  -cmparrS[0].m_Re * val_b1/ temp;
}


// находждение угла сдвига, при котором происходит
// измерений типа огибающей с эллиптоического на гипнрболический
// wchrAlfaOverStitch!=NULL то в папку wchrAlfaOverStitch записываются шейп файлы
// с графиками
bool 	  TParAnt::findAlfaOverSwitch(wchar_t *Fold, double *pvalAlfa, double *pvalb, double *pvalMu)
{
   TParAnt  ParAnt;
   double valMuSingular = 0.,  valBSingular = 0.;

   ParAnt.findEnvelopeSigularPoint (&valMuSingular, &valBSingular);
   double valAlf0 = ParAnt.transformAngToGeneralizedAng ( ParAnt.mAngSdvig  ) ;

   TParAnt  ParAnt1;
   ParAnt1.mLambda = 3.5;
   double valAlfEnd = 1.;//ParAnt1.transformAngToGeneralizedAng ( ParAnt1.mAngSdvig  ) ;
   // метод Эйлера
	// шаг интегрированияф
	double valStep = 0.000001;
	double iNc = fabs(valAlfEnd -  valAlf0)/ valStep;
	valStep =  (valAlfEnd -  valAlf0)/ ((double)iNc);
	bool bChangeSign  = false;
	double valf0 = fncFF(valAlf0, valMuSingular);
	double alfaCur = valAlf0;
	double vmuCur =  valMuSingular;


	const int nBuffRows = iNc ;
	const int nBuffCols = 3;
	double  *parrBuff = new double [nBuffRows  * nBuffCols] ;
	memset(parrBuff, 0, nBuffRows * nBuffCols * sizeof(double));
	const int lenName = 30;
	wchar_t *wcharrFileNames  = new wchar_t [30 * nBuffCols];
	memset(wcharrFileNames, 0, 30 * nBuffCols * sizeof(wchar_t));

	wcscpy( wcharrFileNames, L"AlfSdviga");
	wcscpy( &wcharrFileNames[30], L"MuSingular");

	wcscpy( &wcharrFileNames[60], L"BSingular");
	parrBuff[0] =  valAlf0;
	parrBuff[1] =  valMuSingular;
	parrBuff[2] =  fncBSingular(valAlf0, valMuSingular);

	double *pscaley = new double [nBuffCols] ;




	for (int i =1; i < iNc; i++)
	{
	 double valf1 = 0.;
	 if (!bChangeSign)
	 {
	  valf1 = fncFF(alfaCur, valMuSingular);
	// if( valf1* valf0 <= 0.)
	 if( fabs( fncDiagrSinx_div_x(vmuCur + alfaCur) + fncDiagrSinx_div_x(vmuCur - alfaCur)) <= 0.000001)
	 {
	   bChangeSign = true;
	   *pvalAlfa =   alfaCur;
	   *pvalMu   = vmuCur;
	   *pvalb = fncDiagrSinx_div_x(vmuCur) / ( fncDiagrSinx_div_x(vmuCur + alfaCur) + fncDiagrSinx_div_x(vmuCur - alfaCur));
	 }

	 }
	 alfaCur +=  valStep;
	 vmuCur +=  valStep * valf1;
		if (Fold)
		{
		 parrBuff[ i * nBuffCols] =  alfaCur;
		 parrBuff[ i * nBuffCols + 1] = vmuCur;
		 parrBuff[ i * nBuffCols + 2] =  fncBSingular(alfaCur, vmuCur);

		}

   }




 if( Fold)
 {
		wchar_t Fold1 [400] ={0};
		wcscpy(  Fold1,  Fold);
		wcscat( Fold1, L"\\");
		pscaley[0] = 100.;
		pscaley[1] = 1.;
		pscaley[2] = 1.;
		// wchar_t wchFoldName[] = L"E:\\PROJECTS_C++\\НЛЦ-5П10-03\\Новая папка\\GraphFGr\\";
		for (int i=1 ; i < nBuffCols; i++)
		{

		TYrWriteShapeFile::WriteOneReport(Fold1  // путь к папке
								  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  , nBuffCols // - к-во переменных о корорых накоплена информация в буфере
								  , nBuffRows //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,lenName // максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i  // номер переменной по оси Y
								  ,1. //  масштаб по оси X
								  ,pscaley[i]// масштаб по оси Y
								   ) ;
		}




		wchar_t wchAxesFileName[300] ={0};
		wcscpy(  wchAxesFileName,  Fold);
		wcscat(wchAxesFileName, L"\\AxesArr.shp");
		TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-13., 13.5
		,-13.,13.5, -parrBuff[ 0] * 0.05) ;

		//////////////////////////////////////////////////

  }

		delete  parrBuff;
		delete wcharrFileNames;
		delete pscaley;
}
// правая часть дифференциальног уравнения
double TParAnt::fncFF( double valAlf,  double  valMu)
{
	double temp1 = fncDerivDiagrSinx_div_x(valMu) * (fncDerivDiagrSinx_div_x(valMu + valAlf) - fncDerivDiagrSinx_div_x(valMu - valAlf))
		  - fncDiagrSinx_div_x(valMu) * (fncDeriv2DiagrSinx_div_x(valMu + valAlf) - fncDeriv2DiagrSinx_div_x(valMu - valAlf));

	double temp2 = fncDeriv2DiagrSinx_div_x(valMu)* (fncDiagrSinx_div_x(valMu + valAlf) + fncDiagrSinx_div_x(valMu - valAlf))
		 -  fncDiagrSinx_div_x(valMu) *(fncDeriv2DiagrSinx_div_x(valMu + valAlf) + fncDeriv2DiagrSinx_div_x(valMu - valAlf));
	return -temp1 / temp2;
}

double TParAnt::fncBSingular( double valAlf0,  double  valMu)
{
	return fncDiagrSinx_div_x(valMu) / ( fncDiagrSinx_div_x(valMu +  valAlf0) + fncDiagrSinx_div_x(valMu -  valAlf0));
}
#pragma package(smart_init)


