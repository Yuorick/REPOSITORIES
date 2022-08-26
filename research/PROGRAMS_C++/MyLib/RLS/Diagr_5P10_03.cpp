//---------------------------------------------------------------------------


#pragma hdrstop

#include "Diagr_5P10_03.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "YrWriteShapeFile.h"
#include "URPolyLine.h"
#include "UrPointXY.h"
#include "DiagrSinX.h"
#include "Comp.h"
#include "Gauss.h"
#include "MatrixProccess.h"
#include <dir.h>
#include "URPolygon.h"
#include "Equations.h"

extern const double TET0707;
 // ширина парциальной диаграммы в обобщенном угле  по уровню 0,7
double findPartDiagrWidth07_GeneralizedAng()
{
	return  2.* TET0707;
}

// ширина парциальной диаграммы в обычном угле по уровню 0,7, рад
double findPartDiagrWidth07_Rad(const double VAlLamb,const double VAlApert)
{
	return  2.* asin(TET0707 * VAlLamb/ VAlApert/ M_PI);
}

// построение графиков 2-х парциальных диаграмм и суммарной в обобщенных угловых координатах и реальных
void   createGraphsPartial_and_Sum_Diagrams(wchar_t *wchFoldName1, const double valLamb
	, const double  valAppert , const double  valGenAngSdvig )
{

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
   parrBuff[ i * nBuffCols + 1] = asin ( tet * valLamb/ valAppert/ M_PI);
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
double transformAngToGeneralizedAng (const double valLamb , const double  valAppert
, const double  valAng  )
{
	return  valAppert * M_PI/ valLamb  * sin (valAng);
}

// пересчет обобщенного угла в обычный в рад
double transformGeneralizedAngToAng (const double valLamb , const double  valAppert
, const double  GeneralizedAng  )
{
	return  asin ( GeneralizedAng  * valLamb/ valAppert/ M_PI);
}

// нахождение обощенного угла при котором амплитуда диаграммы равна sqrt(2)/2
double findTet07_For_SumDiagr_5P10_03( const double  valGenAngSdvig)
{
 double tet0 = 1.8;//M_PI/2;

 int i =0;
 double valTemp = fncSumDiagr_5P10_03(  valGenAngSdvig, 0.);
 double a = sqrt(2.)/2. * valTemp;
 for ( i =0; i < 10; i++)
 {
  double del = -(fncSumDiagr_5P10_03(valGenAngSdvig,tet0) - a)/ fncDerivSumDiagr_5P10_03(valGenAngSdvig,tet0);
  tet0 += del;
  if (fabs(del) < 0.0000001) break;

 }
 return tet0;
}

// нахождене первого нуля суммарной диаграммы  обобщ угол
double findFirstZero_For_SumDiagr_5P10_03( const double  valGenAngSdvig)
{
	double valMu = 4.8; // первое пнриближение
	for (int i = 0; i < 100; i++)
	{
	  double valMu1 = valMu - fncSumDiagr_5P10_03(valGenAngSdvig, valMu)/fncDerivSumDiagr_5P10_03(valGenAngSdvig, valMu) ;
	  if (fabs(valMu1 - valMu)< 0.001)
	  {
		return  valMu1;
	  }
	  valMu = valMu1;
	}
	return -10000.0;
}

// нахождене второго нуля суммарной диаграммы  обобщ угол
double findSecondZero_For_SumDiagr_5P10_03( const double  valGenAngSdvig)
{
	double valMu = 8.1; // первое пнриближение
	for (int i = 0; i < 100; i++)
	{
	  double valMu1 = valMu - fncSumDiagr_5P10_03(valGenAngSdvig, valMu)/fncDerivSumDiagr_5P10_03(valGenAngSdvig, valMu) ;
	  if (fabs(valMu1 - valMu)< 0.001)
	  {
		return  valMu1;
	  }
	  valMu = valMu1;
	}
	return -10000.0;
}
double fncSumDiagr_5P10_03( const double  valGenAngSdvig, double tetGeneralized)
{
return fncDiagrSinx_div_x( tetGeneralized - valGenAngSdvig/2.) + fncDiagrSinx_div_x( tetGeneralized + valGenAngSdvig/2.);

}

double fncDerivSumDiagr_5P10_03( const double  valGenAngSdvig, double tetGeneralized)
{
return fncDerivDiagrSinx_div_x( tetGeneralized - valGenAngSdvig/2.) + fncDerivDiagrSinx_div_x( tetGeneralized + valGenAngSdvig/2.);
}


// ширина суммарной диаграммы в обобщенном угле  по уровню 0,7
double findSumDiagrWidth07_GeneralizedAng(const double  valGenAngSdvig)
{
  return 2. *	findTet07_For_SumDiagr_5P10_03(  valGenAngSdvig) ;

}

// ширина  суммарной  диаграммы в обычном угле по уровню 0,7, рад
double findSumDiagrWidth07_Rad(const double  valGenAngSdvig, const double VAlLamb,const double VAlApert)
{
	double valTet07 = findTet07_For_SumDiagr_5P10_03(  valGenAngSdvig) ;
	return  2.* asin(valTet07 * VAlLamb/ VAlApert/ M_PI);
}

// нахождение уровня пересечения двух суммарных сдвинутых диаграмм
double findCrossLevel_For_SumDiagr_5P10_03( const double  valGenAngSdvig)
{
 return  fncSumDiagr_5P10_03(valGenAngSdvig,valGenAngSdvig/2.)/ fncSumDiagr_5P10_03(valGenAngSdvig,0);
}

// имитация замеров 7 первичных диаграмм
// INPUT:
// valLambda
// valAppert
// valAngSdvig  - угол сдвига между соседними диаграммами в рад
// valalfUMTrg
// valalfUMAntp
// cmpKAntp
// OUTPUT:
// cmparrPartS [7]
// cmparrPartSZv [7]
//
//
void ImitateMeasureArrayPartialDiagrams_5P10_03( double valLamb,double valAppert , double valAngSdvig
		, double valalfUMTrg, TComp cmpKTarg, double valalfUMAntp,TComp cmpKAntp, const double VAlNoiseSkz
		,double valAmplFactSig, TComp *cmparrPartS, TComp *cmparrPartSZv)
{
   double valDispNoise = VAlNoiseSkz * VAlNoiseSkz;
   for (int i =0; i < 7; i++)
   {
	  // вычисление угла целей относителбно оси диаграммы с номером i в обобщенных координатах
	  double valDiagrAng =  ((double)(-3 + i))*  valAngSdvig; // угол оси диагр с номером i
	  double valTrg =  valalfUMTrg - valDiagrAng; // угол цели в диаграмме с номером i, рад.
	  double valTrgGen = transformAngToGeneralizedAng ( valLamb , valAppert, valTrg ) ; // угол цели в диаграмме с номером i, обощ угол
	  double valAntp =  valalfUMAntp - valDiagrAng; // угол цели в диаграмме с номером i, рад.
	  double valAntpGen = transformAngToGeneralizedAng ( valLamb , valAppert, valAntp ) ; // угол цели в диаграмме с номером i, обощ угол
	  ///
	  cmparrPartS[i] = (cmpKTarg * TComp(fncDiagrSinx_div_x(valTrgGen), 0.)) + ( cmpKAntp * TComp(fncDiagrSinx_div_x(valAntpGen), 0.) );

	  // диаграмма с направления цели
	  double valDiagrTarg = fncDiagrSinx_div_x(valTrgGen);

	  // диаграмма с навления антипода
	  double valDiagrAntp = fncDiagrSinx_div_x(valAntpGen);

	  // дисперсия действительной части сигнала
	   double valDispRe = ( (valDiagrTarg * cmpKTarg.m_Re * valDiagrTarg * cmpKTarg.m_Re)
						  + (valDiagrAntp * cmpKAntp.m_Re * valDiagrAntp * cmpKAntp.m_Re) )
						   * valAmplFactSig * valAmplFactSig+ valDispNoise;
	   // дисперсия мнипмой части сигнала
	   double valDispIm = ( (valDiagrTarg * cmpKTarg.m_Im * valDiagrTarg * cmpKTarg.m_Im)
						  + (valDiagrAntp * cmpKAntp.m_Im * valDiagrAntp * cmpKAntp.m_Im) )
						   * valAmplFactSig * valAmplFactSig + valDispNoise;


	  cmparrPartSZv[i].m_Re = cmparrPartS[i].m_Re  +  getGauss(0., sqrt(valDispRe));
	  cmparrPartSZv[i].m_Im = cmparrPartS[i].m_Im +  getGauss(0., sqrt(valDispIm));
   }
}

// имитация замеров 6 суммарных  диаграмм
void ImitateMeasureArraySumDiagrams_5P10_03(  TComp *cmparrPartS, TComp *cmparrPartSZv,  TComp *cmparrSumS, TComp *cmparrSumSZv)
{
  for (int i =0; i < 6; i++)
  {
	cmparrSumS[i] =  cmparrPartS[i] + cmparrPartS[i + 1];
	cmparrSumSZv[i] =  cmparrPartSZv[i] + cmparrPartSZv[i + 1];
  }
}

// оценивание обобщенных углов цели и антипода по измерениям 3 суммарных диаграмм
//INPUT:
// wchFoldName1 - путь к папаке с графиками
// cmparrS[3] -  массив замеров
// valGenAngSdvig -  обобщенный угол сдвига диаграмм
//  OUTPUT:
//  *pvalGenTargEps -   угол цели, обобщ коорд
// *pvalGenAntpEps -    угол антипода, обобщ коорд
int EstGenAngsThreeSumDiagr_5P10_03(TComp *cmparrS, double valGenAngSdvig
  ,  double  *pvalGenTargEps,  double  *pvalGenAntpEps, double *pval_b0, double *pval_b1)//, TComp *pcmpKTarg, TComp *pcmpKAntp, double *arrCrlMtrx)
{

  // вычисленние вектора b
  /* double arrS[4] = {0.};
   arrS[0] = cmparrS[0].m_Re;
   arrS[1] = cmparrS[2].m_Re;
   arrS[2] = cmparrS[0].m_Im;
   arrS[3] = cmparrS[2].m_Im;
   double arrs[2] = {0.};
   arrs[0] = cmparrS[1].m_Re;
   arrs[1] = cmparrS[1].m_Im;

   double arr_b[2] = {0.}, arrSInv[4] = {0.};
   bool brez =   InverseMtrx2(arrS, arrSInv);
   MtrxMultMatrx(arrSInv,2, 2, arrs,1, arr_b) ; */
   double arr_b[2] = {0.}  ;
   bool brez = calcVect_b(cmparrS,arr_b) ;
   *pval_b0 = arr_b[0];
   *pval_b1 = arr_b[1];
  // double root1 = findFirstZero_For_SumDiagr_5P10_03( valGenAngSdvig);
   //double root2 = findSecondZero_For_SumDiagr_5P10_03( valGenAngSdvig);

   double arrRoots[2] ={0.};
   int iNUmRoots = findRootsFgr_For_3SumDiagr_5P10_03(  valGenAngSdvig, arr_b,  arrRoots);
   *pvalGenTargEps =  arrRoots[1];
   *pvalGenAntpEps =  arrRoots[0];
   return iNUmRoots ;
}

// нахождение корней функции FGr длшя метода 3 суммарных диаграмм
// возвращает к-во корней
 int findRootsFgr_For_3SumDiagr_5P10_03(  double valGenAngSdvig, double *arr_b,  double *arrRoots)
 {

   int iNUmRoots = 0;

   double valMu0 =  findFirstZero_For_SumDiagr_5P10_03(valGenAngSdvig) + valGenAngSdvig - 0.1;
   double valStep = 0.25;
	const int lenBuff = 2. *valMu0/ valStep ;
   double *parrBuff  = new double [lenBuff];

  double valmu = -valMu0;

  for (int i=0 ; i < lenBuff; i++)
  {
   valmu = -valMu0 + ((double)i) * valStep;
   parrBuff[ i ] = fncFGr_SumDiagr(arr_b, valGenAngSdvig, valmu);

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
	 arrRoots [iNUmRoots ] = findRootMethChord_For_Fgr_For_3SumDiagr_5P10_03(  valGenAngSdvig,arr_b,  valX0, valX1) ;
	 iNUmRoots++;


   }
  }
	delete parrBuff;
  return iNUmRoots;

 }

double 	  findRootMethChord_For_Fgr_For_3SumDiagr_5P10_03(  double valGenAngSdvig
	 , double *arr_b, double  valX0, double valX1)
{
  double valMu0 = valX0 - 0.00001;
  double valMu1 = valX1 + 0.00001;
  double valMuRez =  valMu0 ;
  double eps = 0.001;
  for (int i = 0; i < 100; i++)
  {
	double f0 = fncFGr_SumDiagr(arr_b,  valGenAngSdvig, valMu0 ) ;
	double f1 = fncFGr_SumDiagr(arr_b,  valGenAngSdvig, valMu1 ) ;
	double valMuRez0 = valMu0 - f0 * ( valMu1 - valMu0)/ (f1 - f0);
	if (fabs(valMuRez0 - valMuRez) < eps)
	{
	 return valMuRez0;
	}
	double f2 = fncFGr_SumDiagr(arr_b,  valGenAngSdvig, valMuRez0 ) ;
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
double fncFGr_SumDiagr(double *arr_b, double valGenAngSdvig, double valmu)
{
  return (fncSumDiagr_5P10_03(valGenAngSdvig, valmu) - arr_b[0] * fncSumDiagr_5P10_03(valGenAngSdvig, valmu +valGenAngSdvig)
		 - arr_b[1] * fncSumDiagr_5P10_03(valGenAngSdvig, valmu - valGenAngSdvig));
}



// оценивание обобщенных углов цели и антипода по измерениям 3 суммарных диаграмм
//INPUT:
// wchFoldName1 - путь к папаке с графиками
// cmparrS[6] -  массив замеров
// valAngSdvig -  обобщенный угол сдвига диаграмм
// NumRayTriple-  номер тройки рабочих лучей (нумерация лучей идет снизу начиная с нуля)
//                номер тройки это порядковый номер нижнего луча
//  OUTPUT:
//  *valEstAngTarg -   угол цели, рад
// *valEstAngAntp -    угол антипода, рад
// cmpKTarg - коэф отражения цели
// cmpKAntp - коэф отражения антипода
// arrMtrxCorr  - коррел матрица ошибок измерения угла цели и антипода
int estimateMethThreeSumDiagr_5P10_03(wchar_t *wchFoldName1, double valNoiseSkz
  , double valmAmplFactSig, double valLambda,double valAppert , TComp *cmparrS
, int iNumRayTriple, double valGenAngSdvig
  , double *valEstAngTarg, double *valEstAngAntp, TComp *cmpKTarg , TComp *cmpKAntp , double *arrMtrxCorr
  , double *pval_b0, double *pval_b1 )
{
	// нахождение обобщенных углов цели и антипода
	double  valGenTargEps = -100.,  valGenAntpEps = -100.0;  // обобщ углы цели и антипода
	if (EstGenAngsThreeSumDiagr_5P10_03(  &cmparrS[iNumRayTriple], valGenAngSdvig
	,  &valGenTargEps,  &valGenAntpEps, pval_b0, pval_b1) !=2)
	{
	   return -1;
	}
	//
	// вычисление коэффиц отражения
	TComp cmparrA[4],cmparrAInv[4], cmparrK[2], cmparrSTemp[2];

	cmparrA[0] = TComp(fncSumDiagr_5P10_03(valGenAngSdvig, valGenTargEps + valGenAngSdvig), 0.);
	cmparrA[1] = TComp(fncSumDiagr_5P10_03(valGenAngSdvig, valGenAntpEps + valGenAngSdvig), 0.);
	cmparrA[2] = TComp(fncSumDiagr_5P10_03(valGenAngSdvig, valGenTargEps - valGenAngSdvig), 0.);
	cmparrA[3] = TComp(fncSumDiagr_5P10_03(valGenAngSdvig, valGenAntpEps - valGenAngSdvig), 0.);

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
			double valRadCentreRayPos =  ( -1.5 + ((double)iNumRayTriple))
				* transformGeneralizedAngToAng (valLambda , valAppert , valGenAngSdvig ) ;
		  // вычисление угла цели в радианах в АСК
			double valRadTargEps =  transformGeneralizedAngToAng (valLambda , valAppert , valGenTargEps ) ;// в радианы
			*valEstAngTarg = valRadTargEps + valRadCentreRayPos;
		  // вычисление угла антипода в радианах в АПСК
			double valRadAntpАEps =  transformGeneralizedAngToAng (valLambda , valAppert , valGenAntpEps ) ;// в радианы
			*valEstAngAntp = valRadAntpАEps + valRadCentreRayPos;
	///

	// вычисление коррел матрицы ошибок определения обобщенных углов цели и антипода
	double arrMtrxCorrMu[4] = {0.};
	calcMtrxCorrGenAngs_Meth3SumDiagr(valGenTargEps,  valGenAntpEps
	, valGenAngSdvig,*cmpKTarg,*cmpKAntp,valNoiseSkz, valmAmplFactSig, arrMtrxCorrMu) ;
	///

	// вычисление корреляц матрицы ошибок определения углов цели и антипода в радианаах
	double arrQ[4] = {0.}, arrTemp[4] = {0.};
	double coef = valLambda / valAppert / M_PI;
	arrQ[0] =  coef / sqrt( 1.- coef * valGenTargEps * coef * valGenTargEps);
	arrQ[3] =  coef / sqrt( 1.- coef * valGenAntpEps * coef * valGenAntpEps);
	MtrxMultMatrx(arrQ,2, 2, arrMtrxCorrMu,2, arrTemp)  ;
	MtrxMultMatrxTransp(arrTemp,2, 2, arrQ,2, arrMtrxCorr) ;
	return 1;

}

bool calcMtrxCorrGenAngs_Meth3SumDiagr(double valGenTargEps,  double valGenAntpEps , double valGenAngSdvig
  ,TComp cmpKTarg ,TComp cmpKAntp , double valNoiseSkz,double valmAmplFactSig, double *arrMtrxCorrGenAngs )
{
  // формирование матирицы частных производных dG/dx
  double arrdG_po_dX[36] = {0.};
  arrdG_po_dX[0] =  cmpKTarg.m_Re * fncDerivSumDiagr_5P10_03( valGenAngSdvig, valGenAngSdvig + valGenTargEps);
  arrdG_po_dX[1] =  cmpKAntp.m_Re * fncDerivSumDiagr_5P10_03( valGenAngSdvig, valGenAngSdvig + valGenAntpEps);
  arrdG_po_dX[2] =  fncSumDiagr_5P10_03( valGenAngSdvig, valGenAngSdvig + valGenTargEps);
  arrdG_po_dX[4] =  fncSumDiagr_5P10_03( valGenAngSdvig, valGenAngSdvig + valGenAntpEps);

  arrdG_po_dX[6] =  cmpKTarg.m_Im * fncDerivSumDiagr_5P10_03( valGenAngSdvig, valGenAngSdvig + valGenTargEps);
  arrdG_po_dX[7] =  cmpKAntp.m_Im * fncDerivSumDiagr_5P10_03( valGenAngSdvig, valGenAngSdvig + valGenAntpEps);
  arrdG_po_dX[9] =  fncSumDiagr_5P10_03( valGenAngSdvig, valGenAngSdvig + valGenTargEps);
  arrdG_po_dX[11] =  fncSumDiagr_5P10_03( valGenAngSdvig, valGenAngSdvig + valGenAntpEps);

  arrdG_po_dX[12] =  cmpKTarg.m_Re * fncDerivSumDiagr_5P10_03( valGenAngSdvig,  valGenTargEps);
  arrdG_po_dX[13] =  cmpKAntp.m_Re * fncDerivSumDiagr_5P10_03( valGenAngSdvig,  valGenAntpEps);
  arrdG_po_dX[14] =  fncSumDiagr_5P10_03( valGenAngSdvig,  valGenTargEps);
  arrdG_po_dX[16] =  fncSumDiagr_5P10_03( valGenAngSdvig,  valGenAntpEps);

  arrdG_po_dX[18] =  cmpKTarg.m_Im * fncDerivSumDiagr_5P10_03( valGenAngSdvig,  valGenTargEps);
  arrdG_po_dX[19] =  cmpKAntp.m_Im * fncDerivSumDiagr_5P10_03( valGenAngSdvig,  valGenAntpEps);
  arrdG_po_dX[21] =  fncSumDiagr_5P10_03( valGenAngSdvig,  valGenTargEps);
  arrdG_po_dX[23] =  fncSumDiagr_5P10_03( valGenAngSdvig,  valGenAntpEps);

  arrdG_po_dX[24] =  cmpKTarg.m_Re * fncDerivSumDiagr_5P10_03( valGenAngSdvig, -valGenAngSdvig + valGenTargEps);
  arrdG_po_dX[25] =  cmpKAntp.m_Re * fncDerivSumDiagr_5P10_03( valGenAngSdvig, -valGenAngSdvig + valGenAntpEps);
  arrdG_po_dX[26] =  fncSumDiagr_5P10_03( valGenAngSdvig, -valGenAngSdvig + valGenTargEps);
  arrdG_po_dX[28] =  fncSumDiagr_5P10_03( valGenAngSdvig, -valGenAngSdvig + valGenAntpEps);

  arrdG_po_dX[30] =  cmpKTarg.m_Im * fncDerivSumDiagr_5P10_03( valGenAngSdvig, -valGenAngSdvig + valGenTargEps);
  arrdG_po_dX[31] =  cmpKAntp.m_Im * fncDerivSumDiagr_5P10_03( valGenAngSdvig, -valGenAngSdvig + valGenAntpEps);
  arrdG_po_dX[33] =  fncSumDiagr_5P10_03( valGenAngSdvig, -valGenAngSdvig + valGenTargEps);
  arrdG_po_dX[35] =  fncSumDiagr_5P10_03( valGenAngSdvig, -valGenAngSdvig + valGenAntpEps);
  ///

  // формирование корреляц матрицы измерений
  double arrCorrPartDiagr[64] = {0.} // коррел матрица измерений 4 парциальных диаграмм
  , arrLinTransf[48] ={0.}, arrMtrxCorrSumDiagr[36] ={0.}, arrTEmp[48] = {0.};
  int iNumAnsamble = 1;
  int lenAnsamble = 4;
  createMtrxCorrForMeasuresForAnsambleOfPartDiagrs( valGenTargEps + valGenAngSdvig/2.,   valGenAntpEps + valGenAngSdvig/2.
   , valGenAngSdvig, cmpKTarg , cmpKAntp ,  valNoiseSkz,  valmAmplFactSig,  iNumAnsamble,  lenAnsamble  , arrCorrPartDiagr);


  for (int i =0; i < 6; i++)
  {
	arrLinTransf [ 8 * i + i] = 1.;
	arrLinTransf [ 8 * i + i + 2] = 1.;
  }

  MtrxMultMatrx(arrLinTransf,6, 8, arrCorrPartDiagr,8, arrTEmp) ;
  MtrxMultMatrxTransp(arrTEmp,6, 8, arrLinTransf,6, arrMtrxCorrSumDiagr) ;
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

// INPUT
//  valScaleY  - растяжение пол оси Y
//
//
//
//
//


void createPtesentGraphs_for_SumDiagr( wchar_t *wchrPresntSumDiagrams, double valLambda
	   , double valAppert , double valAngSdvig, int iNumRayTriple
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

   for (int i =0; i < 6; i++)
   {
	wchar_t FileName[400] ={0.};
	wcscpy( FileName, Fold);
	wcscat(FileName, L"\\SumDiagrNo_");
	wchar_t string[10] = {0.};
	_itow(i,string, 10);
	wcscat(FileName, string);
	wcscat(FileName, L".shp");
	double valSmeshenieGoriz =  ( -2.5 + ((double) i)) * valAngSdvig;
	double valSMeshenieVert = 0.;
	createGraphSumDiagr_from_rad(FileName, valLambda
	, valAppert , valAngSdvig, valSmeshenieGoriz, valSMeshenieVert, valScaleY);
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
	double valSmeshenieGoriz = ( -2.5 + ((double)(iNumRayTriple + i))) * valAngSdvig;

	createGraphSumDiagr_from_rad(FileName, valLambda
	, valAppert , valAngSdvig, valSmeshenieGoriz, 0., valScaleY ) ;
   }
   ///

   // построение графика функции FGreece в зависимости от угла рад
   wchar_t Fold1[400] ={0.};
   wcscpy( Fold1, wchrPresntSumDiagrams);
   wcscat(Fold1, L"\\GraphFGreece_from_rad_For_3SumDiagr");
	_wmkdir(Fold1);
	wcscat(Fold1, L"\\");

   // относиьтельный угол цетрального луча тройки в радианах
   double valRadCentreRayPos =  ( -1.5 + ((double)iNumRayTriple)) * valAngSdvig;

   createGraphFGreece_from_rad_For_3SumDiagr(Fold1,  valLambda
	, valAppert ,valAngSdvig,valRadCentreRayPos , valSMeshenieVert, val_b0, val_b1)  ;
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

void createGraphSumDiagr_from_rad(wchar_t *FileName, double valLambda
	,double  valAppert ,double  valAngSdvig,double  valSmeshenieGoriz, double  valSMeshenieVert,  double valScaleY )
{

	double step =0.0001;
	double  valGenAngSdvig = transformAngToGeneralizedAng (valLambda ,  valAppert, valAngSdvig  ) ;
	double valZero =  findFirstZero_For_SumDiagr_5P10_03(   valGenAngSdvig);
	const double VAlDiap = 1.5 *valZero;  // диапахзон графиков по обобщ координате плюс минус
	const int nBuffRows = 2. *VAlDiap/ step ;
	 TURPolyLine pln( 1, nBuffRows) ;



  double tet = -VAlDiap;
  for (int i=0 ; i < nBuffRows; i++)
  {
   tet = -VAlDiap + ((double)i) * step;
   pln.Points[i].X = (transformGeneralizedAngToAng (valLambda ,  valAppert, tet  ) + valSmeshenieGoriz)* 100. ;
   pln.Points[i].Y =fncSumDiagr_5P10_03( valGenAngSdvig, tet) * valScaleY + valSMeshenieVert;

  }

 pln.WriteSetSHPFiles(FileName, &pln, 1);
}

void createGraphFGreece_from_rad_For_3SumDiagr(wchar_t *Fold, double valLambda
	,double  valAppert ,double  valAngSdvig,double  valSmeshenieGoriz
	, double  valSMeshenieVert, double val_b0, double val_b1 )
{
	double arr_b[2] = {0.};
	arr_b[0] =  val_b0 ;
	arr_b[1]  = val_b1;
   // построение графика функции FGreece
   double valGenAngSdvig = transformAngToGeneralizedAng (valLambda , valAppert, valAngSdvig ) ;
   double valMu0 =  findFirstZero_For_SumDiagr_5P10_03(valGenAngSdvig) + valGenAngSdvig - 0.1;
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
   parrBuff[ i * nBuffCols] = (transformGeneralizedAngToAng (valLambda , valAppert, valmu ) + valSmeshenieGoriz) * 100.;
   parrBuff[ i * nBuffCols + 1] = fncFGr_SumDiagr(arr_b, valGenAngSdvig, valmu) + valSMeshenieVert;

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
}


// построение графика функции F греческое для метода 3 суммарных диаграмм в зависимости от обобщенного угла
void createGraphFGreece_For_3SumDiagr_from_GenAng(wchar_t *Fold, double valLambda
	,double  valAppert ,double  valGenAngSdvig,double  valSmeshenieGoriz
	, double  valSMeshenieVert, double val_b0, double val_b1 )
{
	double arr_b[2] = {0.};
	arr_b[0] =  val_b0 ;
	arr_b[1]  = val_b1;
   // построение графика функции FGreece

   double valMu0 =  findFirstZero_For_SumDiagr_5P10_03(valGenAngSdvig) + valGenAngSdvig - 0.1;
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
   parrBuff[ i * nBuffCols + 1] = fncFGr_SumDiagr(arr_b, valGenAngSdvig, valmu) + valSMeshenieVert;

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
}


// построение графико СКЗ оши бок оценивания в звисимости от углового расстояния
// INPUT:
// Fold
// valLambda - раб частота
// valAppert - апертура антиенны
// valAngSdvig - угол сдвига между диаграммами в рад
// valSmeshenieGoriz  - смещение графика по горизонтали
// valSMeshenieVert   - смещение графика по вертикали
// cmpKTarg  - комплексный	 сигнал цели
// cmpKAntp  - комплексный	 сигнал антипода
// valNoiseSkz - скз шума в парциальной диаграмме
// valmAmplFactSig - скз шума коэффиц усиления диаграммы
// INPUT:
// графики зависимости ско ошибок оценивания цели и антипода в зависимости
// от углового расстояния между целью и антиподом
// цельнаходится в центре центрального луча
// антипод приближается к цели
void createGraphsSKZ_from_AngDiffer_For_3SumDiagr(wchar_t *Fold,double valLambda
	,double  valAppert ,double  valAngSdvig,double  valSmeshenieGoriz
	, double  valSMeshenieVert, TComp cmpKTarg, TComp cmpKAntp , double valNoiseSkz,double valmAmplFactSig)
{

   double valGenAngSdvig = transformAngToGeneralizedAng (valLambda , valAppert, valAngSdvig ) ;
  // double valMu0 =  3. * findFirstZero_For_SumDiagr_5P10_03(valGenAngSdvig) + valGenAngSdvig ;
   double valStep = 0.001;
	const int nBuffRows = 3. * findFirstZero_For_SumDiagr_5P10_03(valGenAngSdvig)/ valStep ;
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

  double valRadAngTarg = valAngSdvig /2.;
  double valGenAngTarg =  transformAngToGeneralizedAng (valLambda , valAppert, valRadAngTarg ) ;
  double valGenAngAntp = 0.;
  double   arrMtrxCorr[4] ={0.}, arrMtrxCorrMu[4] = {0.};

  for (int i=0 ; i < nBuffRows; i++)
  {
   valGenAngAntp = valGenAngTarg - ((double)i + 10.) * valStep;


   	// вычисление коррел матрицы ошибок определения обобщенных углов цели и антипода
	calcMtrxCorrGenAngs_Meth3SumDiagr( valGenAngTarg ,  valGenAngAntp , valGenAngSdvig
	  , cmpKTarg, cmpKAntp,valNoiseSkz,  valmAmplFactSig, arrMtrxCorrMu) ;
	///

	// вычисление корреляц матрицы ошибок определения углов цели и антипода в радианаах
	double arrQ[4] = {0.}, arrTemp[4] = {0.};
	double coef = valLambda / valAppert / M_PI;
	arrQ[0] =  coef / sqrt( 1.- coef * valGenAngTarg * coef * valGenAngTarg);
	arrQ[3] =  coef / sqrt( 1.- coef *  valGenAngAntp * coef *  valGenAngAntp);
	MtrxMultMatrx(arrQ,2, 2, arrMtrxCorrMu,2, arrTemp)  ;
	MtrxMultMatrxTransp(arrTemp,2, 2, arrQ,2, arrMtrxCorr) ;


   parrBuff[ i * nBuffCols] = valRadAngTarg - transformGeneralizedAngToAng (valLambda , valAppert,valGenAngAntp);

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
	 createGraphSumDiagr_from_rad(wchAxesFileName,  valLambda
	,  valAppert ,  valAngSdvig,  valSmeshenieGoriz,   valSMeshenieVert,  4. );

		wcscpy(  wchAxesFileName,  Fold);
	wcscat(wchAxesFileName, L"DigarSum0.shp");
	 createGraphSumDiagr_from_rad(wchAxesFileName,  valLambda
	,  valAppert ,  valAngSdvig,  -valAngSdvig,   valSMeshenieVert,  4. );

		wcscpy(  wchAxesFileName,  Fold);
	wcscat(wchAxesFileName, L"DigarSum2.shp");
	 createGraphSumDiagr_from_rad(wchAxesFileName,  valLambda
	,  valAppert ,  valAngSdvig,  valAngSdvig,   valSMeshenieVert,  4. );

	// полигон изображающий истинная цель
	TURPointXY  pntTopRight2(0.05, cmpKTarg.modul()/ 800. * 5.);
	TURPointXY  pntBottomLeft2(- 0.05, 0.);
	TURPolygon plgTargTrue(pntTopRight2, pntBottomLeft2);
	wcscpy( wchAxesFileName, Fold);
	wcscat(wchAxesFileName, L"\\plgTargTrue.shp");
	plgTargTrue.WriteSetSHPFiles(wchAxesFileName ,&plgTargTrue, 1) ;
}

// графики для сравнения суммарной диаграммы м ее аппроксимации полинрмом 4 степени
void createGraphs_for_Compare_Diagrams(wchar_t *Fold , double valLambda
	,double  valAppert ,double  valGenAngSdvig  )
{

	double valStep =0.01;

	double valZero =  findFirstZero_For_SumDiagr_5P10_03(   valGenAngSdvig);
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

	double valF0 = fncSumDiagr_5P10_03( valGenAngSdvig, 0.) ;
	double root1 = findFirstZero_For_SumDiagr_5P10_03( valGenAngSdvig);
	double root2 = findSecondZero_For_SumDiagr_5P10_03( valGenAngSdvig);
   for (int i=0 ; i < nBuffRows; i++)
  {
   double valGenAngCur = VAlDiap - ((double)i ) * valStep;
   parrBuff[ i * nBuffCols] =  valGenAngCur;
   parrBuff[ i * nBuffCols +1] =fncSumDiagr_5P10_03( valGenAngSdvig, valGenAngCur)/valF0 ;
   parrBuff[ i * nBuffCols +2] = ( valGenAngCur * valGenAngCur -  root1 * root1)
		* (valGenAngCur * valGenAngCur -  root2 * root2)/ (root1 * root1 * root2 * root2);

  }

 // double scalex = 100.;
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
}

// графики для сравнения парциальной  диаграммы sinx/x и ее аппроксимации полинрмом 4,6 и 8 степени
void createGraphs_for_Compare_PartDiagrams(wchar_t *Fold  ,double  valGenAngSdvig  )
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
}


// вычисление вектора b
// INPUT:
// cmparrS[3] - массив измерений тройки диаграмм
// OUTPUT:
// arr_b[2] -  коэффициенты уравнения
bool calcVect_b(TComp *cmparrS, double *arr_b)
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
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//////  ВЫЧИСЛИТЕЛЬНЫЕ ФУНКЦИИ ДЛЯ МЕТОДА 3 ПАРЦИАЛЬНЫХ ДИАГРАММ  /////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// оценивание обобщенных углов цели и антипода по измерениям 3 ПАРЦИАЛЬНЫХ диаграмм
//INPUT:
// wchFoldName1 - путь к папаке с графиками
// cmparrS[6] -  массив замеров
// valAngSdvig -  обобщенный угол сдвига диаграмм
// NumRayTriple-  номер тройки рабочих лучей (нумерация лучей идет снизу начиная с нуля)
//                номер тройки это порядковый номер нижнего луча
//  OUTPUT:
//  *valEstAngTarg -   угол цели, рад
// *valEstAngAntp -    угол антипода, рад
// cmpKTarg - коэф отражения цели
// cmpKAntp - коэф отражения антипода
// arrMtrxCorr  - коррел матрица ошибок измерения угла цели и антипода
int estimateMethThreePartDiagr_5P10_03(wchar_t *wchFoldName1, double valNoiseSkz, double valmAmplFactSig
  , double valLambda,double valAppert , TComp *cmparrS
  , int iNumRayTriple, double valGenAngSdvig
  , double *valEstAngTarg, double *valEstAngAntp, TComp *cmpKTarg , TComp *cmpKAntp , double *arrMtrxCorr
  , double *pval_b0, double *pval_b1 )
{
	// нахождение обобщенных углов цели и антипода
	double  valGenTargEps = -100.,  valGenAntpEps = -100.0;  // обобщ углы цели и антипода
	if (EstGenAngsThreePartDiagr_5P10_03(  &cmparrS[iNumRayTriple], valGenAngSdvig
	,  &valGenTargEps,  &valGenAntpEps, pval_b0, pval_b1) !=2)
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
			double valRadCentreRayPos =  ( -2. + ((double)iNumRayTriple))
				* transformGeneralizedAngToAng (valLambda , valAppert , valGenAngSdvig ) ;
		  // вычисление угла цели в радианах в АСК
			double valRadTargEps =  transformGeneralizedAngToAng (valLambda , valAppert , valGenTargEps ) ;// в радианы
			*valEstAngTarg = valRadTargEps + valRadCentreRayPos;
		  // вычисление угла антипода в радианах в АПСК
			double valRadAntpАEps =  transformGeneralizedAngToAng (valLambda , valAppert , valGenAntpEps ) ;// в радианы
			*valEstAngAntp = valRadAntpАEps + valRadCentreRayPos;
	///

	// вычисление коррел матрицы ошибок определения обобщенных углов цели и антипода
	double arrMtrxCorrMu[4] = {0.};
	calcMtrxCorrGenAngs_Meth3PartDiagr(valGenTargEps,  valGenAntpEps , valGenAngSdvig
	  ,*cmpKTarg,*cmpKAntp,valNoiseSkz, valmAmplFactSig, arrMtrxCorrMu) ;
	///

	// вычисление корреляц матрицы ошибок определения углов цели и антипода в радианаах
	double arrQ[4] = {0.}, arrTemp[4] = {0.};
	double coef = valLambda / valAppert / M_PI;
	arrQ[0] =  coef / sqrt( 1.- coef * valGenTargEps * coef * valGenTargEps);
	arrQ[3] =  coef / sqrt( 1.- coef * valGenAntpEps * coef * valGenAntpEps);
	MtrxMultMatrx(arrQ,2, 2, arrMtrxCorrMu,2, arrTemp)  ;
	MtrxMultMatrxTransp(arrTemp,2, 2, arrQ,2, arrMtrxCorr) ;
	return 1;

}


// вычисление коррел матрицы ошибок оценивания обобщенных углов цели и антипода
// INPUT:
// valGenTargEps, valGenAntpEps - оценки обобщенных углов цели и антипода относительно центральной диаграммы тройки
// valGenAngSdvig  -обобщ уголл сдвига диаграмм
// cmpKTarg, cmpKAntp - оценки  комплексных сигналов цели и антипода
// valNoiseSkz - скз внутр шума антеныы
// valAmplFactSig  - скз шума коэффиц усиления (в долях)
// OUTPUT:
// arrMtrxCorrGenAngs[4] - коррел матрица ошибок оценивания углов цели и антоипода
bool calcMtrxCorrGenAngs_Meth3PartDiagr(double valGenTargEps,  double valGenAntpEps , double valGenAngSdvig
  ,TComp cmpKTarg ,TComp cmpKAntp , double valNoiseSkz, double valAmplFactSig, double *arrMtrxCorrGenAngs )
{
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
  int iNumAnsamble = 2;
  createMtrxCorrForMeasuresForAnsambleOfPartDiagrs( valGenTargEps,  valGenAntpEps ,  valGenAngSdvig
  , cmpKTarg , cmpKAntp , valNoiseSkz, valAmplFactSig,  iNumAnsamble, 3
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

//формирование корреляционной матрицы ошибок первичных измерений
// ансамбля первичных диаграмм
// вектор первичных ошибок имеет вид ...S[i].m_Re S , S[i].m_Im...
// INPUT:
// valGenTargEps  - обобщ угол цели
// valGenAntpEps - обобщ угол антипода
// valGenAngSdvig - обобщ угол сдвига диаграмм
// cmpKTarg, cmpKAntp - сигналы(комлексные) цели и антипода
// valNoiseSkz, - скз внутреннего шума антенны
// valmAmplFactSig - скз шума коэффиц усиления
//  iNumAnsamble   - намер первой диаграммы в ансамбле
//          (все диаграммы перенумерованы в порядке возрастания угла сдвига начиная с нулевого номера)
// lenAnsamble  - количество диаграмм в ансамбле
// OTPUT:
// arrMtrxCorr[2 *lenAnsamble *2 *lenAnsamble] - корреляц матрица
void createMtrxCorrForMeasuresForAnsambleOfPartDiagrs(double valGenTargEps,  double valGenAntpEps , double valGenAngSdvig
  ,TComp cmpKTarg ,TComp cmpKAntp , double valNoiseSkz, double valmAmplFactSig, int iNumAnsamble, int lenAnsamble
  , double *arrMtrxCorr )
 {

	memset(arrMtrxCorr, 0, 4 * lenAnsamble * lenAnsamble * sizeof(double));
	double valDispNoise = valNoiseSkz * valNoiseSkz;
  // вычисление коррел матрицы, обусловленной разбросом коэффиц усиления диаграмм
	for (int i = 0; i < lenAnsamble; i++)  // цикл по диаграммам в ансамбле
	{
	   // вычисление угла целей относителбно оси диаграммы с номером i в обобщенных координатах
	  double valDiagrGenAng =  ((double)(-3 + iNumAnsamble + i))*  valGenAngSdvig; // угол оси диагр с номером i
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
						   * valmAmplFactSig * valmAmplFactSig;
	   // дисперсия мнипмой части сигнала
	   double valDispIm = ( (valDiagrTarg * cmpKTarg.m_Im * valDiagrTarg * cmpKTarg.m_Im)
						  + (valDiagrAntp * cmpKAntp.m_Im * valDiagrAntp * cmpKAntp.m_Im) )
						   * valmAmplFactSig * valmAmplFactSig;
	 ///
	  // заполнение диагональных элементов матрицы arrMtrxCorr в соответствующих строках
	  // действительная и комплек5сная части сигнала диаграммы с номером i занимают строки 2*i и 2*i+1 соответственно
	  arrMtrxCorr[2*lenAnsamble * (2 * i) + 2 *i] =  valDispRe  + valDispNoise;
	  arrMtrxCorr[2*lenAnsamble * (2 * i + 1) + 2 *i + 1] =  valDispIm  + valDispNoise;

	}
	return;
 }

// оценивание обобщенных углов цели и антипода по измерениям 3 парциальных диаграмм
//INPUT:
// wchFoldName1 - путь к папаке с графиками
// cmparrS[3] -  массив замеров
// valGenAngSdvig -  обобщенный угол сдвига диаграмм
//  OUTPUT:
//  *pvalGenTargEps -   угол цели, обобщ коорд
// *pvalGenAntpEps -    угол антипода, обобщ коорд
int EstGenAngsThreePartDiagr_5P10_03(TComp *cmparrS, double valGenAngSdvig
  ,  double  *pvalGenTargEps,  double  *pvalGenAntpEps, double *pval_b0, double *pval_b1)//, TComp *pcmpKTarg, TComp *pcmpKAntp, double *arrCrlMtrx)
{

  // вычисленние вектора b
   double arr_b[2] = {0.}  ;
   bool brez = calcVect_b(cmparrS,arr_b) ;
   ///
   *pval_b0 = arr_b[0];
   *pval_b1 = arr_b[1];


   double arrRoots[2] ={0.};
   int iNUmRoots = findRootsFgr_For_3PartDiagr_5P10_03(  valGenAngSdvig, arr_b,  arrRoots);
   *pvalGenTargEps =  arrRoots[1];
   *pvalGenAntpEps =  arrRoots[0];
   return iNUmRoots ;
}








// нахождение корней функции FGr длшя метода 3 суммарных диаграмм
// возвращает к-во корней
 int findRootsFgr_For_3PartDiagr_5P10_03(  double valGenAngSdvig, double *arr_b,  double *arrRoots)
 {

   int iNUmRoots = 0;

   double valMu0 =  M_PI + valGenAngSdvig - 0.1;
   double valStep = 0.25;
	const int lenBuff = 2. *valMu0/ valStep ;
   double *parrBuff  = new double [lenBuff];

  double valmu = -valMu0;

  for (int i=0 ; i < lenBuff; i++)
  {
   valmu = -valMu0 + ((double)i) * valStep;
   parrBuff[ i ] = fncFGr_PartDiagr(arr_b, valGenAngSdvig, valmu);

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
	 arrRoots [iNUmRoots ] = findRootMethChord_For_Fgr_For_3PartDiagr_5P10_03(  valGenAngSdvig,arr_b,  valX0, valX1) ;
	 iNUmRoots++;


   }
  }
	delete parrBuff;
  return iNUmRoots;

 }

double 	  findRootMethChord_For_Fgr_For_3PartDiagr_5P10_03(  double valGenAngSdvig
	 , double *arr_b, double  valX0, double valX1)
{
  double valMu0 = valX0 - 0.00001;
  double valMu1 = valX1 + 0.00001;
  double valMuRez =  valMu0 ;
  double eps = 0.001;
  for (int i = 0; i < 100; i++)
  {
	double f0 = fncFGr_PartDiagr(arr_b,  valGenAngSdvig, valMu0 ) ;
	double f1 = fncFGr_PartDiagr(arr_b,  valGenAngSdvig, valMu1 ) ;
	double valMuRez0 = valMu0 - f0 * ( valMu1 - valMu0)/ (f1 - f0);
	if (fabs(valMuRez0 - valMuRez) < eps)
	{
	 return valMuRez0;
	}
	double f2 = fncFGr_PartDiagr(arr_b,  valGenAngSdvig, valMuRez0 ) ;
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
double fncFGr_PartDiagr(double *arr_b, double valGenAngSdvig, double valmu)
{
  return (fncDiagrSinx_div_x(valmu) - arr_b[0] * fncDiagrSinx_div_x(valmu +valGenAngSdvig)
		 - arr_b[1] * fncDiagrSinx_div_x( valmu - valGenAngSdvig));
}



// построение графика функции F греческое для метода 3 парциальных диаграмм в зависимости от обобщенного угла
void createGraphFGreece_For_3PartDiagr_from_GenAng(wchar_t *Fold, double valLambda
	,double  valAppert ,double  valGenAngSdvig,double  valSmeshenieGoriz
	, double  valSMeshenieVert, double val_b0, double val_b1 )
{
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
   parrBuff[ i * nBuffCols + 1] = fncFGr_PartDiagr(arr_b, valGenAngSdvig, valmu) + valSMeshenieVert;

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
}


// построение основной картинки презентации
// INPUT
//  valScaleY  - растяжение пол оси Y
//
//
//
//
//


void createPtesentGraphs_for_PartDiagr( wchar_t *wchrPresntPartDiagrams, double valLambda
	   , double valAppert , double valAngSdvig, int iNumPartRayTriple
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

   for (int i =0; i < 7; i++)
   {
	wchar_t FileName[400] ={0.};
	wcscpy( FileName, Fold);
	wcscat(FileName, L"\\PartDiagrNo_");
	wchar_t string[10] = {0.};
	_itow(i,string, 10);
	wcscat(FileName, string);
	wcscat(FileName, L".shp");
	double valSmeshenieGoriz =  ( -3. + ((double) i)) * valAngSdvig;
	double valSMeshenieVert = 0.;
	createGraphPartDiagr_from_rad(FileName, valLambda
	, valAppert , valAngSdvig, valSmeshenieGoriz, valSMeshenieVert, valScaleY);
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
	double valSmeshenieGoriz = ( -3. + ((double)(iNumPartRayTriple + i))) * valAngSdvig;

	createGraphPartDiagr_from_rad(FileName, valLambda
	, valAppert , valAngSdvig, valSmeshenieGoriz, 0., valScaleY ) ;
   }
   ///

   // построение графика функции FGreece в зависимости от угла рад
   wchar_t Fold1[400] ={0.};
   wcscpy( Fold1, wchrPresntPartDiagrams);
   wcscat(Fold1, L"\\GraphFGreece_from_rad_For_PartDiagr");
	_wmkdir(Fold1);
	wcscat(Fold1, L"\\");

   // относиьтельный угол цетрального луча тройки в радианах
   double valRadCentreRayPos =  ( -2.+ ((double)iNumPartRayTriple)) * valAngSdvig;

   createGraphFGreece_from_rad_For_3PartDiagr(Fold1,  valLambda
	, valAppert ,valAngSdvig,valRadCentreRayPos , valSMeshenieVert, val_b0, val_b1)  ;
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
void createGraphPartDiagr_from_rad(wchar_t *FileName, double valLambda
	,double  valAppert ,double  valAngSdvig,double  valSmeshenieGoriz, double  valSMeshenieVert,  double valScaleY )
{

	double step =0.0001;
	double  valGenAngSdvig = transformAngToGeneralizedAng (valLambda ,  valAppert, valAngSdvig  ) ;
	double valZero =  M_PI;
	const double VAlDiap = 3. *valZero;  // диапахзон графиков по обобщ координате плюс минус
	const int nBuffRows = 2. *VAlDiap/ step ;
	 TURPolyLine pln( 1, nBuffRows) ;



  double tet = -VAlDiap;
  for (int i=0 ; i < nBuffRows; i++)
  {
   tet = -VAlDiap + ((double)i) * step;
   pln.Points[i].X = (transformGeneralizedAngToAng (valLambda ,  valAppert, tet  ) + valSmeshenieGoriz)* 100. ;
   pln.Points[i].Y =fncDiagrSinx_div_x(tet) * valScaleY + valSMeshenieVert;

  }

 pln.WriteSetSHPFiles(FileName, &pln, 1);
}

// построение графика функции F греческое для метода 3 парциальных диаграмм в зависимости от угла в рад
void createGraphFGreece_from_rad_For_3PartDiagr(wchar_t *Fold, double valLambda
	,double  valAppert ,double  valAngSdvig,double  valSmeshenieGoriz
	, double  valSMeshenieVert, double val_b0, double val_b1 )
{
	double arr_b[2] = {0.};
	arr_b[0] =  val_b0 ;
	arr_b[1]  = val_b1;
   // построение графика функции FGreece
   double valGenAngSdvig = transformAngToGeneralizedAng (valLambda , valAppert, valAngSdvig ) ;
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
   parrBuff[ i * nBuffCols] = (transformGeneralizedAngToAng (valLambda , valAppert, valmu ) + valSmeshenieGoriz) * 100.;
   parrBuff[ i * nBuffCols + 1] = fncFGr_PartDiagr(arr_b, valGenAngSdvig, valmu) + valSMeshenieVert;

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
}

// построение графико СКЗ оши бок оценивания в звисимости от углового расстояния
// INPUT:
// Fold
// valLambda - раб частота
// valAppert - апертура антиенны
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
void createGraphsSKZ_from_AngDiffer_For_3PartDiagr(wchar_t *Fold, double valLambda
	,double  valAppert ,double  valAngSdvig,double  valSmeshenieGoriz
	, double  valSMeshenieVert, TComp cmpKTarg, TComp cmpKAntp , double valNoiseSkz,double valmAmplFactSig)
{

   double valGenAngSdvig = transformAngToGeneralizedAng (valLambda , valAppert, valAngSdvig ) ;
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

  double valRadAngTarg = valAngSdvig /2.;
  double valGenAngTarg =  transformAngToGeneralizedAng (valLambda , valAppert, valRadAngTarg ) ;
  double valGenAngAntp = 0.;
  double   arrMtrxCorr[4] ={0.}, arrMtrxCorrMu[4] = {0.};

  for (int i=0 ; i < nBuffRows; i++)
  {
   valGenAngAntp = valGenAngTarg - ((double)i + 10.) * valStep;


   	// вычисление коррел матрицы ошибок определения обобщенных углов цели и антипода
	calcMtrxCorrGenAngs_Meth3PartDiagr(valGenAngTarg ,  valGenAngAntp , valGenAngSdvig
	  , cmpKTarg, cmpKAntp,valNoiseSkz,  valmAmplFactSig, arrMtrxCorrMu) ;
	///

	// вычисление корреляц матрицы ошибок определения углов цели и антипода в радианаах
	double arrQ[4] = {0.}, arrTemp[4] = {0.};
	double coef = valLambda / valAppert / M_PI;
	arrQ[0] =  coef / sqrt( 1.- coef * valGenAngTarg * coef * valGenAngTarg);
	arrQ[3] =  coef / sqrt( 1.- coef *  valGenAngAntp * coef *  valGenAngAntp);
	MtrxMultMatrx(arrQ,2, 2, arrMtrxCorrMu,2, arrTemp)  ;
	MtrxMultMatrxTransp(arrTemp,2, 2, arrQ,2, arrMtrxCorr) ;


   parrBuff[ i * nBuffCols] = valRadAngTarg - transformGeneralizedAngToAng (valLambda , valAppert,valGenAngAntp);

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
	 createGraphPartDiagr_from_rad(wchAxesFileName,  valLambda
	,  valAppert ,  valAngSdvig,  valSmeshenieGoriz,   valSMeshenieVert,  4. );

		wcscpy(  wchAxesFileName,  Fold);
	wcscat(wchAxesFileName, L"DigarPart0.shp");
	 createGraphPartDiagr_from_rad(wchAxesFileName,  valLambda
	,  valAppert ,  valAngSdvig,  -valAngSdvig,   valSMeshenieVert,  4. );

		wcscpy(  wchAxesFileName,  Fold);
	wcscat(wchAxesFileName, L"DigarPart2.shp");
	 createGraphPartDiagr_from_rad(wchAxesFileName,  valLambda
	,  valAppert ,  valAngSdvig,  valAngSdvig,   valSMeshenieVert,  4. );

	// полигон изображающий истинная цель
	TURPointXY  pntTopRight2(0.05, cmpKTarg.modul()/ 800. * 5.);
	TURPointXY  pntBottomLeft2(- 0.05, 0.);
	TURPolygon plgTargTrue(pntTopRight2, pntBottomLeft2);
	wcscpy( wchAxesFileName, Fold);
	wcscat(wchAxesFileName, L"\\plgTargTrue.shp");
	plgTargTrue.WriteSetSHPFiles(wchAxesFileName ,&plgTargTrue, 1) ;
}



#pragma package(smart_init)
