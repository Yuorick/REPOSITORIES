//---------------------------------------------------------------------------


#pragma hdrstop
#include <string.h>
#include <math.h>
#include "Triple.h"
#include "SincDgr.h"
#include "Comp.h"
#include "DiagrSinX.h"
#include "MatrixProccess.h"
#include "URPolyLine.h"
#include "UrPointXY.h"
#include "YrWriteShapeFile.h"
extern const double TET0707;
extern const double VAL_C ;



TTriple::TTriple()
{

 // длина волны
   mLambda = 3. ;
 // угол сдвига парциальных диаграмм  радианы

 marrRadAngSdvig[0] = 0.025;
 marrRadAngSdvig[1] = 0.028;
 for (int i = 0; i < 3; i++)
 {
   marrDgr[i] = TSincDgr();
   marrDgr[i].mAppert = 98;
 }
}


// Конструктор копирования
TTriple::TTriple (const TTriple &R2)
 {
 // длина волны
   mLambda = R2. mLambda ;
 // угол сдвига парциальных диаграмм  радианы
  memcpy(marrRadAngSdvig, R2.marrRadAngSdvig, sizeof(double) * 3);
 //
  memcpy(marrDgr, R2.marrDgr, sizeof(TSincDgr) * 3);

 }

  // оператор присваивания
  TTriple TTriple::operator=(TTriple  R2)
{
 // длина волны
   mLambda = R2. mLambda ;
 // угол сдвига парциальных диаграмм  радианы
  memcpy(marrRadAngSdvig, R2.marrRadAngSdvig, sizeof(double) * 3);
 //
  memcpy(marrDgr, R2.marrDgr, sizeof(TSincDgr) * 3);
  return *this ;
}

__fastcall TTriple::TTriple(const double Lambda,double *arrRadAngSdvig, TSincDgr *arrDgr)
{
 mLambda = Lambda;
 memcpy( marrRadAngSdvig, arrRadAngSdvig, 2 * sizeof(double));
 memcpy(marrDgr, arrDgr, 3 * sizeof( TSincDgr ));
}

// оценивание обобщенных углов цели и антипода по измерениям 3 ПАРЦИАЛЬНЫХ диаграмм
//INPUT:
// cmparrS[3] -  массив замеров

//  OUTPUT:
//  *valEstAngTarg -   угол цели, рад
// *valEstAngAntp -    угол антипода, рад
// cmpKTarg - коэф отражения цели
// cmpKAntp - коэф отражения антипода
// arrMtrxCorr  - коррел матрица ошибок измерения угла цели и антипода
int TTriple::estimate( TComp *cmparrS, double *valEstAngTarg
, double *valEstAngAntp, TComp *cmpKTarg , TComp *cmpKAntp , double *arrMtrxCorr
  , double *pval_b0, double *pval_b1 )
{
	 int ireturn = -1;
	// нахождение обобщенных углов цели и антипода
	double  valGenTargEps = -100.,  valGenAntpEps = -100.0;  // обобщ углы цели и антипода
   //	if (EstGenAngsThreePartDiagr(  &cmparrS[iNumRayTriple]
   //	,  &valGenTargEps,  &valGenAntpEps, pval_b0, pval_b1) !=2)
		if (calculateGenAngs(  cmparrS,  &valGenTargEps,  &valGenAntpEps, pval_b0, pval_b1) !=2)
	{
	   return -1;
	}
	ireturn = 2;

	// вычисление углоы цели и антипода в рад относительно центрального луча центаральной (2-ой по счету) диаграммы
		  // относиьтельный угол цетрального луча тройки в радианах
	  /*		double valRadCentreRayPos = calcRadAngSdvigaPartDiagr(iNumRayTriple + 1);
			// ( - ((double)mQuantDiagr -1.) /2. + 1. + ((double)iNumRayTriple))				*mAngSdvig ;
		  // вычисление угла цели в радианах в АСК
			const double ValRadTargEps =  transformGeneralizedAngToAng ( valGenTargEps ) ;// в радианы
			*valEstAngTarg = ValRadTargEps + valRadCentreRayPos;
		  // вычисление угла антипода в радианах в АПСК
			const double ValRadAntpАEps =  transformGeneralizedAngToAng (valGenAntpEps ) ;// в радианы
			*valEstAngAntp = ValRadAntpАEps + valRadCentreRayPos; */


			*valEstAngTarg = marrDgr[1].transformGeneralizedAngToAng ( valGenTargEps,mLambda ) ;
			*valEstAngAntp = marrDgr[1].transformGeneralizedAngToAng ( valGenAntpEps,mLambda ) ;
	///

	// вычисление коэффиц отражения
	TComp cmparrA[4],cmparrAInv[4], cmparrK[2], cmparrSTemp[2];

  double temp0 = marrDgr[0].mAppert/ marrDgr[1].mAppert;
  double temp2 = marrDgr[2].mAppert/ marrDgr[1].mAppert;
  double valGenAngSdvig0 =  marrDgr[0].transformAngToGeneralizedAng(marrRadAngSdvig[0], mLambda ) ;
  double valGenAngSdvig2 =  marrDgr[2].transformAngToGeneralizedAng(marrRadAngSdvig[1], mLambda ) ;

	cmparrA[0] = TComp(fncDiagrSinx_div_x(valGenTargEps * temp0 + valGenAngSdvig0), 0.);
	cmparrA[1] = TComp(fncDiagrSinx_div_x (valGenAntpEps * temp0+ valGenAngSdvig0 ), 0.);
	cmparrA[2] = TComp(fncDiagrSinx_div_x (valGenTargEps * temp2 - valGenAngSdvig2 ), 0.);
	cmparrA[3] = TComp(fncDiagrSinx_div_x (valGenAntpEps * temp2- valGenAngSdvig2 ), 0.);

	bool brez =   InverseMtrx2( cmparrA, cmparrAInv);
	if (!brez)
	{
	 return -2;
	}

	cmparrSTemp[0] = cmparrS [0 ];
	cmparrSTemp[1] = cmparrS [ 2 ];
	MtrxMultMatrx(cmparrAInv,2, 2, cmparrSTemp,1, cmparrK);
	*cmpKTarg = cmparrK[0] ;
	*cmpKAntp = cmparrK[1] ;

	///



	// вычисление коррел матрицы ошибок определения обобщенных углов цели и антипода
	double arrMtrxCorrMu[4] = {0.};
	calcMtrxCorrGenAngs(cmparrS,valGenTargEps,  valGenAntpEps
	  ,*cmpKTarg,*cmpKAntp, arrMtrxCorrMu) ;

				///

	// вычисление корреляц матрицы ошибок определения углов цели и антипода в радианаах
	double arrQ[4] = {0.}, arrTemp[4] = {0.};
	double coef = mLambda / marrDgr[1].mAppert / M_PI;
	arrQ[0] =  coef / sqrt( 1.- coef * valGenTargEps * coef * valGenTargEps);
	arrQ[3] =  coef / sqrt( 1.- coef * valGenAntpEps * coef * valGenAntpEps);
	MtrxMultMatrx(arrQ,2, 2, arrMtrxCorrMu,2, arrTemp)  ;
	MtrxMultMatrxTransp(arrTemp,2, 2, arrQ,2, arrMtrxCorr) ;
	///

	return ireturn;
}


// оценивание обобщенных углов цели и антипода по измерениям 3 парциальных диаграмм
//INPUT:
// wchFoldName1 - путь к папаке с графиками
// cmparrS[3] -  массив замеров
//  OUTPUT:
//  *pvalGenTargEps -   угол цели, обобщ коорд
// *pvalGenAntpEps -    угол антипода, обобщ коорд
int TTriple::calculateGenAngs(TComp *cmparrS ,  double  *pvalGenTargEps
   ,  double  *pvalGenAntpEps, double *pval_b0, double *pval_b1)
{

  // вычисленние вектора b
   double arr_b[2] = {0.}  ;
   bool brez = calcVect_b(cmparrS,arr_b) ;
   ///
   *pval_b0 = arr_b[0];
   *pval_b1 = arr_b[1];

   if ((fabs(arr_b[0]) < 0.000000001)&& (fabs(arr_b[1]) < 0.000000001) )
   {
	 *pvalGenTargEps =  M_PI;
	 *pvalGenAntpEps = M_PI;
	 return 1;
   }


   double arrRoots[2] ={0.};
   int iNUmRoots = findRootsFgr(NULL,arr_b,  arrRoots);
   *pvalGenTargEps =  arrRoots[1];
   *pvalGenAntpEps =  arrRoots[0];
   return iNUmRoots ;



}


// нахождение корней функции FGr длшя метода 3 суммарных диаграмм
// возвращает к-во корней
 int TTriple::findRootsFgr(wchar_t *wchFileName,double *arr_b,  double *arrRoots)
 {
	double valGenAngSdvig0 = marrDgr[1].transformAngToGeneralizedAng (marrRadAngSdvig[0], mLambda  ) ;
	double valGenAngSdvig1 = marrDgr[1].transformAngToGeneralizedAng (marrRadAngSdvig[1], mLambda  ) ;

   int iNUmRoots = 0;

  // double valMu0 =  M_PI + max_(valGenAngSdvig0, valGenAngSdvig1) - 0.01;
   double valMu0 =  M_PI + max_(valGenAngSdvig0, valGenAngSdvig1)- 0.7;
   double valStep = 0.25;
	const int lenBuff = 2. *valMu0/ valStep ;
  double *parrBuff  = new double [lenBuff];

  double valmu = -valMu0;

  for (int i=0 ; i < lenBuff; i++)
  {
   valmu = -valMu0 + ((double)i) * valStep;
   parrBuff[ i ] = fncFGr(arr_b, valmu);

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
	 arrRoots [iNUmRoots ] = findRootMethChord_For_Fgr(  arr_b,  valX0, valX1) ;
	 iNUmRoots++;


   }
  }
  if(wchFileName)
  {
	 TURPolyLine pln(1, lenBuff ) ;
	 for (int i=0 ; i < lenBuff; i++)
  {
   pln.Points[i].X = -valMu0 + ((double)i) * valStep;
   pln.Points[i].Y = parrBuff[ i];
   pln.WriteSetSHPFiles(wchFileName,&pln, 1) ;

  }
  }
  delete parrBuff;
  return iNUmRoots;
 }


double TTriple::fncFGr(double *arr_b, double valmu)
{
	double valGenAngSdvig0 = marrDgr[0].transformAngToGeneralizedAng (marrRadAngSdvig[0], mLambda  ) ;
	double valGenAngSdvig2 = marrDgr[2].transformAngToGeneralizedAng (marrRadAngSdvig[1], mLambda  ) ;

  return (fncDiagrSinx_div_x(valmu) - arr_b[0] * fncDiagrSinx_div_x(valmu * marrDgr[0].mAppert/ marrDgr[1].mAppert +valGenAngSdvig0)
		 - arr_b[1] * fncDiagrSinx_div_x( valmu * marrDgr[2].mAppert/ marrDgr[1].mAppert- valGenAngSdvig2));

}


double TTriple::findRootMethChord_For_Fgr(double *arr_b, double  valX0, double valX1)
{
  double valMu0 = valX0 - 0.00001;
  double valMu1 = valX1 + 0.00001;
  double valMuRez =  valMu0 ;
  double eps = 0.001;
  for (int i = 0; i < 100; i++)
  {
	double f0 = fncFGr(arr_b, valMu0 ) ;
	double f1 = fncFGr(arr_b, valMu1 ) ;
	double valMuRez0 = valMu0 - f0 * ( valMu1 - valMu0)/ (f1 - f0);
	if (fabs(valMuRez0 - valMuRez) < eps)
	{
	 return valMuRez0;
	}
	double f2 = fncFGr(arr_b, valMuRez0 ) ;
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
// альтернативным способом через измерения диаграмм
// INPUT:
// valGenTargEps, valGenAntpEps - оценки обобщенных углов цели и антипода относительно центральной диаграммы тройки
// OUTPUT:
// arrMtrxCorrGenAngs[4] - коррел матрица ошибок оценивания углов цели и антоипода
bool TTriple::calcMtrxCorrGenAngs(TComp *cmparrS,  double valGenTargEps,  double valGenAntpEps
  ,TComp cmpKTarg ,TComp cmpKAntp , double *arrMtrxCorrGenAngs )
{


	 // проверки
	 double arr_b[2] ={0.};
	 calcVect_b(cmparrS, arr_b);
	 double val_b0 = arr_b[0];
	 double val_b1 = arr_b[1];

	 // ТЕСТ
	// calcTrueVect_b(double valRadAngTarg, double valRadAngAntp, double *arr_b)

  // формирование матирицы частных производных dB/dS
  // первая строка это градиент функции b0 по S
  // вторая строка это градиент функции b1 по S
  double arrdB_po_dS[12] = {0.};
  calc_dB0_po_dS(cmparrS,   val_b0,  val_b1, arrdB_po_dS);
  calc_dB1_po_dS(cmparrS,   val_b0,  val_b1, &arrdB_po_dS[6]);
  ///

  double arrGenAngs[2] = {0.};
  arrGenAngs[0] =  valGenTargEps;
  arrGenAngs[1] =  valGenAntpEps;
  ///

  // формирование матрицы частных производных углов по измерениям
  double arr_dMu_po_dS[12] = {0.};
  double temp0 = marrDgr[0].mAppert/ marrDgr[1].mAppert;
  double temp2 = marrDgr[2].mAppert/ marrDgr[1].mAppert;
  double valGenAngSdvig0 =  marrDgr[0].transformAngToGeneralizedAng(marrRadAngSdvig[0], mLambda ) ;
  double valGenAngSdvig2 =  marrDgr[2].transformAngToGeneralizedAng(marrRadAngSdvig[1], mLambda ) ;
  for (int i =0; i < 2; i++)
  {

	 double coeff = fncDerivDiagrSinx_div_x(arrGenAngs[i])
				   -  val_b0 *  fncDerivDiagrSinx_div_x(arrGenAngs[i]* temp0 + valGenAngSdvig0)* temp0
				   -  val_b1 *  fncDerivDiagrSinx_div_x(arrGenAngs[i]* temp2 - valGenAngSdvig2)* temp2 ;
	 double arrF[2] = {0.};
	 arrF[0] = -fncDiagrSinx_div_x(arrGenAngs[i] + valGenAngSdvig0) / coeff;
	 arrF[1] = -fncDiagrSinx_div_x(arrGenAngs[i] - valGenAngSdvig2) / coeff;
	// MtrxTranspMultMatrx(arrdB_po_dS,2, 6, arrF, 1, &arr_dMu_po_dS[ i * 6]) ;
	 MtrxMultMatrx(arrF,1, 2, arrdB_po_dS, 6, &arr_dMu_po_dS[ i * 6]) ;
  }

  ///

  // формирование корреляц матрицы измерений
  double arrCorrPartDiagr[36] = {0.};

  createMtrxCorrForMeasures( valGenTargEps,  valGenAntpEps
  , cmpKTarg , cmpKAntp , arrCorrPartDiagr );
  ///

  // вычисление коррел матрицы вектора дельтах
  double  arrtemp0[12] = {0.};
  MtrxMultMatrx(arr_dMu_po_dS,2, 6, arrCorrPartDiagr, 6, arrtemp0) ;
  MtrxMultMatrxTransp(arrtemp0,2, 6, arr_dMu_po_dS,2, arrMtrxCorrGenAngs) ;
  ///
  return true;
}

//формирование корреляционной матрицы ошибок первичных измерений
// ансамбля первичных диаграмм
// вектор первичных ошибок имеет вид ...S[i].m_Re S , S[i].m_Im...
// INPUT:
// valGenTargEps  - обобщ угол цели  в координатах диаграммы с номером 2  (средней)
// valGenAntpEps - обобщ угол антипода  в координатах диаграммы с номером 2  (средней)

// cmpKTarg, cmpKAntp - сигналы(комлексные) цели и антипода

// OTPUT:
// arrMtrxCorr[36] - корреляц матрица
void TTriple::createMtrxCorrForMeasures(double valGenTargEps,  double valGenAntpEps
  ,TComp cmpKTarg ,TComp cmpKAntp, double *arrMtrxCorr )
 {

	double arrRadTargEps[3] = {0.} // угол цели в диаграммах
	, arrRadAntpEps[3] = {0.};  // угол антипода в диаграммах
	arrRadTargEps[1] = marrDgr[1].transformGeneralizedAngToAng(valGenTargEps, mLambda);
	arrRadAntpEps[1] = marrDgr[1].transformGeneralizedAngToAng(valGenAntpEps, mLambda);

	arrRadTargEps[0] = arrRadTargEps[1] + marrRadAngSdvig[0];
	arrRadAntpEps[0] = arrRadAntpEps[1] + marrRadAngSdvig[0];

	arrRadTargEps[2] = arrRadTargEps[1] - marrRadAngSdvig[1];
	arrRadAntpEps[2] = arrRadAntpEps[1] - marrRadAngSdvig[1];

	memset(arrMtrxCorr, 0, 4 * 3* 3 * sizeof(double));

  // вычисление коррел матрицы, обусловленной разбросом коэффиц усиления диаграмм
	for (int i = 0; i < 3; i++)  // цикл по диаграммам в тройке
	{
	  // диаграмма с направления цели
	  double valDiagrTarg = marrDgr[i].fncDiagrFromRad(arrRadTargEps[i], mLambda );

	  // диаграмма с навления антипода
	  double valDiagrAntp = marrDgr[i].fncDiagrFromRad(arrRadAntpEps[i], mLambda );

	  // дисперсия действительной части сигнала
	   double valDispRe = ( (valDiagrTarg * cmpKTarg.m_Re * valDiagrTarg * cmpKTarg.m_Re)
						  + (valDiagrAntp * cmpKAntp.m_Re * valDiagrAntp * cmpKAntp.m_Re) )
						   * marrDgr[i].mAmplFactSig * marrDgr[i].mAmplFactSig;
	   // дисперсия мнипмой части сигнала
	   double valDispIm = ( (valDiagrTarg * cmpKTarg.m_Im * valDiagrTarg * cmpKTarg.m_Im)
						  + (valDiagrAntp * cmpKAntp.m_Im * valDiagrAntp * cmpKAntp.m_Im) )
						   * marrDgr[i].mAmplFactSig * marrDgr[i].mAmplFactSig;
	 ///
	  // заполнение диагональных элементов матрицы arrMtrxCorr в соответствующих строках
	  // действительная и комплек5сная части сигнала диаграммы с номером i занимают строки 2*i и 2*i+1 соответственно
	  arrMtrxCorr[2*3 * (2 * i) + 2 *i] =  valDispRe  + marrDgr[i].mNoiseDisp;
	  arrMtrxCorr[2*3 * (2 * i + 1) + 2 *i + 1] =  valDispIm  + marrDgr[i].mNoiseDisp;

	}
	return;
 }


void calc_dB0_po_dS(TComp *cmparrS, double val_b0, double val_b1, double *arrdB_po_dS)
{


 double temp = cmparrS[0].m_Re * cmparrS[2].m_Im - cmparrS[2].m_Re * cmparrS[0].m_Im;
 arrdB_po_dS[0] = -val_b0 * cmparrS[2].m_Im/ temp;
 arrdB_po_dS[1] =  val_b0 * cmparrS[2].m_Re/ temp;
 arrdB_po_dS[2] =  cmparrS[2].m_Im/ temp;
 arrdB_po_dS[3] =  -cmparrS[2].m_Re/ temp;
 arrdB_po_dS[4] =  -cmparrS[2].m_Im * val_b1/ temp;
 arrdB_po_dS[5] =  cmparrS[2].m_Re * val_b1/ temp;

}


void calc_dB1_po_dS(TComp *cmparrS,  double val_b0, double val_b1, double *arrdB_po_dS)
{
 double temp = cmparrS[0].m_Re * cmparrS[2].m_Im - cmparrS[2].m_Re * cmparrS[0].m_Im;
 arrdB_po_dS[0] = val_b0 * cmparrS[0].m_Im/ temp;
 arrdB_po_dS[1] =  -val_b0 * cmparrS[0].m_Re/ temp;
 arrdB_po_dS[2] = - cmparrS[0].m_Im/ temp;
 arrdB_po_dS[3] =  cmparrS[0].m_Re/ temp;
 arrdB_po_dS[4] =  cmparrS[0].m_Im * val_b1/ temp;
 arrdB_po_dS[5] =  -cmparrS[0].m_Re * val_b1/ temp;
}

// вычисление вектора b
// INPUT:
// cmparrS[3] - массив измерений тройки диаграмм
// OUTPUT:
// arr_b[2] -  коэффициенты уравнения
bool TTriple::calcVect_b(TComp *cmparrS, double *arr_b)
{
 // вычисленние вектора b

 double temp = cmparrS[0].m_Re * cmparrS[2].m_Im - cmparrS[2].m_Re * cmparrS[0].m_Im;
 if (fabs(temp)< 0.00000000001)
 {
   return false; }

 arr_b [0] = (cmparrS[2].m_Im * cmparrS[1].m_Re - cmparrS[2].m_Re* cmparrS[1].m_Im)/ temp;
 arr_b [1] = (-cmparrS[0].m_Im * cmparrS[1].m_Re + cmparrS[0].m_Re * cmparrS[1].m_Im)/ temp;

   return true ;
}


void TTriple::createSystErrGraphs( wchar_t *Fold)
{
	double valStep = 0.1;
	const double VAlDiap =  M_PI;  // диапахзон изменеия угла антипода по обобщенной координате
	const int nBuffRows = VAlDiap/ valStep ;
	const int nBuffCols = 15;
	double  *parrBuff = new double [nBuffRows  * nBuffCols] ;
	memset(parrBuff, 0, nBuffRows * nBuffCols * sizeof(double));

	const int lenName = 30;
	wchar_t *wcharrFileNames  = new wchar_t [30 * nBuffCols];
	memset(wcharrFileNames, 0, 30 * nBuffCols * sizeof(wchar_t));

	wcscpy( wcharrFileNames, L"RadAng");
	wcscpy( &wcharrFileNames[30], L"dTetTarg_po_dDel1");
	wcscpy( &wcharrFileNames[60], L"dTetTarg_po_dDel2");
	wcscpy( &wcharrFileNames[90], L"dTetTarg_po_dDel3");
	wcscpy( &wcharrFileNames[120], L"dTetTarg_po_da1");
	wcscpy( &wcharrFileNames[150], L"dTetTarg_po_da2");
	wcscpy( &wcharrFileNames[180], L"dTetTarg_po_da3");
	wcscpy( &wcharrFileNames[210], L"dTetTarg_po_dLambda");

	wcscpy( &wcharrFileNames[240], L"dTetAntp_po_dDel1");
	wcscpy( &wcharrFileNames[270], L"dTetAntp_po_dDel2");
	wcscpy( &wcharrFileNames[300], L"dTetAntp_po_dDel3");
	wcscpy( &wcharrFileNames[330], L"dTetAntp_po_da1");
	wcscpy( &wcharrFileNames[360], L"dTetAntp_po_da2");
	wcscpy( &wcharrFileNames[390], L"dTetAntp_po_da3");
	wcscpy( &wcharrFileNames[420], L"dTetAntp_po_dLambda");


   double valRadAngTargCur =  0.;
   for (int i=0 ; i < nBuffRows; i++)
  {
   double valGenAngAntpCur = -VAlDiap + ((double)i ) * valStep;
   double valRadAngAntpCur = marrDgr[1].transformGeneralizedAngToAng(valGenAngAntpCur, mLambda) ;
   parrBuff[ i * nBuffCols] =  valRadAngAntpCur ;
   parrBuff[ i * nBuffCols +1] = fnc_dTet_po_dDel1(valRadAngTargCur,valRadAngAntpCur) ;
   parrBuff[ i * nBuffCols +2] = fnc_dTet_po_dDel2(valRadAngTargCur,valRadAngAntpCur) ;
   parrBuff[ i * nBuffCols +3] = fnc_dTet_po_dDel3(valRadAngTargCur,valRadAngAntpCur) ;

   parrBuff[ i * nBuffCols +4] = fnc_dTet_po_da1(valRadAngTargCur,valRadAngAntpCur) ;
   parrBuff[ i * nBuffCols +5] = fnc_dTet_po_da2(valRadAngTargCur,valRadAngAntpCur) ;
   parrBuff[ i * nBuffCols +6] = fnc_dTet_po_da3(valRadAngTargCur,valRadAngAntpCur) ;
   parrBuff[ i * nBuffCols +7] = fnc_dTet_po_dLambda(valRadAngTargCur,valRadAngAntpCur) ;

   parrBuff[ i * nBuffCols +8] = fnc_dTet_po_dDel1(valRadAngAntpCur,valRadAngTargCur) ;
   parrBuff[ i * nBuffCols +9] = fnc_dTet_po_dDel2(valRadAngAntpCur,valRadAngTargCur) ;
   parrBuff[ i * nBuffCols +10] = fnc_dTet_po_dDel3(valRadAngAntpCur,valRadAngTargCur) ;

   parrBuff[ i * nBuffCols +11] = fnc_dTet_po_da1(valRadAngAntpCur,valRadAngTargCur) ;
   parrBuff[ i * nBuffCols +12] = fnc_dTet_po_da2(valRadAngAntpCur,valRadAngTargCur) ;
   parrBuff[ i * nBuffCols +13] = fnc_dTet_po_da3(valRadAngAntpCur,valRadAngTargCur) ;
   parrBuff[ i * nBuffCols +14] = fnc_dTet_po_dLambda(valRadAngAntpCur,valRadAngTargCur) ;
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
  pscaley[7] = 1.;
  pscaley[8] = 1.;
  pscaley[9] = 1.;
  pscaley[10] = 1.;
  pscaley[11] = 1.;
  pscaley[12] = 1.;
  pscaley[13] = 1.;
  pscaley[14] = 1.;
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
								  ,1000. //  масштаб по оси X
								  ,pscaley[i]// масштаб по оси Y
								   ) ;
  }




	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  Fold);
	wcscat(wchAxesFileName, L"AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,0., parrBuff[ 0]  * 1000. + 10.
	 ,-15.,15., 0.25) ;
	 delete pscaley ;
	 delete  parrBuff;
	 delete wcharrFileNames;
}


// частная производная неявно заданной функции решения (корня) уравнения   углов
// по ширине первой (с нулевым номером) диашраммы тройки
double TTriple::fnc_dTet_po_dDel1(double valRadAngTargCur, double valRadAngAntpCur)
{
  double arr_b[2] ={0.};
  calcTrueVect_b(valRadAngTargCur, valRadAngAntpCur, arr_b) ;
  double val_dFGr_po_dTet =  fnc_dFGr_po_dTet(valRadAngTargCur , arr_b);
  double val_dF1_po_dDelta05 = marrDgr[0].fncDerivDgr_po_dDelta05(valRadAngTargCur +marrRadAngSdvig[0], mLambda );
  double temp = arr_b[0] * val_dF1_po_dDelta05/ val_dFGr_po_dTet;


  return arr_b[0] * val_dF1_po_dDelta05/ val_dFGr_po_dTet ;

}


// частная производная неявно заданной функции решения (корня) уравнения   углов
// по ширине второй  (с первым номером) диашраммы тройки
double TTriple::fnc_dTet_po_dDel2(double valRadAngTargCur, double valRadAngAntpCur)
{
  double arr_b[2] ={0.};
  calcTrueVect_b(valRadAngTargCur, valRadAngAntpCur, arr_b) ;
  double val_dFGr_po_dTet =  fnc_dFGr_po_dTet(valRadAngTargCur , arr_b);
  double val_dF2_po_dDelta05 = marrDgr[1].fncDerivDgr_po_dDelta05(valRadAngTargCur, mLambda );

  return - val_dF2_po_dDelta05/val_dFGr_po_dTet ;
}

// частная производная неявно заданной функции решения (корня) уравнения   углов
// по ширине третьей  (со вторым номером) диашраммы тройки
double TTriple::fnc_dTet_po_dDel3(double valRadAngTargCur, double valRadAngAntpCur)
{
  double arr_b[2] ={0.};
  calcTrueVect_b(valRadAngTargCur, valRadAngAntpCur, arr_b) ;
  double val_dFGr_po_dTet =  fnc_dFGr_po_dTet(valRadAngTargCur , arr_b);
  double val_dF3_po_dDelta05 = marrDgr[2].fncDerivDgr_po_dDelta05(valRadAngTargCur - marrRadAngSdvig[1], mLambda );

  return arr_b[1] * val_dF3_po_dDelta05/val_dFGr_po_dTet ;
}


// частная производная функции уравнения углов (FGr)
// по таетт
double TTriple::fnc_dFGr_po_dTet(double valRadTet, double *arr_b)
{
   double val_dF2_po_dTet =  marrDgr[1].fncDerivDgr_po_dTet(valRadTet, mLambda ) ;
   double val_dF1_po_dTet =  marrDgr[0].fncDerivDgr_po_dTet(valRadTet + marrRadAngSdvig[0] , mLambda ) ;
   double val_dF3_po_dTet =  marrDgr[2].fncDerivDgr_po_dTet(valRadTet -marrRadAngSdvig[1], mLambda ) ;
	return val_dF2_po_dTet - arr_b[0] *val_dF1_po_dTet - arr_b[1] *val_dF3_po_dTet  ;

}

// частная производная неявно заданной функции решения (корня) уравнения   углов
// по углу сдвига между перой и второй диаграммами a1
double TTriple::fnc_dTet_po_da1(double valRadAngTargCur, double valRadAngAntpCur)
{
  double arr_b[2] ={0.};
  calcTrueVect_b(valRadAngTargCur, valRadAngAntpCur, arr_b) ;
  double val_dFGr_po_dTet =  fnc_dFGr_po_dTet(valRadAngTargCur , arr_b);
  double val_dF1_po_dTetta = marrDgr[0].fncDerivDgr_po_dTet(valRadAngTargCur +marrRadAngSdvig[0], mLambda );

  return arr_b[0] * val_dF1_po_dTetta/ val_dFGr_po_dTet ;
}

// частная производная неявно заданной функции решения (корня) уравнения   углов
// по углу сдвига  периой и второй диаграммами a2
double TTriple::fnc_dTet_po_da2(double valRadAngTargCur, double valRadAngAntpCur)
{
  double arr_b[2] ={0.};
  calcTrueVect_b(valRadAngTargCur, valRadAngAntpCur, arr_b) ;
  double val_dFGr_po_dTet =  fnc_dFGr_po_dTet(valRadAngTargCur , arr_b);
  double val_dF2_po_dTetta = marrDgr[1].fncDerivDgr_po_dTet(valRadAngTargCur , mLambda );

  return - val_dF2_po_dTetta/ val_dFGr_po_dTet ;
}


// частная производная неявно заданной функции решения (корня) уравнения   углов
// по углу сдвига диаграммы a3
double TTriple::fnc_dTet_po_da3(double valRadAngTargCur, double valRadAngAntpCur)
{
  double arr_b[2] ={0.};
  calcTrueVect_b(valRadAngTargCur, valRadAngAntpCur, arr_b) ;
  double val_dFGr_po_dTet =  fnc_dFGr_po_dTet(valRadAngTargCur , arr_b);
  double val_dF3_po_dTetta = marrDgr[2].fncDerivDgr_po_dTet(valRadAngTargCur - marrRadAngSdvig[1], mLambda );

  return -arr_b[1] * val_dF3_po_dTetta/ val_dFGr_po_dTet ;
}

// частная производная неявно заданной функции решения (корня) уравнения   углов
// по ширине первой (с нулевым номером) диашраммы тройки
double TTriple::fnc_dTet_po_dLambda(double valRadAngTargCur, double valRadAngAntpCur)
{
  double arr_b[2] ={0.};
  calcTrueVect_b(valRadAngTargCur, valRadAngAntpCur, arr_b) ;
  double val_dFGr_po_dTet =  fnc_dFGr_po_dTet(valRadAngTargCur , arr_b);
  double val_dFGr_po_dLambda = fnc_dFGr_po_dLambda(valRadAngTargCur , arr_b);
  return -val_dFGr_po_dLambda/ val_dFGr_po_dTet ;
}



// частная производная функции уравнения углов (FGr)
// по Lambda
double TTriple::fnc_dFGr_po_dLambda(double valRadTet, double *arr_b)
{
   double val_dF2_po_dLambda =  marrDgr[1].fncDerivDgr_po_dLambda(valRadTet, mLambda ) ;
   double val_dF1_po_dLambda =  marrDgr[0].fncDerivDgr_po_dLambda(valRadTet + marrRadAngSdvig[0] , mLambda ) ;
   double val_dF3_po_dLambda =  marrDgr[2].fncDerivDgr_po_dLambda(valRadTet -marrRadAngSdvig[1], mLambda ) ;
	return val_dF2_po_dLambda - arr_b[0] *val_dF1_po_dLambda - arr_b[1] *val_dF3_po_dLambda  ;

}
// вычисление вектора b  через истинные углы цели и антипода
// INPUT:
// valRadAngTarg, valRadAngAntp - углы цели и антипода
// OUTPUT:
// arr_b[2] -  коэффициенты уравнения
bool TTriple::calcTrueVect_b(double valRadAngTarg, double valRadAngAntp, double *arr_b)
{
 // вычисленние вектора b
 double f11 = marrDgr[0].fncDiagrFromRad(valRadAngTarg + marrRadAngSdvig[0], mLambda ) ;
 double f12 = marrDgr[0].fncDiagrFromRad(valRadAngAntp + marrRadAngSdvig[0], mLambda ) ;


 double f21 = marrDgr[1].fncDiagrFromRad(valRadAngTarg , mLambda ) ;
 double f22 = marrDgr[1].fncDiagrFromRad(valRadAngAntp , mLambda ) ;

 double f31 = marrDgr[2].fncDiagrFromRad(valRadAngTarg - marrRadAngSdvig[1], mLambda ) ;
 double f32 = marrDgr[2].fncDiagrFromRad(valRadAngAntp - marrRadAngSdvig[1], mLambda ) ;


 double temp = f12 *f31 - f11*f32;
 if (fabs(temp) < 0.000000000001)
 {
   return false;
 }

 arr_b [0] = (f31 * f22 - f32 * f21)/ temp;
 arr_b [1] = (f21 * f12 - f22 * f11)/ temp;
 double tt =  f22 - arr_b [0]* f12 - arr_b [1] * f32;
 return true ;
}

void TTriple::checkPlaneHypot(wchar_t *wchFileName
  ,TURPolyLine  plnModulGraph0,TURPolyLine  plnArgGraph0,TURPolyLine  plnModulGraph1
  ,TURPolyLine  plnArgGraph1,TURPolyLine  plnModulGraph2,TURPolyLine  plnArgGraph2)
{
 const int IargMax0 = plnModulGraph0.calcGraphArgMax() ;
 const int IargMax1 = IargMax0 + 50;//plnModulGraph1.calcGraphArgMax();
 // создание 2-х базисных векторов
 // векторы создаются в точках максимумоа 0-вой1 и 1-ой диаграмм

   // исходные векторы
   TComp cmpX0[3], cmpX1[3], cmpE0[3], cmpE1W[3], cmpE1[3], cmpZCUr[3], cmpZCUrNorm[3];
   cmpX0[0] = transfTrigonForm(plnModulGraph0.Points[IargMax0].Y, plnArgGraph0.Points[IargMax0].Y);
   cmpX0[1] = transfTrigonForm(plnModulGraph1.Points[IargMax0].Y, plnArgGraph1.Points[IargMax0].Y);
   cmpX0[2] = transfTrigonForm(plnModulGraph2.Points[IargMax0].Y, plnArgGraph2.Points[IargMax0].Y);

   cmpX1[0] = transfTrigonForm(plnModulGraph0.Points[IargMax1].Y, plnArgGraph0.Points[IargMax1].Y);
   cmpX1[1] = transfTrigonForm(plnModulGraph1.Points[IargMax1].Y, plnArgGraph1.Points[IargMax1].Y);
   cmpX1[2] = transfTrigonForm(plnModulGraph2.Points[IargMax1].Y, plnArgGraph2.Points[IargMax1].Y);
   ///

   // нормированные векторы
	MatrxDivideScalar(cmpX0, 1, 3, vectNorm(cmpX0, 3) , cmpE0);
	MatrxDivideScalar(cmpX1, 1, 3, vectNorm(cmpX1, 3) , cmpE1W);
	///

	// ортогонализация вектора cmpEW1
	TComp cmpScal = scalProd(cmpE1W, cmpE0, 3);
	TComp cmparrProj[3], cmparrT1[3];
	MatrxMultScalar(cmpE0, 3, 1, cmpScal ,cmparrProj);
	MtrxMinusMatrx(cmpE1W, cmparrProj,3, 1, cmparrT1);
	MatrxDivideScalar(cmparrT1, 1, 3, vectNorm(cmparrT1, 3) , cmpE1);

///

// построение графика
TComp cmparrT0[3], cmparrT4[3], cmparrT2[3],cmparrT3[3];

   TURPolyLine  plnPerpMod =  plnModulGraph0;
   for (int i = 0; i < plnModulGraph0.NumPoints; i++)
   {
	 cmpZCUr[0] = transfTrigonForm(plnModulGraph0.Points[i].Y, plnArgGraph0.Points[i].Y);
	 cmpZCUr[1] = transfTrigonForm(plnModulGraph1.Points[i].Y, plnArgGraph1.Points[i].Y);
	 cmpZCUr[2] = transfTrigonForm(plnModulGraph2.Points[i].Y, plnArgGraph2.Points[i].Y);

	 MatrxDivideScalar(cmpZCUr, 1, 3, vectNorm(cmpZCUr, 3) , cmpZCUrNorm);

	 TComp cmp0 = scalProd(cmpZCUr, cmpE0, 3);
	 TComp cmp1 = scalProd(cmpZCUr, cmpE1, 3);

	 MatrxMultScalar(cmpE0, 3, 1, cmp0 ,cmparrT0);
	 MatrxMultScalar(cmpE1, 3, 1, cmp1 ,cmparrT4);

	 MtrxMinusMatrx(cmpZCUr, cmparrT0,3, 1, cmparrT2);
	 MtrxMinusMatrx(cmparrT2, cmparrT4,3, 1, cmparrT3);



	 plnPerpMod.Points[i].Y = vectNorm_(cmparrT3, 3) ;

   }


	TURPolyLine::WriteSetSHPFiles(wchFileName,&plnPerpMod, 1) ;




}


#pragma package(smart_init)
