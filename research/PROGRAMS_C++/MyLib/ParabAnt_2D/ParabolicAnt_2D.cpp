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
#include "ParabolicAnt_2D.h"
#include "SincDgr.h";
#include "Triple.h"
#include "MeasStand.h"
#include "PareDgrs.h"
extern const double TET0707;
extern const double VAL_C ;
extern const double VAlCritAngSdvig =  2.928;
extern const double VAL_C = 299792458.;

TParabolicAnt_2D::TParabolicAnt_2D()
{
 // к-во парциальных диагшрамм в веере
  mQuantDiagr = 7;

 // длина волны
   mLambda = 3. ;
 // угол сдвига парциальных диаграмм  радианы

  mparrAngAbsSdvig = new double [mQuantDiagr];
  mparrDgr = new TSincDgr[ mQuantDiagr];
  for (int i = 0; i < mQuantDiagr; i++)
  {
	mparrAngAbsSdvig[i] =  1.5026 *M_PI/ 180.;
	mparrDgr [i] =  TSincDgr();
  }

}

__fastcall TParabolicAnt_2D::~TParabolicAnt_2D()
{
	if(mparrDgr) delete [] mparrDgr ;
	mparrDgr = NULL ;
	if(mparrAngAbsSdvig) //delete mparrAngAbsSdvig ;
	mparrAngAbsSdvig = NULL ;
}

// Конструктор копирования
TParabolicAnt_2D::TParabolicAnt_2D (const TParabolicAnt_2D &R2)
 {
// к-во парциальных диагшрамм в веере
  mQuantDiagr = R2.mQuantDiagr;

 // длина волны
   mLambda = R2. mLambda ;
 // угол сдвига парциальных диаграмм  радианы
  memcpy(mparrAngAbsSdvig, R2.mparrAngAbsSdvig, sizeof(double) * mQuantDiagr);
 //
  memcpy(mparrDgr, R2.mparrDgr, sizeof(TSincDgr) * mQuantDiagr);

 }

  // оператор присваивания
  TParabolicAnt_2D TParabolicAnt_2D::operator=(TParabolicAnt_2D  R2)
{
// к-во парциальных диагшрамм в веере
  mQuantDiagr = R2.mQuantDiagr;

 // длина волны
   mLambda = R2. mLambda ;
 // угол сдвига парциальных диаграмм  радианы
  memcpy(mparrAngAbsSdvig, R2.mparrAngAbsSdvig, sizeof(double) * mQuantDiagr);
 //
  memcpy(mparrDgr, R2.mparrDgr, sizeof(TSincDgr) * mQuantDiagr);
  return *this ;
}


// парам констр  1
 __fastcall TParabolicAnt_2D::TParabolicAnt_2D(const int quantDiagr,const double Lambda
   ,double  *parrAngSdvig, TSincDgr *parrDgr)
 {
	 mLambda = Lambda ;
	 mQuantDiagr = quantDiagr;
	if (mparrAngAbsSdvig != NULL)
	{
   //	delete  parrAngSdvig;
	mparrAngAbsSdvig = NULL;
	}
	if(parrAngSdvig  != NULL)
	{
	   mparrAngAbsSdvig  = new double[quantDiagr];
	   if(mparrAngAbsSdvig == NULL)
	   {
		   ShowMessage(L"Not memory for mparrAngAbsSdvig") ;
		   Abort() ;
	   }
	   memcpy(mparrAngAbsSdvig, parrAngSdvig, mQuantDiagr * sizeof(double));
	 }
	////////////////////////

	if (mparrDgr != NULL)
	{
   //	delete [] mparrDgr;
	mparrDgr = NULL;
	}
	if(parrDgr  != NULL)
	{
	   mparrDgr  = new TSincDgr[quantDiagr];
	   if(mparrDgr == NULL)
	   {
		   ShowMessage(L"Not memory for mparrDgr") ;
		   Abort() ;
	   }
	   memcpy(mparrDgr, parrDgr, mQuantDiagr * sizeof(TSincDgr));
	 }
 }

 // парам констр  1
 __fastcall TParabolicAnt_2D::TParabolicAnt_2D(const int quantDiagr,const double Lambda
   , double *parrAngSdvig,  double *arrRadTet05, const double AmplFactSig,  const double NoiseSkz )
 {
    mLambda = Lambda ;
	 mQuantDiagr = quantDiagr;
	if (mparrAngAbsSdvig != NULL)
	{
   //	delete  parrAngSdvig;
	mparrAngAbsSdvig = NULL;
	}
	if(parrAngSdvig  != NULL)
	{
	   mparrAngAbsSdvig  = new double[quantDiagr];
	   if(mparrAngAbsSdvig == NULL)
	   {
		   ShowMessage(L"Not memory for mparrAngAbsSdvig") ;
		   Abort() ;
	   }
	   memcpy(mparrAngAbsSdvig, parrAngSdvig, mQuantDiagr * sizeof(double));
	 }
	////////////////////////

	if (mparrDgr != NULL)
	{
   //	delete [] mparrDgr;
	mparrDgr = NULL;
	}

	mparrDgr  = new TSincDgr[quantDiagr];
	if(mparrDgr == NULL)
	{
	   ShowMessage(L"Not memory for mparrDgr") ;
	   Abort() ;
	}
	  for (int i = 0; i < quantDiagr; i++)
	  {
		mparrDgr[i] = TSincDgr(arrRadTet05[i], NoiseSkz, AmplFactSig, mLambda);
	  }
 }


 ////////////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////////

// имитация замеров   диаграмм
// INPUT:
// valalfUMTrg - угол цели в рад
// valalfUMAntp - угол антип в рад
// cmpKAntp
// OUTPUT:
// cmparrPartS []
// cmparrPartSZv []
void TParabolicAnt_2D::ImitateMeasureArray( double valalfUMTrg, TComp cmpKTarg, double valalfUMAntp,TComp cmpKAntp
		, TComp *cmparrPartS, TComp *cmparrPartSZv)
{

   for (int i =0; i < mQuantDiagr; i++)
   {
	  // вычисление угла целей относителбно оси диаграммы с номером i в обобщенных координатах
	  double valDiagrAng = mparrAngAbsSdvig[i] ; // угол оси диагр с номером i

	  double valTrg =  valalfUMTrg - valDiagrAng; // угол цели в диаграмме с номером i, рад.
	  double valTargGen = mparrDgr[i].transformAngToGeneralizedAng (  valTrg , mLambda) ; // угол цели в диаграмме с номером i, обощ угол
	  double valAntp =  valalfUMAntp - valDiagrAng; // угол цели в диаграмме с номером i, рад.
	  double valAntpGen = mparrDgr[i].transformAngToGeneralizedAng ( valAntp, mLambda ) ; // угол цели в диаграмме с номером i, обощ угол
	  ///
	   mparrDgr[i].ImitateMeasure(  valTargGen,  cmpKTarg,   valAntpGen, cmpKAntp
		, &cmparrPartS[i] , &cmparrPartSZv[i]) ;
   }

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
int TParabolicAnt_2D::estimateMethThreeDiagr( TComp *cmparrS
  , int iNumRayTriple , double *valEstAngTarg, double *valEstAngAntp, TComp *cmpKTarg , TComp *cmpKAntp , double *arrMtrxCorr
  , double *pval_b0, double *pval_b1 )
{
	 double arrRadAngSdvig[2] = {0.};
	 arrRadAngSdvig[0] = mparrAngAbsSdvig[ iNumRayTriple + 1] - mparrAngAbsSdvig[ iNumRayTriple];
	 arrRadAngSdvig[1] = mparrAngAbsSdvig[ iNumRayTriple + 2] - mparrAngAbsSdvig[ iNumRayTriple+ 1];
	 TTriple Triple(mLambda,arrRadAngSdvig, &mparrDgr[ iNumRayTriple]);
	 if (Triple.estimate( &cmparrS[iNumRayTriple], valEstAngTarg, valEstAngAntp, cmpKTarg , cmpKAntp , arrMtrxCorr
  , pval_b0, pval_b1 ) != 2)
  {
	  return -1;
  }

	 *valEstAngTarg += mparrAngAbsSdvig[ iNumRayTriple + 1];
	 *valEstAngAntp += mparrAngAbsSdvig[ iNumRayTriple + 1];
	  return 0;

 }

 void TParabolicAnt_2D::processMeasureArrFromStandKaurovFile(wchar_t *wchFoldReport, int iNumTriple
	   ,TMeasStand* parrMeas, int iNumMeasures )
{


	 double valAngNormal = (mparrAngAbsSdvig[2] + mparrAngAbsSdvig[3])/ 2.;
	 double valAngOtnosit = mparrAngAbsSdvig[iNumTriple +1] - valAngNormal;
	 double arrRadAngSdvig[2] = {0.};
	 arrRadAngSdvig[0] = mparrAngAbsSdvig[iNumTriple + 1] - mparrAngAbsSdvig[iNumTriple];
	 arrRadAngSdvig[1] = mparrAngAbsSdvig[iNumTriple + 2] - mparrAngAbsSdvig[iNumTriple + 1];
	 TTriple tripleWork(mLambda,arrRadAngSdvig, &mparrDgr[iNumTriple]);

	 TURPolyLine plnTargEst( 10, iNumMeasures) ;
	 TURPolyLine plnAntpEst( 10, iNumMeasures) ;

	 plnTargEst.NumParts =0;
	 plnTargEst.NumPoints = 0;
	 plnTargEst.Parts[0] = 0;

	 plnAntpEst.NumParts =0;
	 plnAntpEst.NumPoints = 0;
	 plnAntpEst.Parts[0] = 0;

	 bool bFindNextPart = true;
	 int iNumParts = 0;
	 int iNumPoints = 0;

	 double valEstAngTarg =0., valEstAngAntp =0.,  arrMtrxCorr [4] = {0.}, pval_b0 =0., pval_b1 =0.;
	 TComp cmpKTarg , cmpKAntp ;


	 for (int i =0; i < iNumMeasures; i++)
	 {
	   if (!isMeasure(parrMeas[i].mpcmparrMeas) )
	   {       // это не замер
		 if (bFindNextPart)
		 {
		   continue;
		 }
		 else
		 {
		  bFindNextPart = true;
		  continue;
		 }

	   }
	   else
	   {  // это замер
		   if(!tripleWork.estimate( parrMeas[i].mpcmparrMeas,  &valEstAngTarg
			 ,&valEstAngAntp, &cmpKTarg , &cmpKAntp ,  arrMtrxCorr ,  &pval_b0,  &pval_b1 ))
		   {
			 bFindNextPart  = true;
			 continue;
		   }

			  double ang = 0.;
			  TComp::angBetveenComps(cmpKTarg, cmpKAntp, &ang)  ;
			  ang = ang * 180./ M_PI;
			  if (bFindNextPart)
			  {
			 iNumParts++;
			 plnTargEst.Parts[iNumParts -1] = iNumPoints ;
			 plnAntpEst.Parts[iNumParts -1] = iNumPoints ;
			 bFindNextPart = false;
			 if (iNumParts > 10)
			 {
			  break;
			 }
			 }


	   plnTargEst.Points[iNumPoints].X =  parrMeas[i].mRadAxeAnt * 1000.;
	   plnTargEst.Points[iNumPoints].Y =  (valEstAngTarg + /*parrMeas[i].mRadAxeAnt + */valAngOtnosit)* 1000.;
	   plnAntpEst.Points[iNumPoints].X =  parrMeas[i].mRadAxeAnt * 1000.;
	   plnAntpEst.Points[iNumPoints].Y =  (valEstAngAntp + /* parrMeas[i].mRadAxeAnt +*/ valAngOtnosit)* 1000.;
	   iNumPoints ++;
	   }
	 }

	 plnTargEst.NumPoints = iNumPoints;
	 plnTargEst.NumParts = iNumParts;
	 plnAntpEst.NumPoints = iNumPoints;
	 plnAntpEst.NumParts = iNumParts;

	wchar_t wchFileName[300] ={0};
	wcscpy(  wchFileName,  wchFoldReport);
	wcscat(wchFileName, L"\\TargEst.shp");
	plnTargEst.WriteSetSHPFiles(wchFileName, &plnTargEst,1 );

	wcscpy(  wchFileName,  wchFoldReport);
	wcscat(wchFileName, L"\\AntpEst.shp");
	plnAntpEst.WriteSetSHPFiles(wchFileName, &plnAntpEst,1 );

	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  wchFoldReport);
	wcscat(wchAxesFileName, L"\\AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-150., 150.
	 ,-150,150, 10.) ;


}



 void TParabolicAnt_2D::processMeasureArrFromStandKaurovFile_(wchar_t *wchFoldReport, int iNumTriple
	   ,TMeasStand* parrMeas, int iNumMeasures )
{


	 double valAngNormal = (mparrAngAbsSdvig[2] + mparrAngAbsSdvig[3])/ 2.;
	 double valAngOtnosit = mparrAngAbsSdvig[iNumTriple +1] - valAngNormal;
	 double arrRadAngSdvig[2] = {0.};
	 arrRadAngSdvig[0] = mparrAngAbsSdvig[iNumTriple + 1] - mparrAngAbsSdvig[iNumTriple];
	 arrRadAngSdvig[1] = mparrAngAbsSdvig[iNumTriple + 2] - mparrAngAbsSdvig[iNumTriple + 1];
	 TTriple tripleWork(mLambda,arrRadAngSdvig, &mparrDgr[iNumTriple]);





	 bool bFindNextPart = true;
	 int iNumParts = 0;
	 int iNumPoints = 0;

	 double valEstAngTarg =0., valEstAngAntp =0.,  arrMtrxCorr [4] = {0.}, pval_b0 =0., pval_b1 =0.;
	 TComp cmpKTarg , cmpKAntp ;
	 TURPointXY *ppntarrTarg  = new TURPointXY  [iNumMeasures];
	 TURPointXY *ppntarrAntp  = new TURPointXY  [iNumMeasures];
	 int iParts [10] = {0};

	 for (int i =0; i < iNumMeasures; i++)
	 {
	   if (i == 28)
	   {
		int iiii = 0;
	   }
	   if (!isMeasure(parrMeas[i].mpcmparrMeas) )
	   {       // это не замер
		 if (bFindNextPart)
		 {
		   continue;
		 }
		 else
		 {
		  bFindNextPart = true;
		  continue;
		 }

	   }
	   else
	   {  // это замер
		   if( tripleWork.estimate( parrMeas[i].mpcmparrMeas,  &valEstAngTarg
			 ,&valEstAngAntp, &cmpKTarg , &cmpKAntp ,  arrMtrxCorr ,  &pval_b0,  &pval_b1 ) != 2)
		   {

			if (i == 50)
			{
			double arr_b[2] = {0.};
			arr_b[0] = pval_b0;
			arr_b[1] = pval_b1;
			double arrRoots [2] ={0.};
			wchar_t wchFileName[300] ={0};
			wcscpy(  wchFileName,  wchFoldReport);
			wcscat(wchFileName, L"\\Fgr.shp");
			tripleWork.findRootsFgr(wchFileName,arr_b, arrRoots) ;
			int iiii = 0;
			}


			 bFindNextPart  = true;
			 continue;
		   }

			  double ang = 0.;
			  TComp::angBetveenComps(cmpKTarg, cmpKAntp, &ang)  ;
			  ang = ang * 180./ M_PI;
			  if (bFindNextPart)
			  {
			 iNumParts++;
			 iParts[iNumParts -1] = iNumPoints ;

			 bFindNextPart = false;
			 if (iNumParts > 10)
			 {
			  break;
			 }
			 }



	   ppntarrTarg [iNumPoints].X =  parrMeas[i].mRadAxeAnt * 1000.;
	   ppntarrTarg [iNumPoints].Y =  (valEstAngTarg + parrMeas[i].mRadAxeAnt + valAngOtnosit)* 1000.;
	   ppntarrAntp[iNumPoints].X =  parrMeas[i].mRadAxeAnt * 1000.;
	   ppntarrAntp[iNumPoints].Y =  (valEstAngAntp +  parrMeas[i].mRadAxeAnt + valAngOtnosit)* 1000.;
	   iNumPoints ++;
	   }
	 }

	  bool bpr = true;
	  while( bpr )
	  {
		bpr = false;
		for (int i = (iNumParts - 1); i >0; i--)
		{
			if ( (iParts [i] - iParts [i -1]) == 1 )
			{
			  bpr = true;
			  for (int j = iParts [i] ; j < iNumPoints -2; j++)
			  {
				ppntarrTarg [j]  =  ppntarrTarg [j +1];
				ppntarrAntp [j]  =  ppntarrAntp [j +1];
			  }

			  for (int j = i ; j < iNumParts -2; j++)
			  {
				iParts [j] = iParts [j + 1] -1;

			  }
			  iNumPoints--;
			  iNumParts--;

			}
		}
	 }

	 TURPolyLine plnTargEst( iNumParts, iNumPoints,iParts,ppntarrTarg) ;
	 TURPolyLine plnAntpEst( iNumParts, iNumPoints,iParts,ppntarrAntp) ;



	wchar_t wchFileName[300] ={0};
	wcscpy(  wchFileName,  wchFoldReport);
	wcscat(wchFileName, L"\\TargEst.shp");
	plnTargEst.WriteSetSHPFiles(wchFileName, &plnTargEst,1 );

	wcscpy(  wchFileName,  wchFoldReport);
	wcscat(wchFileName, L"\\AntpEst.shp");
	plnAntpEst.WriteSetSHPFiles(wchFileName, &plnAntpEst,1 );

	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  wchFoldReport);
	wcscat(wchAxesFileName, L"\\AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-200., 200.
	 ,-150,150, 10.) ;

	 delete [] ppntarrTarg;
	 delete [] ppntarrAntp;

	wcscpy(  wchFileName,  wchFoldReport);
	wcscat(wchFileName, L"\\Digr0.shp");
	double valSdvig  = ( 0.07959 - mparrAngAbsSdvig[iNumTriple] + (mparrAngAbsSdvig[2] + mparrAngAbsSdvig[3]) /2.) * 1000.;
	TURPointXY pntSdvig (valSdvig, -50);
	double scaley = 30.;
	 mparrDgr[iNumTriple].createSHP_Graph(wchFileName, mLambda,  2
	   ,  pntSdvig,1000,  scaley);

		wcscpy(  wchFileName,  wchFoldReport);
	wcscat(wchFileName, L"\\Digr1.shp");
	 valSdvig  = (0.07959 - mparrAngAbsSdvig[iNumTriple+1] + (mparrAngAbsSdvig[2] + mparrAngAbsSdvig[3]) /2.) * 1000.;
	TURPointXY pntSdvig1 (valSdvig, -50);
	 scaley = 30.;
	 mparrDgr[iNumTriple+1].createSHP_Graph(wchFileName, mLambda,  2
	   ,  pntSdvig1,1000,  scaley);

	   	wcscpy(  wchFileName,  wchFoldReport);
	wcscat(wchFileName, L"\\Digr2.shp");
	 valSdvig  = ( 0.07959 - mparrAngAbsSdvig[iNumTriple+ 2] + (mparrAngAbsSdvig[2] + mparrAngAbsSdvig[3]) /2.) * 1000.;
	TURPointXY pntSdvig2 (valSdvig, -50);
	 scaley = 30.;
	 mparrDgr[iNumTriple + 2].createSHP_Graph(wchFileName, mLambda,  2
	   ,  pntSdvig2,1000,  scaley);


		wcscpy(  wchFileName,  wchFoldReport);
	   wcscat(wchFileName, L"\\AxesSdv.shp");
	   TURPointXY pntSdvig3 (0., pntSdvig.Y);
	   TYrWriteShapeFile::CreateShpArrowedAxes(wchFileName,-200,200
	 ,-100.,100.,0.05, pntSdvig3) ;

}

bool isMeasure(TComp *pcmparrMeas)
{
	for (int i = 0; i < 3; i++)
	{
	  if (pcmparrMeas[i].modul() < 0.0000000001)
	  {
		return false;
	  }
	}
	return true;
}


 void TParabolicAnt_2D::processARSM_ForMeasureArrFromStandKaurovFile(wchar_t *wchFoldName1, int iNumPare
		,TMeasStand* parrMeas, int iNumMeasures, double valTargAngApriori )
{

	 double valAngMiddle = (mparrAngAbsSdvig[iNumPare] + mparrAngAbsSdvig[iNumPare + 1]) /2.;
	 double valAngNormal = (mparrAngAbsSdvig[2] + mparrAngAbsSdvig[3])/ 2.;
	 double valAngOtnosit = valAngMiddle  - valAngNormal;



	 TPareDgrs pareDiagrs(mLambda,mparrAngAbsSdvig[iNumPare + 1] - mparrAngAbsSdvig[iNumPare], mparrDgr[iNumPare]);

	wchar_t wchFoldName[300] ={0};
	wcscpy(  wchFoldName,  wchFoldName1);
	wcscat(wchFoldName, L"\\");


	const int nBuffRows = iNumMeasures ;
	const int nBuffCols = 5;
	double  *parrBuff = new double [nBuffRows  * nBuffCols] ;
	memset(parrBuff, 0, nBuffRows * nBuffCols * sizeof(double));

	const int lenName = 30;
	wchar_t *wcharrFileNames  = new wchar_t [30 * nBuffCols];
	memset(wcharrFileNames, 0, 30 * nBuffCols * sizeof(wchar_t));

	wcscpy( wcharrFileNames, L"AntennaAng");
	wcscpy( &wcharrFileNames[30], L"TargAngEst");
	wcscpy( &wcharrFileNames[2 * 30], L"ModulCmpK");

	wcscpy( &wcharrFileNames[3 * 30], L"SKZ");
	wcscpy( &wcharrFileNames[4 * 30], L"Psi2");
   //	wcscpy( &wcharrFileNames[5 * 30], L"SumF");
   //	wcscpy( &wcharrFileNames[6 * 30], L"SumFSdvigMin");
   //	wcscpy( &wcharrFileNames[7 * 30], L"SumFSdvigPlus");
  int iNumRows = 0;
  for (int i=0 ; i < nBuffRows; i++)
  {
   // положение цели отнгосительно РСН пары лучей
   double valRadRSNPosition = valAngOtnosit + parrMeas[i].mRadAxeAnt;
   if (fabs(valRadRSNPosition - valTargAngApriori) > pareDiagrs. mAngSdvig/2.)
   {
	 continue;
   }
   parrBuff[ iNumRows * nBuffCols] = parrMeas[i].mRadAxeAnt * 1000.; //  угол  полворота антенны
   double valEstPartDiaRadAngTarg = 0., valPartDiaDisp = 0., valPsi2 = 0.;
	TComp cmpPartDiaKTarg ;
   pareDiagrs.estimateARSMethPartDiagr( parrMeas[i].mpcmparrMeas,   &valEstPartDiaRadAngTarg,  &cmpPartDiaKTarg
		   ,  &valPartDiaDisp,  &valPsi2);
   double valTargPos = valEstPartDiaRadAngTarg + valRadRSNPosition;
   parrBuff[ iNumRows * nBuffCols + 1] = valTargPos * 1000.;
   parrBuff[ iNumRows * nBuffCols + 2] = cmpPartDiaKTarg.modul();
   parrBuff[ iNumRows * nBuffCols + 3] = sqrt(valPartDiaDisp) * 1000.;
   parrBuff[ iNumRows * nBuffCols + 4] =  valPsi2;
   iNumRows ++;
 }

 // double scalex = 100.;
  double *pscaley = new double [nBuffCols] ;
  pscaley[0] = 1.;
  pscaley[1] = 1.;
  pscaley[2] = 1.;
  pscaley[3] = 1.;
  pscaley[4] = 1.;
  for (int i = 1; i < nBuffCols; i++)
  {
  TYrWriteShapeFile::WriteOneReport(                 wchFoldName  // путь к папке
								  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  , nBuffCols // - к-во переменных о корорых накоплена информация в буфере
								  , iNumRows //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,lenName // максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i  // номер переменной по оси Y
								  ,1. //  масштаб по оси X
								  ,pscaley[i]// масштаб по оси Y
								   ) ;
  }

	delete parrBuff;
	delete pscaley;
  //	wchar_t wchAxesFileName[300] ={0};
  //	wcscpy(  wchAxesFileName,  wchFoldName);
 //	wcscat(wchAxesFileName, L"AxesArr_.shp");
 //	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-150., 150.
 //	 ,-pscaley[0] * 1.1,pscaley[0] * 1.1, 50000.*step) ;

  delete wcharrFileNames;
}
#pragma package(smart_init)
