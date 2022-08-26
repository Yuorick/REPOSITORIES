					//---------------------------------------------------------------------------


#pragma hdrstop

#include "HomingSituation_3D.h"
#include <math.h>
#include <string.h>
#include <dir.h>
#include "UrPointXY.h"
#include "URPolygon.h"
#include "URPolyLine.h"
#include "MatrixProccess.h"
#include <stdlib.h>
#include <string.h>
#include <vcl.h>
#include "Atmosphere.h"

#include "YrWriteShapeFile.h"
const int QUANT_COLS = 11 , QUANT_POINTS_MAX = 20000;




//---------------------------------------------------------------------------
THomingSituation_3D::~THomingSituation_3D()
{

   if ( mpwcharrFoldReport != NULL)
   {
   free( mpwcharrFoldReport );
   mpwcharrFoldReport = NULL ;
   }

}


THomingSituation_3D::THomingSituation_3D( TBombTraj_3D BombTraj_3D //
		, TCntlFuncPar CntlFuncPar
		,long double TCur
		, wchar_t *pwcharrFoldReport // путь у папке с отчетом
		)
{


   mBombTraj_3D = BombTraj_3D ;
   mBombTraj_3D.mTCur = TCur;
   mTCur = TCur;

   mCntlFuncPar = CntlFuncPar;

  // путь у папке с отчетом
   wchar_t  wcharrFoldReport [300] = L"E:\\AMETIST\\PROJECTS_C++\\IMPULS\\FINAL_REPORT";
   mpwcharrFoldReport = NULL ;
   if ((mpwcharrFoldReport = (wchar_t *)malloc(300 * sizeof (wchar_t))) != NULL)
	{
	wcscpy(mpwcharrFoldReport, pwcharrFoldReport);

	}
	else
	{
	ShowMessage(L"Not memory for mpwcharrFoldReport") ;
	Abort() ;
	}



}

//---------------------------------------------------------------------------


// конструктор копирования
 THomingSituation_3D ::THomingSituation_3D (const THomingSituation_3D &R)
 {

	// траектория АУБ
	mBombTraj_3D = R.mBombTraj_3D;
	mCntlFuncPar = R.mCntlFuncPar ;

	mTCur = R.mTCur;
 // путь у папке с отчетом
	mpwcharrFoldReport = NULL;
	if(R.mpwcharrFoldReport != NULL)
	{
	if ((mpwcharrFoldReport = (wchar_t *)malloc(300 * sizeof (wchar_t))) != NULL)
	{
	wcscpy(mpwcharrFoldReport, R.mpwcharrFoldReport);

	}
	else
	{
	ShowMessage(L"Not memory for mpwcharrFoldReport") ;
	Abort() ;
	}

	}

 }

 // оператор присваивания
 THomingSituation_3D THomingSituation_3D::operator=(THomingSituation_3D  R)
 {

	// траектория АУБ
	mBombTraj_3D = R.mBombTraj_3D;
	mCntlFuncPar = R.mCntlFuncPar ;


	mTCur = R.mTCur;
 // путь у папке с отчетом
	mpwcharrFoldReport = NULL;
	if(R.mpwcharrFoldReport != NULL)
	{
	if ((mpwcharrFoldReport = (wchar_t *)malloc(300 * sizeof (wchar_t))) != NULL)
	{
	wcscpy(mpwcharrFoldReport, R.mpwcharrFoldReport);

	}
	else
	{
	ShowMessage(L"Not memory for mpwcharrFoldReport") ;
	Abort() ;
	}

	}

	return *this ;
 }


 // интгрирование уравнений до момента  падения
// OUTPUT:
// valDHoriz - горизонтальтеая дальн точки падения
void THomingSituation_3D::fncMoveBomb_TO_ZeroAlt_AND_ShowGraphs( wchar_t *wcharrPath, long double &valDHoriz )
{
	const double DEL_T = 0.1;
	double *parrBuff;
	double *pscaleY;
	wchar_t *pwcharrFileNames ;
	const int lenName =30 ;// максимальная длина имени переменной
   if (wcharrPath)
   {

	parrBuff = new double [QUANT_COLS * QUANT_POINTS_MAX] ;
	memset (parrBuff, 0, QUANT_COLS * QUANT_POINTS_MAX * sizeof(double)) ;
	pwcharrFileNames = new wchar_t [ QUANT_COLS * lenName] ;
	memset (pwcharrFileNames, 0, QUANT_COLS * lenName* sizeof(wchar_t)) ;

	wcscpy( &pwcharrFileNames[ 0 * 30], L"t");
	wcscpy( &pwcharrFileNames[ 1 * 30], L"X");
	wcscpy( &pwcharrFileNames[ 2* 30],  L"Y");
	wcscpy( &pwcharrFileNames[ 3 * 30], L"PI");
	wcscpy( &pwcharrFileNames[ 4 * 30], L"V");
	wcscpy( &pwcharrFileNames[ 5 * 30], L"Tetta");
	wcscpy( &pwcharrFileNames[ 6 * 30], L"UVert");
	wcscpy( &pwcharrFileNames[ 7 * 30], L"M");
	wcscpy( &pwcharrFileNames[ 8 * 30], L"Cy");
	wcscpy( &pwcharrFileNames[ 9 * 30], L"Z");
	wcscpy( &pwcharrFileNames[10 * 30], L"PSI");

	pscaleY = new double  [QUANT_COLS] ;
	pscaleY[1] = 0.1;
	pscaleY[2] = 0.1;

	pscaleY[3] = 1000.;
	pscaleY[4] = 1.;
	pscaleY[5] = 1000.;
	pscaleY[6] = 1000.;
	pscaleY[7] = 1000.;
	pscaleY[8] = 1000.;
	pscaleY[9] = 0.1;
	pscaleY[10] = 1000.;
}

  int iCirc = 1000. / mBombTraj_3D.mStepInt ;
 // long double valTTemp = mTCur + ((long double)iCirc) * mStepInt ;
  int i = 0 ;
  int iNupPointsOut = 0;
  double valTOut = -DEL_T ;
  for ( i = 0; i < iCirc; i++)
  {

	  if (mTCur > 9.3099999)
	  {
		 int iii = 0;
	  }
	  long double valUGor = 0.,  valUVert = 0.;
	  fncCalcU( mTCur, valUGor, valUVert)  ;
	  mBombTraj_3D.fncEilerStep(valUGor, valUVert, mBombTraj_3D.mStepInt);
	  mTCur = mBombTraj_3D.mTCur;

	  // погтсроение графиов для проверки и отладки
   /*	  TURPointXY PointsLineVelo[2]; // линия мншгновенной сорости УАБ
	  PointsLineVelo[0] = TURPointXY(mBombTraj_3D.marrStrSK_VS [0], mBombTraj_3D.marrStrSK_VS [1]);
	  PointsLineVelo[1] = TURPointXY(mBombTraj_3D.marrStrSK_VS [0] - mBombTraj_3D.marrStrSK_VS [1]/ tan(mBombTraj_3D.marrStrSK_VS [4]), 0.);
	  TURPolyLine lineVelo(PointsLineVelo,2) ;
	  wchar_t wcharrFileName[300]={0};
   wcscpy(wcharrFileName, wcharrPath);
   wcscat(wcharrFileName, L"\\LIneVelo.shp");
   TURPolyLine::WriteSetSHPFiles(wcharrFileName,&lineVelo, 1) ;

   TURPointXY PointTarg(mBombTraj_3D.mBomb.mHomingHead.marrMesureTarg[0],mBombTraj_3D.mBomb.mHomingHead.marrMesureTarg[1])   ;
   wcscpy(wcharrFileName, wcharrPath);
   wcscat(wcharrFileName, L"\\Targ.shp");
   TURPointXY::WriteSetSHPFiles(wcharrFileName,&PointTarg, 1) ;
	*/
	  ///]



	  if (mBombTraj_3D.marrStrSK_VS [1] < 0.) break ;

	 // if (((double)mTCur > (valTOut + DEL_T - 1E-15 ))&& ( iNupPointsOut < QUANT_POINTS_MAX ))
	  if  ( iNupPointsOut < QUANT_POINTS_MAX )
	  {

	  valTOut = (double)mTCur ;

	   if (wcharrPath)
	   {
	   double *p = &parrBuff[iNupPointsOut *  QUANT_COLS] ;
	   p[0] =  (double)mTCur ;
	   p[1] = (double)mBombTraj_3D.marrStrSK_VS [0];
	   p[2] = (double)mBombTraj_3D.marrStrSK_VS [1];
	   p[3] = (double)mBombTraj_3D.marrStrSK_VS [2];
	   p[4] = (double)mBombTraj_3D.marrStrSK_VS [3];
	   p[5] = (double)mBombTraj_3D.marrStrSK_VS [4];
	   long double valUGor = 0., valUVert = 0.;
	   fncCalcU( mTCur, valUGor, valUVert );
	   p[6] = (double)valUVert;
	   long double valTay,  valDerivTay =0., valMach = -1., arrGradMach[5] = {0.};

	   fncCalcNormTemperature(p[2], valTay, valDerivTay) ;
	   mBombTraj_3D.fncCalcMach( valTay,  valDerivTay
	 ,valMach);
	   p[7] = (double)valMach;
	   long double  val_q = 0., arrGrad_q [5] = {0.};
	   mBombTraj_3D.fncCalc_q_and_Grad_q(val_q, arrGrad_q) ;
	   p[8] = mBombTraj_3D.mBomb.fncCy(valMach);
	   p[9] = (double)mBombTraj_3D.marrStrSK_VS [5];
	   p[10] = (double)mBombTraj_3D.marrStrSK_VS [6];

		}

	  iNupPointsOut++;

	  }



  }
   valDHoriz = sqrtl( mBombTraj_3D.marrStrSK_VS [0] * mBombTraj_3D.marrStrSK_VS [0] + mBombTraj_3D.marrStrSK_VS [5] * mBombTraj_3D.marrStrSK_VS [5]);
 if (wcharrPath)
 {
		wchar_t wcharrPath1[300] = {0} ;
		wcscpy(wcharrPath1, wcharrPath );
		wcscat(wcharrPath1, L"\\");

		for (int j = 1; j < QUANT_COLS ; j++)
		{


		TYrWriteShapeFile::WriteOneReport(wcharrPath1  // путь к папке
							  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
							  ,QUANT_COLS // - к-во переменных о корорых накоплена информация в буфере
							  ,iNupPointsOut //  - к-во точек
							  ,pwcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
							  ,lenName // максимальная длина имени переменной
							  ,0  // номер переменной по оси X
							  ,j  // номер переменной по оси Y
							  ,100. //  масштаб по оси X
							  ,pscaleY[j]  // масштаб по оси Y
							   );
		}

		TYrWriteShapeFile::WriteOneReport(wcharrPath1  // путь к папке
							  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
							  ,QUANT_COLS // - к-во переменных о корорых накоплена информация в буфере
							  ,iNupPointsOut //  - к-во точек
							  ,pwcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
							  ,lenName // максимальная длина имени переменной
							  ,1  // номер переменной по оси X
							  ,2  // номер переменной по оси Y
							  ,1.//  масштаб по оси X
							  ,1.  // масштаб по оси Y
							   );

		 TYrWriteShapeFile::WriteOneReport(wcharrPath1  // путь к папке
							  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
							  ,QUANT_COLS // - к-во переменных о корорых накоплена информация в буфере
							  ,iNupPointsOut //  - к-во точек
							  ,pwcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
							  ,lenName // максимальная длина имени переменной
							  ,1  // номер переменной по оси X
							  ,9  // номер переменной по оси Y
							  ,1.//  масштаб по оси X
							  ,1.  // масштаб по оси Y
							   );

		wchar_t wchFileName4[300] = {0} ;
		wcscpy(wchFileName4, wcharrPath );
		wcscat(wchFileName4, L"\\Axes.shp");
		TYrWriteShapeFile::CreateShpAxes(wchFileName4,-40000.,40000,-40000.,40000) ;

		delete parrBuff ;
		delete pwcharrFileNames ;
		delete pscaleY ;
   }

}



 // интгрирование уравнений до момента  падения
// OUTPUT:
// valDHoriz - горизонтальтеая дальн точки падения
void THomingSituation_3D::fncMoveClas_TO_BlastPoint_AND_ShowGraphs( wchar_t *wcharrPath, long double &valDHoriz )
{
	if (mCntlFuncPar.miControlType < 20  )
	{
	 ShowMessage(L"Error fncMoveClas_TO_BlastPoint_AND_ShowGraphs: miControlType < 20");
	 return;
	}
	const double DEL_T = 0.1;
	double *parrBuff;
	double *pscaleY;
	wchar_t *pwcharrFileNames ;
	const int lenName =30 ;// максимальная длина имени переменной
	const int Quant_Cols = QUANT_COLS + 2;
   if (wcharrPath)
   {

	parrBuff = new double [Quant_Cols * QUANT_POINTS_MAX] ;
	memset (parrBuff, 0, Quant_Cols * QUANT_POINTS_MAX * sizeof(double)) ;
	pwcharrFileNames = new wchar_t [ Quant_Cols * lenName] ;
	memset (pwcharrFileNames, 0, Quant_Cols * lenName* sizeof(wchar_t)) ;

	wcscpy( &pwcharrFileNames[ 0 * 30], L"t");
	wcscpy( &pwcharrFileNames[ 1 * 30], L"X");
	wcscpy( &pwcharrFileNames[ 2* 30],  L"Y");
	wcscpy( &pwcharrFileNames[ 3 * 30], L"PI");
	wcscpy( &pwcharrFileNames[ 4 * 30], L"V");
	wcscpy( &pwcharrFileNames[ 5 * 30], L"Tetta");
	wcscpy( &pwcharrFileNames[ 6 * 30], L"AlfAt");
	wcscpy( &pwcharrFileNames[ 7 * 30], L"M");
	wcscpy( &pwcharrFileNames[ 8 * 30], L"Cy");
	wcscpy( &pwcharrFileNames[ 9 * 30], L"CyqSmUmaxMG");
	wcscpy( &pwcharrFileNames[ 10 * 30], L"AngVisir");
	wcscpy( &pwcharrFileNames[ 11 * 30], L"DeltaAng");

	pscaleY = new double  [Quant_Cols] ;
	pscaleY[1] = 0.1;
	pscaleY[2] = 0.1;

	pscaleY[3] = 1000.;
	pscaleY[4] = 1.;
	pscaleY[5] = 1000.;
	pscaleY[6] = 1000.;
	pscaleY[7] = 1000.;
	pscaleY[8] = 1000.;
	pscaleY[9] = 1000.;
	pscaleY[10] = 1000.;
	pscaleY[11] = 1000.;
}

  int iCirc = 1000. / mBombTraj_3D.mStepInt ;
 // long double valTTemp = mTCur + ((long double)iCirc) * mStepInt ;
  int i = 0 ;
  int iNupPointsOut = 0;
  double valTOut = -DEL_T ;
  for ( i = 0; i < iCirc; i++)
  {

	  if (mTCur > 9.3099999)
	  {
		 int iii = 0;
	  }

	  long double valUGor = 0.,  valUVert = 0.;
	  fncCalcU( mTCur, valUGor, valUVert)  ;
	  mBombTraj_3D.fncEilerStep(valUGor, valUVert, mBombTraj_3D.mStepInt);

	  mTCur = mBombTraj_3D.mTCur;


	  // погтсроение графиов для проверки и отладки
   /*	  TURPointXY PointsLineVelo[2]; // линия мншгновенной сорости УАБ
	  PointsLineVelo[0] = TURPointXY(mBombTraj_3D.marrStrSK_VS [0], mBombTraj_3D.marrStrSK_VS [1]);
	  PointsLineVelo[1] = TURPointXY(mBombTraj_3D.marrStrSK_VS [0] - mBombTraj_3D.marrStrSK_VS [1]/ tan(mBombTraj_3D.marrStrSK_VS [4]), 0.);
	  TURPolyLine lineVelo(PointsLineVelo,2) ;
	  wchar_t wcharrFileName[300]={0};
   wcscpy(wcharrFileName, wcharrPath);
   wcscat(wcharrFileName, L"\\LIneVelo.shp");
   TURPolyLine::WriteSetSHPFiles(wcharrFileName,&lineVelo, 1) ;

   TURPointXY PointTarg(mBombTraj_3D.mBomb.mHomingHead.marrMesureTarg[0],mBombTraj_3D.mBomb.mHomingHead.marrMesureTarg[1])   ;
   wcscpy(wcharrFileName, wcharrPath);
   wcscat(wcharrFileName, L"\\Targ.shp");
   TURPointXY::WriteSetSHPFiles(wcharrFileName,&PointTarg, 1) ;
	*/
	  ///]





	  if  ( iNupPointsOut < QUANT_POINTS_MAX )
	  {

	  valTOut = (double)mTCur ;

	   if (wcharrPath)
	   {
	   double *p = &parrBuff[iNupPointsOut *  Quant_Cols] ;
	   p[0] =  (double)mTCur ;
	   p[1] = (double)mBombTraj_3D.marrStrSK_VS [0];
	   p[2] = (double)mBombTraj_3D.marrStrSK_VS [1];
	   p[3] = (double)mBombTraj_3D.marrStrSK_VS [2];
	   p[4] = (double)mBombTraj_3D.marrStrSK_VS [3];
	   p[5] = (double)mBombTraj_3D.marrStrSK_VS [4];
	   long double valUGor = 0., valUVert = 0.;
	   fncCalcU( mTCur, valUGor, valUVert );
	   p[6] = (double)valUVert;

	   long double valTay,  valDerivTay =0., valMach = -1., arrGradMach[5] = {0.};

	   fncCalcNormTemperature(p[2], valTay, valDerivTay) ;
	   mBombTraj_3D.fncCalcMach( valTay,  valDerivTay,valMach);
	   p[7] = (double)valMach;
	   long double  val_q = 0., arrGrad_q [5] = {0.};
	   mBombTraj_3D.fncCalc_q_and_Grad_q(val_q, arrGrad_q) ;
	   p[8] = mBombTraj_3D.mBomb.fncCy(valMach);
	   p[9] = p[8] * val_q * mBombTraj_3D.mBomb.mMidSq * mBombTraj_3D.mBomb.mMaxAt / mBombTraj_3D.mBomb.mMass / G_ZEMLI ;
	  // p[10] = (double)fncAngVisir(); !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ВНИМАНИЕ!!!!
	   p[11] = fabs (p[10] - p[5] + p[6]);

		}

	  iNupPointsOut++;

	  }



  }
   valDHoriz = mBombTraj_3D.marrStrSK_VS [0] ;
 if (wcharrPath)
 {
		wchar_t wcharrPath1[300] = {0} ;
		wcscpy(wcharrPath1, wcharrPath );
		wcscat(wcharrPath1, L"\\");

		for (int j = 1; j < Quant_Cols ; j++)
		{


		TYrWriteShapeFile::WriteOneReport(wcharrPath1  // путь к папке
							  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
							  ,Quant_Cols // - к-во переменных о корорых накоплена информация в буфере
							  ,iNupPointsOut //  - к-во точек
							  ,pwcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
							  ,lenName // максимальная длина имени переменной
							  ,0  // номер переменной по оси X
							  ,j  // номер переменной по оси Y
							  ,100. //  масштаб по оси X
							  ,pscaleY[j]  // масштаб по оси Y
							   );
		}

		TYrWriteShapeFile::WriteOneReport(wcharrPath1  // путь к папке
							  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
							  ,Quant_Cols // - к-во переменных о корорых накоплена информация в буфере
							  ,iNupPointsOut //  - к-во точек
							  ,pwcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
							  ,lenName // максимальная длина имени переменной
							  ,1  // номер переменной по оси X
							  ,2  // номер переменной по оси Y
							  ,1.//  масштаб по оси X
							  ,1.  // масштаб по оси Y
							   );

		wchar_t wchFileName4[300] = {0} ;
		wcscpy(wchFileName4, wcharrPath );
		wcscat(wchFileName4, L"\\Axes.shp");
		TYrWriteShapeFile::CreateShpAxes(wchFileName4,-40000.,40000,-40000.,40000) ;

		delete parrBuff ;
		delete pwcharrFileNames ;
		delete pscaleY ;
   }

}

// вычисление угла атаки (  управления)
// Расчет функции управления УАБ с использованием члена класса TCntlFuncPar mCntlFuncPar
// miControlType = 21 или 22- указывает тип управлкения
// в массивах  miArrParams и  mDblArrParams хранятся парает ры, описывающие оуправляющую фунцию
// В зависимости от типа управления miControlType форируется во внешней програме алгорит управления

 void THomingSituation_3D::fncCalcU( long double TCur,long double &valUGor, long double &valUVert)
 {


	int ind = 0;
	long double valCoef = 0.;
	switch( mCntlFuncPar.miControlType )
	{

		case 0:   // движение происходит по баллистике,  alfaAt = 0
		valUGor = 0.;
		valUVert = 0.;
		return;
		//
		case 21:
		valCoef = mCntlFuncPar.mDblArrParams[ mCntlFuncPar.miArrParams[0]] ;
		valUGor = sqrtl(1. - valCoef * valCoef)  *  mBombTraj_3D.mBomb.mMaxAt ;
		ind =( mCntlFuncPar.miArrParams[0]  % 2 == 0 ?  1 :-1);
		for (int i = 0; i < mCntlFuncPar.miArrParams[0]; i++)
		{
		if(TCur >= mCntlFuncPar.mDblArrParams[i] )
		ind = - ind;
		else
		break;
		}
		valUVert =  valCoef *  mBombTraj_3D.mBomb.mMaxAt * ((long double)ind);
		return  ;
		break;
		case 22:
		valCoef = mCntlFuncPar.mDblArrParams[ mCntlFuncPar.miArrParams[0]] ;
		valUGor = sqrtl(1. - valCoef * valCoef)  *  mBombTraj_3D.mBomb.mMaxAt ;
		ind =(mCntlFuncPar.miArrParams[0]  % 2 == 0 ?  -1 :1);
		for (int i = 0; i < mCntlFuncPar.miArrParams[0]; i++)
		{
		if(TCur >= mCntlFuncPar.mDblArrParams[i] )
		ind = - ind;
		else
		break;
		}
		valUVert =  valCoef *  mBombTraj_3D.mBomb.mMaxAt * ((long double)ind);
		return  ;
		break;

		case 23:
		valCoef = mCntlFuncPar.mDblArrParams[ mCntlFuncPar.miArrParams[0]] ;
		if (fabsl((fabsl(mBombTraj_3D.marrStrSK_VS[5]) - M_PI/2.)) < 0.001)
		{
		 valCoef = 0.9999;
		}
		valUGor = sqrtl(1. - valCoef * valCoef)  *  mBombTraj_3D.mBomb.mMaxAt ;
		ind =( mCntlFuncPar.miArrParams[0]  % 2 == 0 ?  1 :-1);
		for (int i = 0; i < mCntlFuncPar.miArrParams[0]; i++)
		{
		if(TCur >= mCntlFuncPar.mDblArrParams[i] )
		ind = - ind;
		else
		break;
		}
		valUVert =  valCoef *  mBombTraj_3D.mBomb.mMaxAt * ((long double)ind);
		return  ;
		break;

		case 24:
		valCoef = mCntlFuncPar.mDblArrParams[ mCntlFuncPar.miArrParams[0]] ;
		if (fabsl((fabsl(mBombTraj_3D.marrStrSK_VS[5]) - M_PI/2.)) < 0.001) valCoef = 0.9999;
		valUGor = sqrtl(1. - valCoef * valCoef)  *  mBombTraj_3D.mBomb.mMaxAt ;
		ind =(mCntlFuncPar.miArrParams[0]  % 2 == 0 ?  -1 :1);
		for (int i = 0; i < mCntlFuncPar.miArrParams[0]; i++)
		{
		if(TCur >= mCntlFuncPar.mDblArrParams[i] )
		ind = - ind;
		else
		break;
		}
		valUVert =  valCoef *  mBombTraj_3D.mBomb.mMaxAt * ((long double)ind);
		return  ;
		break;

		default:
		break;
	}




 }






// нахождение оптимаального управления для задачи по дальности методом  перебора
// для одной точи переключения
// INPUT:
// valShagT - шаг сетки  по времени
THomingSituation_3D THomingSituation_3D::fncFindOptimalControl_for_Dist_MethodPerebora_1Point(TCntlFuncPar &CntlFuncPar
	   ,const long double valShagT, long double& valDGorizOpt)
{
   THomingSituation_3D HomingSituation_3D  = *this ;
   TCntlFuncPar CntlFuncParTemp =  CntlFuncPar ;
   long double *parrTimes = CntlFuncParTemp.mDblArrParams;
   parrTimes[0] = 0.;
   parrTimes[1] = 0.;
   parrTimes[2] = 1000.;
   CntlFuncParTemp.miArrParams[0] = 3; // к-во точек переключения
   THomingSituation_3D HomingSituation_3DReturn  = *this ;
   long double valDGoroz = 0.;
   HomingSituation_3D.fncMoveBomb_TO_ZeroAlt_AND_ShowGraphs(NULL, valDGoroz) ;

   long double valT = HomingSituation_3D.mBombTraj_3D.mTCur;
   long double valTOpt = 0.;
   valDGorizOpt = valDGoroz;
   long double valTCur = 0., valDGorizCur = 0.;
   int nCircle = valT / valShagT;

   for (int  i=0; i < nCircle; i++)
   {
	THomingSituation_3D HomingSituation_3DCur = *this;// ( mTet0, mV0, mAltit, mX0, mStepInt, miControlType, marrTimes, 3, mBomb) ;
	HomingSituation_3DCur.mCntlFuncPar.miArrParams[0] = 3 ;
	valTCur = ((long double)i) *valShagT ;
	HomingSituation_3DCur.mCntlFuncPar.mDblArrParams[0] = 0.;
	HomingSituation_3DCur.mCntlFuncPar.mDblArrParams[1] =  valTCur;
	HomingSituation_3DCur.mCntlFuncPar.mDblArrParams[2] = 1000.;

   HomingSituation_3DCur.fncMoveBomb_TO_ZeroAlt_AND_ShowGraphs(NULL, valDGorizCur) ;
   if (valDGorizCur > valDGorizOpt )
   {
	valDGorizOpt = valDGorizCur;
	valTOpt =  valTCur ;
	HomingSituation_3DReturn = HomingSituation_3DCur ;
   }
  // else
  // {
  //	  if ((i > 1) && (valDGorizCur < (valDGorizOpt -2.)) )break ;
  // }
 }

  return HomingSituation_3DReturn;
}

int THomingSituation_3D::fncSign (long double a)
{
	return (a >=0.)? 1:-1;
}


// нахождение достижимых точек падения методом перебора
// для однкой точи переключения
// INPUT:
// valShagT - шаг сетки  по времени
void THomingSituation_3D::fncFindSetPointsOfApproachibility(const int nCircleCoef, const long double valShagT
		, TURPointXY **pparrPnt, int &lenArr)
{
  /*
  // Быстрый вариант
   int lenArrCur = 0;
  for (int i =0; i < nCircleCoef; i++)
  {
	 //long double valCoef = ((long double) (1./ nCircleCoef)) * ( (long double)(i + 1.));
	 long double valCoef = cosl( M_PI/ 2.*  ((long double) (1./ nCircleCoef)) * ( (long double)(i )));
	 THomingSituation_3D  HomingSituation_3DTemp = *this;
	 HomingSituation_3DTemp.mCntlFuncPar.mDblArrParams[0] = 0.;
	 HomingSituation_3DTemp.mCntlFuncPar.mDblArrParams[1] = 1000.;
	 HomingSituation_3DTemp.mCntlFuncPar.mDblArrParams[2] = valCoef;
	 HomingSituation_3DTemp.mCntlFuncPar.miArrParams[0] = 2;
	 long double valDGoroz = 0.;
	 HomingSituation_3DTemp.fncMoveBomb_TO_ZeroAlt_AND_ShowGraphs(NULL, valDGoroz) ;
	 long double valTCur = HomingSituation_3DTemp.mTCur;
	 int nCircleTCur =  valTCur / valShagT;

	 if ( (lenArrCur + 2) > lenArr)
	 {
	  lenArr = lenArrCur + 100;
	  *pparrPnt = ( TURPointXY *)realloc (*pparrPnt,sizeof( TURPointXY)* lenArr);

	 }

	 long double valDGorozMin = 10000000., valDGorozMax = 0. ;
	 for (int j = 0; j < nCircleTCur; j++)
	 {
		HomingSituation_3DTemp = *this;
		HomingSituation_3DTemp.mCntlFuncPar.miArrParams[0] = 3;
		HomingSituation_3DTemp.mCntlFuncPar.mDblArrParams[0] = 0.;
		HomingSituation_3DTemp.mCntlFuncPar.mDblArrParams[1] = valTCur - ((long double) j) * valShagT;
		HomingSituation_3DTemp.mCntlFuncPar.mDblArrParams[2] = valTCur;
		HomingSituation_3DTemp.mCntlFuncPar.mDblArrParams[3] = valCoef;
		HomingSituation_3DTemp.fncMoveBomb_TO_ZeroAlt_AND_ShowGraphs(NULL, valDGoroz) ;
		if (valDGoroz > valDGorozMax)
		{
		   valDGorozMax = valDGoroz;
		  (*pparrPnt)[ lenArrCur ].X = (double)HomingSituation_3DTemp.mBombTraj_3D.marrStrSK_VS[0] ;
		  (*pparrPnt)[ lenArrCur ].Y = (double)HomingSituation_3DTemp.mBombTraj_3D.marrStrSK_VS[5] ;

		}
		if (valDGoroz < valDGorozMin)
		{
		   valDGorozMin = valDGoroz;
		  (*pparrPnt)[ lenArrCur +1 ].X = (double)HomingSituation_3DTemp.mBombTraj_3D.marrStrSK_VS[0] ;
		  (*pparrPnt)[ lenArrCur +1 ].Y = (double)HomingSituation_3DTemp.mBombTraj_3D.marrStrSK_VS[5] ;

		}

	 }
	 lenArrCur += 2; ;

  }


	lenArr = 2 * lenArrCur;
	*pparrPnt = ( TURPointXY *)realloc (*pparrPnt,sizeof( TURPointXY)* lenArr);
	for (int i = 0; i < lenArrCur; i++)
	{
	  (*pparrPnt)[ lenArrCur + i ].X =  (*pparrPnt)[ i ].X;
	  (*pparrPnt)[ lenArrCur + i ].Y = -(*pparrPnt)[ i ].Y;
	}


  TURPointXY::fncCleanSetPoints(10., pparrPnt, lenArr);
 */

  //  Точный и медленный вариант
  int lenArrCur = 0;
  for (int i =0; i < nCircleCoef; i++)
  {
	 //long double valCoef = ((long double) (1./ nCircleCoef)) * ( (long double)(i + 1.));
	 long double valCoef = cosl( M_PI/ 2.*  ((long double) (1./ nCircleCoef)) * ( (long double)(i )));
	 THomingSituation_3D  HomingSituation_3DTemp = *this;
	 HomingSituation_3DTemp.mCntlFuncPar.mDblArrParams[0] = 0.;
	 HomingSituation_3DTemp.mCntlFuncPar.mDblArrParams[1] = 1000.;
	 HomingSituation_3DTemp.mCntlFuncPar.mDblArrParams[2] = valCoef;
	 HomingSituation_3DTemp.mCntlFuncPar.miArrParams[0] = 2;
	 long double valDGoroz = 0.;
	 HomingSituation_3DTemp.fncMoveBomb_TO_ZeroAlt_AND_ShowGraphs(NULL, valDGoroz) ;
	 long double valTCur = HomingSituation_3DTemp.mTCur;
	 int nCircleTCur =  valTCur / valShagT;

	 if ( (lenArrCur + nCircleTCur) > lenArr)
	 {
	  lenArr = lenArrCur + nCircleTCur;
	  *pparrPnt = ( TURPointXY *)realloc (*pparrPnt,sizeof( TURPointXY)* lenArr);

	 }
	 for (int j = 0; j < nCircleTCur; j++)
	 {
		HomingSituation_3DTemp = *this;
		HomingSituation_3DTemp.mCntlFuncPar.miArrParams[0] = 3;
		HomingSituation_3DTemp.mCntlFuncPar.mDblArrParams[0] = 0.;
		HomingSituation_3DTemp.mCntlFuncPar.mDblArrParams[1] = valTCur - ((long double) j) * valShagT;
		HomingSituation_3DTemp.mCntlFuncPar.mDblArrParams[2] = valTCur;
		HomingSituation_3DTemp.mCntlFuncPar.mDblArrParams[3] = valCoef;
		HomingSituation_3DTemp.fncMoveBomb_TO_ZeroAlt_AND_ShowGraphs(NULL, valDGoroz) ;
		(*pparrPnt)[ lenArrCur ].X = (double)HomingSituation_3DTemp.mBombTraj_3D.marrStrSK_VS[0] ;
		(*pparrPnt)[ lenArrCur ].Y = (double)HomingSituation_3DTemp.mBombTraj_3D.marrStrSK_VS[5] ;
		lenArrCur++ ;

	 }

  }


	lenArr = 2 * lenArrCur;
	*pparrPnt = ( TURPointXY *)realloc (*pparrPnt,sizeof( TURPointXY)* lenArr);
	for (int i = 0; i < lenArrCur; i++)
	{
	  (*pparrPnt)[ lenArrCur + i ].X =  (*pparrPnt)[ i ].X;
	  (*pparrPnt)[ lenArrCur + i ].Y = -(*pparrPnt)[ i ].Y;
	}


  TURPointXY::fncCleanSetPoints(10., pparrPnt, lenArr);

}

// нахождение достижимых точек падения методом перебора
// для однкой точи переключения
// INPUT:
// valShagT - шаг сетки  по времени
TURPolygon THomingSituation_3D::fncFindSetOfDopPoints(wchar_t *pwchOutFile, const int QuantVar, const long double valStepT )
{
	//  long double valStepInt = 0.01 ;
  //	 TBombTraj_3D BombTraj_3D( mTet0,mV0 ,mY0 ,0.,valStepInt ,valTBegin	);

   //	 const int  iControlType = 23;
   //	 int iArrParams[5] = {0};
   //	 iArrParams[0] = 3;
   //	 long double DblArrParams[20] ={0.};
	// DblArrParams[1] = 13.;
// DblArrParams[2] = 1000.;
   //	 DblArrParams[3] = 0.5;
   //	 TCntlFuncPar CntlFuncPar ( iControlType, iArrParams , DblArrParams);
	 THomingSituation_3D  HomingSituation_3D = *this; //(  BombTraj_3D , CntlFuncPar	,0., mpwchOutFile ) ;

	long double valDHorizPlan = -1.;
  //	HomingSituation_3D.fncMoveBomb_TO_ZeroAlt_AND_ShowGraphs( mpwchOutFile, valDHorizPlan );
  //	LabeledEdit15->Text =  FloatToStr(valDHorizPlan );
	 int lenArr = 1000;
	 TURPointXY *parrPoints = (TURPointXY*)malloc(sizeof(TURPointXY) * lenArr);
	 TURPointXY **pparrPoints = &parrPoints;
	 HomingSituation_3D.fncFindSetPointsOfApproachibility(QuantVar,valStepT,pparrPoints, lenArr) ;
	 wchar_t wchFileName[300] ={0};
	 wcscpy(wchFileName, pwchOutFile );
		wcscat(wchFileName, L"\\Points.shp");
	 TURPointXY::WriteSetSHPFiles(wchFileName, *pparrPoints, lenArr) ;

	 TURPolygon plg =   TURPolygon::Conv(*pparrPoints, lenArr) ;
	 wcscpy(wchFileName, pwchOutFile);
	 wcscat(wchFileName, L"\\Conv.shp");
	 TURPolygon::WriteSetSHPFiles(wchFileName, &plg, 1) ;

   /*	 HomingSituation_3D.mCntlFuncPar.miControlType = 23;
	 lenArr = 1000;
	 *pparrPoints = (TURPointXY*)realloc(*pparrPoints,sizeof(TURPointXY) * lenArr);
	 HomingSituation_3D.fncFindSetPointsOfApproachibility(50,1., pparrPoints, lenArr) ;
	 wcscpy(wchFileName, mpwchOutFile );
	 wcscat(wchFileName, L"\\Points1.shp");
	 TURPointXY::WriteSetSHPFiles(wchFileName, *pparrPoints, lenArr) ; */
	 free(parrPoints);
	return  plg;
}



#pragma package(smart_init)
