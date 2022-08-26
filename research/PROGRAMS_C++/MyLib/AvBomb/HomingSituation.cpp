//---------------------------------------------------------------------------


#pragma hdrstop

#include "HomingSituation.h"
#include <math.h>
#include <dir.h>
#include <string.h>
#include "UrPointXY.h"
#include "URPolygon.h"
#include "URPolyLine.h"
#include "MatrixProccess.h"
#include <stdlib.h>
#include <vcl.h>
#include "Atmosphere.h"
#include "FinalSituation.h"
#include "Sector.h"
 #include "JetPlg.h"
 #include "UABPlg.h"
 #include "FireBallPlg.h"
 #include "VesselPlg.h"
#include "YrWriteShapeFile.h"
#include "BlastSituation.h"
const int QUANT_COLS = 10 , QUANT_POINTS_MAX = 25000;

//---------------------------------------------------------------------------
THomingSituation::~THomingSituation()
{

   if ( mpwcharrFoldReport != NULL)
   {
   free( mpwcharrFoldReport );
   mpwcharrFoldReport = NULL ;
   }

}

 // парам конструктор 1
THomingSituation::THomingSituation(TShipTarg ShipTarg //
		, TBombTraj BombTraj //
		, TCntlFuncPar CntlFuncPar
		,long double TCur
		, wchar_t *pwcharrFoldReport // путь у папке с отчетом
		)
{
   mShipTarg = ShipTarg ;
   mShipTarg.mT = TCur;
   mBombTraj = BombTraj ;
   mBombTraj.mTCur = TCur;
   mTCur = TCur;
   mFilt  = TFilt ( ShipTarg, BombTraj);
   mCntlFuncPar = CntlFuncPar;

  // путь у папке с отчетом

   mpwcharrFoldReport = NULL ;
   if (pwcharrFoldReport != NULL)
   {
		if ((mpwcharrFoldReport = (wchar_t *)malloc(300 * sizeof (wchar_t))) != NULL)
		{
		int i = wcslen(pwcharrFoldReport);
		wcscpy(mpwcharrFoldReport, pwcharrFoldReport);
		int j = wcslen(mpwcharrFoldReport);
		int iii=0;

		}
		else
		{
		ShowMessage(L"Not memory for mpwcharrFoldReport") ;
		Abort() ;
		}
	}

}

//---------------------------------------------------------------------------
// парам конструктор 2
THomingSituation::THomingSituation(TShipTarg ShipTarg //
		, TBombTraj BombTraj //
		, TCntlFuncPar CntlFuncPar
		,TGlonass Glonass
		,long double TCur
		, wchar_t *pwcharrFoldReport // путь у папке с отчетом
		)
{
   mShipTarg = ShipTarg ;
   mShipTarg.mT = TCur;
   mBombTraj = BombTraj ;
   mBombTraj.mTCur = TCur;
   mTCur = TCur;
   mFilt  = TFilt ( ShipTarg, BombTraj);
   mCntlFuncPar = CntlFuncPar;
   mGlonass = Glonass ;

  // путь у папке с отчетом
  // wchar_t  wcharrFoldReport [300] = L"E:\\AMETIST\\PROJECTS_C++\\IMPULS\\FINAL_REPORT";
   mpwcharrFoldReport = NULL ;
   if(pwcharrFoldReport)
   {
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

}


// конструктор копирования
 THomingSituation ::THomingSituation (const THomingSituation &R)
 {
	// корабль цель
	mShipTarg = R.mShipTarg;
	// траектория АУБ
	mBombTraj = R.mBombTraj;
	mCntlFuncPar = R.mCntlFuncPar ;
	mFilt = R.mFilt ;
	mTCur = R.mTCur;
	mGlonass = R.mGlonass;
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
 THomingSituation THomingSituation::operator=(THomingSituation  R)
 {
	// корабль цель
	mShipTarg = R.mShipTarg;
	// траектория АУБ
	mBombTraj = R.mBombTraj;
	mCntlFuncPar = R.mCntlFuncPar ;

	mFilt = R.mFilt ;
	mTCur = R.mTCur;
	mGlonass = R.mGlonass;

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


 //---------------------------------------------------------------------------
//---------------------------------------------------------------------------

// Вычисление вероятности поражения методом стат испытаний на последне участке траектории
// INPUT:
// mpwchStat - имя папаки для вывода графики
// NR - число реализаций
// Bomb - УАБ
// mTStart - момент привязки нач условий движения
// mTet0 -  нач угол бросания  УАБ
//	mV0  - нач скоргсть УАБ
// mAltit - нач высота  УАБ
 //	 mX0 - нач. положение по оси X УАБ
//	 mStepInt шагт интьерирования по времени диф уравнений движения УАБ
// valKNav - коэф навигации
// ShipTarg  - цель
// arrCfW[5] - массив  возмущений в правой части диффер уравнений,  arrWSigSq[0] = arrWSigSq[1] = 0.;  arrWSigSq[0] < 1
//    DelW[i] = arrCfW[i] * W, i = 2,3,4. То есть, это коэффициент показывающий с какой точностью мат модель описывает ускорение
// OUTPUT:
// valProbab - вероятность поражения
void THomingSituation::fncCalcProbability_AND_ShowGraphs( wchar_t *pwchStat, const int NR, TBomb Bomb,long double valTBegin
, long double valTet0, long double valV0, long double valY0, long double valX0, long double valStepInt, long double valKNav
, TShipTarg ShipTarg ,long double *arrCfW, TGlonass Glonass, long double &valProbab )
{


	int nr = 0 ; // ноер текущей реализации метода стат испытаний
	const int  iControlType = 10;
	int iArrParams[5] = {0};
	long double DblArrParams[20] ={0.};
	DblArrParams[0] = valKNav; // КН
	TCntlFuncPar CntlFuncPar ( iControlType, iArrParams , DblArrParams);

  valProbab = 0.;
 for ( nr = 0; nr < NR; nr++)
 {

 // инициализация УАБ

	  // создание АУБ с разбросо параметров
	  TBomb BombReal = Bomb;
	  BombReal.mMidSq += TShipTraj::getGauss(0., 1. ) * BombReal.mMidSq * 0.01;
	  BombReal.mSurfFuzSq += TShipTraj::getGauss(0., 1. ) * BombReal.mSurfFuzSq * 0.01;
	  BombReal.mMass += TShipTraj::getGauss(0., 1. ) * BombReal.mMass * 0.01;
	  BombReal.mLamb_k += TShipTraj::getGauss(0., 1. ) * BombReal.mLamb_k * 0.01;
	  BombReal.mSq_k += TShipTraj::getGauss(0., 1. ) * BombReal.mSq_k * 0.01;
	  TBombTraj BombTraj( valTet0,valV0 ,valY0 ,valX0,valStepInt ,valTBegin, BombReal, arrCfW	);


	  THomingSituation  HomingSituation( ShipTarg , BombTraj , CntlFuncPar, Glonass	,valTBegin, NULL ) ;

	 long double valProbabCur = -1.;
	 long  double valSig_AngVis = -1.;
	 bool breturn = HomingSituation.fncMoveHomeHead_TO_BlastPoint_AND_ShowGraphs_StatMeth( NULL, Bomb, valProbabCur, valSig_AngVis );
	 if ( breturn)
	 {
	   valProbab +=  valProbabCur;
	 }

 }
 valProbab = valProbab/(long double)NR ;

}


// создание картинок привязанных к текущей точке траектории - УАБ, вектора ее скорости, угла AlfaT(простр угла атаки)
// pwchShapeUAB -  путь к директоии с шейп файлами УАБ
// wcharrPath - путь к папке с результирующими графиками
//void THomingSituation:: drawPictures(wchar_t * pwchShapeUAB, wchar_t *pwcharrPath, const double valRastigenie)
void THomingSituation:: drawPictures( wchar_t *pwcharrPath, const double valRastigenie)
{

  // 1. картиеки с УАБ
	wchar_t wchFileName[300] = {0};

	// соэдание путей к результирующему файлу  с УАБ
  wchar_t wchFileNameRez0[300] = {0};
	wcscpy(wchFileNameRez0, pwcharrPath );
	wcscat(wchFileNameRez0, L"\\UABUnionRez00.shp");


	 wchar_t wchFileNameRez10[300] = {0} ;

	wcscpy(wchFileNameRez10, pwcharrPath );
	wcscat(wchFileNameRez10, L"\\UABUnionRez10.shp");



   TURPointXY pntSdvig( mBombTraj.marrStrSK_VS [0], mBombTraj.marrStrSK_VS [1]);

	 drawUAB( wchFileNameRez0,mBombTraj.marrStrSK_VS [4],  pntSdvig,valRastigenie) ;
	 long double valT = (mTCur > 0.01)? mTCur -0.001:mTCur ;
	long double valU = -fncCalcU(valT);
   //	TURPolygon plg1 =  plgBomb.LinTransform(mBombTraj.marrStrSK_VS [4] + M_PI +valU,  pntSdvig, valRastigenie ) ;
   //	plg1.WriteSetSHPFiles(wchFileNameRez10, &plg1,1) ;
	 drawUAB( wchFileNameRez10,mBombTraj.marrStrSK_VS [4] +valU,  pntSdvig,valRastigenie) ;

  // 2. картинки с продольной осью УАБ и касательной к траектории
  double valCoeff = 500.;
  TURPointXY arrPnt[2];
  // точка текущей траектории
  arrPnt[0].X = mBombTraj.marrStrSK_VS[0];
  arrPnt[0].Y = mBombTraj.marrStrSK_VS[1];
  ///

  // точка на касательной к траектории УАБ
   arrPnt[1].X = mBombTraj.marrStrSK_VS[0] + valCoeff * valRastigenie * cos(mBombTraj.marrStrSK_VS[4]);
	arrPnt[1].Y = mBombTraj.marrStrSK_VS[1] + valCoeff * valRastigenie * sin(mBombTraj.marrStrSK_VS[4]);
   ///

   TURPolyLine PlnTang(arrPnt,2) ;
	wchar_t wchFileNamePlnTang[300] = {0}, wchFileNamePlnAxe[300] = {0}, wchFileNameSect0[300] = {0},wchFileNameSect1[300] = {0} ;

	wcscpy(wchFileNamePlnTang, pwcharrPath );
	wcscat(wchFileNamePlnTang, L"\\CasatLine.shp");
	PlnTang.WriteSetSHPFiles(wchFileNamePlnTang, &PlnTang, 1) ;
	///

	// точка на оси УАБ
   arrPnt[1].X = mBombTraj.marrStrSK_VS[0] + valCoeff * valRastigenie * cos(mBombTraj.marrStrSK_VS[4] + valU);
   arrPnt[1].Y = mBombTraj.marrStrSK_VS[1] + valCoeff * valRastigenie * sin(mBombTraj.marrStrSK_VS[4] + valU);
	TURPolyLine PlnAxe(arrPnt,2) ;
   wcscpy(wchFileNamePlnAxe, pwcharrPath );
   wcscat(wchFileNamePlnAxe, L"\\AxeLine.shp");
   PlnTang.WriteSetSHPFiles(wchFileNamePlnAxe, &PlnAxe, 1) ;
   ///

   // секторы для обозначения угла атаки
  TURPolyLine plnSect0 = 	TURPolyLine :: fncCreateSector( arrPnt[0], valCoeff * valRastigenie /2.,
					mBombTraj.marrStrSK_VS[4],mBombTraj.marrStrSK_VS[4] + valU,400) ;
  wcscpy( wchFileNameSect0, pwcharrPath );
   wcscat( wchFileNameSect0, L"\\Sect0.shp");
   plnSect0 .WriteSetSHPFiles( wchFileNameSect0, &plnSect0 , 1) ;

   TURPolyLine plnSect1 = 	TURPolyLine :: fncCreateSector( arrPnt[0], valCoeff * valRastigenie /2.* 0.95,
					mBombTraj.marrStrSK_VS[4],mBombTraj.marrStrSK_VS[4] + valU,400) ;
  wcscpy( wchFileNameSect1, pwcharrPath );
   wcscat( wchFileNameSect1, L"\\Sect1.shp");
   plnSect1 .WriteSetSHPFiles( wchFileNameSect1, &plnSect1 , 1) ;
   ///

   // картинка с сектором обзора координатора ГСН , ось УАБ совпадаеть с вектром скорости
   const double valR = mBombTraj.mBomb.mHomingHead.mR ;
   double valFi0 = mBombTraj.marrStrSK_VS[4] - mBombTraj.mBomb.mHomingHead.mFi/2.;
   double valFi1 = mBombTraj.marrStrSK_VS[4] + mBombTraj.mBomb.mHomingHead.mFi/2.;
   TURPolygon PlgSectorHomeHead0 = TURPolygon::fncCreateSector(arrPnt[0], valR,valFi0,valFi1,400) ;
   wchar_t wchFileNameSectorHomeHead0[300] = {0};
   wcscpy(wchFileNameSectorHomeHead0, pwcharrPath );
   wcscat(wchFileNameSectorHomeHead0, L"\\SectotHH0.shp");
   PlgSectorHomeHead0.WriteSetSHPFiles(wchFileNameSectorHomeHead0, &PlgSectorHomeHead0, 1) ;

   valFi0 = mBombTraj.marrStrSK_VS[4] + valU - mBombTraj.mBomb.mHomingHead.mFi/2.;
   valFi1 = mBombTraj.marrStrSK_VS[4] + valU + mBombTraj.mBomb.mHomingHead.mFi/2.;
   TURPolygon PlgSectorHomeHead1 = TURPolygon::fncCreateSector(arrPnt[0], valR,valFi0,valFi1,400) ;
   wchar_t wchFileNameSectorHomeHead1[300] = {0};
   wcscpy(wchFileNameSectorHomeHead1, pwcharrPath );
   wcscat(wchFileNameSectorHomeHead1, L"\\SectotHH1.shp");
   PlgSectorHomeHead1.WriteSetSHPFiles(wchFileNameSectorHomeHead1, &PlgSectorHomeHead1, 1) ;
   ///

}
 //---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void THomingSituation::createPicturesConusZahv_and_ConusHoming(wchar_t *pwchOutFile)
{
	   wchar_t wcharrPath[300] = {0} ;
		wcscpy(wcharrPath, pwchOutFile );
		wcscat(wcharrPath, L"\\PictConusZahv_And_TPodl");
		_wmkdir(wcharrPath);

		createPicturesConusZahv(wcharrPath);

	   //	createGraph_f(wcharrPath, 270.);
	   createGraph_f(wcharrPath, 270.);

		createPicturesConusHoming(wcharrPath);

		wchar_t wchFileName4[300] = {0} ;
		wcscpy(wchFileName4, wcharrPath );
		wcscat(wchFileName4, L"\\Axes.shp");
		TYrWriteShapeFile::CreateShpAxes(wchFileName4,-80000.,80000,-80000.,80000) ;
		//// построение эллипсов
		wcscpy(wchFileName4, wcharrPath );
		wcscat(wchFileName4, L"\\Ellips");
		_wmkdir( wchFileName4);
		createEllipsGraph(wchFileName4) ;

		//TURPointXY pntCentre( mBombTraj.marrStrSK_VS[0] + mBombTraj.marrStrSK_VS[1]/fabs(tan(mBombTraj.marrStrSK_VS[4])),0.);
	   //	wcscpy(wchFileName4, wcharrPath );
	   //	wcscat(wchFileName4, L"\\PntCentre.shp");
	   //	TURPointXY::WriteSetSHPFiles(wchFileName4, &pntCentre, 1) ;
		// tetta =30
		THomingSituation HomSit = *this;
		double valTet = M_PI /180. * 30. ;
		double valFi = (double)mBombTraj.mBomb.mHomingHead.mFi ;
		double valD  = (double)mBombTraj.mBomb.mHomingHead.mR ;
		double alf = valFi/2.;
		double gam =  valTet -alf ;
		double height = valD * sin(gam);
		double delx = height/ tan( valTet);
		HomSit.mBombTraj = TBombTraj (
	-valTet  // нач угол бросания
	,mBombTraj.marrStrSK_VS[3]  // нач скоргсть
	,height  // нач высота
	 ,mBombTraj.marrStrSK_VS[0] + mBombTraj.marrStrSK_VS[1]/fabs(tan(mBombTraj.marrStrSK_VS[4])) - delx  // нач. положение по оси X
	,mBombTraj.mStepInt  // шагт интьерирования по времени диф уравнений движения
	 ,mBombTraj.mTCur
	,mBombTraj.mBomb
	) ;
	wchar_t pwcharrPath1[300] ={0};
	wcscpy(pwcharrPath1, wcharrPath );
	wcscat(pwcharrPath1, L"\\tetta=30");
	_wmkdir (pwcharrPath1);
	const double valRastigenie = 130.;
	HomSit.drawPictures(pwcharrPath1,  valRastigenie);
	wcscpy(wchFileName4, wcharrPath );
	wcscat(wchFileName4, L"\\Ellips_30");
	_wmkdir( wchFileName4);
	HomSit.createEllipsGraph(wchFileName4) ;

	//
		// tette =15
		HomSit = *this;
		valTet = M_PI /180. * 15. ;
		gam =  valTet -alf ;
		height = valD * sin(gam);
		delx = height/ tan( valTet);
		HomSit.mBombTraj = TBombTraj (
	-valTet  // нач угол бросания
	,mBombTraj.marrStrSK_VS[3]  // нач скоргсть
	,height  // нач высота
	 ,mBombTraj.marrStrSK_VS[0] + mBombTraj.marrStrSK_VS[1]/fabs(tan(mBombTraj.marrStrSK_VS[4])) - delx  // нач. положение по оси X
	,mBombTraj.mStepInt  // шагт интьерирования по времени диф уравнений движения
	 ,mBombTraj.mTCur
	,mBombTraj.mBomb
	) ;

	wcscpy(pwcharrPath1, wcharrPath );
	wcscat(pwcharrPath1, L"\\tetta=15");
	_wmkdir (pwcharrPath1);
	HomSit.drawPictures( pwcharrPath1,  valRastigenie);
	wcscpy(wchFileName4, wcharrPath );
	wcscat(wchFileName4, L"\\Ellips_15");
	_wmkdir(wchFileName4);
	HomSit.createEllipsGraph(wchFileName4) ;
	///

	//
		// tette =60
		HomSit = *this;
		valTet = M_PI /180. * 60. ;
		gam =  valTet -alf ;
		height = valD * sin(gam);
		delx = height/ tan( valTet);
		HomSit.mBombTraj = TBombTraj (
	-valTet  // нач угол бросания
	,mBombTraj.marrStrSK_VS[3]  // нач скоргсть
	,height  // нач высота
	 ,mBombTraj.marrStrSK_VS[0] + mBombTraj.marrStrSK_VS[1]/fabs(tan(mBombTraj.marrStrSK_VS[4])) - delx  // нач. положение по оси X
	,mBombTraj.mStepInt  // шагт интьерирования по времени диф уравнений движения
	 ,mBombTraj.mTCur
	,mBombTraj.mBomb
	) ;
	wcscpy(pwcharrPath1, wcharrPath );
	wcscat(pwcharrPath1, L"\\tetta=60");
	_wmkdir (pwcharrPath1);
	HomSit.drawPictures( pwcharrPath1,  valRastigenie);
	wcscpy(wchFileName4, wcharrPath );
	wcscat(wchFileName4, L"\\Ellips_60");
	_wmkdir(wchFileName4);

	HomSit.createEllipsGraph(wchFileName4) ;
	///

	//
		// tette =90
		HomSit = *this;
		 valTet = M_PI /2. -0.01 ;

		 gam =  valTet -alf ;
		 height = valD * sin(gam);
		delx = height/ tan( valTet);
		HomSit.mBombTraj = TBombTraj (
	-valTet  // нач угол бросания
	,mBombTraj.marrStrSK_VS[3]  // нач скоргсть
	,height  // нач высота
	 ,mBombTraj.marrStrSK_VS[0] + mBombTraj.marrStrSK_VS[1]/fabs(tan(mBombTraj.marrStrSK_VS[4])) - delx  // нач. положение по оси X
	,mBombTraj.mStepInt  // шагт интьерирования по времени диф уравнений движения
	 ,mBombTraj.mTCur
	,mBombTraj.mBomb
	) ;

	wcscpy(pwcharrPath1, wcharrPath );
	wcscat(pwcharrPath1, L"\\tetta=90");
	_wmkdir (pwcharrPath1);

	HomSit.drawPictures( pwcharrPath1,  valRastigenie);
	wcscpy(wchFileName4, wcharrPath );
	wcscat(wchFileName4, L"\\Ellips_90");
	_wmkdir(wchFileName4);

	HomSit.createEllipsGraph(wchFileName4) ;

    	wcscpy(wchFileName4, wcharrPath );
		wcscat(wchFileName4, L"\\VesselCur.shp");
	double valPos = mBombTraj.marrStrSK_VS[0] + mBombTraj.marrStrSK_VS[1]/fabs(tan(mBombTraj.marrStrSK_VS[4])) ;
	double valRastigenieVess = 6.;
	drawVessel( wchFileName4, valPos, mShipTarg.mTraject.marrVectSostGSK[3], valRastigenieVess ) ;

}


void  THomingSituation::createEllipsGraph(wchar_t *wchFoldPath)
{
//эллипс лдля тетта = 45 град
		double valFi = (double)mBombTraj.mBomb.mHomingHead.mFi ;
		double valD  = (double)mBombTraj.mBomb.mHomingHead.mR ;
		double alf = valFi/2.;
		double valTet = fabs( mBombTraj. marrStrSK_VS [4] );
		double a = fnc_a( valTet);
		double b = fnc_b( valTet);
		// центр эллипса
		double x2 = valD * cos( valTet - alf );
		double x1 = valD * sin( valTet - alf )/tan( valTet + alf );

		 TURPointXY pntCentre((x1+x2)/2.+ mBombTraj.marrStrSK_VS[0] ,0.);
		TURPolygon plgEll= TURPolygon::createEllips(a, b , pntCentre, 300 ) ;
		wchar_t wchFileName[300] = {0};
		wcscpy(wchFileName, wchFoldPath );
		wcscat(wchFileName, L"\\Ellips.shp");
		TURPolygon::WriteSetSHPFiles(wchFileName, &plgEll, 1) ;
		wcscpy(wchFileName, wchFoldPath );
		wcscat(wchFileName, L"\\Centre.shp");
		TURPointXY::WriteSetSHPFiles(wchFileName, &pntCentre, 1) ;

		TURPointXY pntAim(  mBombTraj.marrStrSK_VS[0] + mBombTraj.marrStrSK_VS[1]/fabs(tan(mBombTraj.marrStrSK_VS[4])),0.);
		wcscpy(wchFileName, wchFoldPath );
		wcscat(wchFileName, L"\\pntAim.shp");
		TURPointXY::WriteSetSHPFiles(wchFileName, &pntAim, 1) ;

		TURPolygon plgTargCircle = TURPolygon::fncCreateCircle(pntAim, b, 300) ;
		wcscpy(wchFileName, wchFoldPath );
		wcscat(wchFileName, L"\\TargCircle.shp");

		TURPolygon::WriteSetSHPFiles(wchFileName, &plgTargCircle, 1) ;

}
void THomingSituation::createGraph_f(wchar_t *pwchOutFile , double valV)
{
  TURPointXY pnts[99],pnts1[99];
  double valFi = (double)mBombTraj.mBomb.mHomingHead.mFi ;
  double valD  = (double)mBombTraj.mBomb.mHomingHead.mR ;
  double delTetta = (M_PI/2. - valFi)/100.;
  double Cy = mBombTraj.mBomb.fncCy(1.);//    4.5;
   double almax =  mBombTraj.mBomb.mMaxAt;// 0.34;
   double Sm =  mBombTraj.mBomb.mMidSq;// 0.2642;
   double mass =  mBombTraj.mBomb.mMass;//  2200.;
  for (int i = 0; i < 99; i++)
  {
   double valF = 0., valDF =0.;
   double valTet = valFi + delTetta * ((double)(i + 1.));
   fnc_f_and_df(   valV
	 , valTet ,valF,valDF);
   pnts[i].X = valTet * 10. ;
   pnts[i].Y = valF * 10.;
   pnts1[i].X = valTet * 10. ;
   pnts1[i].Y = valDF* 10.;

  }
   wchar_t wcharrFile[300] = {0} ;
   wcscpy(wcharrFile, pwchOutFile );
   wcscat(wcharrFile, L"\\f.shp");
  TURPolyLine pln(pnts,99) ;
  TURPolyLine::WriteSetSHPFiles(wcharrFile, &pln, 1) ;

  wcscpy(wcharrFile, pwchOutFile );
  wcscat(wcharrFile, L"\\df.shp");
  TURPolyLine pln1(pnts1,99) ;
  TURPolyLine::WriteSetSHPFiles(wcharrFile, &pln1, 1) ;


}
void THomingSituation::createPicturesConusZahv(wchar_t *pwchOutFile)
{

  double valFi = (double)mBombTraj.mBomb.mHomingHead.mFi ;
  double valD  = (double)mBombTraj.mBomb.mHomingHead.mR ;

  TURPointXY pnts[99],pnts1[99],pnts2[99];
  double delTetta = (M_PI/2. - valFi)/100.;
  for (int i = 0; i < 99; i++)
  {
   double valTet = valFi + delTetta * ((double)(i + 1.));
   pnts[i].X = valTet * 180./M_PI ;
   pnts[i].Y = fnc_b( valTet)/20.;
   pnts1[i].X = valTet * 1000. ;
   pnts1[i].Y = fnc_a( valTet);
   pnts2[i].X = valTet * 1000. ;
   pnts2[i].Y = fnc_b( valTet);
  }
   wchar_t wcharrFile[300] = {0} ;
   wcscpy(wcharrFile, pwchOutFile );
   wcscat(wcharrFile, L"\\Graph_TPodleta_from_Tetta.shp");
  TURPolyLine pln(pnts,99) ;
  TURPolyLine::WriteSetSHPFiles(wcharrFile, &pln, 1) ;

  wcscpy(wcharrFile, pwchOutFile );
  wcscat(wcharrFile, L"\\a.shp");
  TURPolyLine pln1(pnts1,99) ;
  TURPolyLine::WriteSetSHPFiles(wcharrFile, &pln1, 1) ;

  wcscpy(wcharrFile, pwchOutFile );
  wcscat(wcharrFile, L"\\b.shp");
  TURPolyLine pln2(pnts2,99) ;
  TURPolyLine::WriteSetSHPFiles(wcharrFile, &pln2, 1) ;



}

void THomingSituation::createPicturesConusHoming(wchar_t *pwchOutFile)
{
  double valFi = (double)mBombTraj.mBomb.mHomingHead.mFi ;
  double valD  = (double)mBombTraj.mBomb.mHomingHead.mR ;

  TURPointXY pnts[100];
  double delV = (325. - 50.)/100.;
  for (int i = 0; i < 100; i++)
  {
   double valV = 50. + delV * ((double)i );
   pnts[i].X = valV ;
   pnts[i].Y = fnc_tetta( valV)*180./ M_PI;
  }
   wchar_t wcharrFile[300] = {0} ;
   wcscpy(wcharrFile, pwchOutFile );
   wcscat(wcharrFile, L"\\Graph_TetteGran_from_V.shp");
  TURPolyLine pln(pnts,100) ;
  TURPolyLine::WriteSetSHPFiles(wcharrFile, &pln, 1) ;

}
double     THomingSituation::fnc_b(double valTet)
{
  double valFi = (double)mBombTraj.mBomb.mHomingHead.mFi ;
  double valR  = (double)mBombTraj.mBomb.mHomingHead.mR ;

	double alf = valFi/2.;
	double valt = -(sin(alf)/ tan( valTet) - cos(alf));
   //	double valt1 = (1. -  tan(alf)* tan(alf)/ valt/ valt);
	double valt1 = 1. -  tan(alf)* tan(alf)/ tan( valTet)/ tan( valTet);
	double valx0 = tan(alf)* tan(alf) *  valR * valt/ tan (valTet) / valt1;
	double valbsq = - valt1 * valx0 * valx0 + 2. * valx0 *tan(alf)* tan(alf) /  tan (valTet)* valR * valt
		 +tan(alf)* tan(alf) * valt*valt * valR * valR;
	 if (valbsq <0.)
	 {
	 ShowMessage("Err fnc_b");
	 }

	return sqrt(valbsq);
}


double     THomingSituation::fnc_a( double valTet)
{
  double valFi = (double)mBombTraj.mBomb.mHomingHead.mFi ;
  double valR  = (double)mBombTraj.mBomb.mHomingHead.mR ;

	double alf = valFi/2.;
	double x1 = valR * cos ( valTet - alf);
	double x2 = valR * sin ( valTet - alf) * cos( valTet + alf)/sin( valTet + alf);

	return (x1 -x2)/2.;
}

double     THomingSituation::fnc_tetta( double valV)
{
  double valFi = (double)mBombTraj.mBomb.mHomingHead.mFi ;
  double valD  = (double)mBombTraj.mBomb.mHomingHead.mR ;
  double alf = valFi/2.;

  double valTet0 = M_PI/2.-0.05;

   double Cy = (double)mBombTraj.mBomb.fncCy( 0.);
   double almax = (double)mBombTraj.mBomb.mMaxAt ;
   double Sm = (double)mBombTraj.mBomb.mMidSq;
   double mass = (double)mBombTraj.mBomb.mMass;  for (int i = 0; i < 20; i++)
  {
   double valF = 0., valDF =0. ;
   fnc_f_and_df( valV , valTet0  ,valF, valDF);
   double del = - valF / valDF;
   valTet0 += del;
   if(fabs(del) < 0.001)
   {
	   break;
   }

  }
  return valTet0 ;
}

void     THomingSituation::fnc_f_and_df( double valV
	  ,double valTet0  ,double &valF,double &valDF)
{
  double valFi = (double)mBombTraj.mBomb.mHomingHead.mFi ;
  double valD  = (double)mBombTraj.mBomb.mHomingHead.mR ;
  double alf = valFi/2.;

   double valH =  valD * sin( valTet0 - alf);
   long double valGradTayl = 0., valTayl = 0.;
   fncCalcNormTemperature((long double)valH,valTayl, valGradTayl)  ;
   double pont = fncCalcOtnositP( valH);
   double dpont = -pont *(valGradTayl + G_ZEMLI/ATM_R_UNIVER )/ valTayl;

   double Cy = (double)mBombTraj.mBomb.fncCy( 0.);
   double almax = (double)mBombTraj.mBomb.mMaxAt ;
   double Sm = (double)mBombTraj.mBomb.mMidSq;
   double mass = (double)mBombTraj.mBomb.mMass;

	valF = cos( valTet0) - Cy *  Sm * almax/ mass /G_ZEMLI *  valV * valV/2. * ATM_RoN0 * pont;
	valDF =  -sin(valTet0)- Cy *  Sm * almax/ mass /G_ZEMLI *  valV * valV/2. * ATM_RoN0 *valD * cos(valTet0 - alf) *dpont;
}





long double THomingSituation::fncAngVisir()
{

	 long double arrVS_Ship_GSK[4] = {0.};
	 arrVS_Ship_GSK[0 ] =  mShipTarg.mTraject.marrVectSostGSK[0] ;
	 arrVS_Ship_GSK[1 ] =  mShipTarg.mTraject.marrVectSostGSK[2] ;
	 arrVS_Ship_GSK[2] =   mShipTarg.mTraject.marrVectSostGSK[3] ;
	 arrVS_Ship_GSK[3 ] =  mShipTarg.mTraject.marrVectSostGSK[5] ;


	 long double arrVS_Bomb_GSK[4] = {0.};
	 arrVS_Bomb_GSK[0 ] =  mBombTraj.marrStrSK_VS[0] ;
	 arrVS_Bomb_GSK[1 ] =  mBombTraj.marrStrSK_VS[1] ;
	 arrVS_Bomb_GSK[2] =   mBombTraj.marrStrSK_VS [3] * cosl(mBombTraj.marrStrSK_VS [4]) ;
	 arrVS_Bomb_GSK[3 ] =  mBombTraj.marrStrSK_VS [3] * sinl(mBombTraj.marrStrSK_VS [4]) ;



   long double arrDelS[2] = {0.};
  MtrxMinusMatrx( arrVS_Ship_GSK, arrVS_Bomb_GSK,2, 1, arrDelS);

  long double temp =OuterProduct_2(&(arrVS_Ship_GSK[2]), arrDelS    )/NormVect2(arrDelS)/NormVect2(&(arrVS_Ship_GSK[2]) );
  if (fabs(temp) > 1.)
  {
   //	ShowMessage(L"ERROR2");
   temp = temp / fncSign (temp);;
  }
  temp = asin(temp);
  if (arrVS_Ship_GSK[2]< 0.)
  {
    temp += M_PI;
  }

  //long double t1 = fncCalcDeltaAngVisir();
  return  temp;
 // long double t = atan2l(arrDelS[1],arrDelS[0]);
//  return  t;


}
long double THomingSituation::fncDist()
{

	long double arrVS_Ship_GSK[2] = {0.};
	arrVS_Ship_GSK[0 ] =  mShipTarg.mTraject.marrVectSostGSK[0] ;
	arrVS_Ship_GSK[1 ] =  mShipTarg.mTraject.marrVectSostGSK[2] ;
	long double arrDelS[2] = {0.};
	MtrxMinusMatrx( arrVS_Ship_GSK, mBombTraj.marrStrSK_VS,2, 1, arrDelS);
	return  NormVect2(arrDelS);

}


void  THomingSituation:: ImitateZamer(long double &valRZv,long double &valFiVisZv,long double &valTZv
	 , long double &valSig_BMO_R, long double &valSig_MMO_R,
	   long double &valSig_BMO_Fi, long double &valSig_MMO_Fi)
{
  valSig_BMO_R = mBombTraj.mBomb.mHomingHead.mSigFl_R / sqrtl( mBombTraj.mStepInt/mBombTraj.mBomb.mHomingHead.m_h -1. );
  valSig_MMO_R =  mShipTarg.mTargData.mLTarg/ 6.;

 // valSig_BMO_R = sqrtl( valSig_BMO_R * valSig_BMO_R + valSig_MMO_R  * valSig_MMO_R ); //
 // valSig_MMO_R = 0.;

  valSig_BMO_Fi = mBombTraj.mBomb.mHomingHead.mSigFl_U/ sqrtl( mBombTraj.mStepInt/mBombTraj.mBomb.mHomingHead.m_h -1. );
  valSig_MMO_Fi = (atanl(mBombTraj.marrStrSK_VS [1]/ (mShipTarg.mTraject.marrVectSostGSK[0]
				   - mShipTarg.mTargData.mLTarg/ 2. -mBombTraj.marrStrSK_VS [0]))
				  -atanl(mBombTraj.marrStrSK_VS [1]/ (mShipTarg.mTraject.marrVectSostGSK[0]
				   + mShipTarg.mTargData.mLTarg/ 2. -mBombTraj.marrStrSK_VS [0]))) /6.;// / sqrtl( mBombTraj.mStepInt/mBombTraj.mBomb.mHomingHead.m_h -1. );

  valSig_BMO_Fi = sqrtl(valSig_BMO_Fi * valSig_BMO_Fi +  valSig_MMO_Fi * valSig_MMO_Fi/ 2. );
  valSig_MMO_Fi =  valSig_MMO_Fi/ sqrtl(2.);
  long double d = fncDist();
  long double Fi= fncAngVisir();
  long double errR=  TShipTraj::getGauss(0., sqrtl( valSig_BMO_R * valSig_BMO_R + valSig_MMO_R * valSig_MMO_R) );
  long double errFi=  TShipTraj::getGauss(0., sqrtl(valSig_BMO_Fi*valSig_BMO_Fi +valSig_MMO_Fi*valSig_MMO_Fi) );
  valFiVisZv = Fi + errFi;
  valRZv =  d +errR;
  valTZv = mBombTraj.mTCur;
}



// нахождение  угла бросания обеспечивающего макс дальночть полеты УАБ при постоянной начальной скорости и высоте бросания
// методом перебора и построение графика зависимотсти иакс даьности от угла бросания
// если pwchMaxDistOutFile== NULL, то графики не строятся
// INPUT:
// pwchMaxDistOutFile - путь к папке с отчетом
// valAngStep - шаг по углу, начинается с valAngStep
// quantSteps - к-во шагов по углу
// OUTPUT:
// вся информация об оптимальной траектории содержится в return
THomingSituation THomingSituation::findOptAng_for_MaxDist_AndShowGraphs(wchar_t *pwchMaxDistOutFile, const long double valAngStep, const int quantSteps)
{
	long double valDHorizPlanOpt = -1., valD = -1.;
	THomingSituation HomingSituationReturn = *this;
	THomingSituation HomingSituation1 = HomingSituationReturn ;
	HomingSituation1.fncMoveBomb_TO_ZeroAlt_AND_ShowGraphs( NULL, valDHorizPlanOpt );
	TURPointXY * pPoints = NULL;
	if(pwchMaxDistOutFile)
	{
	 pPoints = new TURPointXY[ quantSteps];
	}
	for (int i = 0; i <  quantSteps; i++)
	{
		long double valTet0 = valAngStep * ((long double)i);
		TBombTraj BombTraj (valTet0,mBombTraj.mV0 ,mBombTraj.mAltit,mBombTraj.mX0,mBombTraj.mStepInt  // шагт интьерирования по времени диф уравнений движения
		,mBombTraj.mTStart,mBombTraj.mBomb);

		THomingSituation HomingSituation(mShipTarg, BombTraj, mCntlFuncPar,mTCur,NULL);

		THomingSituation HomingSituationPlan = HomingSituation;

		bool brez1 = HomingSituationPlan.fncSolveOCP_for_MaxDist_MethNewton(HomingSituation.mCntlFuncPar);

		long double valDHorizPlan = -1.;
			THomingSituation HomingSituation2 = HomingSituation;
		HomingSituation2.fncMoveBomb_TO_ZeroAlt_AND_ShowGraphs( NULL, valDHorizPlan );
		if(pwchMaxDistOutFile)
		{
		  pPoints[i].X = 1000.* valTet0 * 180./M_PI;
		  pPoints[i].Y =  valDHorizPlan ;

        }
		if (valDHorizPlan > valDHorizPlanOpt)
		{
		valDHorizPlanOpt = valDHorizPlan ;
		HomingSituationReturn = HomingSituation;
		}

	}
	if(pwchMaxDistOutFile)
	{
	TURPolyLine  PolyLine( pPoints,quantSteps) ;
	wchar_t wchFileName[300]  = {0};
	wcscpy(wchFileName,pwchMaxDistOutFile );
	wcscat(wchFileName, L"\\MaxDist_fromTetta.shp");

	PolyLine.WriteSetSHPFiles(wchFileName,&PolyLine, 1) ;
	THomingSituation HomingSituation2 = HomingSituationReturn ;
	HomingSituation2.fncMoveBomb_TO_ZeroAlt_AND_ShowGraphs(pwchMaxDistOutFile,valD) ;
	delete []pPoints;
	}
	return  HomingSituationReturn ;
}


// Построениме графика зависимости максимальной дальности от высоты бросания
// Для каждой высоты valHMin + i * valHStep, i=0,..., (quantStep-1) расчитывается максимальнная дальность полета,
// запоминаетися и строится график
// Очевидно, чьто максимальная дальность valMaxDist достигается на аксиальной высоте бросания
// INPUT:
// valV0 - нач скорость бросания
// valStepInt - шаг интегрирован7ия по времени методо Эйлера
// pwchMaxDistOutFile  - директория с графиком
// valHMin   - миниальная высота с которой начинает строится график
// valHStep  - шаг по высоте
//quantStep   - к-во шагов
// OUTPUT:
// valMaxDist  - максимальная дальность полета (на максимальной высоте бросания)
// valTPerecl  - момент переключения управления (угла атаки)
// valAng  - оптимальный угол бросания
void THomingSituation::createGraph_OptMaxDist_fom_Height(long double valV0,long double valStepInt, wchar_t *pwchMaxDistOutFile
	, long double valHMin,long double valHStep,int quantSteps, long double &valMaxDist, long double &valTPerecl, long double &valAng)
{
    const double valTBegin = 0.;
	  const double Bearing = M_PI/ 2., TargCourse = M_PI/ 2.,TargZenitAng = M_PI/ 2., valDist0 = 15000., valH = 0., valV = 20.;
	 TInitTargData InitData( Bearing,  TargCourse,  TargZenitAng , valV
	 ,valDist0, valH, valTBegin);

	 const int quantParts = 1;
	 TPartData arrPartData[1];
	 arrPartData[0].miTypePart = 1;
	 arrPartData[0].mTimePart = 1000.;
	 arrPartData[0].mSigW = 10.;


	 TShipTarg ShipTarg(InitData ,quantParts, arrPartData, NULL);
	  long double  valDHoriz = -1. ;
	  TBomb Bomb;
	  TURPointXY * pPoints =  new TURPointXY[ quantSteps];

	 for (int i = 0; i < quantSteps ; i++)
	 {
		long double valY0 =  valHMin + ((long double)i) * valHStep;
		long double valTet0 = 0., valTBegin = 0.;
		TBombTraj BombTraj (valTet0	,valV0, valY0 ,0., valStepInt, valTBegin	, Bomb 	) ;
		const int  iControlType = 1;
		int iArrParams[5] = {0};
		iArrParams[0] = 3;
		long double DblArrParams[20] ={0.};
		DblArrParams[1] = 13.;
		DblArrParams[2] = 1000.;
		DblArrParams[3] = 1000.;
		TCntlFuncPar CntlFuncPar ( iControlType, iArrParams , DblArrParams);
		THomingSituation  HomingSituation( ShipTarg , BombTraj , CntlFuncPar	,0., NULL ) ;
		THomingSituation HomingSit = HomingSituation.findOptAng_for_MaxDist_AndShowGraphs(NULL, 3./180.*M_PI, 10) ;
		if (i == (quantSteps-1))
		{
			HomingSit.fncMoveBomb_TO_ZeroAlt_AND_ShowGraphs(pwchMaxDistOutFile,valDHoriz) ;
			valMaxDist =  valDHoriz;
			valTPerecl = HomingSit.mCntlFuncPar.mDblArrParams[1];
			valAng =  HomingSit.mBombTraj.mTet0;
		}
		else
		{
		HomingSit.fncMoveBomb_TO_ZeroAlt_AND_ShowGraphs( NULL, valDHoriz );
        }
		pPoints[i].X = valY0;
		pPoints[i].Y =  valDHoriz ;


	 }
	TURPolyLine  PolyLine( pPoints,quantSteps) ;
	wchar_t wchFileName[300]  = {0};
	wcscpy(wchFileName,pwchMaxDistOutFile );
	wcscat(wchFileName, L"\\MaxDist_fromHeight.shp");

	PolyLine.WriteSetSHPFiles(wchFileName,&PolyLine, 1) ;
	delete []pPoints;

}





// Построениме графика зависимости максимальной дальности перехвата цели от высоты бросания
// Для каждой высоты valHMin + i * valHSteps, i=0,..., (quantSteps-1) расчитывается максимальная дальность перехвата,
// запоминаетися и строится график
// считается, что цель может быть перехвачена на некоторой дальности, если найдется такая траектория подлета к точке исходного положения цели ( в момент "0")
// что в момент времени, когда расстояние от УАБ до исходной точки будет равно дальности открытия ГСН и цель при любом возможном движении
// окажется в секторе обзора ГСН
// INPUT:
// valV0 - нач скорость бросания
// valStepInt - шаг интегрирован7ия по времени методом Эйлера
// Bomb - данные УАБ
// QuantOverSwitchPnts - к-во точек перключения
// pwchMaxDistZahvOutFile  - директория с графиком
// valHMin   - миниальная высота с которой начинает строится график
// valHStep  - шаг по высоте
// quantSteps   - к-во шагов
// OUTPUT:
// valHeight
// iTypeOfTraj
// valMaxDistTakeover
// valTetta0
// valTPerecl
// valTZahv
// valVXahv
// valYZahv
// valXZahv
// valTettaZahv

bool THomingSituation::createGraph_MaxDist_For_Takeover_from_Height(long double valV0,long double valStepInt, TBomb Bomb
	 ,wchar_t *pwchMaxDistZahvOutFile, int QuantOverSwitchPnts	, long double valHMin
	 ,long double valHStep,int quantSteps, long double &valHeight,int &iTypeOfTraj
	 , long double &valMaxDistTakeover,long double &valTetta0, long double &valTPerecl
	 , long double &valTZahv, long double &valVXahv, long double &valYZahv,long double & valXZahv
	 , long double &valTettaZahv)
{

	  bool breturn = false;
	  valMaxDistTakeover = -100.;
	  // создание цели. Цель в расчетах не участвует, но создать какую-то
	  // требуется для того, чтобы создать объект  класса THomingSituation
	  double TargDist = 20000.;
	  const double valTBegin = 0.;
	  const double Bearing = M_PI/ 2., TargCourse = M_PI/ 2.,TargZenitAng = M_PI/ 2., valH = 0., valV = 20.;
	  TInitTargData InitData( Bearing,  TargCourse,  TargZenitAng , valV
	 ,TargDist, valH, valTBegin);

	 const int quantParts = 1;
	 TPartData arrPartData[1];
	 arrPartData[0].miTypePart = 1;
	 arrPartData[0].mTimePart = 1000.;
	 arrPartData[0].mSigW = 10.;

	 TShipTarg ShipTarg(InitData ,quantParts, arrPartData, NULL);
	 /// цель создана


	 long double valTet0Temp = 0.;
	 const int  iControlType = 2;
	 int iArrParams[5] = {0};
	 iArrParams[0] = 3;
	 long double DblArrParams[20] ={0.};
	 DblArrParams[1] = 13.;
	 DblArrParams[2] = 1000.;
	 DblArrParams[3] = 1000.;
	 TCntlFuncPar CntlFuncPar ( iControlType, iArrParams , DblArrParams);

	 TURPointXY * pPoints =  new TURPointXY[ quantSteps];
	 int lenLine = 0;
	 long double valMaxDist = -1.;


	 // создание BombTrajTemp

	  TBombTraj BombTrajTemp( valTet0Temp,valV0 ,valHMin ,0.,valStepInt ,valTBegin, Bomb	);
	  ///
	 THomingSituation HomingSituationOpt( ShipTarg , BombTrajTemp , CntlFuncPar	,0., NULL ) ;
	 THomingSituation HomingSituation1( ShipTarg , BombTrajTemp , CntlFuncPar	,0., NULL ) ;
	 THomingSituation HomingSituationReturn ( ShipTarg , BombTrajTemp , CntlFuncPar	,0., NULL ) ;
	 for (int i = 0; i < quantSteps ; i++)
	 {
	  // создание TBombTraj
	  long double valY0 =  valHMin + ((long double)i) * valHStep;
	  TBombTraj BombTraj( valTet0Temp,valV0 ,valY0 ,0.,valStepInt ,valTBegin, Bomb	);
	  ///

	   THomingSituation  HomingSituation( ShipTarg , BombTraj , CntlFuncPar	,0., NULL ) ;
		long double  valDHoriz = -1.;


	THomingSituation  HomingSituation0 = HomingSituation;
	HomingSituation1 = HomingSituation;


	const long double valShagT = 0.1;
	HomingSituationReturn = HomingSituation ;

	long double valMaxDistTakeoverCur = -1. ;
		bool brez = HomingSituation0.fnc_FindMaxDist_For_Takeover_Perebor( QuantOverSwitchPnts, valShagT
	  , HomingSituationReturn, valMaxDistTakeoverCur) ;

		if(brez)
		{
		  pPoints[lenLine].X = valY0;
		  pPoints[lenLine].Y =  valMaxDistTakeoverCur ;
		  lenLine++;
		  if (valMaxDistTakeoverCur > valMaxDistTakeover )
		  {
			valMaxDistTakeover = valMaxDistTakeoverCur ;
			HomingSituationOpt =  HomingSituationReturn;

			HomingSituation1.mCntlFuncPar =  HomingSituationReturn.mCntlFuncPar;
			HomingSituation1.mBombTraj.mTet0= HomingSituationReturn.mBombTraj.mTet0 ;
			HomingSituation1.mBombTraj.fncFillNachalnieUsloviaVS() ;
			HomingSituation1.mShipTarg = HomingSituationReturn.mShipTarg;
			breturn = true;
		// HomingSituation1.fncMoveClass_TO_FixedTime_AND_ShowGraphs( mpwchMaxDistZahvOutFile, HomingSituationReturn.mBombTraj.mTCur);
	  //	HomingSituation1.fncMoveClass_TO_FixedTime_AND_ShowGraphs( mpwchMaxDistZahvOutFile, HomingSituationReturn.mCntlFuncPar.mDblArrParams[HomingSituationReturn.mCntlFuncPar.miArrParams[0] -1] );HomingSituationOpt =  HomingSituationReturn; /// ИСПРАВИТЬ!!!!!
		  }

		}
	 }
	if (breturn)
	{


	TURPolyLine  PolyLine( pPoints,lenLine) ;
	wchar_t wchFileName[300]  = {0};
	wcscpy(wchFileName,pwchMaxDistZahvOutFile );
	wcscat(wchFileName, L"\\MaxDist_fromHeight.shp");
	PolyLine.WriteSetSHPFiles(wchFileName,&PolyLine, 1) ;
	delete []pPoints;


			 HomingSituation1.findTakeoverDist_AND_CreatePictures(pwchMaxDistZahvOutFile , valMaxDistTakeover) ;
			 valHeight =  HomingSituation1.mBombTraj.mAltit;
			 iTypeOfTraj  = HomingSituation1.mCntlFuncPar.miControlType;
			 valTetta0 =  HomingSituation1.mBombTraj.mTet0;
			 valTPerecl = HomingSituation1.mCntlFuncPar.mDblArrParams[1];
			 valTZahv =  HomingSituation1.mTCur;
			 valVXahv =  HomingSituation1.mBombTraj.marrStrSK_VS[3];
			 valYZahv =  HomingSituation1.mBombTraj.marrStrSK_VS[1];
			 valXZahv =  HomingSituation1.mBombTraj.marrStrSK_VS[0];
			 valTettaZahv = HomingSituation1.mBombTraj.marrStrSK_VS[4];
			 return true;
	}

  return false;
}
 // интгрирование уравнений до момента  падения
// OUTPUT:
// valDHoriz - горизонтальтеая дальн точки падения
void THomingSituation::fncMoveBomb_TO_ZeroAlt_AND_ShowGraphs( wchar_t *wcharrPath, long double &valDHoriz )
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
	wcscpy( &pwcharrFileNames[ 6 * 30], L"AlfAt");
	wcscpy( &pwcharrFileNames[ 7 * 30], L"M");
	wcscpy( &pwcharrFileNames[ 8 * 30], L"Cy");
	wcscpy( &pwcharrFileNames[ 9 * 30], L"CyqSmUmaxMG");

	pscaleY = new double  [QUANT_COLS] ;
	pscaleY[1] = 0.1;
	pscaleY[2] = 0.1;

	pscaleY[3] = 1000.;
	pscaleY[4] = 1.;
	pscaleY[5] = 1000.;
	pscaleY[6] = 1000.;
	pscaleY[7] = 1000.;
	pscaleY[8] = 1000.;
	pscaleY[9] = 1000.;
}

  int iCirc = 1000. / mBombTraj.mStepInt ;
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
	  long double valU = fncCalcU( mTCur)  ;
	  mBombTraj.fncEilerStep(valU, mBombTraj.mStepInt);
	  mTCur = mBombTraj.mTCur;





	  if (mBombTraj.marrStrSK_VS [1] < 0.) break ;


	  if  ( iNupPointsOut < QUANT_POINTS_MAX )
	  {

	  valTOut = (double)mTCur ;

	   if (wcharrPath)
	   {
	   double *p = &parrBuff[iNupPointsOut *  QUANT_COLS] ;
	   p[0] =  (double)mTCur ;
	   p[1] = (double)mBombTraj.marrStrSK_VS [0];
	   p[2] = (double)mBombTraj.marrStrSK_VS [1];
	   p[3] = (double)mBombTraj.marrStrSK_VS [2];
	   p[4] = (double)mBombTraj.marrStrSK_VS [3];
	   p[5] = (double)mBombTraj.marrStrSK_VS [4];
	   p[6] = (double)fncCalcU( mTCur);
	   long double valTay,  valDerivTay =0., valMach = -1., arrGradMach[5] = {0.};

	   fncCalcNormTemperature(p[2], valTay, valDerivTay) ;
	   mBombTraj.fncCalcMach_and_GradMach( valTay,  valDerivTay
	 ,valMach, arrGradMach);
	   p[7] = (double)valMach;
	   long double  val_q = 0., arrGrad_q [5] = {0.};
	   mBombTraj.fncCalc_q_and_Grad_q(val_q, arrGrad_q) ;
	   p[8] = mBombTraj.mBomb.fncCy(valMach);
	   p[9] = p[8] * val_q * mBombTraj.mBomb.mMidSq * mBombTraj.mBomb.mMaxAt / mBombTraj.mBomb.mMass / G_ZEMLI ;

		}

	  iNupPointsOut++;

	  }



  }
   valDHoriz = mBombTraj.marrStrSK_VS [0] ;
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

		wchar_t wchFileName4[300] = {0} ;
		wcscpy(wchFileName4, wcharrPath );
		wcscat(wchFileName4, L"\\Axes.shp");
		TYrWriteShapeFile::CreateShpAxes(wchFileName4,-80000.,80000,-80000.,80000) ;

		delete parrBuff ;
		delete pwcharrFileNames ;
		delete pscaleY ;
   }

}
void THomingSituation::fncMoveBomb_TO_ZeroAlt_AND_ShowPictures( wchar_t *wcharrPath,  long double &valDHoriz )
{
	drawJetPictures( wcharrPath);
	THomingSituation HomingSituationCopy0 = *this ;

	wchar_t pwcharrPath[900] = {0};
	for (int i =0; i < 3; i++)
	{
		wcscpy(&pwcharrPath[i *300], wcharrPath) ;
		wcscat(&pwcharrPath[i *300], L"\\Part");
		wchar_t string[5] = {0};
		_itow(i, string, 10);
		wcscat(&pwcharrPath[i *300], string);

	}
	HomingSituationCopy0.fncMoveBomb_TO_ZeroAlt_AND_ShowGraphs(wcharrPath, valDHoriz );
	bool bend = false;
	THomingSituation HomingSituationCur = *this ;
	int  iControlType = mCntlFuncPar.miControlType;
	int  iNotControlType = (iControlType == 1)?2:1;
	int  iControlTypeCur = ( mCntlFuncPar.miArrParams[0]  % 2 == 0 ? iControlType :iNotControlType);

	for (int i = 0; i < (mCntlFuncPar.miArrParams[0] - 1); i++)
	{

	 long double valFixedTime = 0.;
	 long double valFixedEndTime = mCntlFuncPar.mDblArrParams[i + 1];
	 if (mCntlFuncPar.mDblArrParams[i + 1] < HomingSituationCopy0.mTCur)
	 {
	   valFixedTime = ( mCntlFuncPar.mDblArrParams[i + 1] + mCntlFuncPar.mDblArrParams[i ])/ 2.;
	 }
	 else
	 {
	  valFixedTime = (HomingSituationCopy0.mTCur + mCntlFuncPar.mDblArrParams[i ])/ 2.;
	  valFixedEndTime = HomingSituationCopy0.mTCur;
	  bend = true;
	 }
	  _wmkdir( &pwcharrPath[i *300]);

	 THomingSituation HomingSituationCopy1 = *this ;// для рисования картинки УАБ и корабля посередине участка
	 THomingSituation HomingSituationCopy2 = *this ;
	 HomingSituationCopy1.fncMoveClass_TO_FixedTime_AND_ShowGraphs( NULL,  valFixedTime );
	 const double valRastigenie = 3.;

	 wchar_t pwcharrPathPositMiddle[300] = {0},pwcharrPathPositEnd[300] = {0};
	 wcscpy(pwcharrPathPositMiddle ,&pwcharrPath[i *300]) ;
		wcscat(pwcharrPathPositMiddle, L"\\MiddlePoint");
		_wmkdir( pwcharrPathPositMiddle);

	 HomingSituationCopy1.drawPictures( pwcharrPathPositMiddle,valRastigenie);

	  wcscpy(pwcharrPathPositEnd ,&pwcharrPath[i *300]) ;
		wcscat(pwcharrPathPositEnd, L"\\EndPoint");
		_wmkdir( pwcharrPathPositEnd);
	 HomingSituationCur.fncMoveClass_TO_FixedTime_AND_ShowGraphs( &pwcharrPath[i *300],  valFixedEndTime );
	 HomingSituationCur.drawPictures( pwcharrPathPositEnd ,valRastigenie);


	 if (bend ) break;

	 TBombTraj BombTraj  (
	HomingSituationCur.mBombTraj.marrStrSK_VS [4]  // нач угол бросания
	,HomingSituationCur.mBombTraj.marrStrSK_VS [3]  // нач скоргсть
	,HomingSituationCur.mBombTraj.marrStrSK_VS [1]  // нач высота
	 ,HomingSituationCur.mBombTraj.marrStrSK_VS [0]// нач. положение по оси X
	,HomingSituationCur.mBombTraj.mStepInt  // шагт интьерирования по времени диф уравнений движения
	 ,HomingSituationCur.mBombTraj.mTCur
	,mBombTraj.mBomb
	,HomingSituationCur.mBombTraj.marrCfW
	)  ;
	int iArrParams[NUM_INT_PARAMS] ={0};
	iArrParams[0] = 2;
	long double DblArrParams [NUM_LD_PARAMS] = {0.};
	DblArrParams[0] = valFixedEndTime;
	DblArrParams[1] = 1000.;
	iControlTypeCur = (iControlTypeCur == iControlType ?iNotControlType  : iControlType);
	//int icur = iControlType;
   //	iControlType =  iNotControlType;
   //	iNotControlType = icur;
   //	iControlTypeCur = iControlType;

	 TCntlFuncPar CntlFuncPar( iControlTypeCur, iArrParams , DblArrParams);
	 THomingSituation HomingSituationCur1(HomingSituationCur.mShipTarg
		,  BombTraj //
		,  CntlFuncPar
		,mGlonass
		,HomingSituationCur.mTCur
		,NULL // путь у папке с отчетом
		);
	 HomingSituationCur =  HomingSituationCur1;


	}
   *this = HomingSituationCopy0;

}
 // интгрирование уравнений до момента  падения
// OUTPUT:
// valDHoriz - горизонтальтеая дальн точки падения
void THomingSituation::fncMoveClas_TO_BlastPoint_AND_ShowGraphs( wchar_t *wcharrPath, long double &valDHoriz )
{
	if (mCntlFuncPar.miControlType < 10  )
	{
	 ShowMessage(L"Error fncMoveClas_TO_BlastPoint_AND_ShowGraphs: miControlType < 10");
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

  int iCirc = 1000. / mBombTraj.mStepInt ;
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
	  long double valU = fncCalcU( mTCur)  ;
	  mBombTraj.fncEilerStep(valU, mBombTraj.mStepInt);
	  mTCur = mBombTraj.mTCur;
	  mShipTarg.recalcTrajPoint(mTCur);

	  // погтсроение графиов для проверки и отладки
   /*	  TURPointXY PointsLineVelo[2]; // линия мншгновенной сорости УАБ
	  PointsLineVelo[0] = TURPointXY(mBombTraj.marrStrSK_VS [0], mBombTraj.marrStrSK_VS [1]);
	  PointsLineVelo[1] = TURPointXY(mBombTraj.marrStrSK_VS [0] - mBombTraj.marrStrSK_VS [1]/ tan(mBombTraj.marrStrSK_VS [4]), 0.);
	  TURPolyLine lineVelo(PointsLineVelo,2) ;
	  wchar_t wcharrFileName[300]={0};
   wcscpy(wcharrFileName, wcharrPath);
   wcscat(wcharrFileName, L"\\LIneVelo.shp");
   TURPolyLine::WriteSetSHPFiles(wcharrFileName,&lineVelo, 1) ;

   TURPointXY PointTarg(mBombTraj.mBomb.mHomingHead.marrMesureTarg[0],mBombTraj.mBomb.mHomingHead.marrMesureTarg[1])   ;
   wcscpy(wcharrFileName, wcharrPath);
   wcscat(wcharrFileName, L"\\Targ.shp");
   TURPointXY::WriteSetSHPFiles(wcharrFileName,&PointTarg, 1) ;
	*/
	  ///]




	  long double dx = (mBombTraj.marrStrSK_VS[0] - mShipTarg.mTraject.marrVectSostGSK[0]);
	  long double dy = (mBombTraj.marrStrSK_VS[1] - mShipTarg.mTraject.marrVectSostGSK[3]);
	  if (( sqrtl(dx*dx + dy * dy ) < mBombTraj.mBomb.mEbomb.mR) || (mBombTraj.marrStrSK_VS[1] < 0.) )
	  {
	   break ;
	  }
	 // if (((double)mTCur > (valTOut + DEL_T - 1E-15 ))&& ( iNupPointsOut < QUANT_POINTS_MAX ))
	  if  ( iNupPointsOut < QUANT_POINTS_MAX )
	  {

	  valTOut = (double)mTCur ;

	   if (wcharrPath)
	   {
	   double *p = &parrBuff[iNupPointsOut *  Quant_Cols] ;
	   p[0] =  (double)mTCur ;
	   p[1] = (double)mBombTraj.marrStrSK_VS [0];
	   p[2] = (double)mBombTraj.marrStrSK_VS [1];
	   p[3] = (double)mBombTraj.marrStrSK_VS [2];
	   p[4] = (double)mBombTraj.marrStrSK_VS [3];
	   p[5] = (double)mBombTraj.marrStrSK_VS [4];
	   p[6] = (double)fncCalcU( mTCur);
	   long double valTay,  valDerivTay =0., valMach = -1., arrGradMach[5] = {0.};

	   fncCalcNormTemperature(p[2], valTay, valDerivTay) ;
	   mBombTraj.fncCalcMach_and_GradMach( valTay,  valDerivTay
	 ,valMach, arrGradMach);
	   p[7] = (double)valMach;
	   long double  val_q = 0., arrGrad_q [5] = {0.};
	   mBombTraj.fncCalc_q_and_Grad_q(val_q, arrGrad_q) ;
	   p[8] = mBombTraj.mBomb.fncCy(valMach);
	   p[9] = p[8] * val_q * mBombTraj.mBomb.mMidSq * mBombTraj.mBomb.mMaxAt / mBombTraj.mBomb.mMass / G_ZEMLI ;
	   p[10] = (double)fncAngVisir();
	   p[11] = fabs (p[10] - p[5] + p[6]);

		}

	  iNupPointsOut++;

	  }



  }
   valDHoriz = mBombTraj.marrStrSK_VS [0] ;
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
		TYrWriteShapeFile::CreateShpAxes(wchFileName4,-80000.,80000,-80000.,80000) ; ;

		delete parrBuff ;
		delete pwcharrFileNames ;
		delete pscaleY ;
   }

}






// интгрирование уравнений до момента  падения
// wcharrPath - папка уда выводить шейп файлы
// pwchShapeUAB - папка отуда брать шейп файлы с изображением УАБ и молнии
//
// OUTPUT:
// valProb - вероятность поражения
bool THomingSituation::fncMoveHomeHead_TO_BlastPoint_AND_ShowGraphs_StatMeth( wchar_t *wcharrPath
   , TBomb BombModel, long double &valProb,  long double &valSig_AngVis )
{
	bool breturn = false;
	if (mCntlFuncPar.miControlType < 10  )
	{
	 ShowMessage(L"Error fncMoveHomeHead_TO_BlastPoint_AND_ShowGraphs_StatMeth: miControlType < 10");
	 return false;
	}
	const double DEL_T = 0.1;
	double *parrBuff;
	double *pscaleY;
	wchar_t *pwcharrFileNames ;
	const int lenName =30 ;// максимальная длина имени переменной
	const int Quant_Cols = QUANT_COLS + 11;

	double arrDist[39] ={0.};
	double arrSigAlfVisir[39] ={0.};
	double arrSigR[39] ={0.};
	int lenArrSigAlf = 0;
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
	wcscpy( &pwcharrFileNames[ 9 * 30], L"EstAngVis");
	wcscpy( &pwcharrFileNames[ 10 * 30], L"AngVisir");
	wcscpy( &pwcharrFileNames[ 11 * 30], L"DeltaAng");
	wcscpy( &pwcharrFileNames[ 12 * 30], L"ErrorAngVis");
	wcscpy( &pwcharrFileNames[ 13 * 30], L"SigmaDeltaAngVis");
	wcscpy( &pwcharrFileNames[ 14 * 30], L"FiVizZv");
	wcscpy( &pwcharrFileNames[ 15 * 30], L"EstTargX");
	wcscpy( &pwcharrFileNames[ 16 * 30], L"EstTargV");
	wcscpy( &pwcharrFileNames[ 17 * 30], L"TargXZv");
	wcscpy( &pwcharrFileNames[ 18 * 30], L"R");
	wcscpy( &pwcharrFileNames[ 19 * 30], L"RZv");
	wcscpy( &pwcharrFileNames[ 20 * 30], L"EstR");

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
	pscaleY[12] = 1000.;
	pscaleY[13] = 1000.;
	pscaleY[14] = 1000.;
	pscaleY[15] = 1.;
	pscaleY[16] = 1.;
	pscaleY[17] = 1.;
	pscaleY[18] = 1.;
	pscaleY[19] = 1.;
	pscaleY[20] = 1.;
}

  int iCirc = 1000. / mBombTraj.mStepInt ;
 // long double valTTemp = mTCur + ((long double)iCirc) * mStepInt ;
  int i = 0 ;
  int iNupPointsOut = 0;
  double valTOut = -DEL_T ;
  long double valRZv = 0., valFiVisZv  = 0., valTZv  = 0.
	 , valSig_BMO_R  = 0., valSig_MMO_R  = 0., valSig_BMO_Fi  = 0.,valSig_MMO_Fi  = 0.;
  long double valU = 0.;
  for ( i = 0; i < iCirc; i++)
  {
	   ImitateZamer(valRZv, valFiVisZv, valTZv
	   , valSig_BMO_R,  valSig_MMO_R, valSig_BMO_Fi,valSig_MMO_Fi);
	   long double arrBombVSZv[6] ;
	  ImitateZamerGlonass(arrBombVSZv);

	  // valSig_MMO_R = 0. ;

	   long double valEstFiVisPrev = mFilt .marrEstAngVis[0] ;
	   long double sigmAngVis_w = 40. * 350./ valRZv/valRZv /3. ;
	   long double sigmR_w = 1800./  mFilt.marrEstR[0] ;   // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	   long double sigmX_w = 0.01;
	   long double temp_mmo1 = sqrtl(valSig_MMO_R * cosl( valFiVisZv) * valSig_MMO_R * cosl( valFiVisZv)
			+ valRZv *  sinl( valFiVisZv) * valSig_MMO_Fi * valRZv *  sinl( valFiVisZv) * valSig_MMO_Fi);
	   long double sigmX_mmo = sqrtl( mGlonass.mSKOPos * mGlonass.mSKOPos + temp_mmo1 * temp_mmo1);
	   long double sigmX_bmo = sqrtl(valSig_BMO_R * cosl( valFiVisZv) * valSig_BMO_R * cosl( valFiVisZv)
			+ valRZv *  sinl( valFiVisZv) * valSig_BMO_Fi * valRZv *  sinl( valFiVisZv) * valSig_BMO_Fi);
	   sigmX_bmo = sqrtl(sigmX_bmo * sigmX_bmo +  sigmX_mmo * sigmX_mmo);
	   sigmX_mmo = 0.;// sigmX_bmo;
	  // long double valXZv = arrBombVSZv[0] + valRZv * cos (valFiVisZv);
	   long double valXZv = mShipTarg.mTraject.marrVectSostGSK[0] + TShipTraj::getGauss(0., sigmX_bmo );


	   if (i > 0)
	   {
	   mFilt.OneStepGolubev( valFiVisZv, valRZv, valXZv, valTZv
		 , sigmAngVis_w ,valSig_BMO_Fi , valSig_MMO_Fi
		 ,  sigmR_w, valSig_BMO_R ,valSig_MMO_R
		 ,sigmX_w, sigmX_bmo,sigmX_mmo);
	   }

		long double valEstFiVisTochka = 0.;
		if (i ==0)
		{
		  valEstFiVisTochka = mFilt .marrEstAngVis[1] ;
		}
		else
		{
		 valEstFiVisTochka = ( mFilt .marrEstAngVis[0] - valEstFiVisPrev )/ mBombTraj.mStepInt;

		}

	  //	valEstFiVisTochka = mFilt .marrEstAngVis[1] ;


	  TBombTraj BombTrajModel (atanl(arrBombVSZv[3]/arrBombVSZv[2])
	,sqrtl(arrBombVSZv[3] * arrBombVSZv[3] + arrBombVSZv[2] * arrBombVSZv[2])
	,arrBombVSZv[1] ,arrBombVSZv[0]	,mBombTraj.mStepInt ,mTCur ,BombModel );

	TShipTarg ShipTargEst =mShipTarg;
	if (ShipTargEst.mparrBuff != NULL)
	{
	  free(ShipTargEst.mparrBuff );
	  ShipTargEst.mparrBuff = NULL;
	}
	if (ShipTargEst.mpwcharrFoldReport != NULL)
	{
	  free(ShipTargEst.mpwcharrFoldReport );
	  ShipTargEst.mpwcharrFoldReport = NULL;
	}
	ShipTargEst.mLenMemoryAlloc = 0;
	ShipTargEst.mTraject.marrVectSostGSK[0] =  mFilt.marrEstTargX[0];
	ShipTargEst.mTraject.marrVectSostGSK[3] =  mFilt.marrEstTargX[1];

	THomingSituation HomingSituationModel(ShipTargEst, BombTrajModel,mCntlFuncPar ,mTCur, NULL);

	  valU = HomingSituationModel.fncCalcU( mTCur)  ;


	  mBombTraj.fncEilerStep(valU, mBombTraj.mStepInt);
	  mTCur = mBombTraj.mTCur;
	  mShipTarg.recalcTrajPoint(mTCur);

	  long double dx = (mBombTraj.marrStrSK_VS[0] - mShipTarg.mTraject.marrVectSostGSK[0]);
	  long double dy = (mBombTraj.marrStrSK_VS[1] - mShipTarg.mTraject.marrVectSostGSK[3]);
	  if (mBombTraj.marrStrSK_VS[1] < 0.)
	  {
		breturn = false;
		break ;
	  }
	 // if ( sqrtl(dx*dx + dy * dy ) < mBombTraj.mBomb.mEbomb.mR)
	//  {
	//   breturn = true;
	//   break ;
	 // }
	  long double tempy = mFilt.marrEstR[0] * cosl(mFilt.marrEstAngVis[0]) + mShipTarg.mTargData.mLDanger/ 2.;
	  long double distlong = sqrtl(arrBombVSZv[1] * arrBombVSZv[1] + tempy * tempy );
	  if ( distlong + 2. * sqrtl(mFilt.marrK_R[0]) < mBombTraj.mBomb.mEbomb.mR)
	  {
	   breturn = true;
	   break ;
	  }

	 // if (((double)mTCur > (valTOut + DEL_T - 1E-15 ))&& ( iNupPointsOut < QUANT_POINTS_MAX ))
	  if  ( iNupPointsOut < QUANT_POINTS_MAX )
	  {

	  valTOut = (double)mTCur ;

	   if (wcharrPath)
	   {
	   double *p = &parrBuff[iNupPointsOut *  Quant_Cols] ;
	   p[0] =  (double)mTCur ;
	   p[1] = (double)mBombTraj.marrStrSK_VS [0];
	   p[2] = (double)mBombTraj.marrStrSK_VS [1];
	   p[3] = (double)mBombTraj.marrStrSK_VS [2];
	   p[4] = (double)mBombTraj.marrStrSK_VS [3];
	   p[5] = (double)mBombTraj.marrStrSK_VS [4];
	   p[6] = valU;//(double)fncCalcU( mTCur);
	   long double valTay,  valDerivTay =0., valMach = -1., arrGradMach[5] = {0.};

	   fncCalcNormTemperature(p[2], valTay, valDerivTay) ;
	   mBombTraj.fncCalcMach_and_GradMach( valTay,  valDerivTay
	 ,valMach, arrGradMach);
	   p[7] = (double)valMach;
	   long double  val_q = 0., arrGrad_q [5] = {0.};
	   mBombTraj.fncCalc_q_and_Grad_q(val_q, arrGrad_q) ;
	   p[8] = mBombTraj.mBomb.fncCy(valMach);
	   p[9] = (double)mFilt.marrEstAngVis[0] ;
	   //p[8] * val_q * mBombTraj.mBomb.mMidSq * mBombTraj.mBomb.mMaxAt / mBombTraj.mBomb.mMass / G_ZEMLI ;
	   p[10] = (double)fncAngVisir();
	   p[11] = fabs (p[10] - p[5] + p[6]);
	   p[12] = p[9] - p[10];
	   p[13] = sqrt( mFilt.marrK_AngVis[0]) ;
	   p[14] = (double) valFiVisZv ;
	   p[15] = (double)(mFilt.marrEstTargX[0]);
	   p[16] = (double)(mFilt.marrEstTargX[1]);
	   p[17] =  (double)valXZv;
	   p[18] = (double)( fncDist());
	   p[19] = (double)valRZv;
	   p[20] =(double) (mFilt.marrEstR[0]);
	   int iiii = 0;
		}

		if(fncDist() < (1000. -(double (lenArrSigAlf))* 25.))
		{
			arrDist[lenArrSigAlf] =fncDist() ;
			arrSigAlfVisir[lenArrSigAlf] = calcSigAngVisExtrapol();
			arrSigR[lenArrSigAlf] = sqrt(mFilt.marrK_R [0]);
		  //	arrSigR[lenArrSigAlf] = calcSigRExtrapol();
			lenArrSigAlf ++;
		}
	  iNupPointsOut++;
	  }
  }

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
		TYrWriteShapeFile::CreateShpAxes(wchFileName4,-80000.,80000,-80000.,80000) ; ;
		delete parrBuff ;
		delete pwcharrFileNames ;
		delete pscaleY ;
   }
 if(wcharrPath)
 {
   wchar_t pwcharrPathDefeat[300] = {0};
   wcscpy(pwcharrPathDefeat, wcharrPath );
   wcscat(pwcharrPathDefeat, L"\\DefeatSector");
   _wmkdir (pwcharrPathDefeat);
   drawEMB_Blast_Pictures( pwcharrPathDefeat) ;


   wchar_t pwchShapeFileOutVess[300] ={0};
	wcscpy(pwchShapeFileOutVess, pwcharrPathDefeat);
	wcscat(pwchShapeFileOutVess, L"\\VesselBlowOut.shp");
	const double valRastigenie0 =1.;
   drawVessel( pwchShapeFileOutVess, mShipTarg.mTraject.marrVectSostGSK[0], mShipTarg.mTraject.marrVectSostGSK[3],  valRastigenie0 ) ;
   }
  if (breturn)
  {

	// проверка. того, что уязвимая часть цели находится от УАБ на расстоянии меньшем расстоянию поражения

	 long double arrVS_Ship_GSK[2] = {0.};
	 arrVS_Ship_GSK[0 ] =  mShipTarg.mTraject.marrVectSostGSK[0] + mShipTarg.mTargData.mLDanger/ 2.;
	 arrVS_Ship_GSK[1 ] =  mShipTarg.mTraject.marrVectSostGSK[2] ;
	 long double arrDelS[2] = {0.};
	 MtrxMinusMatrx( arrVS_Ship_GSK, mBombTraj.marrStrSK_VS,2, 1, arrDelS);
	 long double dist =NormVect2(arrDelS);
	 if ( dist > mBombTraj.mBomb.mEbomb.mR)
	 {
	   valProb = 0.;
	   return true;
	 }
	 ///
   /*	long double h = mBombTraj.mBomb.mTStabSyst;
	// Экстраполяция
	long double arrKExtr = 0.;

	long double sigm_w = 40. * 350./ mBombTraj.mBomb.mEbomb.mR/ mBombTraj.mBomb.mEbomb.mR /3. ;
	long double sq = sqrt(mFilt.marrK_AngVis [3] - mFilt.marrK_AngVis [1] * mFilt.marrK_AngVis [1] /mFilt.marrK_AngVis [0]) ;
	arrKExtr =  mFilt.marrK_AngVis [0] + 2. * h * mFilt.marrK_AngVis[1] + h * h * mFilt.marrK_AngVis[3]
	 +  sigm_w * sigm_w * h* h * h * h/4. + sigm_w * h * h * h * sq;




	   valSig_AngVis = sqrt(arrKExtr  +2.* mGlonass.mSKOVeloc * mGlonass.mSKOVeloc
			 /mBombTraj.marrStrSK_VS[3]/ mBombTraj.marrStrSK_VS[3]); */

	  valSig_AngVis =  calcSigAngVisExtrapol()  ;
	  const double valPsi1 = calcPsi1();
	  const double valPsi2 = calcPsi2();
	  const double a =  (  valPsi1 - mBombTraj.mBomb.mEbomb.mFi /2.) / valSig_AngVis;
	  const double b =  ( -valPsi2 + mBombTraj.mBomb.mEbomb.mFi /2.) / valSig_AngVis;
	  valProb = fncNormRaspr(  a,  b) ;
  }
  else
  {
   valProb = -1.;
   valSig_AngVis = -1.;
  }

   if(wcharrPath)
 {
   wchar_t pwcharrPath00[300] = {0};
   wcscpy(pwcharrPath00, wcharrPath );
   wcscat(pwcharrPath00, L"\\Graph_P_from_RDefeat_FiDefeat");
   _wmkdir(pwcharrPath00) ;
  createGraphs_P_from_RDefeat_FiDefeat(arrDist, arrSigAlfVisir,  arrSigR,lenArrSigAlf , pwcharrPath00);
  return breturn ;
}
}


// ссоздание графиков зависимости внроятности поражения от дальности сектора поражения и угла сектора диаграммы
void THomingSituation::createGraphs_P_from_RDefeat_FiDefeat(double *parrDist,
			double *parrSigAlfVisir,double *parrSigR, const int lenArrSigAlf , wchar_t *wcharrPath00)
{
  TURPointXY *pPonts0 = new TURPointXY[lenArrSigAlf] ;
  TURPointXY *pPonts1 = new TURPointXY[lenArrSigAlf] ;
  for (int i =0; i < lenArrSigAlf; i++)
  {
	 double FiDefeat = asin(mBombTraj.mBomb.mEbomb.mR/ parrDist[i] * sin(mBombTraj.mBomb.mEbomb.mFi));
	 TBlastSituation BlSitCur(M_PI/ 4., parrDist[i], FiDefeat
	 ,(double)mShipTarg.mTargData.mLDanger ,(double)parrSigR[i] ,(double)parrSigAlfVisir[i]) ;
	 double valDTresh =0., valProbTresh = 0.;
	 BlSitCur.findOptBlastTreshold_and_createProbabGraph(valDTresh 	,valProbTresh, NULL) ;
	 pPonts0[i].X = parrDist[i];
	 pPonts0[i].Y = valProbTresh *100.;
	 pPonts1[i].X = FiDefeat * 180./ M_PI;
	 pPonts1[i].Y = valProbTresh*100.;

  }

  if (wcharrPath00)
   {
   wchar_t pwchFileName0[300] = {0};
   wcscpy(pwchFileName0, wcharrPath00 );
   wcscat(pwchFileName0, L"\\Graph_P_from_RDefeat.shp");
   wchar_t pwchFileName1[300] = {0};
   wcscpy(pwchFileName1, wcharrPath00 );
   wcscat(pwchFileName1, L"\\Graph_P_from_FiDefeat.shp");
	 TURPolyLine polyLine0( pPonts0,lenArrSigAlf)  ;
	 polyLine0.WriteSetSHPFiles(pwchFileName0, &polyLine0,1) ;
	 TURPolyLine polyLine1( pPonts1,lenArrSigAlf)  ;
	 polyLine1.WriteSetSHPFiles(pwchFileName1, &polyLine1,1) ;

   }
   delete []pPonts0 ;
   delete []pPonts1 ;
}
 //(const double AlfVis, const double RDefeat,const double FiDefeat
	// ,const double LDanger ,const double SigR ,const double SigAlf)
// вычисление эустраполированного СКЗ ошибки
// угла визирования

double THomingSituation::calcSigAngVisExtrapol()
{
	long double h = mBombTraj.mBomb.mTStabSyst;
	// Экстраполяция
	long double arrKExtr = 0.;

//long double sigm_w = 40. * 350./ mBombTraj.mBomb.mEbomb.mR/ mBombTraj.mBomb.mEbomb.mR /3. ;
	long double dist = fncDist();
	long double sigm_w = 40. * 350./ dist/ dist/3. ;
	long double sq = sqrt(mFilt.marrK_AngVis [3] - mFilt.marrK_AngVis [1] * mFilt.marrK_AngVis [1] /mFilt.marrK_AngVis [0]) ;
	arrKExtr =  mFilt.marrK_AngVis [0] + 2. * h * mFilt.marrK_AngVis[1] + h * h * mFilt.marrK_AngVis[3]
	 +  sigm_w * sigm_w * h* h * h * h/4. + sigm_w * h * h * h * sq;

	 double  valSig_AngVis = sqrt(arrKExtr  +2.* mGlonass.mSKOVeloc * mGlonass.mSKOVeloc
			 /mBombTraj.marrStrSK_VS[3]/ mBombTraj.marrStrSK_VS[3]);
  return valSig_AngVis;
}

// вычисление эустраполированного СКЗ ошибки
// по дальности

double THomingSituation::calcSigRExtrapol()
{
	long double h = mBombTraj.mBomb.mTStabSyst;
	// Экстраполяция
	long double arrKExtr = 0.;

//
	long double dist = fncDist();
	long double sigm_w = 1800./  dist ;

	long double sq = sqrt(mFilt.marrK_R [3] - mFilt.marrK_R [1] * mFilt.marrK_R [1] /mFilt.marrK_R [0]) ;
	arrKExtr =  mFilt.marrK_R [0] + 2. * h * mFilt.marrK_R[1] + h * h * mFilt.marrK_R[3]
	 +  sigm_w * sigm_w * h* h * h * h/4. + sigm_w * h * h * h * sq;


  return sqrt(arrKExtr);
}

// интгрирование уравнений до момента  падения
// wcharrPath - папка уда выводить шейп файлы
// pwchShapeUAB - папка отуда брать шейп файлы с изображением УАБ и молнии
//
// OUTPUT:
// valProb - вероятность поражения
bool THomingSituation::fncCreateGraph_Prob_from_FiDefeat( wchar_t *wcharrPath )
{
	bool breturn = false;
	if (mCntlFuncPar.miControlType < 10  )
	{
	 ShowMessage(L"Error fncMoveHomeHead_TO_BlastPoint_AND_ShowGraphs_StatMeth: miControlType < 10");
	 return false;
	}
	const double DEL_T = 0.1;
	double *parrBuff;
	double *pscaleY;
	wchar_t *pwcharrFileNames ;
	const int lenName =30 ;// максимальная длина имени переменной
	const int Quant_Cols = 3;
   if (wcharrPath)
   {
	parrBuff = new double [Quant_Cols * QUANT_POINTS_MAX] ;
	memset (parrBuff, 0, Quant_Cols * QUANT_POINTS_MAX * sizeof(double)) ;
	pwcharrFileNames = new wchar_t [ Quant_Cols * lenName] ;
	memset (pwcharrFileNames, 0, Quant_Cols * lenName* sizeof(wchar_t)) ;


	wcscpy( &pwcharrFileNames[ 0], L"FiDefeat");
	wcscpy( &pwcharrFileNames[ 1 * 30], L"RDefeat");
	wcscpy( &pwcharrFileNames[ 2 * 30], L"Prob");


	pscaleY = new double  [Quant_Cols] ;

	pscaleY[0] = 1000.;
	pscaleY[1] = 1.;
	pscaleY[2] = 100.;

}

  int iCirc = 1000. / mBombTraj.mStepInt ;
 // long double valTTemp = mTCur + ((long double)iCirc) * mStepInt ;
  int i = 0 ;
  int iNupPointsOut = 0;
  double valTOut = -DEL_T ;
  long double valRZv = 0., valFiVisZv  = 0., valTZv  = 0.
	 , valSig_BMO_R  = 0., valSig_MMO_R  = 0., valSig_BMO_Fi  = 0.,valSig_MMO_Fi  = 0.;
  long double valU = 0.;
  for ( i = 0; i < iCirc; i++)
  {
	   ImitateZamer(valRZv, valFiVisZv, valTZv
	   , valSig_BMO_R,  valSig_MMO_R, valSig_BMO_Fi,valSig_MMO_Fi);
		valFiVisZv = fncAngVisir();
		valRZv =  fncDist();
        long double arrBombVSZv[4] ;
		arrBombVSZv[0] = mBombTraj.marrStrSK_VS[0] ;
		arrBombVSZv[1] = mBombTraj.marrStrSK_VS[1] ;
		arrBombVSZv[2] = mBombTraj.marrStrSK_VS[3] * cosl(mBombTraj.marrStrSK_VS[4]) ;
		arrBombVSZv[3] = mBombTraj.marrStrSK_VS[3] * sinl(mBombTraj.marrStrSK_VS[4]);

	   long double valEstFiVisPrev = mFilt .marrEstAngVis[0] ;
	   long double sigmAngVis_w = 40. * 350./ valRZv/valRZv /3.* 3. ;
	   long double sigmR_w = 1800./  mFilt.marrEstR[0] ;   // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   long double sigmX_w = 0.4;
	   long double temp_mmo1 = sqrtl(valSig_MMO_R * cosl( valFiVisZv) * valSig_MMO_R * cosl( valFiVisZv)
			+ valRZv *  sinl( valFiVisZv) * valSig_MMO_Fi * valRZv *  sinl( valFiVisZv) * valSig_MMO_Fi);
	   long double sigmX_mmo = sqrtl( mGlonass.mSKOPos * mGlonass.mSKOPos + temp_mmo1 * temp_mmo1);
	   long double sigmX_bmo = sqrtl(valSig_BMO_R * cosl( valFiVisZv) * valSig_BMO_R * cosl( valFiVisZv)
			+ valRZv *  sinl( valFiVisZv) * valSig_BMO_Fi * valRZv *  sinl( valFiVisZv) * valSig_BMO_Fi);

	   long double valXZv = arrBombVSZv[0] + valRZv * cos (valFiVisZv);

	   sigmX_bmo = sqrtl(sigmX_bmo * sigmX_bmo + sigmX_mmo * sigmX_mmo );
	   sigmX_mmo = 0.;


	   mFilt.OneStepGolubev( valFiVisZv, valRZv, valXZv, valTZv
		 , sigmAngVis_w ,valSig_BMO_Fi , valSig_MMO_Fi
		 ,  sigmR_w, valSig_BMO_R ,valSig_MMO_R
		 ,sigmX_w, sigmX_bmo,sigmX_mmo);

		long double valEstFiVisTochka = 0.;
		if (i ==0)
		{
		  valEstFiVisTochka = mFilt .marrEstAngVis[1] ;
		}
		else
		{
		 valEstFiVisTochka = ( mFilt .marrEstAngVis[0] - valEstFiVisPrev )/ mBombTraj.mStepInt;

		}

	  //	valEstFiVisTochka = mFilt .marrEstAngVis[1] ;






	  valU = fncCalcU1(valEstFiVisTochka, mTCur)  ;


	  mBombTraj.fncEilerStep(valU, mBombTraj.mStepInt);
	  mTCur = mBombTraj.mTCur;
	  mShipTarg.recalcTrajPoint(mTCur);






	  if (mBombTraj.marrStrSK_VS[1] < 0.)	  {

		break ;
	  }
	 long double arrVS_Ship_GSK[2] = {0.};
	 arrVS_Ship_GSK[0 ] =  mShipTarg.mTraject.marrVectSostGSK[0] + mShipTarg.mTargData.mLDanger/ 2.;
	 arrVS_Ship_GSK[1 ] =  mShipTarg.mTraject.marrVectSostGSK[2] ;
	 long double arrDelS[2] = {0.};
	 MtrxMinusMatrx( arrVS_Ship_GSK, mBombTraj.marrStrSK_VS,2, 1, arrDelS);
	 long double valRDefeat =NormVect2(arrDelS);

	  if ( valRDefeat < 1500.)
	  {
		if ( valRDefeat < 200.)
		{
		 breturn = true;
		 break ;
		}
		if  ( iNupPointsOut < QUANT_POINTS_MAX )
	  {
	   double *p = &parrBuff[iNupPointsOut *  Quant_Cols] ;
		p[1] =  (double)valRDefeat ;
		double valPowerTot = mBombTraj.mBomb.mEbomb.fncCalcPower ();
		p[0] = 2. * acos(1. - valPowerTot/ 2. / M_PI/(mBombTraj.mBomb.mEbomb.mW*10000.)/ p[1]/p[1] );
	   //	p[2] = calcProb
	   // расчет вероятности
		long double h = mBombTraj.mBomb.mTStabSyst;
	// Экстраполяция
	long double arrKExtr = 0.;

	long double sigm_w = 40. * 350./ valRDefeat/ valRDefeat /3. ;
	long double sq = sqrt(mFilt.marrK_AngVis [3] - mFilt.marrK_AngVis [1] * mFilt.marrK_AngVis [1] /mFilt.marrK_AngVis [0]) ;
	arrKExtr =  mFilt.marrK_AngVis [0] + 2. * h * mFilt.marrK_AngVis[1] + h * h * mFilt.marrK_AngVis[3]
	 +  sigm_w * sigm_w * h* h * h * h/4. + sigm_w * h * h * h * sq;



	   double valFiDefeat = p[0];
	   double valSig_AngVis = sqrt(arrKExtr  +2.* mGlonass.mSKOVeloc * mGlonass.mSKOVeloc
			 /mBombTraj.marrStrSK_VS[3]/ mBombTraj.marrStrSK_VS[3]);
	  const double valPsi1 = calcPsi1();
	  const double valPsi2 = calcPsi2();
	  const double a =  (  valPsi1 - valFiDefeat /2.) / valSig_AngVis;
	  const double b =  ( -valPsi2 + valFiDefeat /2.) / valSig_AngVis;
	  double valProb = fncNormRaspr(  a,  b) ;
	  p[2] =  valProb;
		iNupPointsOut++;
	  }

	  }


  }

 if (wcharrPath)
 {
		wchar_t wcharrPath1[300] = {0} ;
		wcscpy(wcharrPath1, wcharrPath );
		wcscat(wcharrPath1, L"\\");
        TYrWriteShapeFile::WriteOneReport(wcharrPath1  // путь к папке
							  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
							  ,Quant_Cols // - к-во переменных о корорых накоплена информация в буфере
							  ,iNupPointsOut //  - к-во точек
							  ,pwcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
							  ,lenName // максимальная длина имени переменной
							  ,0  // номер переменной по оси X
							  ,2  // номер переменной по оси Y
							  ,1. //  масштаб по оси X
							  ,pscaleY[1]  // масштаб по оси Y
							   );
        TYrWriteShapeFile::WriteOneReport(wcharrPath1  // путь к папке
							  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
							  ,Quant_Cols // - к-во переменных о корорых накоплена информация в буфере
							  ,iNupPointsOut //  - к-во точек
							  ,pwcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
							  ,lenName // максимальная длина имени переменной
							  ,1  // номер переменной по оси X
							  ,2  // номер переменной по оси Y
							  ,1. //  масштаб по оси X
							  ,pscaleY[2]  // масштаб по оси Y
							   );

		TYrWriteShapeFile::WriteOneReport(wcharrPath1  // путь к папке
							  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
							  ,Quant_Cols // - к-во переменных о корорых накоплена информация в буфере
							  ,iNupPointsOut //  - к-во точек
							  ,pwcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
							  ,lenName // максимальная длина имени переменной
							  ,0  // номер переменной по оси X
							  ,1  // номер переменной по оси Y
							  ,1.//  масштаб по оси X
							  ,1.  // масштаб по оси Y
							   );

		wchar_t wchFileName4[300] = {0} ;
		wcscpy(wchFileName4, wcharrPath );
		wcscat(wchFileName4, L"\\Axes.shp");
		TYrWriteShapeFile::CreateShpAxes(wchFileName4,-80000.,80000,-80000.,80000) ; ;
		delete parrBuff ;
		delete pwcharrFileNames ;
		delete pscaleY ;
   }




  return breturn ;
}
// Имитация замера ГЛОНАСС
void THomingSituation::ImitateZamerGlonass(long double  *arrBombVSZv)
{
  arrBombVSZv[0] = mBombTraj.marrStrSK_VS[0] + TShipTraj::getGauss(0., mGlonass.mSKOPos );
  arrBombVSZv[1] = mBombTraj.marrStrSK_VS[1] + TShipTraj::getGauss(0., mGlonass.mSKOPos );
  arrBombVSZv[2] = mBombTraj.marrStrSK_VS[3] * cosl(mBombTraj.marrStrSK_VS[4]) + TShipTraj::getGauss(0., mGlonass.mSKOVeloc );
  arrBombVSZv[3] = mBombTraj.marrStrSK_VS[3] * sinl(mBombTraj.marrStrSK_VS[4]) + TShipTraj::getGauss(0., mGlonass.mSKOVeloc );


}
// вычисление угла атаки (  управления) для метода стат испытаний на последнем участке траектории

//  Если тип управления 10 = miControlType - самонаведение методом пропорциональной навигации
//    mDblArrParams[0] - коэфиициент навигации


 long double THomingSituation::fncCalcU__( long double TCur)
 {

	if( mCntlFuncPar.miControlType != 10 ) return 0.;

	// вычисление нормальной виртуальной температуры  и ее градиента
	long double valTay = 0.,valDerivTay = 0. ;
	fncCalcNormTemperature(mBombTraj.marrStrSK_VS[1],valTay, valDerivTay)  ;
	///

	// вычисление числа Маха и вектора его частных производных
	long double valMach = 0.,arrGradMach[5] ={0.} ;
	mBombTraj.fncCalcMach_and_GradMach(valTay, valDerivTay ,valMach, arrGradMach);
	///
	// вычисление фунции q и вектора ее градиента
	long double val_q = 0., arrGrad_q[5] = {0.} ;
	mBombTraj.fncCalc_q_and_Grad_q( val_q, arrGrad_q);
	///

	// вычисление фунции Cy и вектора ее градиента
	long double val_Cy = 0., arrGrad_Cy[5] = {0.} ;
	mBombTraj.fncCalc_Cy_and_Grad_Cy(valMach,arrGradMach, val_Cy, arrGrad_Cy);
	if(val_Cy < 0.5)
	{
	int iii = 0;
	}
	///
	long double valFiTochka = 0., valAlfTemp = 0.;

	valFiTochka = mFilt.marrEstAngVis[1];
	valAlfTemp = -(G_ZEMLI * cosl(mBombTraj.marrStrSK_VS [4])/ mBombTraj.marrStrSK_VS [3]
	+ mCntlFuncPar.mDblArrParams[0] * valFiTochka) // mCntlFuncPar.mDblArrParams[0] - коэфф навигации
	*  mBombTraj.mBomb.mMass *  mBombTraj.marrStrSK_VS [3] /(val_Cy * val_q * mBombTraj.mBomb.mMidSq );
	if (fabs(valAlfTemp)> mBombTraj.mBomb.mMaxAt )
	{
	valAlfTemp =   mBombTraj.mBomb.mMaxAt *fncSign (valAlfTemp);
	}

	return -valAlfTemp ;

 }

 long double THomingSituation::fncCalcU1(long double valEstFiVisTochka, long double TCur)
 {

	if( mCntlFuncPar.miControlType != 10 ) return 0.;

	// вычисление нормальной виртуальной температуры  и ее градиента
	long double valTay = 0.,valDerivTay = 0. ;
	fncCalcNormTemperature(mBombTraj.marrStrSK_VS[1],valTay, valDerivTay)  ;
	///

	// вычисление числа Маха и вектора его частных производных
	long double valMach = 0.,arrGradMach[5] ={0.} ;
	mBombTraj.fncCalcMach_and_GradMach(valTay, valDerivTay ,valMach, arrGradMach);
	///
	// вычисление фунции q и вектора ее градиента
	long double val_q = 0., arrGrad_q[5] = {0.} ;
	mBombTraj.fncCalc_q_and_Grad_q( val_q, arrGrad_q);
	///

	// вычисление фунции Cy и вектора ее градиента
	long double val_Cy = 0., arrGrad_Cy[5] = {0.} ;
	mBombTraj.fncCalc_Cy_and_Grad_Cy(valMach,arrGradMach, val_Cy, arrGrad_Cy);
	if(val_Cy < 0.5)
	{
	int iii = 0;
	}
	///
	long double  valAlfTemp = 0.;


	valAlfTemp = -(G_ZEMLI * cosl(mBombTraj.marrStrSK_VS [4])/ mBombTraj.marrStrSK_VS [3]
	+ mCntlFuncPar.mDblArrParams[0] * valEstFiVisTochka) // mCntlFuncPar.mDblArrParams[0] - коэфф навигации
	*  mBombTraj.mBomb.mMass *  mBombTraj.marrStrSK_VS [3] /(val_Cy * val_q * mBombTraj.mBomb.mMidSq );
	if (fabs(valAlfTemp)> mBombTraj.mBomb.mMaxAt )
	{
	valAlfTemp =   mBombTraj.mBomb.mMaxAt *fncSign (valAlfTemp);
	}

	return valAlfTemp ;

 }





// вычисление угла атаки (  управления)
// Расчет функции управления УАБ с использованием члена класса TCntlFuncPar mCntlFuncPar
// miControlType - указывает тип управлкения
// в массивах  miArrParams и  mDblArrParams хранятся парает ры, описывающие оуправляющую фунцию
// В зависимости от типа управления miControlType форируется во внешней програме алгорит управления
// 1. Если тип управления 0 = miControlType , то движение происходит по баллистике,  alfaAt = 0
// 2. Если тип управления 1 = miControlType   - на последнем участке  планирование  alfaAt = -  mBomb.mMaxAt
//    При этом, массивы параметров хранят следующую информацию
//    miArrParams[0] = длина массива времен траектории (число переключений +2)
//    mDblArrParams - массив времен
//    mDblArrParams[0] - время начала движениея
//    mDblArrParams[i], i = 1,..., (miArrParams[0] -2) - времена переключений управления
//    mDblArrParams[ miArrParams[0] -1 ]  - конец движения
// 3. Если тип управления 2 = miControlType   - на последнем участке  пикирование  alfaAt =   mBomb.mMaxAt
//    остальное аналогично случаю 2.
// 4. Если тип управления 10 = miControlType - самонаведение методом пропорциональной навигации
//    mDblArrParams[0] - коэфиициент навигации
// 5. Если тип управления 11 = miControlType - самонаведение методом параллельного сближения
// 6. Если тип управления 12 = miControlType - метод оптимального управления по быстродействию на встречу с целью
 long double THomingSituation::fncCalcU( long double TCur)
 {


	int ind = 0;

	switch( mCntlFuncPar.miControlType )
	{
		case 0:   // движение происходит по баллистике,  alfaAt = 0
		return 0.;
		//
		case 1:      //на последнем участке  планирование  alfaAt = -  mBomb.mMaxAt
		ind =( mCntlFuncPar.miArrParams[0]  % 2 == 0 ?  1 :-1);
		for (int i = 0; i < mCntlFuncPar.miArrParams[0]; i++)
		{
		if(TCur >= mCntlFuncPar.mDblArrParams[i] )
		ind = - ind;
		else
		break;
		}
		return  ((long double)ind) * mBombTraj.mBomb.mMaxAt ;
		//
		case 2:         // на последнем участке  пикирование  alfaAt =   mBomb.mMaxAt
		ind =(mCntlFuncPar.miArrParams[0]  % 2 == 0 ?  -1 :1);
		for (int i = 0; i < mCntlFuncPar.miArrParams[0]; i++)
		{
		if(TCur >= mCntlFuncPar.mDblArrParams[i]  )
		ind = - ind;
		else
		break;
		}
		return  ((long double)ind) * mBombTraj.mBomb.mMaxAt ;
		//
		default:
		break;
	}



 // вычисление нормальной виртуальной температуры  и ее градиента
  long double valTay = 0.,valDerivTay = 0. ;
  fncCalcNormTemperature(mBombTraj.marrStrSK_VS[1],valTay, valDerivTay)  ;
  ///

  // вычисление числа Маха и вектора его частных производных
  long double valMach = 0.,arrGradMach[5] ={0.} ;
  mBombTraj.fncCalcMach_and_GradMach(valTay, valDerivTay ,valMach, arrGradMach);
  ///
 // вычисление фунции q и вектора ее градиента
  long double val_q = 0., arrGrad_q[5] = {0.} ;
  mBombTraj.fncCalc_q_and_Grad_q( val_q, arrGrad_q);
  ///

  // вычисление фунции Cy и вектора ее градиента
  long double val_Cy = 0., arrGrad_Cy[5] = {0.} ;
  mBombTraj.fncCalc_Cy_and_Grad_Cy(valMach,arrGradMach, val_Cy, arrGrad_Cy);
  if(val_Cy < 0.5)
  {
      int iii = 0;
  }
  ///
  long double valFiTochka = 0., valAlfTemp = 0.;
  switch( mCntlFuncPar.miControlType )
  {
	  case 10:
		valFiTochka = fncCalcDeltaAngVisir();
		valAlfTemp = -(G_ZEMLI * cosl(mBombTraj.marrStrSK_VS [4])/ mBombTraj.marrStrSK_VS [3]
		 + mCntlFuncPar.mDblArrParams[0] * valFiTochka) // mCntlFuncPar.mDblArrParams[0] - коэфф навигации
		*  mBombTraj.mBomb.mMass *  mBombTraj.marrStrSK_VS [3] /(val_Cy * val_q * mBombTraj.mBomb.mMidSq );
		if (fabs(valAlfTemp)> mBombTraj.mBomb.mMaxAt )
		{
		valAlfTemp =   mBombTraj.mBomb.mMaxAt *fncSign (valAlfTemp);
		}

		return valAlfTemp ;
	  case 11:

	  case 12:

	  default:
	  break;
  }
 return 0.;
 }


 // вычисление скорости изменения угла визирования цели
long double THomingSituation::fncCalcDeltaAngVisir()
{
	 long double arrVS_Ship_GSK[4] = {0.};
	 arrVS_Ship_GSK[0 ] =  mShipTarg.mTraject.marrVectSostGSK[0] ;
	 arrVS_Ship_GSK[1 ] =  mShipTarg.mTraject.marrVectSostGSK[2] ;
	 arrVS_Ship_GSK[2] =   mShipTarg.mTraject.marrVectSostGSK[3] ;
	 arrVS_Ship_GSK[3 ] =  mShipTarg.mTraject.marrVectSostGSK[5] ;


   long double arrVS_Bomb_GSK[4] = {0.};
	 arrVS_Bomb_GSK[0 ] =  mBombTraj.marrStrSK_VS[0] ;
	 arrVS_Bomb_GSK[1 ] =  mBombTraj.marrStrSK_VS[1] ;
	 arrVS_Bomb_GSK[2] =   mBombTraj.marrStrSK_VS [3] * cosl(mBombTraj.marrStrSK_VS [4]) ;
	 arrVS_Bomb_GSK[3 ] =  mBombTraj.marrStrSK_VS [3] * sinl(mBombTraj.marrStrSK_VS [4]) ;



   long double arrDelS[2] = {0.};
  MtrxMinusMatrx( arrVS_Ship_GSK, arrVS_Bomb_GSK,2, 1, arrDelS);


  long double  arrDelV [2] ={0.};
  MtrxMinusMatrx(&(arrVS_Ship_GSK[2]), &(arrVS_Bomb_GSK[2]),2, 1, arrDelV);

  long double temp = (OuterProduct_2(arrDelS , arrDelV)/NormVect2(arrDelS)/NormVect2(arrDelS));
  if (fabs(temp) > 1.)
  {
   //	ShowMessage(L"ERROR3");
   //	int iii = 0;
	temp = temp / fabsl( temp) ;
  }
 // return asinl(OuterProduct_2(arrDelS , arrDelV)/NormVect2(arrDelS)/NormVect2(arrDelS));
 return asinl(temp);


}



// OCP - Optimal Control Problem
// Решение системы нелинейных уравнений отностительно моментов  переключений
// INPUT:
// Tet0 -  угол бросания
// V0  - нач скоргсть
// Altit  - нач высота
//  StepInt  - шаг интьерирования по времени диф уравнений движения
// arrTimes - начльный массив  моментов вреени. arrTimes[0] - начльный момент движения
//         arrTimes[1] - arrTimes[QuantTimes -2]  - моменты переключений, задаются произвольно
//         arrTimes[QuantTimes -1]  - большое число, например, 1000.
// QuantTimes - длина массива моментов времени. = 2 + число переключений
// Bomb - УАБ
 bool THomingSituation::fncSolveOCP_for_MaxDist_MethNewton(TCntlFuncPar &CntlFuncPar)
{
  int iControlType =  CntlFuncPar.miControlType;
  long double *arrTimes = CntlFuncPar.mDblArrParams;
  int QuantTimes = CntlFuncPar.miArrParams[0] ;
  THomingSituation HomingSituation = *this;

  // поиск первого приближения
  long double* arrZ =  new long double[ QuantTimes] ;
  if (!HomingSituation.InitialApprox(arrZ))
  {
	  CntlFuncPar.miArrParams[0] = 2 ; //
	  CntlFuncPar.mDblArrParams[0] = 0.;
	  CntlFuncPar.mDblArrParams[1] = 1000.;
	  return false ;// точек переключения нет
  }
  TCntlFuncPar CntlFuncPar0 = CntlFuncPar;
  CntlFuncPar0.mDblArrParams[1] = arrZ[0];
  CntlFuncPar0.mDblArrParams[2] = arrZ[1];
  CntlFuncPar0.miArrParams[0] = 3;
  ///

  // цикл метода Ньютона
  bool bSuccess = false;
  const long double lEps = 0.001 ;
  long double* arrFgr = new long double[ QuantTimes];
  long double* mtrxH  = new long double[ QuantTimes * QuantTimes];
  long double* mtrxInvH  = new long double[ QuantTimes * QuantTimes];
  long double* arrZt  = new long double[ QuantTimes ];
  long double* arr_dXt  = new long double[ QuantTimes ];
  long double* arr_dXt1  = new long double[ QuantTimes ];

  for (int i = 0; i < 100; i++)
  {
	 THomingSituation  HomingSituation = *this ;
	memcpy( &(HomingSituation.mCntlFuncPar.mDblArrParams[1]), arrZ
		, (HomingSituation.mCntlFuncPar.miArrParams[0] - 1) * sizeof( long double)) ;
	HomingSituation.calc_Fgr_and_H( arrZ,  arrFgr,  mtrxH);
	InverseMtrx( mtrxH, QuantTimes,  mtrxInvH)  ;
	MtrxMultMatrx(mtrxInvH,QuantTimes, QuantTimes, arrFgr, 1, arr_dXt1) ;
	MatrxMultScalar(arr_dXt1, QuantTimes, 1, 0.01,arr_dXt);
	MtrxMinusMatrx(arrZ, arr_dXt,QuantTimes, 1, arrZt);

	memcpy(arrZ, arrZt, QuantTimes * sizeof(long double));
	if(!IsSolutionTrue(arrZ, arrZt))
	{
      CntlFuncPar = CntlFuncPar0;
	  return false ;// точек переключения нет
	}
	long double dF = NormVect(arr_dXt, QuantTimes);
	if (dF < lEps)
	{

	  bSuccess = true;
	  break ;
	}



  }
  if (bSuccess)
  {
   memcpy(&(CntlFuncPar.mDblArrParams[1]),   arrZ, (QuantTimes-1) * sizeof(long double));
  }
  else
  {
	 CntlFuncPar= CntlFuncPar0;
  }
  delete arrZ ;
  delete arrFgr ;
  delete mtrxH ;
  delete mtrxInvH;
  delete arrZt ;
  delete arr_dXt;
  delete arr_dXt1;
  return  bSuccess ;

}

// Проверка того, что решение на следующей итерации подходит
//... НАПИСАТЬ ПОДРОБНЕЕ !!!!

bool THomingSituation::IsSolutionTrue(long double* arrX, long double* arrXt)
{
	if (arrX[0] < 0.) return false;
	if (arrX[1] < 0.) return false;
	if (arrX[0] >= arrX[1]) return false;

	return true ;
}
// поиск первого приближения  моментов переключений
// моменты переклюжчений хранятся в массиве marrTimes
// время падения - mTCur
// OUTPUT:
// arrX - начальное приближение массива переменных (искомых)
// arrX [0] - arrX[QuantPerecluchen -3] - моменты пережлючений
// arrX[QuantPerecluchen -2]   - время падения T
// arrX[QuantPerecluchen -1]  - множитель Лагранжа "лямбда"
 bool THomingSituation:: InitialApprox( long double* arrZ)
 {

   int QuantTimes = mCntlFuncPar.miArrParams[0] ;
   long double *parrTimes = mCntlFuncPar.mDblArrParams;

   if (!(( QuantTimes == 2 )||(  QuantTimes == 3 )))
   {
	  ShowMessage(L"ERROR  !(( mQuantTimes == 3 )||(  mQuantTimes == 4 ))") ;
   }

   THomingSituation HomingSituationTemp = *this ;
   long double valDHoriz = 0.;
   TCntlFuncPar CntlFuncPar;
   THomingSituation HomingSituationTemp1
	   = HomingSituationTemp.fncFindOptimalControl_for_Dist_MethodPerebora_1Point(CntlFuncPar,1.,valDHoriz) ;
   arrZ[0] = HomingSituationTemp1.mCntlFuncPar.mDblArrParams[1];
   arrZ[1] = HomingSituationTemp1.mBombTraj.mTCur;

   if (HomingSituationTemp1.mCntlFuncPar.mDblArrParams[1] == HomingSituationTemp1.mCntlFuncPar.mDblArrParams[0])
   {
	return false;
   }

   THomingSituation HomingSituationTemp2 = HomingSituationTemp1;
   HomingSituationTemp1.fncMoveBomb_TO_ZeroAlt_AND_ShowGraphs(NULL, valDHoriz ) ;
   arrZ[0] = HomingSituationTemp1.mCntlFuncPar.mDblArrParams[1];
   arrZ[1] = HomingSituationTemp1.mBombTraj.mTCur;



   ///
   long double* arrA = new long double [QuantTimes -1];
   long double* arrB = new long double [QuantTimes -1];
   // ищем скорость  в омент падения
   arrA[0] =  HomingSituationTemp1.mBombTraj.marrStrSK_VS [3] * cosl(HomingSituationTemp1.mBombTraj.marrStrSK_VS [4]) ;
   arrB[0] =  HomingSituationTemp1.mBombTraj.marrStrSK_VS [3] * sinl(HomingSituationTemp1.mBombTraj.marrStrSK_VS [4]) ;
   // ищем градиенты по вреенам переключения   &arrA[1], &arrB[1]
   int quantTPerecl = -1;

   HomingSituationTemp2.fncCalc_GradX1_and_GradX2_po_dTay( HomingSituationTemp1.mBombTraj.mTCur,quantTPerecl, &arrA[1], &arrB[1]  );
   ///
   if (quantTPerecl != (QuantTimes - 2))
   {
	// .....
   }
	// нахождение лябда =  arrZ [(mQuantTimes - 1)]
	long double valNormB = NormVect(arrB, quantTPerecl + 1);
	arrZ [QuantTimes - 1] = -ScalProduct(arrA, arrB, quantTPerecl + 1)/ (valNormB * valNormB) ;

	delete arrA ;
	delete arrB ;

	return true;

 }

// вычисленпие правой части системы нелинейных уравнений относительно x
 // и матрицы частных производных  (Fgr - F греческое (круглое) )
 // ДЛЯ НАХОЖДЕНИЯ МТОМЕНТОВ ПЕРЕКЛЮЧЕНИЙ !!!!
 void THomingSituation::calc_Fgr_and_H( long double* arrZ, long double* arrFgr, long double* mtrxH)
 {
   calc_Fgr( arrZ,  arrFgr);
   int QuantTimes = mCntlFuncPar.miArrParams[0] ;
   long double *parrTimes = mCntlFuncPar.mDblArrParams;
   long double* mtrxHTransp = new long double[ QuantTimes * QuantTimes] ;
   long double  *arrTemp0 = new long double[QuantTimes ];
   long double  *arrTemp1 = new long double[QuantTimes ];
   long double  *arrZTemp = new long double[QuantTimes ];
   long double  valStep = 0.1;
   for (int i = 0; i < QuantTimes; i++)
   {
	memcpy(arrZTemp, arrZ, QuantTimes * sizeof (long double)) ;
	arrZTemp [i] +=  valStep;

	calc_Fgr( arrZTemp,  arrTemp0);
	MtrxMinusMatrx(arrTemp0, arrFgr,1, QuantTimes, arrTemp1);
	MatrxMultScalar( arrTemp1, 1, QuantTimes, 1./valStep ,&mtrxHTransp[ i * QuantTimes]);
   }
   MatrTransp(mtrxHTransp,QuantTimes , QuantTimes, mtrxH);
   delete mtrxHTransp ;
   delete arrTemp0;
   delete arrTemp1;
   delete arrZTemp ;
 }
// вычисление правой части систеы нелинейных уравнений нахождения моментов переключений
// задачи на максиум дальности
 void THomingSituation::calc_Fgr( long double* arrZ, long double* arrFgr)
 {
	THomingSituation HomingSitTemp = *this;// BombTrajTemp (mTet0,mV0,mAltit,mX0,mStepInt, miControlType, marrTimes,mQuantTimes ,mBomb );
	int QuantTimes =HomingSitTemp.mCntlFuncPar.miArrParams[0] ;
	long double *parrTimes = HomingSitTemp.mCntlFuncPar.mDblArrParams;
	memcpy(&(parrTimes[1] ),arrZ, (QuantTimes- 1) * sizeof( long double)) ;
	long double  *arrA = new long double[QuantTimes - 1];
	long double  *arrB = new long double[QuantTimes - 1];
	long double  *arrTemp0 = new long double[QuantTimes - 1];

	int quantTPerecl = -1;
	fncCalc_GradX1_and_GradX2_po_dTay__( parrTimes[QuantTimes -1] ,quantTPerecl, arrA, arrB  );
	///
   if (quantTPerecl != (QuantTimes - 2))
   {
	// .....
   }
   arrA [QuantTimes- 2] =  HomingSitTemp.mBombTraj.marrStrSK_VS [3] * cosl(HomingSitTemp.mBombTraj.marrStrSK_VS [4]) ;
   arrB [QuantTimes - 2] =  HomingSitTemp.mBombTraj.marrStrSK_VS [3] * sinl(HomingSitTemp.mBombTraj.marrStrSK_VS [4]) ;
   MatrxMultScalar(arrB, (QuantTimes - 1), 1, arrZ[QuantTimes - 1] ,arrTemp0);
   MtrxSumMatrx(arrTemp0, arrA ,(QuantTimes- 1), 1, arrFgr) ;
   arrFgr[ QuantTimes - 1] =  HomingSitTemp.mBombTraj.marrStrSK_VS [1];
	delete arrA;
	delete arrB ;
	delete arrTemp0;

 }

 // вычисление
void THomingSituation::fncCalc_GradX1_and_GradX2_po_dTay__( long double valTEnd,int &quantTPerecl
	,long double* GradX1_po_dTay,long double* GradX2_po_dTay  )
{
	int QuantTimes = mCntlFuncPar.miArrParams[0] ;
	long double *parrTimes =  mCntlFuncPar.mDblArrParams;
   long double* arrMtrxL = new long double [25 * (QuantTimes - 1)];
   long double* arrMtrxB = new long double [ 5 * (QuantTimes - 1)];
   fncMove_to_EndTime_and_CollectMtrxs(   arrMtrxL,  arrMtrxB, quantTPerecl  ) ;
   long double arrT0[5] ={0.}, arrInv[25] ={0.}, arrRez[5] ={0.};
   int ind = 0;
   ind =( QuantTimes  % 2 == 0 ?  -1 : 1);

   for (int i = 0; i < quantTPerecl; i++)
   {
	 InverseMtrx5(&arrMtrxL[ i * 25], arrInv) ;
	 MtrxMultMatrx(arrInv,5, 5, &arrMtrxB[i * 5],1, arrRez) ;
	 MtrxMultMatrx(&arrMtrxL[ quantTPerecl * 25],5, 5, arrRez,1, arrT0) ;
	// GradX1_po_dTay[i] = arrT0[0]* 2. * ((long double)ind) * mBomb.mMaxAt ;
	// GradX2_po_dTay[i] = arrT0[1]* 2. * ((long double)ind) * mBomb.mMaxAt ;
	 GradX1_po_dTay[i] = arrT0[0]* ((long double)ind)  ;
	 GradX2_po_dTay[i] = arrT0[1]* ((long double)ind)  ;
	 ind = -ind ;
   }
   delete  arrMtrxL;
   delete  arrMtrxB;

}
// вычисление
void THomingSituation::fncCalc_GradX1_and_GradX2_po_dTay( long double valTEnd,int &quantTPerecl
	,long double* GradX1_po_dTay,long double* GradX2_po_dTay  )
{
   int QuantTimes = mCntlFuncPar.miArrParams[0] ;
   THomingSituation HomingSituationTemp = *this ;
   HomingSituationTemp.mCntlFuncPar.mDblArrParams [QuantTimes -1]  =   valTEnd;

   long double* arrMtrxL = new long double [25 * (QuantTimes - 1)];
   long double* arrMtrxB = new long double [ 5 * (QuantTimes - 1)];
   HomingSituationTemp.fncMove_to_EndTime_and_CollectMtrxs(   arrMtrxL,  arrMtrxB, quantTPerecl  ) ;
   long double arrT0[5] ={0.}, arrInv[25] ={0.}, arrRez[5] ={0.};
   int ind = 0;
   ind =( QuantTimes  % 2 == 0 ?  -1 : 1);

   for (int i = 0; i < quantTPerecl; i++)
   {
	 InverseMtrx5(&arrMtrxL[ i * 25], arrInv) ;
	 MtrxMultMatrx(arrInv,5, 5, &arrMtrxB[i * 5],1, arrRez) ;
	 MtrxMultMatrx(&arrMtrxL[ quantTPerecl * 25],5, 5, arrRez,1, arrT0) ;
	 GradX1_po_dTay[i] = arrT0[0]* 2. * ((long double)ind) * mBombTraj.mBomb.mMaxAt ; // ???????!!!!!!!!!!!!
	 GradX2_po_dTay[i] = arrT0[1]* 2. * ((long double)ind) * mBombTraj.mBomb.mMaxAt ; // сравни с  fncCalc_GradX1_and_GradX2_po_dTay__
	 ind = -ind ;                                                                     // как правильно?
   }
   delete  arrMtrxL;
   delete  arrMtrxB;


}
// интгрирование уравнений движения УАБ, матрицы перехода,до момента  падения
// OUTPUT:
// 1. arrMtrxL [ mQuantTimes* 25] - массив матриц перехода из mTStart в
//   моменты переключений marrTimes[i], i = 1,..., mQuantTimes+1
//    последний момент - это момент точки падения
// 2. arrMtrxB [ mQuantTimes* 5]- массив векторов частных производных по
//    управлению AlfaAt в моменты переклюцчений  marrTimes[i], i = 1,..., mQuantTimes - 1
// 3. quantTPerecl - сколько получилось моментов переключений (так как УАБ могла упасть раньше последнего момента переключения)
void THomingSituation::fncMove_to_EndTime_and_CollectMtrxs( long double* arrMtrxL, long double* arrMtrxB, int &quantTPerecl  )
{
  int QuantTimes = mCntlFuncPar.miArrParams[0] ;
  long double *parrTimes =  mCntlFuncPar.mDblArrParams;
  quantTPerecl = 0;
  long double mtrxL [25] = {0.}, mtrxB[5] = {0.};
  for (int i = 0; i < 5; i++) mtrxL [i * 5 + i] = 1. ;


  int iCirc = 1000. / mBombTraj.mStepInt ;

  int i = 0 ;
  for ( i = 0; i < iCirc; i++)
  {
	  long double valU = fncCalcU(mTCur);
	  fncEilerStep_VS_and_MtrxL_and_MtrxB( mBombTraj.mStepInt, mtrxL, mtrxB );
	  if (mTCur > parrTimes[quantTPerecl + 1] )
	  {

		memcpy(&arrMtrxL[quantTPerecl * 25], mtrxL, 25 * sizeof(long double));
		memcpy(&arrMtrxB[quantTPerecl  * 5], mtrxB, 5 * sizeof(long double));

		quantTPerecl++;
	  }
	  if (mTCur > parrTimes[QuantTimes - 1] ) //??????
	  {
		quantTPerecl--;
		break ;
	  }
  }

}

// шаг етода эйлера на время valStepInt
// пересчитывает фазовый вектор - marrStrSK_VS [5]
// матрицу перехода -  mtrxL[25]
// возвращает матрицу частных производных правой части по управлению -  mtrxB [5]
 void THomingSituation::fncEilerStep_VS_and_MtrxL_and_MtrxB( const long double valStepInt
		 ,long double* mtrxL,long double* mtrxB )
 {
   long double arrF[5] ={0.}, arr_dF_po_dx[25] = {0}, arr_dF_po_dU[5] = {0.} ;
   // шаг интегрирования фазового вектора
   // arrF[5] -  правая часть сиситеты диф уравнений
// arr_dF_po_dx[25] -матрицы частных проихводных векторр фуеции правой части по фазовым переенным
// arr_dF_po_dU[5] - вектор частных производныз вектора правой части по переменной управления U (= AlfaP)
   long double valU = fncCalcU( mTCur);
   mBombTraj.fncCalc_F_and_dF_po_dx_and_dF_po_dU(valU , arrF, arr_dF_po_dx,arr_dF_po_dU ) ;
   MatrxMultScalar(arr_dF_po_dU, 1, 5, 2.* fncCalcU( mTCur-0.0001),mtrxB);
  // пересчет фазового вектора:
   long double arrT0[25] = {0.},arrT1[25] = {0.},arrT2[25] = {0.},arrT3[25] = {0.},arrT4[25] = {0.},arrT5[25] = {0.};
   MatrxMultScalar(arrF, 5, 1, valStepInt,arrT0); // f * dt
   MtrxSumMatrx(arrT0, mBombTraj.marrStrSK_VS,5, 1, arrT1) ;
   memcpy(mBombTraj.marrStrSK_VS, arrT1, 5 * sizeof(long double)) ;   ///
  ///

  // пересчет матрицы перехода
  MtrxMultMatrx(arr_dF_po_dx,5, 5, mtrxL,5, arrT2);  // arr_dF_po_dx * mtrxL
  MatrxMultScalar(arrT2, 25, 1, valStepInt,arrT3); // arr_dF_po_dx * dt
  MtrxSumMatrx(arrT3, mtrxL, 25, 1, arrT4) ;
  memcpy(mtrxL, arrT4, 25 * sizeof(long double)) ;   ///
  ///
  mTCur += valStepInt ;

 }

// нахождение оптимаального управления для задачи по дальности методом  перебора
// для одной точи переключения
// INPUT:
// valShagT - шаг сетки  по времени
THomingSituation THomingSituation::fncFindOptimalControl_for_Dist_MethodPerebora_1Point(TCntlFuncPar &CntlFuncPar
	   ,const long double valShagT, long double& valDGorizOpt)
{
   THomingSituation HomingSituation  = *this ;
   TCntlFuncPar CntlFuncParTemp =  CntlFuncPar ;
   long double *parrTimes = CntlFuncParTemp.mDblArrParams;
   parrTimes[0] = 0.;
   parrTimes[1] = 0.;
   parrTimes[2] = 1000.;
   CntlFuncParTemp.miArrParams[0] = 3; // к-во точек переключения
   THomingSituation HomingSituationReturn  = *this ;
   long double valDGoroz = 0.;
   HomingSituation.fncMoveBomb_TO_ZeroAlt_AND_ShowGraphs(NULL, valDGoroz) ;

   long double valT = HomingSituation.mBombTraj.mTCur;
   long double valTOpt = 0.;
   valDGorizOpt = valDGoroz;
   long double valTCur = 0., valDGorizCur = 0.;
   int nCircle = valT / valShagT;

   for (int  i=0; i < nCircle; i++)
   {
	THomingSituation HomingSituationCur = *this;// ( mTet0, mV0, mAltit, mX0, mStepInt, miControlType, marrTimes, 3, mBomb) ;
	HomingSituationCur.mCntlFuncPar.miArrParams[0] = 3 ;
	valTCur = ((long double)i) *valShagT ;
	HomingSituationCur.mCntlFuncPar.mDblArrParams[0] = 0.;
	HomingSituationCur.mCntlFuncPar.mDblArrParams[1] =  valTCur;
	HomingSituationCur.mCntlFuncPar.mDblArrParams[2] = 1000.;

   HomingSituationCur.fncMoveBomb_TO_ZeroAlt_AND_ShowGraphs(NULL, valDGorizCur) ;
   if (valDGorizCur > valDGorizOpt )
   {
	valDGorizOpt = valDGorizCur;
	valTOpt =  valTCur ;
	HomingSituationReturn = HomingSituationCur ;
   }
  // else
  // {
  //	  if ((i > 1) && (valDGorizCur < (valDGorizOpt -2.)) )break ;
  // }
 }

  return HomingSituationReturn;
}

////////////////////////////////////////////////////////////////////////////////////
//////// ФУНКЦИИ НАХОЖДЕНИЕ МАКСИМАЛЬНОЙ ДАЛЬНОСТИ ЗАХВАТА ЦЕЛИ ГОЛОВКОЙ ////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// функция нахождениея асимальной дальности захвата цели
// для заданной высоты бросания и скорости бросания ищется угол бросания и закон
// управления, обеспечивающие аксимальную дальность захвата цели ГСН
// Закон управления  - кусочно-постоянный
//  iQuantPointsOverswitch - к-во переключений
// valShagT - шаг по вреени при переборе точек переключения
// результаты сохраняются в HomingSituationReturn

bool THomingSituation::fnc_FindMaxDist_For_Takeover_Perebor(const int iQuantPointsOverswitch,const long double valShagT
   , THomingSituation &HomingSituationReturn0, long double &valMaxDistTakeover)
{
	valMaxDistTakeover = -1.;
	bool breturn = false;
	const int iNum = 20;
	long double valShagAng = 90./((long double)iNum) / 180. * M_PI ;

	for (int i = 0; i <  iNum-2 ; i++)

	{
	  long double valAng =    valShagAng * ((long double)(i +1));
	  if (valAng  > 85. /180. * M_PI)
	  {
        break;
	  }
	  THomingSituation  HomingSituationCur = *this;
	  THomingSituation  HomingSituationReturn = *this;
	TBombTraj BombTraj ( valAng ,mBombTraj.mV0 ,mBombTraj.mAltit  ,mBombTraj.mX0 	,mBombTraj.mStepInt  ,mBombTraj.mTStart ,mBombTraj.mBomb );
	HomingSituationCur.mBombTraj = BombTraj;
	HomingSituationCur.mCntlFuncPar.miControlType = 1;
	THomingSituation HomingSituationCur1 = HomingSituationCur;
	HomingSituationCur1.mCntlFuncPar.miControlType = 2;
	bool brez = false,brez1 = false;
	long double valDist = -1., valDist1 = -1.;

	switch( iQuantPointsOverswitch)
	{
	case 1:
	brez = HomingSituationCur.fnc_FindMaxDist_For_Takeover_For_FixedTetta0_Perebor_1Point(valShagT, HomingSituationReturn, valDist);
	 if (brez)
	 {
		if (valDist > valMaxDistTakeover)
		{
		 valMaxDistTakeover =valDist ;
		 HomingSituationReturn0 = HomingSituationReturn ;
		}
		breturn = true;
	 }
	  brez1 = HomingSituationCur1.fnc_FindMaxDist_For_Takeover_For_FixedTetta0_Perebor_1Point(valShagT,HomingSituationReturn, valDist);
	   if (brez1)
	   {
		if (valDist > valMaxDistTakeover)
		{
		 valMaxDistTakeover =valDist ;
		 HomingSituationReturn0 = HomingSituationReturn ;
		}
		breturn = true;
	   }

	break;

	case 2:
	 brez = HomingSituationCur.fnc_FindMaxDist_For_Takeover_For_FixedTetta0_Perebor_2Point(valShagT, HomingSituationReturn, valDist);
	  if (brez)
	   {
		if (valDist > valMaxDistTakeover)
		{
		 valMaxDistTakeover =valDist ;
		 HomingSituationReturn0 = HomingSituationReturn ;
		}
		breturn = true;
	 }
	  brez1 = HomingSituationCur1.fnc_FindMaxDist_For_Takeover_For_FixedTetta0_Perebor_2Point(valShagT, HomingSituationReturn, valDist);
      if (brez1)
	   {
		if (valDist > valMaxDistTakeover)
		{
		 valMaxDistTakeover =valDist ;
		 HomingSituationReturn0 = HomingSituationReturn ;
		}
		breturn = true;
	 }
	break;

	case 3:
	brez = HomingSituationCur.fnc_FindMaxDist_For_Takeover_For_FixedTetta0_Perebor_3Point(valShagT, HomingSituationReturn, valDist);
       if (brez1)
	   {
		if (valDist > valMaxDistTakeover)
		{
		 valMaxDistTakeover =valDist ;
		 HomingSituationReturn0 = HomingSituationReturn ;
		}
		breturn = true;
	 }
	  brez1 = HomingSituationCur1.fnc_FindMaxDist_For_Takeover_For_FixedTetta0_Perebor_3Point(valShagT, HomingSituationReturn, valDist1);
	  if (brez1)
	   {
		if (valDist > valMaxDistTakeover)
		{
		 valMaxDistTakeover =valDist ;
		 HomingSituationReturn0 = HomingSituationReturn ;
		}
		breturn = true;
	 }
	break;

	 default:

	return false;
	}


   }
	return breturn ;

}

// функция нахождениея асимальной дальности захвата цели
// для заданной высоты бросания , заданного начального угла бросания Тетта0 и скорости бросания ищется  закон
// управления, обеспечивающие аксимальную дальность захвата цели ГСН
// Закон управления  - кусочно-постоянный с одной точкой переключения
// valShagT - шаг по вреени при переборе точек переключения
// результаты сохраняются в HomingSituationReturn
 bool THomingSituation::fnc_FindMaxDist_For_Takeover_For_FixedTetta0_Perebor_1Point(const long double valShagT
   , THomingSituation &HomingSituationReturn0, long double &valMaxDistTakeover)
{
   bool breturn = false ;
   valMaxDistTakeover =  -1. ;
   THomingSituation HomingSituation  = *this ;
   TCntlFuncPar CntlFuncParTemp =  HomingSituation.mCntlFuncPar ;
   long double *parrTimes = CntlFuncParTemp.mDblArrParams;
   parrTimes[0] = 0.;
   parrTimes[1] = 0.;
   parrTimes[2] = 1000.;
   CntlFuncParTemp.miArrParams[0] = 3; // к-во точек переключения
   THomingSituation HomingSituationReturn  = *this ;
   long double valDGoroz = 0.;
   HomingSituation.fncMoveBomb_TO_ZeroAlt_AND_ShowGraphs(NULL, valDGoroz) ;

   long double valT = HomingSituation.mBombTraj.mTCur;
   long double valTOpt = 0.;

   long double valTCur = 0.;
   int nCircle = valT / valShagT;
   long double	valDistTakeover = -1. ;
   for (int  i=0; i < nCircle; i++)
   {
	THomingSituation HomingSituationCur = *this;// ( mTet0, mV0, mAltit, mX0, mStepInt, miControlType, marrTimes, 3, mBomb) ;
	HomingSituationCur.mCntlFuncPar.miArrParams[0] = 3 ;
	valTCur = ((long double)i) *valShagT ;
	HomingSituationCur.mCntlFuncPar.mDblArrParams[0] = 0.;
	HomingSituationCur.mCntlFuncPar.mDblArrParams[1] =  valTCur;
	HomingSituationCur.mCntlFuncPar.mDblArrParams[2] = 1000.;

   if(HomingSituationCur.findTakeoverDist_AND_ShowGraphs(NULL,valDistTakeover))
   {
	 if (valDistTakeover > valMaxDistTakeover)
	 {
	  HomingSituationCur.mCntlFuncPar.mDblArrParams[2] = HomingSituationCur.mBombTraj.mTCur;
	  memset(HomingSituationCur.mShipTarg.mTraject.marrVectSostGSK_Begin, 0, 9 * sizeof(double));
	  HomingSituationCur.mShipTarg.mTraject.marrVectSostGSK_Begin[0] =  (double)valDistTakeover;
	  valMaxDistTakeover = valDistTakeover ;
	  HomingSituationReturn0 =  HomingSituationCur ;
	  breturn = true;
	 }


   }

 }

  return breturn;
}

// интегрирование уравнений движения УАБ до момента времени , кргда выполняются условия
// осушществления самонаведения и захвата цели на сопровождение ГСН
 bool THomingSituation::findTakeoverDist_AND_ShowGraphs(wchar_t *wcharrPath ,long double &valDistTakeover)
 {
	bool breturn = false;
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
	wcscpy( &pwcharrFileNames[ 6 * 30], L"AlfAt");
	wcscpy( &pwcharrFileNames[ 7 * 30], L"M");
	wcscpy( &pwcharrFileNames[ 8 * 30], L"Cy");
	wcscpy( &pwcharrFileNames[ 9 * 30], L"CyqSmUmaxMG");

	pscaleY = new double  [QUANT_COLS] ;
	pscaleY[1] = 0.1;
	pscaleY[2] = 0.1;

	pscaleY[3] = 1000.;
	pscaleY[4] = 1.;
	pscaleY[5] = 1000.;
	pscaleY[6] = 1000.;
	pscaleY[7] = 1000.;
	pscaleY[8] = 1000.;
	pscaleY[9] = 1000.;
}

  int iCirc = 1000. / mBombTraj.mStepInt ;
 // long double valTTemp = mTCur + ((long double)iCirc) * mStepInt ;
  int i = 0 ;
  int iNupPointsOut = 0;
  double valTOut = -DEL_T ;
  for ( i = 0; i < iCirc; i++)
  {

	//  if (mTCur > 9.3099999)
	  //{
		// int iii = 0;
	//  }
	  long double valU = fncCalcU( mTCur)  ;

	  // шаг интегрования
	  mBombTraj.fncEilerStep(valU, mBombTraj.mStepInt);
	  mTCur = mBombTraj.mTCur;
	  ///


	  if  ( iNupPointsOut < QUANT_POINTS_MAX )
	  {

	  valTOut = (double)mTCur ;

	   if (wcharrPath)
	   {
	   double *p = &parrBuff[iNupPointsOut *  QUANT_COLS] ;
	   p[0] =  (double)mTCur ;
	   p[1] = (double)mBombTraj.marrStrSK_VS [0];
	   p[2] = (double)mBombTraj.marrStrSK_VS [1];
	   p[3] = (double)mBombTraj.marrStrSK_VS [2];
	   p[4] = (double)mBombTraj.marrStrSK_VS [3];
	   p[5] = (double)mBombTraj.marrStrSK_VS [4];
	   p[6] = (double)fncCalcU( mTCur);
	   long double valTay,  valDerivTay =0., valMach = -1., arrGradMach[5] = {0.};

	   fncCalcNormTemperature(p[2], valTay, valDerivTay) ;
	   mBombTraj.fncCalcMach_and_GradMach( valTay,  valDerivTay
	 ,valMach, arrGradMach);
	   p[7] = (double)valMach;
	   long double  val_q = 0., arrGrad_q [5] = {0.};
	   mBombTraj.fncCalc_q_and_Grad_q(val_q, arrGrad_q) ;
	   p[8] = mBombTraj.mBomb.fncCy(valMach);
	   p[9] = p[8] * val_q * mBombTraj.mBomb.mMidSq * mBombTraj.mBomb.mMaxAt / mBombTraj.mBomb.mMass / G_ZEMLI ;

		}

	  iNupPointsOut++;

	  }

	 //  проверка условий захвата
		// высота > 0
	 if ((mBombTraj.marrStrSK_VS[1] <0.)|| (mBombTraj.marrStrSK_VS [4] >= (M_PI/2. - 0.01)))
	 {
		if (parrBuff) delete parrBuff;
		if (pscaleY) delete pscaleY;
		if (pwcharrFileNames) delete pwcharrFileNames;
		return false ;
	 }
	 if (mBombTraj.marrStrSK_VS [4] >=0.) continue ;
	 if (mBombTraj.marrStrSK_VS [1] >= mBombTraj.mBomb.mHomingHead.mR  ) continue ;


	 ///
		// сектор обзора пересекается с линией гтризонта (верхний луч пересекается с горизонтом)  3333
		 if (IsGeom(valDistTakeover))
		 {
		   breturn = IsHoming();
		   break;
		 }



  }
 if (!breturn) return false;


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

		wchar_t wchFileName4[300] = {0} ;
		wcscpy(wchFileName4, wcharrPath );
		wcscat(wchFileName4, L"\\Axes.shp");
		TYrWriteShapeFile::CreateShpAxes(wchFileName4,-80000.,80000,-80000.,80000) ; ;

		delete parrBuff ;
		delete pwcharrFileNames ;
		delete pscaleY ;
   }

	return true ;
}


// интегрирование уравнений движения УАБ до момента времени , кргда выполняются условия
// осушществления самонаведения и захвата цели на сопровождение ГСН
 bool THomingSituation::findTakeoverDist_AND_CreatePictures(wchar_t *wcharrPath
	 , long double &valDistTakeover)
 {
	drawJetPictures( wcharrPath);
	THomingSituation HomingSituationCopy0 = *this ;

	wchar_t pwcharrPath[900] = {0};
	for (int i =0; i < 3; i++)
	{
		wcscpy(&pwcharrPath[i *300], wcharrPath) ;
		wcscat(&pwcharrPath[i *300], L"\\Part");
		wchar_t string[5] = {0};
		_itow(i, string, 10);
		wcscat(&pwcharrPath[i *300], string);

	}
	HomingSituationCopy0.findTakeoverDist_AND_ShowGraphs(wcharrPath ,valDistTakeover);
	bool bend = false;
	THomingSituation HomingSituationCur = *this ;
	int  iControlType = mCntlFuncPar.miControlType;
	int  iNotControlType = (iControlType == 1)?2:1;
	int  iControlTypeCur = ( mCntlFuncPar.miArrParams[0]  % 2 == 0 ? iControlType :iNotControlType);

	for (int i = 0; i < (mCntlFuncPar.miArrParams[0] - 1); i++)
	{

	 long double valFixedTime = 0.;
	 long double valFixedEndTime = mCntlFuncPar.mDblArrParams[i + 1];
	 if (mCntlFuncPar.mDblArrParams[i + 1] < HomingSituationCopy0.mTCur)
	 {
	   valFixedTime = ( mCntlFuncPar.mDblArrParams[i + 1] + mCntlFuncPar.mDblArrParams[i ])/ 2.;
	 }
	 else
	 {
	  valFixedTime = (HomingSituationCopy0.mTCur + mCntlFuncPar.mDblArrParams[i ])/ 2.;
	  valFixedEndTime = HomingSituationCopy0.mTCur;
	  bend = true;
	 }
	  _wmkdir( &pwcharrPath[i *300]);

	 THomingSituation HomingSituationCopy1 = *this ;// для рисования картинки УАБ и корабля посередине участка
	 THomingSituation HomingSituationCopy2 = *this ;
	 HomingSituationCopy1.fncMoveClass_TO_FixedTime_AND_ShowGraphs( NULL,  valFixedTime );
	 const double valRastigenie = 300.;

	 wchar_t pwcharrPathPositMiddle[300] = {0},pwcharrPathPositEnd[300] = {0};
	 wcscpy(pwcharrPathPositMiddle ,&pwcharrPath[i *300]) ;
		wcscat(pwcharrPathPositMiddle, L"\\MiddlePoint");
		_wmkdir( pwcharrPathPositMiddle);

	 HomingSituationCopy1.drawPictures( pwcharrPathPositMiddle,valRastigenie);

	  wcscpy(pwcharrPathPositEnd ,&pwcharrPath[i *300]) ;
		wcscat(pwcharrPathPositEnd, L"\\EndPoint");
		_wmkdir( pwcharrPathPositEnd);
	 HomingSituationCur.fncMoveClass_TO_FixedTime_AND_ShowGraphs( &pwcharrPath[i *300],  valFixedEndTime );
	 HomingSituationCur.drawPictures( pwcharrPathPositEnd ,valRastigenie);


	 if (bend )
	 {
	  // построение графиков для цели
		   // положение цели на момент пеленга
		   wchar_t wchFileName4[300] ={0};
		TURPointXY pntTargPos(valDistTakeover, 0.);
		wcscpy(wchFileName4, wcharrPath );
		wcscat(wchFileName4, L"\\pntTargPosBegin.shp");
		pntTargPos.WriteSetSHPFiles(wchFileName4, &pntTargPos, 1 );
		   ///

		   // отрезок достижимости
		   TURPointXY arrTargPos[2];
		   arrTargPos[0] =  TURPointXY(valDistTakeover - 20.* valFixedEndTime, 0.);
		   arrTargPos[1] =  TURPointXY(valDistTakeover + 20.* valFixedEndTime, 0.);
		   TURPolyLine plnArea(arrTargPos,2) ;
		   wcscpy(wchFileName4, wcharrPath );
		  wcscat(wchFileName4, L"\\plnTargArea.shp");
		  plnArea.WriteSetSHPFiles(wchFileName4, &plnArea, 1 );

		  // построение картинок с кораблем
		 // TShipTraj ShipTrajPlus = mShipTarg.mTraject;
		 // ShipTrajPlus.mTCur =  valFixedEndTime;
		  //ShipTrajPlus.marrVectSostGSK[0] = ShipTrajPlus.marrVectSostGSK_Begin[0] +  20.* valFixedEndTime;
		  double valPos = valDistTakeover +  20.* valFixedEndTime;
		  double valVeloc =  20.;
		 // ShipTrajPlus.marrVectSostGSK[3] = 20.;
		  wcscpy(wchFileName4, wcharrPath );
		  wcscat(wchFileName4, L"\\VesselPlus.shp");
		  //_wmkdir(wchFileName4);

		 const double valRastigenieVess = 4.;
		 drawVessel( wchFileName4,  valPos,  valVeloc,  valRastigenieVess ) ;

		 // ShipTrajPlus.marrVectSostGSK[0] = ShipTrajPlus.marrVectSostGSK_Begin[0] -  20.* valFixedEndTime;
		 // ShipTrajPlus.marrVectSostGSK[3] = -20.;
		 valPos = valDistTakeover - 20.* valFixedEndTime;
		 valVeloc =  -20.;
		  wcscpy(wchFileName4, wcharrPath );
		  wcscat(wchFileName4, L"\\VesselMinus.shp");

		  drawVessel( wchFileName4,  valPos,  valVeloc,  valRastigenieVess ) ;

	  break;
	 }

	 TBombTraj BombTraj  (
	HomingSituationCur.mBombTraj.marrStrSK_VS [4]  // нач угол бросания
	,HomingSituationCur.mBombTraj.marrStrSK_VS [3]  // нач скоргсть
	,HomingSituationCur.mBombTraj.marrStrSK_VS [1]  // нач высота
	 ,HomingSituationCur.mBombTraj.marrStrSK_VS [0]// нач. положение по оси X
	,HomingSituationCur.mBombTraj.mStepInt  // шагт интьерирования по времени диф уравнений движения
	 ,HomingSituationCur.mBombTraj.mTCur
	,mBombTraj.mBomb
	,HomingSituationCur.mBombTraj.marrCfW
	)  ;
	int iArrParams[NUM_INT_PARAMS] ={0};
	iArrParams[0] = 2;
	long double DblArrParams [NUM_LD_PARAMS] = {0.};
	DblArrParams[0] = valFixedEndTime;
	DblArrParams[1] = 1000.;
	iControlTypeCur = (iControlTypeCur == iControlType ?iNotControlType  : iControlType);
	//int icur = iControlType;
   //	iControlType =  iNotControlType;
   //	iNotControlType = icur;
   //	iControlTypeCur = iControlType;

	 TCntlFuncPar CntlFuncPar( iControlTypeCur, iArrParams , DblArrParams);
	 THomingSituation HomingSituationCur1(HomingSituationCur.mShipTarg
		,  BombTraj //
		,  CntlFuncPar
		,mGlonass
		,HomingSituationCur.mTCur
		,NULL // путь у папке с отчетом
		);
	 HomingSituationCur =  HomingSituationCur1;


	}
   *this = HomingSituationCopy0;
	return true ;
}




// создание картиннок с самолетом
void THomingSituation:: drawJetPictures( wchar_t *wcharrPath)
{
  // нахождение радиуса вращениия самолета
  double valR = mBombTraj.mV0 * mBombTraj.mV0 / G_ZEMLI;
  // нахождение центра вращения самолета
  double valTet1 = mBombTraj.mTet0 + M_PI /2.;
  TURPointXY PntCentre( (double)(mBombTraj.mX0 + cos(valTet1) * valR  ), (double)(mBombTraj.mAltit + sin(valTet1) * valR));
  ///

  // создание сектора
  double valFi0 = -M_PI /2.;
  double valFi1 = 0. ;
  TSector Sector(  PntCentre, valR, valFi0,  valFi1 );
  //
  // рисование сектора
  wchar_t wchFileName[300] = {0};
  wchar_t wchFileNameOut[300] = {0};
  wcscpy(wchFileName, wcharrPath );
  wcscat(wchFileName, L"\\JetTraj.shp");
  Sector.ShowMe(wchFileName) ;
  ///

   //считывание файла с центром

 // wcscpy(wchFileName, pwchShapeJet );
 // wcscat(wchFileName, L"\\CentreJet.shp");
 // TURPointXY *pPnt = (TURPointXY*) malloc( sizeof(TURPointXY));
 // TURPointXY **ppPnt = &pPnt;
 // int quantPnt = 1;
 // TURPointXY::ReadSHPFile(wchFileName, ppPnt,  &quantPnt) ;
  ///

  // считывание файла plg
//  wcscpy(wchFileName, pwchShapeJet );
 // wcscat(wchFileName, L"\\plgJet.shp");

  wcscpy(wchFileNameOut, wcharrPath );
  wcscat(wchFileNameOut, L"\\Jet.shp");

 // TURPolygon plgJet = createPictJet();

  TURPointXY pntSdvig( mBombTraj.marrStrSK_VS [0], mBombTraj.marrStrSK_VS [1]);
  const double valRastigenie = 40.;

  drawJet( wchFileNameOut,mBombTraj.marrStrSK_VS [4],  pntSdvig, valRastigenie);

  // построение первой картинки
 // createTransformedPlgShapeFile(wchFileName, wchFileNameOut, mBombTraj.marrStrSK_VS [4]
 //	, (*ppPnt)[0],  pntSdvig, valRastigenie) ;
 // TURPolygon plgJet1=  plgJet.LinTransform(mBombTraj.marrStrSK_VS [4] , pntSdvig, valRastigenie );
 // plgJet1.WriteSetSHPFiles(wchFileNameOut,&plgJet1, 1) ;
  // построение второй картинки
   wcscpy(wchFileNameOut, wcharrPath );
   wcscat(wchFileNameOut, L"\\Jet1.shp"); // куда класть шейп файл

   //  точка на траектории
   double valTet2 = M_PI/ 2.- 0.3;
   pntSdvig.X = PntCentre.X + valR * sin(valTet2);
   pntSdvig.Y = PntCentre.Y - valR * cos(valTet2);



 // TURPolygon plgJet2=  plgJet.LinTransform(valTet2, pntSdvig, valRastigenie );
  //plgJet2.WriteSetSHPFiles(wchFileNameOut,&plgJet2, 1) ;
  drawJet( wchFileNameOut,valTet2,  pntSdvig, valRastigenie);

  //	 free(pPnt);
}


// создание картиннок с подрывом ЭМБ

void THomingSituation:: drawEMB_Blast_Pictures(wchar_t *pwcharrPath)
{
  // создание картинки с сектор поражения



 double arrShipVeloc[2] ={0.};
  arrShipVeloc[0] = mFilt.marrEstTargX[1];
  double arrDelS [2] = {0.};
  arrDelS [1] = - mBombTraj.marrStrSK_VS[1];
  arrDelS [0] = mFilt.marrEstTargX[0] - mBombTraj.marrStrSK_VS[0];
   double valAngVis = atan2(arrDelS[1], arrDelS[0]);
  ///
   TURPointXY arrPnt[2];
  // точка текущей траектории
  arrPnt[0].X = mBombTraj.marrStrSK_VS[0];
  arrPnt[0].Y = mBombTraj.marrStrSK_VS[1];
  ///

  // точка на касательной к траектории УАБ
   arrPnt[1].X = mBombTraj.marrStrSK_VS[0]  + mBombTraj.mBomb.mEbomb.mR * cos(valAngVis);
	arrPnt[1].Y = mBombTraj.marrStrSK_VS[1] + mBombTraj.mBomb.mEbomb.mR * sin(valAngVis);
   ///

   TURPolyLine PlnDirect(arrPnt,2) ;
	wchar_t wchFileNamePlnDirect[300] = {0};
	wcscpy(wchFileNamePlnDirect, pwcharrPath );
	wcscat(wchFileNamePlnDirect, L"\\DirectLine.shp");
	PlnDirect.WriteSetSHPFiles(wchFileNamePlnDirect, &PlnDirect, 1) ;

	 // картинка с сектором обзора координатора ГСН , ось УАБ совпадаеть с вектром скорости

   double valFi0 = valAngVis - mBombTraj.mBomb.mEbomb.mFi/ 2.;
   double valFi1 = valAngVis + mBombTraj.mBomb.mEbomb.mFi/ 2.;
   TURPolygon PlgSectorDefeat = TURPolygon::fncCreateSector(arrPnt[0], mBombTraj.mBomb.mEbomb.mR,valFi0,valFi1,400) ;
   wchar_t wchFileNameSectorDefeat[300] = {0};

   wcscpy(wchFileNameSectorDefeat, pwcharrPath );

   wcscat(wchFileNameSectorDefeat, L"\\DefeatSector.shp");
   PlgSectorDefeat.WriteSetSHPFiles(wchFileNameSectorDefeat, &PlgSectorDefeat, 1) ;
   ///

   // построение молнии

	// путь к результирующему файлу
	 wchar_t pwchShapeFileOut[300] = {0} ;
	 wcscpy(pwchShapeFileOut, pwcharrPath );
	wcscat(pwchShapeFileOut, L"\\Fireball.shp");

	TURPointXY pntSdvig( mBombTraj.marrStrSK_VS [0], mBombTraj.marrStrSK_VS [1]);

	double valRastigenie = 0.3;
 // createTransformedPlgShapeFile(pwchShapeFileInp, pwchShapeFileOut,   -valAng0 +valAngVis
   //	, pntBegin,  pntSdvig,valRastigenie) ;
	drawFireBall( pwchShapeFileOut,valAngVis,  pntSdvig, valRastigenie) ;
	///

	// УАБ

	// соэдание путей к результирующему файлу с УАБ
  wchar_t wchFileNameRez0[300] = {0} ;
	wcscpy(wchFileNameRez0, pwcharrPath );
	wcscat(wchFileNameRez0, L"\\UABUnionRez00.shp");

	valRastigenie = 8.;

	drawUAB( wchFileNameRez0,valAngVis,  pntSdvig, valRastigenie) ;
   ///

}

// создание преобразованного шейп файла
//INPUT:
// pwchShapeInp - путь где лежит исходный шейп файл ...\file.shp
// wcharrPathOut - путь где лежит результирующий шейп файл ...\file.shp
// valAng - угол поворота
// pntCentre - точка центра полигона
// pntSdvig - точка куда перемещается центр полигона
// valRastigenie - коэффициент растяжения
void THomingSituation::createTransformedPlgShapeFile(wchar_t *pwchShapeFileInp, wchar_t *pwchShapeFileOut, const double  valAng
	, const TURPointXY pntCentre, const TURPointXY pntSdvig,const double valRastigenie)
{
	// считывание файла plg

	TURPolygon *pPlg = (TURPolygon*)malloc( sizeof(TURPolygon));
	int quantPlg = 1;
	TURPolygon **ppPlg = & pPlg;
	TURPolygon::ReadSHPFile(pwchShapeFileInp,ppPlg,  & quantPlg) ;
	///
  

	// копировангие полигона
	TURPolygon *pPlgTansf = new TURPolygon [quantPlg];
	memcpy( pPlgTansf, (*ppPlg), quantPlg * sizeof(TURPolygon));
	///

	// преобразование
	for (int i =0; i < quantPlg; i++)
	{
	TURPolygon plgCur = pPlgTansf[i].LinTransform(valAng,pntCentre,  pntSdvig,valRastigenie ) ;
	pPlgTansf[i] =  plgCur;
	}
	///
	 TURPolygon plnUn( pPlgTansf, quantPlg);
	 TURPolygon::WriteSetSHPFiles(pwchShapeFileOut,&plnUn,  1) ;
	// соэдание шейп файла
   //	TURPolygon::WriteSetSHPFiles(pwchShapeFileOut,pPlgTansf,  quantPlg) ;
	free(pPlg);
	delete [] pPlgTansf;
}

// проверка геометрии сектора обзора ГСН
// и проверка того, что область достижимости цели лежтит внутри сектора обзора координатора ГСН
 bool THomingSituation:: IsGeom(long double &valDistTakeover)
 {
 try
 {
	   long double valTet1= mBombTraj.marrStrSK_VS [4] +  mBombTraj.mBomb.mHomingHead.mFi /2. ;
	   if (valTet1 >=0.) return false;
	  // long double valx1 = mBombTraj.marrStrSK_VS [0] - mBombTraj.marrStrSK_VS [1]/ tanl( valTet1);
	  long double valx1 = - mBombTraj.marrStrSK_VS [1]/ tanl( valTet1);
	   long double vald1 = sqrtl(valx1 * valx1 + mBombTraj.marrStrSK_VS [1] * mBombTraj.marrStrSK_VS [1]);
	   bool bgeom1 = (vald1 <= mBombTraj.mBomb.mHomingHead.mR);
	   if (!bgeom1) return false;
	   valx1 += mBombTraj.marrStrSK_VS [0]  ;

	   long double valTet2= mBombTraj.marrStrSK_VS [4] -  mBombTraj.mBomb.mHomingHead.mFi /2. ;
	   if (valTet2 <= -M_PI) return false;
	  // long double valx2 = mBombTraj.marrStrSK_VS [0] - mBombTraj.marrStrSK_VS [1]/ tanl( valTet2);
	   long double valx2 = - mBombTraj.marrStrSK_VS [1]/ tanl( valTet2);
	   long double vald2 = sqrtl(valx2 * valx2 + mBombTraj.marrStrSK_VS [1] * mBombTraj.marrStrSK_VS [1]);
	   bool bgeom2 = (vald2 <= mBombTraj.mBomb.mHomingHead.mR);
	   if (!bgeom2) return false;
	   valx2 += mBombTraj.marrStrSK_VS [0]  ;


	   long double valDistTakeoverTemp =  mBombTraj.marrStrSK_VS [0]  - mBombTraj.marrStrSK_VS [1] / tanl(mBombTraj.marrStrSK_VS [4]);
	   long double valxmax = valDistTakeoverTemp + 20.  * mBombTraj.mTCur;
	   long double valxmin = valDistTakeoverTemp - 20.  * mBombTraj.mTCur;
	   bool barea =  (valxmax < valx1 )&& (valxmin > valx2 ) ;
	   if (!barea) return false;
	   valDistTakeover = valDistTakeoverTemp ;
	   return true;
 }
 catch(...)
 {
     int ii =0;
 }

 }

	// проверка возможности самонаведения
  bool THomingSituation:: IsHoming()
  {

	long double valTay,  valDerivTay =0., valMach = -1., arrGradMach[5] = {0.};
	fncCalcNormTemperature(mBombTraj.marrStrSK_VS [1], valTay, valDerivTay) ;
	mBombTraj.fncCalcMach_and_GradMach( valTay,  valDerivTay
	,valMach, arrGradMach);

	long double  val_q = 0., arrGrad_q [5] = {0.};
	mBombTraj.fncCalc_q_and_Grad_q(val_q, arrGrad_q) ;
	long double temp =  mBombTraj.mBomb.mMaxAt
	* val_q *mBombTraj.mBomb.fncCy(valMach)* mBombTraj.mBomb.mMidSq /  mBombTraj.mBomb.mMass / G_ZEMLI;
	long double temp1 =  fabsl(sinl(mBombTraj.marrStrSK_VS [4]-M_PI/ 2.) ) ;
	return (temp1 <=  temp ) ;
  }



// функция нахождениея асимальной дальности захвата цели
// для заданной высоты бросания , заданного начального угла бросания Тетта0 и скорости бросания ищется  закон
// управления, обеспечивающие аксимальную дальность захвата цели ГСН
	// Закон управления  - кусочно-постоянный с двумя точками переключения
// valShagT - шаг по вреени при переборе точек переключения
// результаты сохраняются в HomingSituationReturn
 bool THomingSituation::fnc_FindMaxDist_For_Takeover_For_FixedTetta0_Perebor_2Point(const long double valShagT
   , THomingSituation &HomingSituationReturn0, long double &valMaxDistTakeover)
{

	valMaxDistTakeover = -1.;
	bool breturn  = false ;
	TCntlFuncPar CntlFuncPar;
	CntlFuncPar.miControlType = mCntlFuncPar.miControlType ;
	CntlFuncPar.miArrParams[0] = 4;
	CntlFuncPar.mDblArrParams[0] = 0.;
	CntlFuncPar.mDblArrParams[1] = 0.;
	CntlFuncPar.mDblArrParams[2] = 0.;
	CntlFuncPar.mDblArrParams[3] = 1000.;

	THomingSituation HomingSituation0(mShipTarg, mBombTraj, CntlFuncPar,mTCur ,NULL);
	THomingSituation HomingSituation1 = HomingSituation0 ;
	long double valDGoriz = -1.;
	HomingSituation1.fncMoveBomb_TO_ZeroAlt_AND_ShowGraphs(NULL, valDGoriz);

   long double valT = HomingSituation1.mBombTraj.mTCur;
   long double valTOpt = 0.;

   long double valTCur = 0., valDGorizCur = 0.;
   int nCircle = valT / valShagT;

   for (int  i=0; i < nCircle; i++)
   {
	THomingSituation HomingSituationCur =  HomingSituation0;
	valTCur = ((long double)i) *valShagT ;
	HomingSituationCur.mCntlFuncPar.mDblArrParams[0] = 0.;
	HomingSituationCur.mCntlFuncPar.mDblArrParams[1] = valTCur;
	HomingSituationCur.mCntlFuncPar.mDblArrParams[0] = 1000.;
	HomingSituationCur.fncMoveClass_TO_FixedTime_AND_ShowGraphs( NULL, HomingSituationCur.mCntlFuncPar.mDblArrParams[1] ) ;


	TBombTraj BombTrajCur1(	mBombTraj.marrStrSK_VS [4] 	,mBombTraj.marrStrSK_VS [3] ,mBombTraj.marrStrSK_VS [1] ,mBombTraj.marrStrSK_VS [0]  // нач. положение по оси X
	,mBombTraj.mStepInt  ,0.,mBombTraj.mBomb );
	TCntlFuncPar CntlFuncParCur;
	CntlFuncParCur.miControlType = mCntlFuncPar.miControlType ;
	CntlFuncParCur.miArrParams[0] = 3;
	CntlFuncParCur.mDblArrParams[0] = 0.;
	CntlFuncParCur.mDblArrParams[1] = 0.;
	CntlFuncParCur.mDblArrParams[2] = 1000.;
	THomingSituation HomingSituationCur1(mShipTarg , BombTrajCur1, CntlFuncParCur,0.,NULL );
	THomingSituation HomingSituationCur2 =  HomingSituationCur1;
	THomingSituation HomingSituationReturn =   HomingSituationCur1;
	long double valMaxDist = -1.;
	  bool brez1 = HomingSituationCur1.fnc_FindMaxDist_For_Takeover_For_FixedTetta0_Perebor_1Point(valShagT, HomingSituationCur2, valMaxDist);
		if (brez1)
		{
		  if (valMaxDist > valMaxDistTakeover)
		  {
		  valMaxDistTakeover = valMaxDist;
		  HomingSituationReturn0 =  HomingSituation0 ;
		  HomingSituationReturn0.mCntlFuncPar.miArrParams[0] = 4;
		  HomingSituationReturn0.mCntlFuncPar.mDblArrParams[0] = 0.;
		  HomingSituationReturn0.mCntlFuncPar.mDblArrParams[1] = valTCur ;
		  HomingSituationReturn0.mCntlFuncPar.mDblArrParams[2] = valTCur + HomingSituationCur2.mCntlFuncPar.mDblArrParams[1];
		  HomingSituationReturn0.mCntlFuncPar.mDblArrParams[3] = valTCur + HomingSituationCur2.mCntlFuncPar.mDblArrParams[2];
		  HomingSituationReturn0.mShipTarg = HomingSituationCur2.mShipTarg ;
		  breturn = true ;
		}

	   }
   }
   return breturn ;
}
// функция нахождениея асимальной дальности захвата цели
// для заданной высоты бросания , заданного начального угла бросания Тетта0 и скорости бросания ищется  закон
// управления, обеспечивающие аксимальную дальность захвата цели ГСН
	// Закон управления  - кусочно-постоянный с тремя точками переключения
// valShagT - шаг по вреени при переборе точек переключения
// результаты сохраняются в HomingSituationReturn
 bool THomingSituation::fnc_FindMaxDist_For_Takeover_For_FixedTetta0_Perebor_3Point(const long double valShagT
   , THomingSituation &HomingSituationReturn0, long double &valMaxDistTakeover)
{
	return true;
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////



// ФУНКЦИИ НАХОЖДЕНИЯ МАКСИМАЛЬНОЙ ПРИЦЕЛЬНОЙ ДАЛЬНОСТИ С УЧЕТОМ УСЛОВИЙ ЗАХВАТА
// нахождение оптимаального управления для задачи по дальности методом  перебора
// для двух точек переключения
// INPUT:
// valShagT - шаг сетки  по времени
TBombTraj THomingSituation::fnc_Check_Possibility_Takeover_Perebor_2Point(const long double valShagT, long double& valDGorizOpt)
{
  /* TBombTraj BombTraj ( mTet0, mV0, mAltit, mX0, mStepInt, miControlType, marrTimes, 4, mBomb) ;

   BombTraj.marrTimes[0] = 0.;
   BombTraj.marrTimes[1] = 0;
   BombTraj.marrTimes[2] = 0;
   BombTraj.marrTimes[3] = 1000.;
   TBombTraj BombReturn = BombTraj ;
   long double valDGoroz = 0.;
   BombTraj.fncMoveClass_TO_ZeroAlt_AND_ShowGraphs(NULL, valDGoroz) ;
   long double valT = BombTraj.mTCur;
   long double valTOpt = 0.;
   valDGorizOpt = valDGoroz;
   long double valTCur = 0., valDGorizCur = 0.;
   int nCircle = valT / valShagT;


   for (int  i=0; i < nCircle; i++)
   {
	TBombTraj BombTrajCur ( mTet0, mV0, mAltit, mX0, mStepInt, miControlType, marrTimes, 4, mBomb) ;
	valTCur = ((long double)i) *valShagT ;
	BombTrajCur.marrTimes[0] = 0.;
	BombTrajCur.marrTimes[1] = valTCur;
	BombTrajCur.marrTimes[2] = 1000.;

	   BombTrajCur.fncMovePhasVector(BombTrajCur.marrTimes[1]);
	   long double arrTimes[20] = {0.} ;


	   TBombTraj  BombTrajCur1(BombTrajCur.marrStrSK_VS [4],BombTrajCur.marrStrSK_VS [3], BombTrajCur.marrStrSK_VS [1]
		 ,BombTrajCur.marrStrSK_VS [0], mStepInt, miControlType, marrTimes, mQuantTimes, BombTrajCur.mBomb) ;
	   TBombTraj TrajCur2 = BombTrajCur1.fncFindOptimalControl_for_Dist_MethodPerebora_1Point(valShagT, valDGorizCur)  ;
	   if (valDGorizCur > valDGorizOpt)
	   {
		  BombReturn.marrTimes[1] =  BombTrajCur.marrTimes[1];
		  BombReturn.marrTimes[2] =  BombTrajCur.marrTimes[1] +  TrajCur2.marrTimes[1];
		  BombReturn.marrTimes[3] =  BombTrajCur.marrTimes[1] +  TrajCur2.marrTimes[2];
		  valDGorizOpt = valDGorizCur ;
	   }
	 //  else
	//   {
	 //	   break;
	 //  }

   }

  return BombReturn ;  */
  return TBombTraj();
}

// ФУНКЦИИ НАХОЖДЕНИЯ МАКСИМАЛЬНОЙ ПРИЦЕЛЬНОЙ ДАЛЬНОСТИ С УЧЕТОМ УСЛОВИЙ ЗАХВАТА


// нахождение оптимаального управления для задачи НАХОЖДЕНИЯ МАКСИМАЛЬНОЙ ПРИЦЕЛЬНОЙ ДАЛЬНОСТИ С УЧЕТОМ УСЛОВИЙ ЗАХВАТА
//  методом  перебора
// для одной  точкит  переключения
// INPUT:
// valShagT - шаг сетки  по времени перебора момента переключения
// valFunnelAng  - угол воронки скорости УАБ для самонаведения
//  valDeltaAimPoint - точность по точке прицеливания
bool THomingSituation::fnc_Check_Possibility_Takeover_For_FixedTraj_Perebor_1Point(const long double valDeltaAimPoint
	 ,const long double valFunnelAng,const long double valShagT, THomingSituation &HomingSituationReturn)
{
   THomingSituation HomingSituation  = *this ;
   TCntlFuncPar CntlFuncParTemp =  HomingSituation.mCntlFuncPar ;
   long double *parrTimes = CntlFuncParTemp.mDblArrParams;
   parrTimes[0] = 0.;
   parrTimes[1] = 0.;
   parrTimes[2] = 1000.;
   CntlFuncParTemp.miArrParams[0] = 3; // к-во точек переключения
   HomingSituationReturn  = *this ;
   long double valDGoroz = 0.;
   HomingSituation.fncMoveBomb_TO_ZeroAlt_AND_ShowGraphs(NULL, valDGoroz) ;

   long double valT = HomingSituation.mBombTraj.mTCur;
   long double valTOpt = 0.;

   long double valTCur = 0.;
   int nCircle = valT / valShagT;

   for (int  i=0; i < nCircle; i++)
   {
	THomingSituation HomingSituationCur = *this;// ( mTet0, mV0, mAltit, mX0, mStepInt, miControlType, marrTimes, 3, mBomb) ;
	HomingSituationCur.mCntlFuncPar.miArrParams[0] = 3 ;
	valTCur = ((long double)i) *valShagT ;
	HomingSituationCur.mCntlFuncPar.mDblArrParams[0] = 0.;
	HomingSituationCur.mCntlFuncPar.mDblArrParams[1] =  valTCur;
	HomingSituationCur.mCntlFuncPar.mDblArrParams[2] = 1000.;

   if(HomingSituationCur.fncMoveBomb_TO_TakeoverPoint_AND_ShowGraphs(NULL,valDeltaAimPoint, valFunnelAng))
   {
	 HomingSituationReturn =  HomingSituationCur ;
	 return true;
   }

 }

  return false;
}

//пернеещение АБ на расстояние захвата
//ЕСли невозможно, то возвращает false
// Если возможно, то true
//
bool THomingSituation::fncMoveBomb_TO_TakeoverPoint_AND_ShowGraphs(wchar_t *wcharrPath
	  ,const long double valDeltaAimPoint,const long double valFunnelAng)
{
	bool breturn = false;
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
	wcscpy( &pwcharrFileNames[ 6 * 30], L"AlfAt");
	wcscpy( &pwcharrFileNames[ 7 * 30], L"M");
	wcscpy( &pwcharrFileNames[ 8 * 30], L"Cy");
	wcscpy( &pwcharrFileNames[ 9 * 30], L"CyqSmUmaxMG");

	pscaleY = new double  [QUANT_COLS] ;
	pscaleY[1] = 0.1;
	pscaleY[2] = 0.1;

	pscaleY[3] = 1000.;
	pscaleY[4] = 1.;
	pscaleY[5] = 1000.;
	pscaleY[6] = 1000.;
	pscaleY[7] = 1000.;
	pscaleY[8] = 1000.;
	pscaleY[9] = 1000.;
}

  int iCirc = 1000. / mBombTraj.mStepInt ;
 // long double valTTemp = mTCur + ((long double)iCirc) * mStepInt ;
  int i = 0 ;
  int iNupPointsOut = 0;
  double valTOut = -DEL_T ;
  for ( i = 0; i < iCirc; i++)
  {

	//  if (mTCur > 9.3099999)
	  //{
		// int iii = 0;
	//  }
	  long double valU = fncCalcU( mTCur)  ;

	  // шаг интегрования
	  mBombTraj.fncEilerStep(valU, mBombTraj.mStepInt);
	  mTCur = mBombTraj.mTCur;
	  ///


	  if  ( iNupPointsOut < QUANT_POINTS_MAX )
	  {

	  valTOut = (double)mTCur ;

	   if (wcharrPath)
	   {
	   double *p = &parrBuff[iNupPointsOut *  QUANT_COLS] ;
	   p[0] =  (double)mTCur ;
	   p[1] = (double)mBombTraj.marrStrSK_VS [0];
	   p[2] = (double)mBombTraj.marrStrSK_VS [1];
	   p[3] = (double)mBombTraj.marrStrSK_VS [2];
	   p[4] = (double)mBombTraj.marrStrSK_VS [3];
	   p[5] = (double)mBombTraj.marrStrSK_VS [4];
	   p[6] = (double)fncCalcU( mTCur);
	   long double valTay,  valDerivTay =0., valMach = -1., arrGradMach[5] = {0.};

	   fncCalcNormTemperature(p[2], valTay, valDerivTay) ;
	   mBombTraj.fncCalcMach_and_GradMach( valTay,  valDerivTay
	 ,valMach, arrGradMach);
	   p[7] = (double)valMach;
	   long double  val_q = 0., arrGrad_q [5] = {0.};
	   mBombTraj.fncCalc_q_and_Grad_q(val_q, arrGrad_q) ;
	   p[8] = mBombTraj.mBomb.fncCy(valMach);
	   p[9] = p[8] * val_q * mBombTraj.mBomb.mMidSq * mBombTraj.mBomb.mMaxAt / mBombTraj.mBomb.mMass / G_ZEMLI ;

		}

	  iNupPointsOut++;

	  }

	 //  проверка условий захвата
	  long double arrDelta[2] ={0.};
	  if (mBombTraj.marrStrSK_VS[1] <0.) return false;


	  MtrxMinusMatrx(mBombTraj.marrStrSK_VS, mShipTarg.mTraject.marrVectSostGSK_Begin,1, 2, arrDelta);
	  long double delta =  NormVect2(arrDelta);
	  if ( delta > mBombTraj.mBomb.mHomingHead.mR)
	  {
	   continue;
	  }
	  else
	  {  if(
			 ( (M_PI / 2. + mBombTraj.marrStrSK_VS[4]) <  valFunnelAng)
		  && ( fabsl(mBombTraj.marrStrSK_VS[0] - mBombTraj.marrStrSK_VS[1]/ tanl(mBombTraj.marrStrSK_VS[4] )
			   - mShipTarg.mTraject.marrVectSostGSK_Begin[0] ) < valDeltaAimPoint)
		  )
		 {
			breturn = true;
		 }
		 else
		 {
		  breturn = false;
		 }
		 break;
	  }

  }
 if (!breturn) return false;


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

		wchar_t wchFileName4[300] = {0} ;
		wcscpy(wchFileName4, wcharrPath );
		wcscat(wchFileName4, L"\\Axes.shp");
		TYrWriteShapeFile::CreateShpAxes(wchFileName4,-80000.,80000,-80000.,80000) ; ;

		delete parrBuff ;
		delete pwcharrFileNames ;
		delete pscaleY ;
   }

	return true ;
}


// нахождение угла бросания для захвата на определенной дальности
//  методом  перебора  для одной  точкит  переключения
//
// INPUT:
// valShagT - шаг сетки  по времени перебора момента переключения
// valFunnelAng  - угол воронки скорости УАБ для самонаведения
//  valDeltaAimPoint - точность по точке прицеливания
bool THomingSituation::fnc_FindPossibleAngle_For_Takeover_Perebor(const int iQuantPointsOverswitch
   , const long double valDeltaAimPoint ,const long double valFunnelAng,const long double valShagT
   , THomingSituation &HomingSituationReturn)
{
	long double valShagAng = 2.5/ 180. * M_PI ;
	const int iNum = 30;
	for (int i = 0; i < 2 * iNum + 1; i++)
	{
	  long double valAng = ((long double)iNum) * valShagAng -   valShagAng * ((long double)i );
	  THomingSituation  HomingSituationCur = *this;
	  TBombTraj BombTraj ( valAng ,mBombTraj.mV0 ,mBombTraj.mAltit  ,mBombTraj.mX0 	,mBombTraj.mStepInt  ,mBombTraj.mTStart ,mBombTraj.mBomb );
	HomingSituationCur.mBombTraj = BombTraj;
	HomingSituationCur.mCntlFuncPar.miControlType = 1;
	THomingSituation HomingSituationCur1 = HomingSituationCur;
	HomingSituationCur1.mCntlFuncPar.miControlType = 2;
	bool brez = false,brez1 = false;

	switch( iQuantPointsOverswitch)
	{
	case 1:
	brez = HomingSituationCur.fnc_Check_Possibility_Takeover_For_FixedTraj_Perebor_1Point( valDeltaAimPoint
	 ,valFunnelAng, valShagT, HomingSituationReturn);
	 if (brez)
	 {
	   return true;
	 }
	  brez1 = HomingSituationCur1.fnc_Check_Possibility_Takeover_For_FixedTraj_Perebor_1Point( valDeltaAimPoint
	 ,valFunnelAng, valShagT, HomingSituationReturn);

	if (brez1)
	 {
	   return true;
	 }
	break;

	case 2:
	 brez = HomingSituationCur.fnc_Check_Possibility_Takeover_For_FixedTraj_Perebor_2Point( valDeltaAimPoint
	 ,valFunnelAng, valShagT, HomingSituationReturn);
	 if (brez)
	 {
	   return true;
	 }
	  brez1 = HomingSituationCur1.fnc_Check_Possibility_Takeover_For_FixedTraj_Perebor_2Point( valDeltaAimPoint
	 ,valFunnelAng, valShagT, HomingSituationReturn);

	if (brez1)
	 {
	   return true;
	 }
	break;

	case 3:
	brez = HomingSituationCur.fnc_Check_Possibility_Takeover_For_FixedTraj_Perebor_3Point( valDeltaAimPoint
	 ,valFunnelAng, valShagT, HomingSituationReturn);
	 if (brez)
	 {
	   return true;
	 }
	  brez1 = HomingSituationCur1.fnc_Check_Possibility_Takeover_For_FixedTraj_Perebor_3Point( valDeltaAimPoint
	 ,valFunnelAng, valShagT, HomingSituationReturn);

	if (brez1)
	 {
	   return true;
	 }
	break;

	 default:

	return false;
	}


   }
	return false ;
}


// интгрирование уравнений до момента  падения
// OUTPUT:
// valFixedTime - конец интеравала интегрирования
void THomingSituation::fncMoveClass_TO_FixedTime_AND_ShowGraphs( wchar_t *wcharrPath, const long double valFixedTime )
{
   int iCirc = (valFixedTime - mTCur) / mBombTraj.mStepInt ;
 	bool breturn = false;
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
	wcscpy( &pwcharrFileNames[ 6 * 30], L"AlfAt");
	wcscpy( &pwcharrFileNames[ 7 * 30], L"M");
	wcscpy( &pwcharrFileNames[ 8 * 30], L"Cy");
	wcscpy( &pwcharrFileNames[ 9 * 30], L"CyqSmUmaxMG");

	pscaleY = new double  [QUANT_COLS] ;
	pscaleY[1] = 0.1;
	pscaleY[2] = 0.1;

	pscaleY[3] = 1000.;
	pscaleY[4] = 1.;
	pscaleY[5] = 1000.;
	pscaleY[6] = 1000.;
	pscaleY[7] = 1000.;
	pscaleY[8] = 1000.;
	pscaleY[9] = 1000.;
}

  int i = 0 ;
  int iNupPointsOut = 0;
  double valTOut = -DEL_T ;
  for ( i = 0; i < iCirc; i++)
  {


	  long double valU = fncCalcU( mTCur)  ;

	  // шаг интегрования
	  mBombTraj.fncEilerStep(valU, mBombTraj.mStepInt);
	  mTCur = mBombTraj.mTCur;
	  ///


	  if  ( iNupPointsOut < QUANT_POINTS_MAX )
	  {

	  valTOut = (double)mTCur ;

	   if (wcharrPath)
	   {
	   double *p = &parrBuff[iNupPointsOut *  QUANT_COLS] ;
	   p[0] =  (double)mTCur ;
	   p[1] = (double)mBombTraj.marrStrSK_VS [0];
	   p[2] = (double)mBombTraj.marrStrSK_VS [1];
	   p[3] = (double)mBombTraj.marrStrSK_VS [2];
	   p[4] = (double)mBombTraj.marrStrSK_VS [3];
	   p[5] = (double)mBombTraj.marrStrSK_VS [4];
	   p[6] = (double)fncCalcU( mTCur);
	   long double valTay,  valDerivTay =0., valMach = -1., arrGradMach[5] = {0.};

	   fncCalcNormTemperature(p[2], valTay, valDerivTay) ;
	   mBombTraj.fncCalcMach_and_GradMach( valTay,  valDerivTay
	 ,valMach, arrGradMach);
	   p[7] = (double)valMach;
	   long double  val_q = 0., arrGrad_q [5] = {0.};
	   mBombTraj.fncCalc_q_and_Grad_q(val_q, arrGrad_q) ;
	   p[8] = mBombTraj.mBomb.fncCy(valMach);
	   p[9] = p[8] * val_q * mBombTraj.mBomb.mMidSq * mBombTraj.mBomb.mMaxAt / mBombTraj.mBomb.mMass / G_ZEMLI ;

		}

	  iNupPointsOut++;

	  }

  }



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

		wchar_t wchFileName4[300] = {0} ;
		wcscpy(wchFileName4, wcharrPath );
		wcscat(wchFileName4, L"\\Axes.shp");
		TYrWriteShapeFile::CreateShpAxes(wchFileName4,-80000.,80000,-80000.,80000) ; ;

		// отрезок достижимости
		   TURPointXY arrOneLine[2];
		   arrOneLine[0] =  TURPointXY(0., pscaleY[7]);
		   arrOneLine[1] =  TURPointXY(40000. , pscaleY[7]);
		   TURPolyLine plnOneLine(arrOneLine, 2) ;
		   wcscpy(wchFileName4, wcharrPath );
		  wcscat(wchFileName4, L"\\plnOneLine.shp");
		  plnOneLine.WriteSetSHPFiles(wchFileName4, &plnOneLine, 1 );
		  ///


		delete parrBuff ;
		delete pwcharrFileNames ;

		// построение графиков для цели
		   // положение цели на момент пеленга
		TURPointXY pntTargPos(mShipTarg.mTraject.marrVectSostGSK_Begin[0], 0.);
		wcscpy(wchFileName4, wcharrPath );
		wcscat(wchFileName4, L"\\pntTargPos.shp");
		pntTargPos.WriteSetSHPFiles(wchFileName4, &pntTargPos, 1 );
		   ///

		   // отрезок достижимости
		   TURPointXY arrTargPos[2];
		   arrTargPos[0] =  TURPointXY(mShipTarg.mTraject.marrVectSostGSK_Begin[0] - 20.* mBombTraj.mTCur, 0.);
		   arrTargPos[1] =  TURPointXY(mShipTarg.mTraject.marrVectSostGSK_Begin[0] + 20.* mBombTraj.mTCur, 0.);
		   TURPolyLine plnArea(arrTargPos,2) ;
		   wcscpy(wchFileName4, wcharrPath );
		  wcscat(wchFileName4, L"\\plnTargArea.shp");
		  plnArea.WriteSetSHPFiles(wchFileName4, &plnArea, 1 );
		  ///

		  // построение графиков для УАБ
		   // график вектора скорости
			 TURPointXY arrVeloBomb[2];
			 double valScale = 1.;
		   arrVeloBomb[0] =  TURPointXY(mBombTraj.marrStrSK_VS[0] ,mBombTraj.marrStrSK_VS[1]);
		   arrVeloBomb[1] =  TURPointXY(mBombTraj.marrStrSK_VS[0] + valScale * mBombTraj.marrStrSK_VS[3] * cos (mBombTraj.marrStrSK_VS [4])
				 , mBombTraj.marrStrSK_VS[1] + valScale * mBombTraj.marrStrSK_VS[3] * sin (mBombTraj.marrStrSK_VS [4]));
		   TURPolyLine plnVeloBomb(arrVeloBomb,2) ;
		   wcscpy(wchFileName4, wcharrPath );
		  wcscat(wchFileName4, L"\\plnVeloBomb.shp");
		  plnVeloBomb.WriteSetSHPFiles(wchFileName4, &plnVeloBomb, 1 );
		  ///


		  // построение сектора захвата
		  const double valFi0 = mBombTraj.marrStrSK_VS [4] + mBombTraj.mBomb.mHomingHead.mFi /2.;
		  const double valFi1 = mBombTraj.marrStrSK_VS [4] - mBombTraj.mBomb.mHomingHead.mFi /2.;
		  TURPolygon plgZahv =  TURPolygon::fncCreateSector(arrVeloBomb[0], mBombTraj.mBomb.mHomingHead.mR, valFi0, valFi1,200) ;
		   wcscpy(wchFileName4, wcharrPath );
		  wcscat(wchFileName4, L"\\plgZahv.shp");
		  plgZahv.WriteSetSHPFiles(wchFileName4, &plgZahv, 1 );

	  delete pscaleY ;

   }



}





// нахождение оптимаального управления для задачи НАХОЖДЕНИЯ МАКСИМАЛЬНОЙ ПРИЦЕЛЬНОЙ ДАЛЬНОСТИ С УЧЕТОМ УСЛОВИЙ ЗАХВАТА
//  методом  перебора
// для 2 точек перключения
// INPUT:
// valShagT - шаг сетки  по времени перебора момента переключения
// valFunnelAng  - угол воронки скорости УАБ для самонаведения
//  valDeltaAimPoint - точность по точке прицеливания
bool THomingSituation::fnc_Check_Possibility_Takeover_For_FixedTraj_Perebor_2Point(const long double valDeltaAimPoint
	 ,const long double valFunnelAng,const long double valShagT, THomingSituation &HomingSituationReturn)
{
	TCntlFuncPar CntlFuncPar;
	CntlFuncPar.miControlType = mCntlFuncPar.miControlType ;
	CntlFuncPar.miArrParams[0] = 4;
	CntlFuncPar.mDblArrParams[0] = 0.;
	CntlFuncPar.mDblArrParams[1] = 0.;
	CntlFuncPar.mDblArrParams[2] = 0.;
	CntlFuncPar.mDblArrParams[3] = 1000.;

	THomingSituation HomingSituation0(mShipTarg, mBombTraj, CntlFuncPar,mTCur ,NULL);
	THomingSituation HomingSituation1 = HomingSituation0 ;
	long double valDGoriz = -1.;
	HomingSituation1.fncMoveBomb_TO_ZeroAlt_AND_ShowGraphs(NULL, valDGoriz);

   long double valT = HomingSituation1.mBombTraj.mTCur;
   long double valTOpt = 0.;

   long double valTCur = 0., valDGorizCur = 0.;
   int nCircle = valT / valShagT;

   for (int  i=0; i < nCircle; i++)
   {
	THomingSituation HomingSituationCur =  HomingSituation0;
	valTCur = ((long double)i) *valShagT ;
	HomingSituationCur.mCntlFuncPar.mDblArrParams[0] = 0.;
	HomingSituationCur.mCntlFuncPar.mDblArrParams[1] = valTCur;
	HomingSituationCur.mCntlFuncPar.mDblArrParams[0] = 1000.;
	HomingSituationCur.fncMoveClass_TO_FixedTime_AND_ShowGraphs( NULL, HomingSituationCur.mCntlFuncPar.mDblArrParams[1] ) ;


	TBombTraj BombTrajCur1(	mBombTraj.marrStrSK_VS [4] 	,mBombTraj.marrStrSK_VS [3] ,mBombTraj.marrStrSK_VS [1] ,mBombTraj.marrStrSK_VS [0]  // нач. положение по оси X
	,mBombTraj.mStepInt  ,0.,mBombTraj.mBomb );
	TCntlFuncPar CntlFuncParCur;
	CntlFuncParCur.miControlType = mCntlFuncPar.miControlType ;
	CntlFuncParCur.miArrParams[0] = 3;
	CntlFuncParCur.mDblArrParams[0] = 0.;
	CntlFuncParCur.mDblArrParams[1] = 0.;
	CntlFuncParCur.mDblArrParams[2] = 1000.;
	THomingSituation HomingSituationCur1(mShipTarg , BombTrajCur1, CntlFuncParCur,0.,NULL );
		THomingSituation HomingSituationCur2 =  HomingSituationCur1;
		HomingSituationReturn =   HomingSituationCur1;
	  bool brez1 = HomingSituationCur1.fnc_Check_Possibility_Takeover_For_FixedTraj_Perebor_1Point
		(valDeltaAimPoint, valFunnelAng, valShagT, HomingSituationCur2);
		if (brez1)
		{
		  HomingSituationReturn =  HomingSituation0 ;
		  HomingSituationReturn.mCntlFuncPar.miArrParams[0] = 4;
		  HomingSituationReturn.mCntlFuncPar.mDblArrParams[0] = 0.;
		  HomingSituationReturn.mCntlFuncPar.mDblArrParams[1] = valTCur ;
		  HomingSituationReturn.mCntlFuncPar.mDblArrParams[2] = valTCur + HomingSituationCur2.mCntlFuncPar.mDblArrParams[1];
		  HomingSituationReturn.mCntlFuncPar.mDblArrParams[3] = valTCur + HomingSituationCur2.mBombTraj.mTCur;
		  return true ;
		}

   }

  return false ;
}


// нахождение оптимаального управления для задачи НАХОЖДЕНИЯ МАКСИМАЛЬНОЙ ПРИЦЕЛЬНОЙ ДАЛЬНОСТИ С УЧЕТОМ УСЛОВИЙ ЗАХВАТА
//  методом  перебора
// для 3 точек перключения
// INPUT:
// valShagT - шаг сетки  по времени перебора момента переключения
// valFunnelAng  - угол воронки скорости УАБ для самонаведения
//  valDeltaAimPoint - точность по точке прицеливания
bool THomingSituation::fnc_Check_Possibility_Takeover_For_FixedTraj_Perebor_3Point(const long double valDeltaAimPoint
	 ,const long double valFunnelAng,const long double valShagT, THomingSituation &HomingSituationReturn)
{
	TCntlFuncPar CntlFuncPar;
	CntlFuncPar.miControlType = mCntlFuncPar.miControlType ;
	CntlFuncPar.miArrParams[0] = 5;
	CntlFuncPar.mDblArrParams[0] = 0.;
	CntlFuncPar.mDblArrParams[1] = 0.;
	CntlFuncPar.mDblArrParams[2] = 0.;
	CntlFuncPar.mDblArrParams[3] = 0.;
	CntlFuncPar.mDblArrParams[4] = 1000.;

	THomingSituation HomingSituation0(mShipTarg, mBombTraj, CntlFuncPar,mTCur ,NULL);
	THomingSituation HomingSituation1 = HomingSituation0 ;
	long double valDGoriz = -1.;
	HomingSituation1.fncMoveBomb_TO_ZeroAlt_AND_ShowGraphs(NULL, valDGoriz);

   long double valT = HomingSituation1.mBombTraj.mTCur;
   long double valTOpt = 0.;

   long double valTCur = 0., valDGorizCur = 0.;
   int nCircle = valT / valShagT;

   for (int  i=0; i < nCircle; i++)
   {
	THomingSituation HomingSituationCur =  HomingSituation0;
	valTCur = ((long double)i) *valShagT ;
	HomingSituationCur.mCntlFuncPar.mDblArrParams[0] = 0.;
	HomingSituationCur.mCntlFuncPar.mDblArrParams[1] = valTCur;
	HomingSituationCur.mCntlFuncPar.mDblArrParams[2] = 0.;
	HomingSituationCur.mCntlFuncPar.mDblArrParams[3] = 1000.;
	HomingSituationCur.fncMoveClass_TO_FixedTime_AND_ShowGraphs( NULL, HomingSituationCur.mCntlFuncPar.mDblArrParams[1] ) ;

	TBombTraj BombTrajCur1(	mBombTraj.marrStrSK_VS [4]	,mBombTraj.marrStrSK_VS [3]  ,mBombTraj.marrStrSK_VS [1]  ,mBombTraj.marrStrSK_VS [0] ,mBombTraj.mStepInt   ,0. ,mBombTraj.mBomb );
	TCntlFuncPar CntlFuncParCur;
	CntlFuncParCur.miControlType = mCntlFuncPar.miControlType ;
	CntlFuncParCur.miArrParams[0] = 4;
	CntlFuncParCur.mDblArrParams[0] = 0.;
	CntlFuncParCur.mDblArrParams[1] = 0.;
	CntlFuncParCur.mDblArrParams[2] = 1000.;
	THomingSituation HomingSituationCur1(mShipTarg //
		, BombTrajCur1 //
		, CntlFuncParCur
		,0.
		,NULL // путь у папке с отчетом
		);
		THomingSituation HomingSituationCur2 =  HomingSituationCur1;
		HomingSituationReturn =   HomingSituationCur1;
	  bool brez1 = HomingSituationCur1.fnc_Check_Possibility_Takeover_For_FixedTraj_Perebor_2Point
		(valDeltaAimPoint, valFunnelAng, valShagT, HomingSituationCur2);
		if (brez1)
		{
		  HomingSituationReturn =  HomingSituation0 ;
		  HomingSituationReturn.mCntlFuncPar.miArrParams[0] = 4;
		  HomingSituationReturn.mCntlFuncPar.mDblArrParams[0] = 0.;
		  HomingSituationReturn.mCntlFuncPar.mDblArrParams[1] = valTCur ;
		  HomingSituationReturn.mCntlFuncPar.mDblArrParams[2] = valTCur + HomingSituationCur2.mCntlFuncPar.mDblArrParams[1];
		  HomingSituationReturn.mCntlFuncPar.mDblArrParams[3] = valTCur + HomingSituationCur2.mCntlFuncPar.mDblArrParams[2];
		  HomingSituationReturn.mCntlFuncPar.mDblArrParams[4] = valTCur + HomingSituationCur2.mBombTraj.mTCur;
		  return true ;
		}

   }

  return false ;
}


long double THomingSituation::calcPsi1()
{
  long double x1 =  mShipTarg.mTraject.marrVectSostGSK[0]  - mBombTraj.marrStrSK_VS [0] + mShipTarg.mTargData. mLDanger /2. ;
  long double y1 =  -mBombTraj.marrStrSK_VS [1] ;
  long double x2 =  mShipTarg.mTraject.marrVectSostGSK[0]  - mBombTraj.marrStrSK_VS [0]  ;
  long double y2 = -mBombTraj.marrStrSK_VS [1] ;
  long double t1 = acosl((x1 * x2 + y1 * y2)/ sqrtl(x1 * x1 + y1 * y1)/ sqrtl(x2 * x2 + y2 * y2));
  return fabsl ( t1);
}
long double THomingSituation::calcPsi2()
{
  long double x1 =  mShipTarg.mTraject.marrVectSostGSK[0]  - mBombTraj.marrStrSK_VS [0] - mShipTarg.mTargData. mLDanger /2. ;
  long double y1 =  -mBombTraj.marrStrSK_VS [1] ;
  long double x2 =  mShipTarg.mTraject.marrVectSostGSK[0]  - mBombTraj.marrStrSK_VS [0]  ;
  long double y2 = -mBombTraj.marrStrSK_VS [1] ;
  return fabsl ( acosl((x1 * x2 + y1 * y2)/ sqrtl(x1 * x1 + y1 * y1)/ sqrtl(x2 * x2 + y2 * y2)));

}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

int THomingSituation::fncSign (long double a)
{
	return (a >=0.)? 1:-1;
}

#pragma package(smart_init)
