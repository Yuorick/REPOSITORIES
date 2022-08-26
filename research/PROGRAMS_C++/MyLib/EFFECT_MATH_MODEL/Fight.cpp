﻿//---------------------------------------------------------------------------


	#pragma hdrstop
	#include <string.h>
	#include <stdlib.h>
	#include <math.h>
	#include <vcl.h>
	#include <dir.h>

	#include "Fight.h"
	#include "CalcCorMatrx.h"
	#include "StabSyst2.h"
	#include "MatrixProccess.h"
	#include "YrWriteShapeFile.h"
	#include "Zamer.h"
	#include "MyShellTraj.h"
	#include "Equations.h"
	#include "Comp.h"
	#include "NeighbourhoodAppPoint.h"
	#include "ProbabilityTheory.h"
	#include "YrRastr.h"
	#include "Gauss.h"
	#include "CoastTargNeibourhood.h"
	#include "OutPutAeroShot.h"

   	#define COEFF_BOEGOTOVNOSTY_SNARIADA 0.9
  //	#define COEFF_BOEGOTOVNOSTY_SNARIADA 1.


 extern const double NODATA;
 #define VALKMMO 0. // часть ошибки измерения положения антенны, которая уходит в ММО
 #define SIGW 0.1//33.0    // скз шума в движении цели, исползуемое в алгоритме фильтрации сопровождения цели
 const int QUANT_COLS_BUFF = 67;
 // тип фыильтра. если = 1 то оптимальный, если = 0, то штатный КФ, если -1 то КФ с пересчетом
 const int  TYPE_OF_FILTER = 1;
//---------------------------------------------------------------------------
// размерность вектора разбросов параметров движения снаряда

TFight::TFight()
{

	// корабль
	mVessel =TVessel();
	//  цель
	mTarget = TTarget();
	// фильтр сопровождения
	mTraceFlt = TTraceFlt();
	// эталонный сигнал
	TEtalonSign mEtalonSign = TEtalonSign();
	//  внешняя среда (веиер, качка)
	TEnvironment mEnvironment =TEnvironment();

	// Искусственный промах
   //	mMissSimulated =0.;

	// Искусственный Nedolet
   //	mPereletSimulated =0.;

	// текущее время
	double mT;
	mpwcharrFoldReport = NULL ;
	mparrBuff   = NULL ;
	mQuantPntReport = 0;
	mLenMemoryAlloc = 0 ;

}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
 TFight::~TFight()
{
   if (mparrBuff != NULL)
   {
   free( mparrBuff)  ;
   mparrBuff = NULL ;
   }
   if ( mpwcharrFoldReport != NULL)
   {
   free( mpwcharrFoldReport );
   mpwcharrFoldReport = NULL ;
   }


}


// конструктор копирования
 TFight ::TFight (const TFight &R)
 {
	mVessel  =  R.mVessel ;
	mTarget = R.mTarget;
	mT = R.mT;
	mTraceFlt = R.mTraceFlt;
	mEnvironment = R.mEnvironment;
	mTExtrap = R.mTExtrap;
	mTargDesBet = R.mTargDesBet;
	mTargDesEps = R. mTargDesEps;
	mTargDesR = R.mTargDesR ;
	mSigTargDesBet = R.mSigTargDesBet;
	mSigTargDesEps =R.mSigTargDesEps;
	mSigTargDesR = R.mSigTargDesR ;
	memcpy( marrVSExtrap_KGSK, R.marrVSExtrap_KGSK, 6 * sizeof(double));
	memcpy( marrSigExtrapPosit_KGSK, marrSigExtrapPosit_KGSK, 3 * sizeof(double));
	mSigModulV_KGSK = mSigModulV_KGSK ;
	mEtalonSign = R.mEtalonSign;
	mAntCoeff = R.mAntCoeff;
   //	mMissSimulated = R.mMissSimulated;
   //	mPereletSimulated = R.mPereletSimulated;

			// для отчета
		mLenMemoryAlloc = R.mLenMemoryAlloc ;
		mQuantPntReport = R.mQuantPntReport ;

		mparrBuff  = NULL;

		if(R.mparrBuff  != NULL)
		{
		if((mparrBuff  = (double*)malloc( QUANT_COLS_BUFF * R.mLenMemoryAlloc * sizeof(double))) != NULL)
		{
		memcpy( mparrBuff ,R.mparrBuff , QUANT_COLS_BUFF * R.mLenMemoryAlloc * sizeof(double));
		}
		else
		{
		ShowMessage(L"Not memory for mparrBuff ") ;
		Abort() ;
		}

		}

		//
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
 TFight &TFight::operator=(const TFight  &R)
 {
	mVessel  =  R.mVessel ;
	mTarget = R.mTarget;
	mT = R.mT;
	mTraceFlt = R.mTraceFlt;
	mEnvironment = R.mEnvironment;
	mTExtrap = R.mTExtrap;
	mTargDesBet = R.mTargDesBet;
	mTargDesEps = R. mTargDesEps;
	mTargDesR = R.mTargDesR ;
	mSigTargDesBet = R.mSigTargDesBet;
	mSigTargDesEps =R.mSigTargDesEps;
	mSigTargDesR = R.mSigTargDesR ;
	memcpy( marrVSExtrap_KGSK, R.marrVSExtrap_KGSK, 6 * sizeof(double));
	memcpy( marrSigExtrapPosit_KGSK, marrSigExtrapPosit_KGSK, 3 * sizeof(double));
	mSigModulV_KGSK = mSigModulV_KGSK ;
	mEtalonSign = R.mEtalonSign;
	mAntCoeff = R.mAntCoeff;
  //	mMissSimulated = R.mMissSimulated;
   //	mPereletSimulated = R.mPereletSimulated;



			// для отчета
		mLenMemoryAlloc = R.mLenMemoryAlloc ;
		mQuantPntReport = R.mQuantPntReport ;

		mparrBuff  = NULL;

		if(R.mparrBuff  != NULL)
		{
		if((mparrBuff  = (double*)malloc( QUANT_COLS_BUFF * R.mLenMemoryAlloc * sizeof(double))) != NULL)
		{
		memcpy( mparrBuff ,R.mparrBuff , QUANT_COLS_BUFF * R.mLenMemoryAlloc * sizeof(double));
		}
		else
		{
		ShowMessage(L"Not memory for mparrBuff ") ;
		Abort() ;
		}

		}

		//
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
/*
	// парам конструктор  1
 TFight::TFight (const TVessel Vessel, const TTarget Target
		,const double TExtrap, wchar_t *pwcharrFoldReport)
 {
		mEtalonSign = TEtalonSign();

		mVessel  = Vessel ;
		mTarget = Target ;
		if( mTarget.mpwcharrFoldReport != NULL)
		{
		mTarget.mQuantPntReport = 0 ;
		}
		if( mVessel.mpwcharrFoldReport != NULL)
		{
		mVessel.mQuantPntReport = 0 ;

		}
		// инициализация фильтра
		//  1.Время привязки оценок
		mTraceFlt.mTf = mTarget.mTraject.mTCur ;

		// 2.Оенка вектора состояния цели в КГСК на момент mTarget.mT
		double arrTemp [ 6 ]= {0.} ;
		memcpy(arrTemp , mTarget.mTraject.marrVectSostGSK, 6 * sizeof(double)) ;
		MtrxMinusMatrx(arrTemp ,  mVessel.marrVectSost,3, 1, mTraceFlt.marrVS_KGSK);

		// 3.Оценка вектора состояния в АСК на момент mTf
		RecalcVect_KGSK_INTO_ASK_True (mTraceFlt.marrVS_KGSK,mTraceFlt.marrVS_ASK,6 ) ;


		//  4.Объединенная корреляционная матрица ошибок оценивания в КГСК

		memset(mTraceFlt.marrK_KGSK, 0, 36 * sizeof(double)) ;
		for (int i = 0; i < 3; i++)
		{
		mTraceFlt.marrK_KGSK[ 6 * i + i ] = 1000000. ;
		mTraceFlt.marrK_KGSK[ 6 * i + i + 3] = 100000. ;
		mTraceFlt.marrK_KGSK[ 6 * (i + 3) + i] = 100000. ;
		mTraceFlt.marrK_KGSK[ 6 * (i + 3) + i + 3] = 100000. ;

		}

		//	5.Объединенная корреляционная матрица ошибок оценивания в ASK
		memcpy( mTraceFlt.marrK_ASK,mTraceFlt.marrK_KGSK,36 * sizeof(double));

		// 6. вектор углов
		mTraceFlt.marrMu[0] = mVessel.mSins.mEstQ;
		mTraceFlt.marrMu[1] = mVessel.mSins.mEstPsi ;
		mTraceFlt.marrMu[2] = mVessel.mSins.mEstTet ;
		mTraceFlt.marrMu[3] = mVessel.mDriver.mRealBet ;
		mTraceFlt.marrMu[4] = mVessel.mDriver.mRealEps ;

		//  7.матрица перехода из АСК в ПСК на момент mTf

		calcMatr_ASK_v_PSK( &mTraceFlt.marrMu[3],mTraceFlt.marrF_ASK_V_PSK) ;


		// 8.расширенная матрица перехода из АСК в ПСК на момент mTf
		createExtendMtrx(mTraceFlt.marrF_ASK_V_PSK, mTraceFlt.marrFExt_ASK_V_PSK) ;

		//  9.матрица перехода из ПСК в КГСК на момент mTf
		calcMatr_PSK_v_KGSK(mTraceFlt.marrMu, mTraceFlt.marrF_PSK_V_KGSK) ;


		// 10. расширенная матрица перехода из  ПСК в КГСК мна момент mTf
		createExtendMtrx(mTraceFlt.marrF_PSK_V_KGSK, mTraceFlt.marrFExt_PSK_V_KGSK) ;

		// 11. матрица перехода из  АСК в КГСК мна момент mTf
		MtrxMultMatrx(mTraceFlt.marrF_PSK_V_KGSK,3, 3, mTraceFlt.marrF_ASK_V_PSK,3, mTraceFlt.marrF_ASK_V_KGSK) ;

		// 12. расширенная матрица перехода из  АСК в КГСК мна момент mTf
		createExtendMtrx(mTraceFlt.marrF_ASK_V_KGSK, mTraceFlt.marrFExt_ASK_V_KGSK) ;
		// 11. признак того, что фильтр прошел инициализацию
		mTraceFlt.mbInit = true ;

		// ИНИЦИАЛИЗАЦИЯ ПЕРВИЧНЫХ ЦЕЛЕУКЫАЗАНИЙ
		mT =  mTarget.mTraject.mTCur ;
		mTExtrap = TExtrap;
		// дяля  ычисление истинного положения цели в ПСК - АС
		double  arrVSExtrap_PSK_AS [6] = {0};

// расчет первичных целеуказаний (на основе внешних данных о положении цели)
	getTargDes_PSK_True( TExtrap, mTargDesBet, mTargDesEps, mTargDesR, arrVSExtrap_PSK_AS ) ;




	// для формироывания буфера памяти для отчета
	mpwcharrFoldReport = NULL ;
	mparrBuff   = NULL ;
	mQuantPntReport = 0 ;
	mLenMemoryAlloc = 0 ;
	if(pwcharrFoldReport != NULL)
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

	mLenMemoryAlloc = 2000;
	if((mparrBuff  = (double*)malloc( QUANT_COLS_BUFF * mLenMemoryAlloc * sizeof(double))) != NULL)
	{

	memset( mparrBuff,  0,  QUANT_COLS_BUFF * mLenMemoryAlloc * sizeof(double));
	}
	else
	{
	ShowMessage(L"Not memory for mparrBuff ") ;
	Abort() ;
	}

	}
 }
 */

 /*
	// парам конструктор 2
 TFight::TFight (const TVessel Vessel, const TTarget Target
		,const double TExtrap, TEtalonSign ETalonSignal, wchar_t *pwcharrFoldReport)
 {
	 mEtalonSign = ETalonSignal;

	 mVessel  = Vessel ;
   mTarget = Target ;
   if( mTarget.mpwcharrFoldReport != NULL)
   {
	mTarget.mQuantPntReport = 0 ;

   }
	if( mVessel.mpwcharrFoldReport != NULL)
   {
	mVessel.mQuantPntReport = 0 ;

   }
   // инициализация фильтра
	//  1.Время привязки оценок
	mTraceFlt.mTf = mTarget.mTraject.mTCur ;

	// 2.Оенка вектора состояния цели в КГСК на момент mTarget.mT
	double arrTemp [ 6 ]= {0.} ;
	memcpy(arrTemp , mTarget.mTraject.marrVectSostGSK, 6 * sizeof(double)) ;
	MtrxMinusMatrx(arrTemp ,  mVessel.marrVectSost,3, 1, mTraceFlt.marrVS_KGSK);

// 3.Оценка вектора состояния в АСК на момент mTf
	 RecalcVect_KGSK_INTO_ASK_True (mTraceFlt.marrVS_KGSK,mTraceFlt.marrVS_ASK,6 ) ;


//  4.Объединенная корреляционная матрица ошибок оценивания в КГСК

	memset(mTraceFlt.marrK_KGSK, 0, 36 * sizeof(double)) ;
	for (int i = 0; i < 3; i++)
	{
	mTraceFlt.marrK_KGSK[ 6 * i + i ] = 1000000. ;
	mTraceFlt.marrK_KGSK[ 6 * i + i + 3] = 100000. ;
	mTraceFlt.marrK_KGSK[ 6 * (i + 3) + i] = 100000. ;
	mTraceFlt.marrK_KGSK[ 6 * (i + 3) + i + 3] = 100000. ;

	}

   //	5.Объединенная корреляционная матрица ошибок оценивания в ASK
	memcpy( mTraceFlt.marrK_ASK,mTraceFlt.marrK_KGSK,36 * sizeof(double));

	// 6. вектор углов
	mTraceFlt.marrMu[0] = mVessel.mSins.mEstQ;
	mTraceFlt.marrMu[1] = mVessel.mSins.mEstPsi ;
	mTraceFlt.marrMu[2] = mVessel.mSins.mEstTet ;
	mTraceFlt.marrMu[3] = mVessel.mDriver.mRealBet ;
	mTraceFlt.marrMu[4] = mVessel.mDriver.mRealEps ;

	//  7.матрица перехода из АСК в ПСК на момент mTf

	calcMatr_ASK_v_PSK( &mTraceFlt.marrMu[3],mTraceFlt.marrF_ASK_V_PSK) ;


 // 8.расширенная матрица перехода из АСК в ПСК на момент mTf
	   createExtendMtrx(mTraceFlt.marrF_ASK_V_PSK, mTraceFlt.marrFExt_ASK_V_PSK) ;

//  9.матрица перехода из ПСК в КГСК на момент mTf
	calcMatr_PSK_v_KGSK(mTraceFlt.marrMu, mTraceFlt.marrF_PSK_V_KGSK) ;


 // 10. расширенная матрица перехода из  ПСК в КГСК мна момент mTf
	 createExtendMtrx(mTraceFlt.marrF_PSK_V_KGSK, mTraceFlt.marrFExt_PSK_V_KGSK) ;

 // 11. матрица перехода из  АСК в КГСК мна момент mTf
	 MtrxMultMatrx(mTraceFlt.marrF_PSK_V_KGSK,3, 3, mTraceFlt.marrF_ASK_V_PSK,3, mTraceFlt.marrF_ASK_V_KGSK) ;

	 // 12. расширенная матрица перехода из  АСК в КГСК мна момент mTf
	 createExtendMtrx(mTraceFlt.marrF_ASK_V_KGSK, mTraceFlt.marrFExt_ASK_V_KGSK) ;
 // 11. признак того, что фильтр прошел инициализацию
	 mTraceFlt.mbInit = true ;

 // ИНИЦИАЛИЗАЦИЯ ПЕРВИЧНЫХ ЦЕЛЕУКЫАЗАНИЙ
	mT =  mTarget.mTraject.mTCur ;
    mTExtrap = TExtrap;
	// дяля  ычисление истинного положения цели в ПСК - АС
			double  arrVSExtrap_PSK_AS [6] = {0};

// расчет первичных целеуказаний (на основе внешних данных о положении цели)
	getTargDes_PSK_True( TExtrap, mTargDesBet, mTargDesEps, mTargDesR, arrVSExtrap_PSK_AS ) ;
//	mVessel.mRadar.mRzv = TFight::Norm3(mTraceFlt.marrVS_KGSK ); //  замер дальности
//	mVessel.mRadar.mVzv = 0;  // замер по углу V (курс)
//	mVessel.mRadar.mUzv = 0;

//	mVessel.mMeasurer.takeRadarData();



	// для формироывания буфера памяти для отчета
	mpwcharrFoldReport = NULL ;
	mparrBuff   = NULL ;
	mQuantPntReport = 0 ;
	mLenMemoryAlloc = 0 ;
	if(pwcharrFoldReport != NULL)
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

	mLenMemoryAlloc = 2000;
	if((mparrBuff  = (double*)malloc( QUANT_COLS_BUFF * mLenMemoryAlloc * sizeof(double))) != NULL)
	{

	memset( mparrBuff,  0,  QUANT_COLS_BUFF * mLenMemoryAlloc * sizeof(double));
	}
	else
	{
	ShowMessage(L"Not memory for mparrBuff ") ;
	Abort() ;
	}

	}
 }
 */

 /*
	// парам конс 3
 TFight::TFight (const TVessel Vessel, const TTarget Target
		,const double TExtrap, TEtalonSign ETalonSignal
		, TEnvironment Environment,wchar_t *pwcharrFoldReport)
 {
	mEtalonSign = ETalonSignal;

	mVessel  = Vessel ;
	mTarget = Target ;
   mEnvironment = Environment;
   if( mTarget.mpwcharrFoldReport != NULL)
   {
	mTarget.mQuantPntReport = 0 ;

   }
	if( mVessel.mpwcharrFoldReport != NULL)
   {
	mVessel.mQuantPntReport = 0 ;

   }
   // инициализация фильтра
	//  1.Время привязки оценок
	mTraceFlt.mTf = mTarget.mTraject.mTCur ;

	// 2.Оенка вектора состояния цели в КГСК на момент mTarget.mT
	double arrTemp [ 6 ]= {0.} ;
	memcpy(arrTemp , mTarget.mTraject.marrVectSostGSK, 6 * sizeof(double)) ;
	MtrxMinusMatrx(arrTemp ,  mVessel.marrVectSost,3, 1, mTraceFlt.marrVS_KGSK);

// 3.Оценка вектора состояния в АСК на момент mTf
	 RecalcVect_KGSK_INTO_ASK_True (mTraceFlt.marrVS_KGSK,mTraceFlt.marrVS_ASK,6 ) ;


//  4.Объединенная корреляционная матрица ошибок оценивания в КГСК

	memset(mTraceFlt.marrK_KGSK, 0, 36 * sizeof(double)) ;
	for (int i = 0; i < 3; i++)
	{
	mTraceFlt.marrK_KGSK[ 6 * i + i ] = 1000000. ;
	mTraceFlt.marrK_KGSK[ 6 * i + i + 3] = 100000. ;
	mTraceFlt.marrK_KGSK[ 6 * (i + 3) + i] = 100000. ;
	mTraceFlt.marrK_KGSK[ 6 * (i + 3) + i + 3] = 100000. ;

	}

   //	5.Объединенная корреляционная матрица ошибок оценивания в ASK
	memcpy( mTraceFlt.marrK_ASK,mTraceFlt.marrK_KGSK,36 * sizeof(double));

	// 6. вектор углов
	mTraceFlt.marrMu[0] = mVessel.mSins.mEstQ;
	mTraceFlt.marrMu[1] = mVessel.mSins.mEstPsi ;
	mTraceFlt.marrMu[2] = mVessel.mSins.mEstTet ;
	mTraceFlt.marrMu[3] = mVessel.mDriver.mRealBet ;
	mTraceFlt.marrMu[4] = mVessel.mDriver.mRealEps ;

	//  7.матрица перехода из АСК в ПСК на момент mTf

	calcMatr_ASK_v_PSK( &mTraceFlt.marrMu[3],mTraceFlt.marrF_ASK_V_PSK) ;


 // 8.расширенная матрица перехода из АСК в ПСК на момент mTf
	   createExtendMtrx(mTraceFlt.marrF_ASK_V_PSK, mTraceFlt.marrFExt_ASK_V_PSK) ;

//  9.матрица перехода из ПСК в КГСК на момент mTf
	calcMatr_PSK_v_KGSK(mTraceFlt.marrMu, mTraceFlt.marrF_PSK_V_KGSK) ;


 // 10. расширенная матрица перехода из  ПСК в КГСК мна момент mTf
	 createExtendMtrx(mTraceFlt.marrF_PSK_V_KGSK, mTraceFlt.marrFExt_PSK_V_KGSK) ;

 // 11. матрица перехода из  АСК в КГСК мна момент mTf
	 MtrxMultMatrx(mTraceFlt.marrF_PSK_V_KGSK,3, 3, mTraceFlt.marrF_ASK_V_PSK,3, mTraceFlt.marrF_ASK_V_KGSK) ;

	 // 12. расширенная матрица перехода из  АСК в КГСК мна момент mTf
	 createExtendMtrx(mTraceFlt.marrF_ASK_V_KGSK, mTraceFlt.marrFExt_ASK_V_KGSK) ;
 // 11. признак того, что фильтр прошел инициализацию
	 mTraceFlt.mbInit = true ;

 // ИНИЦИАЛИЗАЦИЯ ПЕРВИЧНЫХ ЦЕЛЕУКЫАЗАНИЙ
	mT =  mTarget.mTraject.mTCur ;
    mTExtrap = TExtrap;
	// дяля  ычисление истинного положения цели в ПСК - АС
			double  arrVSExtrap_PSK_AS [6] = {0};

// расчет первичных целеуказаний (на основе внешних данных о положении цели)
 	getTargDes_PSK_True( TExtrap, mTargDesBet, mTargDesEps, mTargDesR, arrVSExtrap_PSK_AS ) ;




	// для формироывания буфера памяти для отчета
	mpwcharrFoldReport = NULL ;
	mparrBuff   = NULL ;
	mQuantPntReport = 0 ;
	mLenMemoryAlloc = 0 ;
	if(pwcharrFoldReport != NULL)
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

	mLenMemoryAlloc = 2000;
	if((mparrBuff  = (double*)malloc( QUANT_COLS_BUFF * mLenMemoryAlloc * sizeof(double))) != NULL)
	{

	memset( mparrBuff,  0,  QUANT_COLS_BUFF * mLenMemoryAlloc * sizeof(double));
	}
	else
	{
	ShowMessage(L"Not memory for mparrBuff ") ;
	Abort() ;
	}

	}
 }
 */

 	// парам конс 3
 TFight::TFight (const TVessel Vessel, const TTarget Target
		,const double TExtrap, TEtalonSign ETalonSignal
		, TEnvironment Environment,wchar_t *pwcharrFoldReport)
 {
	mAntCoeff = 1.;
	mEtalonSign = ETalonSignal;

	mVessel  = Vessel ;
	mTarget = Target ;
	mEnvironment = Environment;
   //	mMissSimulated = 0.;
   //	mPereletSimulated = 0.;

   if( mTarget.mpwcharrFoldReport != NULL)
   {
	mTarget.mQuantPntReport = 0 ;

   }
	if( mVessel.mpwcharrFoldReport != NULL)
   {
	mVessel.mQuantPntReport = 0 ;

   }
	 // инициализация фильтра
	 	// 6. вектор углов
	mTraceFlt.marrMu[0] = mVessel.mSins.mEstQ;
	mTraceFlt.marrMu[1] = mVessel.mSins.mEstPsi ;
	mTraceFlt.marrMu[2] = mVessel.mSins.mEstTet ;
	mTraceFlt.marrMu[3] = mVessel.mDriver.mRealBet ;
	mTraceFlt.marrMu[4] = mVessel.mDriver.mRealEps ;
	//  1.Время привязки оценок
	mTraceFlt.mTf = mTarget.mTraject.mTCur ;

	// 2.Оенка вектора состояния цели в КГСК на момент mTarget.mT
	double arrTemp [ 6 ]= {0.} ;
	memcpy(arrTemp , mTarget.mTraject.marrVectSostGSK, 6 * sizeof(double)) ;
	MtrxMinusMatrx(arrTemp ,  mVessel.marrVectSost,3, 1, mTraceFlt.marrVS_KGSK);

// 3.Оценка вектора состояния в АСК на момент mTf
	 RecalcVect_KGSK_INTO_ASK_True (mTraceFlt.marrVS_KGSK,mTraceFlt.marrVS_ASK,6 ) ;




	//  7.матрица перехода из АСК в ПСК на момент mTf

	calcMatr_ASK_v_PSK( &mTraceFlt.marrMu[3],mTraceFlt.marrF_ASK_V_PSK) ;


 // 8.расширенная матрица перехода из АСК в ПСК на момент mTf
	   createExtendMtrx(mTraceFlt.marrF_ASK_V_PSK, mTraceFlt.marrFExt_ASK_V_PSK) ;

//  9.матрица перехода из ПСК в КГСК на момент mTf
	calcMatr_PSK_v_KGSK(mTraceFlt.marrMu, mTraceFlt.marrF_PSK_V_KGSK) ;


 // 10. расширенная матрица перехода из  ПСК в КГСК мна момент mTf
	 createExtendMtrx(mTraceFlt.marrF_PSK_V_KGSK, mTraceFlt.marrFExt_PSK_V_KGSK) ;

 // 11. матрица перехода из  АСК в КГСК мна момент mTf
	 MtrxMultMatrx(mTraceFlt.marrF_PSK_V_KGSK,3, 3, mTraceFlt.marrF_ASK_V_PSK,3, mTraceFlt.marrF_ASK_V_KGSK) ;

	 // 12. расширенная матрица перехода из  АСК в КГСК мна момент mTf
	 createExtendMtrx(mTraceFlt.marrF_ASK_V_KGSK, mTraceFlt.marrFExt_ASK_V_KGSK) ;

	 //  4.Объединенная корреляционная матрица ошибок оценивания в АСК
		 // вывод коорреляционной матрицы ошибок оцениывания дальности на стационый режим
		 TStabSyst2 StabSyst(  mVessel.mControlSyst.mFiltT,mTarget.mTraject.marrSigW[1], 0.
	,mVessel.mFar_2D.mDistSKZ ) ;
	double  arrKStab [4] = {1000000., 1000000., 1000000., 10000000.
	};
	double arrPStab [2] = {0.5,0.5
	};
	StabSyst.CalcCorrMatrGolubev(arrKStab,arrPStab
		,NULL, NULL, NULL
		,NULL,NULL) ;



	memset(mTraceFlt.marrK_ASK, 0, 36 * sizeof(double)) ;
	for (int i = 0; i < 3; i++)
	{
	mTraceFlt.marrK_ASK[ 6 * i + i ] = 40000. ;
	mTraceFlt.marrK_ASK[ 6 * i + i + 3] = 5000. ;
	mTraceFlt.marrK_ASK[ 6 * (i + 3) + i] = 5000. ;
	mTraceFlt.marrK_ASK[ 6 * (i + 3) + i + 3] = 2500. ;
	}
	mTraceFlt.marrK_ASK[ 6  + 1 ] =  arrKStab[0];
	mTraceFlt.marrK_ASK[ 6  + 1 + 3] = arrKStab[1];
	mTraceFlt.marrK_ASK[ 6 * (1 + 3) + 1] = arrKStab[1];
	mTraceFlt.marrK_ASK[ 6 * (1 + 3) + 1 + 3] = arrKStab[3];

   //	5.Объединенная корреляционная матрица ошибок оценивания в ASK

	double arrTEmp0[36] = {0.};
	MtrxMultMatrx(mTraceFlt.marrFExt_ASK_V_KGSK, 6,6, mTraceFlt.marrK_ASK, 6, arrTEmp0);
	MtrxMultMatrxTransp(arrTEmp0,6, 6, mTraceFlt.marrFExt_ASK_V_KGSK,6, mTraceFlt.marrK_KGSK) ;


 // 11. признак того, что фильтр прошел инициализацию
	 mTraceFlt.mbInit = true ;



 // ИНИЦИАЛИЗАЦИЯ ПЕРВИЧНЫХ ЦЕЛЕУКЫАЗАНИЙ
	mT =  mTarget.mTraject.mTCur ;
    mTExtrap = TExtrap;
	// дяля  ычисление истинного положения цели в ПСК - АС
			double  arrVSExtrap_PSK_AS [6] = {0};

// расчет первичных целеуказаний (на основе внешних данных о положении цели)
 	getTargDes_PSK_True( TExtrap, mTargDesBet, mTargDesEps, mTargDesR, arrVSExtrap_PSK_AS ) ;




	// для формироывания буфера памяти для отчета
	mpwcharrFoldReport = NULL ;
	mparrBuff   = NULL ;
	mQuantPntReport = 0 ;
	mLenMemoryAlloc = 0 ;
	if(pwcharrFoldReport != NULL)
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

	mLenMemoryAlloc = 2000;
	if((mparrBuff  = (double*)malloc( QUANT_COLS_BUFF * mLenMemoryAlloc * sizeof(double))) != NULL)
	{

	memset( mparrBuff,  0,  QUANT_COLS_BUFF * mLenMemoryAlloc * sizeof(double));
	}
	else
	{
	ShowMessage(L"Not memory for mparrBuff ") ;
	Abort() ;
	}

	}
 }


	// парам конс 4
 TFight::TFight (const TVessel Vessel, const TTarget Target
		,const double TExtrap, TEtalonSign ETalonSignal
		, TEnvironment Environment, const double VAlAntCoeff, wchar_t *pwcharrFoldReport)
 {
	mAntCoeff =  VAlAntCoeff;
	mEtalonSign = ETalonSignal;

	mVessel  = Vessel ;
	mTarget = Target ;
	mEnvironment = Environment;
   //	mMissSimulated = 0.;
   //	mPereletSimulated = 0.;

   if( mTarget.mpwcharrFoldReport != NULL)
   {
	mTarget.mQuantPntReport = 0 ;

   }
	if( mVessel.mpwcharrFoldReport != NULL)
   {
	mVessel.mQuantPntReport = 0 ;

   }
	 // инициализация фильтра
	 	// 6. вектор углов
	mTraceFlt.marrMu[0] = mVessel.mSins.mEstQ;
	mTraceFlt.marrMu[1] = mVessel.mSins.mEstPsi ;
	mTraceFlt.marrMu[2] = mVessel.mSins.mEstTet ;
	mTraceFlt.marrMu[3] = mVessel.mDriver.mRealBet ;
	mTraceFlt.marrMu[4] = mVessel.mDriver.mRealEps ;
	//  1.Время привязки оценок
	mTraceFlt.mTf = mTarget.mTraject.mTCur ;

	// 2.Оенка вектора состояния цели в КГСК на момент mTarget.mT
	double arrTemp [ 6 ]= {0.} ;
	memcpy(arrTemp , mTarget.mTraject.marrVectSostGSK, 6 * sizeof(double)) ;
	MtrxMinusMatrx(arrTemp ,  mVessel.marrVectSost,3, 1, mTraceFlt.marrVS_KGSK);

// 3.Оценка вектора состояния в АСК на момент mTf
	 RecalcVect_KGSK_INTO_ASK_True (mTraceFlt.marrVS_KGSK,mTraceFlt.marrVS_ASK,6 ) ;




	//  7.матрица перехода из АСК в ПСК на момент mTf

	calcMatr_ASK_v_PSK( &mTraceFlt.marrMu[3],mTraceFlt.marrF_ASK_V_PSK) ;


 // 8.расширенная матрица перехода из АСК в ПСК на момент mTf
	   createExtendMtrx(mTraceFlt.marrF_ASK_V_PSK, mTraceFlt.marrFExt_ASK_V_PSK) ;

//  9.матрица перехода из ПСК в КГСК на момент mTf
	calcMatr_PSK_v_KGSK(mTraceFlt.marrMu, mTraceFlt.marrF_PSK_V_KGSK) ;


 // 10. расширенная матрица перехода из  ПСК в КГСК мна момент mTf
	 createExtendMtrx(mTraceFlt.marrF_PSK_V_KGSK, mTraceFlt.marrFExt_PSK_V_KGSK) ;

 // 11. матрица перехода из  АСК в КГСК мна момент mTf
	 MtrxMultMatrx(mTraceFlt.marrF_PSK_V_KGSK,3, 3, mTraceFlt.marrF_ASK_V_PSK,3, mTraceFlt.marrF_ASK_V_KGSK) ;

	 // 12. расширенная матрица перехода из  АСК в КГСК мна момент mTf
	 createExtendMtrx(mTraceFlt.marrF_ASK_V_KGSK, mTraceFlt.marrFExt_ASK_V_KGSK) ;

	 //  4.Объединенная корреляционная матрица ошибок оценивания в АСК
		 // вывод коорреляционной матрицы ошибок оцениывания дальности на стационый режим
		 TStabSyst2 StabSyst(  mVessel.mControlSyst.mFiltT,mTarget.mTraject.marrSigW[1], 0.
	,mVessel.mFar_2D.mDistSKZ ) ;
	double  arrKStab [4] = {1000000., 1000000., 1000000., 10000000.
	};
	double arrPStab [2] = {0.5,0.5
	};
	StabSyst.CalcCorrMatrGolubev(arrKStab,arrPStab
		,NULL, NULL, NULL
		,NULL,NULL) ;



	memset(mTraceFlt.marrK_ASK, 0, 36 * sizeof(double)) ;
	for (int i = 0; i < 3; i++)
	{
	mTraceFlt.marrK_ASK[ 6 * i + i ] = 40000. ;
	mTraceFlt.marrK_ASK[ 6 * i + i + 3] = 5000. ;
	mTraceFlt.marrK_ASK[ 6 * (i + 3) + i] = 5000. ;
	mTraceFlt.marrK_ASK[ 6 * (i + 3) + i + 3] = 2500. ;
	}
	mTraceFlt.marrK_ASK[ 6  + 1 ] =  arrKStab[0];
	mTraceFlt.marrK_ASK[ 6  + 1 + 3] = arrKStab[1];
	mTraceFlt.marrK_ASK[ 6 * (1 + 3) + 1] = arrKStab[1];
	mTraceFlt.marrK_ASK[ 6 * (1 + 3) + 1 + 3] = arrKStab[3];

   //	5.Объединенная корреляционная матрица ошибок оценивания в ASK

	double arrTEmp0[36] = {0.};
	MtrxMultMatrx(mTraceFlt.marrFExt_ASK_V_KGSK, 6,6, mTraceFlt.marrK_ASK, 6, arrTEmp0);
	MtrxMultMatrxTransp(arrTEmp0,6, 6, mTraceFlt.marrFExt_ASK_V_KGSK,6, mTraceFlt.marrK_KGSK) ;


 // 11. признак того, что фильтр прошел инициализацию
	 mTraceFlt.mbInit = true ;



 // ИНИЦИАЛИЗАЦИЯ ПЕРВИЧНЫХ ЦЕЛЕУКЫАЗАНИЙ
	mT =  mTarget.mTraject.mTCur ;
    mTExtrap = TExtrap;
	// дяля  ычисление истинного положения цели в ПСК - АС
			double  arrVSExtrap_PSK_AS [6] = {0};

// расчет первичных целеуказаний (на основе внешних данных о положении цели)
 	getTargDes_PSK_True( TExtrap, mTargDesBet, mTargDesEps, mTargDesR, arrVSExtrap_PSK_AS ) ;




	// для формироывания буфера памяти для отчета
	mpwcharrFoldReport = NULL ;
	mparrBuff   = NULL ;
	mQuantPntReport = 0 ;
	mLenMemoryAlloc = 0 ;
	if(pwcharrFoldReport != NULL)
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

	mLenMemoryAlloc = 2000;
	if((mparrBuff  = (double*)malloc( QUANT_COLS_BUFF * mLenMemoryAlloc * sizeof(double))) != NULL)
	{

	memset( mparrBuff,  0,  QUANT_COLS_BUFF * mLenMemoryAlloc * sizeof(double));
	}
	else
	{
	ShowMessage(L"Not memory for mparrBuff ") ;
	Abort() ;
	}

	}
 }
  /*
	// парам конс 5
 TFight::TFight (const TVessel Vessel, const TTarget Target
		,const double TExtrap, TEtalonSign ETalonSignal
		, TEnvironment Environment,const double  VAlMissSimulated,const double  VAlNedoletSimulated
		,const double VAlAntCoeff, wchar_t *pwcharrFoldReport)
 {
	mAntCoeff =  VAlAntCoeff;
	mEtalonSign = ETalonSignal;

	mVessel  = Vessel ;
	mTarget = Target ;
	mEnvironment = Environment;
	mMissSimulated = VAlMissSimulated;
	mPereletSimulated = VAlNedoletSimulated;

   if( mTarget.mpwcharrFoldReport != NULL)
   {
	mTarget.mQuantPntReport = 0 ;

   }
	if( mVessel.mpwcharrFoldReport != NULL)
   {
	mVessel.mQuantPntReport = 0 ;

   }
	 // инициализация фильтра
	 	// 6. вектор углов
	mTraceFlt.marrMu[0] = mVessel.mSins.mEstQ;
	mTraceFlt.marrMu[1] = mVessel.mSins.mEstPsi ;
	mTraceFlt.marrMu[2] = mVessel.mSins.mEstTet ;
	mTraceFlt.marrMu[3] = mVessel.mDriver.mRealBet ;
	mTraceFlt.marrMu[4] = mVessel.mDriver.mRealEps ;
	//  1.Время привязки оценок
	mTraceFlt.mTf = mTarget.mTraject.mTCur ;

	// 2.Оенка вектора состояния цели в КГСК на момент mTarget.mT
	double arrTemp [ 6 ]= {0.} ;
	memcpy(arrTemp , mTarget.mTraject.marrVectSostGSK, 6 * sizeof(double)) ;
	MtrxMinusMatrx(arrTemp ,  mVessel.marrVectSost,3, 1, mTraceFlt.marrVS_KGSK);

// 3.Оценка вектора состояния в АСК на момент mTf
	 RecalcVect_KGSK_INTO_ASK_True (mTraceFlt.marrVS_KGSK,mTraceFlt.marrVS_ASK,6 ) ;




	//  7.матрица перехода из АСК в ПСК на момент mTf

	calcMatr_ASK_v_PSK( &mTraceFlt.marrMu[3],mTraceFlt.marrF_ASK_V_PSK) ;


 // 8.расширенная матрица перехода из АСК в ПСК на момент mTf
	   createExtendMtrx(mTraceFlt.marrF_ASK_V_PSK, mTraceFlt.marrFExt_ASK_V_PSK) ;

//  9.матрица перехода из ПСК в КГСК на момент mTf
	calcMatr_PSK_v_KGSK(mTraceFlt.marrMu, mTraceFlt.marrF_PSK_V_KGSK) ;


 // 10. расширенная матрица перехода из  ПСК в КГСК мна момент mTf
	 createExtendMtrx(mTraceFlt.marrF_PSK_V_KGSK, mTraceFlt.marrFExt_PSK_V_KGSK) ;

 // 11. матрица перехода из  АСК в КГСК мна момент mTf
	 MtrxMultMatrx(mTraceFlt.marrF_PSK_V_KGSK,3, 3, mTraceFlt.marrF_ASK_V_PSK,3, mTraceFlt.marrF_ASK_V_KGSK) ;

	 // 12. расширенная матрица перехода из  АСК в КГСК мна момент mTf
	 createExtendMtrx(mTraceFlt.marrF_ASK_V_KGSK, mTraceFlt.marrFExt_ASK_V_KGSK) ;

	 //  4.Объединенная корреляционная матрица ошибок оценивания в АСК
		 // вывод коорреляционной матрицы ошибок оцениывания дальности на стационый режим
		 TStabSyst2 StabSyst(  mVessel.mControlSyst.mFiltT,mTarget.mTraject.marrSigW[1], 0.
	,mVessel.mFar_2D.mDistSKZ ) ;
	double  arrKStab [4] = {1000000., 1000000., 1000000., 10000000.
	};
	double arrPStab [2] = {0.5,0.5
	};
	StabSyst.CalcCorrMatrGolubev(arrKStab,arrPStab
		,NULL, NULL, NULL
		,NULL,NULL) ;



	memset(mTraceFlt.marrK_ASK, 0, 36 * sizeof(double)) ;
	for (int i = 0; i < 3; i++)
	{
	mTraceFlt.marrK_ASK[ 6 * i + i ] = 40000. ;
	mTraceFlt.marrK_ASK[ 6 * i + i + 3] = 5000. ;
	mTraceFlt.marrK_ASK[ 6 * (i + 3) + i] = 5000. ;
	mTraceFlt.marrK_ASK[ 6 * (i + 3) + i + 3] = 2500. ;
	}
	mTraceFlt.marrK_ASK[ 6  + 1 ] =  arrKStab[0];
	mTraceFlt.marrK_ASK[ 6  + 1 + 3] = arrKStab[1];
	mTraceFlt.marrK_ASK[ 6 * (1 + 3) + 1] = arrKStab[1];
	mTraceFlt.marrK_ASK[ 6 * (1 + 3) + 1 + 3] = arrKStab[3];

   //	5.Объединенная корреляционная матрица ошибок оценивания в ASK

	double arrTEmp0[36] = {0.};
	MtrxMultMatrx(mTraceFlt.marrFExt_ASK_V_KGSK, 6,6, mTraceFlt.marrK_ASK, 6, arrTEmp0);
	MtrxMultMatrxTransp(arrTEmp0,6, 6, mTraceFlt.marrFExt_ASK_V_KGSK,6, mTraceFlt.marrK_KGSK) ;


 // 11. признак того, что фильтр прошел инициализацию
	 mTraceFlt.mbInit = true ;



 // ИНИЦИАЛИЗАЦИЯ ПЕРВИЧНЫХ ЦЕЛЕУКЫАЗАНИЙ
	mT =  mTarget.mTraject.mTCur ;
    mTExtrap = TExtrap;
	// дяля  ычисление истинного положения цели в ПСК - АС
			double  arrVSExtrap_PSK_AS [6] = {0};

// расчет первичных целеуказаний (на основе внешних данных о положении цели)
 	getTargDes_PSK_True( TExtrap, mTargDesBet, mTargDesEps, mTargDesR, arrVSExtrap_PSK_AS ) ;




	// для формироывания буфера памяти для отчета
	mpwcharrFoldReport = NULL ;
	mparrBuff   = NULL ;
	mQuantPntReport = 0 ;
	mLenMemoryAlloc = 0 ;
	if(pwcharrFoldReport != NULL)
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

	mLenMemoryAlloc = 2000;
	if((mparrBuff  = (double*)malloc( QUANT_COLS_BUFF * mLenMemoryAlloc * sizeof(double))) != NULL)
	{

	memset( mparrBuff,  0,  QUANT_COLS_BUFF * mLenMemoryAlloc * sizeof(double));
	}
	else
	{
	ShowMessage(L"Not memory for mparrBuff ") ;
	Abort() ;
	}

	}
 }

 */
// перечсчет вектора сотояния из KGSK в ASK   истинную
// есдли  lenarrKGSK == 6 , то пересчитывается положение и скорость
// если   lenarrKGSK == 3   , тот пересчитывается только положение
// на вход подается вектиор состоящий из 3 или 6 координат.
// первые 3 координаты представляют из себя  положение точки в КГС
// последние 3 координаты представляют из себя скорость точки в КГСК
 void  __fastcall TFight::RecalcVect_KGSK_INTO_ASK_True (double *arrKGSK,double *arrASK,int lenarrKGSK )
 {
    // персчет вектора положения в АСК
	// создание вектора углов ориентации коррабля
	double arrMu[5] = {0} ;
/*	arrMu[0] = mVessel.mQ;
	arrMu[1] = mVessel.mPsi;
	arrMu[2] = mVessel.mTet ;
	arrMu[3] = mVessel.mDriver.mRealBet ;
	arrMu[4] = mVessel.mDriver.mRealEps ; */
	mVessel.calcVectTrueDeckAngles_For_Far(arrMu);
	// создание матирицы перехода из ПСК в КГСК
	arrMu[3] = 0 ;
	arrMu[4] = 0;
	double matrPereh_PSK_V_KGSK[9] = {0} ;
	calcMatr_ASK_v_KGSK(arrMu, matrPereh_PSK_V_KGSK);
	// создание вектора положения с ПСК-ЦТ
		double arrPosPSK_CT[3] = {0} ;
		MtrxTranspMultMatrx(matrPereh_PSK_V_KGSK,3, 3, arrKGSK,1, arrPosPSK_CT) ;
	// создание вектора положение в ПСК-АС
		double arrPosPSK_AS[3] = {0} ;
		for (int i =0; i < 3; i++) arrPosPSK_AS[i] =  arrPosPSK_CT [i] - mVessel.marrParral[i] ;



	// создание матрицы перехода из  АСК  в ПСК
	mVessel.calcVectTrueDeckAngles_For_Far(arrMu);
	arrMu[0] = 0;
	arrMu[1] = 0;
	arrMu[2] = 0 ;
 //	arrMu[3] = mVessel.mDriver.mRealBet ;
 //	arrMu[4] = mVessel.mDriver.mRealEps ;
	double matrPereh_ASK_V_PSK[9] = {0} ;
	calcMatr_ASK_v_KGSK(arrMu, matrPereh_ASK_V_PSK) ;
	// создание вектора положения в истинной линейной АСК
	MtrxTranspMultMatrx( matrPereh_ASK_V_PSK,3, 3, arrPosPSK_AS,1, arrASK) ;
	// пересчет скорости
	if (lenarrKGSK == 6)
	{
		// создание марицы перехода из АСК в КГСК
		double matrPereh_ASK_V_KGSK[9] = {0.} ;
		//arrMu[0] = mVessel.mQ;
		//arrMu[1] = mVessel.mPsi;
		//	arrMu[2] = mVessel.mTet ;
		mVessel.calcVectTrueDeckAngles_For_Far(arrMu);
		calcMatr_ASK_v_KGSK(arrMu, matrPereh_ASK_V_KGSK);
		MtrxTranspMultMatrx( matrPereh_ASK_V_KGSK,3, 3, &arrKGSK[3],1, &arrASK[3]) ;
	}
 }
// перечсчет вектора сотояния из ПСК-АС в ГСК
// есдли  lenarrVS == 6 , то пересчитывается положение и скорость
// если   lenarrVS == 3   , тот пересчитывается только положение
// на вход подается вектиор состоящий из 3 или 6 координат.
// первые 3 координаты представляют из себя  положение точки в ПСК-АС
// последние 3 координаты представляют из себя скорость точки в ПСК-АС
 void  __fastcall TFight::RecalcVect_PSK_AS_INTO_GSK (double *arrPSK_AS,double *arrGSK,int lenarrVS )
 {
   double arr0[3] = {0.},arr1[3] = {0.};

	//1.  персчет вектора положения

	 MtrxSumMatrx(arrPSK_AS, mVessel.marrParral  ,3, 1, arr0) ;
	 MtrxMultMatrx(mTraceFlt.marrF_PSK_V_KGSK,3, 3, arr0,1, arr1) ;
	 MtrxSumMatrx(arr1, mVessel.marrVectSost  ,3, 1, arrGSK) ;
	 //2. пересчет скорости
	  if (lenarrVS == 6) MtrxMultMatrx(mTraceFlt.marrF_PSK_V_KGSK,3, 3, &arrPSK_AS[3],1, &arrGSK[3]) ;

}

//  пересчет вектора из системы коорд ГСК в ПСК-АС
// есдли  lenarrVS == 6 , то пересчитывается положение и скорость
// если   lenarrVS == 3   , тот пересчитывается только положение
// на вход подается вектиор состоящий из 3 или 6 координат.
// первые 3 координаты представляют из себя  положение точки в ГСК
// последние 3 координаты представляют из себя скорость точки в ГСК

void  __fastcall TFight::RecalcVect_GSK_INTO_PSK_AS_True (double *arrGSK,double *arrPSK_AS,int lenarrKGSK )
{
   // персчет вектора состояния по положениею из ГСК в ПСК -ЦТ
	double  arrTemp [ 3 ] ={0.}  ;
	RecalcVect_GSK_INTO_PSK_CT_True (arrGSK,arrPSK_AS, lenarrKGSK ) ;

	// пересчет  вектора состояния по положениею из  ПСК -ЦТ  в  ПСК -АС
   MtrxMinusMatrx(arrPSK_AS, mVessel.marrParral,3, 1, arrTemp ) ;
   memcpy(arrPSK_AS, arrTemp, 3 * sizeof(double));
   return ;

}

//  пересчет вектора из истинной системы коорд ГСК в истинную ПСК-ЦТ
// есдли  lenarrVS == 6 , то пересчитывается положение и скорость
// если   lenarrVS == 3   , тот пересчитывается только положение
// на вход подается вектиор состоящий из 3 или 6 координат.
// первые 3 координаты представляют из себя  положение точки в ГСК
// последние 3 координаты представляют из себя скорость точки в ГСК

void  __fastcall TFight::RecalcVect_GSK_INTO_PSK_CT_True (double *arrGSK,double *arrPSK_CT,int lenarrKGSK )
{
   // персчет вектора состояния по положениею из ГСК в KGSK
   double arrKGSK[3] = {0};
   MtrxMinusMatrx(arrGSK, mVessel.marrVectSost,3, 1, arrKGSK ) ;
   // пересчет вектора состояния из КГСК в ПСК-ЦТ
	 double arrMu[5] = {0.} , matrPereh_PSK_V_KGSK [9] = {0.} ;
	// arrMu[0] =  mVessel.mQ ;
	// arrMu[1] =  mVessel.mPsi ;
	// arrMu[2] =  mVessel.mTet ;
	 mVessel.calcVectTrueDeckAngles_For_Far(arrMu);
   calcMatr_PSK_v_KGSK(arrMu, matrPereh_PSK_V_KGSK) ;
   MtrxTranspMultMatrx(matrPereh_PSK_V_KGSK,3, 3, arrKGSK,1, arrPSK_CT) ;

   if (lenarrKGSK == 6 ) MtrxTranspMultMatrx(matrPereh_PSK_V_KGSK,3, 3, &arrGSK [3],1, &arrPSK_CT[3]) ;

}

/*
void  __fastcall TFight::RecalcVect_GSK_INTO_PSK_AS_True (double *arrGSK,double *arrPSK_AS,int lenarrKGSK )
{
   // персчет вектора состояния по положениею из ГСК в KGSK
   double arrKGSK[6] = {0}, arrTemp0 [ 6 ] = {0};
   memcpy( arrTemp0, arrGSK , lenarrKGSK * sizeof(double)) ;
   MtrxMinusMatrx(arrTemp0, mVessel.marrVectSost,3, 1, arrKGSK ) ;
   // пересчет вектора состояния из КГСК в ПСК-ЦТ
   double arrPSK_CT [6] = {0} ;
   MtrxTranspMultMatrx(mTraceFlt.marrF_PSK_V_KGSK,3, 3, arrKGSK,1, arrPSK_CT) ;
   // персчет вектора состояния по положениею из ПСК-ЦТ в ПСК-АС
   double arrPSK_AS_Temp[6] = {0} ;
   MtrxMinusMatrx(arrPSK_CT, mVessel.marrParral,3, 1, arrPSK_AS_Temp ) ;

   if (lenarrKGSK == 6 ) MtrxTranspMultMatrx(mTraceFlt.marrF_PSK_V_KGSK,3, 3, &arrKGSK [3],1, &arrPSK_AS_Temp[3]) ;
	memcpy( arrPSK_AS, arrPSK_AS_Temp, lenarrKGSK * sizeof(double)) ;


}
  */
// расчет углов целеуазаний valTargDesEps и  valTargDesBet на момент времени valTExtrap в истинной системе координат ПСК(И)
	 void  __fastcall TFight::getTargDes_PSK_True(const double valTExtrap
   ,double &valTargDesBet,double &valTargDesEps,double &valTargDesR, double *arrVSExtrap_PSK)
{
	 // пресчет оценок фильра из АСК в ПСК
	 double arrSEst_PSK_AS[6] = {0.} ;
	 MtrxMultMatrx(mTraceFlt.marrFExt_ASK_V_PSK  ,6, 6, mTraceFlt.marrVS_ASK,1, arrSEst_PSK_AS);
	 // линейеый прогноз углов Q, Psi, Tet
	 double arrMu [5] = {0.} ;

	 arrMu [0] =  mVessel.mSins.mEstQ + (valTExtrap - mTraceFlt.mTf) * mVessel.mSins.mEstVQ ;
	 arrMu [1] =  mVessel.mSins.mEstPsi + (valTExtrap - mTraceFlt.mTf) * mVessel.mSins.mEstVPsi ;
	 arrMu [2] =  mVessel.mSins.mEstTet + (valTExtrap - mTraceFlt.mTf) * mVessel.mSins.mEstVTet ;

	 // нахождение матрицы перехода из ПСК спрогнозированной в КГСК - arrMPereh_PSKExtr_V_KGSK_Extend
	 double arrMPereh_PSKExtr_V_KGSK[9]={0.},  arrMPereh_PSKExtr_V_KGSK_Extend[36]={0.};
	 calcMatr_ASK_v_KGSK(arrMu, arrMPereh_PSKExtr_V_KGSK) ;
	 createExtendMtrx(arrMPereh_PSKExtr_V_KGSK, arrMPereh_PSKExtr_V_KGSK_Extend);

	 // расчет матрицы FT (t + h)* L(h) * F(t)
	 double arrL6[36] = {0}, arr0[36] ={0},arr1[36]={0}, arrVect0[6]= {0},arrVect1[6]= {0},arrVect2[6]= {0} ;
	 CreateMatrL6 (valTExtrap - mTraceFlt.mTf, arrL6);
	 MtrxTranspMultMatrx(arrMPereh_PSKExtr_V_KGSK_Extend,6, 6, arrL6,6, arr0);
	 MtrxMultMatrx(arr0,6, 6, mTraceFlt.marrFExt_PSK_V_KGSK  ,6, arr1) ;

	 // рапсчет  FT* L * F * S
	 MtrxMultMatrx(arr1,6, 6, arrSEst_PSK_AS,1, arrVect0)  ;
	 // расчет f
	 calcVect_f(valTExtrap - mTraceFlt.mTf, arrMu [0]
	   ,mVessel.mSins.mEstVVess ,  mVessel.mSins.mEstVH * (valTExtrap - mTraceFlt.mTf) , arrVect1) ;

	 //
	  MtrxTranspMultMatrx(arrMPereh_PSKExtr_V_KGSK_Extend,6, 6, arrVect1,1, arrVect2)  ;

	// результат  экстраполяции в ПСК
	MtrxSumMatrx(arrVect0, arrVect2,6, 1, arrVSExtrap_PSK) ;

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 recalcCoord_INTO_Spherical(arrVSExtrap_PSK, valTargDesR, valTargDesBet, valTargDesEps) ;

}

// расчет углов целеуазаний valTargDesEps и  valTargDesBet на момент времени valTExtrap в истинной системе координат ПСК(И)
// INPUT:
// arrVS_KGSK - вектор состояния в КГСК на текущий момоент воемени
// valDelTExtrap - время на которое считаются целеуазания (время экстраполяции - текущее время))
// OUTPUT:
// valTargDesBet - целеуказание по углу Bet
// valTargDesEps  - целеуказание по углу Eps
// valTargDesR    - целеуказание по дальности
// arrVSExtrap_PSK - оценка вектора состояния в ПСК на момент  valTExtrap
   void  __fastcall TFight::getTargDes_From_KGSK_PSK_AS_True(const double valDelTExtrap, double *arrVS_KGSK
   ,double &valTargDesBet,double &valTargDesEps,double &valTargDesR, double *arrVSExtrap_PSK)
{
	  // экстраполяция вектора состояния цели в КГСК на момент  valTExtrap
	   double arrVS_KGSK_Extr[9] = {0.} ; // экстрап векторор сост в КГСК
	   memcpy( arrVS_KGSK_Extr, arrVS_KGSK, 6 * sizeof(double)) ;
	   for (int i = 0; i < 3; i++) arrVS_KGSK_Extr[i] +=  valDelTExtrap *  arrVS_KGSK_Extr[i + 3]  ;


	 // линейеый прогноз углов Q, Psi, Tet
	 double arrMu [5] = {0.} ;
	 arrMu [0] =  mVessel.mSins.mEstQ   + valDelTExtrap * mVessel.mSins.mEstVQ ;
	 arrMu [1] =  mVessel.mSins.mEstPsi + valDelTExtrap * mVessel.mSins.mEstVPsi ;
	 arrMu [2] =  mVessel.mSins.mEstTet + valDelTExtrap * mVessel.mSins.mEstVTet ;

	 // нахождение матрицы перехода из ПСК спрогнозированной в КГСК - arrMPereh_PSKExtr_V_KGSK_Extend
	 double arrMPereh_PSKExtr_V_KGSK[9]={0.} ; //,  arrMPereh_PSKExtr_V_KGSK_Extend[36]={0.};
	 calcMatr_ASK_v_KGSK(arrMu, arrMPereh_PSKExtr_V_KGSK) ;
// createExtendMtrx(arrMPereh_PSKExtr_V_KGSK, arrMPereh_PSKExtr_V_KGSK_Extend);

	 // пересчет вектора состояния в спрогнозированную ПСК
	 double   arrVect1[6]= {0},arrVect2[6]= {0} ;

	// MtrxTranspMultMatrx(arrMPereh_PSKExtr_V_KGSK_Extend,6, 6, arrVS_KGSK_Extr,6, arrVSExtrap_PSK);
	  MtrxTranspMultMatrx(arrMPereh_PSKExtr_V_KGSK,3, 3, arrVS_KGSK_Extr,1, arrVSExtrap_PSK);


	 // расчет f  (поправки на ход корабля)
	   calcVect_f(valDelTExtrap, arrMu [0]
	   ,mVessel.mSins.mEstVVess ,  mVessel.mSins.mEstVH * valDelTExtrap  , arrVect1) ;

	 //
	 // MtrxTranspMultMatrx(arrMPereh_PSKExtr_V_KGSK_Extend,6, 6, arrVect1,1, arrVect2)  ;
	  MtrxTranspMultMatrx(arrMPereh_PSKExtr_V_KGSK,3, 3, arrVect1,1, arrVect2)  ;

	// результат  экстраполяции в ПСК
	arrVSExtrap_PSK [0] +=  (arrVect2 [0] - mVessel.marrParral[0]);
	arrVSExtrap_PSK [1] +=  (arrVect2 [1] - mVessel.marrParral[1]);
	arrVSExtrap_PSK [2] +=  (arrVect2 [2] - mVessel.marrParral[2]);

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 recalcCoord_INTO_Spherical(arrVSExtrap_PSK, valTargDesR, valTargDesBet, valTargDesEps) ;

	 return ;
}

// пересчет векотра  из кажущейся АСК в кажушуюся ГСК
void  __fastcall TFight::RecalcVect_ASK_Seem_INTO_GSK (double *arrASK,double *arrGSK,int lenarrASK )
{
	double arrMu[5] ={0} ;
	double matrPereh_ASK_V_PSK [9] = {0}, matrPereh_PSK_V_KGSK[9] = {0};
	double arrT0[6] = {0},arrT1[6] = {0},arrT2[6] = {0} ;
	arrMu[0] = mVessel.mSins.mEstQ;
	arrMu[1] = mVessel.mSins.mEstPsi ;
	arrMu[2] = mVessel.mSins.mEstTet ;
	arrMu[3] = mVessel.mDriver.mEstBet ;
	arrMu[4] = mVessel.mDriver.mEstEps ;
// 1. пересчет вектора положения
	// переход в ПСК-АС
	calcMatr_ASK_v_PSK( &arrMu[3], matrPereh_ASK_V_PSK) ;
	MtrxMultMatrx(matrPereh_ASK_V_PSK,3,3, arrASK,1, arrT0) ;
	///
	// переъход в ПСК-ЦТ
	MtrxSumMatrx(arrT0, mVessel.marrParral, 3, 1, arrT1 ) ;

	///
	// переход в КГСК
	calcMatr_PSK_v_KGSK(arrMu, matrPereh_PSK_V_KGSK) ;
	MtrxMultMatrx(matrPereh_PSK_V_KGSK,3,3, arrT1, 1, arrT2) ;

	///

	// переход в ГСК
	MtrxSumMatrx(arrT2, mVessel.marrEstVectSost, 3, 1, arrGSK ) ;
	if (lenarrASK < 6)return ;
// 2. пересчет по скорости
	double  matrPereh_ASK_V_GSK[9];
	MtrxMultMatrx(matrPereh_PSK_V_KGSK,3,3, matrPereh_ASK_V_PSK,3,matrPereh_ASK_V_GSK) ;
	MtrxMultMatrx(matrPereh_ASK_V_GSK,3,3, &arrASK[3], 1, & arrGSK [3]) ;
	if (lenarrASK < 9)return ;
// 3. пересчет по ускорению
	MtrxMultMatrx(matrPereh_ASK_V_GSK,3,3, &arrASK[6], 1, & arrGSK [6]) ;

}
// пересчет векотра  из кажущейся АСК в кажушуюся КГСК
void  __fastcall TFight::RecalcVect_ASK_Seem_INTO_KGSK_Seem (double *arrASK,double *arrKGSK,int lenarrASK )
{
	 double arrPSK_CT_Seem[9] = {0.}, arrKGSK_Seem [9] = {0.} ;
	 RecalcVect_ASK_Seem_INTO_PSK_CT_Seem (arrASK,arrPSK_CT_Seem,lenarrASK );

	 RecalcVect_PSK_CT_Seem_INTO_KGSK_Seem (arrPSK_CT_Seem,arrKGSK_Seem,lenarrASK );

	 memcpy( arrKGSK, arrKGSK_Seem, lenarrASK * sizeof(double)) ;

	 return ;
}

// пересчет векотра  из кажущейся ПСК-ЦТ в кажушуюся КГСК
void  __fastcall TFight::RecalcVect_PSK_CT_Seem_INTO_KGSK_Seem (double *arrPSK,double *arrKGSK,int lenarrPSK )
{
	double arrMu[5] ={0} ;
	double  matrPereh_PSK_V_KGSK[9] = {0};
	arrMu[0] = mVessel.mSins.mEstQ;
	arrMu[1] = mVessel.mSins.mEstPsi ;
	arrMu[2] = mVessel.mSins.mEstTet ;
	arrMu[3] = mVessel.mDriver.mEstBet ;
	arrMu[4] = mVessel.mDriver.mEstEps ;

	///
	// переход в КГСК
	calcMatr_PSK_v_KGSK(arrMu, matrPereh_PSK_V_KGSK) ;
	MtrxMultMatrx(matrPereh_PSK_V_KGSK,3,3, arrPSK, 1, arrKGSK) ;

	///


	if (lenarrPSK < 6)return ;
//  пересчет по скорости


	MtrxMultMatrx(matrPereh_PSK_V_KGSK, 3, 3, &arrPSK[3], 1, & arrKGSK [3]) ;
	if (lenarrPSK < 9)return ;
//  пересчет по ускорению
	MtrxMultMatrx(matrPereh_PSK_V_KGSK, 3, 3, &arrPSK[6], 1, & arrKGSK [6]) ;

	return ;
}



// пересчет векотра  из кажущейся АСК в кажушуюся ПСК-ЦТ
void  __fastcall TFight::RecalcVect_ASK_Seem_INTO_PSK_CT_Seem (double *arrASK,double *arrPSK_Seem,int lenarrASK )
{
	double arrMu[5] ={0} ;
	double matrPereh_ASK_V_PSK [9] = {0};
	double arrT0[6] = {0};
	arrMu[0] = mVessel.mSins.mEstQ;
	arrMu[1] = mVessel.mSins.mEstPsi ;
	arrMu[2] = mVessel.mSins.mEstTet ;
	arrMu[3] = mVessel.mDriver.mEstBet ;
	arrMu[4] = mVessel.mDriver.mEstEps ;
// 1. пересчет вектора положения
	// переход в ПСК-АС
	calcMatr_ASK_v_PSK( &arrMu[3], matrPereh_ASK_V_PSK) ;
	MtrxMultMatrx(matrPereh_ASK_V_PSK,3,3, arrASK,1, arrT0) ;
	///
	// переъход в ПСК-ЦТ
	MtrxSumMatrx(arrT0, mVessel.marrParral, 3, 1, arrPSK_Seem ) ;
	if (lenarrASK < 6)return ;

// 2. пересчет по скорости
	MtrxMultMatrx(matrPereh_ASK_V_PSK,3,3, &arrASK[3],1, &arrPSK_Seem[3]) ;
	if (lenarrASK < 9)return ;

// 3. пересчет по ускорению
	MtrxMultMatrx(matrPereh_ASK_V_PSK,3,3, &arrASK[6], 1, &arrPSK_Seem[6]) ;

	return ;

}

// пересчет векотра  из кажущейся АСК в кажушуюся ПСК-АС
void  __fastcall TFight::RecalcVect_ASK_Seem_INTO_PSK_AS_Seem (double *arrASK,double *arrPSK_Seem,int lenarrASK )
{

	double arrMu[5] ={0} ;
	double matrPereh_ASK_V_PSK [9] = {0};
	arrMu[0] = mVessel.mSins.mEstQ;
	arrMu[1] = mVessel.mSins.mEstPsi ;
	arrMu[2] = mVessel.mSins.mEstTet ;
	arrMu[3] = mVessel.mDriver.mEstBet ;
	arrMu[4] = mVessel.mDriver.mEstEps ;
	// 1. переход в ПСК-АС
	calcMatr_ASK_v_PSK( &arrMu[3], matrPereh_ASK_V_PSK) ;
	MtrxMultMatrx(matrPereh_ASK_V_PSK,3,3, arrASK,1, arrPSK_Seem) ;
	///

	if (lenarrASK < 6)return ;

// 2. пересчет по скорости
	MtrxMultMatrx(matrPereh_ASK_V_PSK,3,3, &arrASK[3],1, &arrPSK_Seem[3]) ;
	if (lenarrASK < 9)return ;

// 3. пересчет по ускорению
	MtrxMultMatrx(matrPereh_ASK_V_PSK,3,3, &arrASK[6], 1, &arrPSK_Seem[6]) ;

	return ;

}

 /*
// пересчет замера   из кажущейся АСфСК в кажущуюся КГСК
void  __fastcall TFight::RecalcZamer_INTO_KGSK (const double valVZv
		  ,const double valRZv, const double valUZv, double *valX_KGSKZv
		  ,double *valY_KGSKZv,double *valZ_KGSKZv)
{
	double arrMu[5] ={0} ;
	double matrPereh_ASK_V_PSK [9] = {0}, matrPereh_PSK_V_KGSK[9] = {0};
	double arrT0[6] = {0},arrT1[6] = {0},arrT2[6] = {0} ;
	double arrASK [3] ={0.} ;
	arrASK [0] =  valRZv * valVZv ;
	arrASK [1] =  valRZv ;
	arrASK [2] =  valRZv * valUZv ;
	arrMu[0] = mVessel.mSins.mEstQ;
	arrMu[1] = mVessel.mSins.mEstPsi ;
	arrMu[2] = mVessel.mSins.mEstTet ;
	arrMu[3] = mVessel.mDriver.mEstBet ;
	arrMu[4] = mVessel.mDriver.mEstEps ;
// 1. пересчет вектора положения
	// переход в ПСК-АС
	calcMatr_ASK_v_PSK( &arrMu[3], matrPereh_ASK_V_PSK) ;
	MtrxMultMatrx(matrPereh_ASK_V_PSK,3,3, arrASK,1, arrT0) ;
	///
	// переъход в ПСК-ЦТ
	MtrxSumMatrx(arrT0, mVessel.marrParral, 3, 1, arrT1 ) ;

	///
	// переход в КГСК
	calcMatr_PSK_v_KGSK(arrMu, matrPereh_PSK_V_KGSK) ;
	MtrxMultMatrx(matrPereh_PSK_V_KGSK,3,3, arrT1, 1, arrT2) ;

	///
	*valX_KGSKZv = arrT2[0] ;
	*valY_KGSKZv = arrT2[1] ;
	*valZ_KGSKZv = arrT2[2] ;

}
	*/

	/*
// расчет корреляционной матрицы ошибок флуктуацинного происх и ошибок привода
// расчет корреляционной матрицы ошибок оценивания в ГСК
// учитываются флуктуационные составляющие и составляющие,
// вызванные ошибками определения углов ориентации привода BET и Eps
// ош ибки определения углов ориентации корабля не учитываются
 void  __fastcall TFight::RealcCorMtrx_Bet_Eps_VRU_ASK_Seem_INTO_GSK (double *arrCorrMtrx_GSK )
{
   double ar_dF_po_dBet_sq [9],ar_dF_po_dEps_sq [9], matrPereh_ASK_V_KGSK[9],arrMu[5] ;
	   double  arrDT0[9] = {0}, arrDT1[9] = {0}, arrDT2[9] = {0}
				,arrDT3[9] = {0},arrDT4[9] = {0},arrDT5[9] = {0};
	arrMu[0] = mVessel.mSins.mEstQ;
	arrMu[1] = mVessel.mSins.mEstPsi ;
	arrMu[2] = mVessel.mSins.mEstTet ;
	arrMu[3] = mVessel.mDriver.mEstBet ;
	arrMu[4] = mVessel.mDriver.mEstEps ;
	   calc_dF_po_dBet_sq(mVessel.mDriver.mRealBet,  mVessel.mDriver.mRealEps,ar_dF_po_dBet_sq ) ;
	   calc_dF_po_dEps_sq(mVessel.mDriver.mRealBet,  mVessel.mDriver.mRealEps,ar_dF_po_dEps_sq ) ;
	   calcMatr_ASK_v_KGSK(arrMu,matrPereh_ASK_V_KGSK) ;
	   arrDT0[0] = mVessel.mMeasurer.mSigV * mVessel.mMeasurer.mSigV * mVessel.mMeasurer.mRzv * mVessel.mMeasurer.mRzv ;
	   arrDT0[4] = mVessel.mMeasurer.mSigR * mVessel.mMeasurer.mSigR ;
	   arrDT0[8] = mVessel.mMeasurer.mSigU * mVessel.mMeasurer.mSigU * mVessel.mMeasurer.mRzv * mVessel.mMeasurer.mRzv ;
	   MtrxMultMatrx(matrPereh_ASK_V_KGSK,3, 3, arrDT0,3, arrDT1) ;
	   MtrxMultMatrxTransp(arrDT1,3, 3, matrPereh_ASK_V_KGSK, 3, arrDT2) ;

	   double temp0 =  mVessel.mMeasurer.mRzv * mVessel.mMeasurer.mRzv
	   * mVessel.mDriver.mSigBet* mVessel.mDriver.mSigBet ;
	   MatrxMultScalar(ar_dF_po_dBet_sq , 3, 3, temp0,arrDT3);
	   MtrxSumMatrx(arrDT2, arrDT3,3,3, arrDT4) ;

	   double temp1 =  mVessel.mMeasurer.mRzv * mVessel.mMeasurer.mRzv
	   * mVessel.mDriver.mSigEps* mVessel.mDriver.mSigEps ;
	   MatrxMultScalar(ar_dF_po_dEps_sq , 3, 3, temp1,arrDT5);
	   MtrxSumMatrx(arrDT4, arrDT5,3,3, arrCorrMtrx_GSK) ;
 }
 */
 /*
 // расчет результирующей корреляционной матрицы ошибок, включающней ошибки
 // флуктуацинного происх ,  ошибки привода и ошибки определония качек в ГСК
// все ошибки складываютя под корнем
 void  __fastcall TFight::SummarizedCorMtrx_ErrMes_In_GSK (double *arrCorrMtrx_GSK )
{
   double ar_dF_po_dBet_sq [9],ar_dF_po_dEps_sq [9], matrPereh_ASK_V_KGSK[9],arrMu[5] ;
	   double  arrDT0[9] = {0}, arrDT1[9] = {0}, arrDT2[9] = {0}
				,arrDT3[9] = {0},arrDT4[9] = {0},arrDT5[9] = {0},arrDT51[9] = {0};
	arrMu[0] = mVessel.mSins.mEstQ;
	arrMu[1] = mVessel.mSins.mEstPsi ;
	arrMu[2] = mVessel.mSins.mEstTet ;
	arrMu[3] = mVessel.mDriver.mEstBet ;
	arrMu[4] = mVessel.mDriver.mEstEps ;
	   calc_dF_po_dBet_sq(mVessel.mDriver.mRealBet,  mVessel.mDriver.mRealEps,ar_dF_po_dBet_sq ) ;
	   calc_dF_po_dEps_sq(mVessel.mDriver.mRealBet,  mVessel.mDriver.mRealEps,ar_dF_po_dEps_sq ) ;
	   calcMatr_ASK_v_KGSK(arrMu,matrPereh_ASK_V_KGSK) ;
	   arrDT0[0] = mVessel.mMeasurer.mSigV * mVessel.mMeasurer.mSigV * mVessel.mMeasurer.mRzv * mVessel.mMeasurer.mRzv ;
	   arrDT0[4] = mVessel.mMeasurer.mSigR * mVessel.mMeasurer.mSigR ;
	   arrDT0[8] = mVessel.mMeasurer.mSigU * mVessel.mMeasurer.mSigU * mVessel.mMeasurer.mRzv * mVessel.mMeasurer.mRzv ;
	   MtrxMultMatrx(matrPereh_ASK_V_KGSK,3, 3, arrDT0,3, arrDT1) ;
	   MtrxMultMatrxTransp(arrDT1,3, 3, matrPereh_ASK_V_KGSK, 3, arrDT2) ;

	   double temp0 =  mVessel.mMeasurer.mRzv * mVessel.mMeasurer.mRzv
	   * mVessel.mDriver.mSigBet* mVessel.mDriver.mSigBet ;
	   MatrxMultScalar(ar_dF_po_dBet_sq , 3, 3, temp0,arrDT3);
	   MtrxSumMatrx(arrDT2, arrDT3,3,3, arrDT4) ;

	   double temp1 =  mVessel.mMeasurer.mRzv * mVessel.mMeasurer.mRzv
	   * mVessel.mDriver.mSigEps* mVessel.mDriver.mSigEps ;
	   MatrxMultScalar(ar_dF_po_dEps_sq , 3, 3, temp1,arrDT5);
	   MtrxSumMatrx(arrDT4, arrDT5,3,3, arrDT51) ;

	   // расчет матриц, вызванных ошибками определения качек
	   double ar_dF_po_dTet_sq [9] = {0},ar_dF_po_dQ_sq [9] = {0},ar_dF_po_dPsi_sq [9] = {0}, arrDT6[9] = {0}
		  , arrDT7[9]  = {0}, arrDT8[9] = {0}, arrDT9[9] = {0}, arrDT10[9]  = {0};
	   calc_dF_po_dQ_sq  (mVessel.mDriver.mEstBet,  mVessel.mDriver.mEstEps, ar_dF_po_dQ_sq) ;
	   calc_dF_po_dPsi_sq(mVessel.mDriver.mEstBet,  mVessel.mDriver.mEstEps, ar_dF_po_dPsi_sq) ;
	   calc_dF_po_dTet_sq(mVessel.mDriver.mEstBet,  mVessel.mDriver.mEstEps, ar_dF_po_dTet_sq) ;

	   double temp2 = mVessel.mMeasurer.mRzv * mVessel.mMeasurer.mRzv
	   * mVessel.mSins.mMaxSig_Q * mVessel.mSins.mMaxSig_Q  ;
	   double temp3 = mVessel.mMeasurer.mRzv * mVessel.mMeasurer.mRzv
	   * mVessel.mSins.mMaxSig_Psi* mVessel.mSins.mMaxSig_Psi ;
	   double temp4 = mVessel.mMeasurer.mRzv * mVessel.mMeasurer.mRzv
	   * mVessel.mSins.mMaxSig_Tet* mVessel.mSins.mMaxSig_Tet  ;

	   MatrxMultScalar(ar_dF_po_dQ_sq   , 3, 3, temp2,arrDT6);
	   MatrxMultScalar(ar_dF_po_dPsi_sq , 3, 3, temp3,arrDT7);
	   MatrxMultScalar(ar_dF_po_dTet_sq , 3, 3, temp4,arrDT8);

	   MtrxSumMatrx(arrDT51, arrDT6,3,3, arrDT9) ;
	   MtrxSumMatrx(arrDT9, arrDT7,3,3, arrDT10) ;
	   MtrxSumMatrx(arrDT10, arrDT8,3,3, arrCorrMtrx_GSK) ;
 }
  */

// расчет вектора f
void  __fastcall TFight::calcVect_f(const double h, const double valQ
	   , const double valV, const double valH, double *arr_f)
{
   memset(arr_f, 0, 6 * sizeof(double)) ;
   arr_f[0] = - h *valV * sin(valQ);
   arr_f[1] = - h *valV * cos(valQ);
   arr_f[2] = -valH ;

}

void  __fastcall TFight::CreateMatrL6 (const double h, double *arrL)
{
	for (int i = 0; i < 3; i++)
	{
	   arrL[i       * 6 +i    ] = 1;
	   arrL[i       * 6 +i + 3] = h;
	   arrL[(i + 3) * 6 +i + 3] = 1;

	}
}


// расчет вектора истинного положения в антенной сферической системе координат
void __fastcall TFight::SimulatePositionASphSK(double &valR,double &valV,double &valU)
{
   // расчет вектора положения цели в КГСК
  double arrPosKGSK[3] ={0.};
  for (int i = 0; i < 3; i++) arrPosKGSK[i] = mTarget.mTraject.marrVectSostGSK[i] - mVessel.marrVectSost[i] ;
  //

  // персчет вектора положения в АСК
	double arrPosASK[3] = {0} ;
  RecalcVect_KGSK_INTO_ASK_True (arrPosKGSK,arrPosASK,3 );
  //  пресчет вектора положения цели из линейной АСК в АСфСК

	Recalc_ASK_INTO_ASphSK(arrPosASK,valR,valV,valU) ;




}
// расчет вектора истинного положения в антенной сферичксекой системе координат
// arrPosASK - вектор положения в ASK
void __fastcall TFight::Recalc_ASK_INTO_ASphSK(double *arrPosASK, double &valR,double &valV,double &valU)
{
	valR = sqrt( arrPosASK[0] * arrPosASK[0] + arrPosASK[1] * arrPosASK[1] + arrPosASK[2] * arrPosASK[2]) ;
	double temp =  arrPosASK[0] / valR;
	if (fabs (temp) > 0.9999999999)
	{
	temp =  0.9999999999 * SIGN_D(temp);
	}
 valV  = asin(temp ) ;

	temp =  arrPosASK[2] / valR;
	if (fabs (temp) > 0.9999999999)
	{
	temp =  0.9999999999 * SIGN_D(temp);
	}

 valU  = asin( temp ) ;
}
// расчет положения цели в КГСК
void __fastcall TFight::SimulatePositionKGSK(double *arrPosKGSK)
{
  for (int i = 0; i < 3; i++) arrPosKGSK[i] = mTarget.mTraject.marrVectSostGSK[i] - mVessel.marrVectSost[i] ;

}

/*
// шаг фильтрации сопровождения цели
// INPUT:
// valTExtrap - время расчета целеуказаний
// OUTPUT:
// valTargDesBet,valTargDesEps ,valTargDesR  - целеуказания в сферической ПСК намоммент  valTExtrap
// SigTargDesBet, SigTargDesEps, SigTargDesR  - СКЗ ошибок целеуказаний
// arrVSExtrap_GSK - экстраполированный вектор состояния на момент valTExtrap, персчитанный в ГСК
// arrSigExtrapPosit_GSK[0] - СКЗ ошибки оцениывания экстраполированного вектора состояния на момент valTExtrap по координате X в ГСК
// arrSigExtrapPosit_GSK[1] - СКЗ ошибки оцениывания экстраполированного вектора состояния на момент valTExtrap по координате Y в ГСК
// arrSigExtrapPosit_GSK[2] - СКЗ ошибки оцениывания экстраполированного вектора состояния на момент valTExtrap по координате Z (H) в ГСК
// SigModulV_GSK - СКЗ ошибки оценивания экстраполированного значения модуля скорости на момент valTExtrap
void __fastcall TFight::TrackFiltration(TZamer InpASKZamer, const double valTExtrap,const int TypeOfFilter
		 ,double &valTargDesBet,double &valTargDesEps ,double &valTargDesR
		 , double & SigTargDesBet,double & SigTargDesEps,double & SigTargDesR
		 , double *arrVSExtrap_GSK, double *arrSigExtrapPosit_GSK, double &SigModulV_GSK )
{
	try
	{
	const double valT = InpASKZamer.mT ;// время привязки замера
		 // 1. Экстраполчяция оценок
	double arrMu[5] = {0.} ; // вектор обновленных на момент  valT углов

	arrMu[0] = mVessel.mSins.mEstQ;
	arrMu[1] = mVessel.mSins.mEstPsi ;
	arrMu[2] = mVessel.mSins.mEstTet ;
	arrMu[3] = mVessel.mDriver.mEstBet ;
	arrMu[4] = mVessel.mDriver.mEstEps ;
		 // 1. расчет матриц переходла на момент valT
		 double arrF_ASK_INTO_KGSK[9] = {0}
			   ,arrFBig_ASK_INTO_KGSK[36] = {0}
			   ,arrF_ASK_INTO_PSK[9]= {0}
			   ,arrFBig_ASK_INTO_PSK[36] = {0}
			   ,arrF_PSK_INTO_KGSK[9] ={0}
			   ,arrFBig_PSK_INTO_KGSK [36] = {0} ;

		calcMatr_ASK_v_PSK(&arrMu[3],arrF_ASK_INTO_PSK) ;
		calcMatr_PSK_v_KGSK(arrMu,arrF_PSK_INTO_KGSK) ;
		MtrxMultMatrx(arrF_PSK_INTO_KGSK,3, 3, arrF_ASK_INTO_PSK,3, arrF_ASK_INTO_KGSK) ;

		createExtendMtrx(arrF_ASK_INTO_KGSK, arrFBig_ASK_INTO_KGSK);
		createExtendMtrx(arrF_ASK_INTO_PSK, arrFBig_ASK_INTO_PSK);
		createExtendMtrx(arrF_PSK_INTO_KGSK, arrFBig_PSK_INTO_KGSK);
		///
		// 2. расчет экстраполированного вектора состояния цели на момент valT
		 double valSigW = 0 ;  //  скз шума в движении цели
		 const double h = valT - mTraceFlt.mTf;
		// расчет FT(t + h)  * F(t)

		double arr1[36] = {0},arr2[36] = {0},arr3[9] = {0}
		,arr4[9] = {0},arr5[6] = {0.},arr6[6] = {0.},arr7[6] = {0.};
		   MtrxTranspMultMatrx(arrFBig_ASK_INTO_KGSK,6, 6, mTraceFlt.marrFExt_ASK_V_KGSK ,6, arr1 ) ;// изм. 22ю09 !!!!
		///

	   // расчет FT(t + h) *  FT(t) * S(t)

		MtrxMultMatrx(arr1,6, 6, mTraceFlt.marrVS_ASK  ,1, arr2);
		///
		 // экстраполяция
		 for (int i =0; i < 3; i++)  arr2[i ] += h * arr2[  i + 3];
		  ///

		//расчет слагаемого с вектором параллакса
		MtrxMinusMatrx(mTraceFlt.marrF_PSK_V_KGSK  , arrF_PSK_INTO_KGSK,3, 3, arr3) ;
		MtrxTranspMultMatrx(arrF_ASK_INTO_KGSK,3, 3, arr3,3, arr4 ) ;
		MtrxMultMatrx(arr4, 3, 3, mVessel.marrParral, 1, arr5);
		///
		// расчет слагаемого с вектором  f
		double arr_f[6] = {0.} ;
		calcVect_f( h, mVessel.mSins.mEstQ, mVessel.mSins.mEstVVess, mVessel.mSins.mEstVH * h, arr_f) ;
		//	calcVect_f( h, mVessel.mSins.mEstQ, mVessel.mSins.mEstVVess, mVessel.mSins.mEstH , arr_f) ;
		MtrxTranspMultMatrx(arrF_ASK_INTO_KGSK,3, 3, arr_f,1, arr6 ) ;
		///
		// расчет экстраполированного вектора состояния цели на момент valT
		MtrxSumMatrx(arr2, arr5,6, 1, arr7) ;
		MtrxSumMatrx(arr7, arr6,6, 1, mTraceFlt.marrVS_ASK) ;

		// 3. трасформация коррел матрицы ошибок оценивания из АСК(t) в АСК(valT)
		double arrK[36] = {0.},arrKT0[36] = {0.} ;
		MtrxMultMatrx(arr1,6, 6, mTraceFlt.marrK_ASK ,6, arrKT0) ;
		MtrxMultMatrxTransp(arrKT0,6, 6, arr1,6, mTraceFlt.marrK_ASK) ;

		valSigW = (TypeOfFilter == -1)?sqrt(3.1):SIGW;
		///


		///
		// 4. дезинтеграция  корреляц матрицы по осям  на 3  матрицы 2х2
			 double arrKAxes[12] ={0.} ;
			 DesintegrateMatrK(mTraceFlt.marrK_ASK ,arrKAxes ) ;
		 ///
		 // 5. создание массивов параметров фильтрации по каждой из осей координат
		   double  arrSigBMO[3] = {0.} // массив скз БМО измерений
				  ,arrSigMMO[3] = {0.}  // массив скз ММО измерений
				  ,arrrSig_Add[3] ={0.} ; // массив дополнительного шума в канале положения цели
		arrSigBMO [0] = InpASKZamer.marrMeas[1] * sqrt( mVessel.mDriver.mSigBet * mVessel.mDriver.mSigBet * cos( mVessel.mDriver.mEstEps)* cos( mVessel.mDriver.mEstEps)
				 * (1. - VALKMMO) + InpASKZamer.marrCorr[0]) ;
		arrSigBMO [1] = sqrt(InpASKZamer.marrCorr[4]);
		arrSigBMO [2] = InpASKZamer.marrMeas[1] * sqrt( mVessel.mDriver.mSigEps * mVessel.mDriver.mSigEps * (1. - VALKMMO)
						 + InpASKZamer.marrCorr[9]) ;
						 //
		arrSigMMO [0] = InpASKZamer.marrMeas[1] * cos( mVessel.mDriver.mEstEps) * mVessel.mDriver.mSigBet * sqrt ( VALKMMO ) ;
		arrSigMMO [1] = 0 ;
		arrSigMMO [2] = InpASKZamer.marrMeas[1] * mVessel.mDriver.mSigEps * sqrt(VALKMMO ) ;
					  //

					   //
	   // расчет добавки  в правую часть уравн движения
		arrrSig_Add [0] = InpASKZamer.marrMeas[1] * sqrt ( cos( mVessel.mDriver.mEstEps) * cos( mVessel.mDriver.mEstEps) * mVessel.mSins.mSig_dQdt* mVessel.mSins.mSig_dQdt
						 + sin(mVessel.mDriver.mEstBet) * sin(mVessel.mDriver.mEstBet)*sin( mVessel.mDriver.mEstEps)*sin( mVessel.mDriver.mEstEps)* mVessel.mSins.mSig_dPsidt*mVessel.mSins.mSig_dPsidt
						 + cos(mVessel.mDriver.mEstBet) * cos(mVessel.mDriver.mEstBet)*sin( mVessel.mDriver.mEstEps)*sin( mVessel.mDriver.mEstEps)* mVessel.mSins.mSig_dTetdt* mVessel.mSins.mSig_dTetdt)
						 * sqrt(h) ;
		arrrSig_Add [1] =  0.;
		arrrSig_Add [2] = InpASKZamer.marrMeas[1] * sqrt ( cos(mVessel.mDriver.mEstBet) * cos(mVessel.mDriver.mEstBet)* mVessel.mSins.mSig_dPsidt*mVessel.mSins.mSig_dPsidt
						  + sin(mVessel.mDriver.mEstBet) * sin(mVessel.mDriver.mEstBet)* mVessel.mSins.mSig_dTetdt* mVessel.mSins.mSig_dTetdt)
						  * sqrt(h) ;
		 ///
		 // 6. Создание вектора замера в линейной АСК (valT)
		 double arrYZv_ASK[3] = {0.} ;
		 arrYZv_ASK[0] = InpASKZamer.marrMeas[0] * InpASKZamer.marrMeas[1];
		 arrYZv_ASK[1] = InpASKZamer.marrMeas[1] ;
		 arrYZv_ASK[2] = InpASKZamer.marrMeas[2] * InpASKZamer.marrMeas[1];


		 // 7. фильтрация  по осям и экстраполяция матрицы в АСК(valT) на момент valTExtrap
		double arrKExtrAxes[12]= {0.},arrKExtr[36]= {0.}, arrP[2],arrKOut[4] = {0.} ;
		switch(TypeOfFilter)
		{
		case -1:

			for (int i = 0; i < 3; i++)
			{
			TStabSyst2 Sys(h, valSigW , arrSigMMO [i], arrSigBMO [i ],arrrSig_Add[i]);

			Sys.OneStepFltrKlm_Real(&arrKAxes[ i * 4],arrKOut,arrP);
			memcpy(&arrKAxes[ i * 4], arrKOut, 4 * sizeof(double)) ;

			mTraceFlt.marrVS_ASK[     i ] +=   arrP [ 0] * (arrYZv_ASK[ i ] - mTraceFlt.marrVS_ASK[ i ]) ;
			mTraceFlt.marrVS_ASK[  i + 3] +=   arrP [ 1] * (arrYZv_ASK[ i ] - mTraceFlt.marrVS_ASK[ i ]) ;

			Sys.ExtrapolStepFltrKlm_Real(valTExtrap -  valT, &arrKAxes [i * 4], &arrKExtrAxes [ i * 4] ) ;

			}
		break;

		case 1:
			for (int i = 0; i < 3; i++)
			{
			TStabSyst2 Sys(h, valSigW , arrSigMMO [i], arrSigBMO [i ],arrrSig_Add[i]);

			Sys.OneStepGolubev_flAdd(&arrKAxes[ i * 4],arrP);

			double del =  (arrYZv_ASK[ i    ] - mTraceFlt.marrVS_ASK[ i ]) ;
			mTraceFlt.marrVS_ASK [ i    ] +=   arrP [ 0] * del ;
			mTraceFlt.marrVS_ASK [ i + 3] +=   arrP [ 1] * del ;
			Sys.ExtrapolStepGolubev_flAdd(valTExtrap -  valT, &arrKAxes [i * 4], &arrKExtrAxes [ i * 4] ) ;
			}
		break;
		default:
		ShowMessage(L"TYPE_OF_FILTER was appointed mistakenly" );
		break;
		}

		// 8. интеграция обединенной корреляционной матрицы 6х6  из 3-х матриц 2х2 по осям координат

		 IntegrateMatrK(arrKAxes , mTraceFlt.marrK_ASK ) ;
		 IntegrateMatrK(arrKExtrAxes , arrKExtr ) ;

		 // 9. запоминание рекурсивной информации
		 mTraceFlt.mTf = valT ;
		 memcpy(mTraceFlt.marrF_ASK_V_PSK, arrF_ASK_INTO_PSK, 9 * sizeof(double)) ;
		 memcpy(mTraceFlt.marrFExt_ASK_V_PSK, arrFBig_ASK_INTO_PSK, 36 * sizeof(double)) ;
		 memcpy(mTraceFlt.marrF_PSK_V_KGSK, arrF_PSK_INTO_KGSK, 9 * sizeof(double)) ;
		 memcpy(mTraceFlt.marrFExt_PSK_V_KGSK, arrFBig_PSK_INTO_KGSK, 36 * sizeof(double)) ;
		 memcpy(mTraceFlt.marrF_ASK_V_KGSK, arrF_ASK_INTO_KGSK, 9 * sizeof(double)) ;
		 memcpy(mTraceFlt.marrFExt_ASK_V_KGSK, arrFBig_ASK_INTO_KGSK, 36 * sizeof(double)) ;

		 memcpy(mTraceFlt.marrMu, arrMu, 5 * sizeof(double)) ;

			// пересчет оценки вектора состояния в КГСК
			// MtrxMultMatrx(mTraceFlt.marrFExt_ASK_V_KGSK,6, 6, mTraceFlt.marrVS_ASK,1, mTraceFlt.marrVS_KGSK) ;
			 RecalcVect_ASK_Seem_INTO_KGSK_Seem (mTraceFlt.marrVS_ASK , mTraceFlt.marrVS_KGSK,6 )  ;
			 for (int i = 0; i < 3; i++) mTraceFlt.marrVS_KGSK[  3 + i] -=  mVessel.marrEstVectSost [ 3 + i] ;



			 // пересчет коррел матрицы ошибок из АСК в КГСК
			 double   arrTempr0[36] = {0.} ;;
			 MtrxMultMatrx(mTraceFlt.marrFExt_ASK_V_KGSK,6, 6, mTraceFlt.marrK_ASK,6,  arrTempr0) ;
			 MtrxMultMatrxTransp(arrTempr0,6, 6, mTraceFlt.marrFExt_ASK_V_KGSK,6, mTraceFlt.marrK_KGSK ) ;
		 // 10. расчет целеуказаний и вектора состояния в ПСК  на момент  valTExtrap
		 double arrVSExtrap_PSK[6] = {0.} ;
		 getTargDes_PSK_True(valTExtrap, valTargDesBet, valTargDesEps, valTargDesR, arrVSExtrap_PSK);
		 ///
		 // 11. расчет точности экстраполяции
			 // 11.1 расчет точности целеуказаний в ПСК(И)
				 // 11.1.1 пересчет маторцы  arrKExtr из АСК в ПСК
			 double arrTT0[36] = {0}, arrKExtr_PSK[36] = {0} ;
			 MtrxMultMatrx(mTraceFlt.marrFExt_ASK_V_PSK,6, 6, arrKExtr,6, arrTT0 ) ;

			 MtrxMultMatrxTransp(arrTT0 ,6, 6,mTraceFlt.marrFExt_ASK_V_PSK,6, arrKExtr_PSK) ;

			 //11.1.2 расчет точности экстраполяции в сферич сист коорд
			  // создание коррел матрицы ошибок оценивания в ПСК экстраполированного вектора положения цели
			  double arrKPosExtrap_PSK [9] ={0} ;
			  for (int i = 0; i < 3;  i++)
			  for (int j = 0; j  < 3; j++)
			  {
				arrKPosExtrap_PSK [3 * i + j]   =  arrKExtr_PSK [ 6 * i + j ] ;
			  }
			  CalcToleranceSpherical(arrVSExtrap_PSK , arrKPosExtrap_PSK, SigTargDesBet, SigTargDesEps, SigTargDesR  ) ;


			 // 11.2 расчет точности оценивания модуля скорости(в любой СК)
			 SigModulV_GSK = 0 ;
			 for (int i = 0; i < 3; i++)
			 {
			  SigModulV_GSK +=  mTraceFlt.marrVS_ASK [3 +i] * mTraceFlt.marrVS_ASK [3 + i] * arrKExtr [ (3 + i)* 6 + 3 + i];
			 }

			 double valVTemp =  Norm3(& mTraceFlt.marrVS_ASK[3]) ;
			 if( valVTemp > 100.)
			 {
			 SigModulV_GSK = sqrt(SigModulV_GSK)/ valVTemp ;
			 }
			 else
			 {
			   SigModulV_GSK = sqrt(arrKExtr [ (3 + 0)* 6 + 3 + 0] + arrKExtr [ (3 + 1)* 6 + 3 + 1] + arrKExtr [ (3 + 2)* 6 + 3 + 2]) ;
             }
			 // 11.3 персчет вектора целеуказаний в ГСК
			 RecalcVect_PSK_AS_INTO_GSK (arrVSExtrap_PSK, arrVSExtrap_GSK, 6 );

			 // 11.4 персчет корреляц матрицы ошибок экстраполяции  в ГСК (== КГСК)
			 double arrT0[36] = {0.},arrT1[36] = {0.};
			 MtrxMultMatrx(mTraceFlt.marrF_ASK_V_KGSK,6, 6, arrKExtr,6, arrT0) ;
			 MtrxMultMatrxTransp(arrT0,6, 6, mTraceFlt.marrF_ASK_V_KGSK,6, arrT1) ;
			 // 11.5 расчет СКЗ ошибок экстраполяции в ГСК (= КГСК)
			 arrSigExtrapPosit_GSK[0] =  sqrt(arrT1[0]) ;
			 arrSigExtrapPosit_GSK[1] =  sqrt(arrT1[6 + 1]) ;
			 arrSigExtrapPosit_GSK[2] =  sqrt(arrT1[6 * 2 + 2]) ;
  }
  catch(...)
  {
		ShowMessage(L"Error in TFight::TrackFiltration");
  }

}

*/


// шаг фильтрации сопровождения цели
// INPUT:
// valTExtrap - время расчета целеуказаний
// OUTPUT:
// valTargDesBet,valTargDesEps ,valTargDesR  - целеуказания в сферической ПСК намоммент  valTExtrap
// SigTargDesBet, SigTargDesEps, SigTargDesR  - СКЗ ошибок целеуказаний
// arrVSExtrap_GSK - экстраполированный вектор состояния на момент valTExtrap, персчитанный в ГСК
// arrSigExtrapPosit_GSK[0] - СКЗ ошибки оцениывания экстраполированного вектора состояния на момент valTExtrap по координате X в ГСК
// arrSigExtrapPosit_GSK[1] - СКЗ ошибки оцениывания экстраполированного вектора состояния на момент valTExtrap по координате Y в ГСК
// arrSigExtrapPosit_GSK[2] - СКЗ ошибки оцениывания экстраполированного вектора состояния на момент valTExtrap по координате Z (H) в ГСК
// SigModulV_GSK - СКЗ ошибки оценивания экстраполированного значения модуля скорости на момент valTExtrap
void __fastcall TFight::TrackFiltration(TZamer InpASKZamer, const double valTExtrap ,const int TypeOfFilter
	,double &valTargDesBet,double &valTargDesEps ,double &valTargDesR
	, double & SigTargDesBet,double & SigTargDesEps,double & SigTargDesR
	, double *arrVSExtrap_GSK, double *arrSigExtrapPosit_GSK, double &SigModulV_GSK )
{
	try
	{
	const double valT = InpASKZamer.mT ;// время привязки замера
		 // 1. Экстраполчяция оценок
	double arrMu[5] = {0.} ; // вектор обновленных на момент  valT углов

	arrMu[0] = mVessel.mSins.mEstQ;
	arrMu[1] = mVessel.mSins.mEstPsi ;
	arrMu[2] = mVessel.mSins.mEstTet ;
	arrMu[3] = mVessel.mDriver.mEstBet ;
	arrMu[4] = mVessel.mDriver.mEstEps ;
		 // 1. расчет матриц переходла на момент valT
		 double arrF_ASK_INTO_KGSK[9] = {0}
			   ,arrFBig_ASK_INTO_KGSK[36] = {0}
				 ,arrF_ASK_INTO_PSK[9]= {0}
			   ,arrFBig_ASK_INTO_PSK[36] = {0}
			   ,arrF_PSK_INTO_KGSK[9] ={0}
			   ,arrFBig_PSK_INTO_KGSK [36] = {0} ;

		calcMatr_ASK_v_PSK(&arrMu[3],arrF_ASK_INTO_PSK) ;
		calcMatr_PSK_v_KGSK(arrMu,arrF_PSK_INTO_KGSK) ;
		MtrxMultMatrx(arrF_PSK_INTO_KGSK,3, 3, arrF_ASK_INTO_PSK,3, arrF_ASK_INTO_KGSK) ;

		createExtendMtrx(arrF_ASK_INTO_KGSK, arrFBig_ASK_INTO_KGSK);
		createExtendMtrx(arrF_ASK_INTO_PSK, arrFBig_ASK_INTO_PSK);
		createExtendMtrx(arrF_PSK_INTO_KGSK, arrFBig_PSK_INTO_KGSK);
		///
		// 2. расчет экстраполированного вектора состояния цели на момент valT
	   //	 double valSigW = 0 ;  //  скз шума в движении цели
		 const double h = valT - mTraceFlt.mTf;
		// расчет FT(t + h)  * F(t)

		double arr1[36] = {0},arr2[36] = {0},arr3[9] = {0}
		,arr4[9] = {0},arr5[6] = {0.},arr6[6] = {0.},arr7[6] = {0.};
		   MtrxTranspMultMatrx(arrFBig_ASK_INTO_KGSK,6, 6, mTraceFlt.marrFExt_ASK_V_KGSK ,6, arr1 ) ;// изм. 22ю09 !!!!
		///

	   // расчет FT(t + h) *  FT(t) * S(t)

		MtrxMultMatrx(arr1,6, 6, mTraceFlt.marrVS_ASK  ,1, arr2);
		///
		 // экстраполяция
		 for (int i =0; i < 3; i++)  arr2[i ] += h * arr2[  i + 3];
		  ///

		//расчет слагаемого с вектором параллакса
		MtrxMinusMatrx(mTraceFlt.marrF_PSK_V_KGSK  , arrF_PSK_INTO_KGSK,3, 3, arr3) ;
		MtrxTranspMultMatrx(arrF_ASK_INTO_KGSK,3, 3, arr3,3, arr4 ) ;
		MtrxMultMatrx(arr4, 3, 3, mVessel.marrParral, 1, arr5);
		///
		// расчет слагаемого с вектором  f
		double arr_f[6] = {0.} ;
		calcVect_f( h, mVessel.mSins.mEstQ, mVessel.mSins.mEstVVess, mVessel.mSins.mEstVH * h, arr_f) ;
		//	calcVect_f( h, mVessel.mSins.mEstQ, mVessel.mSins.mEstVVess, mVessel.mSins.mEstH , arr_f) ;
		MtrxTranspMultMatrx(arrF_ASK_INTO_KGSK,3, 3, arr_f,1, arr6 ) ;
		///
		// расчет экстраполированного вектора состояния цели на момент valT
		MtrxSumMatrx(arr2, arr5,6, 1, arr7) ;
		MtrxSumMatrx(arr7, arr6,6, 1, mTraceFlt.marrVS_ASK) ;

		// 3. трасформация коррел матрицы ошибок оценивания из АСК(t) в АСК(valT)
		double arrKT0[36] = {0.} ;
		MtrxMultMatrx(arr1,6, 6, mTraceFlt.marrK_ASK ,6, arrKT0) ;
		MtrxMultMatrxTransp(arrKT0,6, 6, arr1,6, mTraceFlt.marrK_ASK) ;

	  //	valSigW = (TypeOfFilter == -1)?sqrt(3.1):SIGW;
		///


		///
		// 4. дезинтеграция  корреляц матрицы по осям  на 3  матрицы 2х2
			 double arrKAxes[12] ={0.} ;
			 DesintegrateMatrK(mTraceFlt.marrK_ASK ,arrKAxes ) ;
		 ///

		 // 5. создание массивов параметров фильтрации по каждой из осей координат
		   double  arrSigBMO[3] = {0.} // массив скз БМО измерений
				  ,arrSigMMO[3] = {0.}  // массив скз ММО измерений
				  ,arrrSig_Add[3] ={0.} ; // массив дополнительного шума в канале положения цели
		arrSigBMO [0] = InpASKZamer.marrMeas[1]   * sqrt( mVessel.mDriver.mSigBet * mVessel.mDriver.mSigBet * cos( mVessel.mDriver.mEstEps)* cos( mVessel.mDriver.mEstEps)
				 * (1. - VALKMMO) + InpASKZamer.marrCorr[0]  ) ;
		arrSigBMO [1] = sqrt(InpASKZamer.marrCorr[4] );
		arrSigBMO [2] = InpASKZamer.marrMeas[1] * sqrt( mVessel.mDriver.mSigEps * mVessel.mDriver.mSigEps * (1. - VALKMMO)
						 + InpASKZamer.marrCorr[8]) ;
					   //
		arrSigMMO [0] = InpASKZamer.marrMeas[1] * cos( mVessel.mDriver.mEstEps) * mVessel.mDriver.mSigBet * sqrt ( VALKMMO ) ;
		arrSigMMO [1] = 0 ;
		arrSigMMO [2] = InpASKZamer.marrMeas[1] * mVessel.mDriver.mSigEps * sqrt(VALKMMO ) ;
					  //

					   //
	   // расчет добавки  в правую часть уравн движения
			arrrSig_Add [0] = InpASKZamer.marrMeas[1] * sqrt ( cos( mVessel.mDriver.mEstEps) * cos( mVessel.mDriver.mEstEps) * mVessel.mSins.mSig_dQdt* mVessel.mSins.mSig_dQdt
						 + sin(mVessel.mDriver.mEstBet) * sin(mVessel.mDriver.mEstBet)*sin( mVessel.mDriver.mEstEps)*sin( mVessel.mDriver.mEstEps)* mVessel.mSins.mSig_dPsidt*mVessel.mSins.mSig_dPsidt
						 + cos(mVessel.mDriver.mEstBet) * cos(mVessel.mDriver.mEstBet)*sin( mVessel.mDriver.mEstEps)*sin( mVessel.mDriver.mEstEps)* mVessel.mSins.mSig_dTetdt* mVessel.mSins.mSig_dTetdt)
						 * sqrt(h) ;
		arrrSig_Add [1] =  0.;
		arrrSig_Add [2] = InpASKZamer.marrMeas[1] * sqrt ( cos(mVessel.mDriver.mEstBet) * cos(mVessel.mDriver.mEstBet)* mVessel.mSins.mSig_dPsidt*mVessel.mSins.mSig_dPsidt
						  + sin(mVessel.mDriver.mEstBet) * sin(mVessel.mDriver.mEstBet)* mVessel.mSins.mSig_dTetdt* mVessel.mSins.mSig_dTetdt)
						  * sqrt(h) ;
		 ///
		 // 6. Создание вектора замера в линейной АСК (valT)
		 double arrYZv_ASK[3] = {0.} ;
		 arrYZv_ASK[0] = InpASKZamer.marrMeas[0] * InpASKZamer.marrMeas[1];
		 arrYZv_ASK[1] = InpASKZamer.marrMeas[1] ;
		 arrYZv_ASK[2] = InpASKZamer.marrMeas[2] * InpASKZamer.marrMeas[1];



		 // 7. фильтрация  по осям и экстраполяция матрицы в АСК(valT) на момент valTExtrap
		double arrKExtrAxes[12]= {0.},arrKExtr[36]= {0.}, arrP[2],arrKOut[4] = {0.} ;
		switch(TypeOfFilter)
		{
		case -1:

			for (int i = 0; i < 3; i++)
			{
			TStabSyst2 Sys(h, mTarget.mTraject.marrSigW[0] , arrSigMMO [i], arrSigBMO [i ],arrrSig_Add[i]);

			Sys.OneStepFltrKlm_Real(&arrKAxes[ i * 4],arrKOut,arrP);
			memcpy(&arrKAxes[ i * 4], arrKOut, 4 * sizeof(double)) ;

			mTraceFlt.marrVS_ASK[     i ] +=   arrP [ 0] * (arrYZv_ASK[ i ] - mTraceFlt.marrVS_ASK[ i ]) ;
			mTraceFlt.marrVS_ASK[  i + 3] +=   arrP [ 1] * (arrYZv_ASK[ i ] - mTraceFlt.marrVS_ASK[ i ]) ;

			Sys.ExtrapolStepFltrKlm_Real(valTExtrap -  valT, &arrKAxes [i * 4], &arrKExtrAxes [ i * 4] ) ;

			}
		break;

		case 1:
			for (int i = 0; i < 3; i++)
			{
			TStabSyst2 Sys(h, mTarget.mTraject.marrSigW[0], arrSigMMO [i], arrSigBMO [i ],arrrSig_Add[i]);

			Sys.OneStepGolubev_flAdd(&arrKAxes[ i * 4],arrP);

			double del =  (arrYZv_ASK[ i    ] - mTraceFlt.marrVS_ASK[ i ]) ;
			mTraceFlt.marrVS_ASK [ i    ] +=   arrP [ 0] * del ;
			mTraceFlt.marrVS_ASK [ i + 3] +=   arrP [ 1] * del ;
			Sys.ExtrapolStepGolubev_flAdd(valTExtrap -  valT, &arrKAxes [i * 4], &arrKExtrAxes [ i * 4] ) ;
			}
		break;
		default:
		ShowMessage(L"TYPE_OF_FILTER was appointed mistakenly" );
		break;
		}

		// 8. интеграция обединенной корреляционной матрицы 6х6  из 3-х матриц 2х2 по осям координат

		 IntegrateMatrK(arrKAxes , mTraceFlt.marrK_ASK ) ;
		 IntegrateMatrK(arrKExtrAxes , arrKExtr ) ;

		 // 9. запоминание рекурсивной информации
		 mTraceFlt.mTf = valT ;
		 memcpy(mTraceFlt.marrF_ASK_V_PSK, arrF_ASK_INTO_PSK, 9 * sizeof(double)) ;
		 memcpy(mTraceFlt.marrFExt_ASK_V_PSK, arrFBig_ASK_INTO_PSK, 36 * sizeof(double)) ;
		 memcpy(mTraceFlt.marrF_PSK_V_KGSK, arrF_PSK_INTO_KGSK, 9 * sizeof(double)) ;
		 memcpy(mTraceFlt.marrFExt_PSK_V_KGSK, arrFBig_PSK_INTO_KGSK, 36 * sizeof(double)) ;
		 memcpy(mTraceFlt.marrF_ASK_V_KGSK, arrF_ASK_INTO_KGSK, 9 * sizeof(double)) ;
		 memcpy(mTraceFlt.marrFExt_ASK_V_KGSK, arrFBig_ASK_INTO_KGSK, 36 * sizeof(double)) ;

		 memcpy(mTraceFlt.marrMu, arrMu, 5 * sizeof(double)) ;

			// пересчет оценки вектора состояния в КГСК
			// MtrxMultMatrx(mTraceFlt.marrFExt_ASK_V_KGSK,6, 6, mTraceFlt.marrVS_ASK,1, mTraceFlt.marrVS_KGSK) ;
			 RecalcVect_ASK_Seem_INTO_KGSK_Seem (mTraceFlt.marrVS_ASK , mTraceFlt.marrVS_KGSK,6 )  ;
			 for (int i = 0; i < 3; i++) mTraceFlt.marrVS_KGSK[  3 + i] -=  mVessel.marrEstVectSost [ 3 + i] ;



			 // пересчет коррел матрицы ошибок из АСК в КГСК
			 double   arrTempr0[36] = {0.} ;;
			 MtrxMultMatrx(mTraceFlt.marrFExt_ASK_V_KGSK,6, 6, mTraceFlt.marrK_ASK,6,  arrTempr0) ;
			 MtrxMultMatrxTransp(arrTempr0,6, 6, mTraceFlt.marrFExt_ASK_V_KGSK,6, mTraceFlt.marrK_KGSK ) ;


		 // 10. расчет целеуказаний и вектора состояния в ПСК  на момент  valTExtrap
		 double arrVSExtrap_PSK[6] = {0.} ;
		 getTargDes_PSK_True(valTExtrap, valTargDesBet, valTargDesEps, valTargDesR, arrVSExtrap_PSK);
		 ///
		 // 11. расчет точности экстраполяции
			 // 11.1 расчет точности целеуказаний в ПСК(И)
				 // 11.1.1 пересчет маторцы  arrKExtr из АСК в ПСК
			 double arrTT0[36] = {0}, arrKExtr_PSK[36] = {0} ;
			 MtrxMultMatrx(mTraceFlt.marrFExt_ASK_V_PSK,6, 6, arrKExtr,6, arrTT0 ) ;

			 MtrxMultMatrxTransp(arrTT0 ,6, 6,mTraceFlt.marrFExt_ASK_V_PSK,6, arrKExtr_PSK) ;

			 //11.1.2 расчет точности экстраполяции в сферич сист коорд
			  // создание коррел матрицы ошибок оценивания в ПСК экстраполированного вектора положения цели
			  double arrKPosExtrap_PSK [9] ={0} ;
			  for (int i = 0; i < 3;  i++)
			  for (int j = 0; j  < 3; j++)
			  {
				arrKPosExtrap_PSK [3 * i + j]   =  arrKExtr_PSK [ 6 * i + j ] ;
			  }
			  CalcToleranceSpherical(arrVSExtrap_PSK , arrKPosExtrap_PSK, SigTargDesBet, SigTargDesEps, SigTargDesR  ) ;


			 // 11.2 расчет точности оценивания модуля скорости(в любой СК)
			 SigModulV_GSK = 0 ;
			 for (int i = 0; i < 3; i++)
			 {
			  SigModulV_GSK +=  mTraceFlt.marrVS_ASK [3 +i] * mTraceFlt.marrVS_ASK [3 + i] * arrKExtr [ (3 + i)* 6 + 3 + i];
			 }

			 double valVTemp =  Norm3(& mTraceFlt.marrVS_ASK[3]) ;
			 if( valVTemp > 100.)
			 {
			 SigModulV_GSK = sqrt(SigModulV_GSK)/ valVTemp ;
			 }
			 else
			 {
			   SigModulV_GSK = sqrt(arrKExtr [ (3 + 0)* 6 + 3 + 0] + arrKExtr [ (3 + 1)* 6 + 3 + 1] + arrKExtr [ (3 + 2)* 6 + 3 + 2]) ;
             }
			 // 11.3 персчет вектора целеуказаний в ГСК
			 RecalcVect_PSK_AS_INTO_GSK (arrVSExtrap_PSK, arrVSExtrap_GSK, 6 );

			 // 11.4 персчет корреляц матрицы ошибок экстраполяции  в ГСК (== КГСК)
			 double arrT0[36] = {0.},arrT1[36] = {0.};
			 MtrxMultMatrx(mTraceFlt.marrF_ASK_V_KGSK,6, 6, arrKExtr,6, arrT0) ;
			 MtrxMultMatrxTransp(arrT0,6, 6, mTraceFlt.marrF_ASK_V_KGSK,6, arrT1) ;
			 // 11.5 расчет СКЗ ошибок экстраполяции в ГСК (= КГСК)
			 arrSigExtrapPosit_GSK[0] =  sqrt(arrT1[0]) ;
			 arrSigExtrapPosit_GSK[1] =  sqrt(arrT1[6 + 1]) ;
			 arrSigExtrapPosit_GSK[2] =  sqrt(arrT1[6 * 2 + 2]) ;
  }
  catch(...)
  {
		ShowMessage(L"Error in TFight::TrackFiltration");
  }

}
// нахождение точности определения координат в сферической СК
// Задан вектор arrSInp в прямоугольной СК с ошибками,
// которые характеризуются корреляц матрицей  arrKInp
// Требуется определеить СКЗ ошибок расчетва сферических коодинат
// R, Eps, Bet , где
// R*R = arrSInp[0] *  arrSInp[0] + arrSInp[1] *  arrSInp[1] + arrSInp[2] *  arrSInp[2]
// Eps = asin ( arrSInp[2] / R)
//Bet = asin ( arrSInp[0] / sqrt ( arrSInp[0] *  arrSInp[0] + arrSInp[1] *  arrSInp[1])

void  __fastcall TFight::CalcToleranceSpherical(double *arrSInp , double *arrKInp
		, double &sigB, double &sigE, double &sigR  )
{
	double valR = Norm3(arrSInp) ;
	double valr = sqrt( arrSInp[0] *  arrSInp[0] + arrSInp[1]*  arrSInp[1]) ;
  // 1. расчет точности опрделения R  = f1(x,y,z)
	//1.1 нахождение вуектора  градиента функции f1 в точке arrSInp
	double arrGradF1[3] ={0} ;
	arrGradF1 [0] = arrSInp[ 0 ] / valR ;
	arrGradF1 [1] = arrSInp[ 1 ] / valR ;
	arrGradF1 [2] = arrSInp[ 2 ] / valR ;
	// 1.2 расчет точности
	sigR = sqrt(CalcDispLinTransf ( arrGradF1, arrKInp ) );
	///

  // 2. расчет точности лпределения Bet = f2(x,y,z)
	 // 2.1 назхождение вектира градиента функции f2 вы точке  arrSInp
		double arrGradF2[3] ={0} ;
		arrGradF2 [0] =   arrSInp[1] / valr/ valr ;
		arrGradF2[ 1] = - arrSInp[ 0 ] / valr/ valr ;
		arrGradF2 [2] =   0 ;
	  // 2.2 расчет точности
		sigB = sqrt(CalcDispLinTransf ( arrGradF2, arrKInp ) );
   ///

	 // 3. расчет точности лпределения Eps = f2(x,y,z)
	 // 3.1 назхождение вектира градиента функции f3 вы точке  arrSInp
		double arrGradF3[3] ={0} ;
		arrGradF3 [0] = - arrSInp[2] *arrSInp[0] / valr/ valR / valR ;
		arrGradF3[ 1] = - arrSInp[2] *arrSInp[1] / valr/ valR / valR ;
		arrGradF3 [2] =  valr/ valR / valR ;
	  // 1.2 расчет точности
		sigE = sqrt(CalcDispLinTransf ( arrGradF3, arrKInp ) );

}
double  __fastcall TFight::CalcDispLinTransf (double * arrGradF, double *arrKInp )
{
    double sum = 0 ;
	for (int i = 0; i < 3; i++)
	for (int j = 0; j < 3; j++)
	{
	   sum +=  arrGradF[ i ] * arrGradF[ j ] *  arrKInp [ 3 * i + j ] ;
	}
	return sum ;
}
double  __fastcall TFight::Norm3(double *parr)
{
	return  (double)sqrt(parr[0]  * parr[0] + parr[1]  * parr[1] +parr[2]  * parr[2] );
}
// Разбиение корреляционной матрицы 6х6 на 3 матрицы 2х2
// INPUT: arrInp - исходная матрица [36]
// OUTPUT: arrOut[12]
// в матрице arrOut хранятся последовательно 3 матрицы 2х2
 void  __fastcall TFight::DesintegrateMatrK(double *arrInp , double *arrOut )
 {
   for (int i = 0; i < 3; i++)
   {
	 arrOut[ i * 4    ] = arrInp [ 6 * i + i          ] ;
	 arrOut[ i * 4 + 1] = arrInp [ 6 * i + i       + 3] ;
	 arrOut[ i * 4 + 2] = arrOut[ i * 4 + 1] ;
	 arrOut[ i * 4 + 3] = arrInp [ 6 * (i + 3) + i + 3] ;
   }
 }
// Создание объединенной  корреляционной матрицы 6х6 из 3 матриц 2х2
// INPUT: arrInp[12] - исходная матрица, состоящая из 3-х матирц 2х2
// OUTPUT: arrOut[36]
 void  __fastcall TFight::IntegrateMatrK(double *arrInp , double *arrOut )
 {
   memset( arrOut, 0, 36 * sizeof(double)) ;
	 for (int i = 0; i < 3; i++)
   {
	  arrOut [ 6 * i       + i    ] =  arrInp[ i * 4    ];
	  arrOut [ 6 * i       + i + 3]  = arrInp[ i * 4 + 1];
	  arrOut [ 6 * (i + 3) + i    ]  = arrInp[ i * 4 + 2];
	  arrOut [ 6 * (i + 3) + i + 3]  = arrInp[ i * 4 + 3];
   }
 }

// перечсчет вектора сотояния из KGSK в ПСК   кажущуюся
// есдли  lenarrKGSK == 6 , то пересчитывается положение и скорость
// если   lenarrKGSK == 3   , то пересчитывается только положение
// если   lenarrKGSK == 9   , тт пересчитывается только положение и скорость и ускоренике
// на вход подается вектиор состоящий из 3 , 6 или 9 координат.
// первые 3 координаты представляют из себя  положение точки в КГС
// последние 3 координаты представляют из себя скорость точки в КГСК
 void  __fastcall TFight::RecalcVect_KGSK_INTO_PSK_AS_Seem (double *arrKGSK,double *arrPSK,int lenarrKGSK )
 {
    // персчет вектора положения в АСК
	// создание вектора углов ориентации коррабля
	double arrMu[5] = {0} ;
	arrMu[0] = mVessel.mSins.mEstQ;
	arrMu[1] = mVessel.mSins.mEstPsi ;
	arrMu[2] = mVessel.mSins.mEstTet ;
	arrMu[3] = mVessel.mDriver.mEstBet ;
	arrMu[4] = mVessel.mDriver.mEstEps ;
	// создание матирицы перехода из   КГСК в ПСК
	arrMu[3] = 0 ;
	arrMu[4] = 0;
	double matrPereh_PSK_V_KGSK[9] = {0} ;
	calcMatr_PSK_v_KGSK( arrMu, matrPereh_PSK_V_KGSK) ;
	// создание вектора положения с ПСК-ЦТ
	   //	double arrPosPSK_CT[3] = {0} ;
		MtrxTranspMultMatrx(matrPereh_PSK_V_KGSK,3, 3, arrKGSK,1, arrPSK) ;
	// создание вектора положение в ПСК-АС

		for (int i =0; i < 3; i++) arrPSK[i] =  arrPSK [i] - mVessel.marrParral[i] ;

	// пересчет скорости
   if (lenarrKGSK == 6) MtrxTranspMultMatrx( matrPereh_PSK_V_KGSK,3, 3, &arrKGSK[3],1, &arrPSK[3]) ;
   if (lenarrKGSK == 9) MtrxTranspMultMatrx( matrPereh_PSK_V_KGSK,3, 3, &arrKGSK[6],1, &arrPSK[6]) ;
   return ;
 }
 // перечсчет вектора сотояния из  ПСК   в кажущуюся  АСК
// есдли  lenarrKGSK == 6 , то пересчитывается положение и скорость
// если   lenarrKGSK == 3   , то пересчитывается только положение
// если   lenarrKGSK == 9   , тт пересчитывается только положение и скорость и ускоренике
// на вход подается вектиор состоящий из 3 , 6 или 9 координат.
// первые 3 координаты представляют из себя  положение точки в КГС
// последние 3 координаты представляют из себя скорость точки в КГСК
 void  __fastcall TFight::RecalcVect_PSK_AS_Seem_INTO_ASK_Seem (double *arrPSK,double *arrASK,int lenarrPSK )
 {
    // персчет вектора положения в АСК
	// создание вектора углов ориентации коррабля
	double arrMu[2] = {0} ;
	arrMu[0] = mVessel.mDriver.mEstBet ;
	arrMu[1] = mVessel.mDriver.mEstEps ;

	// создание матрицы перехода из  АСК  в ПСК
	double  matrPereh_ASK_V_PSK[9] = {0.} ;

	calcMatr_ASK_v_PSK(arrMu, matrPereh_ASK_V_PSK) ;
	// создание вектора положения в истинной линейной АСК
	MtrxTranspMultMatrx( matrPereh_ASK_V_PSK,3, 3, arrPSK,1, arrASK) ;
	// пересчет скорости
   if (lenarrPSK == 6) MtrxTranspMultMatrx( matrPereh_ASK_V_PSK,3, 3, &arrPSK[3], 1, &arrASK[3] ) ;
   if (lenarrPSK == 9) MtrxTranspMultMatrx( matrPereh_ASK_V_PSK,3, 3, &arrPSK[6], 1, &arrASK[6]) ;
   return ;
 }

 // перечсчет вектора сотояния из KGSK в АСК  кажущуюся
// есдли  lenarrKGSK == 6 , то пересчитывается положение и скорость
// если   lenarrKGSK == 3   , то пересчитывается только положение
// если   lenarrKGSK == 9   , тт пересчитывается только положение и скорость и ускоренике
// на вход подается вектиор состоящий из 3 , 6 или 9 координат.
// первые 3 координаты представляют из себя  положение точки в КГС
// последние 3 координаты представляют из себя скорость точки в КГСК
 void  __fastcall TFight::RecalcVect_KGSK_INTO_ASK_Seem (double *arrKGSK,double *arrASK,int lenarrKGSK )
 {
	 double arrTemp[9] = {0.} ;
	 RecalcVect_KGSK_INTO_PSK_AS_Seem (arrKGSK,arrTemp, lenarrKGSK ) ;
	 RecalcVect_PSK_AS_Seem_INTO_ASK_Seem (arrTemp, arrASK, lenarrKGSK )  ;
	 return ;

 }


// занесение информации в массивы для отчета
 void TFight::updateReportData(TZamer InpASKZamer)
 {
	 if(mpwcharrFoldReport == NULL) return ;
	 if (mQuantPntReport ==  mLenMemoryAlloc)
	 {
	   mLenMemoryAlloc += 2000 ;
	   mparrBuff = (double *)realloc(mparrBuff,QUANT_COLS_BUFF * mLenMemoryAlloc * sizeof(double)) ;

	 }

	   int num0 =  mQuantPntReport * QUANT_COLS_BUFF ;

	   mparrBuff [num0] =     mT ;

	   // вектор состояния цели истинный в КГСК с 1 по 6
	   double *parrTargKGSK =  &mparrBuff [num0 + 1] ;
	   memcpy(parrTargKGSK, mTarget.mTraject.marrVectSostGSK, 6 * sizeof(double)) ;
	   for (int i = 0; i < 6; i++) parrTargKGSK [i] -=  mVessel.marrVectSost[i];

	   ///

		 // вектор состояния цели в истинной ПСК-АС с 7 по 12 позиции

		 RecalcVect_GSK_INTO_PSK_AS_True (mTarget.mTraject.marrVectSostGSK,&mparrBuff [num0 +1 + 6],6 )  ;
	   //
	   // положение цели в сферич системем коорд ПСК  c 13 по 15 поз. R, Bet, Eps
		 recalcCoord_INTO_Spherical(&mparrBuff [num0 +1 + 6], mparrBuff [num0 +13], mparrBuff [num0 +14], mparrBuff [num0 +15]) ;
	   ///

	   //  вектор состояния цели в линейной истинной АСК  c 16 по 21 поз

		 RecalcVect_KGSK_INTO_ASK_True ( &mparrBuff [num0 +1], &mparrBuff [num0 +16], 6 )  ;
	   ///
	   // положение цели в АСфСК(R,V,U)   с 22 по 24 поз
		 Recalc_ASK_INTO_ASphSK(&mparrBuff [num0 +16], mparrBuff [num0 +22],mparrBuff [num0 +23],mparrBuff [num0 +24]) ;
	   ///

	   // информация по фильтрации

	   	/////////////////////////////////////////////////////////////////////////////

			 // занесение информ в массивы по фильтрации  в АСК

	// оценка веукторc 25а сост в АСК c 25 по 30 поз
	double *parrEstASK = &mparrBuff [num0 +25]  ;
	memcpy(parrEstASK, mTraceFlt.marrVS_ASK, 6 * sizeof(double)) ;
	///
	// ошибки фильтрации в АСК  c 31 по 36 поз
  //	double *parrErrASK =  mparrBuff [num0 + 31]  ;
  //	for (int i = 0; i < 6; i++)
  //	 {
  //	  parrErrASK [i] =   parrEstASK[i] - arrTargASK[i] ;
  //	}
	// cкз ошибок фильтрации  c 31 gj 36 поз
	double *parrSigASK =  &mparrBuff [num0 + 31]  ;
   parrSigASK [0 ] =  sqrt(mTraceFlt.marrK_ASK [0]);
   parrSigASK [1 ] =  sqrt(mTraceFlt.marrK_ASK [1 * 6 + 1]);
   parrSigASK [2 ] =  sqrt(mTraceFlt.marrK_ASK [2 * 6 + 2]);
   parrSigASK [3 ] =  sqrt(mTraceFlt.marrK_ASK [3 * 6 + 3]);
   parrSigASK [4 ] =  sqrt(mTraceFlt.marrK_ASK [4 * 6 + 4]);
   parrSigASK [5 ] =  sqrt(mTraceFlt.marrK_ASK [5 * 6 + 5]);
   ///


	// оценка вектора состояния в КГСК c 37 по 42 поз
	double *parrEstKGSK = &mparrBuff [num0 + 37]  ;
	memcpy(parrEstKGSK, mTraceFlt.marrVS_KGSK, 6 * sizeof(double)) ;
	///

	// оценка векитра состояния в ПСК-АС  с 43 по 48 поз
	double *parrEstPSK = &mparrBuff [num0 + 43]  ;
	RecalcVect_ASK_Seem_INTO_PSK_AS_Seem(mTraceFlt.marrVS_ASK, parrEstPSK, 6 ) ;
	///

	//
	// целеуказания по Bet
	  mparrBuff [num0 + 49] =  mTargDesBet ;
	  // целеуказания по Eps
	  mparrBuff [num0 + 50] =  mTargDesEps ;
	  //cкз по Bet
		mparrBuff [num0 + 51] =  mSigTargDesBet ;
	  //cкз по Eps
	  mparrBuff [num0 + 52] =  mSigTargDesEps ;

	  // замер в линейной АСК
	  double *parrMeasASK =  & mparrBuff [num0 + 53] ;
		mparrBuff [num0 + 53] = InpASKZamer.marrMeas[0] * InpASKZamer.marrMeas[1] ;
		mparrBuff [num0 + 54] = InpASKZamer.marrMeas[1];
		mparrBuff [num0 + 55] = InpASKZamer.marrMeas[2] * InpASKZamer.marrMeas[1] ;
	  ///
	  // замер в ПСК-АС
	   double *parrMeasPSK_AS =  & mparrBuff [num0 + 56] ;
		 RecalcVect_ASK_Seem_INTO_PSK_AS_Seem (parrMeasASK,parrMeasPSK_AS,3 ) ;
	   ///
	   //замер в КГСК
		double *parrMeasKGSK =  & mparrBuff [num0 + 59] ;
		double arrTemp[3] = {0.} ;
		MtrxSumMatrx(parrMeasPSK_AS, mVessel.marrParral, 3, 1, arrTemp ) ;
		RecalcVect_PSK_CT_Seem_INTO_KGSK_Seem (arrTemp, parrMeasKGSK, 3 ) ;
		///
		// замер в ГСК
		double *parrMeasGSK =  & mparrBuff [num0 + 62] ;
		 MtrxSumMatrx(parrMeasKGSK, mVessel.marrEstVectSost, 3, 1, parrMeasGSK ) ;

		 // целеуказание по R
		 mparrBuff [num0 + 65] = mTargDesR ;
		 // скз по R
		 mparrBuff [num0 + 66] = mSigTargDesR ;
		mQuantPntReport++ ;

 }

// публикация  отчета
 void TFight::WriteReport()
 {
	 if (mpwcharrFoldReport == NULL )return ;




	 wchar_t *wcharrFileNames = new wchar_t[QUANT_COLS_BUFF * 30];
	 memset(wcharrFileNames, 0,  QUANT_COLS_BUFF * 30 * sizeof( wchar_t)) ;

	 wcscpy( &wcharrFileNames[ 0 * 30], L"t");
	 wcscpy( &wcharrFileNames[ 1 * 30], L"TrgKGSK_Real_X");
	 wcscpy( &wcharrFileNames[ 2* 30],  L"TrgKGSK_Real_Y");
	 wcscpy( &wcharrFileNames[ 3 * 30], L"TrgKGSK_Real_Z");
	 wcscpy( &wcharrFileNames[ 4 * 30], L"TrgKGSK_Real_VX");
	 wcscpy( &wcharrFileNames[ 5 * 30], L"TrgKGSK_Real_VY");
	 wcscpy( &wcharrFileNames[ 6 * 30], L"TrgKGSK_Real_VZ");

	 wcscpy( &wcharrFileNames[ 7 * 30], L"TrgPSK_Real_X");
	 wcscpy( &wcharrFileNames[ 8 * 30], L"TrgPSK_Real_Y");
	 wcscpy( &wcharrFileNames[ 9 * 30], L"TrgPSK_Real_Z");
	 wcscpy( &wcharrFileNames[10 * 30], L"TrgPSK_Real_VX");
	 wcscpy( &wcharrFileNames[11 * 30], L"TrgPSK_Real_VY");
	 wcscpy( &wcharrFileNames[12 * 30], L"TrgPSK_Real_VZ");

	 wcscpy( &wcharrFileNames[13 * 30], L"Trg_Real_R");
	 wcscpy( &wcharrFileNames[14 * 30], L"Trg_Real_Bet");
	 wcscpy( &wcharrFileNames[15 * 30], L"Trg_Real_Eps");

	 wcscpy( &wcharrFileNames[16 * 30], L"TrgASK_Real_X");
	 wcscpy( &wcharrFileNames[17 * 30], L"TrgASK_Real_Y");
	 wcscpy( &wcharrFileNames[18 * 30], L"TrgASK_Real_Z");
	 wcscpy( &wcharrFileNames[19 * 30], L"TrgASK_Real_VX");
	 wcscpy( &wcharrFileNames[20 * 30], L"TrgASK_Real_VY");
	 wcscpy( &wcharrFileNames[21 * 30], L"TrgASK_Real_VZ");

	 wcscpy( &wcharrFileNames[22 * 30], L"TrgASK_Real_R");
	 wcscpy( &wcharrFileNames[23 * 30], L"TrgASK_Real_V");
	 wcscpy( &wcharrFileNames[24 * 30], L"TrgASK_Real_U");

	 wcscpy( &wcharrFileNames[25 * 30], L"TrgASK_Est_X");
	 wcscpy( &wcharrFileNames[26 * 30], L"TrgASK_Est_Y");
	 wcscpy( &wcharrFileNames[27 * 30], L"TrgASK_Est_Z");
	 wcscpy( &wcharrFileNames[28 * 30], L"TrgASK_Est_VX");
	 wcscpy( &wcharrFileNames[29 * 30], L"TrgASK_Est_VY");
	 wcscpy( &wcharrFileNames[30 * 30], L"TrgASK_Est_VZ");

	 wcscpy( &wcharrFileNames[31 * 30], L"SigASK_X");
	 wcscpy( &wcharrFileNames[32 * 30], L"SigASK_Y");
	 wcscpy( &wcharrFileNames[33 * 30], L"SigASK_Z");
	 wcscpy( &wcharrFileNames[34 * 30], L"SigASK_VX");
	 wcscpy( &wcharrFileNames[35 * 30], L"SigASK_VY");
	 wcscpy( &wcharrFileNames[36 * 30], L"SigASK_VZ");

	 wcscpy( &wcharrFileNames[37 * 30], L"TrgKGSK_Est_X");
	 wcscpy( &wcharrFileNames[38 * 30], L"TrgKGSK_Est_Y");
	 wcscpy( &wcharrFileNames[39 * 30], L"TrgKGSK_Est_Z");
	 wcscpy( &wcharrFileNames[40 * 30], L"TrgKGSK_Est_VX");
	 wcscpy( &wcharrFileNames[41 * 30], L"TrgKGSK_Est_VY");
	 wcscpy( &wcharrFileNames[42 * 30], L"TrgKGSK_Est_VZ");

	 wcscpy( &wcharrFileNames[43 * 30], L"TrgPSK_Est_X");
	 wcscpy( &wcharrFileNames[44 * 30], L"TrgPSK_Est_Y");
	 wcscpy( &wcharrFileNames[45 * 30], L"TrgPSK_Est_Z");
	 wcscpy( &wcharrFileNames[46 * 30], L"TrgPSK_Est_VX");
	 wcscpy( &wcharrFileNames[47 * 30], L"TrgPSK_Est_VY");
	 wcscpy( &wcharrFileNames[48 * 30], L"TrgPSK_Est_VZ");

	 wcscpy( &wcharrFileNames[49 * 30], L"CelUkBet");
	 wcscpy( &wcharrFileNames[50 * 30], L"CelUkEps");
	 wcscpy( &wcharrFileNames[51 * 30], L"SigCelUkBet");
	 wcscpy( &wcharrFileNames[52 * 30], L"SigCelUkEps");

	 wcscpy( &wcharrFileNames[53 * 30], L"X_Zv_ASK");
	 wcscpy( &wcharrFileNames[54 * 30], L"Y_Zv_ASK");
	 wcscpy( &wcharrFileNames[55 * 30], L"Z_Zv_ASK");

	 wcscpy( &wcharrFileNames[56 * 30], L"X_Zv_PSK");
	 wcscpy( &wcharrFileNames[57 * 30], L"Y_Zv_PSK");
	 wcscpy( &wcharrFileNames[58 * 30], L"Z_Zv_PSK");

	 wcscpy( &wcharrFileNames[59 * 30], L"X_Zv_KGSK");
	 wcscpy( &wcharrFileNames[60 * 30], L"Y_Zv_KGSK");
	 wcscpy( &wcharrFileNames[61 * 30], L"Z_Zv_KGSK");

	 wcscpy( &wcharrFileNames[62 * 30], L"X_Zv_GSK");
	 wcscpy( &wcharrFileNames[63 * 30], L"Y_Zv_GSK");
	 wcscpy( &wcharrFileNames[64 * 30], L"Z_Zv_GSK");

	 wcscpy( &wcharrFileNames[65 * 30], L"CelUkR");
	 wcscpy( &wcharrFileNames[66 * 30], L"SigCelUkR");

	 double *pscaleY = new double  [QUANT_COLS_BUFF] ;
	 for (int i = 0; i < QUANT_COLS_BUFF; i++) pscaleY[i] = 1;


	 pscaleY[14] = 100.;
	 pscaleY[15] = 100;
	 pscaleY[23] = 1000;
	 pscaleY[24] = 1000.;
	 pscaleY[49] = 1000.;
	 pscaleY[50] = 1000.;
	 pscaleY[51] = 1000.;
	 pscaleY[52] = 1000.;

		// 1. создание папки FightReport
		wchar_t wcharFightReportFold [300] = {0} ;
	 wcscpy(wcharFightReportFold, mpwcharrFoldReport);
	wcscat(wcharFightReportFold, L"\\FightReport");
	_wmkdir(wcharFightReportFold);


	 //1. выод  координат цели в АСК
	 wchar_t wcharrPath0 [300] = {0} ;
	 wcscpy(wcharrPath0, wcharFightReportFold);
	wcscat(wcharrPath0, L"\\ASK_Inf");
	_wmkdir(wcharrPath0);
	wcscat(wcharrPath0, L"\\");

	 for (int i = 16; i < 37; i++)
	 {
	 TYrWriteShapeFile::WriteOneReport(wcharrPath0  // путь к папке
								  ,mparrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,QUANT_COLS_BUFF // - к-во переменных о корорых накоплена информация в буфере
								  ,mQuantPntReport  //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30// максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i // номер переменной по оси Y
								  ,100 //  масштаб по оси X
								  ,pscaleY[i]  // масштаб по оси Y
								   )  ;
	 }

	  double *parrTemp = new double [ 7 * mQuantPntReport] ;
	 for (int i = 0; i < mQuantPntReport; i++)
	 {
	   int num0 = i * QUANT_COLS_BUFF ;
	   MtrxMinusMatrx(&mparrBuff [num0 +25], &mparrBuff [num0 +16], 6, 1, &parrTemp[ 7 * i]);
	   parrTemp[ 7 * i + 6 ] = mparrBuff [num0]; // время
	 }
	 wchar_t wcharrNames[30 * 7] = {0};
	  wcscpy( &wcharrNames[0 * 30], L"ErrAsk_X");
	  wcscpy( &wcharrNames[1 * 30], L"ErrAsk_Y");
	  wcscpy( &wcharrNames[2 * 30], L"ErrAsk_Z");
	  wcscpy( &wcharrNames[3 * 30], L"ErrAsk_VX");
	  wcscpy( &wcharrNames[4 * 30], L"ErrAsk_VY");
	  wcscpy( &wcharrNames[5 * 30], L"ErrAsk_VZ");
	  wcscpy( &wcharrNames[6 * 30], L"t");
	  for (int i =0; i < 6; i++)
	  {
		TYrWriteShapeFile::WriteOneReport(wcharrPath0  // путь к папке
								  ,parrTemp // массив с информацией - матрица nBuffRows x nBuffCols
								  ,7  // - к-во переменных о корорых накоплена информация в буфере
								  ,mQuantPntReport  //  - к-во точек
								  ,wcharrNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30// максимальная длина имени переменной
								  ,6  // номер переменной по оси X от 0 до 6
								  ,i // номер переменной по оси Y
								  ,100 //  масштаб по оси X
								  ,1  // масштаб по оси Y
								   )  ;

	  }
	wchar_t wchFileName [300] = {0} ;
	wcscpy(wchFileName, wcharrPath0 );
	wcscat(wchFileName, L"Axes.shp");
	TYrWriteShapeFile::CreateShpAxes(wchFileName,-40000.,40000,-40000.,40000) ;
	////
	//2.  ВЫВОД ИНФОРМ В КГСК
	 wchar_t wcharrPath1 [300] = {0} ;
	 wcscpy(wcharrPath1, wcharFightReportFold);

	wcscat(wcharrPath1, L"\\KGSK_Inf");
	_wmkdir(wcharrPath1);
	wcscat(wcharrPath1, L"\\");
	  for (int i =0; i < 6; i++)
	  {
		TYrWriteShapeFile::WriteOneReport(wcharrPath1  // путь к папке
								  ,mparrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,QUANT_COLS_BUFF // - к-во переменных о корорых накоплена информация в буфере
								  ,mQuantPntReport  //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30// максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i + 1 // номер переменной по оси Y
								  ,100 //  масштаб по оси X
								  ,1  // масштаб по оси Y
								   )  ;
	   TYrWriteShapeFile::WriteOneReport(wcharrPath1  // путь к папке
								  ,mparrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,QUANT_COLS_BUFF // - к-во переменных о корорых накоплена информация в буфере
								  ,mQuantPntReport  //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30// максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i + 37// номер переменной по оси Y
								  ,100 //  масштаб по оси X
								  ,1  // масштаб по оси Y
								   )  ;
	  }

	  for (int i = 0; i < mQuantPntReport; i++)
	 {
	   int num0 = i * QUANT_COLS_BUFF ;
	   MtrxMinusMatrx(&mparrBuff [num0 +1], &mparrBuff [num0 +37], 6, 1, &parrTemp[ 7 * i]);
	   parrTemp[ 7 * i + 6 ] = mparrBuff [num0]; // время
	 }
	  memset( wcharrNames, 0, sizeof(wchar_t) *30 * 7);
	  wcscpy( &wcharrNames[0 * 30], L"ErrKGSK_X");
	  wcscpy( &wcharrNames[1 * 30], L"ErrKGSK_Y");
	  wcscpy( &wcharrNames[2 * 30], L"ErrKGSK_Z");
	  wcscpy( &wcharrNames[3 * 30], L"ErrKGSK_VX");
	  wcscpy( &wcharrNames[4 * 30], L"ErrKGSK_VY");
	  wcscpy( &wcharrNames[5 * 30], L"ErrKGSK_VZ");
	  wcscpy( &wcharrNames[6 * 30], L"t");
	  for (int i =0; i < 6; i++)
	  {
		TYrWriteShapeFile::WriteOneReport(wcharrPath1 // путь к папке
								  ,parrTemp // массив с информацией - матрица nBuffRows x nBuffCols
								  ,7  // - к-во переменных о корорых накоплена информация в буфере
								  ,mQuantPntReport  //  - к-во точек
								  ,wcharrNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30// максимальная длина имени переменной
								  ,6  // номер переменной по оси X от 0 до 6
								  ,i // номер переменной по оси Y
								  ,100 //  масштаб по оси X
								  ,1  // масштаб по оси Y
								   )  ;

	  }

	wchar_t wchFileName1 [300] = {0} ;
	wcscpy(wchFileName1, wcharrPath1 );
	wcscat(wchFileName1, L"Axes.shp");
	TYrWriteShapeFile::CreateShpAxes(wchFileName1,-40000.,40000,-40000.,40000) ;
	  ///

		//3.  ВЫВОД ИНФОРМ В PSK
 wchar_t wcharrPath2 [300] = {0} ;
	 wcscpy(wcharrPath2, wcharFightReportFold);

	wcscat(wcharrPath2, L"\\PSK_Inf");
	_wmkdir(wcharrPath2);
	wcscat(wcharrPath2, L"\\");
	for (int i =0; i < 6; i++)
	  {
		TYrWriteShapeFile::WriteOneReport(wcharrPath2  // путь к папке
								  ,mparrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,QUANT_COLS_BUFF // - к-во переменных о корорых накоплена информация в буфере
								  ,mQuantPntReport  //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30// максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i + 7 // номер переменной по оси Y
								  ,100 //  масштаб по оси X
								  ,1  // масштаб по оси Y
								   )  ;
	   TYrWriteShapeFile::WriteOneReport(wcharrPath2  // путь к папке
								  ,mparrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,QUANT_COLS_BUFF // - к-во переменных о корорых накоплена информация в буфере
								  ,mQuantPntReport  //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30// максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i + 43// номер переменной по оси Y
								  ,100 //  масштаб по оси X
								  ,1  // масштаб по оси Y
								   )  ;
	  }

	  for (int i = 0; i < mQuantPntReport; i++)
	 {
	   int num0 = i * QUANT_COLS_BUFF ;
	   MtrxMinusMatrx(&mparrBuff [num0 +7], &mparrBuff [num0 +43], 6, 1, &parrTemp[ 7 * i]);
	   parrTemp[ 7 * i + 6 ] = mparrBuff [num0]; // время
	 }
	  memset( wcharrNames, 0, sizeof(wchar_t) *30 * 7);
	  wcscpy( &wcharrNames[0 * 30], L"ErrPSK_X");
	  wcscpy( &wcharrNames[1 * 30], L"ErrPSK_Y");
	  wcscpy( &wcharrNames[2 * 30], L"ErrPSK_Z");
	  wcscpy( &wcharrNames[3 * 30], L"ErrPSK_VX");
	  wcscpy( &wcharrNames[4 * 30], L"ErrPSK_VY");
	  wcscpy( &wcharrNames[5 * 30], L"ErrPSK_VZ");
	  wcscpy( &wcharrNames[6 * 30], L"t");
	  for (int i =0; i < 6; i++)
	  {
		TYrWriteShapeFile::WriteOneReport(wcharrPath2 // путь к папке
								  ,parrTemp // массив с информацией - матрица nBuffRows x nBuffCols
								  ,7  // - к-во переменных о корорых накоплена информация в буфере
								  ,mQuantPntReport  //  - к-во точек
								  ,wcharrNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30// максимальная длина имени переменной
								  ,6  // номер переменной по оси X от 0 до 6
								  ,i // номер переменной по оси Y
								  ,100 //  масштаб по оси X
								  ,1  // масштаб по оси Y
								   )  ;

	  }
	  wchar_t wchFileName2 [300] = {0} ;
	wcscpy(wchFileName2, wcharrPath2 );
	wcscat(wchFileName2, L"Axes.shp");
	TYrWriteShapeFile::CreateShpAxes(wchFileName2,-40000.,40000,-40000.,40000) ;
	  ///

			//3.  ВЫВОД ИНФОРМ В целеуазаниям
	wchar_t wcharrPath3 [300] = {0} ;
	 wcscpy(wcharrPath3, wcharFightReportFold);
	// wcscat(wcharrPath3, L"FightReport\\CELEUKAZANIA\\");
	wcscat(wcharrPath3, L"\\CELEUKAZANIA");
		_wmkdir(wcharrPath3);
  wcscat(wcharrPath3, L"\\");
	 pscaleY[0] = 1000.;
	 pscaleY[1] = 1000000;
	 pscaleY[2] = 1000000;
	for (int i =0; i < 3; i++) // истинные V, R, U
	  {
		TYrWriteShapeFile::WriteOneReport(wcharrPath3  // путь к папке
								  ,mparrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,QUANT_COLS_BUFF // - к-во переменных о корорых накоплена информация в буфере
								  ,mQuantPntReport  //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30// максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i + 13 // номер переменной по оси Y
								  ,100 //  масштаб по оси X
								  ,1.//pscaleY[i]  // масштаб по оси Y
								   )  ;

	  }


	 pscaleY[0] = 1.;
	 pscaleY[1] = 1.;
	 pscaleY[2] = 1000000;
	 pscaleY[3] = 1000000;
	  for (int i =0; i < 4; i++) // целеуказания по  V,  U и скз целеуказаний по V, U
	  {
		TYrWriteShapeFile::WriteOneReport(wcharrPath3  // путь к папке
								  ,mparrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,QUANT_COLS_BUFF // - к-во переменных о корорых накоплена информация в буфере
								  ,mQuantPntReport  //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30// максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i + 49 // номер переменной по оси Y
								  ,100 //  масштаб по оси X
								  ,pscaleY[i]  // масштаб по оси Y
								   )  ;

	  }


	 pscaleY[0] = 1.;
	 pscaleY[1] = 1000.;
	   for (int i =0; i < 2; i++) // целеуказания по  R и скз целеуказаний по R
	  {
		TYrWriteShapeFile::WriteOneReport(wcharrPath3  // путь к папке
								  ,mparrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,QUANT_COLS_BUFF // - к-во переменных о корорых накоплена информация в буфере
								  ,mQuantPntReport  //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30// максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i + 65 // номер переменной по оси Y
								  ,100 //  масштаб по оси X
								  ,pscaleY[i]  // масштаб по оси Y
								   )  ;

	  }

	  for (int i = 0; i < ( mQuantPntReport -1 ) ; i++) // ошибки целеуказаний
	 {
	   int num0 = (i + 1) * QUANT_COLS_BUFF ;
	   int num1 = i * QUANT_COLS_BUFF ;
	   MtrxMinusMatrx(&mparrBuff [num0 + 14], &mparrBuff [num1 +49], 2, 1, &parrTemp[ 7 * i]);

	   double arrBet[3] = {0.},arrModulBet[3] = {0.};
	   arrBet[0] = parrTemp[ 7 * i];
	   arrBet[1] = parrTemp[ 7 * i] + 2.* M_PI;
	   arrBet[2] = parrTemp[ 7 * i] - 2.* M_PI;
	   for(int ii = 0; ii < 3;ii++)
	   {
		arrModulBet[ii] = fabs( arrBet[ii] ) ;
       }
	   int NumArgMin = -1;
	   parrTemp[ 7 * i] = arrBet[ NumArgMin ];

	   parrTemp[ 7 * i + 3 ] =  mparrBuff [num0 + 16] - mparrBuff [num0 + 65] ;


	 }
	  memset( wcharrNames, 0, sizeof(wchar_t) *30 * 7);
	  wcscpy( &wcharrNames[0 * 30], L"DelCelUkBet");
	  wcscpy( &wcharrNames[1 * 30], L"DelCelUkEps");
	  wcscpy( &wcharrNames[2 * 30], L"DelCelUkR");
	  wcscpy( &wcharrNames[6 * 30], L"t");
	 pscaleY[0] = 1000000;
	 pscaleY[1] = 1000000;
	 pscaleY[2] = 1000.;

	  for (int i =0; i < 3; i++)
	  {
		TYrWriteShapeFile::WriteOneReport(wcharrPath3 // путь к папке
								  ,parrTemp // массив с информацией - матрица nBuffRows x nBuffCols
								  ,7  // - к-во переменных о корорых накоплена информация в буфере
								  ,mQuantPntReport-1  //  - к-во точек
								  ,wcharrNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30// максимальная длина имени переменной
								  ,6  // номер переменной по оси X от 0 до 4
								  ,i  // номер переменной по оси Y
								  ,100 //  масштаб по оси X
								   ,pscaleY[i]  // масштаб по оси Y
								   )  ;

		}

			wchar_t wchFileName3  [300] = {0} ;
	 wcscpy(wchFileName3 , wcharFightReportFold);
	// wcscat(wcharrPath3, L"FightReport\\CELEUKAZANIA\\");
	wcscat(wchFileName3 , L"\\Axes.shp");
	TYrWriteShapeFile::CreateShpAxes(wchFileName3,-40000.,40000,-40000.,40000) ;
	  ///

	delete pscaleY ;
	delete parrTemp ;

	// ВЫВОД ИНФОРМАЦМИИ ПО ЗАМЕРАМ В РАЗНЫХ СК
	 wchar_t wcharrPath4 [300] = {0} ;
	 wcscpy(wcharrPath4, wcharFightReportFold);
	 wcscat(wcharrPath4, L"\\Measures");
	 _wmkdir(wcharrPath4);
	 wcscat(wcharrPath4, L"\\");
	 for (int i = 53; i < 65; i++)
	 {
	 TYrWriteShapeFile::WriteOneReport(wcharrPath4  // путь к папке
								  ,mparrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,QUANT_COLS_BUFF // - к-во переменных о корорых накоплена информация в буфере
								  ,mQuantPntReport  //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30// максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i // номер переменной по оси Y
								  ,100 //  масштаб по оси X
								  ,1 // масштаб по оси Y
								   )  ;
	 }
	 wchar_t wchFileName4 [300] = {0} ;
	wcscpy(wchFileName4, wcharrPath4 );
	wcscat(wchFileName4, L"Axes.shp");
	TYrWriteShapeFile::CreateShpAxes(wchFileName4,-40000.,40000,-40000.,40000) ;

 }


 // публикация  отчета  перегруженная
 void TFight::WriteReport(wchar_t *pwcharrPath)
 {
	 if (pwcharrPath == NULL )return ;
	 wchar_t wcharrPath [300] = {0} ;
	 wcscpy(wcharrPath, pwcharrPath);
	 wcscat(wcharrPath, L"\\");






	 wchar_t *wcharrFileNames = new wchar_t[QUANT_COLS_BUFF * 30];
	 memset(wcharrFileNames, 0,  QUANT_COLS_BUFF * 30 * sizeof( wchar_t)) ;

	 wcscpy( &wcharrFileNames[ 0 * 30], L"t");
	 wcscpy( &wcharrFileNames[ 1 * 30], L"TrgKGSK_Real_X");
	 wcscpy( &wcharrFileNames[ 2* 30],  L"TrgKGSK_Real_Y");
	 wcscpy( &wcharrFileNames[ 3 * 30], L"TrgKGSK_Real_Z");
	 wcscpy( &wcharrFileNames[ 4 * 30], L"TrgKGSK_Real_VX");
	 wcscpy( &wcharrFileNames[ 5 * 30], L"TrgKGSK_Real_VY");
	 wcscpy( &wcharrFileNames[ 6 * 30], L"TrgKGSK_Real_VZ");

	 wcscpy( &wcharrFileNames[ 7 * 30], L"TrgPSK_Real_X");
	 wcscpy( &wcharrFileNames[ 8 * 30], L"TrgPSK_Real_Y");
	 wcscpy( &wcharrFileNames[ 9 * 30], L"TrgPSK_Real_Z");
	 wcscpy( &wcharrFileNames[10 * 30], L"TrgPSK_Real_VX");
	 wcscpy( &wcharrFileNames[11 * 30], L"TrgPSK_Real_VY");
	 wcscpy( &wcharrFileNames[12 * 30], L"TrgPSK_Real_VZ");

	 wcscpy( &wcharrFileNames[13 * 30], L"Trg_Real_R");
	 wcscpy( &wcharrFileNames[14 * 30], L"Trg_Real_Bet");
	 wcscpy( &wcharrFileNames[15 * 30], L"Trg_Real_Eps");

	 wcscpy( &wcharrFileNames[16 * 30], L"TrgASK_Real_X");
	 wcscpy( &wcharrFileNames[17 * 30], L"TrgASK_Real_Y");
	 wcscpy( &wcharrFileNames[18 * 30], L"TrgASK_Real_Z");
	 wcscpy( &wcharrFileNames[19 * 30], L"TrgASK_Real_VX");
	 wcscpy( &wcharrFileNames[20 * 30], L"TrgASK_Real_VY");
	 wcscpy( &wcharrFileNames[21 * 30], L"TrgASK_Real_VZ");

	 wcscpy( &wcharrFileNames[22 * 30], L"TrgASK_Real_R");
	 wcscpy( &wcharrFileNames[23 * 30], L"TrgASK_Real_V");
	 wcscpy( &wcharrFileNames[24 * 30], L"TrgASK_Real_U");

	 wcscpy( &wcharrFileNames[25 * 30], L"TrgASK_Est_X");
	 wcscpy( &wcharrFileNames[26 * 30], L"TrgASK_Est_Y");
	 wcscpy( &wcharrFileNames[27 * 30], L"TrgASK_Est_Z");
	 wcscpy( &wcharrFileNames[28 * 30], L"TrgASK_Est_VX");
	 wcscpy( &wcharrFileNames[29 * 30], L"TrgASK_Est_VY");
	 wcscpy( &wcharrFileNames[30 * 30], L"TrgASK_Est_VZ");

	 wcscpy( &wcharrFileNames[31 * 30], L"SigASK_X");
	 wcscpy( &wcharrFileNames[32 * 30], L"SigASK_Y");
	 wcscpy( &wcharrFileNames[33 * 30], L"SigASK_Z");
	 wcscpy( &wcharrFileNames[34 * 30], L"SigASK_VX");
	 wcscpy( &wcharrFileNames[35 * 30], L"SigASK_VY");
	 wcscpy( &wcharrFileNames[36 * 30], L"SigASK_VZ");

	 wcscpy( &wcharrFileNames[37 * 30], L"TrgKGSK_Est_X");
	 wcscpy( &wcharrFileNames[38 * 30], L"TrgKGSK_Est_Y");
	 wcscpy( &wcharrFileNames[39 * 30], L"TrgKGSK_Est_Z");
	 wcscpy( &wcharrFileNames[40 * 30], L"TrgKGSK_Est_VX");
	 wcscpy( &wcharrFileNames[41 * 30], L"TrgKGSK_Est_VY");
	 wcscpy( &wcharrFileNames[42 * 30], L"TrgKGSK_Est_VZ");

	 wcscpy( &wcharrFileNames[43 * 30], L"TrgPSK_Est_X");
	 wcscpy( &wcharrFileNames[44 * 30], L"TrgPSK_Est_Y");
	 wcscpy( &wcharrFileNames[45 * 30], L"TrgPSK_Est_Z");
	 wcscpy( &wcharrFileNames[46 * 30], L"TrgPSK_Est_VX");
	 wcscpy( &wcharrFileNames[47 * 30], L"TrgPSK_Est_VY");
	 wcscpy( &wcharrFileNames[48 * 30], L"TrgPSK_Est_VZ");

	 wcscpy( &wcharrFileNames[49 * 30], L"CelUkBet");
	 wcscpy( &wcharrFileNames[50 * 30], L"CelUkEps");
	 wcscpy( &wcharrFileNames[51 * 30], L"SigCelUkBet");
	 wcscpy( &wcharrFileNames[52 * 30], L"SigCelUkEps");

	 wcscpy( &wcharrFileNames[53 * 30], L"X_Zv_ASK");
	 wcscpy( &wcharrFileNames[54 * 30], L"Y_Zv_ASK");
	 wcscpy( &wcharrFileNames[55 * 30], L"Z_Zv_ASK");

	 wcscpy( &wcharrFileNames[56 * 30], L"X_Zv_PSK");
	 wcscpy( &wcharrFileNames[57 * 30], L"Y_Zv_PSK");
	 wcscpy( &wcharrFileNames[58 * 30], L"Z_Zv_PSK");

	 wcscpy( &wcharrFileNames[59 * 30], L"X_Zv_KGSK");
	 wcscpy( &wcharrFileNames[60 * 30], L"Y_Zv_KGSK");
	 wcscpy( &wcharrFileNames[61 * 30], L"Z_Zv_KGSK");

	 wcscpy( &wcharrFileNames[62 * 30], L"X_Zv_GSK");
	 wcscpy( &wcharrFileNames[63 * 30], L"Y_Zv_GSK");
	 wcscpy( &wcharrFileNames[64 * 30], L"Z_Zv_GSK");

	 wcscpy( &wcharrFileNames[65 * 30], L"CelUkR");
	 wcscpy( &wcharrFileNames[66 * 30], L"SigCelUkR");

	 double *pscaleY = new double  [QUANT_COLS_BUFF] ;
	 for (int i = 0; i < QUANT_COLS_BUFF; i++) pscaleY[i] = 1;


	 pscaleY[14] = 100.;
	 pscaleY[15] = 100;
	 pscaleY[23] = 1000;
	 pscaleY[24] = 1000.;
	 pscaleY[49] = 1000.;
	 pscaleY[50] = 1000.;
	 pscaleY[51] = 1000.;
	 pscaleY[52] = 1000.;



	 //1. выод  координат цели в АСК
	 wchar_t wcharrPath0 [300] = {0} ;
	 wcscpy(wcharrPath0, pwcharrPath);
	 wcscat(wcharrPath0, L"ASK_Inf\\");
	 _wmkdir(wcharrPath0);

	 for (int i = 16; i < 37; i++)
	 {
	 TYrWriteShapeFile::WriteOneReport(wcharrPath0  // путь к папке
								  ,mparrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,QUANT_COLS_BUFF // - к-во переменных о корорых накоплена информация в буфере
								  ,mQuantPntReport  //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30// максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i // номер переменной по оси Y
								  ,100 //  масштаб по оси X
								  ,pscaleY[i]  // масштаб по оси Y
								   )  ;
	 }

	  double *parrTemp = new double [ 7 * mQuantPntReport] ;
	 for (int i = 0; i < mQuantPntReport; i++)
	 {
	   int num0 = i * QUANT_COLS_BUFF ;
	   MtrxMinusMatrx(&mparrBuff [num0 +25], &mparrBuff [num0 +16], 6, 1, &parrTemp[ 7 * i]);
	   parrTemp[ 7 * i + 6 ] = mparrBuff [num0]; // время
	 }
	 wchar_t wcharrNames[30 * 7] = {0};
	  wcscpy( &wcharrNames[0 * 30], L"ErrAsk_X");
	  wcscpy( &wcharrNames[1 * 30], L"ErrAsk_Y");
	  wcscpy( &wcharrNames[2 * 30], L"ErrAsk_Z");
	  wcscpy( &wcharrNames[3 * 30], L"ErrAsk_VX");
	  wcscpy( &wcharrNames[4 * 30], L"ErrAsk_VY");
	  wcscpy( &wcharrNames[5 * 30], L"ErrAsk_VZ");
	  wcscpy( &wcharrNames[6 * 30], L"t");
	  for (int i =0; i < 6; i++)
	  {
		TYrWriteShapeFile::WriteOneReport(wcharrPath0  // путь к папке
								  ,parrTemp // массив с информацией - матрица nBuffRows x nBuffCols
								  ,7  // - к-во переменных о корорых накоплена информация в буфере
								  ,mQuantPntReport  //  - к-во точек
								  ,wcharrNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30// максимальная длина имени переменной
								  ,6  // номер переменной по оси X от 0 до 6
								  ,i // номер переменной по оси Y
								  ,100 //  масштаб по оси X
								  ,1  // масштаб по оси Y
								   )  ;

	  }
	wchar_t wchFileName [300] = {0} ;
	wcscpy(wchFileName, wcharrPath0 );
	wcscat(wchFileName, L"Axes.shp");
	TYrWriteShapeFile::CreateShpAxes(wchFileName,-40000.,40000,-40000.,40000) ;
	////
	//2.  ВЫВОД ИНФОРМ В КГСК

	 wchar_t wcharrPath1 [300] = {0} ;
	 wcscpy(wcharrPath1, pwcharrPath);
	 wcscat(wcharrPath1, L"KGSK_Inf\\");
	 _wmkdir(wcharrPath1);
	  for (int i =0; i < 6; i++)
	  {
		TYrWriteShapeFile::WriteOneReport(wcharrPath1  // путь к папке
								  ,mparrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,QUANT_COLS_BUFF // - к-во переменных о корорых накоплена информация в буфере
								  ,mQuantPntReport  //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30// максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i + 1 // номер переменной по оси Y
								  ,100 //  масштаб по оси X
								  ,1  // масштаб по оси Y
								   )  ;
	   TYrWriteShapeFile::WriteOneReport(wcharrPath1  // путь к папке
								  ,mparrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,QUANT_COLS_BUFF // - к-во переменных о корорых накоплена информация в буфере
								  ,mQuantPntReport  //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30// максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i + 37// номер переменной по оси Y
								  ,100 //  масштаб по оси X
								  ,1  // масштаб по оси Y
								   )  ;
	  }

	  for (int i = 0; i < mQuantPntReport; i++)
	 {
	   int num0 = i * QUANT_COLS_BUFF ;
	   MtrxMinusMatrx(&mparrBuff [num0 +1], &mparrBuff [num0 +37], 6, 1, &parrTemp[ 7 * i]);
	   parrTemp[ 7 * i + 6 ] = mparrBuff [num0]; // время
	 }
	  memset( wcharrNames, 0, sizeof(wchar_t) *30 * 7);
	  wcscpy( &wcharrNames[0 * 30], L"ErrKGSK_X");
	  wcscpy( &wcharrNames[1 * 30], L"ErrKGSK_Y");
	  wcscpy( &wcharrNames[2 * 30], L"ErrKGSK_Z");
	  wcscpy( &wcharrNames[3 * 30], L"ErrKGSK_VX");
	  wcscpy( &wcharrNames[4 * 30], L"ErrKGSK_VY");
	  wcscpy( &wcharrNames[5 * 30], L"ErrKGSK_VZ");
	  wcscpy( &wcharrNames[6 * 30], L"t");
	  for (int i =0; i < 6; i++)
	  {
		TYrWriteShapeFile::WriteOneReport(wcharrPath1 // путь к папке
								  ,parrTemp // массив с информацией - матрица nBuffRows x nBuffCols
								  ,7  // - к-во переменных о корорых накоплена информация в буфере
								  ,mQuantPntReport  //  - к-во точек
								  ,wcharrNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30// максимальная длина имени переменной
								  ,6  // номер переменной по оси X от 0 до 6
								  ,i // номер переменной по оси Y
								  ,100 //  масштаб по оси X
								  ,1  // масштаб по оси Y
								   )  ;

	  }

	wchar_t wchFileName1 [300] = {0} ;
	wcscpy(wchFileName1, wcharrPath1 );
	wcscat(wchFileName1, L"Axes.shp");
	TYrWriteShapeFile::CreateShpAxes(wchFileName1,-40000.,40000,-40000.,40000) ;
	  ///

		//3.  ВЫВОД ИНФОРМ В PSK

	 wchar_t wcharrPath2 [300] = {0} ;
	 wcscpy(wcharrPath2, pwcharrPath);
	 wcscat(wcharrPath2, L"KGSK_Inf\\");
	 _wmkdir(wcharrPath2);
	for (int i =0; i < 6; i++)
	  {
		TYrWriteShapeFile::WriteOneReport(wcharrPath2  // путь к папке
								  ,mparrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,QUANT_COLS_BUFF // - к-во переменных о корорых накоплена информация в буфере
								  ,mQuantPntReport  //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30// максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i + 7 // номер переменной по оси Y
								  ,100 //  масштаб по оси X
								  ,1  // масштаб по оси Y
								   )  ;
	   TYrWriteShapeFile::WriteOneReport(wcharrPath2  // путь к папке
								  ,mparrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,QUANT_COLS_BUFF // - к-во переменных о корорых накоплена информация в буфере
								  ,mQuantPntReport  //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30// максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i + 43// номер переменной по оси Y
								  ,100 //  масштаб по оси X
								  ,1  // масштаб по оси Y
								   )  ;
	  }

	  for (int i = 0; i < mQuantPntReport; i++)
	 {
	   int num0 = i * QUANT_COLS_BUFF ;
	   MtrxMinusMatrx(&mparrBuff [num0 +7], &mparrBuff [num0 +43], 6, 1, &parrTemp[ 7 * i]);
	   parrTemp[ 7 * i + 6 ] = mparrBuff [num0]; // время
	 }
	  memset( wcharrNames, 0, sizeof(wchar_t) *30 * 7);
	  wcscpy( &wcharrNames[0 * 30], L"ErrPSK_X");
	  wcscpy( &wcharrNames[1 * 30], L"ErrPSK_Y");
	  wcscpy( &wcharrNames[2 * 30], L"ErrPSK_Z");
	  wcscpy( &wcharrNames[3 * 30], L"ErrPSK_VX");
	  wcscpy( &wcharrNames[4 * 30], L"ErrPSK_VY");
	  wcscpy( &wcharrNames[5 * 30], L"ErrPSK_VZ");
	  wcscpy( &wcharrNames[6 * 30], L"t");
	  for (int i =0; i < 6; i++)
	  {
		TYrWriteShapeFile::WriteOneReport(wcharrPath2 // путь к папке
								  ,parrTemp // массив с информацией - матрица nBuffRows x nBuffCols
								  ,7  // - к-во переменных о корорых накоплена информация в буфере
								  ,mQuantPntReport  //  - к-во точек
								  ,wcharrNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30// максимальная длина имени переменной
								  ,6  // номер переменной по оси X от 0 до 6
								  ,i // номер переменной по оси Y
								  ,100 //  масштаб по оси X
								  ,1  // масштаб по оси Y
								   )  ;

	  }
	  wchar_t wchFileName2 [300] = {0} ;
	wcscpy(wchFileName2, wcharrPath2 );
	wcscat(wchFileName2, L"Axes.shp");
	TYrWriteShapeFile::CreateShpAxes(wchFileName2,-40000.,40000,-40000.,40000) ;
	  ///

			//3.  ВЫВОД ИНФОРМ В целеуазаниям

	  wchar_t wcharrPath3 [300] = {0} ;
	 wcscpy(wcharrPath3, pwcharrPath);
	 wcscat(wcharrPath3, L"CELEUKAZANIA\\");
	 _wmkdir(wcharrPath3);
	 pscaleY[0] = 1000.;
	 pscaleY[1] = 1000000;
	 pscaleY[2] = 1000000;
	for (int i =0; i < 3; i++) // истинные V, R, U
	  {
		TYrWriteShapeFile::WriteOneReport(wcharrPath3  // путь к папке
								  ,mparrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,QUANT_COLS_BUFF // - к-во переменных о корорых накоплена информация в буфере
								  ,mQuantPntReport  //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30// максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i + 13 // номер переменной по оси Y
								  ,100 //  масштаб по оси X
								  ,1.//pscaleY[i]  // масштаб по оси Y
								   )  ;

	  }


	 pscaleY[0] = 1.;
	 pscaleY[1] = 1.;
	 pscaleY[2] = 1000000;
	 pscaleY[3] = 1000000;
	  for (int i =0; i < 4; i++) // целеуказания по  V,  U и скз целеуказаний по V, U
	  {
		TYrWriteShapeFile::WriteOneReport(wcharrPath3  // путь к папке
								  ,mparrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,QUANT_COLS_BUFF // - к-во переменных о корорых накоплена информация в буфере
								  ,mQuantPntReport  //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30// максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i + 49 // номер переменной по оси Y
								  ,100 //  масштаб по оси X
								  ,pscaleY[i]  // масштаб по оси Y
								   )  ;

	  }


	 pscaleY[0] = 1.;
	 pscaleY[1] = 1000.;
	   for (int i =0; i < 2; i++) // целеуказания по  R и скз целеуказаний по R
	  {
		TYrWriteShapeFile::WriteOneReport(wcharrPath3  // путь к папке
								  ,mparrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,QUANT_COLS_BUFF // - к-во переменных о корорых накоплена информация в буфере
								  ,mQuantPntReport  //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30// максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i + 65 // номер переменной по оси Y
								  ,100 //  масштаб по оси X
								  ,pscaleY[i]  // масштаб по оси Y
								   )  ;

	  }

	  for (int i = 0; i < ( mQuantPntReport -1 ) ; i++) // ошибки целеуказаний
	 {
	   int num0 = (i + 1) * QUANT_COLS_BUFF ;
	   int num1 = i * QUANT_COLS_BUFF ;
	   MtrxMinusMatrx(&mparrBuff [num0 + 14], &mparrBuff [num1 +49], 2, 1, &parrTemp[ 7 * i]);

	   double arrBet[3] = {0.},arrModulBet[3] = {0.};
	   arrBet[0] = parrTemp[ 7 * i];
	   arrBet[1] = parrTemp[ 7 * i] + 2.* M_PI;
	   arrBet[2] = parrTemp[ 7 * i] - 2.* M_PI;
	   for(int ii = 0; ii < 3;ii++)
	   {
		arrModulBet[ii] = fabs( arrBet[ii] ) ;
       }
	   int NumArgMin = -1;

	   parrTemp[ 7 * i] = arrBet[ NumArgMin ];

	   parrTemp[ 7 * i + 3 ] =  mparrBuff [num0 + 16] - mparrBuff [num0 + 65] ;

	 }
	  memset( wcharrNames, 0, sizeof(wchar_t) *30 * 7);
	  wcscpy( &wcharrNames[0 * 30], L"DelCelUkBet");
	  wcscpy( &wcharrNames[1 * 30], L"DelCelUkEps");
	  wcscpy( &wcharrNames[2 * 30], L"DelCelUkR");
	  wcscpy( &wcharrNames[6 * 30], L"t");
	 pscaleY[0] = 1000000;
	 pscaleY[1] = 1000000;
	 pscaleY[2] = 1000.;

	  for (int i =0; i < 3; i++)
	  {
		TYrWriteShapeFile::WriteOneReport(wcharrPath3 // путь к папке
								  ,parrTemp // массив с информацией - матрица nBuffRows x nBuffCols
								  ,7  // - к-во переменных о корорых накоплена информация в буфере
								  ,mQuantPntReport-1  //  - к-во точек
								  ,wcharrNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30// максимальная длина имени переменной
								  ,6  // номер переменной по оси X от 0 до 4
								  ,i  // номер переменной по оси Y
								  ,100 //  масштаб по оси X
								   ,pscaleY[i]  // масштаб по оси Y
								   )  ;

	  }
	  wchar_t wchFileName3 [300] = {0} ;
	wcscpy(wchFileName3, wcharrPath3 );
	wcscat(wchFileName3, L"Axes.shp");
	TYrWriteShapeFile::CreateShpAxes(wchFileName3,-40000.,40000,-40000.,40000) ;
	  ///

	delete pscaleY ;
	delete parrTemp ;

	// ВЫВОД ИНФОРМАЦМИИ ПО ЗАМЕРАМ В РАЗНЫХ СК


	  wchar_t wcharrPath4 [300] = {0} ;
	 wcscpy(wcharrPath4, pwcharrPath);
	 wcscat(wcharrPath4, L"Measures\\");
	 _wmkdir(wcharrPath4);
	 for (int i = 53; i < 65; i++)
	 {
	 TYrWriteShapeFile::WriteOneReport(wcharrPath4  // путь к папке
								  ,mparrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,QUANT_COLS_BUFF // - к-во переменных о корорых накоплена информация в буфере
								  ,mQuantPntReport  //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30// максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i // номер переменной по оси Y
								  ,100 //  масштаб по оси X
								  ,1 // масштаб по оси Y
								   )  ;
	 }
	 wchar_t wchFileName4 [300] = {0} ;
	wcscpy(wchFileName4, wcharrPath4 );
	wcscat(wchFileName4, L"Axes.shp");
	TYrWriteShapeFile::CreateShpAxes(wchFileName4,-40000.,40000,-40000.,40000) ;

 }

bool TFight::recalcFight( const double valT,const double valTExtrap)
 {
		// пересчет положения корабля на момент valT
		 if (fabs(valT - mTExtrap) > 0.000000001)
		 {
		 return false ;
		 }

		///

		mVessel.recalcVess(valT,mTargDesEps,mTargDesBet);


		mTarget.recalcTrajPoint(valT);
		///

		// расчет истинных координат цели в истинной линейной АСК
		double arrVSTargTrue_ASK[3] = {0.} , valTargTrue_R = 0, valTargTrue_V = 0., valTargTrue_U = 0.;
		double arrVSAntpTrue_ASK[3] = {0.} , valAntpTrue_R = 0, valAntpTrue_V = 0., valAntpTrue_U = 0.;
		double arrVSTargTrue_KGSK [3] = {0.}, arrVSAntpTrue_KGSK [3] = {0.} ;
		MtrxMinusMatrx(mTarget.mTraject.marrVectSostGSK, mVessel.marrVectSost,3, 1,arrVSTargTrue_KGSK) ;

		RecalcVect_KGSK_INTO_ASK_True (arrVSTargTrue_KGSK,arrVSTargTrue_ASK,3 ) ;
		TFight::Recalc_ASK_INTO_ASphSK(arrVSTargTrue_ASK, valTargTrue_R,valTargTrue_V, valTargTrue_U ) ;
		///

		// расчет истинных координат антипода в истинной линейной АСК
		memcpy(arrVSAntpTrue_KGSK, arrVSTargTrue_KGSK, 3 * sizeof(double));
		arrVSAntpTrue_KGSK[2] = -arrVSAntpTrue_KGSK[2];
		RecalcVect_KGSK_INTO_ASK_True (arrVSAntpTrue_KGSK,arrVSAntpTrue_ASK,3 ) ;
		TFight::Recalc_ASK_INTO_ASphSK(arrVSAntpTrue_ASK, valAntpTrue_R,valAntpTrue_V, valAntpTrue_U ) ;
		///

		TZamer CurZamer;
		if((mTarget.menumTargetType == GARPUN_V300) ||(mTarget.menumTargetType == GARPUN_V700)||(mTarget.menumTargetType == GARPUN_V1000)
		 ||(mTarget.menumTargetType ==JET_F16 )||(mTarget.menumTargetType == JET_A_10A ))
		 {
			mVessel.createMeasure(mAntCoeff,valTargTrue_R,valTargTrue_V, valTargTrue_U
			, valT,mTargDesEps,mTargDesBet,mEtalonSign, mTarget.mTargEPR
			, mVessel.mTransmitAnt.mPowerPrd, mVessel.mTransmitAnt.mKYPrd, &CurZamer ) ;
		 }

		 if(
		    (mTarget.menumTargetType == DESTROYER) ||(mTarget.menumTargetType ==  MOTORBOAT )
		 ||	(mTarget.menumTargetType == OPEN_MANPOWER_LIE)
		 ||(mTarget.menumTargetType == OPEN_MANPOWER_STAND)
		 ||(mTarget.menumTargetType == BULLET_PROOF_LIE)
		 ||(mTarget.menumTargetType == BULLET_PROOF_STAND)
		 ||(mTarget.menumTargetType == COVERED_MANPOWER_ENTRENCH)
		 ||(mTarget.menumTargetType == COVERED_MANPOWER_TRENCH)
		 ||(mTarget.menumTargetType ==  MANPOWER_ARMOURED_CARRIER)
		 ||(mTarget.menumTargetType ==  MANPOWER_CAR)
		 ||(mTarget.menumTargetType ==  PLATOON_POINT)
		 ||(mTarget.menumTargetType ==  GROUP_POINT_COAST)
		 )

		 {
			mVessel.createMeasure_ForSingleTarg(valTargTrue_R,valTargTrue_V, valTargTrue_U
			, valT,mTargDesEps,mTargDesBet,mEtalonSign, mTarget.mTargEPR
			, mVessel.mTransmitAnt.mPowerPrd, mVessel.mTransmitAnt.mKYPrd, &CurZamer ) ;
		 }
			 ///

		//  фильтр сопроводжджения




		TrackFiltration( CurZamer, valTExtrap ,TYPE_OF_FILTER
		,mTargDesBet,mTargDesEps ,mTargDesR
		, mSigTargDesBet,mSigTargDesEps,mSigTargDesR
		, marrVSExtrap_KGSK, marrSigExtrapPosit_KGSK, mSigModulV_KGSK );
		  ///



		// вычисление истинного положения цели в ПСК - АС

	// RecalcVect_GSK_INTO_PSK_AS_True (mTarget.mTraject.marrVectSostGSK, arrVSPSK_AS, 3 ) ;
	 mT =  valT ;
	 mTExtrap = valTExtrap;
	 updateReportData(CurZamer) ;
	 return true ;
 }

// пересчет координат из прямоугольной СК в сферическую CК
//  arrInp - координаты исходной точки
// valR - даольность
// valBet - угол отсчитанный по часовой стрелке от оси OY до напроравления на проекцию точки на плоскость OXY
// valEps - угол между радиус-вектором точки и горизонтальной плоскостью
void TFight::recalcCoord_INTO_Spherical(double *arrInp, double &valR, double &valBet, double &valEps)
{
	valR = sqrt ( arrInp[0] * arrInp[0] +arrInp[1] * arrInp[1] +arrInp[2] * arrInp[2] ) ;

	if (fabs(arrInp[2] / valR) > 0.99999)
	{
		 ShowMessage(L"ERROR");

		 return;
	}


	double temp =  arrInp[2] / valR;
	if (fabs(temp) > 0.99999999)
	{
		temp = SIGN_D(temp) *  0.99999999;
	}
	valEps  = asin(temp ) ;

	temp =  arrInp[0]/ sqrt (  arrInp[0] * arrInp[0] +arrInp[1] * arrInp[1]);
	if (fabs(temp) > 0.99999999)
	{
		temp = SIGN_D(temp) *  0.99999999;
	}

	valBet  = asin(temp) ;
	if (arrInp[1] < 0.)
	 {
		if (arrInp[0] < 0.)
		{
		  valBet = -M_PI - valBet ;
		}
		else
		{
		   valBet = M_PI - valBet ;
        }
	 }
}


// нахождение момента первого выстрела по цели
// VAlDistFireBegin - дальность точки встречи дальнего рубежа
// VAlDelT - шаг решения
// начиная с момента t=0 с шагом  VAlDelT
//  решается задача о точке встречи для идельного движения цели и корабля
//  как только дальность становится меньше  VAlDistFireBegin
//  время начала движения изделия возвращается
bool TFight::calcFirstShotTime(const double VAlDistFireBegin,  double *pvalFireBeginTime)
{

	TFight FightCur = *this;
	FightCur.mVessel.mMaxQ = 0.; // корабль двигается по идеальной прямой
	FightCur.mTarget.mTraject.marrSigW[0] = 0.;// цель двигается по идеальной прямой



	int iC = 1000./ mVessel.mControlSyst.mFiltT;
	double  valKGSKEps =0.,  valKGSKBet = 0.;
	for (int i = 0; i < iC; i++)
	{

		double valTCur = ((double)i) * mVessel.mControlSyst.mFiltT;
		//double valTCur = ((double)i) * 0.5;
		double arrVectAppointmentPoint[3] = {0.};
		FightCur.mVessel.recalcVess(valTCur, 0., 0.) ;
		FightCur.mTarget.recalcTrajPoint(valTCur) ;
		//
		double arrTargVS_KGSK0[6] ={0.};
		MtrxMinusMatrx(FightCur.mTarget.mTraject.marrVectSostGSK, FightCur.mVessel.marrVectSost  ,1, 6, arrTargVS_KGSK0);
		///

		//
		double arrPositionAY_KGSK[3] ={0.};
		FightCur.mVessel.calcAY_Position(arrPositionAY_KGSK);
		///



		double    valTFlight =0., valMiss = 0.;



		double *pvalEps0 = NULL, *pvalBet0 = NULL;
		double valKGSKEps1 = valKGSKEps -0.1;
		if (valKGSKEps1< 0.)
		{
		valKGSKEps1 = 0.002;
		}

		if ( i > 0)
		{
		pvalEps0 = &valKGSKEps1 ;
		pvalBet0 = &valKGSKBet;
		}

		if(!FightCur.calcAppointmentPoint(pvalEps0,pvalBet0, &FightCur.mVessel.marrVectSost[3], arrTargVS_KGSK0 , arrPositionAY_KGSK
		,&valKGSKEps,    &valKGSKBet, &valTFlight, arrVectAppointmentPoint, &valMiss))
		{
		return false;
		}

		double valDistCur = Norm3( arrVectAppointmentPoint);

		if (valDistCur <=  VAlDistFireBegin)
		{
		*pvalFireBeginTime  = valTCur;
		return true;
		}

	}
	return false;
}



// нахождение момента первого выстрела по цели
// VAlDistFireBegin - дальность точки встречи дальнего рубежа
// VAlDelT - шаг решения
// начиная с момента t=0 с шагом  VAlDelT
//  решается задача о точке встречи для идельного движения цели и корабля
//  как только дальность становится меньше  VAlDistFireBegin
//  время начала движения изделия возвращается
bool TFight::calcFirstShotTime_(const double VAlDistFireBegin,  double *pvalFireBeginTime,  double *pvalMostRemoteAppPointDist)
{
	double valTCur = 0.;
	double val_Dist = 0., val_DerivDist = 0.;
	double valD0 = calcDistAppointmentPoint( 0.)  ;

	if (valD0 < VAlDistFireBegin)
	{
	*pvalFireBeginTime = 1.;
	*pvalMostRemoteAppPointDist = valD0;
	return true;
	}

	for (int i =0 ; i < 100; i++)
	{
		calcDist_and_DerivDist_AppointmentPoint( valTCur, &val_Dist,&val_DerivDist);
		double valDelt = (val_Dist -VAlDistFireBegin) / val_DerivDist;
		valTCur = valTCur -  0.2 *valDelt;

		if (fabs(valDelt) < 0.1)
		{
			*pvalFireBeginTime = valTCur;
			if (valTCur < 1.)
			{
			valTCur = 1.;
			}
			*pvalMostRemoteAppPointDist = val_Dist;
			return true;
		}
	}
	return false ;
}



// нахождение момента  выстрела
// VAlDistFireBegin - дальность точки встречи дальнего рубежа
// VAlDelT - шаг решения
// начиная с момента t=0 с шагом  VAlDelT
//  решается задача о точке встречи для идельного движения цели и корабля
//  как только дальность становится меньше  VAlDistFireBegin
//  время начала движения изделия возвращается
bool TFight::calcShotTime(const double VAlAppointmentPointDist,  double *pvalFireTime)
{
	double valTCur = 0.;
	double val_Dist = 0., val_DerivDist = 0.;
	double valD0 = calcDistAppointmentPoint( 0.)  ;

	if (valD0 < VAlAppointmentPointDist)
	{

	return false;
	}

	for (int i =0 ; i < 100; i++)
	{

		calcDist_and_DerivDist_AppointmentPoint( valTCur, &val_Dist,&val_DerivDist);
		double valDelt = (val_Dist -VAlAppointmentPointDist) / val_DerivDist;
		valTCur = valTCur -  0.2 *valDelt;

		if (fabs(valDelt) < 0.1)
		{
			*pvalFireTime = valTCur;
			return true;
		}
	}
	return false ;
}

void TFight::calcDist_and_DerivDist_AppointmentPoint( const double VAlT, double *pval_Dist, double *pval_DerivDist)
{
  double valDelt = 0.1;
 *pval_Dist =  calcDistAppointmentPoint(  VAlT) ;

 int jjj = 0;
 double valDist1 = calcDistAppointmentPoint(  VAlT + valDelt) ;
 *pval_DerivDist = (valDist1 - (*pval_Dist))/ valDelt;
}
double TFight::calcDistAppointmentPoint( const double VAlT)
{
	TFight FightCur = *this;
	FightCur.mVessel.mMaxQ = 0.; // корабль двигается по идеальной прямой
	FightCur.mTarget.mTraject.marrSigW[0] = 0.;// цель двигается по идеальной прямой
	FightCur.mTarget.mTraject.marrSigW[1] = 0.;
	FightCur.mTarget.mTraject.marrSigW[2] = 0.;
	FightCur.shift(VAlT);
	///

	// экстраполяция вектиора состояния цели на время  (valTCurrentFight  - FightCur.mT) вперед
	double arrTargExtrapVS_GSK[9] = {0.};
	FightCur.mTarget.mTraject.extrapolateTargVS(VAlT - FightCur.mT, arrTargExtrapVS_GSK);
	///

	// экстраполяция вектиора состояния корабля на время  (valTCurrentFight - FightCur.mT) вперед
	double arrVessExtrapVS_GSK[9] ={0.};
	FightCur.mVessel.extrapolateTrueVS_GSK(VAlT - FightCur.mT, arrVessExtrapVS_GSK);
	///


	// вычисление истинного вект ора состояния цели в КГСК на момомент выстрела
	double arrTargVS_KGSK0[6] ={0.};
	MtrxMinusMatrx(arrTargExtrapVS_GSK, arrVessExtrapVS_GSK,1, 6, arrTargVS_KGSK0);
	///

	// вычисление вектора положения АУ в КГСК
	double arrPositionAY_KGSK[3] = {0.};
	FightCur.mVessel.calcAY_Position(arrPositionAY_KGSK);
	///
	double valKGSKEps = 0., valKGSKBet = 0., valTFlight = 0., valMiss = -1.;
	double arrVectAppointmentPointGSK [9] ={0.};
	FightCur.calcAppointmentPoint(NULL, NULL, &arrVessExtrapVS_GSK[3]
	,arrTargVS_KGSK0, arrPositionAY_KGSK
	, &valKGSKEps, &valKGSKBet, &valTFlight, arrVectAppointmentPointGSK, &valMiss) ;
	// экстраполяция вектиора состояния корабля на полетное время  valTFlight вперед
	double arrVessExtrapVS_GSK_0[9] ={0.};
	MatrxMultScalar(&arrVessExtrapVS_GSK[3], 3, 1, valTFlight,arrVessExtrapVS_GSK_0);

	///
	// положение точки встречи в КГСК на момент встречи
	double arrT0[3] ={0.};
	MtrxMinusMatrx(arrVessExtrapVS_GSK_0, arrVectAppointmentPointGSK,1, 3, arrT0);

	double valDist = Norm3(arrT0);
	return valDist;


}



// вычисление вероятности ОСНОВНАЯ ФУНКЦИЯ !!!!!!!!!!!!
// сначала ищется точка встречи на бли жний огневой рубеж
//  INPUT:
//  *pvalDistFireBegin  - дальний огн рубеж
//    уточняется в процессе работы функции, т.к захват цели может быть меньше этого рубежа
//  VAlDistFireFinish  - ближний огн рубеж
//  VAlHAntenna  - высота антеннвы
// OUTPUT:
// *pvalDistBeginSopr   - дальность начала сопровождения
// *pvalSigE, *pvalSigQ  - СКЗ ошибок угловых измерений на начало сопровождения
// *pvalFireBeginTime  - момент первого выстрела
// *piQuantShots  - к-во выстрелов
//  arrProbab[ILenArr]   - массив вероятностей попадания каждым из выстрелов
//  arrSKZPromach[ILenArr]  - массив СКЗ промахов
// ILenArr   - длина зарезервированной памыти под массивы
// *pvalProb  - вероятность поражения без Гладковского
// *pvalGladkProb  - вероятность поражения с Гладковским
//
//
bool TFight::calcSuccessProbAero(double *pvalDistFireBegin, const double VAlDistFireFinish0
   ,const double VAlHAntenna, double *pvalDistBeginSopr, double *pvalSigE, double *pvalSigQ, double *pvalFireBeginTime
   , int *piQuantShots,double *arrProbab, double *arrSKZPromach,double * arrSKZNedolet, double *arrSKZ_GSK_Z
   , double *arrCorMtrxCartinSK, double *arrDist, TURPolygon *plgArrProjection
	, const double ILenArr, double *pvalProb, double *pvalGladkProb)
{
		// нахождение начала сопровождения
	 bool bLat = true;

	*pvalDistBeginSopr =  mVessel.mFar_2D.calc_TwoTargsZahvatDist(mAntCoeff,mTarget.mTraject.marrVectSostGSK_Begin[2]
	, mTarget.mTargEPR, mVessel.mTransmitAnt.mPowerPrd,  mVessel.mTransmitAnt.mKYPrd, mEtalonSign, VAlHAntenna
	,&bLat  ,pvalSigE ,pvalSigQ  ) ;

	if ((*pvalDistBeginSopr)  < 0.)
	{  // нет разрешения
		(*pvalDistBeginSopr)  = 0;
		*pvalProb = 0.;
		*pvalGladkProb = 0.0;

		return false;  //
	}


	double valDist0 = Norm3(mTarget.mTraject.marrVectSostGSK_Begin);

	  double arrVOtn[3] = {0.};

		arrVOtn[0]= mTarget.mTraject.marrVectSostGSK_Begin[3] - mVessel.mVVess * sin(mVessel.mQ0) ;
		arrVOtn[1]= mTarget.mTraject.marrVectSostGSK_Begin[4] - mVessel.mVVess * cos(mVessel.mQ0) ;
		arrVOtn[2]= mTarget.mTraject.marrVectSostGSK_Begin[5];
		double valModVOtn = Norm3(arrVOtn);

		// выяснение ближнего огневого рубежа
		double valScal = ScalProduct(arrVOtn , mTarget.mTraject.marrVectSostGSK_Begin, 3) ;
		double valDistMin =  sqrt(valDist0 * valDist0 - valScal * valScal / valModVOtn/valModVOtn) + 100.;
		const double VAlDistFireFinish = (valDistMin <  VAlDistFireFinish0)? VAlDistFireFinish0: valDistMin;

		///
	if ((*pvalDistBeginSopr) < valDist0 )
	{
	  double val_TBeginSopr = calcT_BeginSopr( *pvalDistBeginSopr );
		// положение цели в КГСК в момент входа в зону сопровождения
		mTarget.mTraject.marrVectSostGSK_Begin[0] = mTarget.mTraject.marrVectSostGSK_Begin[0]  +  arrVOtn[0]* val_TBeginSopr;
		mTarget.mTraject.marrVectSostGSK_Begin[1] = mTarget.mTraject.marrVectSostGSK_Begin[1]  +  arrVOtn[1]* val_TBeginSopr;
		mTarget.mTraject.marrVectSostGSK_Begin[2] = mTarget.mTraject.marrVectSostGSK_Begin[2]  +  arrVOtn[2]* val_TBeginSopr;
		memcpy(mTarget.mTraject.marrVectSostGSK, mTarget.mTraject.marrVectSostGSK_Begin, 9 * sizeof(double));

		valDist0 = Norm3(mTarget.mTraject.marrVectSostGSK);
	}

	if ((*pvalDistFireBegin) > (valDist0 -  valModVOtn))
	{
		*pvalDistFireBegin =  valDist0 -  valModVOtn;

		if( (*pvalDistFireBegin) < VAlDistFireFinish)
		{ // дальность первого выстрела меньше ближнего рубежа
		(*pvalDistBeginSopr)   = 0;
		*pvalProb = 0.;
		*pvalGladkProb = 0.0;


		 return false;
		}
	}





	///////////////////////////////////////////////

	TFight FightCur0 (mVessel, mTarget,mVessel.mControlSyst.mFiltT,mEtalonSign
	,  mEnvironment, mAntCoeff, NULL ) ;
	*this =  FightCur0;

		///
   double valFireFinishTime = -1.;
	 if(!calcShotTime(VAlDistFireFinish0 ,  &valFireFinishTime))
	 {
		 return false;
	 }

	double valMostRemoteAppPointDist = -1.;
	if(!calcFirstShotTime_((*pvalDistFireBegin) ,  pvalFireBeginTime, &valMostRemoteAppPointDist))
	 {
	   return false;
	 }

	FightCur0.mVessel.mMaxQ = 0.; // корабль двигается по идеальной прямой

 //	memset( FightCur.mTarget.mTraject.marrSigW, 0, 3 * sizeof(double));  // цель двигается по идеальной прямой

	*piQuantShots =0;

	double valKGSKEps = 0., valKGSKBet  = 0., valKGSKEpsPrev = 0., valKGSKBetPrev  = 0.;
	double arrVectAppointmentPointGSK[6] = {0.};
	double valMiss = -1., valTFlight = -1.;
	double valTCurrentShot = -1.;

	for (int i = 0; i < ILenArr; i++)
	{
	 // для отладки
	 if (3 == i) {
		 int iii = 0;
	 }
	 ///



	TFight FightCur = FightCur0;
	valTCurrentShot = valFireFinishTime - (double(i)) *  FightCur.mVessel.mArtComplex.mArtCannon.mRateOfFire;
	// передвижение члена класса FightCur на время  valTCurrentFight
	// поскольку передвижение совершается дискретными шагами кратными
	// темпу фильтрации, то текущее время FightCur вообще гворя может быть меньше  valTCurrentFight
	FightCur.shift(valTCurrentShot);
	///

	// экстраполяция вектиора состояния цели на время  (valTCurrentFight  - FightCur.mT) вперед
	double arrTargExtrapVS_GSK[9] = {0.};
	FightCur.mTarget.mTraject.extrapolateTargVS(valTCurrentShot - FightCur.mT, arrTargExtrapVS_GSK);
	///

	// экстраполяция вектиора состояния корабля на время  (valTCurrentFight - FightCur.mT) вперед
	double arrVessExtrapVS_GSK[9] ={0.};
	FightCur.mVessel.extrapolateTrueVS_GSK(valTCurrentShot - FightCur.mT, arrVessExtrapVS_GSK);
	///


	// вычисление истинного вект ора состояния цели в КГСК на момомент выстрела
	double arrTargVS_KGSK0[6] ={0.};
	MtrxMinusMatrx(arrTargExtrapVS_GSK, arrVessExtrapVS_GSK,1, 6, arrTargVS_KGSK0);
	///

	// вычисление вектора положения АУ в КГСК
	double arrPositionAY_KGSK[3] = {0.};
	FightCur.mVessel.calcAY_Position(arrPositionAY_KGSK);
	///

	// вычисление точки встречи
	if (i ==0)
	{
	FightCur.calcAppointmentPoint(NULL, NULL, &arrVessExtrapVS_GSK[3]
	,arrTargVS_KGSK0, arrPositionAY_KGSK
	, &valKGSKEps, &valKGSKBet, &valTFlight, arrVectAppointmentPointGSK, &valMiss) ;
	valKGSKEpsPrev =  valKGSKEps -0.05;
	valKGSKBetPrev = valKGSKBet;

	}
	else
	{
	FightCur.calcAppointmentPoint(&valKGSKEpsPrev, &valKGSKBetPrev, &arrVessExtrapVS_GSK[3]
	,arrTargVS_KGSK0, arrPositionAY_KGSK
	, &valKGSKEps, &valKGSKBet, &valTFlight, arrVectAppointmentPointGSK, &valMiss) ;
	valKGSKEpsPrev =  valKGSKEps -0.05;
	valKGSKBetPrev = valKGSKBet;

	}


	// экстраполяция вектиора состояния корабля на полетное время  valTFlight вперед
	double arrVessExtrapVS_GSK_0[9] ={0.};
	MatrxMultScalar(&arrVessExtrapVS_GSK[3], 3, 1, valTFlight,arrVessExtrapVS_GSK_0);

	///
	// положение точки встречи в КГСК на момент встречи
	double arrT0[3] ={0.};
	MtrxMinusMatrx(arrVectAppointmentPointGSK, arrVessExtrapVS_GSK_0,1, 3, arrT0);
	arrDist [i] = Norm3(arrT0);
	if(arrDist [i] >  valMostRemoteAppPointDist)
	{
	break;
	}

	if (i > 0)
	{
		 if (arrDist [i] < arrDist [i-1])
		 {
			// break;
		 }
	}

	double valProb =0.,  valDispMiss = 0., valDispNedolet;


	calcProbability_For_Fixed_AppointmentPoint_AirTargs(valTCurrentShot, valTFlight, valKGSKEps, valKGSKBet
		,&valProb,  &valDispMiss, &valDispNedolet, &arrCorMtrxCartinSK[4 *(*piQuantShots)]
		  , &arrSKZ_GSK_Z[(*piQuantShots)],  &plgArrProjection[(*piQuantShots)]);
	 arrProbab[i] = valProb;
	 arrSKZPromach[i] =  sqrt(valDispMiss);
	 arrSKZNedolet[i] =  sqrt(valDispNedolet);

	(*piQuantShots)++;

	}

	flipArray(arrSKZPromach, *piQuantShots);
	flipArray(arrSKZNedolet, *piQuantShots);
	flipArray(arrProbab, *piQuantShots);
	flipArray(arrDist, *piQuantShots);
	flipArray(arrSKZ_GSK_Z, *piQuantShots);
	flipTwoDimArray(arrSKZ_GSK_Z, *piQuantShots, 4);
	flipTwoDimArray(arrCorMtrxCartinSK, *piQuantShots, 4);
	TURPolygon:: flipArray(plgArrProjection, *piQuantShots) ;

	calcAeroTargRezultProbability(arrProbab, arrDist, *piQuantShots, pvalProb, pvalGladkProb);

	return true;
}


// ПЕРЕГРУЖЕННАЯ!!!!!
// вычисление вероятности ОСНОВНАЯ ФУНКЦИЯ !!!!!!!!!!!!
// сначала ищется точка встречи на бли жний огневой рубеж
//  INPUT:
//  *pvalDistFireBegin  - дальний огн рубеж
//    уточняется в процессе работы функции, т.к захват цели может быть меньше этого рубежа
//  VAlDistFireFinish  - ближний огн рубеж
//  VAlHAntenna  - высота антеннвы
// OUTPUT:
// *pvalDistBeginSopr   - дальность начала сопровождения
// *pvalSigE, *pvalSigQ  - СКЗ ошибок угловых измерений на начало сопровождения
// *pvalFireBeginTime  - момент первого выстрела
// *piQuantShots  - к-во выстрелов
//  arrProbab[ILenArr]   - массив вероятностей попадания каждым из выстрелов
//  arrSKZPromach[ILenArr]  - массив СКЗ промахов
// ILenArr   - длина зарезервированной памыти под массивы
// *pvalProb  - вероятность поражения без Гладковского
// *pvalGladkProb  - вероятность поражения с Гладковским
//
//
bool TFight::calcSuccessProbAero(double *pvalDistFireBegin, const double VAlDistFireFinish0
   ,const double VAlHAntenna, double *pvalDistBeginSopr, double *pvalSigE, double *pvalSigQ, double *pvalFireBeginTime
   , int *piQuantShots,double *arrProbab, double *arrSKZPromach,double * arrSKZNedolet, double *arrSKZ_GSK_Z
   , double *arrCorMtrxCartinSK, double *arrDist, TURPolygon *plgArrProjection
	, const double ILenArr, double *pvalProb, double *pvalGladkProb
	, TOutPutAeroShot *parrOutPutAeroShot)
{
		// нахождение начала сопровождения
	 bool bLat = true;

	*pvalDistBeginSopr =  mVessel.mFar_2D.calc_TwoTargsZahvatDist(mAntCoeff,mTarget.mTraject.marrVectSostGSK_Begin[2]
	, mTarget.mTargEPR, mVessel.mTransmitAnt.mPowerPrd,  mVessel.mTransmitAnt.mKYPrd, mEtalonSign, VAlHAntenna
	,&bLat  ,pvalSigE ,pvalSigQ  ) ;

	if ((*pvalDistBeginSopr)  < 0.)
	{  // нет разрешения
		(*pvalDistBeginSopr)  = 0;
		*pvalProb = 0.;
		*pvalGladkProb = 0.0;

		return false;  //
	}


	double valDist0 = Norm3(mTarget.mTraject.marrVectSostGSK_Begin);

	  double arrVOtn[3] = {0.};

		arrVOtn[0]= mTarget.mTraject.marrVectSostGSK_Begin[3] - mVessel.mVVess * sin(mVessel.mQ0) ;
		arrVOtn[1]= mTarget.mTraject.marrVectSostGSK_Begin[4] - mVessel.mVVess * cos(mVessel.mQ0) ;
		arrVOtn[2]= mTarget.mTraject.marrVectSostGSK_Begin[5];
		double valModVOtn = Norm3(arrVOtn);

		// выяснение ближнего огневого рубежа
		double valScal = ScalProduct(arrVOtn , mTarget.mTraject.marrVectSostGSK_Begin, 3) ;
		double valDistMin =  sqrt(valDist0 * valDist0 - valScal * valScal / valModVOtn/valModVOtn) + 100.;
		const double VAlDistFireFinish = (valDistMin <  VAlDistFireFinish0)? VAlDistFireFinish0: valDistMin;

		///
	if ((*pvalDistBeginSopr) < valDist0 )
	{
	  double val_TBeginSopr = calcT_BeginSopr( *pvalDistBeginSopr );
		// положение цели в КГСК в момент входа в зону сопровождения
		mTarget.mTraject.marrVectSostGSK_Begin[0] = mTarget.mTraject.marrVectSostGSK_Begin[0]  +  arrVOtn[0]* val_TBeginSopr;
		mTarget.mTraject.marrVectSostGSK_Begin[1] = mTarget.mTraject.marrVectSostGSK_Begin[1]  +  arrVOtn[1]* val_TBeginSopr;
		mTarget.mTraject.marrVectSostGSK_Begin[2] = mTarget.mTraject.marrVectSostGSK_Begin[2]  +  arrVOtn[2]* val_TBeginSopr;
		memcpy(mTarget.mTraject.marrVectSostGSK, mTarget.mTraject.marrVectSostGSK_Begin, 9 * sizeof(double));

		valDist0 = Norm3(mTarget.mTraject.marrVectSostGSK);
	}

	if ((*pvalDistFireBegin) > (valDist0 -  valModVOtn))
	{
		*pvalDistFireBegin =  valDist0 -  valModVOtn;

		if( (*pvalDistFireBegin) < VAlDistFireFinish)
		{ // дальность первого выстрела меньше ближнего рубежа
		(*pvalDistBeginSopr)   = 0;
		*pvalProb = 0.;
		*pvalGladkProb = 0.0;


		 return false;
		}
	}





	///////////////////////////////////////////////


  //	TFight FightCur0 (mVessel, mTarget,mVessel.mControlSyst.mFiltT,mEtalonSign,mEnvironment
   //	,mMissSimulated,mPereletSimulated, mAntCoeff, NULL ) ;

		TFight FightCur0 (mVessel, mTarget,mVessel.mControlSyst.mFiltT,mEtalonSign,mEnvironment
	, mAntCoeff, NULL ) ;
	*this =  FightCur0;

		///
   double valFireFinishTime = -1.;
	 if(!calcShotTime(VAlDistFireFinish0 ,  &valFireFinishTime))
	 {
		 return false;
	 }

	double valMostRemoteAppPointDist = -1.;
	if(!calcFirstShotTime_((*pvalDistFireBegin) ,  pvalFireBeginTime, &valMostRemoteAppPointDist))
	 {
	   return false;
	 }

	FightCur0.mVessel.mMaxQ = 0.; // корабль двигается по идеальной прямой

 //	memset( FightCur.mTarget.mTraject.marrSigW, 0, 3 * sizeof(double));  // цель двигается по идеальной прямой

	*piQuantShots =0;

	double valKGSKEps = 0., valKGSKBet  = 0., valKGSKEpsPrev = 0., valKGSKBetPrev  = 0.;
	double arrVectAppointmentPointGSK[6] = {0.};
	double valMiss = -1., valTFlight = -1.;
	double valTCurrentShot = -1.;

	for (int i = 0; i < ILenArr; i++)
	{
	if(54 == i)
	{
	int hhhhi=0;;
	}

	TFight FightCur = FightCur0;
	valTCurrentShot = valFireFinishTime - (double(i)) *  FightCur.mVessel.mArtComplex.mArtCannon.mRateOfFire;
	parrOutPutAeroShot[(*piQuantShots)].mTshot =  valTCurrentShot; //**
	// передвижение члена класса FightCur на время  valTCurrentFight
	// поскольку передвижение совершается дискретными шагами кратными
	// темпу фильтрации, то текущее время FightCur вообще гворя может быть меньше  valTCurrentFight
	FightCur.shift(valTCurrentShot);
	///

	// экстраполяция вектиора состояния цели на время  (valTCurrentFight  - FightCur.mT) вперед
	double arrTargExtrapVS_GSK[9] = {0.};
	FightCur.mTarget.mTraject.extrapolateTargVS(valTCurrentShot - FightCur.mT, arrTargExtrapVS_GSK);
	memcpy(parrOutPutAeroShot[(*piQuantShots)].marrVS_GSK_Targ
	, FightCur.mTarget.mTraject.marrVectSostGSK, 6 * sizeof(double));  //**
	///

	// экстраполяция вектиора состояния корабля на время  (valTCurrentFight - FightCur.mT) вперед
	double arrVessExtrapVS_GSK[9] ={0.};
	FightCur.mVessel.extrapolateTrueVS_GSK(valTCurrentShot - FightCur.mT, arrVessExtrapVS_GSK);
	memcpy(parrOutPutAeroShot[(*piQuantShots)].marrVS_GSK_Vessel
	, FightCur.mVessel.marrVectSost, 6 * sizeof(double)); //**
	///


	// вычисление истинного вект ора состояния цели в КГСК на момомент выстрела
	double arrTargVS_KGSK0[6] ={0.};
	MtrxMinusMatrx(arrTargExtrapVS_GSK, arrVessExtrapVS_GSK,1, 6, arrTargVS_KGSK0);
	///

	// вычисление вектора положения АУ в КГСК
	double arrPositionAY_KGSK[3] = {0.};
	FightCur.mVessel.calcAY_Position(arrPositionAY_KGSK);
	///

	// вычисление точки встречи
	if (i ==0)
	{
	FightCur.calcAppointmentPoint(NULL, NULL, &arrVessExtrapVS_GSK[3]
	,arrTargVS_KGSK0, arrPositionAY_KGSK
	, &valKGSKEps, &valKGSKBet, &valTFlight, arrVectAppointmentPointGSK, &valMiss) ;
	valKGSKEpsPrev =  valKGSKEps -0.001;
	valKGSKBetPrev = valKGSKBet;

	}
	else
	{
	FightCur.calcAppointmentPoint(&valKGSKEpsPrev, &valKGSKBetPrev, &arrVessExtrapVS_GSK[3]
	,arrTargVS_KGSK0, arrPositionAY_KGSK
	, &valKGSKEps, &valKGSKBet, &valTFlight, arrVectAppointmentPointGSK, &valMiss) ;
	valKGSKEpsPrev =  valKGSKEps -0.001;
	valKGSKBetPrev = valKGSKBet;

	}

	parrOutPutAeroShot[(*piQuantShots)].mTfly  =valTFlight; //**
	parrOutPutAeroShot[(*piQuantShots)].mEpsKGSK =  valKGSKEps; //**
	parrOutPutAeroShot[(*piQuantShots)].mBetKGSK = valKGSKBet; //**
	// экстраполяция вектиора состояния корабля на полетное время  valTFlight вперед
	double arrVessExtrapVS_GSK_0[9] ={0.};
	MatrxMultScalar(&arrVessExtrapVS_GSK[3], 3, 1, valTFlight,arrVessExtrapVS_GSK_0);

	///
	// положение точки встречи в КГСК на момент встречи
	double arrT0[3] ={0.};
	MtrxMinusMatrx(arrVectAppointmentPointGSK, arrVessExtrapVS_GSK_0,1, 3, arrT0);
	arrDist [i] = Norm3(arrT0);
	if(arrDist [i] >  valMostRemoteAppPointDist)
	{
	break;
	}

	if (i > 0)
	{
		 if (arrDist [i] < arrDist [i-1])
		 {
			// break;
		 }
	}

	double valProb =0.,  valDispMiss = 2., valDispNedolet=0.1;


	calcProbability_For_Fixed_AppointmentPoint_AirTargs(valTCurrentShot, valTFlight, valKGSKEps, valKGSKBet
		,&valProb,  &valDispMiss, &valDispNedolet, &arrCorMtrxCartinSK[4 *(*piQuantShots)]
		, &arrSKZ_GSK_Z[(*piQuantShots)],  &plgArrProjection[(*piQuantShots)]
		, parrOutPutAeroShot[(*piQuantShots)].marrKTarg
	   ,parrOutPutAeroShot[(*piQuantShots)].marrKShell, parrOutPutAeroShot[(*piQuantShots)].marrMiss);   // ** **
	 arrProbab[i] = valProb;
	 arrSKZPromach[i] =  sqrt(valDispMiss);
	 arrSKZNedolet[i] =  sqrt(valDispNedolet);

	(*piQuantShots)++;

	}

	flipArray(arrSKZPromach, *piQuantShots);
	flipArray(arrSKZNedolet, *piQuantShots);
	flipArray(arrProbab, *piQuantShots);
	flipArray(arrDist, *piQuantShots);
	flipArray(arrSKZ_GSK_Z, *piQuantShots);
	flipTwoDimArray(arrSKZ_GSK_Z, *piQuantShots, 4);
	flipTwoDimArray(arrCorMtrxCartinSK, *piQuantShots, 4);
	if (mVessel.mShellBody.mDetonator.mEnumDetonatorType == CONTACT)
	{
		TURPolygon:: flipArray(plgArrProjection, *piQuantShots) ;
	}
	TOutPutAeroShot:: flipArray( parrOutPutAeroShot, *piQuantShots) ;

	calcAeroTargRezultProbability(arrProbab, arrDist, *piQuantShots, pvalProb, pvalGladkProb);

	return true;
}



	// нахождение точки встречи
	// задано текущий вектор состояния корабля и положение АУ в ГСК
	// задан текущий вектор состояния цели в ГСК
	// вычисляются ПУГН и ПУВН, находится точка встречи и подлетное время
	//INPUT:
	// arrTargVS_KGSK[6] - положение   цели
	// arrPositionAY_KGSK - вектор положения АУ в КГСК
	// *pEps0 - начальное значение УМ
	// *pBet0 - начальное значение КУ
	// Если pEps0 = NULL и  *pBet0 = NULL, то начальные значения вычисляются внутри функции
	//  OUTPUT:
	// valKGSKEps, valKGSKBet - углфы наведения в КГСК (двигающейся с кораблем)
	// valTFlight - подлетное вренмя
	// arrVectAppointmentPoint  -координаты точки встречи в КГСК зафиксированной на момент выстрела
	// valMiss - промах
	// возвращает true если точка встречи есть, и false в противнгом случае
	// НЕ ЗАБЫТЬ про положение АУ!!!!
bool TFight::calcAppointmentPoint(double *pEps0, double *pBet0,double *arrVessVelocity_GSK
	 , double *arrTargVS_KGSK0,double *arrPositionAY_KGSK
 , double * valKGSKEps, double * valKGSKBet, double * valTFlight, double *arrVectAppointmentPointGSK, double *valMiss)
{

	double dR = 2.;//5.; //0.25;//
	if (mTarget.menumTargetType == OPEN_MANPOWER_LIE)
	{
     dR = 50.;
	}

	//	double valF0 = 0., F0prev = 0.; // текущее и предыдущее значение целевой ф-ии
	double valEps0 = 0., valBet0 = 0.;

	// если начальные приближения углов заданы, то берутся они, если неьт, то вычисляются
	if (pEps0)
	{
	valEps0 =  *pEps0;
	valBet0 =  *pBet0;
	}
	else
	{
	fncInit_Eps0_Bet0(arrTargVS_KGSK0, mVessel.mShellBody.mV0, valEps0, valBet0)  ;
	}
	///

	// arrVSTargGSK0 - вектор состояния цели в неподвижной КГСК (сдвинутая ГСК)начальгый
	double arrTargVS_GSK0[6] ={0.};
	memcpy(arrTargVS_GSK0,arrTargVS_KGSK0, 6 * sizeof(double));
	MtrxSumMatrx(arrVessVelocity_GSK, &arrTargVS_GSK0[3],1, 3, &arrTargVS_GSK0[3]) ;
	///

	double valF0 = 0;
	double valF0prev = 1000000000.0;
	double val_dtInt = 0.01;
	bool bShag = false;
	//double cosQu, sinQu;
	//double Fi0prev, Qu0prev;
	bool breturn = false;
	int i =0;
	TMyShellTraj ShellTraj;

	double arrTargGSKPosCur[3] ={0.};// положение цели в точке встречи текущейж
	for ( i =0; i < 100; i++)
	{
		// valEps0,  valBet0 - углы наведения в КГСК !!!
		ShellTraj = TMyShellTraj(arrVessVelocity_GSK, mVessel.mShellBody,  valEps0,  valBet0 );

		///

		// нахождения углов наведения в ГСК
		double valEpsGSK0 = 0.,  valGSKBet0 = 0.;
		transform_EpsBetKGSK_to_EpsBetGSK(arrVessVelocity_GSK, mVessel.mShellBody.mV0,  valEps0,  valBet0
		,  &valEpsGSK0,  &valGSKBet0);
		///

		// пересчет вектора состояния цели в ССК
		double arrTargVS_SSK0[6] = {0.};
		ShellTraj.transform_xyzGSK_To_xyzSSK(6, arrTargVS_GSK0, arrTargVS_SSK0) ;
		///

		valF0 =  ShellTraj.calcPointMissMinimum(mEnvironment, arrTargVS_SSK0, val_dtInt);




		if (((valF0 < 20)  || (fabs(valF0 - valF0prev) < 0.1)) && ( !bShag)&& ( i >0))
		{
		bShag= true;
		val_dtInt=   0.0001;
		valF0prev = 10000.;//F0;
		continue;
		}



		if ((valF0 < dR)&& bShag)
		{
		breturn = true;
		break;
		}
		// вычисление вектора положения цели в момент  ShellTraj.mTCur
		double arrTargSSKPosCur[3] ={0.}, arrTemp[3]= {0.};
		MatrxMultScalar(&arrTargVS_SSK0[3], 1, 3, ShellTraj.mTCur, arrTemp);
		MtrxSumMatrx(arrTargVS_SSK0, arrTemp,1, 3, arrTargSSKPosCur) ;
		ShellTraj.transform_xyzSSK_To_xyzGSK( 3, arrTargSSKPosCur, arrTargGSKPosCur) ;
		///

		// пересчет вектора положения снаряда в ГСК
		double arrShellGSKPosCur[3] ={0.};
		ShellTraj.transform_xyzSSK_To_xyzGSK( 3, ShellTraj.marrStrSK_VS, arrShellGSKPosCur) ;
		///

		double temp =  arrShellGSKPosCur[2] / Norm3( arrShellGSKPosCur);
		if (fabs(temp) > 0.999999999)
		{
			temp = SIGN_D(temp) * 0.999999999;
		}
		double valShellEps = asin(temp);


		temp =  arrTargGSKPosCur[2]/ Norm3( arrTargGSKPosCur);
		if (fabs(temp) > 0.999999999)
		{
			temp = SIGN_D(temp) * 0.999999999;
		}
		double valTargEps = asin(temp);

		// dFi= acos((x1*x2+z1*z2)/(sqrt(x1*x1+z1*z1)*sqrt(x2*x2+z2*z2)));
		double val_dEps =	valTargEps -valShellEps;
		valEpsGSK0 = valEpsGSK0 + val_dEps * 0.7;
		//valEpsGSK0 = valEpsGSK0 + val_dEps * 0.4;
		//	 valEps0 = valEps0 + val_dEps * 0.8;

		double val_dBet = atan2(arrTargGSKPosCur[0],arrTargGSKPosCur[1])-atan2(arrShellGSKPosCur[0],arrShellGSKPosCur[1]);
		valGSKBet0 = valGSKBet0 + val_dBet *0.9;
		//valBet0= valBet0 + val_dBet *0.9;
		if(valGSKBet0 >= 2. * M_PI)
		{
		valGSKBet0 = valGSKBet0 - 2. * M_PI;
		}
		if(valGSKBet0 <= -2. * M_PI)
		{
		valGSKBet0 = valGSKBet0 + 2. * M_PI;
		}

		//calc_EpsKGSK0_BetKGSK0(arrVessVelocity_GSK, mVessel.mShellBody,  valEpsGSK0,  valGSKBet0
	   //	,  &valEps0,  &valBet0);
		transform_EpsBetGSK_to_EpsBetKGSK(arrVessVelocity_GSK, mVessel.mShellBody.mV0,  valEpsGSK0,  valGSKBet0
		,  &valEps0,  &valBet0);

		valF0prev = valF0;
	}


	if (!breturn)
	{
	return false;
	}

	if(valBet0 < 0.0)
	{
	 valBet0 +=  2.0*M_PI ;
	 }

	if(valBet0 >=  2.0*M_PI)
	{
	 valBet0 -=  2.0*M_PI ;
	 }

	*valKGSKEps = valEps0;
	*valKGSKBet = valBet0;
	*valTFlight =ShellTraj.mTCur;
	memcpy(arrVectAppointmentPointGSK,arrTargGSKPosCur, 3 * sizeof(double));
	*valMiss = valF0;

	return true;
}

//-----------------------------------------------------------------------

	// нахождение точки встречи   ПЕРЕГРУЖЕННАЯ!!!!!
	// задано текущий вектор состояния корабля и положение АУ в ГСК
	// задан текущий вектор состояния цели в ГСК
	// вычисляются ПУГН и ПУВН, находится точка встречи и подлетное время
	//INPUT:
	// arrTargVS_KGSK[6] - положение   цели
	// arrPositionAY_KGSK - вектор положения АУ в КГСК
	// *pEps0 - начальное значение УМ
	// *pBet0 - начальное значение КУ
	// Если pEps0 = NULL и  *pBet0 = NULL, то начальные значения вычисляются внутри функции
	//  OUTPUT:
	// valKGSKEps, valKGSKBet - углфы наведения в КГСК (двигающейся с кораблем)
	// valTFlight - подлетное вренмя
	// arrVectAppointmentPoint  -координаты точки встречи в КГСК зафиксированной на момент выстрела
	// valMiss - промах
	// возвращает true если точка встречи есть, и false в противнгом случае
	// НЕ ЗАБЫТЬ про положение АУ!!!!
	//arrShellVeloAppointmPoint_GSK [3] - скорость снаряда в точке встречи
	// arrTargVeloAppointmPoint_GSK[3] - скорость цели в точке встречи
	bool TFight::calcAppointmentPoint(double *pEps0, double *pBet0,double *arrVessVelocity_GSK
	, double *arrTargVS_KGSK0,double *arrPositionAY_KGSK
	, double * valKGSKEps, double * valKGSKBet, double * valTFlight, double *arrVectAppointmentPointGSK, double *valMiss
	, double *arrShellVeloAppointmPoint_GSK,  double *arrTargVeloAppointmPoint_GSK)
{
	double dR = 5.;
	switch(mTarget.menumTargetType)
	{
	case OPEN_MANPOWER_LIE:

	case OPEN_MANPOWER_STAND:

	case BULLET_PROOF_LIE:

	case BULLET_PROOF_STAND:

	case COVERED_MANPOWER_ENTRENCH:

	case COVERED_MANPOWER_TRENCH:

	case MANPOWER_ARMOURED_CARRIER:

	case MANPOWER_CAR:
	dR = 50.;
	break;
	default:
	break;
	}

	//	double valF0 = 0., F0prev = 0.; // текущее и предыдущее значение целевой ф-ии
	double valEps0 = 0., valBet0 = 0.;

	// если начальные приближения углов заданы, то берутся они, если неьт, то вычисляются
	if (pEps0)
	{
	valEps0 =  *pEps0;
	valBet0 =  *pBet0;
	}
	else
	{
	fncInit_Eps0_Bet0(arrTargVS_KGSK0, mVessel.mShellBody.mV0, valEps0, valBet0)  ;
	valEps0 = valEps0/2.;
	}
	///

	// arrVSTargGSK0 - вектор состояния цели в неподвижной КГСК (сдвинутая ГСК)начальгый
	double arrTargVS_GSK0[6] ={0.};
	memcpy(arrTargVS_GSK0,arrTargVS_KGSK0, 6 * sizeof(double));
	MtrxSumMatrx(arrVessVelocity_GSK, &arrTargVS_GSK0[3],1, 3, &arrTargVS_GSK0[3]) ;
	///

	double valF0 = 0;
	double valF0prev = 1000000000.0;
	double val_dtInt = 0.0005;
	bool bShag = false;
	//double cosQu, sinQu;
	//double Fi0prev, Qu0prev;
	bool breturn = false;
	int i =0;
	TMyShellTraj ShellTraj;

	double arrTargGSKPosCur[3] ={0.};// положение цели в точке встречи текущейж
	for ( i =0; i < 60; i++)
	{

		// valEps0,  valBet0 - углы наведения в КГСК !!!
		ShellTraj = TMyShellTraj(arrVessVelocity_GSK, mVessel.mShellBody,  valEps0,  valBet0 );
		///

		// нахождения углов наведения в ГСК
		double valEpsGSK0 = 0.,  valGSKBet0 = 0.;
		transform_EpsBetKGSK_to_EpsBetGSK(arrVessVelocity_GSK, mVessel.mShellBody.mV0,  valEps0,  valBet0
		,  &valEpsGSK0,  &valGSKBet0);
		///

		// пересчет вектора состояния цели в ССК
		double arrTargVS_SSK0[6] = {0.};
		ShellTraj.transform_xyzGSK_To_xyzSSK(6, arrTargVS_GSK0, arrTargVS_SSK0) ;
		///

		valF0 =  ShellTraj.calcPointMissMinimum(mEnvironment, arrTargVS_SSK0, val_dtInt);
		if (((valF0 < 20)  || (fabs(valF0 - valF0prev) < 0.1)) && ( !bShag))
		{
		bShag= true;
		val_dtInt= 0.0001; //0.001;//
		valF0prev = 10000.;//F0;
		continue;
		}



		if ((valF0 < dR)&& bShag)
		{
		breturn = true;
		break;
		}
		// вычисление вектора положения цели в момент  ShellTraj.mTCur
		double arrTargSSKPosCur[3] ={0.}, arrTemp[3]= {0.};
		MatrxMultScalar(&arrTargVS_SSK0[3], 1, 3, ShellTraj.mTCur, arrTemp);
		MtrxSumMatrx(arrTargVS_SSK0, arrTemp,1, 3, arrTargSSKPosCur) ;
		ShellTraj.transform_xyzSSK_To_xyzGSK( 3, arrTargSSKPosCur, arrTargGSKPosCur) ;
		///

		// пересчет вектора положения снаряда в ГСК
		double arrShellGSKPosCur[3] ={0.};
		ShellTraj.transform_xyzSSK_To_xyzGSK( 3, ShellTraj.marrStrSK_VS, arrShellGSKPosCur) ;
		///

		double temp =  arrShellGSKPosCur[2] / Norm3( arrShellGSKPosCur);
		if (fabs(temp) > 0.999999999)
		{
			temp = SIGN_D(temp) * 0.999999999;
		}
		double valShellEps = asin(temp);


		temp =  arrTargGSKPosCur[2]/ Norm3( arrTargGSKPosCur);
		if (fabs(temp) > 0.999999999)
		{
			temp = SIGN_D(temp) * 0.999999999;
		}
		double valTargEps = asin(temp);


		// dFi= acos((x1*x2+z1*z2)/(sqrt(x1*x1+z1*z1)*sqrt(x2*x2+z2*z2)));
		double val_dEps =	valTargEps -valShellEps;
		valEpsGSK0 = valEpsGSK0 + val_dEps * 1.2;

		//	 valEps0 = valEps0 + val_dEps * 0.8;

		double val_dBet = atan2(arrTargGSKPosCur[0],arrTargGSKPosCur[1])-atan2(arrShellGSKPosCur[0],arrShellGSKPosCur[1]);
		valGSKBet0 = valGSKBet0 + val_dBet *1.2;
		//valBet0= valBet0 + val_dBet *0.9;
		if(valGSKBet0 >= 2. * M_PI)
		{
		valGSKBet0 = valGSKBet0 - 2. * M_PI;
		}
		if(valGSKBet0 <= -2. * M_PI)
		{
		valGSKBet0 = valGSKBet0 + 2. * M_PI;
		}

	   //	calc_EpsKGSK0_BetKGSK0(arrVessVelocity_GSK, mVessel.mShellBody,  valEpsGSK0,  valGSKBet0
	   //	,  &valEps0,  &valBet0);

		transform_EpsBetGSK_to_EpsBetKGSK(arrVessVelocity_GSK, mVessel.mShellBody.mV0,  valEpsGSK0,  valGSKBet0
		,  &valEps0,  &valBet0);



		valF0prev = valF0;
	}


	if (!breturn)
	{
	return false;
	}

	if(valBet0 < 0.0) valBet0 +=  2.0*M_PI ;
	if(valBet0 >=  2.0*M_PI) valBet0 -=  2.0*M_PI ;

	*valKGSKEps = valEps0;
	*valKGSKBet = valBet0;
	*valTFlight =ShellTraj.mTCur;
	memcpy(arrVectAppointmentPointGSK,arrTargGSKPosCur, 3 * sizeof(double));
	*valMiss = valF0;

	memcpy(arrTargVeloAppointmPoint_GSK, &arrTargVS_GSK0[3], 3 * sizeof(double));

	 double arrVS_PrStSK[6] = {0.};
	ShellTraj.fncCalcVS_v_PrStSK(  arrVS_PrStSK  ) ;
	ShellTraj.transform_xyzSSK_To_xyzGSK( 3, &arrVS_PrStSK[3], arrShellVeloAppointmPoint_GSK) ;
	return true;
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

	// нахождение точки встречи   методом Ньютона
	// задано текущий вектор состояния корабля и положение АУ в ГСК
	// задан текущий вектор состояния цели в ГСК
	// вычисляются ПУГН и ПУВН, находится точка встречи и подлетное время
	//INPUT:
	// arrTargVS_KGSK[6] - положение   цели
	// arrPositionAY_KGSK - вектор положения АУ в КГСК
	// *pEps0 - начальное значение УМ
	// *pBet0 - начальное значение КУ
	// Если pEps0 = NULL и  *pBet0 = NULL, то начальные значения вычисляются внутри функции
	//  OUTPUT:
	// valKGSKEps, valKGSKBet - углфы наведения в КГСК (двигающейся с кораблем)
	// valTFlight - подлетное вренмя
	// arrVectAppointmentPoint  -координаты точки встречи в КГСК зафиксированной на момент выстрела
	// valMiss - промах
	// возвращает true если точка встречи есть, и false в противнгом случае
	// НЕ ЗАБЫТЬ про положение АУ!!!!
	//arrShellVeloAppointmPoint_GSK [3] - скорость снаряда в точке встречи
	// arrTargVeloAppointmPoint_GSK[3] - скорость цели в точке встречи
	bool TFight::findAppointmentPoint_NewtonMeth(double *arrVessVelocity_GSK
	, double *arrTargVS_KGSK0,double *arrPositionAY_KGSK
	, double * valKGSKEps, double * valKGSKBet, double * valTFlight, double *arrVectAppointmentPointGSK, double *valMiss
	, double *arrShellVeloAppointmPoint_GSK,  double *arrTargVeloAppointmPoint_GSK)
{
	double dR = 5;
	switch(mTarget.menumTargetType)
	{
	case OPEN_MANPOWER_LIE:

	case OPEN_MANPOWER_STAND:

	case BULLET_PROOF_LIE:

	case BULLET_PROOF_STAND:

	case COVERED_MANPOWER_ENTRENCH:

	case COVERED_MANPOWER_TRENCH:

	case MANPOWER_ARMOURED_CARRIER:

	case MANPOWER_CAR:
	dR = 50.;
	break;
	default:
	break;
	}


	// если начальные приближения углов заданы, то берутся они, если неьт, то вычисляются
	double valEps0 = 0., valBet0 = 0.;
	fncInit_Eps0_Bet0(arrTargVS_KGSK0, mVessel.mShellBody.mV0, valEps0, valBet0)  ;

	///

	// arrVSTargGSK0 - вектор состояния цели в неподвижной КГСК (сдвинутая ГСК)начальгый
	double arrTargVS_GSK0[6] ={0.};
	memcpy(arrTargVS_GSK0,arrTargVS_KGSK0, 6 * sizeof(double));
	MtrxSumMatrx(arrVessVelocity_GSK, &arrTargVS_GSK0[3],1, 3, &arrTargVS_GSK0[3]) ;
	///

	double valF0 = 0;
 //	double valF0prev = 1000000000.0;
	double val_dtInt = 0.00005;
 //	bool bShag = false;
	//double cosQu, sinQu;
	//double Fi0prev, Qu0prev;
	bool breturn = false;
	int i =0;
	TMyShellTraj ShellTraj;

	double arrTargGSKPosCur[3] ={0.};// положение цели в точке встречи текущейж
	for ( i =0; i < 60; i++)
	{

		// valEps0,  valBet0 - углы наведения в КГСК !!!
		ShellTraj = TMyShellTraj(arrVessVelocity_GSK, mVessel.mShellBody,  valEps0,  valBet0 );
		///

		// нахождения углов наведения в ГСК
		double valEpsGSK0 = 0.,  valGSKBet0 = 0.;
		transform_EpsBetKGSK_to_EpsBetGSK(arrVessVelocity_GSK, mVessel.mShellBody.mV0,  valEps0,  valBet0
		,  &valEpsGSK0,  &valGSKBet0);
		///

		// пересчет вектора состояния цели в ССК
		double arrTargVS_SSK0[6] = {0.};
		ShellTraj.transform_xyzGSK_To_xyzSSK(6, arrTargVS_GSK0, arrTargVS_SSK0) ;
		///

		valF0 =  ShellTraj.calcPointMissMinimum(mEnvironment, arrTargVS_SSK0, val_dtInt);
	/*	if (((valF0 < 20)  || (fabs(valF0 - valF0prev) < 0.1)) && ( !bShag))
		{
		bShag= true;
		val_dtInt= 0.0001; //0.001;//
		valF0prev = 10000.;//F0;
		continue;
		}  */

    	// вычисление вектора положения цели в момент  ShellTraj.mTCur
		double arrTargSSKPosCur[3] ={0.}, arrTemp[3]= {0.};
		MatrxMultScalar(&arrTargVS_SSK0[3], 1, 3, ShellTraj.mTCur, arrTemp);
		MtrxSumMatrx(arrTargVS_SSK0, arrTemp,1, 3, arrTargSSKPosCur) ;
		ShellTraj.transform_xyzSSK_To_xyzGSK( 3, arrTargSSKPosCur, arrTargGSKPosCur) ;
		///

		// пересчет вектора положения снаряда в ГСК
		double arrShellGSKPosCur[3] ={0.};
		ShellTraj.transform_xyzSSK_To_xyzGSK( 3, ShellTraj.marrStrSK_VS, arrShellGSKPosCur) ;
		///
		double arrDElta[3] = {0.};
		MtrxMinusMatrx(arrShellGSKPosCur, arrTargGSKPosCur,3, 1, arrDElta);



		double valDelta = 0.001;
		double arrF_SSKT[9] = {0.};
		double arrVectPartialDeriv [8] = {0.};
		ShellTraj.fncCalcVectPartialDeriv(mEnvironment, val_dtInt, 3,   valDelta,  arrVectPartialDeriv );
		memcpy(arrF_SSKT, arrVectPartialDeriv, 3 * sizeof(double));
		ShellTraj.fncCalcVectPartialDeriv(mEnvironment, val_dtInt, 7,   valDelta,  arrVectPartialDeriv );
		memcpy(&arrF_SSKT[3], arrVectPartialDeriv, 3 * sizeof(double));


		ShellTraj.fncCalcVS_v_PrStSK(   arrVectPartialDeriv  );
		memcpy(&arrF_SSKT[6], &arrVectPartialDeriv[3], 3 * sizeof(double));

		double arrF_SSK[9] = {0.};
		MatrTransp(arrF_SSKT, 3, 3, arrF_SSK);

		double arrF_GSK[9] = {0.}, arrMtrxTransformOut[9] = {0.};
		ShellTraj.createMtrxTransform_xyzSSK_To_xyzGSK( 3,  arrMtrxTransformOut);
		MtrxMultMatrx(arrMtrxTransformOut, 3, 3,  arrF_SSK,3, arrF_GSK) ;



		arrF_GSK [2] -=  arrTargVS_GSK0[3];
		arrF_GSK [5] -=  arrTargVS_GSK0[4];
		arrF_GSK [8] -=  arrTargVS_GSK0[5];

		double arrFInv [9] = {0.};
		InverseMtrx3(arrF_GSK, arrFInv);

		double ar_x[3] ={0.}, ar_x0[3] ={0.};
		ar_x[0] = valBet0 ;
		ar_x[1] = valEps0 ;
		ar_x[2] = ShellTraj.mTCur;

		double arrtEmp[3] = {0.};
		MtrxMultMatrx(arrFInv,3, 3, arrDElta,1, arrtEmp) ;
		MtrxMinusMatrx(ar_x , arrtEmp,3, 1, ar_x0);
		valBet0 = ar_x0[0];
		valEps0 = ar_x0[1];



		if (valF0 < dR)
		{
		breturn = true;
		break;
		}



	}


	if (!breturn)
	{
	return false;
	}

	if(valBet0 < 0.0) valBet0 +=  2.0*M_PI ;
	if(valBet0 >=  2.0*M_PI) valBet0 -=  2.0*M_PI ;

	*valKGSKEps = valEps0;
	*valKGSKBet = valBet0;
	*valTFlight =ShellTraj.mTCur;
	memcpy(arrVectAppointmentPointGSK,arrTargGSKPosCur, 3 * sizeof(double));
	*valMiss = valF0;

	memcpy(arrTargVeloAppointmPoint_GSK, &arrTargVS_GSK0[3], 3 * sizeof(double));

	 double arrVS_PrStSK[6] = {0.};
	ShellTraj.fncCalcVS_v_PrStSK(  arrVS_PrStSK  ) ;
	ShellTraj.transform_xyzSSK_To_xyzGSK( 3, &arrVS_PrStSK[3], arrShellVeloAppointmPoint_GSK) ;
	return true;
}

//-----------------------------------------------------------------------
// передвигание всего класса на время  VAlTime > tFilt
void TFight::shift(double VAlTime)
{
	int iFiltSteps =  (VAlTime - mT)/ mVessel.mControlSyst.mFiltT;
	double valT0 = mT;
	for (int i = 0; i < iFiltSteps ; i++)
	{
	double valTCur = valT0 + (double(i + 1.) ) * mVessel.mControlSyst.mFiltT ;

	recalcFight( valTCur,valTCur + mVessel.mControlSyst.mFiltT);
	}

}

//------------------------------------------------------------------------------------------------------------------
// Вычисление вероятности поражения в заданной точке траектории
// INPUT:
// VAlTCurrentShot - время выстрела от начвала движения
// VAlTFlight   - полетное время
// VAlKGSKEps, VAlKGSKBet - углы наведения в КГСК
// OUTPUT:
// *pvalProb   - вероятность
// *pvalDispMiss - дисперсия промаха
//
//
void TFight::calcProbability_For_Fixed_AppointmentPoint_AirTargs(const double VAlTCurrentShot
, const double VAlTFlight ,const double VAlKGSKEps, const double VAlKGSKBet,double *pvalProb
,double  *pvalDispMiss,double *pvalDispNedolet
	 ,double *arrCorMtrxCartinSK, double *pvalDisp_GSK_Z, TURPolygon *pPlgProjection)
{
   TFight FightCur = *this;
	FightCur.mVessel.mMaxQ = 0.; // корабль двигается по идеальной прямой

 //	memset( FightCur.mTarget.mTraject.marrSigW, 0, 3 * sizeof(double)); // цель двигается по идеальной прямой





	// передвижение члена класса FightCur на время  VAlTCurrentShot
	// поскольку передвижение совершается дискретными шагами кратными
	// темпу фильтрации, то текущее время FightCur вообще гворя может быть меньше  VAlTCurrentShot
	//double val0 = FightCur.mTraceFlt.marrK_KGSK[0];

	FightCur.shift(VAlTCurrentShot);

	///


	// экстраполяция вектиора состояния цели на время  (VAlTCurrentShot  - FightCur.mT) вперед
	double arrTargExtrapVS_GSK[9] = {0.};// это вектор состояния цели на момент выстрела  VAlTCurrentShot
	FightCur.mTarget.mTraject.extrapolateTargVS(VAlTCurrentShot - FightCur.mT, arrTargExtrapVS_GSK);
	///

	// экстраполяция вектиора состояния корабля на время  (valTCurrentFight - FightCur.mT) вперед
	double arrVessExtrapVS_GSK[9] ={0.};
	FightCur.mVessel.extrapolateTrueVS_GSK(VAlTCurrentShot - FightCur.mT, arrVessExtrapVS_GSK);
	///


// формирование диагональной матрицы LEN_ARR_SCATTERS x LEN_ARR_SCATTERS дисперсий возмущений
//  - 8 тначаольных условий, коэфф формы, масса, ветер по X, ветер по Y
	double arrMtrxShellDisp [LEN_ARR_SCATTERS * LEN_ARR_SCATTERS] = {0.};
	fillShellVozmDispMatr (false, VAlKGSKEps,  VAlKGSKBet, arrMtrxShellDisp);
	///

	// формирование матрицы частных производных фазового вектора снаряда в точке встречи по начальным условиям Пси, Пи, Тетта, Масса, Cx? Cz, 3 параметра ветра - 9 штук
	// массе, коэффициенту формы, горизонтальному ветру по осям X и Z в ССК , коэффициенту формы по оси Z- всего 13 штук
	double arrStrSK_Jacobian[8 * LEN_ARR_SCATTERS] = {0.}, arrShellScatteringsCorMtrx_GSK[36] = {0.};
	double arrShellVS_GSK[6] = {0.}; // ВС снаряда в ГСК
	TMyShellTraj ShellTraj (&arrVessExtrapVS_GSK[3], mVessel.mShellBody,  VAlKGSKEps,VAlKGSKBet );

	const double VAlStepInt = 0.0005;
	double arrShellScatteringsCorMtrxPos_SSK[9] ={0.};

	// вычисление корреляционной матрицы разброса вектора состояния снаряда
// в ГСК
// INPUT:
// Environment -внешняя среда
// VAlStepInt - шаг интегрирования
// VAlTFlight - полетное время
// arrMtrxShellDisp  - диагональная матрица диспервсий разбросов параметров
//OUTPUT:
// arrStrSK_Jacobian [LEN_ARR_SCATTERS * 8] - матрица частных производных фазового вектора
// снаряда в траектоореной сиситеме координат по параметрам
// arrShellScatteringsCorMtarx_GSK [ 6 * 6] - корреляционная матрица в ГСК
// arrShellScatteringsCorMtrxPos_SSK [3*3] - корреляцилнная матрица ошибак рассеяния положения в ССК
	ShellTraj.calc_VS_GSK_And_ScatteringsCorrMatrx_GSK (mEnvironment, VAlStepInt
	,VAlTFlight,arrMtrxShellDisp, arrStrSK_Jacobian, arrShellScatteringsCorMtrx_GSK, arrShellVS_GSK, arrShellScatteringsCorMtrxPos_SSK);

	///
	// вычисление корреляционной матрицы разбросов вектора состояния цели на момент  VAlTCurrentShot + VAlTFlight
	double arrTargExtrpScatteringsCorMtrx_GSK[36] = {0.}, arrTargExtrpVS_KGSK[6] = {0.};
	 FightCur.extrapolateTrueTargVS_And_ScatteringCorrMatrx_KGSK( VAlTFlight
		 , arrTargExtrpScatteringsCorMtrx_GSK, arrTargExtrpVS_KGSK);

	//

	// формирование корреляционной матрицы разбросов вект ора промаха
	double arrVectMiss[6] = {0.}, arrCorrMatrxVectMiss [36] ={0.};
	MtrxSumMatrx( arrTargExtrpScatteringsCorMtrx_GSK, arrShellScatteringsCorMtrx_GSK,6, 6, arrCorrMatrxVectMiss) ;
	*pvalDisp_GSK_Z = sqrt(arrCorrMatrxVectMiss[ 6 * 2 + 2]);
	MtrxMinusMatrx(&arrShellVS_GSK[3], &arrTargExtrapVS_GSK [3], 3, 1, &arrVectMiss[3]);
 //	memcpy(&arrVectMiss[3], &arrTargExtrpVS_KGSK[3], 3 * sizeof(double));
	///


		// вычисление дальности точки встречи
	double arrTrueTargPosKGSK[3] = {0.};
	calcTargPosKGSK(VAlTCurrentShot + VAlTFlight,arrTrueTargPosKGSK);
	double valDistAppPoint = Norm3(arrTrueTargPosKGSK);
	//
	///

	// вычисление вероятности
	double valMissSimulated =0.;
	double valPereletSimulated = mVessel.mControlSyst.mTblPereletSimulated.calcValue(valDistAppPoint);
      TNeighbourhoodAppPoint NeighbourhoodAppPoint(arrVectMiss, arrCorrMatrxVectMiss,  mTarget
	 , mVessel.mShellBody, valDistAppPoint, mVessel.mArtComplex.mArtCannon.mSigRabT, VAlTFlight
	 ,valMissSimulated,valPereletSimulated);

	 NeighbourhoodAppPoint.calcDestructionProb_For_AirTargs(pvalProb,pvalDispMiss
	 ,arrCorMtrxCartinSK,pvalDispNedolet,  pPlgProjection);

}


//------------------------------------------------------------------------------------------------------------------
// Вычисление вероятности поражения в заданной точке траектории  ПЕРЕГРУЖЕННАЯ!!!
// INPUT:
// VAlTCurrentShot - время выстрела от начвала движения
// VAlTFlight   - полетное время
// VAlKGSKEps, VAlKGSKBet - углы наведения в КГСК
// OUTPUT:
// *pvalProb   - вероятность
// *pvalDispMiss - дисперсия промаха
//
//
void TFight::calcProbability_For_Fixed_AppointmentPoint_AirTargs(const double VAlTCurrentShot
 ,const double VAlTFlight,const double VAlKGSKEps, const double VAlKGSKBet,double *pvalProb
 ,double  *pvalDispMiss,double *pvalDispNedolet,double *arrCorMtrxCartinSK
 ,double *pvalDisp_GSK_Z, TURPolygon *pPlgProjection, double *arrKTarg, double *arrKShell
 , double *arrMiss)
{
	TFight FightCur = *this;
	FightCur.mVessel.mMaxQ = 0.; // корабль двигается по идеальной прямой

 //	memset( FightCur.mTarget.mTraject.marrSigW, 0, 3 * sizeof(double)); // цель двигается по идеальной прямой





	// передвижение члена класса FightCur на время  VAlTCurrentShot
	// поскольку передвижение совершается дискретными шагами кратными
	// темпу фильтрации, то текущее время FightCur вообще гворя может быть меньше  VAlTCurrentShot


	FightCur.shift(VAlTCurrentShot);

	///

	// осреднение корреляцуионной матрицы ошибок оценивания вектора состояния цели
	// по начальным фазам палубных углов
	double arrKTarg_GSK_Middle[36]={0.};
	double arrTargExtrpScatteringsCorMtrx_GSK[36] = {0.}, arrTargExtrpVS_KGSK[6] = {0.};
	 FightCur.extrapolateTrueTargVS_And_ScatteringCorrMatrx_KGSK(VAlTCurrentShot + VAlTFlight
		 , arrTargExtrpScatteringsCorMtrx_GSK, arrTargExtrpVS_KGSK);
	int nC0 = 100;
	  for (int j = 0; j < nC0; j++)
	  {
		TFight FightTemp = *this;
		FightTemp.mVessel.mMaxQ = 0.; // корабль двигается по идеальной прямой
		for (int k =0; k < 4; k++)

		{
		FightTemp.mVessel.marrDelt [k] = getRand01() * 2. * M_PI;
		}

		FightTemp.shift(VAlTCurrentShot);

		// вычисление корреляционной матрицы разбросов вектора состояния цели на момент  VAlTCurrentShot + VAlTFlight
	double arrTargExtrpScatteringsCorMtrx_GSK0[36] = {0.}, arrTargExtrpVS_KGSK0[6] = {0.};
	 FightTemp.extrapolateTrueTargVS_And_ScatteringCorrMatrx_KGSK(VAlTCurrentShot + VAlTFlight
		 , arrTargExtrpScatteringsCorMtrx_GSK0, arrTargExtrpVS_KGSK0);
		 MtrxSumMatrx(arrKTarg_GSK_Middle, arrTargExtrpScatteringsCorMtrx_GSK0,6, 6, arrKTarg_GSK_Middle) ;

	  }
	  MatrxDivideScalar(arrKTarg_GSK_Middle, 6, 6, ((double)nC0),arrTargExtrpScatteringsCorMtrx_GSK);
	 memcpy(arrKTarg, arrTargExtrpScatteringsCorMtrx_GSK , 36 * sizeof(double));
	 ///

	// экстраполяция вектиора состояния цели на время  (VAlTCurrentShot  - FightCur.mT) вперед
	double arrTargExtrapVS_GSK[9] = {0.};// это вектор состояния цели на момент выстрела  VAlTCurrentShot
	FightCur.mTarget.mTraject.extrapolateTargVS(VAlTCurrentShot - FightCur.mT, arrTargExtrapVS_GSK);
	///

	// экстраполяция вектиора состояния корабля на время  (valTCurrentFight - FightCur.mT) вперед
	double arrVessExtrapVS_GSK[9] ={0.};
	FightCur.mVessel.extrapolateTrueVS_GSK(VAlTCurrentShot - FightCur.mT, arrVessExtrapVS_GSK);
	///


// формирование диагональной матрицы LEN_ARR_SCATTERS x LEN_ARR_SCATTERS дисперсий возмущений
//  - 8 тначаольных условий, коэфф формы, масса, ветер по X, ветер по Y
	double arrMtrxShellDisp [LEN_ARR_SCATTERS * LEN_ARR_SCATTERS] = {0.};
	fillShellVozmDispMatr (false, VAlKGSKEps,  VAlKGSKBet, arrMtrxShellDisp);
	///

	// формирование матрицы частных производных фазового вектора снаряда в точке встречи по начальным условиям Пси, Пи, Тетта, Масса, Cx? Cz, 3 параметра ветра - 9 штук
	// массе, коэффициенту формы, горизонтальному ветру по осям X и Z в ССК , коэффициенту формы по оси Z- всего 13 штук
	double arrStrSK_Jacobian[8 * LEN_ARR_SCATTERS] = {0.}, arrShellScatteringsCorMtrx_GSK[36] = {0.};
	double arrShellVS_GSK[6] = {0.}; // ВС снаряда в ГСК
	TMyShellTraj ShellTraj (&arrVessExtrapVS_GSK[3], mVessel.mShellBody,  VAlKGSKEps,VAlKGSKBet );

	//const double VAlStepInt = 0.0005;
	const double VAlStepInt = 0.0001;
	double arrShellScatteringsCorMtrxPos_SSK[9] ={0.};

	// вычисление корреляционной матрицы разброса вектора состояния снаряда
// в ГСК
// INPUT:
// Environment -внешняя среда
// VAlStepInt - шаг интегрирования
// VAlTFlight - полетное время
// arrMtrxShellDisp  - диагональная матрица диспервсий разбросов параметров
//OUTPUT:
// arrStrSK_Jacobian [LEN_ARR_SCATTERS * 8] - матрица частных производных фазового вектора
// снаряда в траектоореной сиситеме координат по параметрам
// arrShellScatteringsCorMtarx_GSK [ 6 * 6] - корреляционная матрица в ГСК
// arrShellScatteringsCorMtrxPos_SSK [3*3] - корреляцилнная матрица ошибак рассеяния положения в ССК
	ShellTraj.calc_VS_GSK_And_ScatteringsCorrMatrx_GSK (mEnvironment, VAlStepInt
	,VAlTFlight,arrMtrxShellDisp, arrStrSK_Jacobian, arrShellScatteringsCorMtrx_GSK
	, arrShellVS_GSK, arrShellScatteringsCorMtrxPos_SSK);

	memcpy(arrKShell, arrShellScatteringsCorMtrx_GSK , 36 * sizeof(double));

	///

	/*
	// вычисление корреляционной матрицы разбросов вектора состояния цели на момент  VAlTCurrentShot + VAlTFlight
	double arrTargExtrpScatteringsCorMtrx_GSK[36] = {0.}, arrTargExtrpVS_KGSK[6] = {0.};
	 FightCur.extrapolateTrueTargVS_And_ScatteringCorrMatrx_KGSK(VAlTCurrentShot + VAlTFlight
		 , arrTargExtrpScatteringsCorMtrx_GSK, arrTargExtrpVS_KGSK);
   // ВНИМАНИЕ\ !!!
	memcpy(arrKTarg, arrTargExtrpScatteringsCorMtrx_GSK , 36 * sizeof(double));
	 */
		// вычисление дальности точки встречи
	double arrTrueTargPosKGSK[3] = {0.};
	calcTargPosKGSK(VAlTCurrentShot + VAlTFlight,arrTrueTargPosKGSK);
	double valDistAppPoint = Norm3(arrTrueTargPosKGSK);
	//

	// формирование корреляционной матрицы разбросов вект ора промаха
	double arrVectMiss[6] = {0.}, arrCorrMatrxVectMiss [36] ={0.};
	MtrxSumMatrx( arrTargExtrpScatteringsCorMtrx_GSK, arrShellScatteringsCorMtrx_GSK,6, 6, arrCorrMatrxVectMiss) ;
	*pvalDisp_GSK_Z = sqrt(arrCorrMatrxVectMiss[ 6 * 2 + 2]);
	MtrxMinusMatrx(&arrShellVS_GSK[3], &arrTargExtrapVS_GSK [3], 3, 1, &arrVectMiss[3]);
	memcpy(arrMiss, arrVectMiss, 6 * sizeof(double));
 //	memcpy(&arrVectMiss[3], &arrTargExtrpVS_KGSK[3], 3 * sizeof(double));
	///




	// вычисление вероятности
 //TNeighbourhoodAppPoint::calcDestructionProb_For_AirTargs(valDistAppPoint, arrVectMiss, arrCorrMatrxVectMiss
   //	 , mTarget, mVessel.mShellBody, pvalProb, pvalDispMiss);

	 double valMissSimulated = 0.;
	double valPereletSimulated = mVessel.mControlSyst.mTblPereletSimulated.calcValue(valDistAppPoint);
	 TNeighbourhoodAppPoint NeighbourhoodAppPoint(arrVectMiss, arrCorrMatrxVectMiss,  mTarget
	 , mVessel.mShellBody, valDistAppPoint, mVessel.mArtComplex.mArtCannon.mSigRabT, VAlTFlight
	 ,valMissSimulated,valPereletSimulated);

	 NeighbourhoodAppPoint.calcDestructionProb_For_AirTargs(pvalProb,pvalDispMiss
	 ,arrCorMtrxCartinSK,pvalDispNedolet,  pPlgProjection);

}
//------------------------------------------------------------------------------------------------------------------

/*

//------------------------------------------------------------------------------------------------------------------
// Вычисление вероятности поражения в заданной точке траектории  ПЕРЕГРУЖЕННАЯ!!!
// INPUT:
// VAlTCurrentShot - время выстрела от начвала движения
// VAlTFlight   - полетное время
// VAlKGSKEps, VAlKGSKBet - углы наведения в КГСК
// OUTPUT:
// *pvalProb   - вероятность
// *pvalDispMiss - дисперсия промаха
//
//
void TFight::calcProbability_For_Fixed_AppointmentPoint_AirTargs(const double VAlTCurrentShot
 ,const double VAlTFlight,const double VAlKGSKEps, const double VAlKGSKBet,double *pvalProb
 ,double  *pvalDispMiss,double *pvalDispNedolet,double *arrCorMtrxCartinSK
 ,double *pvalDisp_GSK_Z, TURPolygon *pPlgProjection, double *arrKTarg, double *arrKShell
 , double *arrMiss)
{
	TFight FightCur = *this;
	FightCur.mVessel.mMaxQ = 0.; // корабль двигается по идеальной прямой

 //	memset( FightCur.mTarget.mTraject.marrSigW, 0, 3 * sizeof(double)); // цель двигается по идеальной прямой





	// передвижение члена класса FightCur на время  VAlTCurrentShot
	// поскольку передвижение совершается дискретными шагами кратными
	// темпу фильтрации, то текущее время FightCur вообще гворя может быть меньше  VAlTCurrentShot
   //	double val0 = FightCur.mTraceFlt.marrK_KGSK[0];

	FightCur.shift(VAlTCurrentShot);

	///


	// экстраполяция вектиора состояния цели на время  (VAlTCurrentShot  - FightCur.mT) вперед
	double arrTargExtrapVS_GSK[9] = {0.};// это вектор состояния цели на момент выстрела  VAlTCurrentShot
	FightCur.mTarget.mTraject.extrapolateTargVS(VAlTCurrentShot - FightCur.mT, arrTargExtrapVS_GSK);
	///

	// экстраполяция вектиора состояния корабля на время  (valTCurrentFight - FightCur.mT) вперед
	double arrVessExtrapVS_GSK[9] ={0.};
	FightCur.mVessel.extrapolateTrueVS_GSK(VAlTCurrentShot - FightCur.mT, arrVessExtrapVS_GSK);
	///


// формирование диагональной матрицы LEN_ARR_SCATTERS x LEN_ARR_SCATTERS дисперсий возмущений
//  - 8 тначаольных условий, коэфф формы, масса, ветер по X, ветер по Y
	double arrMtrxShellDisp [LEN_ARR_SCATTERS * LEN_ARR_SCATTERS] = {0.};
	fillShellVozmDispMatr (false, VAlKGSKEps,  VAlKGSKBet, arrMtrxShellDisp);
	///

	// формирование матрицы частных производных фазового вектора снаряда в точке встречи по начальным условиям Пси, Пи, Тетта, Масса, Cx? Cz, 3 параметра ветра - 9 штук
	// массе, коэффициенту формы, горизонтальному ветру по осям X и Z в ССК , коэффициенту формы по оси Z- всего 13 штук
	double arrStrSK_Jacobian[8 * LEN_ARR_SCATTERS] = {0.}, arrShellScatteringsCorMtrx_GSK[36] = {0.};
	double arrShellVS_GSK[6] = {0.}; // ВС снаряда в ГСК
	TMyShellTraj ShellTraj (&arrVessExtrapVS_GSK[3], mVessel.mShellBody,  VAlKGSKEps,VAlKGSKBet );

	//const double VAlStepInt = 0.0005;
	const double VAlStepInt = 0.0001;
	double arrShellScatteringsCorMtrxPos_SSK[9] ={0.};

	// вычисление корреляционной матрицы разброса вектора состояния снаряда
// в ГСК
// INPUT:
// Environment -внешняя среда
// VAlStepInt - шаг интегрирования
// VAlTFlight - полетное время
// arrMtrxShellDisp  - диагональная матрица диспервсий разбросов параметров
//OUTPUT:
// arrStrSK_Jacobian [LEN_ARR_SCATTERS * 8] - матрица частных производных фазового вектора
// снаряда в траектоореной сиситеме координат по параметрам
// arrShellScatteringsCorMtarx_GSK [ 6 * 6] - корреляционная матрица в ГСК
// arrShellScatteringsCorMtrxPos_SSK [3*3] - корреляцилнная матрица ошибак рассеяния положения в ССК
	ShellTraj.calc_VS_GSK_And_ScatteringsCorrMatrx_GSK (mEnvironment, VAlStepInt
	,VAlTFlight,arrMtrxShellDisp, arrStrSK_Jacobian, arrShellScatteringsCorMtrx_GSK
	, arrShellVS_GSK, arrShellScatteringsCorMtrxPos_SSK);

	memcpy(arrKShell, arrShellScatteringsCorMtrx_GSK , 36 * sizeof(double));

	///
	// вычисление корреляционной матрицы разбросов вектора состояния цели на момент  VAlTCurrentShot + VAlTFlight
	double arrTargExtrpScatteringsCorMtrx_GSK[36] = {0.}, arrTargExtrpVS_KGSK[6] = {0.};
	 FightCur.extrapolateTrueTargVS_And_ScatteringCorrMatrx_KGSK(VAlTCurrentShot + VAlTFlight
  // FightCur.extrapolateTrueTargVS_And_ScatteringCorrMatrx_KGSK( VAlTFlight
		 , arrTargExtrpScatteringsCorMtrx_GSK, arrTargExtrpVS_KGSK);
   // ВНИМАНИЕ\ !!!
	memcpy(arrKTarg, arrTargExtrpScatteringsCorMtrx_GSK , 36 * sizeof(double));

		// вычисление дальности точки встречи
	double arrTrueTargPosKGSK[3] = {0.};
	calcTargPosKGSK(VAlTCurrentShot + VAlTFlight,arrTrueTargPosKGSK);
	double valDistAppPoint = Norm3(arrTrueTargPosKGSK);
	//

	// формирование корреляционной матрицы разбросов вект ора промаха
	double arrVectMiss[6] = {0.}, arrCorrMatrxVectMiss [36] ={0.};
	MtrxSumMatrx( arrTargExtrpScatteringsCorMtrx_GSK, arrShellScatteringsCorMtrx_GSK,6, 6, arrCorrMatrxVectMiss) ;
	*pvalDisp_GSK_Z = sqrt(arrCorrMatrxVectMiss[ 6 * 2 + 2]);
	MtrxMinusMatrx(&arrShellVS_GSK[3], &arrTargExtrapVS_GSK [3], 3, 1, &arrVectMiss[3]);
	memcpy(arrMiss, arrVectMiss, 6 * sizeof(double));
 //	memcpy(&arrVectMiss[3], &arrTargExtrpVS_KGSK[3], 3 * sizeof(double));
	///




	// вычисление вероятности
 //TNeighbourhoodAppPoint::calcDestructionProb_For_AirTargs(valDistAppPoint, arrVectMiss, arrCorrMatrxVectMiss
   //	 , mTarget, mVessel.mShellBody, pvalProb, pvalDispMiss);

	 double valMissSimulated = 0.;
	double valPereletSimulated = mVessel.mControlSyst.mTblPereletSimulated.calcValue(valDistAppPoint);
	 TNeighbourhoodAppPoint NeighbourhoodAppPoint(arrVectMiss, arrCorrMatrxVectMiss,  mTarget
	 , mVessel.mShellBody, valDistAppPoint, mVessel.mArtComplex.mArtCannon.mSigRabT, VAlTFlight
	 ,valMissSimulated,valPereletSimulated);

	 NeighbourhoodAppPoint.calcDestructionProb_For_AirTargs(pvalProb,pvalDispMiss
	 ,arrCorMtrxCartinSK,pvalDispNedolet,  pPlgProjection);

}
//------------------------------------------------------------------------------------------------------------------
*/
/*
// Вычисление вероятности поражения в заданной точке траектории
// INPUT:
// VAlTCurrentShot - время выстрела от начвала движения
// VAlTFlight   - полетное время
// VAlKGSKEps, VAlKGSKBet - углы наведения в КГСК
// OUTPUT:
// *pvalProb   - вероятность
// *pvalDispMiss - дисперсия промаха
//
//
void TFight::calcProbability_For_Fixed_AppointmentPoint_SeaTarg(const double VAlTCurrentShot, const double VAlTFlight
	 ,const double VAlKGSKEps, const double VAlKGSKBet, double *pvalProb, double *pvalDispMiss)
{
	double arrCoMtrx_1_and_3_Groups[36] = {	0.}, arrCoMtrx_2_Group[36] = {0.}
				, arrCorrMatrxVectMiss [36 ]= {0.}, arrShellVS_GSK [6] = {0.};
		calcCorMtrx_First_And_Third_Group(VAlKGSKEps, VAlKGSKBet, VAlTFlight,arrCoMtrx_1_and_3_Groups, arrShellVS_GSK);
		const double VAlObservTime = 10.;
		calcCorMtrx_Second_Group(VAlObservTime, VAlTFlight , arrCoMtrx_2_Group);
		MtrxSumMatrx(arrCoMtrx_1_and_3_Groups, arrCoMtrx_2_Group, 6, 6, arrCorrMatrxVectMiss);

			double arrVectMiss[6] = {0.};

		MtrxMinusMatrx(&arrShellVS_GSK[3], &(mTarget.mTraject.marrVectSostGSK_Begin[3]), 3, 1, &arrVectMiss[3]);

	// экстраполяция вектиора состояния цели на время  (VAlTCurrentShot  - FightCur.mT) вперед
	double arrTargExtrapVS_GSK[9] = {0.};// это вектор состояния цели на момент выстрела  VAlTCurrentShot
	FightCur.mTarget.mTraject.extrapolateTargVS(VAlTCurrentShot - FightCur.mT, arrTargExtrapVS_GSK);
	///

	// экстраполяция вектиора состояния корабля на время  (valTCurrentFight - FightCur.mT) вперед
	double arrVessExtrapVS_GSK[9] ={0.};
	FightCur.mVessel.extrapolateTrueVS_GSK(VAlTCurrentShot - FightCur.mT, arrVessExtrapVS_GSK);
	///


// формирование диагональной матрицы LEN_ARR_SCATTERS x LEN_ARR_SCATTERS дисперсий возмущений - 8 тначаольных условий, коэфф формы, масса, ветер по X, ветер по Y
	double arrMtrxShellDisp [LEN_ARR_SCATTERS * LEN_ARR_SCATTERS] = {0.};
	fillShellVozmDispMatr (VAlKGSKEps,  VAlKGSKBet, arrMtrxShellDisp);
	///

	// формирование матрицы частных производных фазового вектора снаряда в точке встречи по начальным условиям Пси, Пи, Тетта, Масса, Cx? Cz, 3 параметра ветра - 9 штук
	// массе, коэффициенту формы, горизонтальному ветру по осям X и Z в ССК , коэффициенту формы по оси Z- всего 13 штук
	double arrStrSK_Jacobian[8 * LEN_ARR_SCATTERS] = {0.}, arrShellScatteringsCorMtrx_GSK[36] = {0.};
	double arrShellVS_GSK[6] = {0.}; // ВС снаряда в ГСК
	TMyShellTraj ShellTraj (&arrVessExtrapVS_GSK[3], mVessel.mShellBody,  VAlKGSKEps,VAlKGSKBet );

	const double VAlStepInt = 0.0005;
	double arrShellScatteringsCorMtrxPos_SSK[9] ={0.};
	ShellTraj.calc_VS_GSK_And_ScatteringsCorrMatrx_GSK (mEnvironment, VAlStepInt
	,VAlTFlight,arrMtrxShellDisp, arrStrSK_Jacobian, arrShellScatteringsCorMtrx_GSK, arrShellVS_GSK, arrShellScatteringsCorMtrxPos_SSK);

	///
	// вычисление корреляционной матрицы разбросов вектора состояния цели на момент  VAlTCurrentShot + VAlTFlight
	double arrTargExtrpScatteringsCorMtrx_GSK[36] = {0.}, arrTargExtrpVS_KGSK[6] = {0.};
	 FightCur.extrapolateTrueTargVS_And_ScatteringCorrMatrx_KGSK(VAlTCurrentShot + VAlTFlight
		 , arrTargExtrpScatteringsCorMtrx_GSK, arrTargExtrpVS_KGSK);
	//

	// формирование корреляционной матрицы разбросов вект ора промаха
	double arrVectMiss[6] = {0.}, arrCorrMatrxVectMiss [36] ={0.};
	MtrxSumMatrx( arrTargExtrpScatteringsCorMtrx_GSK, arrShellScatteringsCorMtrx_GSK,6, 6, arrCorrMatrxVectMiss) ;
		MtrxMinusMatrx(&arrShellVS_GSK[3], &arrTargExtrapVS_GSK [3], 3, 1, &arrVectMiss[3]);
 //	memcpy(&arrVectMiss[3], &arrTargExtrpVS_KGSK[3], 3 * sizeof(double));
	///

	// вычисление дальности точки встречи
   //	double arrAppPoint[3] ={0.};
	 //	MtrxMinusMatrx(arrTargExtrpVS_GSK, arrShellVS_GSK, 1, 3, arrAppPoint);

	 //	double valDistAppPoint = Norm3(arrAppPoint);
	const double VAlDistAppPoint = Norm3(mTarget.mTraject.marrVectSostGSK_Begin);
	///

	// вычисление вероятности
	TNeighbourhoodAppPoint NeighbourhoodAppPoint(arrVectMiss, arrCorrMatrxVectMiss, mTarget, mVessel.mShellBody, VAlDistAppPoint);
	NeighbourhoodAppPoint.calcProbDirectHitting_SeaTarg(pvalProb, pvalDispMiss);

}
*/
//------------------------------------------------------------------------------------------------------------------
void TFight::extrapolateTrueTargVS_And_ScatteringCorrMatrx_KGSK(const double VAlTExtrp
			 , double *arrTargExtrpScatteringsCorMtrx_GSK,  double *arrTargExtrpVS_KGSK)
{
	 //корреляционная матрица  mTraceFlt.marrK_KGSK не учитвает ошибок измерений палубных углов, так как
	 // рассчитывалась для нужд сопровождения
	 // надо добавить разбросы палубных углов
	 // 1. вычисление вектора состояния цели в ПСК
	 double arrTargVS_PSK[9] = {0.};
	 MtrxMultMatrx(mTraceFlt.marrFExt_ASK_V_PSK,6, 6, mTraceFlt.marrVS_ASK,1, arrTargVS_PSK) ;
	 ///

	 // 2. исление дисперсий ошибок измерения углов с учетом задержки СИНС
	 double valDispDEltaExtrIS_Psi = 0., valDispDEltaExtrIS_Tet = 0., valDispDEltaExtrIS_Q = 0.;
	 double  valDispDEltaExtrDef_Psi = 0., valDispDEltaExtrDef_Tet = 0.;
	 double valDispDEltaExtr_Psi = 0., valDispDEltaExtr_Tet = 0., valDispDEltaExtr_Q = 0.;
	 mVessel.calcDispExtrapDeckAngles_SINS(mVessel.mControlSyst.mSinsDelayT   , &valDispDEltaExtrIS_Psi
	 ,&valDispDEltaExtrIS_Tet, &valDispDEltaExtrIS_Q);

	mVessel.calcDispExtrapDeckAngles_Deform(&valDispDEltaExtrDef_Psi,&valDispDEltaExtrDef_Tet);
	 valDispDEltaExtr_Psi =  valDispDEltaExtrIS_Psi +  valDispDEltaExtrDef_Psi;
	 valDispDEltaExtr_Tet =  valDispDEltaExtrIS_Tet +  valDispDEltaExtrDef_Tet;
	 valDispDEltaExtr_Q = valDispDEltaExtrIS_Q;

	 // 3. Вычисление корреляционной матрицы разбросов вектора положения и скорости
	 double arrCorrPosition[9] = {0.}, arrCorrVelocity[9] = {0.}, arrT00[3] = {0.};
	 calcCorMatrx_PSK_KGSK(mVessel.mQ0,0.,0., valDispDEltaExtr_Q, valDispDEltaExtr_Psi
	 , valDispDEltaExtr_Tet, arrTargVS_PSK, arrT00, arrCorrPosition);

	 calcCorMatrx_PSK_KGSK(mVessel.mQ0,0.,0., valDispDEltaExtr_Q, valDispDEltaExtr_Psi
	 , valDispDEltaExtr_Tet, &arrTargVS_PSK[3], arrT00, arrCorrVelocity);

	 // 4. Формирование результирующей корр матрицы 6х6
	 double arrMtrxCor_Q_Psi_Tet[36]= {0.};
	 for (int i =0; i < 3; i++)
	 {
	   for (int j = 0; j < 3; j++)
	   {
		 arrMtrxCor_Q_Psi_Tet [ 6 * i + j] =  arrCorrPosition[ 3 * i + j];
		 arrMtrxCor_Q_Psi_Tet [ 6 * ( 3 + i) + 3 + j] =  arrCorrVelocity[ 3 * i + j];
	   }
	 }
	 // 5.
	 double arrTarg_KGSK [36] = {0.};
	 MtrxSumMatrx(arrTarg_KGSK, mTraceFlt.marrK_KGSK,6, 6, arrTarg_KGSK) ;
		///////////////////////////////////////////////////////////////////////
	 const double VAlDeltaTExtr = VAlTExtrp - mTraceFlt.mTf;
	 double arrL[36] ={0.}, arrTemp [36] ={0.}, arrTemp0 [36] ={0.};
	 CreateMatrL6 (VAlDeltaTExtr, arrL) ;
	 MtrxMultMatrx(arrL, 6, 6, arrTarg_KGSK, 6, arrTemp) ;
	 MtrxMultMatrxTransp(arrTemp,6, 6, arrL,6, arrTemp0) ;

	double arrF [12] ={0.}, arrFluctAdd[36] = {0.};
	arrF [0] = mTarget.mTraject.marrSigW[0] * mTarget.mTraject.marrSigW[0] * VAlDeltaTExtr * VAlDeltaTExtr * VAlDeltaTExtr /3. ;
	arrF [ 1] =  mTarget.mTraject.marrSigW[0] * mTarget.mTraject.marrSigW[0] * VAlDeltaTExtr * VAlDeltaTExtr  /2. ;
	arrF [ 2] =arrF [ 1] ;
	arrF [ 3] = mTarget.mTraject.marrSigW[0] * mTarget.mTraject.marrSigW[0] * VAlDeltaTExtr ;
	memcpy(&arrF[4], arrF, 4 * sizeof(double));
	memcpy(&arrF[8], arrF, 4 * sizeof(double));
	IntegrateMatrK(arrF , arrFluctAdd )   ;
	MtrxSumMatrx(arrTemp0, arrFluctAdd,6, 6, arrTargExtrpScatteringsCorMtrx_GSK) ;


   MtrxMultMatrx(arrL, 6, 6, mTraceFlt.marrVS_KGSK, 1, arrTargExtrpVS_KGSK) ;
   MtrxSumMatrx(arrMtrxCor_Q_Psi_Tet, arrTargExtrpScatteringsCorMtrx_GSK,6, 6, arrTargExtrpScatteringsCorMtrx_GSK) ;   //!!!!!! 09.04.2018
}

//------------------------------------------------------------------------------------------------------------------
 /*
// формирование 10х10 диагональной матрицы возмущений для снаряда
void TFight::fillShellVozmDispMatr (const double VAlKGSKEps, const double VAlKGSKBet, double *arrMtrxShellDisp)
{
	double  valDispRzvT = TProbabilityTheory::calcDispRavnomern(mVessel.mControlSyst.mRzvT);
	double valTExtrap =  mVessel.mControlSyst.mSinsDelayT + sqrt(valDispRzvT);
	// фазовый  вектор  в стрельбовой СК:

// 0. marrStrSK_VS [3]-  угол Пси

// 1. marrStrSK_VS [5]-  относит плотность атмосферы Пи



// 2.  marrStrSK_VS [7]-  угол Тетта
// 3. marrStrSK_VS [6]-  путевая скорость
// дополнительные разбросы:
// 4.  масса
// 5.  коэфф формы Cx
// 6. коэф формы по оси Z
// 7. модуль горизонтальной скорости ветра
// 8. направление горизонтального ветра
// 9. модуль вертикального ветра
 double valDispDeltaPsi = 0., valDispDeltaTet = 0., valDispDeltaQ = 0. ;
 double valDispDEltaExtrIS_Psi = 0., valDispDEltaExtrIS_Tet = 0., valDispDEltaExtrIS_Q = 0. ; // дисперсии ошибок экстраполяции углов СИНС
 double valDispDEltaExtrDeform_Psi = 0., valDispDEltaExtrDeform_Tet = 0. ; // дисперсии ошибок вызванных деформацикй корпуса корабля
 // вычисление дисперсий ишибок экстраполяцмии палубных углов СИНС
 mVessel.calcDispExtrapDeckAngles_SINS(valTExtrap, &valDispDEltaExtrIS_Psi, &valDispDEltaExtrIS_Tet, &valDispDEltaExtrIS_Q);
 mVessel.calcDispExtrapDeckAngles_Deform(valTExtrap, &valDispDEltaExtrDeform_Psi, &valDispDEltaExtrDeform_Tet);

 valDispDeltaPsi  =  valDispDEltaExtrIS_Psi + valDispDEltaExtrDeform_Psi ;
 valDispDeltaTet  =  valDispDEltaExtrIS_Tet + valDispDEltaExtrDeform_Tet ;
 valDispDeltaQ  =  valDispDEltaExtrIS_Q ;


 // вычисление дисперсий ошибок оринетации углов наведения Eps и Bet без учета кучности
 double valDispEps = 0., valDispBet = 0.;
 calcCorMtrx_OshibokOtrabotki_EpsBet(VAlKGSKEps, VAlKGSKBet,valDispDeltaPsi, valDispDeltaTet, valDispDeltaQ, &valDispBet , &valDispEps);
 valDispBet += mVessel.mArtComplex.mArtCannon.mAngGroupedFire * mVessel.mArtComplex.mArtCannon.mAngGroupedFire;
 valDispEps += mVessel.mArtComplex.mArtCannon.mAngGroupedFire * mVessel.mArtComplex.mArtCannon.mAngGroupedFire;

 arrMtrxShellDisp[0] =  valDispBet + mVessel.mArtComplex.mSigU * mVessel.mArtComplex.mSigU; //путевой угол
 arrMtrxShellDisp[ 2 * LEN_ARR_SCATTERS + 2] =  valDispEps + mVessel.mArtComplex.mSigU * mVessel.mArtComplex.mSigU; // угол наклона траетктории


 arrMtrxShellDisp[  LEN_ARR_SCATTERS +1 ] = 0.0001;// плотн атмосфкеры

 arrMtrxShellDisp[ 3 * LEN_ARR_SCATTERS + 3] = mVessel.mShellBody.mDispV0; // путевая скорость
 arrMtrxShellDisp[ 4 * LEN_ARR_SCATTERS + 4] = mVessel.mShellBody.mDispMass0; // масса
 arrMtrxShellDisp[ 5 * LEN_ARR_SCATTERS + 5] = mVessel.mShellBody.mDispCx; // коэф Cx
 arrMtrxShellDisp[ 6 * LEN_ARR_SCATTERS + 6] = mVessel.mShellBody.mDispCz;//клэф Cz
 arrMtrxShellDisp[ 7 * LEN_ARR_SCATTERS + 7] =  (mEnvironment.mWind_V * 0.1 ) * (mEnvironment.mWind_V * 0.1); // гориз ветер модуль
 arrMtrxShellDisp[ 8 * LEN_ARR_SCATTERS + 8] = (5. / 180. * M_PI) * (5. / 180. * M_PI);
 arrMtrxShellDisp[ 9 * LEN_ARR_SCATTERS + 9] =  (mEnvironment.mWind_VertV * 0.1 ) * (mEnvironment.mWind_VertV * 0.1); // вертик ветер модуль


}
*/

//-------------------------------------------------------------------------

// формирование 10х10 диагональной матрицы возмущений для снаряда
// bCOrrection - признак учета сиситематической ошибки (метео и разброса атмосферы)
void TFight::fillShellVozmDispMatr (const bool bCOrrection, const double VAlKGSKEps
	 , const double VAlKGSKBet, double *arrMtrxShellDisp)
{
	double valTExtrap =  mVessel.mControlSyst.mSinsDelayT/2.
	 +  mVessel.mArtComplex.mArtCannon.mRabT;


	// фазовый  вектор  в стрельбовой СК:

// 0. marrStrSK_VS [3]-  угол Пси

// 1. marrStrSK_VS [5]-  относит плотность атмосферы Пи



// 2.  marrStrSK_VS [7]-  угол Тетта
// 3. marrStrSK_VS [6]-  путевая скорость
// дополнительные разбросы:
// 4.  масса
// 5.  коэфф формы Cx
// 6. коэф формы по оси Z
// 7. модуль горизонтальной скорости ветра
// 8. направление горизонтального ветра
// 9. модуль вертикального ветра
 memset(arrMtrxShellDisp, 0,  sizeof(double) * LEN_ARR_SCATTERS * LEN_ARR_SCATTERS);

 double valDispDeltaPsiExtr = 0., valDispDeltaTetExtr = 0., valDispDeltaQExtr = 0. ;
 double valDispDEltaExtrIS_Psi = 0., valDispDEltaExtrIS_Tet = 0.
 , valDispDEltaExtrIS_Q = 0. ; // дисперсии ошибок экстраполяции углов СИНС
 double valDispDEltaExtrDeform_Psi = 0., valDispDEltaExtrDeform_Tet = 0. ; // дисперсии ошибок вызванных деформацикй корпуса корабля
 // вычисление дисперсий ишибок экстраполяцмии палубных углов СИНС
 mVessel.calcDispExtrapDeckAngles_SINS(valTExtrap, &valDispDEltaExtrIS_Psi, &valDispDEltaExtrIS_Tet, &valDispDEltaExtrIS_Q);
 mVessel.calcDispExtrapDeckAngles_Deform( &valDispDEltaExtrDeform_Psi, &valDispDEltaExtrDeform_Tet);

 valDispDeltaPsiExtr  =  valDispDEltaExtrIS_Psi + valDispDEltaExtrDeform_Psi ;
 valDispDeltaTetExtr  =  valDispDEltaExtrIS_Tet + valDispDEltaExtrDeform_Tet ;
 valDispDeltaQExtr  =  valDispDEltaExtrIS_Q ;


 // вычисление дисперсий ошибок оринетации углов наведения Eps и Bet без учета кучности
 double valDispEps = 0., valDispBet = 0.;

 double VAlPsi = 5./ 180. * M_PI; // среднее значение угла корм крена
 double VAlTet = 10./ 180. * M_PI; // среднее значение угла бортового крена
 double arrK_DelBetDelEpsGSK[4] = {0.};
 calcCorMtrx_OshibokOtrabotki_EpsBet(VAlKGSKEps, VAlKGSKBet
  ,valDispDeltaPsiExtr,valDispDeltaTetExtr, valDispDeltaQExtr
  ,VAlPsi ,VAlTet,mVessel.mQ0, arrK_DelBetDelEpsGSK);
  valDispBet = arrK_DelBetDelEpsGSK [0];
  valDispEps = arrK_DelBetDelEpsGSK [3];


 // добавление ошибок кучности (технического рассеяния АУ)
 valDispBet += mVessel.mArtComplex.mArtCannon.mAngGroupedFire * mVessel.mArtComplex.mArtCannon.mAngGroupedFire;
 valDispEps += mVessel.mArtComplex.mArtCannon.mAngGroupedFire * mVessel.mArtComplex.mArtCannon.mAngGroupedFire;

 arrMtrxShellDisp[0] =  valDispBet + mVessel.mArtComplex.mSigU * mVessel.mArtComplex.mSigU; //путевой угол
 arrMtrxShellDisp[ 2 * LEN_ARR_SCATTERS + 2] =  valDispEps + mVessel.mArtComplex.mSigU * mVessel.mArtComplex.mSigU; // угол наклона траетктории
 arrMtrxShellDisp[ 3 * LEN_ARR_SCATTERS + 3] = mVessel.mShellBody.mDispV0; // путевая скорость
 arrMtrxShellDisp[ 4 * LEN_ARR_SCATTERS + 4] = mVessel.mShellBody.mDispMass0; // масса
 arrMtrxShellDisp[ 5 * LEN_ARR_SCATTERS + 5] = mVessel.mShellBody.mDispCx; // коэф Cx
 arrMtrxShellDisp[ 6 * LEN_ARR_SCATTERS + 6] = mVessel.mShellBody.mDispCz;//клэф Cz
	if (bCOrrection)
	{
		return;
	}



 arrMtrxShellDisp[ 7 * LEN_ARR_SCATTERS + 7] =  (mEnvironment.mWind_V * 0.1 ) * (mEnvironment.mWind_V * 0.1); // гориз ветер модуль
 arrMtrxShellDisp[ 8 * LEN_ARR_SCATTERS + 8] = (5. / 180. * M_PI) * (5. / 180. * M_PI);
 arrMtrxShellDisp[9 * LEN_ARR_SCATTERS + 9] =  (mEnvironment.mWind_VertV * 0.1 ) * (mEnvironment.mWind_VertV * 0.1); // вертик ветер модуль
 arrMtrxShellDisp[  LEN_ARR_SCATTERS +1 ] = 0.0001;// плотн атмосфкеры


}

//-----------------------------------------------------------------------------------------------------------------
// заданы углы наведения в КГСК
// надо вычислить ошибки отработки этих углов в ГСК,
// вызванные погрешностями СИНС и деформациями корпуса корабля
// INPUT:
// VAlKGSKEps,  VAlKGSKBet  - углы наведения
// VAlDispDEltaPsi,  VAlDispDEltaTet,  VAlDispDEltaQ - дисперсии палогрешностей определения палубных углов
// VAlPsi,  VAlTet,  VAlQ  - палубные углы
// OUTPUT:
// arrK_DelBetDelEpsGSK - корреляционная матриуца ошибок отработки
// вектора углов наведения  в ГСК AlfT = || Bet  EPS||
//
//-----------------------------------------------------------------------------------------------------------------
// заданы углы наведения в КГСК
// надо вычислить ошибки отработки этих углов в ГСК,
// вызванные погрешностями СИНС и деформациями корпуса корабля
// INPUT:
// VAlKGSKEps,  VAlKGSKBet  - углы наведения
// VAlDispDEltaPsi,  VAlDispDEltaTet,  VAlDispDEltaQ - дисперсии палогрешностей определения палубных углов
// VAlPsi,  VAlTet,  VAlQ  - палубные углы
// OUTPUT:
// arrK_DelBetDelEpsGSK - корреляционная матриуца ошибок отработки
// вектора углов наведения  в ГСК AlfT = || Bet  EPS||
//
void TFight::calcCorMtrx_OshibokOtrabotki_EpsBet(const double VAlKGSKEps, const double VAlKGSKBet
	 ,const double VAlDispDEltaPsi, const double VAlDispDEltaTet, const double VAlDispDEltaQ
	 ,const double VAlPsi, const double VAlTet, const double VAlQ
	, double *arrK_DelBetDelEpsGSK)
{

	// 1. вычислекние углов прицеливания в ГСК
	double valGSKEps = 0., valGSKBet = 0.;
	transform_EpsBetKGSK_to_EpsBetGSK(&(mVessel.marrVectSost[3]), mVessel.mShellBody.mV0
	 ,  VAlKGSKEps, VAlKGSKBet, &valGSKEps, &valGSKBet);

	 double valKGSKEps0 = 0., valKGSKBet0 = 0., valk = 0.;
	 transform_EpsBetGSK_to_EpsBetKGSK(&(mVessel.marrVectSost[3]), mVessel.mShellBody.mV0
   ,valGSKEps,valGSKBet, &valKGSKEps0, &valKGSKBet0, &valk);
	///
   // 2. вектор направляющих косинусов в ГСК
   double arrx[3] ={0.};
   arrx[0]= cos(valGSKEps)* sin(valGSKBet);
   arrx[1]= cos(valGSKEps)* cos(valGSKBet) ;
   arrx[2]= sin(valGSKEps);
   ///

   // вычисление матрицы dF_po_dx И обратной
   double arr_dF_po_dx[9] = {0.}, arr_dF_po_dx_Inv[9] = {0.}, arr_dk_po_dx[3] = {0.}
   , arrT1[3] = {0.}, arrT2[9] = {0.}, arrT3[9] = {0.};
   double val_ro = (valk- 1./mVessel.mShellBody.mV0)/
   (valk - ScalProduct(arrx , &(mVessel.marrVectSost[3]),3)/mVessel.mShellBody.mV0)
   /mVessel.mShellBody.mV0;

   MatrxMultScalar(arrx, 1, 3, val_ro,arrT1);
   MtrxMultMatrxTransp(arrT1,3, 1, &(mVessel.marrVectSost[3]),1, arrT2) ;
   for (int i =0; i < 3; i++)
   {
	 arrT3[i*3 +i] =valk;
   }
   MtrxSumMatrx(arrT2, arrT3,3, 3, arr_dF_po_dx) ;
   InverseMtrx3(arr_dF_po_dx, arr_dF_po_dx_Inv);
   // 3.вычисление корреляционной матрицы рассеяния вектора направляющих косинусов в ГСК,
   // вызванного погрешностями СИНС по скорости
   double arrQ_KGSK[9] = {0.}, arrMtrxK_from_V_Sins[9] = {0.};
   double  arr_vvess[3] = {0.}, arr_kx[3] = {0.}
   , arrT4[3] = {0.}, arrT5[3] = {0.}, arrT6[9] = {0.}, arrT7[9] = {0.};
   MatrxMultScalar(&(mVessel.marrVectSost[3]), 1, 3, 1./mVessel.mShellBody.mV0,arr_vvess);

   MatrxMultScalar(arrx, 1, 3, valk,arr_kx);
   MtrxSumMatrx(arr_kx, arr_vvess,3, 3,  arrT4) ;

   double temp1 = 1./(valk - ScalProduct(arrx , arr_vvess,3)) ;
   MatrxMultScalar(arrx, 1, 3, temp1,arrT5);

   MtrxMultMatrxTransp(arrT5,3, 1, arrT4,3, arrT6) ;
   for (int i=0; i < 9; i++)
   {
	 arrT6[i] +=1.;
   }
   MtrxMultMatrx(arr_dF_po_dx_Inv, 3, 3, arrT6,3, arrT7);
   MatrxMultScalar(arrT7, 1, 9, -1./(mVessel.mShellBody.mV0 * valk),arrQ_KGSK);



		   // матрица дисперсий погрешностей ИС по компонентам скорости
   double arrK_Sins_V[9] = {0.};
   arrK_Sins_V[0] = mVessel.mSins.mSigV * sin(mVessel.mQ0) * mVessel.mSins.mSigV * sin(mVessel.mQ0);
   arrK_Sins_V[4] = mVessel.mSins.mSigV * cos(mVessel.mQ0) * mVessel.mSins.mSigV * cos(mVessel.mQ0);
   arrK_Sins_V[8] = mVessel.mSins.mSig_VH * mVessel.mSins.mSig_VH;

	MtrxMultMatrx_MultMatrxTransp(arrQ_KGSK, arrK_Sins_V, arrQ_KGSK,3, arrMtrxK_from_V_Sins);
	///

	// вычисление корреляционной матрицы рассеяния вектора направляющих косинусов в ГСК,
	// выхванного погрешэкстраполяции ностями по палцбным углам
	 // 4. вектор направляющих косинусов в КГСК
   double arry[3] ={0.};
   arry[0]= cos(VAlKGSKEps)* sin(VAlKGSKBet);
   arry[1]= cos(VAlKGSKEps)* cos(VAlKGSKBet) ;
   arry[2]= sin(VAlKGSKEps);
   ///

   // вычисление матрицы Q ПСК-КГСК
   double arrQ_PSK_KGSK[9] = {0.},arrQ_PSK_GSK[9] = {0.};
   arrQ_PSK_KGSK[0] = -arry[1];
   arrQ_PSK_KGSK[1] = -arry[2] * sin(VAlQ);
   arrQ_PSK_KGSK[2] = -arry[1]* sin(VAlPsi) - arry[2] * cos(VAlQ)* cos(VAlPsi);
   arrQ_PSK_KGSK[3] =  arry[0];
   arrQ_PSK_KGSK[4] = -arry[2] * cos(VAlQ);
   arrQ_PSK_KGSK[5] =  arry[0]* sin(VAlPsi) + arry[2] * sin(VAlQ)* cos(VAlPsi);
   arrQ_PSK_KGSK[6] =  0.;
   arrQ_PSK_KGSK[7] =  arry[0] * sin(VAlQ) +arry[1] * cos(VAlQ);
   arrQ_PSK_KGSK[8] =  arry[0] * cos(VAlQ)*cos(VAlPsi)- arry[1] * sin(VAlQ) *cos(VAlPsi);
	 ///
   // вычисление матрицы Q ПСК-ГСК
   MtrxMultMatrx(arr_dF_po_dx_Inv, 3, 3, arrQ_PSK_KGSK,3, arrQ_PSK_GSK);
   ///
	  // корреляционная матрица погрешностей СИНС по углам
	  // в векторе параметров углы следуют в последовательности Q, PSI, Tet
	  double arrK_Sins_Q_Psi_Tet[9] = {0.};
	  arrK_Sins_Q_Psi_Tet[0] = VAlDispDEltaQ;
	  arrK_Sins_Q_Psi_Tet[4] = VAlDispDEltaPsi;
	  arrK_Sins_Q_Psi_Tet[8] = VAlDispDEltaTet;

	  double arrMtrxK_fromDeckAngs[9] = {0.};
	  MtrxMultMatrx_MultMatrxTransp(arrQ_PSK_GSK, arrK_Sins_Q_Psi_Tet, arrQ_PSK_GSK,3, arrMtrxK_fromDeckAngs);
	  ///

	  // вуычисление суммарной колррел матрицы разбросов вектора направляющих
	  // косинусов в ГСК
	  double arrMtrxK_x[9] = {0.};
	  MtrxSumMatrx(arrMtrxK_fromDeckAngs, arrMtrxK_from_V_Sins,3, 3, arrMtrxK_x) ;
	  ///

	  double arrMtrx_dfi_po_dx[6] = {0.};
	  double temp = arrx[0] * arrx[0] + arrx[1] * arrx[1];
	  arrMtrx_dfi_po_dx[0] =  arrx[1]/ temp;
	  arrMtrx_dfi_po_dx[1] = -arrx[0]/ temp;
	  arrMtrx_dfi_po_dx[5] = 1./ sqrt(temp);

	  double arrTemp[6] = {0.};
	  MtrxMultMatrx(arrMtrx_dfi_po_dx,2, 3, arrMtrxK_x,3, arrTemp) ;
	  MtrxMultMatrxTransp(arrTemp,2, 3, arrMtrx_dfi_po_dx,2, arrK_DelBetDelEpsGSK) ;

}
//------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------

//
void transformVectNaprCos_To_EpsBet(double *arr_g, double* pvalEps0, double* pvalBet0)
{
 double valv = sqrt( arr_g[0]* arr_g[0] +  arr_g[1]* arr_g[1] );
 if ( valv < 0.000000001)
 {
   *pvalBet0 = 0;
	 *pvalEps0 = M_PI /2.;
 }

double temp =  arr_g[2];
if (fabs (temp) > 0.9999999999)
{
 temp =  0.9999999999 * SIGN_D(temp);
}
 *pvalEps0 = asin(temp);
 *pvalBet0 = atan2( arr_g[0], arr_g[1] );
}

//----------------------------------------------------------------------------------------------
// пересчет углов наведения из КГСК в ГСК
void TFight::transform_EpsBetKGSK_to_EpsBetGSK(double *arrVessVelocity_GSK, const double VAlShellV0
	 ,const double  VAlEpsKGSK0,const double    VAlKGSKBet0
			 , double *valGSKEps0,double *valGSKBet0)
{
  // вектор начальной скорости в ГСК
		double arrGSKV0[3] = {0.};
	  arrGSKV0[0] =  arrVessVelocity_GSK[0]   + VAlShellV0 * cos(VAlEpsKGSK0) * sin( VAlKGSKBet0);
		arrGSKV0[1] =  arrVessVelocity_GSK[1] + VAlShellV0 * cos(VAlEpsKGSK0) * cos( VAlKGSKBet0);
		arrGSKV0[2] =  arrVessVelocity_GSK[2] + VAlShellV0 * sin(VAlEpsKGSK0) ;

		double valV = Norm3( arrGSKV0) ;
		double temp =  arrGSKV0[2]/ valV;
			if (fabs (temp ) > 0.99999999999)
			{
			temp = 0.999999999* SIGN_D (temp);
			}

		*valGSKEps0 = asin( temp);
		*valGSKBet0 = atan2(arrGSKV0[0],arrGSKV0[1]);
}
// ---------------------------------------------------------------------------

// Вычисление результирующей вероятности поражения возд цели нп основе вычисленного массива вероятностей
// поражения в каждой точке встречи
// INPUT:
// arrProbab[IQuantShots] -  массив вероятностей
// IQuantShots - длина этого массива (к-во выстрелов)
//  OUTPUT:
// arrGladkProbab  - вероятности поражения с учетом кривой Гладковского нарастающим итогом
// *pvalProb  - вероятность поражения безучета кривой Гладковского
// *pvalGladkProb  - вероятность с учетом кривой галадковского
//
void TFight::calcAeroTargRezultProbability(double *arrProbab, double *arrDist, const int IQuantShots
	, double *pvalProb, double *pvalGladkProb)
{

// калибр 30 прямое попадание 2 раза  воздушные цели
	 if (mVessel.mShellBody.mEnumShellType == CALIBRO_30    )
	 {
			double valMult0 = 1.;
			for (int i = 0; i < IQuantShots; i++)
			{
				valMult0 *= (1. -arrProbab[i] * COEFF_BOEGOTOVNOSTY_SNARIADA);
			}

			double valSum = 0.;
			for(int j =0 ;j < IQuantShots; j++)
			{
			 double valMult1 = 1.;
				for (int i = 0; i < IQuantShots; i++)
				{
					 if ( i == j)
					 {
						 continue;
					 }
					 valMult1 *= (1. -arrProbab[i]* COEFF_BOEGOTOVNOSTY_SNARIADA);
				}
				valSum +=  arrProbab[j]* COEFF_BOEGOTOVNOSTY_SNARIADA * valMult1 ;

			}
			*pvalProb = 1. -  valMult0  - valSum;
			return;
	 }
///




// калибр 130, 100, 76 воздушные цели,  детонатор AR32A фугас  или контактный
if (
			 ((mVessel.mShellBody.mEnumShellType == CALIBRO_SHTAT_130    )
			 ||(mVessel.mShellBody.mEnumShellType == CALIBRO_TARAN_130    )
			 || (mVessel.mShellBody.mEnumShellType == CALIBRO_100    )
			 ||(mVessel.mShellBody.mEnumShellType == CALIBRO_76_SHTAT    )
			 ||(mVessel.mShellBody.mEnumShellType == CALIBRO_76_BARRIER   ) )
				&&
			 ((mVessel.mShellBody.mDetonator.mEnumDetonatorType == AR32A )
			 ||(mVessel.mShellBody.mDetonator.mEnumDetonatorType == CONTACT )
			 ||(mVessel.mShellBody.mDetonator.mEnumDetonatorType == DVM )
			 ||(mVessel.mShellBody.mDetonator.mEnumDetonatorType == AR51_LM )
			 ||(mVessel.mShellBody.mDetonator.mEnumDetonatorType == BARRIER ))
 )
 {



				 double  valPGladk = 0.;
				 double valMult  =1.;
				 for (int i = 0; i < IQuantShots; i++)
				 {
					 double valWeightGladk = mTarget.mplnGlagkovsky.LinearValueApprox(arrDist[i]);
					 valPGladk += valMult * arrProbab[i] * COEFF_BOEGOTOVNOSTY_SNARIADA * valWeightGladk;
					 valMult *= (1. - arrProbab[i] * COEFF_BOEGOTOVNOSTY_SNARIADA );
				 }
				 *pvalProb = 1. - valMult;
				 *pvalGladkProb = valPGladk;
				 return;

 }

		 ///
   if (
			( (mVessel.mShellBody.mEnumShellType == CALIBRO_SHTAT_130)
			||(mVessel.mShellBody.mEnumShellType == CALIBRO_TARAN_130    ))
			 &&
			 (mVessel.mShellBody.mDetonator.mEnumDetonatorType == MFIVU )
			 )

			 {



				 double  valPGladk = 0.;
				 double valMult  =1.;
				 for (int i = 0; i < IQuantShots; i++)
				 {
					 double valWeightGladk = mTarget.mplnGlagkovsky.LinearValueApprox(arrDist[i]);
					 valPGladk += valMult * arrProbab[i] * COEFF_BOEGOTOVNOSTY_SNARIADA * valWeightGladk;
					 valMult *= (1. - arrProbab[i] * COEFF_BOEGOTOVNOSTY_SNARIADA );
				 }
				 *pvalProb = 1. - valMult;
				 *pvalGladkProb = valPGladk;
				 return;

				}

}


// нахождение момента времени вхождения в зону сопровождения
double TFight::calcT_BeginSopr(const double VAlDistBeginSopr)
{
	   // нахождение точки дальнего рубежа в КГСК

		double arrVOtn[3] = {0.};

		arrVOtn[0]= mTarget.mTraject.marrVectSostGSK[3] - mVessel.mVVess * sin(mVessel.mQ0) ;
		arrVOtn[1]= mTarget.mTraject.marrVectSostGSK[4] - mVessel.mVVess * cos(mVessel.mQ0) ;
		arrVOtn[2]= mTarget.mTraject.marrVectSostGSK[5];
		double vala = Norm3 (arrVOtn) * Norm3 (arrVOtn);
		double valb = 2. * ScalProduct(mTarget.mTraject.marrVectSostGSK , arrVOtn, 3) ;
		double valc =   Norm3 (mTarget.mTraject.marrVectSostGSK) * Norm3 (mTarget.mTraject.marrVectSostGSK) - VAlDistBeginSopr * VAlDistBeginSopr;
		TComp cmpx1, cmpx2;
		int irez00 =  SolvEq2(vala,valb,valc, cmpx1, cmpx2);
		if (irez00 >2)
		{
			 return -1;
		}
		double valTBegin = MIN__(cmpx1.m_Re, cmpx2.m_Re);

		return  valTBegin;
}



// ---------------------------------------------------------------------------
// OUTPUT:
// *pvalProb   - вероятность поражения очередью
// *valProb0  - вероятность успеха одного испытания
// *pvalSKZPromach  - СКЗ промаха
// arrVectMiss_GSK[6]  - средний вектор промаха в точке встречи
// arrCoMtrx_1_and_3_Groups [36] - коррел. матрица разброса вектора промаха в точке встречи
bool TFight::calcSuccessProbSeaTarg(const int QuantShells,const double VAlHAntenna
	   ,double *pvalSKZPromach, double *arrCoMtrx_1_and_3_Groups, double *arrCoMtrx_2_Group
		,double *pvalKGSKEps, double *pvalKGSKBet, double *pvalTFlight
		,double *arrVectMiss_GSK, double *pvalProb,  double *pvalProb0)
{

			// 1. Решение задачи о точке встречи
  // вычисление вектора положения АУ в КГСК
	double arrPositionAY_KGSK[3] = {0.};
	mVessel.calcAY_Position(arrPositionAY_KGSK);
	double valKGSKEps = 0., valKGSKBet  = 0.;
	double arrVectAppointmentPointGSK[6] = {0.};
	double valMiss = -1., valTFlight = -1.;

	// вычисление истинного вектора состояния цели в КГСК на момомент выстрела
	double arrTargVS_KGSK0[6] ={0.};
	MtrxMinusMatrx(mTarget.mTraject.marrVectSostGSK_Begin, mVessel.marrVectSost,1, 6, arrTargVS_KGSK0);
	///
	double arrShellVeloAppointmPoint_GSK [3] ={0.},  arrTargVeloAppointmPoint_GSK[3] ={0.};
	if (!calcAppointmentPoint(NULL, NULL, &(mVessel.marrVectSost[3])
	,arrTargVS_KGSK0, arrPositionAY_KGSK
	, &valKGSKEps, &valKGSKBet, &valTFlight, arrVectAppointmentPointGSK, &valMiss
	, arrShellVeloAppointmPoint_GSK,  arrTargVeloAppointmPoint_GSK))
	{
		ShowMessage(L"ERROR_ TFight::calcSuccessProbSeaTarg");
		return false;
	}

 /*	if (!findAppointmentPoint_NewtonMeth(&(mVessel.marrVectSost[3])
	,arrTargVS_KGSK0, arrPositionAY_KGSK
	, &valKGSKEps, &valKGSKBet, &valTFlight, arrVectAppointmentPointGSK, &valMiss
	, arrShellVeloAppointmPoint_GSK,  arrTargVeloAppointmPoint_GSK))
	{
		ShowMessage(L"ERROR_ TFight::calcSuccessProbSeaTarg");
		return false;
	} */


 /////////////////////////////////////////////////////////////////////////
	//////////////  ЭТО ДЛЯ ОТЛАДКИ!!! СВЕРХУ РАСКОММЕНТИРОВАТЬ А ЭТО УДАЛИТЬ !!!!///////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////

//valKGSKEps =	0.35873669492725 ;
//valKGSKBet =	6.27370551180827 ;
//valTFlight =	43.3708000003162 ;
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////



	double valDispMiss = -1.;


	//вычисление вероятности прямого попадания

		double arrCorrMatrxVectMiss_GSK[36] = {	0.}, arrShellVS_GSK [6] = {0.};
		calcCorMtrx_First_And_Third_Group( true, valKGSKEps, valKGSKBet, valTFlight,arrCoMtrx_1_and_3_Groups, arrShellVS_GSK);
		const double VAlObservTime = 10.;
		calcCorMtrx_Second_Group(VAlObservTime, valTFlight , arrCoMtrx_2_Group);
		MtrxSumMatrx(arrCoMtrx_1_and_3_Groups, arrCoMtrx_2_Group, 6, 6, arrCorrMatrxVectMiss_GSK);


		memset(arrVectMiss_GSK, 0, 6 * sizeof(double));
		MtrxMinusMatrx(&arrShellVS_GSK[3], &(mTarget.mTraject.marrVectSostGSK_Begin[3]), 3, 1, &arrVectMiss_GSK[3]);

		const double VAlDistAppPoint = Norm3(mTarget.mTraject.marrVectSostGSK_Begin);
		///

		// вычисление вероятности
		TNeighbourhoodAppPoint NeighbourhoodAppPoint(arrVectMiss_GSK, arrCorrMatrxVectMiss_GSK , mTarget
		, mVessel.mShellBody, VAlDistAppPoint, mVessel.mArtComplex.mArtCannon.mSigRabT);
	//	NeighbourhoodAppPoint.calcProbDirectHitting_SeaTarg(pvalProb0, &valDispMiss);   // тут была ошибка в этой функции 15 мая 2018

		// OUTPUT:
// arrCorMtrxCartinSK[4] - коррел матрица вектора промаха в картинной плоскости
// *pPlgProjection - составной полигон, являющийся проекций полигонов, описывающих геометрию цели, на картинную плоскость
// *pvalDispMiss - СКЗ промаха
// Возвращает вероятнлость прямого попадания
double arrCorMtrxCartinSK[4] = {0.};
TURPolygon plgProjection ;
double valDispMiss1;
*pvalProb0 =  NeighbourhoodAppPoint.__calcProbDirectHitting__( &valDispMiss
, arrCorMtrxCartinSK, &plgProjection);
 ///



	 *pvalSKZPromach = sqrt( valDispMiss);
	 if ((*pvalProb0) < 0.0000000001)
	 {
		 *pvalProb = 0.  ;
		 return true;
	 }

	 if ((*pvalProb0) > 0.99999999999)
	 {
		 *pvalProb = 1.  ;
		 return true;
	 }
//double TNeighbourhoodAppPoint::calcRezultProbForSeaTarg( const int QUantShells
//, const double VAlProb0)

	 *pvalProb = NeighbourhoodAppPoint.calcRezultProbForSeaTarg( QuantShells
	 , *pvalProb0);

	return true;
}
// ---------------------------------------------------------------------------



// ВЫчисление вероятности поражения берег цели
//INPUT:
// QUantShells  - к-во выстрелов
// bCOrrection  - признак пристрелки
// QAantIspit - к-во испытаний метода монте-карло
// VAlKillingRange - радиус приведенной площади
// OUTPUT:
// *pvalSKZPromach  - скз промаха (не нужно, это атттавизм)
// arrScatterCorMtrx_1_and_3_Groups[4] - коррел матрица рзабросов точки падения, 1 и 3 группы
// arrScatterCorMtrx_2_Group[4] - - коррел матрица рзабросов точки падения, 2 группа
// *pvalProb - вероятность теоретическая
// EnumTypeOfControlAlgorithm  - тип алгоритма управления  (Fight.h)
// ppntArrAimingPoints - массив точек прицеливания
//*piQuantAimingPoints -длтна массива точек прицеливания
// piarrRepeatAimingPoints - массив повторений точек прицеливания (целочисленный)
// *pvalProb_Monte_Carlo - вероятность поражения, вычисленная методом монте-карло
bool TFight::calcSuccessProbCoast(const int QUantShells,const  bool bCOrrection
		,const int QAantIspit, const double VAlKillingRange, double *pvalSKZPromach
	, double *arrScatterCorMtrx_1_and_3_Groups, double *arrScatterCorMtrx_2_Group, double *arrShellVeloAppointmPoint_GSK, double *pvalProb
		,enumTypeOfControlAlgorithm EnumTypeOfControlAlgorithm, TURPointXY *ppntArrAimingPoints, int *piarrRepeatAimingPoints
		, int* piQuantAimingPoints, double *pvalProb_Monte_Carlo)

{
		memset(piarrRepeatAimingPoints, 0, QUantShells * sizeof(int));
		memset(ppntArrAimingPoints , 0, QUantShells * sizeof(TURPointXY));

	// 1. Решение задачи о точке встречи
	// вычисление вектора положения АУ в КГСК
	double arrPositionAY_KGSK[3] = {0.};
	mVessel.calcAY_Position(arrPositionAY_KGSK);
	double valKGSKEps = 0., valKGSKBet  = 0.;
	double arrVectAppointmentPointGSK[6] = {0.}, arrTargVeloAppointmPoint_GSK[3] ={0.};
	double valMiss = -1., valTFlight = -1.;

	// вычисление истинного вект ора состояния цели в КГСК на момомент выстрела
	double arrTargVS_KGSK0[6] ={0.};
	MtrxMinusMatrx(mTarget.mTraject.marrVectSostGSK_Begin, mVessel.marrVectSost,1, 6, arrTargVS_KGSK0);
	///
	if (!calcAppointmentPoint(NULL, NULL, &(mVessel.marrVectSost[3])
	,arrTargVS_KGSK0, arrPositionAY_KGSK
	, &valKGSKEps, &valKGSKBet, &valTFlight, arrVectAppointmentPointGSK, &valMiss
	, arrShellVeloAppointmPoint_GSK, arrTargVeloAppointmPoint_GSK))
	{
		ShowMessage(L"ERROR_ TFight::calcSuccessProbCoast");
		return false;
	}

	///

	// 2. Вычисление эллипсоида рассеяния
	 // 2.1 формирование диагональной матрицы LEN_ARR_SCATTERS x LEN_ARR_SCATTERS дисперсий возмущений - 8 тначаольных условий, коэфф формы, масса, ветер по X, ветер по Y
	double arrMtrxShellDisp [LEN_ARR_SCATTERS * LEN_ARR_SCATTERS] = {0.};
	fillShellVozmDispMatr (bCOrrection, valKGSKEps,  valKGSKBet, arrMtrxShellDisp);
		///

	// 2.2 формирование матрицы частных производных фазового вектора снаряда в точке встречи по начальным условиям Пси, Пи, Тетта, Масса, Cx? Cz, 3 параметра ветра - 9 штук
		 // массе, коэффициенту формы, горизонтальному ветру по осям X и Z в ССК , коэффициенту формы по оси Z- всего 13 штук
	double arrStrSK_Jacobian[8 * LEN_ARR_SCATTERS] = {0.};//, arrShellScatteringsCorMtrx_GSK[36] = {0.};
	double arrShellVS_GSK[6] = {0.}; // ВС снаряда в ГСК
	TMyShellTraj ShellTraj (&(mVessel.marrVectSost[3]), mVessel.mShellBody, valKGSKEps,  valKGSKBet );

	const double VAlStepInt = 0.0005;
	double arrShellScatteringsCorMtrxPos_SSK[9] = {0.};
	//вычисление вероятности прямого попадания

		double arrCorrMatrxVectMiss_GSK[36] = {	0.}, arrShellVS_GSKTemp [6] = {0.}
			 ,arrCoMtrx_1_and_3_Groups[36] = {0.},  arrCoMtrx_2_Group[36] = {0.};
		calcCorMtrx_First_And_Third_Group( bCOrrection, valKGSKEps, valKGSKBet, valTFlight,arrCoMtrx_1_and_3_Groups, arrShellVS_GSKTemp);
		const double VAlObservTime = 10.;
		calcCorMtrx_Second_Group(VAlObservTime, valTFlight , arrCoMtrx_2_Group);

		arrScatterCorMtrx_1_and_3_Groups[0] =  arrCoMtrx_1_and_3_Groups[0];
		arrScatterCorMtrx_1_and_3_Groups[1] =  arrCoMtrx_1_and_3_Groups[1];
		arrScatterCorMtrx_1_and_3_Groups[2] =  arrScatterCorMtrx_1_and_3_Groups[1];
		arrScatterCorMtrx_1_and_3_Groups[3] =  arrCoMtrx_1_and_3_Groups[7];

		arrScatterCorMtrx_2_Group[0] =  arrCoMtrx_2_Group[0];
		arrScatterCorMtrx_2_Group[1] =  arrCoMtrx_2_Group[1];
		arrScatterCorMtrx_2_Group[2] =  arrScatterCorMtrx_2_Group[1];
		arrScatterCorMtrx_2_Group[3] =  arrCoMtrx_2_Group[7];

		//	Формирование матрицы эллипса рассеяния точки падения в плоскости OXY ГСК
		 double arrElK[4] = {0.};
		MtrxSumMatrx(arrScatterCorMtrx_1_and_3_Groups, arrScatterCorMtrx_2_Group, 2, 2, arrElK);

	///
	// формирование корреляционной матрицы разбросов вект ора промаха
	double arrVectMiss[6] = {0.};
	memcpy( &arrVectMiss[3], &arrVectAppointmentPointGSK [3], 3 * sizeof(double));

	double valDist =  Norm3( mTarget.mTraject.marrVectSostGSK_Begin) ;

		   double valTarfCellSize = 15.;
		   if(
			 (mTarget.menumTargetType == COVERED_MANPOWER_ENTRENCH)
		   ||(mTarget.menumTargetType == COVERED_MANPOWER_TRENCH)
		   ||(mTarget.menumTargetType == MANPOWER_ARMOURED_CARRIER)
		   ||(mTarget.menumTargetType == MANPOWER_CAR)
		   ||(mTarget.menumTargetType == PLATOON_POINT)
			 )
		   {
			 valTarfCellSize = 0.25;
		   }
			const double VAlTargCellSize = valTarfCellSize;
			const double VAlAimCellSize  = 15.;
			// вычисление угла поворота (ориентации) цели в ГСК
			// считается от оси X положительное напргоавление против часовой стрелки
			double valRotateAng = atan2(mTarget.mTraject.marrVectSostGSK_Begin[4], mTarget.mTraject.marrVectSostGSK_Begin[3]);
			///
			// поворот полигона цели на угол  valRotateAng
			const TURPointXY pntSdvig(0.,0.);
			const double valRastigenie  = 1.;
			TURPolygon  plgTargGSK = mTarget.mpArrPlanePolygon[0].mPolygon.LinTransform(valRotateAng,pntSdvig, valRastigenie ) ;
			TURFigure *pFigure ;
			TURPolygon *pPolygon =  &plgTargGSK;
			pFigure = (TURFigure*)(pPolygon);
				double valObj = -1;
	 if (EnumTypeOfControlAlgorithm == OPTIMAL)
	 { // НАДО БУДЕТ ПОТОМ ПЕРЕДЕЛАТЬ !!!!!!

			TCoastTargNeighbourhood  CoastTargNeighbourhood ( pFigure ,
			VAlTargCellSize ,  VAlAimCellSize , arrElK,
			VAlKillingRange, QUantShells ) ;



			CoastTargNeighbourhood.calcOptimalArray_Of_AimPoints( ppntArrAimingPoints
			,piQuantAimingPoints, piarrRepeatAimingPoints, &valObj) ;
			*pvalProb = ( ((double)(CoastTargNeighbourhood.mMultPntTarg.NumPoints)) -  valObj)/ ((double)(CoastTargNeighbourhood.mMultPntTarg.NumPoints));

	 }

	 if (EnumTypeOfControlAlgorithm == STANDART)
	 {
		 TCoastTargNeighbourhood::calcAimingPoints_For_OpenManPower_MUS(arrElK, VAlKillingRange, QUantShells
				,mTarget.mpArrPlanePolygon[0].mPolygon,valDist,  ppntArrAimingPoints, piarrRepeatAimingPoints
				,piQuantAimingPoints);

		for (int i = 0; i < (*piQuantAimingPoints); i++)
		{
			ppntArrAimingPoints[i] = ppntArrAimingPoints[i].LinTransform(valRotateAng,pntSdvig, valRastigenie ) ;
		}
		TCoastTargNeighbourhood CoastTargNeighbourhood(pFigure , VAlTargCellSize ,VAlAimCellSize, ppntArrAimingPoints
		 , *piQuantAimingPoints,arrScatterCorMtrx_1_and_3_Groups
			, arrScatterCorMtrx_2_Group,VAlKillingRange, QUantShells );
			double *pX = new double  [*piQuantAimingPoints];
			for (int i = 0; i < (*piQuantAimingPoints); i++)
			{
			 pX[i] = (double)piarrRepeatAimingPoints[i];
			}
		 //	valObj = CoastTargNeighbourhood.calcEfficiencyOfStrategy( pX);
		  	valObj = CoastTargNeighbourhood.calcEfficiencyOfStrategy_With_SystMtrx_Var1( pX);
		 // valObj = 0.;



			delete []pX;

			*pvalProb = ( ((double)(CoastTargNeighbourhood.mMultPntTarg.NumPoints )) -  valObj)
		  	/ ((double)(CoastTargNeighbourhood.mMultPntTarg.NumPoints));
			*pvalProb_Monte_Carlo = CoastTargNeighbourhood.estimateStrategy_MonteCarlo(QAantIspit, piarrRepeatAimingPoints) ;
	 }

return true;
}


// ---------------------------------------------------------------------------
// INPUT:
// bCOrrection -  признак коррекции (пристрелки)
// VAlKGSKEps, VAlKGSKBet -  углы наведения
// VAlTFlight - полеьтное время
// OUTPUT:
// arrCoMtrxGSK_1_and_3_Groups[36] - кррел матрица разбросов
// arrShellVS_GSK [6] - вектор состояния снар в ГСК

void TFight::calcCorMtrx_First_And_Third_Group(const bool bCOrrection ,const double VAlKGSKEps
	,  const double VAlKGSKBet, const double VAlTFlight
	, double* arrCorMtrxGSK_1_and_3_Groups, double *arrShellVS_GSK)
{

	 // 1 формирование диагональной матрицы LEN_ARR_SCATTERS x LEN_ARR_SCATTERS дисперсий возмущений - 8 тначаольных условий, коэфф формы, масса, ветер по X, ветер по Y
	double arrMtrxShellDisp [LEN_ARR_SCATTERS * LEN_ARR_SCATTERS] = {0.};

	fillShellVozmDispMatr (bCOrrection, VAlKGSKEps,  VAlKGSKBet, arrMtrxShellDisp);
		///

	// 2 формирование матрицы частных производных фазового вектора снаряда в точке встречи по начальным условиям Пси, Пи, Тетта, Масса, Cx? Cz, 3 параметра ветра - 9 штук
		 // массе, коэффициенту формы, горизонтальному ветру по осям X и Z в ССК , коэффициенту формы по оси Z- всего 13 штук
	double arrStrSK_Jacobian[8 * LEN_ARR_SCATTERS] = {0.};
	TMyShellTraj ShellTraj (&(mVessel.marrVectSost[3]), mVessel.mShellBody, VAlKGSKEps,  VAlKGSKBet );

	const double VAlStepInt = 0.0005;
	double arrShellScatteringsCorMtrxPos_SSK[9] = {0.};

// вычисление корреляционной матрицы разброса вектора состояния снаряда
// в ГСК
// INPUT:
// Environment -внешняя среда
// VAlStepInt - шаг интегрирования
// VAlTFlight - полетное время
// arrMtrxShellDisp  - диагональная матрица диспервсий разбросов параметров
//OUTPUT:
// arrStrSK_Jacobian [LEN_ARR_SCATTERS * 8] - матрица частных производных фазового вектора
// снаряда в траектоореной сиситеме координат по параметрам
// arrShellScatteringsCorMtarx_GSK [ 6 * 6] - корреляционная матрица в ГСК
// arrShellScatteringsCorMtrxPos_SSK [3*3] - корреляцилнная матрица ошибак рассеяния положения в ССК
// arrShellVS_GSK[6] - вектор состояния снаряда в ГСК
	ShellTraj.calc_VS_GSK_And_ScatteringsCorrMatrx_GSK (mEnvironment, VAlStepInt
	,VAlTFlight,arrMtrxShellDisp, arrStrSK_Jacobian, arrCorMtrxGSK_1_and_3_Groups, arrShellVS_GSK, arrShellScatteringsCorMtrxPos_SSK);
}

//---------------------------------------------------------------------------------

// arrCorMtrxGSK_2_Group[36]
void TFight::calcCorMtrx_Second_Group(const double VAlObservTime, const double VAlExtrapTime,  double* arrCorMtrxGSK_2_Group)
{
		TFight FightCur = *this;
		if(
		(mTarget.menumTargetType == OPEN_MANPOWER_LIE)
		 ||(mTarget.menumTargetType == OPEN_MANPOWER_STAND)
		 ||(mTarget.menumTargetType == BULLET_PROOF_LIE)
		 ||(mTarget.menumTargetType == BULLET_PROOF_STAND)
		 ||(mTarget.menumTargetType == COVERED_MANPOWER_ENTRENCH)
		 ||(mTarget.menumTargetType == COVERED_MANPOWER_TRENCH)
		 ||(mTarget.menumTargetType ==  MANPOWER_ARMOURED_CARRIER)
		 ||(mTarget.menumTargetType ==  MANPOWER_CAR)
		 ||(mTarget.menumTargetType ==  PLATOON_POINT)
		 ||(mTarget.menumTargetType ==  GROUP_POINT_COAST)
		 )
		 {
			FightCur.mTarget.mTargEPR = 1000.;
			FightCur.mVessel.mTransmitAnt.mPowerPrd = 10000000.;
			FightCur.mVessel.mTransmitAnt.mKYPrd = 10000.;
			memset(FightCur.mTarget.mTraject.marrSigW, 0., 3 * sizeof(double));
			memset(&(FightCur.mTarget.mTraject.marrVectSostGSK_Begin[3]), 0., 3 * sizeof(double));
			memset(&(FightCur.mTarget.mTraject.marrVectSostGSK[3]), 0., 3 * sizeof(double));
		 }
	FightCur.mVessel.mMaxQ = 0.; // корабль двигается по идеальной прямой

	memset( FightCur.mTarget.mTraject.marrSigW, 0, 3 * sizeof(double)); // цель двигается по идеальной прямой





	// передвижение члена класса FightCur на время  valTCurrentFight
	// поскольку передвижение совершается дискретными шагами кратными
	// темпу фильтрации, то текущее время FightCur вообще гворя может быть меньше  valTCurrentFight
	FightCur.shift(VAlObservTime);
	///



	// вычисление корреляционной матрицы разбросов вектора состояния цели на момент  VAlTCurrentShot + VAlTFlight
	double  arrTargExtrpVS_KGSK[6] = {0.};
	 FightCur.extrapolateTrueTargVS_And_ScatteringCorrMatrx_KGSK(VAlObservTime + VAlExtrapTime
		 , arrCorMtrxGSK_2_Group, arrTargExtrpVS_KGSK);

}

//---------------------------------------------------------------------------------
// приближенное решене РЗВ
//INPUT:
// arrVSTarg_KGSK0  - вектор состояния цели
// val_Shell_V0 - модуль скорости снаряда
// OUTPUT:
// Eps0, Bet0 - углы
//
void fncInit_Eps0_Bet0(double *arrVSTarg_KGSK0, double val_Shell_V0, double& Eps0, double& Bet0)
 {
	double Xc, Yc, Hc, Vcx, Vcy, Vch;
	double A, B, C, D;
	double S1, S2, T0;
   //	double Vsnx, Vsny, Vsnh, Vsn;

	// Вычисление начальных значений углов наведения Fi0, Qu0
	Xc =   arrVSTarg_KGSK0[0];
	Yc =   arrVSTarg_KGSK0[1];
	Hc =   arrVSTarg_KGSK0[2];
	Vcx =  arrVSTarg_KGSK0[3];
	Vcy =  arrVSTarg_KGSK0[4];
	Vch =  arrVSTarg_KGSK0[5];

	A = Xc * Xc + Yc * Yc + Hc * Hc;
	B = 2 * (Xc * Vcx + Yc * Vcy + Hc * Vch);
	C = Vcx * Vcx + Vcy * Vcy + Vch * Vch - val_Shell_V0 * val_Shell_V0;
 //	A /= 10000.0;
 //	B /= 10000.0;
 //	C /= 10000.0;
 A /= A;
 B /= A;
 C /= A;

 //	double temp = B * B - 4 * A * C;
 //	if (temp < 0.00001))
 //	{
 //	 temp = 0.00001;
 //	}

	D = sqrt(B * B - 4 * A * C);

   //	S1 = (-B - D) / (2.0 * A);
	S2 = (-B + D) / (2.0 * A);

	T0 = 1.0 / S2;

   //	Vsnx = Xc * S2 + Vcx;
  //	Vsny = Yc * S2 + Vcy;
  //	Vsnh = Hc * S2 + Vch;
   //	Vsn = sqrt(Vsnx * Vsnx + Vsny * Vsny + Vsnh * Vsnh);

	// Начальные значения углов бросания Fi0 Qu0 для решения задачи встречи
	Xc = Xc + Vcx * T0;
	Yc = Yc + Vcy * T0;
	Hc = Hc + Vch * T0;
	D = sqrt(Xc * Xc + Yc * Yc + Hc * Hc);

			double temp =  Hc / D;
			if (fabs (temp ) > 0.99999999999)
			{
			temp = 0.999999999* SIGN_D (temp);
			}

	Eps0 = asin(temp);
	D = sqrt(Xc * Xc + Yc * Yc);

	temp =  Yc / D;
	if (fabs (temp ) > 0.99999999999)
	{
	temp = 0.999999999* SIGN_D (temp);
	}

	if (Xc >= 0)
	{

		Bet0 = acos(temp);
	}
	else
	{
		Bet0 = 2.0 * M_PI - acos(temp);
	}
}

//---------------------------------------------------
//-------------------------------------------------------------------------

// формирование 10х10 диагональной матрицы возмущений для снаряда
// при стрельбе с земли

void TFight::fillShellVozmDispMatr_ShootingEarth (double *arrMtrxShellDisp)
{



// 0. marrStrSK_VS [3]-  угол Пси

// 1. marrStrSK_VS [5]-  относит плотность атмосферы Пи
// 2.  marrStrSK_VS [7]-  угол Тетта
// 3. marrStrSK_VS [6]-  путевая скорость
// дополнительные разбросы:
// 4.  масса
// 5.  коэфф формы Cx
// 6. коэф формы по оси Z
// 7. модуль горизонтальной скорости ветра
// 8. направление горизонтального ветра
// 9. модуль вертикального ветра
 memset(arrMtrxShellDisp, 0,  sizeof(double) * LEN_ARR_SCATTERS * LEN_ARR_SCATTERS);





 arrMtrxShellDisp[0] =  0.0007 * 0.0007 ; //путевой угол      это для 76 кал
 arrMtrxShellDisp[ 2 * LEN_ARR_SCATTERS + 2] =  0.0007 * 0.0007 ; // угол наклона траетктории  это для 76 кал
 arrMtrxShellDisp[ 3 * LEN_ARR_SCATTERS + 3] = mVessel.mShellBody.mDispV0; // путевая скорость
 arrMtrxShellDisp[ 4 * LEN_ARR_SCATTERS + 4] = mVessel.mShellBody.mDispMass0; // масса
 arrMtrxShellDisp[ 5 * LEN_ARR_SCATTERS + 5] = mVessel.mShellBody.mDispCx; // коэф Cx
 arrMtrxShellDisp[ 6 * LEN_ARR_SCATTERS + 6] = mVessel.mShellBody.mDispCz;//клэф Cz




 //arrMtrxShellDisp[ 7 * LEN_ARR_SCATTERS + 7] =  (mEnvironment.mWind_V * 0.1 ) * (mEnvironment.mWind_V * 0.1); // гориз ветер модуль
 //arrMtrxShellDisp[ 8 * LEN_ARR_SCATTERS + 8] = (5. / 180. * M_PI) * (5. / 180. * M_PI);
 //arrMtrxShellDisp[9 * LEN_ARR_SCATTERS + 9] =  (mEnvironment.mWind_VertV * 0.1 ) * (mEnvironment.mWind_VertV * 0.1); // вертик ветер модуль
 arrMtrxShellDisp[  LEN_ARR_SCATTERS +1 ] = 0.0001;// плотн атмосфкеры


}



//------------------------------------------------------------------------------------------------------------------
// пересчет углов наведения из ГСК в КГСК
void TFight::transform_EpsBetGSK_to_EpsBetKGSK(double *arrVessVelocity_GSK, const double VAlShellV0
   ,const double  VAlEpsGSK0,const double    VAlGSKBet0
			 , double *valEps0,double *valBet0)
{
  // вектор направляющих косинусов скорости в КГСК (искомый)
  double arr_g[3] = {0.};
  ///
  // модуль скорости корабля
  double valVessV0 = Norm3( arrVessVelocity_GSK) ;
  ///

  double arr_f[3] ={0.};// вектор направляющих косинусов
  arr_f[0] = cos(VAlEpsGSK0)* sin (VAlGSKBet0);
  arr_f[1] = cos(VAlEpsGSK0)* cos (VAlGSKBet0);
  arr_f[2] = sin(VAlEpsGSK0);
  ///
 // квадратное уравнение относительно k
 double a = 1. ;
		  double b = -2. * ScalProduct(arr_f , arrVessVelocity_GSK, 3) ;
		  double c = valVessV0 * valVessV0 - VAlShellV0 * VAlShellV0;
		  TComp x1,x2 ;
			// Решение квадраьного улавнения a*x*x + b*x + c =0
			// возвращает
			// 0 - 2 действительных некраьных корня
			// 1 - 2 действительных кратных корня
			// 2- имеется по крайней мере один нулевой корень a!= 0, c=0
			// 3 - 2 комплексно сопряженных корня
			// 4 -  1 действительный корень (a =0)
			// 5  - несовместность a=b=0, c!=0
			// 6 - тождество )a=b=c=0
		  int irez =  SolvEq2( a, b, c, x1, x2);

		  double arrx[2] = {0.};
		  arrx[0] = x1.m_Re;
		  arrx[1] = x2.m_Re;
		  for (int i =0; i < 2; i++)
		  {
		  double arrT0[3] = {0.}, arrT1[3] = {0.}, arrT2[3] = {0.};
		  MatrxMultScalar(arr_f, 1, 3, arrx[i],arrT0);
		  MtrxMinusMatrx(arrT0, arrVessVelocity_GSK,1, 3, arrT1);
		  MatrxMultScalar(arrT1, 1, 3, 1./VAlShellV0, arr_g);
		  if (arr_g[2] >= 0.)
		  {
			break;
		  }
		  }
	   double t = sqrt(arr_g[0] * arr_g[0] + arr_g[1] * arr_g[1]);
	   if (fabs(t)> 1.)
	   {
		 t = 0.999999 * SIGN_D(t);
	   }
	  *valEps0 =  acos(t );

	  t =  arr_g[0]/ t;
	  if (fabs(t)> 1.)
	   {
		 t = 0.999999 * SIGN_D(t);
	   }

	  //*valBet0 =  asin(t );
	  *valBet0 =  atan2(arr_g[0], arr_g[1]);


	return;

  }

//------------------------------------------------------------------------------------------------------------------
// пересчет углов наведения из ГСК в КГСК    ПЕРЕГРУЖЕННАЯ!
void TFight::transform_EpsBetGSK_to_EpsBetKGSK(double *arrVessVelocity_GSK
, const double VAlShellV0 ,const double  VAlEpsGSK0,const double    VAlGSKBet0
			 , double *valEps0,double *valBet0,double *valk)
{
  // вектор направляющих косинусов скорости в КГСК (искомый)
  double arr_g[3] = {0.};
  ///
  // модуль скорости корабля
  double valVessV0 = Norm3( arrVessVelocity_GSK) ;
  ///

  double arr_f[3] ={0.};// вектор направляющих косинусов
  arr_f[0] = cos(VAlEpsGSK0)* sin (VAlGSKBet0);
  arr_f[1] = cos(VAlEpsGSK0)* cos (VAlGSKBet0);
  arr_f[2] = sin(VAlEpsGSK0);
  ///
 // квадратное уравнение относительно k
		double a = 1. ;
		double b = -2. * ScalProduct(arr_f , arrVessVelocity_GSK, 3)/VAlShellV0 ;
		double c = valVessV0 * valVessV0 / (VAlShellV0 * VAlShellV0) -1.;
		*valk = -b/2. + sqrt(b/2.* b/2. -c);

		double arrT0[3] = {0.}, arrT1[3] = {0.}, arrT2[3] = {0.};
		  MatrxMultScalar(arr_f, 1, 3, *valk,arrT0);
		  MatrxMultScalar(arrVessVelocity_GSK, 1, 3, 1./VAlShellV0,arrT1);
		  MtrxMinusMatrx(arrT0, arrT1,1, 3, arr_g);


	   double t = sqrt(arr_g[0] * arr_g[0] + arr_g[1] * arr_g[1]);
	   if (fabs(t)> 1.)
	   {
		 t = 0.999999 * SIGN_D(t);
	   }
	  *valEps0 =  acos(t );

	  t =  arr_g[0]/ t;
	  if (fabs(t)> 1.)
	   {
		 t = 0.999999 * SIGN_D(t);
	   }

	  *valBet0 =  asin(t );


	return;

  }

void TFight::calcTargPosKGSK(const double valTExtr, double *arrTrueTargPosKGSK)
{
  TVessel VesselCur = mVessel;
  double arrVSVess[9] = {0.};
  VesselCur.VSProlong(valTExtr, arrVSVess);


  TTarget TargetCur = mTarget;
  TargetCur.recalcTrajPoint(valTExtr + mTarget.mTraject.mTCur);


  MtrxMinusMatrx(TargetCur.mTraject.marrVectSostGSK, arrVSVess,1, 3, arrTrueTargPosKGSK);
}

#pragma package(smart_init)
