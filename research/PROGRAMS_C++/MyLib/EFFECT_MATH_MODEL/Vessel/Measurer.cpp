//---------------------------------------------------------------------------


#pragma hdrstop

#include "Measurer.h"

#include <vcl.h>

#include "Traject.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
 #include "YrWriteShapeFile.h"

const int QUANT_COLS_BUFF_MEASURER = 7 ;

TMeasurer::TMeasurer()
{
	mSigV =  0.00021;
	mSigU = 0.00042;
	mSigR = 15 ;

	mRzv  = 0;; //  замер дальности
	mVzv  = 0;  // замер по углу V (курс)
	mUzv  = 0; // замер по углу U (угол места)
	mDelR = 0;// ошибка измерения дальности
	mDelV = 0; // ошибка ищзмерения угла V
	mDelU = 0;
	mT = 0. ;
	mbOEChannel = false ;
	mRadar = TRadar();
	mOEChannel = TOEChannel();
	mpwcharrFoldReport = NULL ;
	mparrBuff   = NULL ;
	mQuantPntReport = 0;
	mLenMemoryAlloc = 0 ;
	mAlf = 1 ;

}

//---------------------------------------------------------------------------
  __fastcall TMeasurer::~TMeasurer()
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
 TMeasurer ::TMeasurer (const TMeasurer &R)
 {
	mSigV  =  R.mSigV ;
	mSigU = R.mSigU;
	mSigR = R.mSigR ;

	mT = R.mT ;
	mRzv =  R.mRzv ;
	mVzv = R.mVzv ;
	mUzv = R.mUzv ;
	mDelR = R.mDelR ;
	mDelV = R.mDelV ;
	mDelU = R.mDelU ;
	mRadar = R.mRadar ;
	mOEChannel = R.mOEChannel ;
	mbOEChannel =  R.mbOEChannel ;
	mAlf = R.mAlf;

			// для отчета
		mLenMemoryAlloc = R.mLenMemoryAlloc ;
		mQuantPntReport = R.mQuantPntReport ;

		mparrBuff  = NULL;

		if(R.mparrBuff  != NULL)
		{
		if((mparrBuff  = (double*)malloc( QUANT_COLS_BUFF_MEASURER * R.mLenMemoryAlloc * sizeof(double))) != NULL)
		{
		memcpy( mparrBuff ,R.mparrBuff , QUANT_COLS_BUFF_MEASURER * R.mLenMemoryAlloc * sizeof(double));
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
 TMeasurer TMeasurer::operator=(TMeasurer  R)
 {
	mSigV  =  R.mSigV ;
	mSigU = R.mSigU;
	mSigR = R.mSigR ;

	mT = R.mT ;
	mRzv =  R.mRzv ;
	mVzv = R.mVzv ;
	mUzv = R.mUzv ;
	mDelR = R.mDelR ;
	mDelV = R.mDelV ;
	mDelU = R.mDelU ;
	mRadar = R.mRadar ;
	mOEChannel = R.mOEChannel ;
	mbOEChannel =  R.mbOEChannel ;

			// для отчета
		mLenMemoryAlloc = R.mLenMemoryAlloc ;
		mQuantPntReport = R.mQuantPntReport ;

		mparrBuff  = NULL;

		if(R.mparrBuff  != NULL)
		{
		if((mparrBuff  = (double*)malloc( QUANT_COLS_BUFF_MEASURER * R.mLenMemoryAlloc * sizeof(double))) != NULL)
		{
		memcpy( mparrBuff ,R.mparrBuff , QUANT_COLS_BUFF_MEASURER * R.mLenMemoryAlloc * sizeof(double));
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

  // парам конструктор1
 TMeasurer::TMeasurer (const double SigV, const double SigU, const double SigR ,const double T,const double h, wchar_t *pwcharrFoldReport)
 {
	mbOEChannel =  false ;
	mOEChannel = TOEChannel() ;
	mRadar = TRadar ( SigV,  SigU,  SigR ,T, h,  pwcharrFoldReport);
		  // дл отчета
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
	if((mparrBuff  = (double*)malloc(QUANT_COLS_BUFF_MEASURER * mLenMemoryAlloc * sizeof(double))) != NULL)
	{

	memset( mparrBuff,  0, QUANT_COLS_BUFF_MEASURER * mLenMemoryAlloc * sizeof(double));
	}
	else
	{
	ShowMessage(L"Not memory for mparrBuff ") ;
	Abort() ;
	}

	}

 }
  // парам конструктор 2
 TMeasurer::TMeasurer (const double SigV, const double SigU, const double SigR ,const double T,const double  Rzv
	,const double Vzv, const double Uzv ,const double DelR, const double DelV,const double DelU
	,const bool bOEChannel , const TRadar Radar, const TOEChannel OEChannel, wchar_t *pwcharrFoldReport)
 {
	mSigV  = SigV ;
	mSigU = SigU ;
	mSigR = SigR;

	mT = T;
	mRzv = Rzv;; //  замер дальности
	mVzv = Vzv;  // замер по углу V (курс)
	mUzv = Uzv; // замер по углу U (угол места)
	mDelR = DelR;// ошибка измерения дальности
	mDelV = DelV; // ошибка ищзмерения угла V
	mDelU = DelU;
	mbOEChannel = bOEChannel ;
	mRadar = Radar ;
	mOEChannel = OEChannel ;
		  // дл отчета
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
	if((mparrBuff  = (double*)malloc(QUANT_COLS_BUFF_MEASURER * mLenMemoryAlloc * sizeof(double))) != NULL)
	{

	memset( mparrBuff,  0, QUANT_COLS_BUFF_MEASURER * mLenMemoryAlloc * sizeof(double));
	}
	else
	{
	ShowMessage(L"Not memory for mparrBuff ") ;
	Abort() ;
	}

	}

 }

  // парам конструк   3
  TMeasurer::TMeasurer(const double SigV, const double SigU, const double SigR ,const double T,const double h,const double  Rzv
	,const double Vzv, const double Uzv ,const double DelR, const double DelV,const double DelU , wchar_t *pwcharrFoldReport)
 {
	mbOEChannel =  false ;
	mOEChannel = TOEChannel() ;
	mRadar = TRadar ( SigV,  SigU,  SigR ,T, h, Rzv
	, Vzv,  Uzv , DelR,  DelV, DelU , pwcharrFoldReport);
		  // дл отчета
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
	if((mparrBuff  = (double*)malloc(QUANT_COLS_BUFF_MEASURER * mLenMemoryAlloc * sizeof(double))) != NULL)
	{

	memset( mparrBuff,  0, QUANT_COLS_BUFF_MEASURER * mLenMemoryAlloc * sizeof(double));
	}
	else
	{
	ShowMessage(L"Not memory for mparrBuff ") ;
	Abort() ;
	}

	}

 }

 // парам конструктор 4
 TMeasurer::TMeasurer ( const TRadar Radar, const TOEChannel OEChannel, wchar_t *pwcharrFoldReport)
{
	mbOEChannel =  true ;
	mOEChannel = OEChannel;
	mRadar = Radar;
		  // дл отчета
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
	if((mparrBuff  = (double*)malloc(QUANT_COLS_BUFF_MEASURER * mLenMemoryAlloc * sizeof(double))) != NULL)
	{

	memset( mparrBuff,  0, QUANT_COLS_BUFF_MEASURER * mLenMemoryAlloc * sizeof(double));
	}
	else
	{
	ShowMessage(L"Not memory for mparrBuff ") ;
	Abort() ;
	}

	}
	}
 // формирование мзмерения и ошибок измерения
 // valR - Истинн дальность
 // valV - истинн угол V (курс)
 // valU - исттинн угол U (угол места)
 void TMeasurer::createMeasure (const double valR,const double valV, const double valU, const double valT)
 {

	mRadar.createMeasure( valR, valV, valU,  valT);
	if( mbOEChannel && ( ((valT - mOEChannel.mT) >  mOEChannel.mh) )
		 &&   (fabs(mRadar.mVzv) < mOEChannel.mDiagWidth / 2.) && (fabs(mRadar.mUzv) < mOEChannel.mDiagWidth / 2.))
	{

		mOEChannel.createMeasure( valR, valV, valU,  valT)  ;
		uniteDatas() ;
		mAlf = sqrt(0.5);


	}
	else
	{
	   takeRadarData();
	   mAlf = 1 ;
    }


	mT = valT ;
	updateReportData();

 }


 // занесение информации в массивы для отчета
 void TMeasurer::updateReportData()
 {
       if(mpwcharrFoldReport == NULL) return ;
	 if (mQuantPntReport ==  mLenMemoryAlloc)
	 {
	   mLenMemoryAlloc += 2000 ;
	   mparrBuff = (double *)realloc(mparrBuff,QUANT_COLS_BUFF_MEASURER * mLenMemoryAlloc * sizeof(double)) ;

	 }

	   int num0 =  mQuantPntReport * QUANT_COLS_BUFF_MEASURER;


	   mparrBuff [num0] =     mT ;

	   mparrBuff [num0 +1]  =  mRzv;

	   mparrBuff [num0 +2]  =  mVzv ;

	   mparrBuff [num0 +3]  =  mUzv ;

	   mparrBuff [num0 +4]  = mDelR  ;

	   mparrBuff [num0 +5]  = mDelV ;

	   mparrBuff [num0 +6]  = mDelU  ;

		mQuantPntReport++ ;

 }

// публикация  отчета
 void TMeasurer::WriteReport()
 {
	 if (mpwcharrFoldReport == NULL )return ;
	 const int mColsBuff =  QUANT_COLS_BUFF_MEASURER ;
	 wchar_t wcharrPath [300] = {0} ;
	 wcscpy(wcharrPath, mpwcharrFoldReport);
	 wcscat(wcharrPath, L"MeasurerReport\\");


	 wchar_t wcharrFileNames[QUANT_COLS_BUFF_MEASURER * 30] = {0};
	 wcscpy( &wcharrFileNames[ 0 * 30], L"t");
	 wcscpy( &wcharrFileNames[ 1 * 30], L"Rzv");
	 wcscpy( &wcharrFileNames[ 2* 30],  L"Vzv");
	 wcscpy( &wcharrFileNames[ 3 * 30], L"Uzv");
	 wcscpy( &wcharrFileNames[ 4 * 30], L"DelR");
	 wcscpy( &wcharrFileNames[ 5 * 30], L"DelV");
	 wcscpy( &wcharrFileNames[ 6 * 30], L"DelU");

	 double *pscaleY = new double  [QUANT_COLS_BUFF_MEASURER] ;
	 pscaleY[1] = 1.;
	 pscaleY[2] = 1000;
	 pscaleY[3] = 1000.;
	 pscaleY[4] = 1.;
	 pscaleY[5] = 1000.;
	 pscaleY[6] = 1000.;



	 for (int i = 1; i < QUANT_COLS_BUFF_MEASURER; i++)
	 {
	 TYrWriteShapeFile::WriteOneReport(wcharrPath  // путь к папке
								  ,mparrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,QUANT_COLS_BUFF_MEASURER // - к-во переменных о корорых накоплена информация в буфере
								  ,mQuantPntReport  //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30// максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i // номер переменной по оси Y
								  ,100 //  масштаб по оси X
								  ,pscaleY[i]  // масштаб по оси Y
								   )  ;
	 }

	wchar_t wchFileName [300] = {0} ;
	wcscpy(wchFileName, wcharrPath );
	wcscat(wchFileName, L"Axes.shp");
	TYrWriteShapeFile::CreateShpAxes(wchFileName,-40000.,40000,-40000.,40000) ;
	delete pscaleY ;
 }


 // публикация  отчета   перегруженный
 void TMeasurer::WriteReport(wchar_t *pwcharrPath)
 {
	 if (pwcharrPath == NULL )return ;
	 wchar_t wcharrPath [300] = {0} ;
	 wcscpy(wcharrPath, pwcharrPath);
	 wcscat(wcharrPath, L"\\");

	 wchar_t wcharrFileNames[QUANT_COLS_BUFF_MEASURER * 30] = {0};
	 wcscpy( &wcharrFileNames[ 0 * 30], L"t");
	 wcscpy( &wcharrFileNames[ 1 * 30], L"Rzv");
	 wcscpy( &wcharrFileNames[ 2* 30],  L"Vzv");
	 wcscpy( &wcharrFileNames[ 3 * 30], L"Uzv");
	 wcscpy( &wcharrFileNames[ 4 * 30], L"DelR");
	 wcscpy( &wcharrFileNames[ 5 * 30], L"DelV");
	 wcscpy( &wcharrFileNames[ 6 * 30], L"DelU");

	 double *pscaleY = new double  [QUANT_COLS_BUFF_MEASURER] ;
	 pscaleY[1] = 1.;
	 pscaleY[2] = 1000;
	 pscaleY[3] = 1000.;
	 pscaleY[4] = 1.;
	 pscaleY[5] = 1000.;
	 pscaleY[6] = 1000.;



	 for (int i = 1; i < QUANT_COLS_BUFF_MEASURER; i++)
	 {
	 TYrWriteShapeFile::WriteOneReport(wcharrPath  // путь к папке
								  ,mparrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,QUANT_COLS_BUFF_MEASURER // - к-во переменных о корорых накоплена информация в буфере
								  ,mQuantPntReport  //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30// максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i // номер переменной по оси Y
								  ,100 //  масштаб по оси X
								  ,pscaleY[i]  // масштаб по оси Y
								   )  ;
	 }

	wchar_t wchFileName [300] = {0} ;
	wcscpy(wchFileName, wcharrPath );
	wcscat(wchFileName, L"Axes.shp");
	TYrWriteShapeFile::CreateShpAxes(wchFileName,-40000.,40000,-40000.,40000) ;
	delete pscaleY ;
 }

 //
void TMeasurer::takeRadarData()
{
	mDelV = mRadar.mDelV; // ошибка ищзмерения угла V
	mDelU = mRadar.mDelU;


	mVzv = mRadar.mVzv;  // замер по углу V (курс)
	mUzv = mRadar.mUzv; // замер по углу U (угол места)


	mRzv = mRadar.mRzv ;
	mDelR = mRadar.mDelR ;
	mT = mRadar.mT;

	mSigV = mRadar.mSigV ;
	mSigU = mRadar.mSigU ;
	mSigR = mRadar.mSigR ;
}

void TMeasurer::uniteDatas()
{
  uniteTwoMeasures(mRadar.mVzv  , mRadar.mSigV, mRadar.mDelV
			 ,mOEChannel.mVzv  , mOEChannel.mSigV, mOEChannel.mDelV
		   , &mVzv, &mSigV, & mDelV) ;
  uniteTwoMeasures(mRadar.mRzv  , mRadar.mSigR, mRadar.mDelR,
			 mOEChannel.mRzv  , mOEChannel.mSigR, mOEChannel.mDelR
		   , &mRzv, &mSigR, & mDelR) ;
  uniteTwoMeasures(mRadar.mUzv  , mRadar.mSigU, mRadar.mDelU,
			 mOEChannel.mUzv  , mOEChannel.mSigU, mOEChannel.mDelU
		   , &mUzv, &mSigU, & mDelU) ;
  mT = mRadar.mT;


}

void TMeasurer::uniteTwoMeasures(const double valMeasure0, const double valSig0, const double valDel0
  ,const double valMeasure1, const double valSig1,  const double valDel1
  , double *pMeasureRez,  double *pSigRez,  double *pDelRez)
{
   const double valDisp1 = valSig1 * valSig1 ;
   const double valDisp0 = valSig0 * valSig0 ;
   const double alf = valDisp1 / (valDisp0 + valDisp1 ) ;
  *pMeasureRez = alf * valMeasure0 + (1 - alf) * valMeasure1 ;
  *pSigRez =  sqrt(valDisp0 * valDisp1 / (valDisp0 + valDisp1 ) );
  *pDelRez =  alf * valDel0 + (1 - alf) * valDel1 ;
}
#pragma package(smart_init)
