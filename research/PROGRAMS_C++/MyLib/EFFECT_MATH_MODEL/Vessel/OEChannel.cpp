//---------------------------------------------------------------------------


#pragma hdrstop

#include "OEChannel.h"

#include <vcl.h>

#include "Traject.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
 #include "YrWriteShapeFile.h"
  #include "Gauss.h"
const int QUANT_COLS_BUFF_OEChannel = 7 ;
 extern const double PI;
TOEChannel::TOEChannel()
{
	mSigV =  0.00021;
	mSigU = 0.00042;
	mSigR = 5./3. ;
	mh = 0.1 ;
	mDiagWidth = 1. * PI/3000 ;
	mRzv  = 0;; //  замер дальности
	mVzv  = 0;  // замер по углу V (курс)
	mUzv  = 0; // замер по углу U (угол места)
	mDelR = 0;// ошибка измерения дальности
	mDelV = 0; // ошибка ищзмерения угла V
	mDelU = 0;
	mT = 0. ;
	mpwcharrFoldReport = NULL ;
	mparrBuff   = NULL ;
	mQuantPntReport = 0;
	mLenMemoryAlloc = 0 ;

}

//---------------------------------------------------------------------------
  __fastcall TOEChannel::~TOEChannel()
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
 TOEChannel ::TOEChannel (const TOEChannel &R)
 {
	mSigV  =  R.mSigV ;
	mSigU = R.mSigU;
	mSigR = R.mSigR ;
	mh   = R.mh ;
	mT = R.mT ;
	mRzv =  R.mRzv ;
	mVzv = R.mVzv ;
	mUzv = R.mUzv ;
	mDelR = R.mDelR ;
	mDelV = R.mDelV ;
	mDelU = R.mDelU ;
	mDiagWidth = R.mDiagWidth ;

			// для отчета
		mLenMemoryAlloc = R.mLenMemoryAlloc ;
		mQuantPntReport = R.mQuantPntReport ;

		mparrBuff  = NULL;

		if(R.mparrBuff  != NULL)
		{
		if((mparrBuff  = (double*)malloc( QUANT_COLS_BUFF_OEChannel * R.mLenMemoryAlloc * sizeof(double))) != NULL)
		{
		memcpy( mparrBuff ,R.mparrBuff , QUANT_COLS_BUFF_OEChannel * R.mLenMemoryAlloc * sizeof(double));
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
 TOEChannel TOEChannel::operator=(TOEChannel  R)
 {
	mSigV  =  R.mSigV ;
	mSigU = R.mSigU;
	mSigR = R.mSigR ;
	mh   = R.mh ;
	mT = R.mT ;
	mRzv =  R.mRzv ;
	mVzv = R.mVzv ;
	mUzv = R.mUzv ;
	mDelR = R.mDelR ;
	mDelV = R.mDelV ;
	mDelU = R.mDelU ;
	mDiagWidth = R.mDiagWidth ;

			// для отчета
		mLenMemoryAlloc = R.mLenMemoryAlloc ;
		mQuantPntReport = R.mQuantPntReport ;

		mparrBuff  = NULL;

		if(R.mparrBuff  != NULL)
		{
		if((mparrBuff  = (double*)malloc( QUANT_COLS_BUFF_OEChannel * R.mLenMemoryAlloc * sizeof(double))) != NULL)
		{
		memcpy( mparrBuff ,R.mparrBuff , QUANT_COLS_BUFF_OEChannel * R.mLenMemoryAlloc * sizeof(double));
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
 TOEChannel::TOEChannel (const double SigV, const double SigU, const double SigR ,const double T,const double h
 , const double DiagWidth, wchar_t *pwcharrFoldReport)
 {
	mSigV  = SigV ;
	mSigU = SigU ;
	mSigR = SigR;

	mT = T;
	mh = h ;
	mRzv = 10000; //  замер дальности
	mVzv = 0;  // замер по углу V (курс)
	mUzv = 0; // замер по углу U (угол места)
	mDelR = 0;// ошибка измерения дальности
	mDelV = 0; // ошибка ищзмерения угла V
	mDelU = 0;
	mDiagWidth = DiagWidth ;

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
	if((mparrBuff  = (double*)malloc( QUANT_COLS_BUFF_OEChannel * mLenMemoryAlloc * sizeof(double))) != NULL)
	{

	memset( mparrBuff,  0,  QUANT_COLS_BUFF_OEChannel * mLenMemoryAlloc * sizeof(double));
	}
	else
	{
	ShowMessage(L"Not memory for mparrBuff ") ;
	Abort() ;
	}

	}

 }
  // парам конструктор 2
 TOEChannel::TOEChannel (const double SigV, const double SigU, const double SigR ,const double T,const double h,const double  Rzv
	,const double Vzv, const double Uzv ,const double DelR, const double DelV,const double DelU
	, const double DiagWidth  , wchar_t *pwcharrFoldReport)
 {
	mSigV  = SigV ;
	mSigU = SigU ;
	mSigR = SigR;
	mh   = h ;
	mT = T;
	mRzv = Rzv;; //  замер дальности
	mVzv = Vzv;  // замер по углу V (курс)
	mUzv = Uzv; // замер по углу U (угол места)
	mDelR = DelR;// ошибка измерения дальности
	mDelV = DelV; // ошибка ищзмерения угла V
	mDelU = DelU;
	mDiagWidth = DiagWidth ;
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
	if((mparrBuff  = (double*)malloc(QUANT_COLS_BUFF_OEChannel * mLenMemoryAlloc * sizeof(double))) != NULL)
	{

	memset( mparrBuff,  0, QUANT_COLS_BUFF_OEChannel * mLenMemoryAlloc * sizeof(double));
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
 void TOEChannel::createMeasure (const double valR,const double valV, const double valU, const double valT)
 {

	mDelV = getGauss(0, mSigV ); // ошибка ищзмерения угла V
	mDelU = getGauss(0, mSigU );


	mVzv = valV + mDelV;  // замер по углу V (курс)
	mUzv = valU + mDelU; // замер по углу U (угол места)

	double vR = valR + getGauss(0, mSigR );
	int iR = vR / mSigR ;
	mRzv = ((double)iR) * mSigR + mSigR/ 2. ;
	mDelR =  mRzv  - valR ;
	mT = valT ;
	updateReportData();


 }


 // занесение информации в массивы для отчета
 void TOEChannel::updateReportData()
 {
       if(mpwcharrFoldReport == NULL) return ;
	 if (mQuantPntReport ==  mLenMemoryAlloc)
	 {
	   mLenMemoryAlloc += 2000 ;
	   mparrBuff = (double *)realloc(mparrBuff,QUANT_COLS_BUFF_OEChannel * mLenMemoryAlloc * sizeof(double)) ;

	 }

	   int num0 =  mQuantPntReport * QUANT_COLS_BUFF_OEChannel;


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
 void TOEChannel::WriteReport()
 {
	 if (mpwcharrFoldReport == NULL )return ;
	 const int mColsBuff =  QUANT_COLS_BUFF_OEChannel ;
	 wchar_t wcharrPath [300] = {0} ;
	 wcscpy(wcharrPath, mpwcharrFoldReport);
	 wcscat(wcharrPath, L"OEChannelReport\\");


	 wchar_t wcharrFileNames[QUANT_COLS_BUFF_OEChannel * 30] = {0};
	 wcscpy( &wcharrFileNames[ 0 * 30], L"t");
	 wcscpy( &wcharrFileNames[ 1 * 30], L"Rzv");
	 wcscpy( &wcharrFileNames[ 2* 30],  L"Vzv");
	 wcscpy( &wcharrFileNames[ 3 * 30], L"Uzv");
	 wcscpy( &wcharrFileNames[ 4 * 30], L"DelR");
	 wcscpy( &wcharrFileNames[ 5 * 30], L"DelV");
	 wcscpy( &wcharrFileNames[ 6 * 30], L"DelU");

	 double *pscaleY = new double  [QUANT_COLS_BUFF_OEChannel] ;
	 pscaleY[1] = 1.;
	 pscaleY[2] = 1000;
	 pscaleY[3] = 1000.;
	 pscaleY[4] = 1.;
	 pscaleY[5] = 1000.;
	 pscaleY[6] = 1000.;



	 for (int i = 1; i < QUANT_COLS_BUFF_OEChannel; i++)
	 {
	 TYrWriteShapeFile::WriteOneReport(wcharrPath  // путь к папке
								  ,mparrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,QUANT_COLS_BUFF_OEChannel // - к-во переменных о корорых накоплена информация в буфере
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


 // публикация  отчета   перегруженная
 void TOEChannel::WriteReport(wchar_t *pwcharrPath)
 {
	 if (pwcharrPath == NULL )return ;
	 wchar_t wcharrPath [300] = {0} ;
	 wcscpy(wcharrPath, pwcharrPath);
	 wcscat(wcharrPath, L"\\");


	 wchar_t wcharrFileNames[QUANT_COLS_BUFF_OEChannel * 30] = {0};
	 wcscpy( &wcharrFileNames[ 0 * 30], L"t");
	 wcscpy( &wcharrFileNames[ 1 * 30], L"Rzv");
	 wcscpy( &wcharrFileNames[ 2* 30],  L"Vzv");
	 wcscpy( &wcharrFileNames[ 3 * 30], L"Uzv");
	 wcscpy( &wcharrFileNames[ 4 * 30], L"DelR");
	 wcscpy( &wcharrFileNames[ 5 * 30], L"DelV");
	 wcscpy( &wcharrFileNames[ 6 * 30], L"DelU");

	 double *pscaleY = new double  [QUANT_COLS_BUFF_OEChannel] ;
	 pscaleY[1] = 1.;
	 pscaleY[2] = 1000;
	 pscaleY[3] = 1000.;
	 pscaleY[4] = 1.;
	 pscaleY[5] = 1000.;
	 pscaleY[6] = 1000.;



	 for (int i = 1; i < QUANT_COLS_BUFF_OEChannel; i++)
	 {
	 TYrWriteShapeFile::WriteOneReport(wcharrPath  // путь к папке
								  ,mparrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,QUANT_COLS_BUFF_OEChannel // - к-во переменных о корорых накоплена информация в буфере
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



#pragma package(smart_init)
