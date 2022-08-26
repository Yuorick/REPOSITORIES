//---------------------------------------------------------------------------


#pragma hdrstop
#include <vcl.h>
#include "DriverMech.h"
#include "Traject.h"
#include <stdlib.h>
#include <dir.h>
#include <string.h>
 #include "YrWriteShapeFile.h"
  #include "Gauss.h"

const int QUANT_COLS_BUFF_DRIVER =  7  ;
TDriverMech::TDriverMech()
{
	mSigBet  =  0.000021 ;
	mSigEps  = 0.000021 ;
	mSigDrBet = 0.003141;
	mSigDrEps = 0.003141;
	mpwcharrFoldReport = NULL ;
	mparrBuff   = NULL ;
	mQuantPntReport = 0;
	mLenMemoryAlloc = 0 ;

}

//---------------------------------------------------------------------------
 __fastcall TDriverMech::~TDriverMech()
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


// ����������� �����������
 TDriverMech ::TDriverMech (const TDriverMech &R)
 {
	mSigBet  =  R.mSigBet ;
	mSigEps = R.mSigEps;
	mSigDrBet = R.mSigDrBet ;
	mSigDrEps = R.mSigDrEps ;
	mTDr = R.mTDr; // ������� ����� �������� ����������
	mEstEps = R.mEstEps;// ������(���������) ���� Eps
	mEstBet = R.mEstBet; // ������(���������) ���� Bet
	mRealEps = R.mRealEps;// �������� ���� Eps
	mRealBet = R.mRealBet; // �������� ���� Bet
	mDelEps = R.mDelEps;// ������ ��  Eps
	mDelBet = R.mDelBet;
		// ��� ������
	mLenMemoryAlloc = R.mLenMemoryAlloc ;
	mQuantPntReport = R.mQuantPntReport ;

	////
	mparrBuff  = NULL;

	if(R.mparrBuff  != NULL)
	{
		if((mparrBuff  = (double*)malloc( 7 * R.mLenMemoryAlloc * sizeof(double))) != NULL)
		{
		memcpy( mparrBuff ,R.mparrBuff , 7 * R.mLenMemoryAlloc * sizeof(double));
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

 // �������� ������������
 TDriverMech &TDriverMech::operator=(const TDriverMech  &R)
 {
	mSigBet  =  R.mSigBet ;
	mSigEps = R.mSigEps;
	mSigDrBet = R.mSigDrBet ;
	mSigDrEps = R.mSigDrEps ;
	mTDr = R.mTDr; // ������� ����� �������� ����������
	mEstEps = R.mEstEps;// ������(���������) ���� Eps
	mEstBet = R.mEstBet; // ������(���������) ���� Bet
	mRealEps = R.mRealEps;// �������� ���� Eps
	mRealBet = R.mRealBet; // �������� ���� Bet
	mDelEps = R.mDelEps;// ������ ��  Eps
	mDelBet = R.mDelBet;
		// ��� ������
	mLenMemoryAlloc = R.mLenMemoryAlloc ;
	mQuantPntReport = R.mQuantPntReport ;

	////
	mparrBuff  = NULL;

	if(R.mparrBuff  != NULL)
	{
		if((mparrBuff  = (double*)malloc( QUANT_COLS_BUFF_DRIVER * R.mLenMemoryAlloc * sizeof(double))) != NULL)
		{
		memcpy( mparrBuff ,R.mparrBuff , QUANT_COLS_BUFF_DRIVER * R.mLenMemoryAlloc * sizeof(double));
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
  // ����� ����������� 1
 TDriverMech  :: TDriverMech ( const double SigBet, const double SigEps
	, const double SigDrBet, const double SigDrEps, wchar_t *pwcharrFoldReport)
 {
	mSigDrEps = SigDrEps ;
	mSigDrBet = SigDrBet ;
	mSigEps = SigEps;
	mSigBet = SigBet ;

	  // �� ������
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
	if((mparrBuff  = (double*)malloc( QUANT_COLS_BUFF_DRIVER * mLenMemoryAlloc * sizeof(double))) != NULL)
	{

	memset( mparrBuff,  0,  QUANT_COLS_BUFF_DRIVER * mLenMemoryAlloc * sizeof(double));
	}
	else
	{
	ShowMessage(L"Not memory for mparrBuff ") ;
	Abort() ;
	}

	}

 }

	 // ����� ����������� 2
 TDriverMech :: TDriverMech ( const double DriverSigBet // �������� ��������� ���� Bet �������
								,const double DriverSigEps // �������� ��������� ���� Eps  ������� (���� �����)
								,const double DriverDynamicSigBet // �������� ��������� ���� �����  �������
								,const double DriverDynamicSigEps // ��������  ������� ��������� ���� �����
								,const double TDr  // ������� ����� �������� ����������
								,const double VAlEstEps // ������(���������) ���� Eps
								,const double VAlEstBet  // ������(���������) ���� Bet
								,const double VAlRealEps // �������� ���� Eps
								,const double VAlRealBet  // �������� ���� Bet
								,const double VAlDelEps // ������ ��  Eps
								,const double VAlDelBet  // ������ �� Bet
								,wchar_t *pwcharrFoldReport)
 {
		mSigDrEps = DriverDynamicSigEps ;
		mSigDrBet = DriverDynamicSigBet ;
		mSigEps = DriverSigEps;
		mSigBet = DriverSigBet;
		mEstEps = VAlEstEps ;// ������(���������) ���� Eps
		mEstBet =  VAlEstBet; // ������(���������) ���� Bet
		mRealEps = VAlRealEps;// �������� ���� Eps
		mRealBet = VAlRealBet; // �������� ���� Bet
		mDelEps = VAlDelEps;// ������ ��  Eps
		mDelBet = VAlDelBet; // ������ �� Bet
		mTDr = TDr;
	  // �� ������
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
	if((mparrBuff  = (double*)malloc( QUANT_COLS_BUFF_DRIVER * mLenMemoryAlloc * sizeof(double))) != NULL)
	{

	memset( mparrBuff,  0,  QUANT_COLS_BUFF_DRIVER * mLenMemoryAlloc * sizeof(double));
	}
	else
	{
	ShowMessage(L"Not memory for mparrBuff ") ;
	Abort() ;
	}

	}

 }


 //--------------------------------------------------
 // ������������� ������� �� ���������� ����
 void TDriverMech :: init (	 const double TDr  // ������� ����� �������� ����������
							,const double VAlEstEps // ������(���������) ���� Eps
							,const double VAlEstBet  // ������(���������) ���� Bet
							,const double VAlRealEps // �������� ���� Eps
							,const double VAlRealBet  // �������� ���� Bet
							,const double VAlDelEps // ������ ��  Eps
							,const double VAlDelBet  // ������ �� Bet
								)
 {

		mEstEps = VAlEstEps ;// ������(���������) ���� Eps
		mEstBet =  VAlEstBet; // ������(���������) ���� Bet
		mRealEps = VAlRealEps;// �������� ���� Eps
		mRealBet = VAlRealBet; // �������� ���� Bet
		mDelEps = VAlDelEps;// ������ ��  Eps
		mDelBet = VAlDelBet; // ������ �� Bet
		mTDr = TDr;


 }

 // ������������� �������
 // INPUT
 // valEps0, valBet0 - ������������
 void TDriverMech::recalcDriver(const double valT,const double valEps0,const double valBet0 )
 {
   if ( valT < mTDr )
	{
		ShowMessage(L"������ � �������� �������� �������") ;
		return ;
	 }
	 mRealEps = valEps0 + getGauss(0, mSigDrEps ); // �������� ���� Eps
	 mRealBet = valBet0 + getGauss(0, mSigDrBet );  // �������� ���� Bet
	 mDelEps = 0.;//getGauss(0, mSigEps );// ������ ��  Eps
	 mDelBet = 0.;//getGauss(0, mSigBet ); // ������ �� Bet
   mEstEps = mRealEps + mDelEps ;
   mEstBet = mRealBet + mDelBet ;
	mTDr =  valT ;
	updateReportData();
 }


 // ��������� ���������� � ������� ��� ������
 void TDriverMech::updateReportData()
 {
     if(mpwcharrFoldReport == NULL) return ;
	 if (mQuantPntReport ==  mLenMemoryAlloc)
	 {
	   mLenMemoryAlloc += 2000 ;
	   mparrBuff = (double *)realloc(mparrBuff,QUANT_COLS_BUFF_DRIVER * mLenMemoryAlloc * sizeof(double)) ;

	 }

	   int num0 =  mQuantPntReport * QUANT_COLS_BUFF_DRIVER;

	   mparrBuff [num0] =     mTDr ;

	   mparrBuff [num0 +1]  =  mEstEps;

	   mparrBuff [num0 +2]  =  mEstBet ;

	   mparrBuff [num0 +3]  =  mRealEps ;

	   mparrBuff [num0 +4]  = mRealBet  ;

	   mparrBuff [num0 +5]  = mDelEps ;

	   mparrBuff [num0 +6]  = mDelBet  ;

		mQuantPntReport++ ;

 }

// ����������  ������
 void TDriverMech::WriteReport()
 {
	 if (mpwcharrFoldReport == NULL )return ;

	 wchar_t wcharrPath [300] = {0} ;
	 wcscpy(wcharrPath, mpwcharrFoldReport);
	 wcscat(wcharrPath, L"\\DriverReport");
		_wmkdir(wcharrPath);
			 wcscat(wcharrPath, L"\\");


	 wchar_t wcharrFileNames[QUANT_COLS_BUFF_DRIVER * 30] = {0};
	 wcscpy( &wcharrFileNames[ 0 * 30], L"t");
	 wcscpy( &wcharrFileNames[ 1 * 30], L"EstEps");
	 wcscpy( &wcharrFileNames[ 2* 30],  L"EstBet");
	 wcscpy( &wcharrFileNames[ 3 * 30], L"RealEps ");
	 wcscpy( &wcharrFileNames[ 4 * 30], L"RealBet");
	 wcscpy( &wcharrFileNames[ 5 * 30], L"DelEps");
	 wcscpy( &wcharrFileNames[ 6 * 30], L"DelBet");



	 for (int i = 1; i < QUANT_COLS_BUFF_DRIVER; i++)
	 {
	 TYrWriteShapeFile::WriteOneReport(wcharrPath  // ���� � �����
								  ,mparrBuff // ������ � ����������� - ������� nBuffRows x nBuffCols
								  ,QUANT_COLS_BUFF_DRIVER // - �-�� ���������� � ������� ��������� ���������� � ������
								  ,mQuantPntReport  //  - �-�� �����
								  ,wcharrFileNames //������� � �������� ���������� - ������� nBuffCols x lenName
								  ,30// ������������ ����� ����� ����������
								  ,0  // ����� ���������� �� ��� X
								  ,i // ����� ���������� �� ��� Y
								  ,100 //  ������� �� ��� Y
								  ,10000  // ������� �� ��� X
								   )  ;
	 }

	wchar_t wchFileName [300] = {0} ;
	wcscpy(wchFileName, wcharrPath );
	wcscat(wchFileName, L"Axes.shp");
	TYrWriteShapeFile::CreateShpAxes(wchFileName,-40000.,40000,-40000.,40000) ;
 }


 // ����������  ������  �������������
 void TDriverMech::WriteReport(wchar_t *pwcharrPath)
 {



	 if (pwcharrPath == NULL )return ;
	 wchar_t wcharrPath [300] = {0} ;
	 wcscpy(wcharrPath, pwcharrPath);
	 wcscat(wcharrPath, L"\\");

	 wchar_t wcharrFileNames[QUANT_COLS_BUFF_DRIVER * 30] = {0};
	 wcscpy( &wcharrFileNames[ 0 * 30], L"t");
	 wcscpy( &wcharrFileNames[ 1 * 30], L"EstEps");
	 wcscpy( &wcharrFileNames[ 2* 30],  L"EstBet");
	 wcscpy( &wcharrFileNames[ 3 * 30], L"RealEps ");
	 wcscpy( &wcharrFileNames[ 4 * 30], L"RealBet");
	 wcscpy( &wcharrFileNames[ 5 * 30], L"DelEps");
	 wcscpy( &wcharrFileNames[ 6 * 30], L"DelBet");



	 for (int i = 1; i < QUANT_COLS_BUFF_DRIVER; i++)
	 {
	 TYrWriteShapeFile::WriteOneReport(wcharrPath  // ���� � �����
								  ,mparrBuff // ������ � ����������� - ������� nBuffRows x nBuffCols
								  ,QUANT_COLS_BUFF_DRIVER // - �-�� ���������� � ������� ��������� ���������� � ������
								  ,mQuantPntReport  //  - �-�� �����
								  ,wcharrFileNames //������� � �������� ���������� - ������� nBuffCols x lenName
								  ,30// ������������ ����� ����� ����������
								  ,0  // ����� ���������� �� ��� X
								  ,i // ����� ���������� �� ��� Y
								  ,100 //  ������� �� ��� Y
								  ,10000  // ������� �� ��� X
								   )  ;
	 }

	wchar_t wchFileName [300] = {0} ;
	wcscpy(wchFileName, wcharrPath );
	wcscat(wchFileName, L"Axes.shp");
	TYrWriteShapeFile::CreateShpAxes(wchFileName,-40000.,40000,-40000.,40000) ;
 }

#pragma package(smart_init)