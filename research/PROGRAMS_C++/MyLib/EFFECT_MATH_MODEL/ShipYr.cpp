//---------------------------------------------------------------------------


#pragma hdrstop

#include "ShipYr.h"
#include <vcl.h>


#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "MatrixProccess.h"
#include "CalcCorMatrx.h"
#include "YrWriteShapeFile.h"

TShipYr::TShipYr()
{
	// ���������
	memset(marrParral ,0, 3 * sizeof(double)) ;//  ������ ����������


	mTShip = 0.; // ����� �������� ����������� ����������
	memset(marrEstVectSost, 0 , 9 * sizeof(double)); // ������ ����������( ���������, �������� � ��������� � ��� )
	mQ = 0. ;// ���� �����
	mPsi = 0. ; // ���� ������� �����
	mTet = 0. ; // ���� �������� �����
	mVQ = 0. ; //�������� ��������� ���� �����
	mVPsi = 0. ; // �������� ��������� ���� ������� �����
	mVTet = 0.; // �������� ��������� ���� �������� �����

  //	mpwcharrFoldReport = NULL ;
   //	mparrBuff   = NULL ;
   //	mQuantPntReport = 0;
   //	mLenMemoryAlloc = 0 ;

}

//---------------------------------------------------------------------------
 __fastcall TShipYr::~TShipYr()
{
 /*  if (mparrBuff != NULL)
   {
   free( mparrBuff)  ;
   mparrBuff = NULL ;
   }
   if ( mpwcharrFoldReport != NULL)
   {
   free( mpwcharrFoldReport );
   mpwcharrFoldReport = NULL ;
   }
   */
}

// ����������� �����������
 TShipYr ::TShipYr (const TShipYr &R)
 {

	memcpy(marrParral , R.marrParral, 3 * sizeof(double)) ;//  ������ ����������


	// ���������� ��������

	mTShip = R.mTShip; // ����� �������� ����������� ����������

	memcpy(marrEstVectSost,R.marrEstVectSost , 9 * sizeof(double)); // ������  ��������� � ��� )

	mQ = R.mQ ;
	mPsi = R.mPsi ; // ���� ������� �����
	mTet = R.mTet ; // ���� �������� �����
	mVQ = R.mVQ ; //�������� ��������� ���� �����
	mVPsi = R.mVPsi ; // �������� ��������� ���� ������� �����
	mVTet = R.mVTet; // �������� ��������� ���� �������� ������

 }
 // �������� ������������
 TShipYr TShipYr::operator=(TShipYr  R)
 {

	memcpy(marrParral , R.marrParral, 3 * sizeof(double)) ;//  ������ ����������


	// ���������� ��������

	mTShip = R.mTShip; // ����� �������� ����������� ����������

	memcpy(marrEstVectSost,R.marrEstVectSost , 9 * sizeof(double)); // ������  ��������� � ��� )

	mQ = R.mQ ;
	mPsi = R.mPsi ; // ���� ������� �����
	mTet = R.mTet ; // ���� �������� �����
	mVQ = R.mVQ ; //�������� ��������� ���� �����
	mVPsi = R.mVPsi ; // �������� ��������� ���� ������� �����
	mVTet = R.mVTet; // �������� ��������� ���� �������� ������

	return *this ;
 }

  // ����� �����������1
 TShipYr::TShipYr ( double *arrParral,const  double TSh
		 ,const  double Q ,const  double Psi
		 ,const  double Tet,const  double VQ,const  double VPsi, const double VTet
		 , const double VShip, const double ZShip, const double ZVShip )
 {

	memcpy(marrParral, arrParral , 3 * sizeof(double)) ;
	mTShip  =TSh ;
	mQ = Q ;
	mPsi =Psi;
	mTet  =Tet ;
	mVQ  = VQ ;
	mVPsi = VPsi ;
	mVTet = VTet;
	memset(marrEstVectSost, 0, 9 * sizeof(double)) ;
	marrEstVectSost[ 2] =  ZShip;
	marrEstVectSost[ 3] =  VShip * sin (mQ ) ;
	marrEstVectSost[ 4] =  VShip * cos (mQ ) ;
	marrEstVectSost[ 5] =  ZVShip ;

 }

  // ����� ����������� 2
 TShipYr::TShipYr ( double *arrParral)
 {

	memcpy(marrParral, arrParral , 3 * sizeof(double)) ;
	mTShip  = 0. ;
	mQ =  0. ;
	mPsi  = 0. ;
	mTet  = 0. ;
	mVQ  = 0. ;
	mVPsi = 0. ;
	mVTet = 0. ;
	memset(marrEstVectSost, 0, 9 * sizeof(double)) ;

 }

void TShipYr::recalcVS(const  double valT
		 ,const  double Q ,const  double Psi
		 ,const  double Tet,const  double VQ,const  double VPsi, const double VTet
		 , const double VShip, const double ZShip, const double ZVShip )
{
	double h =  valT -  mTShip ;
	if (h < 0.)return ;

	marrEstVectSost[0] += marrEstVectSost[3] * h;
	marrEstVectSost[1] += marrEstVectSost[4] * h;
	marrEstVectSost[2] = ZShip ;
	marrEstVectSost[3] = VShip * sin( Q) ;
	marrEstVectSost[4] = VShip * cos( Q) ;
	marrEstVectSost[5] = ZVShip ;
	mTShip  = valT ;
	mQ = Q ;
	mPsi =Psi;
	mTet  =Tet ;
	mVQ  = VQ ;
	mVPsi = VPsi ;
	mVTet = VTet;


}
// �������� ������������� ��������� ������� - �� ����� �����
// INPUT:
// valT - ������� �������������
// OUTPUT:
// arrVSExtrap - ������������������ ������ ���������, 6-�� ������
void TShipYr::LinExtrap(const  double valT  , TShipYr &ShipYrExtr )
{
	double h =  valT -  mTShip ;
	ShipYrExtr = *this ;
	ShipYrExtr.marrEstVectSost[0] +=  ShipYrExtr.marrEstVectSost[3] * h;
	ShipYrExtr.marrEstVectSost[1] +=  ShipYrExtr.marrEstVectSost[4] * h;
	ShipYrExtr.marrEstVectSost[2] +=  ShipYrExtr.marrEstVectSost[5] * h;
	ShipYrExtr.mQ += ShipYrExtr.mVQ * h;
	ShipYrExtr.mPsi += ShipYrExtr.mVPsi* h;
	ShipYrExtr.mTet += ShipYrExtr.mVTet * h;

}

 /*
  // ����� �����������3
 TShipYr::TShipYr (const double Bearing, const double TargCourse
	, const double TargZenitAng,  const double V, const double H ,
	const double R,const double valT, wchar_t *pwcharrFoldReport)
 {
	*this =  TShipYr() ;
	 mTShip =  valT ;
	 mSins.mTSins = valT ;
	 mDriver.mTDr = valT ;
   //	mTraceFlt = TTraceFlt( Bearing, TargCourse
   //	, TargZenitAng,   V,  H , R, valT);
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
	if((mparrBuff  = (double*)malloc( QUANT_COLS_BUFF_ShipYr * mLenMemoryAlloc * sizeof(double))) != NULL)
	{

	memset( mparrBuff,  0,  QUANT_COLS_BUFF_ShipYr * mLenMemoryAlloc * sizeof(double));
	}
	else
	{
	ShowMessage(L"Not memory for mparrBuff ") ;
	Abort() ;
	}

	}


 }

   // ����� ����������� 4
 TShipYr::TShipYr (const TSins Sins, const TMeasurer Measurer, const TDriverMech Driver , const TEnvironment Environment
		 ,const double Width,const double Length, double *arrPar,const  double MaxQ ,const  double T_Q
		 ,const double MaxPsi,const double T_Psi ,const  double MaxTet
		 ,const double T_Tet,const double MaxVert, const double Q0,const double VVess
		 , double *arrDelt ,const TInitTargData InitTargData, wchar_t *pwcharrFoldReport)
{

	const	double valCoeffEnv = ((double) (Environment.mBallWave))/9.;
	mWidth  =Width ;
	mLength =Length ;
	memcpy(marrParral, arrPar , 3 * sizeof(double)) ;
	mMaxQ =MaxQ * valCoeffEnv ;
	mT_Q  =T_Q;
	mMaxPsi =MaxPsi * valCoeffEnv ;
	mT_Psi = T_Psi ;
	mMaxTet =MaxTet * valCoeffEnv ;
	mT_Tet =T_Tet ;
	mMaxVert = MaxVert * valCoeffEnv ;

	mQ0 =Q0 ;
	mVVess = VVess ;
	mTShip  = InitTargData.mT ;
	memcpy(marrDelt, arrDelt, 4*sizeof(double)) ;

	memset(marrVectSost,0, 9 * sizeof(double)) ;
	mQ = mQ0 +  mMaxQ * cos(marrDelt[0]);
	mVQ =  2. * PI/mT_Q * mMaxQ * sin(marrDelt[0]);

	mPsi = mMaxPsi * cos(marrDelt[1]);
	mVPsi =  2. * PI/mT_Psi * mMaxPsi * sin( marrDelt[1]);

	mTet =  mMaxTet * cos( marrDelt[2]);
	mVTet =  -2. * PI/mT_Tet * mMaxTet * sin( marrDelt[2]);
	// ������ ��������� ������� � ��� ��������
	 marrVectSost[0] = 0;
	 marrVectSost[1] = 0;
	 marrVectSost[2] = mMaxVert * sin( marrDelt[3]) ;
	 marrVectSost[3] = mVVess * sin(mQ) ;
	 marrVectSost[4] = mVVess * cos(mQ) ;
	 marrVectSost[5] = 2. * PI/ mT_Psi * mMaxVert * cos( marrDelt[3]) ;

	 ///
	 //  2.����� ������� ��������� ���� � ���� �� ������ InitTargData.mT

  const double valeps = asin(InitTargData.mH / InitTargData.mR) ;

   double arrVectSostGSK[6] ={0.}, arrVectSost_KGSK [6] = {0} ;
   arrVectSostGSK[0] = InitTargData.mR * cos(valeps) * sin(InitTargData.mBearing) ;
   arrVectSostGSK[1] = InitTargData.mR * cos(valeps) * cos(InitTargData.mBearing) ;
   arrVectSostGSK[2] = InitTargData.mH ;
   arrVectSostGSK[3] = InitTargData.mV * sin(InitTargData.mTargZenitAng) * sin(InitTargData.mTargCourse) ;
   arrVectSostGSK[4] = InitTargData.mV * sin(InitTargData.mTargZenitAng) * cos(InitTargData.mTargCourse) ;
   arrVectSostGSK[5] = InitTargData.mV * cos(InitTargData.mTargZenitAng) ;
   MtrxMinusMatrx(arrVectSostGSK,  marrVectSost,3, 1, arrVectSost_KGSK);

	 	// ������������� SINS

		mSins =  Sins;
		mSins.mTSins = InitTargData.mT ;
		mSins.mDelQ     = 0 ;
		mSins.mDelVQ    = 0;
		mSins.mDelPsi   = 0 ;
		mSins.mDelVPsi  = 0 ;
		mSins.mDelTet   = 0;
		mSins.mDelVTet  = 0 ;
		mSins.mDelH     = 0 ;
		mSins.mDelVH    = 0 ;
		mSins.mDelVVess = 0 ;
		// 20.������ ������� ���� ����� :
		mSins.mEstQ =mQ + mSins.mDelQ;
		// 21.������� �������� ���������� ���� �����:
		mSins.mEstVQ = mVQ + mSins.mDelVQ;
		// 22.������ ���� ������� ����� :
		mSins.mEstPsi = mPsi + mSins.mDelPsi;
		// 23.������ �������� ���������� ���� ������� �����:
		mSins.mEstVPsi = mVPsi + mSins.mDelVPsi;
		// 24.������ ���� ��������� ����� :
		mSins.mEstTet = mTet + mSins.mDelTet;
		// 25. ������ �������� ���������� ���� ��������� �����:
		mSins.mEstVTet  = mVTet + mSins.mDelVTet;
		// 26.������ ������ ������� ������� :
		mSins.mEstH = marrVectSost[2] + mSins.mDelH ;
		// 27.������ �������� ���������  ������ ������� ������� :
		mSins.mEstVH = marrVectSost[5] + mSins.mDelVH;
		// 28.������ �������� �������:
		mSins.mEstVVess = mVVess + mSins.mDelVVess;

		 ///

	// ������������� �������
	// �������� ������� �����
	double arrMu[5] ={0} ;
	arrMu[0] = mQ;
	arrMu[1] = mPsi ;
	arrMu[2] = mTet ;
	double matrPereh_PSK_V_KGSK [9] = {0.}, arrVS_PSK[3] = {0.} ; ;
	calcMatr_ASK_v_KGSK(arrMu,matrPereh_PSK_V_KGSK) ;
	MtrxTranspMultMatrx(matrPereh_PSK_V_KGSK,3, 3, arrVectSost_KGSK,1, arrVS_PSK) ;

		double rr = 0., vRealEps = 0.,vRealBet = 0. ;
		recalcCoord_INTO_Spherical(arrVS_PSK, rr, vRealBet, vRealEps) ;


	mDriver = Driver ;
	mDriver.mTDr =  InitTargData.mT ;
	mDriver.mRealBet =  vRealBet ;
	mDriver.mRealEps = vRealEps ;
	mDriver.mDelEps  = 0 ;
	mDriver.mDelBet = 0 ;
	mDriver.mEstEps =  vRealEps + mDriver.mDelEps ;
	mDriver.mEstBet =  vRealBet + mDriver.mDelBet ;

	// ������������� ����������
	mMeasurer = Measurer ;
	mMeasurer.createMeasure(rr, 0.,  0., InitTargData.mT );
	 // ������ ��������� ������� � ��� - ���������

	 marrEstVectSost[0] = 0;
	 marrEstVectSost[1] = 0;
	 marrEstVectSost[2] = mSins.mEstH ;
	 marrEstVectSost[3] = mSins.mEstVVess * sin( mSins.mEstQ) ;
	 marrEstVectSost[4] = mSins.mEstVVess * cos( mSins.mEstQ) ;
	 marrEstVectSost[5] = mSins.mEstVH ;

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
	if((mparrBuff  = (double*)malloc( QUANT_COLS_BUFF_ShipYr * mLenMemoryAlloc * sizeof(double))) != NULL)
	{

	memset( mparrBuff,  0,  QUANT_COLS_BUFF_ShipYr * mLenMemoryAlloc * sizeof(double));
	}
	else
	{
	ShowMessage(L"Not memory for mparrBuff ") ;
	Abort() ;
	}

	}

}

// �������� ��������� �� ������������� �� � ����������� C�
//  arrInp - ���������� �������� �����
// valR - ����������
// valBet - ���� ����������� �� ������� ������� �� ��� OY �� ������������� �� �������� ����� �� ��������� OXY
// valEps - ���� ����� ������-�������� ����� � �������������� ����������
void TShipYr::recalcCoord_INTO_Spherical(double *arrInp, double &valR, double &valBet, double &valEps)
{
	valR = sqrt ( arrInp[0] * arrInp[0] +arrInp[1] * arrInp[1] +arrInp[2] * arrInp[2] ) ;
	valEps  = asin(arrInp[2] / valR ) ;
	valBet  = asin(arrInp[0]/ sqrt (  arrInp[0] * arrInp[0] +arrInp[1] * arrInp[1])) ;
	if (arrInp[1] < 0.)
	 {
		if (arrInp[0] < 0.)
		{
		  valBet = -PI - valBet ;
		}
		else
		{
           valBet = PI - valBet ;
        }
	 }


}

 void TShipYr::recalcVess(const double valT)
 {
	 const double h = valT - mTShip ;
	 if (h < 0)
	 {
		ShowMessage(L"������ � �������� �������� �������") ;
		return ;
	 }
	 calcAngles(valT) ;
	// �������� ������� ��������� ������� � ��� ���������
	 marrVectSost[0] += marrVectSost[3]*h;
	 marrVectSost[1] += marrVectSost[4]*h;
	 marrVectSost[2] = mMaxVert * sin(2. * PI / mT_Psi*valT + marrDelt[3]) ;
	 marrVectSost[3] = mVVess * sin(mQ) ;
	 marrVectSost[4] = mVVess * cos(mQ) ;
	 marrVectSost[5] = 2. * PI/ mT_Psi * mMaxVert * cos(2. * PI/ mT_Psi*valT + marrDelt[3]) ;
	 mTShip = valT ;

	 ///

	 // �������� ������� ��������� - ������
	 marrEstVectSost[0] += marrEstVectSost[3] * h;
	 marrEstVectSost[1] += marrEstVectSost[4] * h;
	 marrEstVectSost[2] = mSins.mEstH ;
	 marrEstVectSost[3] = mSins.mEstVVess * sin( mSins.mEstQ) ;
	 marrEstVectSost[4] = mSins.mEstVVess * cos( mSins.mEstQ) ;
	 marrEstVectSost[5] = mSins.mEstVH ;
	 updateReportData();
	 ///



 }

 void TShipYr::calcAngles(const double valT)
 {
   mQ = mQ0 +  mMaxQ * cos(2. * PI/mT_Q * valT -  marrDelt[0]);
   mVQ =  -2. * PI/mT_Q * mMaxQ * sin(2. * PI/mT_Q * valT -  marrDelt[0]);

   mPsi = mMaxPsi * cos(2. * PI/mT_Psi * valT -  marrDelt[1]);
   mVPsi =  -2. * PI/mT_Psi * mMaxPsi * sin(2. * PI/mT_Psi * valT -  marrDelt[1]);

   mTet =  mMaxTet * cos(2. * PI/mT_Tet * valT -  marrDelt[2]);
   mVTet =  -2. * PI/mT_Tet * mMaxTet * sin(2. * PI/mT_Tet * valT -  marrDelt[2]);

 }

// ��������� ���������� � ������� ��� ������
 void TShipYr::updateReportData()
 {
     if(mpwcharrFoldReport == NULL) return ;
	 if (mQuantPntReport ==  mLenMemoryAlloc)
	 {
	   mLenMemoryAlloc += 2000 ;
	   mparrBuff = (double *)realloc(mparrBuff,QUANT_COLS_BUFF_ShipYr * mLenMemoryAlloc * sizeof(double)) ;

	 }

	   int num0 =  mQuantPntReport * QUANT_COLS_BUFF_ShipYr;

	   mparrBuff [num0] =     mTShip ;

	   mparrBuff [num0 +1]  =  marrVectSost [0];  // arrShipGSK_X

	   mparrBuff [num0 +2]  =  marrVectSost [1] ; // arrShipGSK_Y

	   mparrBuff [num0 +3]  =  marrVectSost [2] ; // arrShipGSK_Z

	   mparrBuff [num0 +4]  = marrVectSost [3] ;   // arrShipGSK_VX

	   mparrBuff [num0 +5]  = marrVectSost [4];  // arrShipGSK_VY

	   mparrBuff [num0 +6]  = marrVectSost [5]  ; // arrShipGSK_VZ

		mparrBuff [num0 + 7]  =  mVVess;

	   mparrBuff [num0 + 8]  =  mQ; //

	   mparrBuff [num0 + 9]  =  mPsi;

	   mparrBuff [num0 + 10]  = mTet;

	   mparrBuff [num0 + 11]  = mVQ;

	   mparrBuff [num0 + 12]  = mVPsi;

	   mparrBuff [num0 + 13]  = mVTet;

	   mparrBuff [num0 + 14]  =  marrEstVectSost [0];  // arrShipGSK_X

	   mparrBuff [num0 + 15]  =  marrEstVectSost [1] ; // arrShipGSK_Y

	   mparrBuff [num0 + 16]  =  marrEstVectSost [2] ; // arrShipGSK_Z

	   mparrBuff [num0 + 17]  =  marrEstVectSost [3] ;   // arrShipGSK_VX

	   mparrBuff [num0 + 18]  =  marrEstVectSost [4];  // arrShipGSK_VY

	   mparrBuff [num0 + 19]  =  marrEstVectSost [5]  ; // arrShipGSK_VZ

		mQuantPntReport++ ;

 }

// ����������  ������
 void TShipYr::WriteReport()
 {
	 if (mpwcharrFoldReport == NULL )return ;
	 const int mColsBuff =  QUANT_COLS_BUFF_ShipYr ;
	 wchar_t wcharrPath [300] = {0} ;
	 wcscpy(wcharrPath, mpwcharrFoldReport);
	 wcscat(wcharrPath, L"ShipYrReport\\");


	 wchar_t wcharrFileNames[QUANT_COLS_BUFF_ShipYr * 30] = {0};

	 wcscpy( &wcharrFileNames[ 0 * 30], L"t");
	 wcscpy( &wcharrFileNames[ 1 * 30], L"ShipGSK_Real_X");
	 wcscpy( &wcharrFileNames[ 2* 30],  L"ShipGSK_Real_Y");
	 wcscpy( &wcharrFileNames[ 3 * 30], L"ShipGSK_Real_Z");
	 wcscpy( &wcharrFileNames[ 4 * 30], L"ShipGSK_Real_VX");
	 wcscpy( &wcharrFileNames[ 5 * 30], L"ShipGSK_Real_VY");
	 wcscpy( &wcharrFileNames[ 6 * 30], L"ShipGSK_Real_VZ");
	 wcscpy( &wcharrFileNames[ 7* 30], L"VVess_Real");
	 wcscpy( &wcharrFileNames[ 8 * 30], L"Q_Real");
	 wcscpy( &wcharrFileNames[ 9 * 30], L"Psi_Real");
	 wcscpy( &wcharrFileNames[10 * 30], L"Tet_Real");
	 wcscpy( &wcharrFileNames[11 * 30], L"VQ_Real");
	 wcscpy( &wcharrFileNames[12 * 30], L"VPsi_Real");
	 wcscpy( &wcharrFileNames[13 * 30], L"VTet_Real");


	 wcscpy( &wcharrFileNames[14 * 30], L"ShipGSK_Est_X");
	 wcscpy( &wcharrFileNames[15 * 30], L"ShipGSK_Est_Y");
	 wcscpy( &wcharrFileNames[16 * 30], L"ShipGSK_Est_Z");
	 wcscpy( &wcharrFileNames[17* 30],  L"ShipGSK_Est_VX");
	 wcscpy( &wcharrFileNames[18 * 30], L"ShipGSK_Est_VY");
	 wcscpy( &wcharrFileNames[19 * 30], L"ShipGSK_Est_VZ");

	 double *pscaleY = new double  [QUANT_COLS_BUFF_ShipYr] ;
	 for (int i = 0; i < QUANT_COLS_BUFF_ShipYr; i++)  pscaleY[i] = 1.;




	 pscaleY[3]  = 100.;
	 pscaleY[6]  = 100.;
	 pscaleY[9]  = 10000;
	 pscaleY[10] = 10000.;
	 pscaleY[11] = 10000;
	 pscaleY[12] = 10000;
	 pscaleY[13] = 10000.;



	 for (int i = 1; i < QUANT_COLS_BUFF_ShipYr; i++)
	 {
	 TYrWriteShapeFile::WriteOneReport(wcharrPath  // ���� � �����
								  ,mparrBuff // ������ � ����������� - ������� nBuffRows x nBuffCols
								  ,QUANT_COLS_BUFF_ShipYr // - �-�� ���������� � ������� ��������� ���������� � ������
								  ,mQuantPntReport  //  - �-�� �����
								  ,wcharrFileNames //������� � �������� ���������� - ������� nBuffCols x lenName
								  ,30// ������������ ����� ����� ����������
								  ,0  // ����� ���������� �� ��� X
								  ,i // ����� ���������� �� ��� Y
								  ,100 //  ������� �� ��� X
								  ,pscaleY[i]  // ������� �� ��� Y
								   )  ;
	 }

	 // ���������� � ��������� XY
	 TYrWriteShapeFile::WriteOneReport(wcharrPath  // ���� � �����
								  ,mparrBuff // ������ � ����������� - ������� nBuffRows x nBuffCols
								  ,QUANT_COLS_BUFF_ShipYr // - �-�� ���������� � ������� ��������� ���������� � ������
								  ,mQuantPntReport  //  - �-�� �����
								  ,wcharrFileNames //������� � �������� ���������� - ������� nBuffCols x lenName
								  ,30// ������������ ����� ����� ����������
								  ,1  // ����� ���������� �� ��� X
								  ,2 // ����� ���������� �� ��� Y
								  ,1 //  ������� �� ��� X
								  ,1 // ������� �� ��� Y
								   )  ;
   ///

	wchar_t wchFileName [300] = {0} ;
	wcscpy(wchFileName, wcharrPath );
	wcscat(wchFileName, L"Axes.shp");
	TYrWriteShapeFile::CreateShpAxes(wchFileName,-40000.,40000,-40000.,40000) ;
	delete pscaleY ;
 }

 // ������������� ������� �� �����  valT
 // SINS, Driver , Measurer �� ���������������!!!
  void TShipYr::Move(const double valT,const double valStep)
 {
	  int iCirc = (valT - mTShip)/ valStep -1;
  // double valTTemp = mTShip + ((double)iCirc) *valStep  ;

   for (int i = 0; i < iCirc; i++)
   {
	 recalcVess(mTShip + valStep) ;

   }

   recalcVess(valT ) ;


 }

// ������������� ������� ��������� ������� ������ TShipYr �� �����  valTExtr
// INPUT:
// valTExtr - ������ �������, �� ��������� ��������� ���������� �������������
// OUTPUT:
// arrVSVessExtrap - ������������������ ������ ���������
void TShipYr::VSProlong(const double valTExtr, double *arrVSVessExtrap)
{

  // ������������� ������� ��������� ������� �� ������  valT

  double h = valTExtr - mTShip ;
  memcpy ( arrVSVessExtrap, marrEstVectSost, 6 * sizeof(double)) ;
  const double valModVessV = sqrt(marrEstVectSost[3] * marrEstVectSost[3] + marrEstVectSost[4] * marrEstVectSost[4]) ;
  // ������ ����� �������
   double valEst_Q = mSins.mEstQ  ;
   double valEst_dQ_Po_dt = mSins.mEstVQ ;

  const double valQExtr = valEst_Q + valEst_dQ_Po_dt * h ;
  arrVSVessExtrap [2] +=  arrVSVessExtrap [5] * h ;

  if (fabs(valEst_dQ_Po_dt) < 1E-15)
  {
  arrVSVessExtrap [0] +=  valModVessV *  sin (valEst_Q ) * h ;
  arrVSVessExtrap [1] +=  valModVessV *  cos (valEst_Q ) * h ;
  }
  else
  {
  arrVSVessExtrap [0] -=  valModVessV * ( cos (valQExtr ) - cos (valEst_Q) )/ valEst_dQ_Po_dt ; //      ( sin (valEst_Q) *  h +  valEst_dQ_Po_dt * cos (valEst_Q) *  h *h / 2.);
  arrVSVessExtrap [1] +=  valModVessV * ( sin (valQExtr ) - sin (valEst_Q) )/ valEst_dQ_Po_dt ; // ( cos (valEst_Q) *  h  -  valEst_dQ_Po_dt * sin (valEst_Q) *  h *h / 2.);
  arrVSVessExtrap [3] =  valModVessV * sin (valEst_Q + valEst_dQ_Po_dt * h) ;
  arrVSVessExtrap [4] =  valModVessV * cos (valEst_Q + valEst_dQ_Po_dt * h) ;
  }
}


// �������� ������   �� ��������� ����� � ��������� ����  � �������� ��������
//  OUTPUT:
// arrZam_GSK - ������ ������ � ���� - 3-� ������ ������
// arrZam_GSK - ������ ������� ������, 9-�� ������ ������
//
void  TShipYr::GetZamer_IN_KGSK (double *pTZam, double *arrZam_KGSK, double *arrKZam_KGSK )

{
	const double valRZv  = mMeasurer.mRzv ;
	const double  valVZv = mMeasurer.mVzv ;
	const double  valUZv = mMeasurer.mUzv ;
	double arrMu[5] ={0} ;
	double matrPereh_ASK_V_PSK [9] = {0}, matrPereh_PSK_V_KGSK[9] = {0};
	double arrT0[6] = {0},arrT1[6] = {0} ;
	double arrASK [3] ={0.} ;
	arrASK [0] =  valRZv * valVZv ;
	arrASK [1] =  valRZv ;
	arrASK [2] =  valRZv * valUZv ;
	arrMu[0] = mSins.mEstQ;
	arrMu[1] = mSins.mEstPsi ;
	arrMu[2] = mSins.mEstTet ;
	arrMu[3] = mDriver.mEstBet ;
	arrMu[4] = mDriver.mEstEps ;
// 1. �������� ������� ���������
	// ������� � ���-��
	calcMatr_ASK_v_PSK( &arrMu[3], matrPereh_ASK_V_PSK) ;
	MtrxMultMatrx(matrPereh_ASK_V_PSK,3,3, arrASK,1, arrT0) ;
	///
	// �������� � ���-��
	MtrxSumMatrx(arrT0, marrParral, 3, 1, arrT1 ) ;

	///
	// ������� � ����
	calcMatr_PSK_v_KGSK(arrMu, matrPereh_PSK_V_KGSK) ;
	MtrxMultMatrx(matrPereh_PSK_V_KGSK,3,3, arrT1, 1, arrZam_KGSK) ;

	///

	calcSummarizedCorMtrx_ErrMes_In_GSK (arrKZam_KGSK) ;
	*pTZam = mMeasurer.mT ;

}

// �������� ������   �� ��������� ����� � ��������� ���  � �������� ��������
//  OUTPUT:
// arrZam_GSK - ������ ������ � ��� - 3-� ������ ������
// arrKZam_GSK - ������ ������� ������, 9-�� ������ ������
//
void  TShipYr::GetZamer_IN_GSK (double *pTZam, double *arrZam_GSK, double *arrKZam_GSK )

{
  GetZamer_IN_KGSK (pTZam, arrZam_GSK, arrKZam_GSK);
  for (int i = 0; i < 3; i++) arrZam_GSK [i] +=  marrEstVectSost [i] ;



}

// ������ �������������� �������������� ������� ������ ���������� ������, ����������� ������
 // �������������� ������ ,  ������ ������� � ������ ����������� ����� � ���
// ��� ������ ����������� ��� ������
 void  TShipYr::calcSummarizedCorMtrx_ErrMes_In_GSK (double *arrCorrMtrx_GSK )
{
   double ar_dF_po_dBet_sq [9],ar_dF_po_dEps_sq [9], matrPereh_ASK_V_KGSK[9],arrMu[5] ;
	   double  arrDT0[9] = {0}, arrDT1[9] = {0}, arrDT2[9] = {0}
				,arrDT3[9] = {0},arrDT4[9] = {0},arrDT5[9] = {0},arrDT51[9] = {0};
	arrMu[0] = mSins.mEstQ;
	arrMu[1] = mSins.mEstPsi ;
	arrMu[2] = mSins.mEstTet ;
	arrMu[3] = mDriver.mEstBet ;
	arrMu[4] = mDriver.mEstEps ;
	   calc_dF_po_dBet_sq(mDriver.mRealBet,  mDriver.mRealEps,ar_dF_po_dBet_sq ) ;
	   calc_dF_po_dEps_sq(mDriver.mRealBet,  mDriver.mRealEps,ar_dF_po_dEps_sq ) ;
	   calcMatr_ASK_v_KGSK(arrMu,matrPereh_ASK_V_KGSK) ;
	   arrDT0[0] = mMeasurer.mSigV * mMeasurer.mSigV * mMeasurer.mRzv * mMeasurer.mRzv ;
	   arrDT0[4] = mMeasurer.mSigR * mMeasurer.mSigR ;
	   arrDT0[8] = mMeasurer.mSigU * mMeasurer.mSigU * mMeasurer.mRzv * mMeasurer.mRzv ;
	   MtrxMultMatrx(matrPereh_ASK_V_KGSK,3, 3, arrDT0,3, arrDT1) ;
	   MtrxMultMatrxTransp(arrDT1,3, 3, matrPereh_ASK_V_KGSK, 3, arrDT2) ;

	   double temp0 =  mMeasurer.mRzv * mMeasurer.mRzv
	   * mDriver.mSigBet* mDriver.mSigBet ;
	   MatrxMultScalar(ar_dF_po_dBet_sq , 3, 3, temp0,arrDT3);
	   MtrxSumMatrx(arrDT2, arrDT3,3,3, arrDT4) ;

	   double temp1 =  mMeasurer.mRzv * mMeasurer.mRzv
	   * mDriver.mSigEps* mDriver.mSigEps ;
	   MatrxMultScalar(ar_dF_po_dEps_sq , 3, 3, temp1,arrDT5);
	   MtrxSumMatrx(arrDT4, arrDT5,3,3, arrDT51) ;

	   // ������ ������, ��������� �������� ����������� �����
	   double ar_dF_po_dTet_sq [9] = {0},ar_dF_po_dQ_sq [9] = {0},ar_dF_po_dPsi_sq [9] = {0}, arrDT6[9] = {0}
		  , arrDT7[9]  = {0}, arrDT8[9] = {0}, arrDT9[9] = {0}, arrDT10[9]  = {0};
	   calc_dF_po_dQ_sq  (mDriver.mEstBet,  mDriver.mEstEps, ar_dF_po_dQ_sq) ;
	   calc_dF_po_dPsi_sq(mDriver.mEstBet,  mDriver.mEstEps, ar_dF_po_dPsi_sq) ;
	   calc_dF_po_dTet_sq(mDriver.mEstBet,  mDriver.mEstEps, ar_dF_po_dTet_sq) ;

	   double temp2 = mMeasurer.mRzv * mMeasurer.mRzv
	   * mSins.mMaxSig_Q * mSins.mMaxSig_Q  ;
	   double temp3 = mMeasurer.mRzv * mMeasurer.mRzv
	   * mSins.mMaxSig_Psi* mSins.mMaxSig_Psi ;
	   double temp4 = mMeasurer.mRzv * mMeasurer.mRzv
	   * mSins.mMaxSig_Tet* mSins.mMaxSig_Tet  ;

	   MatrxMultScalar(ar_dF_po_dQ_sq   , 3, 3, temp2,arrDT6);
	   MatrxMultScalar(ar_dF_po_dPsi_sq , 3, 3, temp3,arrDT7);
	   MatrxMultScalar(ar_dF_po_dTet_sq , 3, 3, temp4,arrDT8);

	   MtrxSumMatrx(arrDT51, arrDT6,3,3, arrDT9) ;
	   MtrxSumMatrx(arrDT9, arrDT7,3,3, arrDT10) ;
	   MtrxSumMatrx(arrDT10, arrDT8,3,3, arrCorrMtrx_GSK) ;
 }

 */

#pragma package(smart_init)
