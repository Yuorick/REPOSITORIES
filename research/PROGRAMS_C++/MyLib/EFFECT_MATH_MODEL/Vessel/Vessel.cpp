//---------------------------------------------------------------------------


#pragma hdrstop
#include <vcl.h>

#include "Vessel.h"
#include <string.h>
#include <stdlib.h>
#include <dir.h>
#include <math.h>
#include "MatrixProccess.h"
#include "CalcCorMatrx.h"
 #include "YrWriteShapeFile.h"
 #include "Traject.h"
 #include "Gauss.h"
 #include "Zamer.h"
 #include "EtalonSign.h"

#define KKK 0.
//#define KKK 1.

const int QUANT_COLS_BUFF_VESSEL = 20;
//---------------------------------------------------------------------------

TVessel::TVessel()
{
	mSins  =  TSins();
	mFar_2D = TFar_2D();
	mTransmitAnt = TTransmitAnt();
	mDriver = TDriverMech() ;
	mShellBody = TShellBody();
	mArtComplex = TArtComplex();
	memset(marrArtParral ,0, 3 * sizeof(double)) ;//  ������ ���������� ��

  //	mTraceFlt = TTraceFlt() ;
	// ���������
	mWidth = 0; // ������(�)
	mLength = 0; // �����(�)
	memset(marrParral ,0, 3 * sizeof(double)) ;//  ������ ����������
	mMaxQ = 3./180.*M_PI; // ������������ ���� ���������� �� ������������ �����(��������� ���� ��������)
	mT_Q = 18.; // ������ ��������
	mMaxPsi = 3./180.*M_PI;// ������������ ���� ������� �����(���������)
	mT_Psi = 12; // ������ ������� �����
	mMaxTet = 12./180.*M_PI; //������������ ���� ��������� �����(���������)
	mT_Tet = 6; // ������ �������� �����
	mMaxVert = 1. ;

	// ���������� ��������
	mQ0 =0. ; // ����������� ����
	mVVess = 0.514 * 20 ;// �������� ������� ������ 20 �����
	mTVess = 0.; // ����� �������� ����������� ����������
	memset(marrVectSost, 0 , 9 * sizeof(double)); // ������ ����������( ���������, �������� � ��������� � ��� )
	memset(marrEstVectSost, 0 , 9 * sizeof(double)); // ������ ����������( ���������, �������� � ��������� � ��� )
	marrVectSost[4] = mVVess ;
	mQ = 0. ;// ���� �����
	mPsi = 0. ; // ���� ������� �����
	mTet = 0. ; // ���� �������� �����
	mVQ = 0. ; //�������� ��������� ���� �����
	mVPsi = 0. ; // �������� ��������� ���� ������� �����
	mVTet = 0.; // �������� ��������� ���� �������� ������

	memset(marrDelt, 0, 4 * sizeof(double)) ;

	mpwcharrFoldReport = NULL ;
	mparrBuff   = NULL ;
	mQuantPntReport = 0;
	mLenMemoryAlloc = 0 ;


	//��������� ��������� ������ ������� � ��� �� 100 �
	mAmp_AftFlexure = 0.;
	// ������ ��������� ��������� ������
	mT_AftFlexure = 4.;
	//��������� ��������� ������ ������� � ��� �� 100 �
	mAmp_BoardFlexure = 0.;
	// ������ ��������� ��������� ������
	mT_BoardFlexure = 2.;

	// ���� ��������� ���� ��������� ��������� ������ ������� �������
	mPhase0_AftFlexure = 0.;
		// ���� ��������� ���� ��������� ��������� ������ ������� �������
	mPhase0_BoardFlexure = 0.;

	// ���� ���������, �
	 mControlSyst = TControlSyst();


}

//---------------------------------------------------------------------------
 __fastcall TVessel::~TVessel()
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
 TVessel ::TVessel (const TVessel &R)
 {
	mSins  =  R.mSins ;
	mFar_2D = R.mFar_2D;
	mTransmitAnt = R.mTransmitAnt;
	mDriver = R.mDriver ;
	mShellBody = R.mShellBody ;

	mArtComplex = R.mArtComplex;
	memcpy(marrArtParral , R.marrArtParral, 3 * sizeof(double)) ;//  ������ ����������

	mWidth = R.mWidth; // ������(�)
	mLength = R.mLength; // �����(�)
	memcpy(marrParral , R.marrParral, 3 * sizeof(double)) ;//  ������ ����������
	mMaxQ = R.mMaxQ; // ������������ ���� ���������� �� ������������ �����(��������� ���� ��������)
	mT_Q = R.mT_Q; // ������ ��������
	mMaxPsi = R.mMaxPsi;// ������������ ���� ������� �����(���������)
	mT_Psi = R.mT_Psi; // ������ ������� �����
	mMaxTet = R.mMaxTet; //������������ ���� ��������� �����(���������)
	mT_Tet = R.mT_Tet; // ������ �������� �����
	mMaxVert = R.mMaxVert ;

	// ���������� ��������
	mQ0 =R.mQ0 ; // ����������� ����
	mVVess = R.mVVess ;
	mTVess = R.mTVess; // ����� �������� ����������� ����������
	memcpy(marrVectSost,R.marrVectSost , 9 * sizeof(double)); // ������  ��������� � ��� )
	memcpy(marrEstVectSost,R.marrEstVectSost , 9 * sizeof(double)); // ������  ��������� � ��� )

	mQ = R.mQ ;
	mPsi = R.mPsi ; // ���� ������� �����
	mTet = R.mTet ; // ���� �������� �����
	mVQ = R.mVQ ; //�������� ��������� ���� �����
	mVPsi = R.mVPsi ; // �������� ��������� ���� ������� �����
	mVTet = R.mVTet; // �������� ��������� ���� �������� ������

		//��������� ��������� ������ ������� � ��� �� 100 �
	mAmp_AftFlexure = R.mAmp_AftFlexure;
	// ������ ��������� ��������� ������
	mT_AftFlexure = R.mT_AftFlexure ;
	//��������� ��������� ������ ������� � ��� �� 100 �
	mAmp_BoardFlexure = R.mAmp_BoardFlexure;
	// ������ ��������� ��������� ������
	mT_BoardFlexure = R.mT_BoardFlexure;
	// ���� ��������� ���� ��������� ��������� ������ ������� �������
	mPhase0_AftFlexure = R.mPhase0_AftFlexure;
		// ���� ��������� ���� ��������� ��������� ������ ������� �������
	mPhase0_BoardFlexure = R.mPhase0_BoardFlexure;

	mControlSyst = R.mControlSyst;

	memcpy(marrDelt, R.marrDelt, 4 * sizeof(double)) ;

			// ��� ������
		mLenMemoryAlloc = R.mLenMemoryAlloc ;
		mQuantPntReport = R.mQuantPntReport ;

		mparrBuff  = NULL;

		if(R.mparrBuff  != NULL)
		{
		if((mparrBuff  = (double*)malloc( QUANT_COLS_BUFF_VESSEL * R.mLenMemoryAlloc * sizeof(double))) != NULL)
		{
		memcpy( mparrBuff ,R.mparrBuff , QUANT_COLS_BUFF_VESSEL * R.mLenMemoryAlloc * sizeof(double));
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
 TVessel &TVessel::operator=(const TVessel  &R)
 {
	mSins  =  R.mSins ;
	mFar_2D = R.mFar_2D;
	mTransmitAnt = R.mTransmitAnt;
	mDriver = R.mDriver ;
	mShellBody = R.mShellBody ;

	mArtComplex = R.mArtComplex;
	memcpy(marrArtParral , R.marrArtParral, 3 * sizeof(double)) ;//  ������ ����������

	mWidth = R.mWidth; // ������(�)
	mLength = R.mLength; // �����(�)
	memcpy(marrParral , R.marrParral, 3 * sizeof(double)) ;//  ������ ����������
	mMaxQ = R.mMaxQ; // ������������ ���� ���������� �� ������������ �����(��������� ���� ��������)
	mT_Q = R.mT_Q; // ������ ��������
	mMaxPsi = R.mMaxPsi;// ������������ ���� ������� �����(���������)
	mT_Psi = R.mT_Psi; // ������ ������� �����
	mMaxTet = R.mMaxTet; //������������ ���� ��������� �����(���������)
	mT_Tet = R.mT_Tet; // ������ �������� �����
	mMaxVert = R.mMaxVert ;

	// ���������� ��������
	mQ0 =R.mQ0 ; // ����������� ����
	mVVess = R.mVVess ;
	mTVess = R.mTVess; // ����� �������� ����������� ����������
	memcpy(marrVectSost,R.marrVectSost , 9 * sizeof(double)); // ������  ��������� � ��� )
	memcpy(marrEstVectSost,R.marrEstVectSost , 9 * sizeof(double)); // ������  ��������� � ��� )

	mQ = R.mQ ;
	mPsi = R.mPsi ; // ���� ������� �����
	mTet = R.mTet ; // ���� �������� �����
	mVQ = R.mVQ ; //�������� ��������� ���� �����
	mVPsi = R.mVPsi ; // �������� ��������� ���� ������� �����
	mVTet = R.mVTet; // �������� ��������� ���� �������� ������

		//��������� ��������� ������ ������� � ��� �� 100 �
	mAmp_AftFlexure = R.mAmp_AftFlexure;
	// ������ ��������� ��������� ������
	mT_AftFlexure = R.mT_AftFlexure ;
	//��������� ��������� ������ ������� � ��� �� 100 �
	mAmp_BoardFlexure = R.mAmp_BoardFlexure;
	// ������ ��������� ��������� ������
	mT_BoardFlexure = R.mT_BoardFlexure;
	// ���� ��������� ���� ��������� ��������� ������ ������� �������
	mPhase0_AftFlexure = R.mPhase0_AftFlexure;
		// ���� ��������� ���� ��������� ��������� ������ ������� �������
	mPhase0_BoardFlexure = R.mPhase0_BoardFlexure;

	mControlSyst = R.mControlSyst;

	memcpy(marrDelt, R.marrDelt, 4 * sizeof(double)) ;

			// ��� ������
		mLenMemoryAlloc = R.mLenMemoryAlloc ;
		mQuantPntReport = R.mQuantPntReport ;

		mparrBuff  = NULL;

		if(R.mparrBuff  != NULL)
		{
		if((mparrBuff  = (double*)malloc( QUANT_COLS_BUFF_VESSEL * R.mLenMemoryAlloc * sizeof(double))) != NULL)
		{
		memcpy( mparrBuff ,R.mparrBuff , QUANT_COLS_BUFF_VESSEL * R.mLenMemoryAlloc * sizeof(double));
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

  // ����� �����������1
 TVessel::TVessel (const TSins Sins, const TFar_2D  Far_2D, const TDriverMech Driver
		 ,const double Width,const double Length, double *arrPar,const  double MaxQ ,const  double T_Q
		 ,const double MaxPsi,const double T_Psi ,const  double MaxTet
		 ,const double T_Tet,const double MaxVert, const double Q0,const double VVess,const  double TVess
		 , double *arrVectSost, double *arrEstVectSost,const  double Q ,const  double Psi
		 ,const  double Tet,const  double VQ,const  double VPsi, const double VTet,  wchar_t *pwcharrFoldReport)
 {
	mSins  = Sins ;
	mFar_2D = Far_2D ;
	mDriver = Driver;

	mWidth  =Width ;
	mLength =Length ;
	memcpy(marrParral, arrPar , 3 * sizeof(double)) ;
	mMaxQ =MaxQ  ;
	mT_Q  =T_Q;
	mMaxPsi =MaxPsi ;
	mT_Psi = T_Psi ;
	mMaxTet =MaxTet ;
	mT_Tet =T_Tet ;
	mMaxVert = MaxVert ;

	mQ0 =Q0 ;
	mVVess = VVess ;
	mTVess  =TVess ;
	memcpy(marrVectSost,arrVectSost, 9 * sizeof(double)) ;
	memcpy(marrEstVectSost,arrEstVectSost, 9 * sizeof(double)) ;
	mQ = Q ;
	mPsi =Psi;
	mTet  =Tet ;
	mVQ  =VQ ;
	mVPsi =VPsi ;
	mVTet = VTet;

	for (int i = 0 ; i  <4; i++)
	{
	 marrDelt [i] = getRand01() * 2. * M_PI;
	}

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
	if((mparrBuff  = (double*)malloc( QUANT_COLS_BUFF_VESSEL * mLenMemoryAlloc * sizeof(double))) != NULL)
	{

	memset( mparrBuff,  0,  QUANT_COLS_BUFF_VESSEL * mLenMemoryAlloc * sizeof(double));
	}
	else
	{
	ShowMessage(L"Not memory for mparrBuff ") ;
	Abort() ;
	}

	}
	mArtComplex = TArtComplex();
	memset(marrArtParral ,0, 3 * sizeof(double)) ;//  ������ ���������� ��


 }


   // ����� ����������� 4
 TVessel::TVessel ( 	const TFar_2D  Far_2D,
							const double DriverSigBet // �������� ��������� ���� Bet �������
								,const double DriverSigEps // �������� ��������� ���� Eps  ������� (���� �����)
								,const double DriverDynamicSigBet // �������� ��������� ���� �����  �������
								,const double DriverDynamicSigEps // ��������  ������� ��������� ���� �����
								,const double MaxSig_Q, const double MaxSig_Psi, const double MaxSig_Tet      // ����
								,const double MaxSig_dQdt, const double MaxSig_dPsidt,const double MaxSig_dTetdt   // ����
								,const double MaxSig_H,const double MaxSig_VH,const double K1,const double SigV   // ����
								,const TEnvironment Environment
								,const double Width,const double Length, double *arrPar,const  double MaxQ ,const  double T_Q
								,const double MaxPsi,const double T_Psi ,const  double MaxTet
								,const double T_Tet,const double MaxVert, const double Q0,const double VVess
								,const TInitTargData InitTargData
								,const double MaxAmp_AftFlexure,const double T_AftFlexure,const double MaxAmp_BoardFlexure
								,const double T_BoardFlexure,wchar_t *pwcharrFoldReport)
{
	 // ������������� ����������
 	mFar_2D = Far_2D ;

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

	mAmp_AftFlexure = MaxAmp_AftFlexure * valCoeffEnv ;
	mT_AftFlexure = T_AftFlexure;
	mAmp_BoardFlexure = MaxAmp_BoardFlexure * valCoeffEnv ;
	mT_BoardFlexure = T_BoardFlexure;
	mPhase0_AftFlexure = KKK * getRand01() * 2. * M_PI;
	mPhase0_BoardFlexure = KKK * getRand01() * 2. * M_PI;


	mQ0 =Q0 ;
	mVVess = VVess ;
	mTVess  = InitTargData.mT ;
	for (int i = 0 ; i < 4; i++)
	{
	 marrDelt [i] = KKK *getRand01() * 2. * M_PI;
	}

	memset(marrVectSost,0, 9 * sizeof(double)) ;
	mQ = mQ0 +  mMaxQ * cos(marrDelt[0]);
	mVQ =  2. * M_PI/mT_Q * mMaxQ * sin(marrDelt[0]);

	mPsi = mMaxPsi * cos(marrDelt[1]);
	mVPsi =  2. * M_PI/mT_Psi * mMaxPsi * sin( marrDelt[1]);

	mTet =  mMaxTet * cos( marrDelt[2]);
	mVTet =  -2. * M_PI/mT_Tet * mMaxTet * sin( marrDelt[2]);
	// ������ ��������� ������� � ��� ��������
	 marrVectSost[0] = 0;
	 marrVectSost[1] = 0;
	 marrVectSost[2] = mMaxVert * sin( marrDelt[3]) ;
	 marrVectSost[3] = mVVess * sin(mQ) ;
	 marrVectSost[4] = mVVess * cos(mQ) ;
	 marrVectSost[5] = 2. * M_PI/ mT_Psi * mMaxVert * cos( marrDelt[3]) ;

	 ///
	 //  2.����� ������� ��������� ���� � ���� �� ������ InitTargData.mT
	 if (fabs(InitTargData.mH / InitTargData.mR) > 1.)
{
	ShowMessage(L"ERROR");
	int iii = 0;


}
	const double valeps = arcSin(InitTargData.mH / InitTargData.mR) ;

   double arrVectSostGSK[6] ={0.}, arrVectSost_KGSK [6] = {0} ;
   arrVectSostGSK[0] = InitTargData.mR * cos(valeps) * sin(InitTargData.mBearing) ;
   arrVectSostGSK[1] = InitTargData.mR * cos(valeps) * cos(InitTargData.mBearing) ;
   arrVectSostGSK[2] = InitTargData.mH ;
   arrVectSostGSK[3] = InitTargData.mV * sin(InitTargData.mTargZenitAng) * sin(InitTargData.mTargCourse) ;
   arrVectSostGSK[4] = InitTargData.mV * sin(InitTargData.mTargZenitAng) * cos(InitTargData.mTargCourse) ;
   arrVectSostGSK[5] = InitTargData.mV * cos(InitTargData.mTargZenitAng) ;
   MtrxMinusMatrx(arrVectSostGSK,  marrVectSost,3, 1, arrVectSost_KGSK);

	 	// ������������� SINS
		 mSins = TSins ( Environment
	,    MaxSig_Q,     MaxSig_Psi,     MaxSig_Tet
	,    MaxSig_dQdt,     MaxSig_dPsidt,    MaxSig_dTetdt
	,    MaxSig_H,    MaxSig_VH,    K1,    SigV
	,     mTVess     // �����
	,    mVVess   // ��� ��������
	,     mQ  // �������� ���� �����
	,     mPsi  // �������� ���� ������� �����
	,     mTet  // �������� ���� �������� �����
	,     mVQ  //�������� �������� ��������� ���� �����
	,     mVPsi   //�������� �������� ��������� ���� ������� �����
	,     mVTet   //�������� �������� ��������� ���� �������� ������
	,    marrVectSost[2]
	,    marrVectSost[5]
	,    mT_Q
	,     mT_Psi
	,    mT_Tet
	,   marrDelt
	,pwcharrFoldReport) ;

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



	const double VAlDelEps  = 0 ;
	const double VAlDelBet = 0 ;
	const double VAlEstEps =  vRealEps + VAlDelEps  ;
	const double VAlEstBet =  vRealBet + VAlDelBet ;
		 // ����� ����������� 2
 mDriver = TDriverMech (      DriverSigBet // �������� ��������� ���� Bet �������
								,    DriverSigEps // �������� ��������� ���� Eps  ������� (���� �����)
								,    DriverDynamicSigBet // �������� ��������� ���� �����  �������
								,    DriverDynamicSigEps // ��������  ������� ��������� ���� �����
								,InitTargData.mT  // ������� ����� �������� ����������
								,    VAlEstEps // ������(���������) ���� Eps
								,    VAlEstBet  // ������(���������) ���� Bet
								,     vRealEps // �������� ���� Eps
								,    vRealBet  // �������� ���� Bet
								,    VAlDelEps // ������ ��  Eps
								,    VAlDelBet  // ������ �� Bet
								,   pwcharrFoldReport) ;


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
	if((mparrBuff  = (double*)malloc( QUANT_COLS_BUFF_VESSEL * mLenMemoryAlloc * sizeof(double))) != NULL)
	{

	memset( mparrBuff,  0,  QUANT_COLS_BUFF_VESSEL * mLenMemoryAlloc * sizeof(double));
	}
	else
	{
	ShowMessage(L"Not memory for mparrBuff ") ;
	Abort() ;
	}

	}
	mArtComplex = TArtComplex();
	memset(marrArtParral ,0, 3 * sizeof(double)) ;//  ������ ���������� ��

}



   // ����� ����������� 5
 TVessel::TVessel ( const TShellBody ShellBody
								  ,const TFar_2D  Far_2D
								  ,const double DriverSigBet // �������� ��������� ���� Bet �������
								,const double DriverSigEps // �������� ��������� ���� Eps  ������� (���� �����)
								,const double DriverDynamicSigBet // �������� ��������� ���� �����  �������
								,const double DriverDynamicSigEps // ��������  ������� ��������� ���� �����
								,const double MaxSig_Q, const double MaxSig_Psi, const double MaxSig_Tet      // ����
								,const double MaxSig_dQdt, const double MaxSig_dPsidt,const double MaxSig_dTetdt   // ����
								,const double MaxSig_H,const double MaxSig_VH,const double K1,const double SigV   // ����
								,const TEnvironment Environment
								,const double Width,const double Length, double *arrPar,const  double MaxQ ,const  double T_Q
								,const double MaxPsi,const double T_Psi ,const  double MaxTet
								,const double T_Tet,const double MaxVert, const double Q0,const double VVess
								,const TInitTargData InitTargData
								,const double MaxAmp_AftFlexure,const double T_AftFlexure,const double MaxAmp_BoardFlexure
								,const double T_BoardFlexure,wchar_t *pwcharrFoldReport)
{
	mShellBody = ShellBody;
	 // ������������� ����������
 	mFar_2D = Far_2D ;

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

	mAmp_AftFlexure = MaxAmp_AftFlexure * valCoeffEnv ;
	mT_AftFlexure = T_AftFlexure;
	mAmp_BoardFlexure = MaxAmp_BoardFlexure * valCoeffEnv ;
	mT_BoardFlexure = T_BoardFlexure;
	mPhase0_AftFlexure = KKK * getRand01() * 2. * M_PI;
	mPhase0_BoardFlexure = KKK * getRand01() * 2. * M_PI;


	mQ0 =Q0 ;
	mVVess = VVess ;
	mTVess  = InitTargData.mT ;
	for (int i = 0 ; i < 4; i++)
	{
	 marrDelt [i] = KKK * getRand01() * 2. * M_PI;
	}

	memset(marrVectSost,0, 9 * sizeof(double)) ;
	mQ = mQ0 +  mMaxQ * cos(marrDelt[0]);
	mVQ =  2. * M_PI/mT_Q * mMaxQ * sin(marrDelt[0]);

	mPsi = mMaxPsi * cos(marrDelt[1]);
	mVPsi =  2. * M_PI/mT_Psi * mMaxPsi * sin( marrDelt[1]);

	mTet =  mMaxTet * cos( marrDelt[2]);
	mVTet =  -2. * M_PI/mT_Tet * mMaxTet * sin( marrDelt[2]);
	// ������ ��������� ������� � ��� ��������
	 marrVectSost[0] = 0;
	 marrVectSost[1] = 0;
	 marrVectSost[2] = mMaxVert * sin( marrDelt[3]) ;
	 marrVectSost[3] = mVVess * sin(mQ) ;
	 marrVectSost[4] = mVVess * cos(mQ) ;
	 marrVectSost[5] = 2. * M_PI/ mT_Psi * mMaxVert * cos( marrDelt[3]) ;

	 ///
	 //  2.����� ������� ��������� ���� � ���� �� ������ InitTargData.mT

   	 if (fabs(InitTargData.mH / InitTargData.mR) > 1.)
{
	ShowMessage(L"ERROR");
	int iii = 0;


}


	const double valeps = arcSin(InitTargData.mH / InitTargData.mR) ;

   double arrVectSostGSK[6] ={0.}, arrVectSost_KGSK [6] = {0} ;
   arrVectSostGSK[0] = InitTargData.mR * cos(valeps) * sin(InitTargData.mBearing) ;
   arrVectSostGSK[1] = InitTargData.mR * cos(valeps) * cos(InitTargData.mBearing) ;
   arrVectSostGSK[2] = InitTargData.mH ;
   arrVectSostGSK[3] = InitTargData.mV * sin(InitTargData.mTargZenitAng) * sin(InitTargData.mTargCourse) ;
   arrVectSostGSK[4] = InitTargData.mV * sin(InitTargData.mTargZenitAng) * cos(InitTargData.mTargCourse) ;
   arrVectSostGSK[5] = InitTargData.mV * cos(InitTargData.mTargZenitAng) ;
   MtrxMinusMatrx(arrVectSostGSK,  marrVectSost,3, 1, arrVectSost_KGSK);

	 	// ������������� SINS
		 mSins = TSins ( Environment
	,    MaxSig_Q,     MaxSig_Psi,     MaxSig_Tet
	,    MaxSig_dQdt,     MaxSig_dPsidt,    MaxSig_dTetdt
	,    MaxSig_H,    MaxSig_VH,    K1,    SigV
	,     mTVess     // �����
	,    mVVess   // ��� ��������
	,     mQ  // �������� ���� �����
	,     mPsi  // �������� ���� ������� �����
	,     mTet  // �������� ���� �������� �����
	,     mVQ  //�������� �������� ��������� ���� �����
	,     mVPsi   //�������� �������� ��������� ���� ������� �����
	,     mVTet   //�������� �������� ��������� ���� �������� ������
	,    marrVectSost[2]
	,    marrVectSost[5]
	,    mT_Q
	,     mT_Psi
	,    mT_Tet
	,   marrDelt
	,pwcharrFoldReport) ;

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



	const double VAlDelEps  = 0 ;
	const double VAlDelBet = 0 ;
	const double VAlEstEps =  vRealEps + VAlDelEps  ;
	const double VAlEstBet =  vRealBet + VAlDelBet ;
		 // ����� ����������� 2
 mDriver = TDriverMech (      DriverSigBet // �������� ��������� ���� Bet �������
								,    DriverSigEps // �������� ��������� ���� Eps  ������� (���� �����)
								,    DriverDynamicSigBet // �������� ��������� ���� �����  �������
								,    DriverDynamicSigEps // ��������  ������� ��������� ���� �����
								,InitTargData.mT  // ������� ����� �������� ����������
								,    VAlEstEps // ������(���������) ���� Eps
								,    VAlEstBet  // ������(���������) ���� Bet
								,     vRealEps // �������� ���� Eps
								,    vRealBet  // �������� ���� Bet
								,    VAlDelEps // ������ ��  Eps
								,    VAlDelBet  // ������ �� Bet
								,   pwcharrFoldReport) ;


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
	if((mparrBuff  = (double*)malloc( QUANT_COLS_BUFF_VESSEL * mLenMemoryAlloc * sizeof(double))) != NULL)
	{

	memset( mparrBuff,  0,  QUANT_COLS_BUFF_VESSEL * mLenMemoryAlloc * sizeof(double));
	}
	else
	{
	ShowMessage(L"Not memory for mparrBuff ") ;
	Abort() ;
	}

	}

	mArtComplex = TArtComplex();
	memset(marrArtParral ,0, 3 * sizeof(double)) ;//  ������ ���������� ��

}




	 // ����� ����������� 6
 TVessel::TVessel  ( const TShellBody ShellBody
			  ,const TFar_2D  Far_2D
			  ,const TTransmitAnt TransmitAnt
			  ,const double DriverSigBet // �������� ��������� ���� Bet �������
			,const double DriverSigEps // �������� ��������� ���� Eps  ������� (���� �����)
			,const double DriverDynamicSigBet // �������� ��������� ���� �����  �������
			,const double DriverDynamicSigEps // ��������  ������� ��������� ���� �����
			,const double MaxSig_Q, const double MaxSig_Psi, const double MaxSig_Tet      // ����
			,const double MaxSig_dQdt, const double MaxSig_dPsidt,const double MaxSig_dTetdt   // ����
			,const double MaxSig_H,const double MaxSig_VH,const double K1,const double SigV   // ����
			,const TEnvironment Environment
			,const double Width,const double Length, double *arrPar,const  double MaxQ ,const  double T_Q
			,const double MaxPsi,const double T_Psi ,const  double MaxTet
			,const double T_Tet,const double MaxVert, const double Q0,const double VVess
			,const TInitTargData InitTargData
			,const double MaxAmp_AftFlexure,const double T_AftFlexure,const double MaxAmp_BoardFlexure
			,const double T_BoardFlexure, const TControlSyst ControlSyst, double *arrArtPar,TArtComplex ArtComplex
			, wchar_t *pwcharrFoldReport)
{
  mControlSyst = ControlSyst;
	mShellBody = ShellBody;
	 // ������������� ����������
	mFar_2D = Far_2D ;

	mTransmitAnt = TransmitAnt;

	memcpy(marrArtParral, arrArtPar , 3 * sizeof(double));
	mArtComplex = ArtComplex;

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

	mAmp_AftFlexure = MaxAmp_AftFlexure * valCoeffEnv ;
	mT_AftFlexure = T_AftFlexure;
	mAmp_BoardFlexure = MaxAmp_BoardFlexure * valCoeffEnv ;
	mT_BoardFlexure = T_BoardFlexure;
	mPhase0_AftFlexure = KKK * getRand01() * 2. * M_PI;
	mPhase0_BoardFlexure = KKK * getRand01() * 2. * M_PI;


	mQ0 =Q0 ;
	mVVess = VVess ;
	mTVess  = InitTargData.mT ;
	for (int i = 0 ; i < 4; i++)
	{
	 marrDelt [i] =KKK * getRand01() * 2. * M_PI;
	}

	memset(marrVectSost,0, 9 * sizeof(double)) ;
	mQ = mQ0 +  mMaxQ * cos(marrDelt[0]);
	mVQ =  2. * M_PI/mT_Q * mMaxQ * sin(marrDelt[0]);

	mPsi = mMaxPsi * cos(marrDelt[1]);
	mVPsi =  2. * M_PI/mT_Psi * mMaxPsi * sin( marrDelt[1]);

	mTet =  mMaxTet * cos( marrDelt[2]);
	mVTet =  -2. * M_PI/mT_Tet * mMaxTet * sin( marrDelt[2]);
	// ������ ��������� ������� � ��� ��������
	 marrVectSost[0] = 0;
	 marrVectSost[1] = 0;
	 marrVectSost[2] = mMaxVert * sin( marrDelt[3]) ;
	 marrVectSost[3] = mVVess * sin(mQ) ;
	 marrVectSost[4] = mVVess * cos(mQ) ;
	 marrVectSost[5] = 2. * M_PI/ mT_Psi * mMaxVert * cos( marrDelt[3]) ;

	 ///
	 //  2.����� ������� ��������� ���� � ���� �� ������ InitTargData.mT

   	 if (fabs(InitTargData.mH / InitTargData.mR) > 1.)
{
	ShowMessage(L"ERROR");
	int iii = 0;


}


	const double valeps = arcSin(InitTargData.mH / InitTargData.mR) ;

   double arrVectSostGSK[6] ={0.}, arrVectSost_KGSK [6] = {0} ;
   arrVectSostGSK[0] = InitTargData.mR * cos(valeps) * sin(InitTargData.mBearing) ;
   arrVectSostGSK[1] = InitTargData.mR * cos(valeps) * cos(InitTargData.mBearing) ;
   arrVectSostGSK[2] = InitTargData.mH ;
   arrVectSostGSK[3] = InitTargData.mV * sin(InitTargData.mTargZenitAng) * sin(InitTargData.mTargCourse) ;
   arrVectSostGSK[4] = InitTargData.mV * sin(InitTargData.mTargZenitAng) * cos(InitTargData.mTargCourse) ;
   arrVectSostGSK[5] = InitTargData.mV * cos(InitTargData.mTargZenitAng) ;
   MtrxMinusMatrx(arrVectSostGSK,  marrVectSost,3, 1, arrVectSost_KGSK);

		// ������������� SINS
		 mSins = TSins ( Environment
	,    MaxSig_Q,     MaxSig_Psi,     MaxSig_Tet
	,    MaxSig_dQdt,     MaxSig_dPsidt,    MaxSig_dTetdt
	,    MaxSig_H,    MaxSig_VH,    K1,    SigV
	,     mTVess     // �����
	,    mVVess   // ��� ��������
	,     mQ  // �������� ���� �����
	,     mPsi  // �������� ���� ������� �����
	,     mTet  // �������� ���� �������� �����
	,     mVQ  //�������� �������� ��������� ���� �����
	,     mVPsi   //�������� �������� ��������� ���� ������� �����
	,     mVTet   //�������� �������� ��������� ���� �������� ������
	,    marrVectSost[2]
	,    marrVectSost[5]
	,    mT_Q
	,     mT_Psi
	,    mT_Tet
	,   marrDelt
	,pwcharrFoldReport) ;

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



	const double VAlDelEps  = 0 ;
	const double VAlDelBet = 0 ;
	const double VAlEstEps =  vRealEps + VAlDelEps  ;
	const double VAlEstBet =  vRealBet + VAlDelBet ;
		 // ����� ����������� 2
 mDriver = TDriverMech (      DriverSigBet // �������� ��������� ���� Bet �������
								,    DriverSigEps // �������� ��������� ���� Eps  ������� (���� �����)
								,    DriverDynamicSigBet // �������� ��������� ���� �����  �������
								,    DriverDynamicSigEps // ��������  ������� ��������� ���� �����
								,InitTargData.mT  // ������� ����� �������� ����������
								,    VAlEstEps // ������(���������) ���� Eps
								,    VAlEstBet  // ������(���������) ���� Bet
								,     vRealEps // �������� ���� Eps
								,    vRealBet  // �������� ���� Bet
								,    VAlDelEps // ������ ��  Eps
								,    VAlDelBet  // ������ �� Bet
								,   pwcharrFoldReport) ;


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
	if((mparrBuff  = (double*)malloc( QUANT_COLS_BUFF_VESSEL * mLenMemoryAlloc * sizeof(double))) != NULL)
	{

	memset( mparrBuff,  0,  QUANT_COLS_BUFF_VESSEL * mLenMemoryAlloc * sizeof(double));
	}
	else
	{
	ShowMessage(L"Not memory for mparrBuff ") ;
	Abort() ;
	}

	}
}



	 // ����� ����������� 7
 TVessel::TVessel  (
		const TFar_2D  Far_2D
		,const TTransmitAnt TransmitAnt
		,TDriverMech Driver
		,TSins Sins
		,const TEnvironment Environment
		,const double Width,const double Length, double *arrPar
		,const  double MaxQ ,const  double T_Q
		,const double MaxPsi,const double T_Psi ,const  double MaxTet
		,const double T_Tet,const double MaxVert, const double Q0,const double VVess
		,const double MaxAmp_AftFlexure,const double T_AftFlexure,const double MaxAmp_BoardFlexure
		,const double T_BoardFlexure, const TControlSyst ControlSyst, double *arrArtPar
		, wchar_t *pwcharrFoldReport)
{
  mControlSyst = ControlSyst;

	 // ������������� ����������
	mFar_2D = Far_2D ;

	mTransmitAnt = TransmitAnt;

	memcpy(marrArtParral, arrArtPar , 3 * sizeof(double));


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

	mAmp_AftFlexure = MaxAmp_AftFlexure * valCoeffEnv ;
	mT_AftFlexure = T_AftFlexure;
	mAmp_BoardFlexure = MaxAmp_BoardFlexure * valCoeffEnv ;
	mT_BoardFlexure = T_BoardFlexure;
	mPhase0_AftFlexure = KKK * getRand01() * 2. * M_PI;
	mPhase0_BoardFlexure =  KKK *getRand01() * 2. * M_PI;


	mQ0 =Q0 ;
	mVVess = VVess ;

	for (int i = 0 ; i < 4; i++)
	{
	 marrDelt [i] =KKK * getRand01() * 2. * M_PI;
	}

	memset(marrVectSost,0, 9 * sizeof(double)) ;
	mQ = mQ0 +  mMaxQ * cos(marrDelt[0]);
	mVQ =  2. * M_PI/mT_Q * mMaxQ * sin(marrDelt[0]);

	mPsi = mMaxPsi * cos(marrDelt[1]);
	mVPsi =  2. * M_PI/mT_Psi * mMaxPsi * sin( marrDelt[1]);

	mTet =  mMaxTet * cos( marrDelt[2]);
	mVTet =  -2. * M_PI/mT_Tet * mMaxTet * sin( marrDelt[2]);
	// ������ ��������� ������� � ��� ��������
	 marrVectSost[0] = 0;
	 marrVectSost[1] = 0;
	 marrVectSost[2] = mMaxVert * sin( marrDelt[3]) ;
	 marrVectSost[3] = mVVess * sin(mQ) ;
	 marrVectSost[4] = mVVess * cos(mQ) ;
	 marrVectSost[5] = 2. * M_PI/ mT_Psi * mMaxVert * cos( marrDelt[3]) ;


	mSins = Sins;


		 // ����� ����������� 2
 mDriver = Driver;



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
	if((mparrBuff  = (double*)malloc( QUANT_COLS_BUFF_VESSEL * mLenMemoryAlloc * sizeof(double))) != NULL)
	{

	memset( mparrBuff,  0,  QUANT_COLS_BUFF_VESSEL * mLenMemoryAlloc * sizeof(double));
	}
	else
	{
	ShowMessage(L"Not memory for mparrBuff ") ;
	Abort() ;
	}

	}
}

// ������������� ������� �� ���������� ���� � ��
void TVessel::initVess  ( const TShellBody ShellBody
,const TInitTargData InitTargData
,const TControlSyst ControlSyst, TArtComplex ArtComplex
)
{

	mShellBody = ShellBody;
	mArtComplex = ArtComplex;

	mPhase0_AftFlexure =  getRand01() * 2. * M_PI;
	mPhase0_BoardFlexure =  getRand01() * 2. * M_PI;
	mControlSyst = ControlSyst;


	mTVess  = InitTargData.mT ;
	for (int i = 0 ; i < 4; i++)
	{
	 marrDelt [i] = getRand01() * 2. * M_PI;
	}

	// ��� ������� ����� ��������� �� ������� �� ������������ � �������� � ������������ ��������
	marrDelt[0] = 0.26078;
	marrDelt[1] = 1.0377;
	marrDelt[2] = 5.1238;
	marrDelt[3] = 4.30735;
	///

	memset(marrVectSost,0, 9 * sizeof(double)) ;
	mQ = mQ0 +  mMaxQ * cos(marrDelt[0]);
	mVQ =  2. * M_PI/mT_Q * mMaxQ * sin(marrDelt[0]);

	mPsi = mMaxPsi * cos(marrDelt[1]);
	mVPsi =  2. * M_PI/mT_Psi * mMaxPsi * sin( marrDelt[1]);

	mTet =  mMaxTet * cos( marrDelt[2]);
	mVTet =  -2. * M_PI/mT_Tet * mMaxTet * sin( marrDelt[2]);
	// ������ ��������� ������� � ��� ��������
	 marrVectSost[0] = 0;
	 marrVectSost[1] = 0;
	 marrVectSost[2] = mMaxVert * sin( marrDelt[3]) ;
	 marrVectSost[3] = mVVess * sin(mQ) ;
	 marrVectSost[4] = mVVess * cos(mQ) ;
	 marrVectSost[5] = 2. * M_PI/ mT_Psi * mMaxVert * cos( marrDelt[3]) ;

	 ///
	 //  2.����� ������� ��������� ���� � ���� �� ������ InitTargData.mT

   	 if (fabs(InitTargData.mH / InitTargData.mR) > 1.)
{
	ShowMessage(L"ERROR");
	int iii = 0;


}


	const double valeps = arcSin(InitTargData.mH / InitTargData.mR) ;

   double arrVectSostGSK[6] ={0.}, arrVectSost_KGSK [6] = {0} ;
   arrVectSostGSK[0] = InitTargData.mR * cos(valeps) * sin(InitTargData.mBearing) ;
   arrVectSostGSK[1] = InitTargData.mR * cos(valeps) * cos(InitTargData.mBearing) ;
   arrVectSostGSK[2] = InitTargData.mH ;
   arrVectSostGSK[3] = InitTargData.mV * sin(InitTargData.mTargZenitAng) * sin(InitTargData.mTargCourse) ;
   arrVectSostGSK[4] = InitTargData.mV * sin(InitTargData.mTargZenitAng) * cos(InitTargData.mTargCourse) ;
   arrVectSostGSK[5] = InitTargData.mV * cos(InitTargData.mTargZenitAng) ;
   MtrxMinusMatrx(arrVectSostGSK,  marrVectSost,3, 1, arrVectSost_KGSK);

		// ������������� SINS
		 mSins.fillValues_Delts_and_Ests(InitTargData.mT
		 , mVVess
		 ,	mQ
		 ,mPsi
		 ,	mTet
		 ,mVQ
	,mVPsi
	, mVTet
	,marrVectSost[2]
	,marrVectSost[5]
	,mT_Q
	,mT_Psi
	,mT_Tet
	, marrDelt) ;




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



	const double VAlDelEps  = 0 ;
	const double VAlDelBet = 0 ;
	const double VAlEstEps =  vRealEps + VAlDelEps  ;
	const double VAlEstBet =  vRealBet + VAlDelBet ;

		mDriver.init (InitTargData.mT  // ������� ����� �������� ����������
								,    VAlEstEps // ������(���������) ���� Eps
								,    VAlEstBet  // ������(���������) ���� Bet
								,     vRealEps // �������� ���� Eps
								,    vRealBet  // �������� ���� Bet
								,    VAlDelEps // ������ ��  Eps
								,    VAlDelBet  // ������ �� Bet
					) ;


	 // ������ ��������� ������� � ��� - ���������

	 marrEstVectSost[0] = 0;
	 marrEstVectSost[1] = 0;
	 marrEstVectSost[2] = mSins.mEstH ;
	 marrEstVectSost[3] = mSins.mEstVVess * sin( mSins.mEstQ) ;
	 marrEstVectSost[4] = mSins.mEstVVess * cos( mSins.mEstQ) ;
	 marrEstVectSost[5] = mSins.mEstVH ;



}


double TVessel::calcAmpAftFlexure(const double VAly)
{
	return VAly * mAmp_AftFlexure / 100.;
}

double TVessel::calcAmpBoardFlexure(const double VAly)
{
	return VAly * mAmp_BoardFlexure / 100.;
}


//���������� �������� ����� ������� ������ � ����� ��������� ���
// OUTPUT:
// arrMu [5] - ������ �����
//Q, Psi, Tet, Bet, Eps
//
 void TVessel::calcVectTrueDeckAngles_For_Far(double *arrMu)
 {
	 double  valVQ =0., valVPsi =0., valVTet = 0.;
	calcDeckAngles(mTVess, marrParral,&arrMu[0], &valVQ,&arrMu[1], &valVPsi,&arrMu[2] ,&valVTet   );
	arrMu[3] = mDriver.mRealBet ;
	arrMu[4] = mDriver.mRealEps ;
 }

// �������� ��������� �� ������������� ���� � ����������� C�
//  arrInp - ���������� �������� �����
// valR - ����������
// valBet - ���� ����������� �� ������� ������� �� ��� OY �� ������������� �� �������� ����� �� ��������� OXY
// valEps - ���� ����� ������-�������� ����� � �������������� ����������
void TVessel::recalcCoord_INTO_Spherical(double *arrInp, double &valR, double &valBet, double &valEps)
{
	valR = sqrt ( arrInp[0] * arrInp[0] +arrInp[1] * arrInp[1] +arrInp[2] * arrInp[2] ) ;

	if(fabs (arrInp[2] / valR) > 1.)
	{
		ShowMessage(L"ERROR");
		int iii = 0;
		return;
	}

	valEps  = arcSin(arrInp[2] / valR ) ;



	if(fabs (arrInp[0]/ sqrt (  arrInp[0] * arrInp[0] +arrInp[1] * arrInp[1])) > 1.)
	{
		ShowMessage(L"ERROR");
		int iii = 0;
		return;
	}
	valBet  = arcSin(arrInp[0]/ sqrt (  arrInp[0] * arrInp[0] +arrInp[1] * arrInp[1])) ;
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

 void TVessel::recalcVess(const double valT, const double VAlTargDesEps, const double VAlTargDesBet)
 {
	 const double h = valT - mTVess ;
	 if (h < 0)
	 {
		ShowMessage(L"������ � �������� �������� �������") ;
		return ;
	 }
	 calcCentreDeckAngles(valT) ;
	// �������� ������� ��������� ������� � ��� ���������
	 marrVectSost[0] += marrVectSost[3]*h;
	 marrVectSost[1] += marrVectSost[4]*h;
	 marrVectSost[2] = mMaxVert * sin(2. * M_PI / mT_Psi*valT + marrDelt[3]) ;
	 marrVectSost[3] = mVVess * sin(mQ) ;
	 marrVectSost[4] = mVVess * cos(mQ) ;
	 marrVectSost[5] = 2. * M_PI/ mT_Psi * mMaxVert * cos(2. * M_PI/ mT_Psi*valT + marrDelt[3]) ;
	 mTVess = valT ;

	 ///

	 // �������� ������� ��������� - ������
	 marrEstVectSost[0] += marrEstVectSost[3] * h;
	 marrEstVectSost[1] += marrEstVectSost[4] * h;
	 marrEstVectSost[2] = mSins.mEstH ;
	 marrEstVectSost[3] = mSins.mEstVVess * sin( mSins.mEstQ) ;
	 marrEstVectSost[4] = mSins.mEstVVess * cos( mSins.mEstQ) ;
	 marrEstVectSost[5] = mSins.mEstVH ;

		 mSins.recalcSins_v0(valT
		,mVVess //VVess
		,mQ //	const double Q
		,mPsi// const double Psi
		,mTet// Tet
		,mVQ //  / VQ
		, mVPsi //   VPsi
		,mVTet//  VTet
		,marrVectSost[2] //    H
		,marrVectSost[5] // VH
		,mT_Q// c  T_Q
		,mT_Psi //   T_Psi
		,mT_Tet //   T_Tet
		,marrDelt// double *arrDelt
							   ) ;
		 // ������� Driver
		//  valTargDesEps,valTargDesBet - ���� ������������
		mDriver.recalcDriver(valT,VAlTargDesEps,VAlTargDesBet);
		updateReportData();  ///
 }

 // ���������� ����� ���������� ������ � ������ ������� � ������ valT
 void TVessel::calcCentreDeckAngles(const double valT)
 {
   mQ = mQ0 +  mMaxQ * cos(2. * M_PI/mT_Q * valT -  marrDelt[0]);
   mVQ =  -2. * M_PI/mT_Q * mMaxQ * sin(2. * M_PI/mT_Q * valT -  marrDelt[0]);

	 mPsi = mMaxPsi * cos(2. * M_PI/mT_Psi * valT -  marrDelt[1]); // �������
   mVPsi =  -2. * M_PI/mT_Psi * mMaxPsi * sin(2. * M_PI/mT_Psi * valT -  marrDelt[1]);

   mTet =  mMaxTet * cos(2. * M_PI/mT_Tet * valT -  marrDelt[2]); // ��������
   mVTet =  -2. * M_PI/mT_Tet * mMaxTet * sin(2. * M_PI/mT_Tet * valT -  marrDelt[2]);

 }

 // ���������� ����� ���������� ������ � ������ ������� � ������ valT
 void TVessel::calcDeckAngles(const double valT, double *arrPointPositionPSK
	,double *pvalQ, double *pvalVQ,double *pvalPsi,double *pvalVPsi
	,double *pvalTet ,double *pvalVTet)
 {
		calcCentreDeckAngles(valT);
	 *pvalQ = mQ ;
	 *pvalVQ =mVQ ;
	 //
	 double valTemp = mAmp_AftFlexure/ 100.* arrPointPositionPSK[1];
	 *pvalPsi = mPsi + valTemp * cos(2. * M_PI/mT_AftFlexure * valT -  mPhase0_AftFlexure);
	 *pvalVPsi = mVPsi - valTemp * 2.* M_PI/mT_AftFlexure * sin(2. * M_PI/mT_AftFlexure * valT -  mPhase0_AftFlexure);
		valTemp = mAmp_BoardFlexure/ 100.* arrPointPositionPSK[1];
	 *pvalTet = mTet + valTemp * cos(2. * M_PI/mT_BoardFlexure * valT -  mPhase0_BoardFlexure);
	 *pvalVTet = mVTet - valTemp * 2.* M_PI/mT_BoardFlexure * sin(2. * M_PI/mT_BoardFlexure * valT -  mPhase0_BoardFlexure);

 }


// ��������� ���������� � ������� ��� ������
 void TVessel::updateReportData()
 {
     if(mpwcharrFoldReport == NULL) return ;
	 if (mQuantPntReport ==  mLenMemoryAlloc)
	 {
	   mLenMemoryAlloc += 2000 ;
	   mparrBuff = (double *)realloc(mparrBuff,QUANT_COLS_BUFF_VESSEL * mLenMemoryAlloc * sizeof(double)) ;

	 }

	   int num0 =  mQuantPntReport * QUANT_COLS_BUFF_VESSEL;

	   mparrBuff [num0] =     mTVess ;

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
 void TVessel::WriteReport()
 {
	 if (mpwcharrFoldReport == NULL )return ;

	 wchar_t wcharrPath [300] = {0} ;
	 wcscpy(wcharrPath, mpwcharrFoldReport);
	 wcscat(wcharrPath, L"\\VesselReport");

		_wmkdir(wcharrPath);
		 wcscat(wcharrPath, L"\\");

	 wchar_t wcharrFileNames[QUANT_COLS_BUFF_VESSEL * 30] = {0};

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

	 double *pscaleY = new double  [QUANT_COLS_BUFF_VESSEL] ;
	 for (int i = 0; i < QUANT_COLS_BUFF_VESSEL; i++)  pscaleY[i] = 1.;




	 pscaleY[3]  = 100.;
	 pscaleY[6]  = 100.;
	 pscaleY[9]  = 10000;
	 pscaleY[10] = 10000.;
	 pscaleY[11] = 10000;
	 pscaleY[12] = 10000;
	 pscaleY[13] = 10000.;



	 for (int i = 1; i < QUANT_COLS_BUFF_VESSEL; i++)
	 {
	 TYrWriteShapeFile::WriteOneReport(wcharrPath  // ���� � �����
								  ,mparrBuff // ������ � ����������� - ������� nBuffRows x nBuffCols
								  ,QUANT_COLS_BUFF_VESSEL // - �-�� ���������� � ������� ��������� ���������� � ������
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
								  ,QUANT_COLS_BUFF_VESSEL // - �-�� ���������� � ������� ��������� ���������� � ������
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
	delete []pscaleY ;
 }


 // ����������  ������  �������������
 void TVessel::WriteReport(wchar_t *pwcharrPath)
 {
	 if (pwcharrPath == NULL )return ;
	 wchar_t wcharrPath [300] = {0} ;
	 wcscpy(wcharrPath, pwcharrPath);
	 wcscat(wcharrPath, L"\\");


	 wchar_t wcharrFileNames[QUANT_COLS_BUFF_VESSEL * 30] = {0};

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

	 double *pscaleY = new double  [QUANT_COLS_BUFF_VESSEL] ;
	 for (int i = 0; i < QUANT_COLS_BUFF_VESSEL; i++)  pscaleY[i] = 1.;




	 pscaleY[3]  = 100.;
	 pscaleY[6]  = 100.;
	 pscaleY[9]  = 10000;
	 pscaleY[10] = 10000.;
	 pscaleY[11] = 10000;
	 pscaleY[12] = 10000;
	 pscaleY[13] = 10000.;



	 for (int i = 1; i < QUANT_COLS_BUFF_VESSEL; i++)
	 {
	 TYrWriteShapeFile::WriteOneReport(wcharrPath  // ���� � �����
								  ,mparrBuff // ������ � ����������� - ������� nBuffRows x nBuffCols
								  ,QUANT_COLS_BUFF_VESSEL // - �-�� ���������� � ������� ��������� ���������� � ������
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
								  ,QUANT_COLS_BUFF_VESSEL // - �-�� ���������� � ������� ��������� ���������� � ������
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
	delete []pscaleY ;
 }

 // ������������� ������� �� �����  valT
 // SINS, Driver , Measurer �� ���������������!!!
  void TVessel::Move(const double valT,const double valStep)
 {
	  int iCirc = (valT - mTVess)/ valStep -1;
  // double valTTemp = mTVess + ((double)iCirc) *valStep  ;

   for (int i = 0; i < iCirc; i++)
   {
	 recalcVess(mTVess + valStep,0.,0.) ;

   }

   recalcVess(valT, 0.,0. ) ;


 }

// ������������� ������� ��������� ������� ������ TVessel �� �����  valTExtr
// INPUT:
// valTExtr - ������ �������, �� ��������� ��������� ���������� �������������
// OUTPUT:
// arrVSVessExtrap - ������������������ ������ ���������
void TVessel::VSProlong(const double valTExtr, double *arrVSVessExtrap)
{

  // ������������� ������� ��������� ������� �� ������  valT

  double h = valTExtr - mTVess ;
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
// INPUT:
//InpZamer - ������� ����� � ���
// InpZamer.marrMeas [0] - ���� V
// InpZamer.marrMeas [1] - ��������� R
// InpZamer.marrMeas[2] - ���� U
//  OUTPUT:
// arrZam_GSK - ������ ������ � ���� - 3-� ������ ������
// arrZam_GSK - ������ ������� ������, 9-�� ������ ������
//
void  TVessel::GetZamer_IN_KGSK (TZamer InpASKZamer, TZamer *pOutKGSKZamer )

{
	const double valRZv  = InpASKZamer.marrMeas [1] ;
	const double  valVZv = InpASKZamer.marrMeas [0] ;
	const double  valUZv = InpASKZamer.marrMeas[2] ;
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
	MtrxMultMatrx(matrPereh_PSK_V_KGSK,3,3, arrT1, 1, (*pOutKGSKZamer ).marrMeas) ;

	///

	calcSummarizedCorMtrx_ErrMes_In_GSK (InpASKZamer,(*pOutKGSKZamer ).marrCorr) ;
	(*pOutKGSKZamer ).mT= InpASKZamer.mT ;

}

// �������� ������   �� ��������� ����� � ��������� ���  � �������� ��������
//  OUTPUT:
// arrZam_GSK - ������ ������ � ��� - 3-� ������ ������
// arrKZam_GSK - ������ ������� ������, 9-�� ������ ������
//
void  TVessel::GetZamer_IN_GSK (TZamer InpASKZamer,  TZamer *pOutGSKZamer )

{
	GetZamer_IN_KGSK (InpASKZamer, pOutGSKZamer);
  for (int i = 0; i < 3; i++) (*pOutGSKZamer).marrMeas[i] +=  marrEstVectSost [i] ;



}

// ������ �������������� �������������� ������� ������ ���������� ������, ����������� ������
 // �������������� ������ ,  ������ ������� � ������ ����������� ����� � ���  � ������� ������� �������
// ��� ������ ����������� ��� ������
// INPUT:
//InpZamer - ������� ����� � ���
// InpZamer.marrMeas [0] - ���� V
// InpZamer.marrMeas [1] - ��������� R
// InpZamer.marrMeas[2] - ���� U
 void  TVessel::calcSummarizedCorMtrx_ErrMes_In_GSK (TZamer InpASKZamer,double *arrCorrMtrx_GSK )
{
   double ar_dF_po_dBet_sq [9],ar_dF_po_dEps_sq [9], matrPereh_ASK_V_KGSK[9],arrMu[5] ;
	   double  arrDT0[9] = {0}, arrDT1[9] = {0}, arrDT2[9] = {0}
				,arrDT3[9] = {0},arrDT4[9] = {0},arrDT5[9] = {0},arrDT51[9] = {0};
	arrMu[0] = mSins.mEstQ;
	arrMu[1] = mSins.mEstPsi ;
	arrMu[2] = mSins.mEstTet ;
	arrMu[3] = mDriver.mEstBet ;
	arrMu[4] = mDriver.mEstEps ;

	double valDispV = InpASKZamer.marrCorr[0];
	double valDispR = InpASKZamer.marrCorr[4];
	double valDispU = InpASKZamer.marrCorr[8];
		 calc_dF_po_dBet_sq(mDriver.mEstBet,  mDriver.mEstEps,ar_dF_po_dBet_sq ) ;
	   calc_dF_po_dEps_sq(mDriver.mEstBet,  mDriver.mEstEps,ar_dF_po_dEps_sq ) ;
	   calcMatr_ASK_v_KGSK(arrMu,matrPereh_ASK_V_KGSK) ;
		 arrDT0[0] = valDispV * InpASKZamer.marrMeas[1]* InpASKZamer.marrMeas[1] ;
		 arrDT0[4] = valDispR ;
		 arrDT0[8] = valDispU * InpASKZamer.marrMeas[1]* InpASKZamer.marrMeas[1] ;
	   MtrxMultMatrx(matrPereh_ASK_V_KGSK,3, 3, arrDT0,3, arrDT1) ;
	   MtrxMultMatrxTransp(arrDT1,3, 3, matrPereh_ASK_V_KGSK, 3, arrDT2) ;

		 double temp0 =  InpASKZamer.marrMeas[1]* InpASKZamer.marrMeas[1]
		 * mDriver.mSigBet* mDriver.mSigBet ;
	   MatrxMultScalar(ar_dF_po_dBet_sq , 3, 3, temp0,arrDT3);
	   MtrxSumMatrx(arrDT2, arrDT3,3,3, arrDT4) ;

		 double temp1 =  InpASKZamer.marrMeas[1]* InpASKZamer.marrMeas[1]
	   * mDriver.mSigEps* mDriver.mSigEps ;
	   MatrxMultScalar(ar_dF_po_dEps_sq , 3, 3, temp1,arrDT5);
	   MtrxSumMatrx(arrDT4, arrDT5,3,3, arrDT51) ;

	   // ������ ������, ��������� �������� ����������� �����
	   double ar_dF_po_dTet_sq [9] = {0},ar_dF_po_dQ_sq [9] = {0},ar_dF_po_dPsi_sq [9] = {0}, arrDT6[9] = {0}
		  , arrDT7[9]  = {0}, arrDT8[9] = {0}, arrDT9[9] = {0}, arrDT10[9]  = {0};
	   calc_dF_po_dQ_sq  (mDriver.mEstBet,  mDriver.mEstEps, ar_dF_po_dQ_sq) ;
	   calc_dF_po_dPsi_sq(mDriver.mEstBet,  mDriver.mEstEps, ar_dF_po_dPsi_sq) ;
	   calc_dF_po_dTet_sq(mDriver.mEstBet,  mDriver.mEstEps, ar_dF_po_dTet_sq) ;

		 double temp2 = InpASKZamer.marrMeas[1]* InpASKZamer.marrMeas[1]
		 * mSins.mSig_Q * mSins.mSig_Q  ;
		 double temp5 = calcAmpAftFlexure(marrParral[1]);
	   double temp3 = InpASKZamer.marrMeas[1]* InpASKZamer.marrMeas[1]
		 *( mSins.mSig_Psi* mSins.mSig_Psi + temp5 * temp5/2.);
		 double temp6 = calcAmpBoardFlexure(marrParral[1]);
	   double temp4 = InpASKZamer.marrMeas[1]* InpASKZamer.marrMeas[1]
		 * (mSins.mSig_Tet* mSins.mSig_Tet + temp6 * temp6/2.);

	   MatrxMultScalar(ar_dF_po_dQ_sq   , 3, 3, temp2,arrDT6);
	   MatrxMultScalar(ar_dF_po_dPsi_sq , 3, 3, temp3,arrDT7);
	   MatrxMultScalar(ar_dF_po_dTet_sq , 3, 3, temp4,arrDT8);

	   MtrxSumMatrx(arrDT51, arrDT6,3,3, arrDT9) ;
	   MtrxSumMatrx(arrDT9, arrDT7,3,3, arrDT10) ;
	   MtrxSumMatrx(arrDT10, arrDT8,3,3, arrCorrMtrx_GSK) ;
 }


 // �������� ������ � ��� c ������ ������ ��������
void TVessel::createMeasure(const double VAlCoefAntp,const double VAlTrue_R,const double VAlTrue_V, const double VAlTrue_U
	, const double VAlT, const double VAlTargDesEps, const double VAlTargDesBet
	 , TEtalonSign EtalonSign , const double VAlTargEPR
	, const double VAlPowerPrd, const double VAlKYPrd, TZamer *pOutASKZamer )
{
	// ���������� ������ ������� ��� ������� ������������ � �������� ������ �������
	double arrFarKGSK[3] ={0.};
	RecalcVect_PSK_CT_True_INTO_KGSK_True (marrParral
	,VAlTargDesEps,  VAlTargDesBet,arrFarKGSK,3) ;

	///
	double arrMeas[3] = {0.};
	arrMeas[0] = VAlTrue_V;
	arrMeas[1] = VAlTrue_R ;
	arrMeas[2] = VAlTrue_U;
	double arrCorr [9] = {0.};
	TZamer InpASKZamer( arrMeas,arrCorr, 0.) ;
	TZamer OutKGSKZamer;
	GetZamer_IN_KGSK( InpASKZamer, &OutKGSKZamer);
	 ///

	 double valAntpPhaze = 0.;
	 double valSigmaEps = -1., valSigmaBet = -1.0;
		mFar_2D.calc_SKZ_LAT(VAlCoefAntp, VAlTrue_R ,OutKGSKZamer.marrMeas[2]
		,arrFarKGSK[2],  VAlTargDesEps, VAlTargEPR
		,EtalonSign ,VAlPowerPrd,  VAlKYPrd,  valAntpPhaze
		,&valSigmaEps ,&valSigmaBet);

	  //valSigmaEps = 0.0015;///!!!!! ����� ������ !!!!!
	 (*pOutASKZamer).marrMeas[0] = VAlTrue_V + KKK* getGauss(0., valSigmaBet );
	 (*pOutASKZamer).marrMeas[1] = VAlTrue_R + KKK* getGauss(0., 10. );
	 (*pOutASKZamer).marrMeas[2] = VAlTrue_U + KKK* getGauss(0., valSigmaEps );
	 memset((*pOutASKZamer).marrCorr, 0, 9 * sizeof(double));
	 (*pOutASKZamer).marrCorr[0] = valSigmaBet * valSigmaBet;
	 (*pOutASKZamer).marrCorr[4] = mFar_2D.mDistSKZ * mFar_2D.mDistSKZ ;
	 (*pOutASKZamer).marrCorr[8] = valSigmaEps * valSigmaEps;
	 (*pOutASKZamer).mT = VAlT;
	 ///

//	mFar_2D.createMeasure_Targ_Plus_Measure_Antp(const double VAlTrue_R,const double VAlTrue_V, const double VAlTrue_U
 //	, const double VAlT,
}

/*
// �������� ������ � ��� c ������ ������ ��������
void TVessel::createMeasure(const double VAlTrue_R,const double VAlTrue_V, const double VAlTrue_U
	, const double VAlT, const double VAlTargDesEps, const double VAlTargDesBet
	 , TEtalonSign EtalonSign , const double VAlTargEPR
	, const double VAlPowerPrd, const double VAlKYPrd, TZamer *pOutASKZamer )
{
	// ���������� ������ ������� ��� ������� ������������ � �������� ������ �������
	double arrFarKGSK[3] ={0.};
	RecalcVect_PSK_CT_True_INTO_KGSK_True (marrParral
	,VAlTargDesEps,  VAlTargDesBet,arrFarKGSK,3) ;

	///
	double arrMeas[3] = {0.};
	arrMeas[0] = VAlTrue_V;
	arrMeas[1] = VAlTrue_R ;
	arrMeas[2] = VAlTrue_U;
	double arrCorr [9] = {0.};
	TZamer InpASKZamer( arrMeas,arrCorr, 0.) ;
	TZamer OutKGSKZamer;
	GetZamer_IN_KGSK( InpASKZamer, &OutKGSKZamer);
	 ///

	 double valAntpPhaze = 0.;
	 double valSigmaEps = -1., valSigmaBet = -1.0;
		mFar_2D.calc_SKZ_LAT(VAlTrue_R ,OutKGSKZamer.marrMeas[2]
		,arrFarKGSK[2],  VAlTargDesEps, VAlTargEPR
		,EtalonSign ,VAlPowerPrd,  VAlKYPrd,  valAntpPhaze
		,&valSigmaEps ,&valSigmaBet);


	 (*pOutASKZamer).marrMeas[0] = VAlTrue_V + getGauss(0., valSigmaBet );
	 (*pOutASKZamer).marrMeas[1] = VAlTrue_R + getGauss(0., 10. );
	 (*pOutASKZamer).marrMeas[2] = VAlTrue_U + getGauss(0., valSigmaEps );
	 memset((*pOutASKZamer).marrCorr, 0, 9 * sizeof(double));
	 (*pOutASKZamer).marrCorr[0] = valSigmaBet * valSigmaBet;
	 (*pOutASKZamer).marrCorr[4] = 100.;
	 (*pOutASKZamer).marrCorr[8] = valSigmaEps * valSigmaEps;
	 (*pOutASKZamer).mT = VAlT;
	 ///

//	mFar_2D.createMeasure_Targ_Plus_Measure_Antp(const double VAlTrue_R,const double VAlTrue_V, const double VAlTrue_U
 //	, const double VAlT,
}
*/


void TVessel::createMeasure_ForSingleTarg(const double VAlTrue_R,const double VAlTrue_V, const double VAlTrue_U
	, const double VAlT, const double VAlTargDesEps, const double VAlTargDesBet
	 , TEtalonSign EtalonSign , const double VAlTargEPR
	, const double VAlPowerPrd, const double VAlKYPrd, TZamer *pOutASKZamer )
{
	// ���������� ������ ������� ��� ������� ������������ � �������� ������ �������
	double arrFarKGSK[3] ={0.};
	RecalcVect_PSK_CT_True_INTO_KGSK_True (marrParral
	,VAlTargDesEps,  VAlTargDesBet,arrFarKGSK,3) ;

	double arrMeas[3] = {0.};
	arrMeas[0] = VAlTrue_V;
	arrMeas[1] = VAlTrue_R ;
	arrMeas[2] = VAlTrue_U;
	double arrCorr [9] = {0.};
	TZamer InpASKZamer( arrMeas,arrCorr, 0.) ;
	TZamer OutKGSKZamer;
	GetZamer_IN_KGSK( InpASKZamer, &OutKGSKZamer);
	 ///
	  double valSigmaEps = -1., valSigmaBet = -1.0;
	mFar_2D.calc_SKZ_SingleTarg(VAlTrue_R, VAlTargDesEps,  VAlTargEPR
  ,  EtalonSign ,VAlPowerPrd,  VAlKYPrd ,&valSigmaEps ,&valSigmaBet);

	 (*pOutASKZamer).marrMeas[0] = VAlTrue_V + getGauss(0., valSigmaBet );
	 (*pOutASKZamer).marrMeas[1] = VAlTrue_R + getGauss(0., 10. );
	 (*pOutASKZamer).marrMeas[2] = VAlTrue_U + getGauss(0., valSigmaEps );
	 memset((*pOutASKZamer).marrCorr, 0, 9 * sizeof(double));
	 (*pOutASKZamer).marrCorr[0] = valSigmaBet * valSigmaBet;
	 (*pOutASKZamer).marrCorr[4] = 100.;
	 (*pOutASKZamer).marrCorr[8] = valSigmaEps * valSigmaEps;
	 (*pOutASKZamer).mT = VAlT;
	 ///

//	mFar_2D.createMeasure_Targ_Plus_Measure_Antp(const double VAlTrue_R,const double VAlTrue_V, const double VAlTrue_U
 //	, const double VAlT,
}


// �������� �������  �� ��������  ���-�� � �������� ����
void  __fastcall TVessel::RecalcVect_PSK_CT_True_INTO_KGSK_True (double *arrPSK
, const double VAlTargDesEps, const double VAlTargDesBet,double *arrKGSK,int lenarrPSK )
{
	double arrMu[5] ={0} ;
	double  matrPereh_PSK_V_KGSK[9] = {0};
	arrMu[0] = mQ;
	arrMu[1] = mPsi ;
	arrMu[2] = mTet ;
	arrMu[3] = VAlTargDesBet ;
	arrMu[4] = VAlTargDesEps ;

	///
	// ������� � ����
	calcMatr_PSK_v_KGSK(arrMu, matrPereh_PSK_V_KGSK) ;
	MtrxMultMatrx(matrPereh_PSK_V_KGSK,3,3, arrPSK, 1, arrKGSK) ;

	///


	if (lenarrPSK < 6)return ;
//  �������� �� ��������


	MtrxMultMatrx(matrPereh_PSK_V_KGSK, 3, 3, &arrPSK[3], 1, & arrKGSK [3]) ;
	if (lenarrPSK < 9)return ;
//  �������� �� ���������
	MtrxMultMatrx(matrPereh_PSK_V_KGSK, 3, 3, &arrPSK[6], 1, & arrKGSK [6]) ;

	return ;
}


//����������� ���������� �� � ������� ������� � ����
// ��� ��������!!!!!
void  __fastcall TVessel::calcAY_Position(double *arr_AY_Position_KGSK)
{
	memset(arr_AY_Position_KGSK, 0, sizeof(double) * 3);
}

// �������� ������������� ������� ��������� �� �����  VAlTExtr �������
 void TVessel::extrapolateTrueVS_GSK(const double VAlTExtr, double *arrVessExtrapVS_GSK)
 {
	 memcpy(arrVessExtrapVS_GSK, marrVectSost, 9 * sizeof(double));
	 double arrT[3] = {0.};
	 MatrxMultScalar(&marrVectSost[3], 1, 3, VAlTExtr, arrT);
	 MtrxSumMatrx(arrT, marrVectSost,3, 1, arrVessExtrapVS_GSK) ;
 }


 //--------------------------------------------------------------
 // ���������� ��������� ������  �������� �������������� �������� ����� ����
 // INPUT:
 // VAlTExtrap - ����� �������������
 // OUTPUT:
 // pvalDispDEltaExtrIS_Psi, *pvalDispDEltaExtrIS_Tet , *pvalDispDEltaExtrIS_Q - ���������
 void TVessel::calcDispExtrapDeckAngles_SINS(const double VAlTExtrap, double *pvalDispDEltaExtrIS_Psi
	 , double *pvalDispDEltaExtrIS_Tet, double *pvalDispDEltaExtrIS_Q)
{
	*pvalDispDEltaExtrIS_Psi = mSins.mSig_Psi * mSins.mSig_Psi +  VAlTExtrap *VAlTExtrap * mSins.mSig_dPsidt * mSins.mSig_dPsidt
		+ VAlTExtrap * VAlTExtrap *VAlTExtrap * VAlTExtrap/4. * mMaxPsi* mMaxPsi * ( 2. * M_PI / mT_Psi)* ( 2. * M_PI / mT_Psi)
		* ( 2. * M_PI / mT_Psi)* ( 2. * M_PI / mT_Psi)/ 2.;

	*pvalDispDEltaExtrIS_Tet = mSins.mSig_Tet * mSins.mSig_Tet +  VAlTExtrap *VAlTExtrap * mSins.mSig_dTetdt * mSins.mSig_dTetdt
		+ VAlTExtrap * VAlTExtrap *VAlTExtrap * VAlTExtrap/4. * mMaxTet* mMaxTet * ( 2. * M_PI / mT_Tet)* ( 2. * M_PI / mT_Tet)
		* ( 2. * M_PI / mT_Tet)* ( 2. * M_PI / mT_Tet)/ 2.;

  *pvalDispDEltaExtrIS_Q = mSins.mSig_Q * mSins.mSig_Q +  VAlTExtrap * VAlTExtrap * mSins.mSig_dQdt * mSins.mSig_dQdt
		+ VAlTExtrap * VAlTExtrap *VAlTExtrap * VAlTExtrap/4. * mMaxQ* mMaxQ * ( 2. * M_PI / mT_Q)* ( 2. * M_PI / mT_Q)
		* ( 2. * M_PI / mT_Q)* ( 2. * M_PI / mT_Q)/ 2.;
}
 //--------------------------------------------------------------
 // ���������� ��������� ������    �������� �����, ��������� ����������� ������� �������

 // OUTPUT:
 // pvalDispDEltaExtrIS_Psi, *pvalDispDEltaExtrIS_Tet , *pvalDispDEltaExtrIS_Q - ���������
 void TVessel::calcDispExtrapDeckAngles_Deform( double *pvalDispDEltaExtrDeform_Psi
	 , double *pvalDispDEltaExtrDeform_Tet)
{
	double valaPsi = mAmp_AftFlexure /100. * marrArtParral[1];
	*pvalDispDEltaExtrDeform_Psi = valaPsi * valaPsi /2.;

	double valbTet = mAmp_BoardFlexure /100. * marrArtParral[1];
	*pvalDispDEltaExtrDeform_Tet = valbTet * valbTet /2.;
}

//-------------------------------------
double TVessel::arcSin(const double x)
{
double y = x;
if (fabs(y) > 0.99999999)
{
y = 0.99999999 * SIGNUM0(y);
}
return asin(y);
}
//-----------------------------------------
double TVessel::SIGNUM0(const double y)
{
 return (y >=0.)?1.:-1.;
}


	 #pragma package(smart_init)
