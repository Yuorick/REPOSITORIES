//---------------------------------------------------------------------------

#include <vcl.h>
#include <dir.h>
#pragma hdrstop

#include "FormAeroTargs.h"


#include "YrString.h"
#include <stdio.h>
#include <stdlib.h>
#include "AM_2D.h"
#include "Far_2D.h"
#include "Comp.h"
#include "Faceta.h"
#include "Far.h"
#include "Diagrams.h"
#include "Gauss.h"

#include "MatrixProccess.h"


#include "Target.h"
#include "YrWrite.h"

#include "Sins.h"
#include "DriverMech.h"

#include "Measurer.h"
#include "YrWriteShapeFile.h"

#include "InitTargData.h"
#include "CalcCorMatrx.h"

#include "Environment.h"
#include "Equations.h"
#include "MyShellTraj.h"

#include "SincDgr.h"
#include "URPointZ.h"
#include "ShellBody.h"
#include "ControlSyst.h"
#include "Line.h"
#include "UrPointXY.h"

#define MAX_QUANT_AM 1600
static int I_ENTER_COUNT = 0;
extern const double VAL_C;
#define ILenArr  1500 // ����������� ��������� ����� ���������
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TForm3 *Form3;
//---------------------------------------------------------------------------
__fastcall TForm3::TForm3(TComponent* Owner)
	: TForm(Owner)
{
	//ComboBox3->ItemIndex =2;
		mWSkz = 0.005;
		LabeledEdit27->Text = double (int (mWSkz * 1000. ))/ 1000.;
	   //	ComboBox1->ItemIndex =0;
	   //	ComboBox2->ItemIndex = 0;
		LabeledEdit19->Text = 10.;
		LabeledEdit20->Text = 5.8;

		create5P10();
		Button7->Visible = false;
		Edit3->Visible = false;
		mbCalculated = false;

}
//---------------------------------------------------------------------------
void __fastcall TForm3::Button1Click(TObject *Sender)
{
Close();
}
//---------------------------------------------------------------------------


void __fastcall TForm3::fncInputData()
{
		Button7->Visible = false;
		Edit3->Visible = false;
		mbCalculated = false;
		mAntCoeff = StrTo_Dbl_(LabeledEdit5->Text);
		// ����
		// ���� �������, ����

		mBearing0 = StrTo_Dbl_(LabeledEdit23->Text) * M_PI / 180.;
		// ������, �
		mElev0 = StrTo_Dbl_(LabeledEdit22->Text);

		// ���������, �
		mDist0 = StrTo_Dbl_(LabeledEdit24->Text);

		//
		mTargZenitAng0 = StrTo_Dbl_(LabeledEdit26->Text) * M_PI / 180. + M_PI/2.;
		// ���� �����
		mTargCourse0 = StrTo_Dbl_(LabeledEdit29->Text) * M_PI / 180.;

		///




		//mRateOfFire = StrTo_Dbl_(LabeledEdit49->Text) ;
		mAUDelayT = StrTo_Dbl_(LabeledEdit50->Text)/ 1000.;


		mSigDrivAY_U  =  StrTo_Dbl_(LabeledEdit52->Text) / 1000.;
		mSigDrivAY_dU_po_dt  =  StrTo_Dbl_(LabeledEdit53->Text) / 1000.;
		mWSkz = StrTo_Dbl_(LabeledEdit27->Text) ;

	 //	mFireBegin =  StrTo_Dbl_(LabeledEdit19->Text) * 1000. ;
	//}
 //	else
 //	{
		// ���������������
		// mWidthDgr =  StrTo_Dbl_(LabeledEdit3->Text) * M_PI/3000./2. ;
	//}

		switch(ComboBox3->ItemIndex)   // ����
	{
	case 0:
	mEnumTargType = GARPUN_V300;
	mVelocity0 = 300.;
	mTargEPR = 0.1;
   //	mWSkz = 0.002;
	break;
	case 1:
	mEnumTargType = GARPUN_V700;
	mVelocity0 = 700.;
	mTargEPR = 0.1;
   //	mWSkz = 0.01;
	break ;

	case 2:
	mEnumTargType = 	PLANE ;
	mVelocity0 = 300.;
	mTargEPR = 1.;
  //	mWSkz = 0.01;
	break;

	}
	LabeledEdit25->Text = int (mVelocity0 ) ;
	LabeledEdit28->Text = double (int (mTargEPR * 1000. ))/ 1000.;
 //	LabeledEdit27->Text = double (int (mWSkz * 1000. ))/ 1000.;
	 //	LabeledEdit22->Text = mElev0 ;
   // ��. ���������� ShellBody (mEnumShellType). ShellBody ���������� ����� ��������� �����������
	switch( ComboBox1->ItemIndex) // ��� ��
		{
			case 0:
			mEnumShellType = CALIBRO_130;
			menumCannonType = A_192M;
			LabeledEdit2->Text = 130;
			break;

			case 1:
			mEnumShellType = CALIBRO_100;
			menumCannonType = A_190_01;
			LabeledEdit2->Text = 100;
			break;

			case 2:
			 mEnumShellType = CALIBRO_30;
			 menumCannonType = AK_630;
			 ComboBox2->ItemIndex = 2;
			 LabeledEdit2->Text = 30;
			break;

			case 3:
			 mEnumShellType = CALIBRO_30;
			 menumCannonType = AK_630_2M;
			 ComboBox2->ItemIndex = 2;
			 LabeledEdit2->Text = 30;
			break;

			case 4:
			 menumCannonType = AK_176; // ��� ������� �� ���������!!! ����� ���� ��� ������� ��� ������!!!!
			 LabeledEdit2->Text = 76;
			break;

			default:
			mEnumShellType = CALIBRO_UNKNOWN;
			menumCannonType = CANNON_UNKNOWN;
			break;
		}



      // ���������� ����� � 76 ���������
	  if (menumCannonType == AK_176 )
	  {
	   switch( ComboBox4->ItemIndex)
	   {
		   case 4:
		   mEnumShellType = CALIBRO_76_SHTAT;
		   break;


		   case 5:
		   mEnumShellType = CALIBRO_76_BARRIER;
		   break;

		   default:
		   break;
       }

	  }


	if( (mEnumTargType == GARPUN_V700) || (mEnumTargType == GARPUN_V300) )
	{
		if ((menumCannonType == A_192M) || (menumCannonType == A_190_01)|| (menumCannonType == AK_176))
		{
			mFireBegin = 6050.;
		 //	mFireFinish = 800.;
			mFireFinish = 1000.;
		}
		if ( (menumCannonType == AK_630) || (menumCannonType == AK_630_2M) )
		{
			mFireBegin = 1650.;
			mFireFinish = 500.;
			mDetonatorType = CONTACT;
			ComboBox2->ItemIndex = 2;
		}
	}

	if( mEnumTargType == PLANE)
	{
		if ((menumCannonType == A_192M) || (menumCannonType == A_190_01)|| (menumCannonType == AK_176))
		{
		 	mFireBegin = 10000.;
			mFireFinish = 5800.;
		}

		if ( (menumCannonType == AK_630) || (menumCannonType == AK_630_2M) )
		{
			mFireBegin = 1650.;
			mFireFinish = 500.;
		}
	}
	

	///


		switch(ComboBox2->ItemIndex)
		{
			case 0:
			mDetonatorType = AR32A;
			break;

			case 1:
			mDetonatorType = MFIVU;
			break;

			case 2:
			mDetonatorType = CONTACT;
			break;

			case 3:
			mDetonatorType = DVM;
			break;

			case 4:
			mDetonatorType = AR31A;
			break;

			case 5:
			mDetonatorType = BARRIER;
			break;

			default:
			mDetonatorType = DETONATOR_UNKNOWN;
			break;
		}


		///
		double arrDetonatorParams [LEN_DOUBLE_ARR_DETONATOR_PARAMS] = {0.};
		int iarrDetonatorParams [LEN_DOUBLE_ARR_DETONATOR_PARAMS] = {0};

		arrDetonatorParams[0] =  StrTo_Dbl_(LabeledEdit3->Text)/ 1000.; // �������
		TDetonator Detonator (mDetonatorType,arrDetonatorParams, iarrDetonatorParams);

		mShellBody = TShellBody(mEnumShellType, Detonator ) ;
		mArtCannon =  TArtCannon (menumCannonType, mAUDelayT )  ;
		mArtComplex  = TArtComplex (mArtCannon, 	mSigDrivAY_U , mSigDrivAY_dU_po_dt )  ;

		

	LabeledEdit19->Text = ((double)(int(mFireBegin /1000.* 100.)))/ 100.;
	LabeledEdit20->Text = ((double)(int(mFireFinish /1000.* 100.)))/ 100.;
	LabeledEdit49->Text =  mArtCannon.mRateOfFire;

	mQuantShells = StrToInt(LabeledEdit10->Text);

	mInitTargData = TInitTargData(mBearing0, mTargCourse0, mTargZenitAng0,
		mVelocity0, mDist0, mElev0, 0.);
   const double TCur0 = 0.;
	double arrWSkz[3] = {0.};
	arrWSkz [0] = mWSkz;
	arrWSkz [1] = mWSkz;
	arrWSkz [2] = mWSkz;
	TTraject Traj0(TCur0,  arrWSkz, mInitTargData );
	TTarget Target(Traj0,mEnumTargType, mTargEPR, NULL) ;
	mVessel  =TVessel(        mShellBody, mFar_2D , mTransmitAnt
								, mDriverSigBet, mDriverSigEps, mDriverDynamicSigBet, mDriverDynamicSigEps
								, mMaxSig_Q,  mMaxSig_Psi,  mMaxSig_Tet, mMaxSig_dQdt,  mMaxSig_dPsidt, mMaxSig_dTetdt   // ����
								, mMaxSig_H, mMaxSig_VH, mK1, mSigV, mEnvironment
								, mVesselWidth , mVesselLength, marrFarParallacs, mMaxQ ,mT_Q
								, mMaxPsi, mT_Psi , mMaxTet, mT_Tet, mMaxVert,  mQ0, mVVess, mInitTargData
								, mMaxAmp_AftFlexure, mT_AftFlexure, mMaxAmp_BoardFlexure, mT_BoardFlexure
								,  mControlSyst,marrArtParral, mArtComplex, NULL ) ;
	mFight = TFight( mVessel, Target ,mVessel.mControlSyst.mFiltT,mEtalonSign,  mEnvironment, mAntCoeff,NULL ) ;

	 bool bLat = true;
	 double valSigE = -1.,valSigQ = -1.;
	 double valZAhvatDist  = mDist0;

	   double valHAntenna = 20.;
		 valZAhvatDist = mFight.mVessel.mFar_2D.calc_TwoTargsZahvatDist(mElev0, mTargEPR
		,mVessel. mTransmitAnt.mPowerPrd, mVessel.mTransmitAnt.mKYPrd, mEtalonSign , valHAntenna
	  , &bLat  , &valSigE ,&valSigQ  ) ;




	int ib = (valSigE * 100000.);
	LabeledEdit33->Text = ((double) ib) /100.;

	ib =  valZAhvatDist ;
	LabeledEdit31->Text = ib;

	TFar Far(mVessel.mFar_2D, true);
   const double VAlTetta0 = findDiagrWidth(Far.mFaceta.m_d  ,Far.m_D, Far.mFaceta.m_n
  , Far.m_N ,mVessel.mFar_2D.mLambda);
  ib = (VAlTetta0 * 100000.);
	LabeledEdit34->Text = ((double) ib) /100.;

	LabeledEdit51->Text = double( int(mVessel.mArtComplex.mArtCannon.mAngGroupedFire * 100000.) )/ 100.;

	char arrch[2] = {0};
	String s_22 = arrch;
	LabeledEdit31->Text =s_22;
	LabeledEdit33->Text =s_22;
	LabeledEdit34->Text =s_22;
	LabeledEdit45->Text =s_22;
	LabeledEdit46->Text =s_22;
	LabeledEdit47->Text =s_22;

}






void __fastcall TForm3::Button5Click(TObject *Sender)
{
   OpenDialog1->Filter = L"����� � ��������� (*.csv)|*.csv" ;

	if (OpenDialog1->Execute())
	{
	mpwchOutFileAppointmentPoints =  (OpenDialog1->FileName).w_str();

	}
	Edit1->Text =  mpwchOutFileAppointmentPoints;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void __fastcall TForm3::Button4Click(TObject *Sender)
{
   OpenDialog1->Filter = L"����� � ��������� (*.shp)|*.shp" ;

	if (OpenDialog1->Execute())
	{
	mpwchOutFile0 =  (OpenDialog1->FileName).w_str();

	}
	Edit2->Text =  mpwchOutFile0;
}
//---------------------------------------------------------------------------

void __fastcall TForm3::Button2Click(TObject *Sender)
{


			// �������� ���� � ����� � �������
	if (!mpwchOutFileAppointmentPoints)
	{
	 if (I_ENTER_COUNT > 4) Close();


	 ShowMessage(L"������� ���� � ����� � ���������") ;
	 I_ENTER_COUNT++;
	 return;
	}


	 fncInputData() ;

	int iNumCols = 9;
  int iNumReservedRows = 1000;
  int iLenName = 30;
  double *parrBuff = new double [iNumCols * iNumReservedRows ];

  wchar_t *pwcharrColNames = new  wchar_t[ iNumCols * iLenName];
  memset (pwcharrColNames, 0, iNumCols * iLenName * sizeof(wchar_t));


	wcscpy(pwcharrColNames, L"VISOTA");
	wcscpy(&pwcharrColNames[iLenName], L"DALN. HOR.");
	wcscpy(&pwcharrColNames[2 *iLenName], L"V CELI");
	wcscpy(&pwcharrColNames[3 *iLenName], L"V IZDEL.");
	wcscpy(&pwcharrColNames[4 *iLenName], L"UGOL XI, RAD");

	wcscpy(&pwcharrColNames[5 *iLenName], L"K11");
	wcscpy(&pwcharrColNames[6 *iLenName], L"K12");
	wcscpy(&pwcharrColNames[7*iLenName], L"K21");
	wcscpy(&pwcharrColNames[8 *iLenName], L"K22");





  // ���������� ������� ������� �������� �� ����
// VAlDistFireBegin - ��������� ����� ������� �������� ������
// VAlDelT - ��� �������
// ������� � ������� t=0 � �����  VAlDelT
//  �������� ������ � ����� ������� ��� ��������� �������� ���� � �������
//  ��� ������ ��������� ���������� ������  VAlDistFireBegin
//  ����� ������ �������� ������� ������������
	double valFireBeginTime = -1.;
	double valMostRemoteAppPointDist = -1.;
	if(!mFight.calcFirstShotTime_(mFireBegin,  &valFireBeginTime, &valMostRemoteAppPointDist))
  {
	  return;
  }
//  valFireBeginTime = 35.91;
	///////////////////////////////////////////////
	TFight FightCur = mFight;
	FightCur.mVessel.mShellBody.mDetonator.mEnumDetonatorType = CONTACT;
	FightCur.mVessel.mMaxQ = 0.; // ������� ��������� �� ��������� ������

 //	memset( FightCur.mTarget.mTraject.marrSigW, 0, 3 * sizeof(double));  // ���� ��������� �� ��������� ������

	int iQuantShots =0;

	double valKGSKEps = 0., valKGSKBet  = 0., valKGSKEpsPrev = 0., valKGSKBetPrev  = 0.;
	double arrVectAppointmentPointGSK[6] = {0.};
	double valMiss = -1., valTFlight = -1.;
	double valTCurrentShot = 0. ;


	for (int i = 0; i < iNumReservedRows; i++)
	{

	if (9 == i) {
    int jjj =0;
	}
	valTCurrentShot = valFireBeginTime + (double(i)) *  FightCur.mVessel.mArtComplex.mArtCannon.mRateOfFire;
	// ������������ ����� ������ FightCur �� �����  valTCurrentFight
	// ��������� ������������ ����������� ����������� ������ ��������
	// ����� ����������, �� ������� ����� FightCur ������ ����� ����� ���� ������  valTCurrentFight
	FightCur.shift(valTCurrentShot);
	///

	// ������������� �������� ��������� ���� �� �����  (valTCurrentFight  - FightCur.mT) ������
	double arrTargExtrapVS_GSK[9] = {0.};
	FightCur.mTarget.mTraject.extrapolateTargVS(valTCurrentShot - FightCur.mT, arrTargExtrapVS_GSK);
	///

	// ������������� �������� ��������� ������� �� �����  (valTCurrentFight - FightCur.mT) ������
	double arrVessExtrapVS_GSK[9] ={0.};
	FightCur.mVessel.extrapolateTrueVS_GSK(valTCurrentShot - FightCur.mT, arrVessExtrapVS_GSK);
	///


	// ���������� ��������� ���� ��� ��������� ���� � ���� �� �������� ��������
	double arrTargVS_KGSK0[6] ={0.};
	MtrxMinusMatrx(arrTargExtrapVS_GSK, arrVessExtrapVS_GSK,1, 6, arrTargVS_KGSK0);
	///

	// ���������� ������� ��������� �� � ����
	double arrPositionAY_KGSK[3] = {0.};
	FightCur.mVessel.calcAY_Position(arrPositionAY_KGSK);
	///

	double arrShellVeloAppointmPoint_GSK[3] = {0.}, arrTargVeloAppointmPoint_GSK [3]= {0.};
	// ���������� ����� �������
	if (i ==0)
	{
	FightCur.calcAppointmentPoint(NULL, NULL, &arrVessExtrapVS_GSK[3]
	,arrTargVS_KGSK0, arrPositionAY_KGSK
	, &valKGSKEps, &valKGSKBet, &valTFlight, arrVectAppointmentPointGSK, &valMiss
	, arrShellVeloAppointmPoint_GSK, arrTargVeloAppointmPoint_GSK) ;
	valKGSKEpsPrev =  valKGSKEps -0.05;
	valKGSKBetPrev = valKGSKBet;

	}
	else
	{
	FightCur.calcAppointmentPoint(&valKGSKEpsPrev, &valKGSKBetPrev, &arrVessExtrapVS_GSK[3]
	,arrTargVS_KGSK0, arrPositionAY_KGSK
	, &valKGSKEps, &valKGSKBet, &valTFlight, arrVectAppointmentPointGSK, &valMiss
	, arrShellVeloAppointmPoint_GSK, arrTargVeloAppointmPoint_GSK) ;
	valKGSKEpsPrev =  valKGSKEps -0.05;
	valKGSKBetPrev = valKGSKBet;
	}

	// ������������� �������� ��������� ������� �� �������� �����  valTFlight ������
	double arrVessExtrapVS_GSK_0[9] ={0.};
	MatrxMultScalar(&arrVessExtrapVS_GSK[3], 3, 1, valTFlight,arrVessExtrapVS_GSK_0);

	///
	// ��������� ����� ������� � ���� �� ������ �������
	double arrT0[3] ={0.};
	MtrxMinusMatrx(arrVectAppointmentPointGSK , arrVessExtrapVS_GSK_0,1, 3, arrT0);
	double valDist  = Norm3(arrT0);
	if(valDist < mFireFinish)
	{
	break;
	}

	double valProb = -1, valDispMiss = -1.,valDispNedolet = -1., valSKZ_GSK_Z = -1.0;
	double arrCorMtrxCartinSK[4] ={0.};
	TURPolygon  plgProjection ;
	FightCur.calcProbability_For_Fixed_AppointmentPoint_AirTargs(valTCurrentShot, valTFlight, valKGSKEps, valKGSKBet
		,&valProb,  &valDispMiss, &valDispNedolet,arrCorMtrxCartinSK, &valSKZ_GSK_Z,  &plgProjection);

	iQuantShots++;

	int ii = int (FightCur.mTarget.mTraject.marrVectSostGSK_Begin[2] * 100.) ;
	parrBuff [i * iNumCols] = (double (ii)) /100.;

	ii = int( sqrt(arrT0[0] * arrT0[0] + arrT0[1] * arrT0[1]) * 100.);
	parrBuff [i * iNumCols + 1] = (double (ii)) /100.;

	ii = int(Norm3(arrTargVeloAppointmPoint_GSK) * 100.);
	parrBuff [i * iNumCols + 2] =  (double (ii)) /100.;

	ii = int( Norm3(arrShellVeloAppointmPoint_GSK) * 100.);
	parrBuff [i * iNumCols + 3]  =  (double (ii)) /100.;


	double temp0 =  ScalProduct(arrTargVeloAppointmPoint_GSK , arrShellVeloAppointmPoint_GSK, 3);
	double temp1 = parrBuff [i * iNumCols + 2]*parrBuff [i * iNumCols + 3];
	double temp2 = temp0/ temp1;
	 if (temp2 < -1.)
	{
		temp2 = -0.999999999999;
	}
	if (temp2 > 1.)
	{
		temp2 = 0.999999999999;
	}
	ii = int( acos( temp2) * 1000.);
 // 	ii = int( acos(ScalProduct(arrTargVeloAppointmPoint_GSK , arrShellVeloAppointmPoint_GSK, 3)
 //	 /parrBuff [i * iNumCols + 2]/parrBuff [i * iNumCols + 3]) * 1000.);
		parrBuff [i * iNumCols + 4] =   (double (ii)) /1000.;

	 ii = int( arrCorMtrxCartinSK[0] * 100.);
	 parrBuff [i * iNumCols + 5] =   (double (ii)) /100.;

	 ii = int( arrCorMtrxCartinSK[1] * 100.);
	 parrBuff [i * iNumCols + 6] =   (double (ii)) /100.;

	 parrBuff [i * iNumCols + 7] = parrBuff [i * iNumCols + 6];

	 ii = int( arrCorMtrxCartinSK[3] * 100.);
	 parrBuff [i * iNumCols + 8]=   (double (ii)) /100.;
	}

 wchar_t *pwcharrRowNames = new wchar_t [iQuantShots * iLenName];
 memset(pwcharrRowNames, 0, iQuantShots * iLenName * sizeof(wchar_t));
 for (int i = 0; i < iQuantShots; i++)
 {
	 swprintf (&pwcharrRowNames[ i *iLenName] , L"%d",i + 1);
 }

	 // ������ ������� ���������� �  CSV ����
// INPUT:
// FileName
// parrBuff[ iNumRows * iNumCols] - ������ � �����������
// iNumRows- �-�� ����� �������
// iNumCols - �-�� �������� �������
// pwcharrRowNames[ iLenName* iNumRows] - ����� ����� �������
// pwcharrColNames [ iLenName* iNumCols] - ����� �������� �������
 TYrWrite::WriteMassiveInFIleSCV(mpwchOutFileAppointmentPoints,parrBuff,  iQuantShots,  iNumCols
							 ,pwcharrRowNames,pwcharrColNames, iLenName) ;

 delete parrBuff;
 delete pwcharrColNames;
 delete pwcharrRowNames;


}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void __fastcall TForm3::Button3Click(TObject *Sender)
{

					// �������� ���� � ����� � �������
	if (!mpwchOutFile0)
	{
	 if (I_ENTER_COUNT > 4) Close();


	 ShowMessage(L"������� ���� � ����� � ���������") ;
	 I_ENTER_COUNT++;
	 return;
	}
	wcscpy(mwchOutFold,mpwchOutFile0);
	wchar_t *pwchr = wcsrchr(mwchOutFold, L'\\');
	pwchr[0] = 0;
	   //	String wchFoldName = mwchOutFold;

	 fncInputData() ;


	 double valProb_Gladk = -1.0,  valProb = -1.0;
	 double valSigE = -1., valSigQ = -1., valDistBeginSopr = -1;;

	 // ������ �������������

	 double valFireBeginTime = -1.;
	 mQuantShots = -1;


	 memset(mparrProbab, 0, ILenArr * sizeof(double));

	 memset(mparrSKZPromach, 0, ILenArr * sizeof(double));

	 memset(mparrSKZNedolet, 0, ILenArr * sizeof(double));

	 double valHAntenna = 20.;

	 memset(mparrSKZ_GSK_Z, 0, ILenArr * sizeof(double));

	 memset(mparrCorMtrxCartinSK , 0, ILenArr * 4 *sizeof(double));
	 memset (mparrDist, 0, ILenArr * 4 *sizeof(double));

 	 memset(mpPlgArrProjection , 0, ILenArr *sizeof(TURPolygon));
		 if(!	mFight.calcSuccessProbAero(&mFireBegin,mFireFinish
			,valHAntenna, &valDistBeginSopr, &valSigE, &valSigQ,&valFireBeginTime
			,&mQuantShots , mparrProbab,  mparrSKZPromach, mparrSKZNedolet, mparrSKZ_GSK_Z
			,mparrCorMtrxCartinSK, mparrDist,  mpPlgArrProjection
			,ILenArr, &valProb, &valProb_Gladk)
		 )
	 {
		LabeledEdit46->Text = 0;
		LabeledEdit45->Text = 0;
		return;
     }

	 //
	int ia = valDistBeginSopr ;
	LabeledEdit31->Text = ia;

	ia = (valSigE * 100000);
	 LabeledEdit33->Text  =(( double)ia)/ 100.;
	///
		ia = (valProb * 100);
		LabeledEdit46->Text  =(( double)ia)/ 100.;

		ia = (valProb_Gladk * 100);
		LabeledEdit47->Text  =(( double)ia)/ 100.;

		LabeledEdit45->Text =  mQuantShots;

		ia = mparrDist[mQuantShots -1];
		LabeledEdit1->Text =  ia;

		ia = mparrDist[0];
		LabeledEdit4->Text =  ia;

		Button7->Visible = true;
		Edit3->Visible = true;
		mbCalculated = true;

		// ������ ����������� �� ��������� �� ����� �������
		wchar_t wchFoldName[300] ={0}, wchFileName[300] ={0};
   //	, wchFileNamePolygonEll1[300] ={0}, wchFileNamePolygonEll2[300] ={0}
  //	, wchFileNamePolygonEll3[300] ={0};
	wcscpy(  wchFoldName,  mwchOutFold);
	wcscat(wchFoldName, L"\\");

	wcscpy(  wchFileName,  wchFoldName);
	wcscat(wchFileName, L"Probab.shp");
	double scalex = 0.1, scaley = 100.;
		TYrWriteShapeFile::CreateShpFile(wchFileName, mparrProbab, mparrDist
	 ,mQuantShots,scalex , scaley);
	 ///


	 // ������ SKZ ������� �� ��������� �� ����� �������
	wcscpy(  wchFileName,  wchFoldName);
	wcscat(wchFileName, L"SKZ_Promach.shp");
	scalex = 0.1;
	scaley = 1.;
		TYrWriteShapeFile::CreateShpFile(wchFileName, mparrSKZPromach, mparrDist
	 ,mQuantShots,scalex , scaley);
	 ///

	 // ������ SKZ �������� �� ��������� �� ����� �������
	wcscpy(  wchFileName,  wchFoldName);
	wcscat(wchFileName, L"SKZ_Nedolet.shp");
	scalex = 0.1;
	scaley = 1.;
	TYrWriteShapeFile::CreateShpFile(wchFileName, mparrSKZNedolet, mparrDist
	 ,mQuantShots,scalex , scaley);
	 ///
		wcscpy(  wchFileName,  wchFoldName);
	wcscat(wchFileName, L"Axes.shp");
	TYrWriteShapeFile::CreateShpAxes(wchFileName,-10.,10000.
	 ,-1.,1000.);


	///

	// ���������� ������� ����������� ��������� � ����������� �� �-�� ��������� ���������
	// ������ ������������ ��������� ��������� ��������� ���������, ��� ������������� ������ ��������� ����� ����
	// �� �� ����  ����������� ���� ������. �������, ��������� ������� ������� �-�� ����������� ��������� ���������

	double *parrProbAccumulated  = new double  [mQuantShots];
	memset(parrProbAccumulated, 0, sizeof(double) * mQuantShots);
	double *parrProbGlagkAccumulated  = new double  [mQuantShots];
	memset(parrProbGlagkAccumulated, 0, sizeof(double) * mQuantShots);

	for (int i = 0; i < mQuantShots ; i++)
	{
		 mFight.calcAeroTargRezultProbability(&mparrProbab[i], &mparrDist[i], mQuantShots - i
	, &parrProbAccumulated[i],&parrProbGlagkAccumulated [i]) ;
	}

	wcscpy(  wchFoldName,  mwchOutFold);
	wcscat(wchFoldName, L"\\");

	wcscpy(  wchFileName,  wchFoldName);
	wcscat(wchFileName, L"ProbAccumulated.shp");
  scaley = 100.;
	TYrWriteShapeFile::CreateShpFile(wchFileName, parrProbAccumulated, mparrDist
	 ,mQuantShots,scalex , scaley);

	wcscpy(  wchFoldName,  mwchOutFold);
	wcscat(wchFoldName, L"\\");

	wcscpy(  wchFileName,  wchFoldName);
	wcscat(wchFileName, L"ProbGlagkAccumulated.shp");

	TYrWriteShapeFile::CreateShpFile(wchFileName, parrProbGlagkAccumulated, mparrDist
	 ,mQuantShots,scalex , scaley);
	// �������� ���������� ���������
	TURPolyLine plnLinDiagr = TURPolyLine::createLineDiagram(
		 mparrDist , parrProbGlagkAccumulated
		 , mQuantShots);
	plnLinDiagr.stretchDiagrAlongXY(scalex, scaley) ;
	wcscpy(  wchFoldName,  mwchOutFold);
	wcscat(wchFoldName, L"\\");

	wcscpy(  wchFileName,  wchFoldName);
	wcscat(wchFileName, L"LineDiagr.shp");
	plnLinDiagr.WriteSetSHPFiles(wchFileName, &plnLinDiagr, 1);
	delete parrProbAccumulated;
	delete parrProbGlagkAccumulated;

	


		return;


}
//---------------------------------------------------------------------------




void __fastcall TForm3::ComboBox3Change(TObject *Sender)
{
fncInputData();
}
//---------------------------------------------------------------------------

void __fastcall TForm3::ComboBox1Change(TObject *Sender)
{
fncInputData();
}
//---------------------------------------------------------------------------

void __fastcall TForm3::ComboBox2Change(TObject *Sender)
{
fncInputData();
}
//---------------------------------------------------------------------------





void __fastcall TForm3::LabeledEdit22Change(TObject *Sender)
{

 fncInputData();
}
//---------------------------------------------------------------------------

void __fastcall TForm3::Panel7Click(TObject *Sender)
{
fncInputData();
}
//---------------------------------------------------------------------------

void __fastcall TForm3::Panel10Click(TObject *Sender)
{
//fncInputData();
}
//---------------------------------------------------------------------------

// ---------------------------------------------------------------------------

void __fastcall TForm3::create5P10()
{
   	//���������
	 double valEtalonAmp = 500.;
	//���������
	 double valEtalonDist = 12000.;
	//���
	 double valEtalonAPR = 1.;
	//��� ����� ���� ��������� ��������� 5�10
	 double valNoiseSKZ_5P10 = 18.735;
	// ��� �������� ������� �������� ��������� ��������� 5�10
	 double valEtalonSigAmplFact_5P10 = 0.01;
	//
	// �������� �� ��������
	 double valEtalonPowerPrd = 4000.;
	// �� �� ��������
	 double valEtalonKYPrd = 840.;
	// �������� �� �����
	 double valEtalonKYPriem = 5200.;

	 mEtalonSign = TEtalonSign( valEtalonAmp, valEtalonDist,  valEtalonAPR,
			 valNoiseSKZ_5P10, valEtalonSigAmplFact_5P10, valEtalonPowerPrd,  valEtalonKYPrd
			,  valEtalonKYPriem);

	///

	// �������
	// ������ �� ��
	// �-�� ����������� �� �����������
	 int iNumEmitCols = 8;
	// �-�� ����������� �� ���������
	 int iNumEmitRows = 8;
		// ����� �����
		double valLambda = VAL_C / 8.* 100. / 1000000000.;

			// ���������� ����� ������������
		double valdEmitCol = valLambda * 0.55;
		double valdEmitRow = valLambda * 0.55;


	// ������ �� ���
	// �-�� ��  �� �����������
	 int iNumAMCols = 8;
	// �-�� �� �� ���������
	 int iNumAMRows = 8;
		// ���������� ����� �� �� ���������
	double valdAMRow = ((double)iNumEmitRows) *  valdEmitRow;
	// ���������� ����� �� �� �����������
	double valdAMCol = ((double)iNumEmitCols) *  valdEmitCol;
	// ��� ���� � ��������� ���������
	bool barrAM [5000] = {0};
	 for (int i = 0; i < 5000; i++)
	{
		barrAM[i] = true;
	}


//	double arrArtParral[3] ={0.},  arrFarParallacs[3] ={0.};



	// �������� �� ��������
	 double valPowerPrd =4000.;
	// �� �� ��������
	 double valKYPrd = 3000.;
	 mTransmitAnt.mKYPrd  =valKYPrd;
	 mTransmitAnt.mPowerPrd = valPowerPrd;


	///

const double SigEmitNoise =  valNoiseSKZ_5P10 / sqrt(((double)(8 * 16 * 28)));
	const double SigEmitAmplFact = valEtalonSigAmplFact_5P10 / sqrt
		(((double)(8 * 16 * 28)));
 TAM_2D AM_2D(iNumEmitCols, iNumEmitRows,  valdEmitCol,  valdEmitRow,
		SigEmitNoise, SigEmitAmplFact);

	mFar_2D = TFar_2D(iNumAMCols, iNumAMRows,  valLambda,  valdAMCol,  valdAMRow,
		AM_2D, barrAM);

	for (int i = 0; i < iNumAMCols * iNumAMRows; i++) {

		if (!barrAM[i]) {
			mFar_2D.mpAm2D[i].mOtklCoefUs = 1.;
			mFar_2D.mpAm2D[i].mSigEmitNoise = 0.;
		}
	}


	//-------------------------------------------------------------------
	// ��� ������� ������ ���� (����� �����)
	mSigSins = 0.00041;

	// ��� ������� ������ ������ �������� ����� (����� �����)
	mSig_d_po_dt_Sins = 0.00116;

	// ������� ���
	 mVesselWidth = 0.; // ������(�)
	 mVesselLength = 0. ;

	mMaxQ =     3./180.*M_PI; /// ������������ ���� ���������� �� ������������ �����(��������� ���� ��������)
	mT_Q = 18.; // ������ ��������
	mMaxPsi =      3./180.*M_PI;// ������������ ���� ������� �����(���������)
	mT_Psi = 12; // ������ ������� �����
	mMaxTet =      12./180.*M_PI; //������������ ���� ��������� �����(���������)
	mT_Tet = 6; // ������ �������� �����
	mMaxVert =     1. ;

	// ���������� ��������  ������� ������
	 mQ0 = 0. ; // ����������� ����
	 mVVess = 20. * 0.514 ;// �������� ������� ������ 20 �����
	// double arrDelt[4] = {0.};


	mMaxAmp_AftFlexure  = 1. * M_PI/180.;
	// ������ ��������� ��������� ������
	mT_AftFlexure = 4.;
	//������������ ��������� ��������� ������ ������� � ��� �� 100 �
	mMaxAmp_BoardFlexure =  1. * M_PI/180.;
	// ������ ��������� ��������� ������
	mT_BoardFlexure = 2.;

		 // 3.1 �������� ����
	mMaxSig_Q      =     mSigSins ; //0.000582;
	mMaxSig_Psi    =      mSigSins ; //0.00145;
	mMaxSig_Tet    =    mSigSins ; // 0.00145;
	mMaxSig_dQdt   =      mSig_d_po_dt_Sins ;
	mMaxSig_dPsidt =      mSig_d_po_dt_Sins ;
	mMaxSig_dTetdt =      mSig_d_po_dt_Sins ;
	mK1         = 0.01 ;
	mSigV       =      0.2 * sqrt(2.) ;
	mSigH       =      0.1 ;
	mMaxSig_H =     0.1 ;
	mMaxSig_VH =     0.05 ;

    	// 3.3 �������� �������
	mDriverSigBet  =      0.00021 ;
	mDriverSigEps  =      0.00021 ;
	mDriverDynamicSigBet =      0.0003141;
	mDriverDynamicSigEps =      0.0003141;

		// �������� ����� �����������
	mMeasT = 0.02;

	// �������� ����
	mSinsDelayT = 0.02;
  		// ���� ������� ���
	mRzvT = 0.00001;
	mControlSyst = TControlSyst(mMeasT, mSinsDelayT,mRzvT );

}

void __fastcall TForm3::Button6Click(TObject *Sender)
{
	// �������� ���� � ����� � �������
	if (!mpwchOutFile0)
	{
	 if (I_ENTER_COUNT > 4) Close();


	 ShowMessage(L"������� ���� � ����� � ���������") ;
	 I_ENTER_COUNT++;
	 return;
	}
	wcscpy(mwchOutFold,mpwchOutFile0);
	wchar_t *pwchr = wcsrchr(mwchOutFold, L'\\');
	pwchr[0] = 0;
		String wchFoldName = mwchOutFold;

	 fncInputData() ;
	 TURPolygon arrPlg0 [13], plg[3] ;


	 switch(mFight.mTarget.menumTargetType)
	 {
		 case GARPUN_V300:
		 case GARPUN_V700:

		 for (int i = 0; i < 13; i++)
		 {
			arrPlg0 [ i] = mFight.mTarget.mpArrPlanePolygon[i].mPolygon;
		 }
		 plg[0] = TURPolygon( arrPlg0, 13);
		 plg[1] = plg[0];
		 plg[2] = mFight.mTarget.mpArrPlanePolygon[13].mPolygon;
		 break;
		 case PLANE:
		 plg[0] =	mFight.mTarget.mpArrPlanePolygon[1].mPolygon;
		 plg[1] = mFight.mTarget.mpArrPlanePolygon[0].mPolygon;
		 plg[2] = mFight.mTarget.mpArrPlanePolygon[2].mPolygon;
		 break;
		 default:
		 break;
   }


		 plg[0].calcBoundBox()  ;
		 TURPointXY pntSdvig(-plg[0].Box[2]- 1.,- plg[0].Box[1]+1.);
		 TURPolygon plgTEmp = plg[0].SdvigTransform( pntSdvig );
		 plg[0] =plgTEmp ;

		 plg[1].calcBoundBox()  ;
		 pntSdvig = TURPointXY (-plg[1].Box[2]- 1.,- plg[1].Box[3]-1.);
		 plgTEmp = plg[1].SdvigTransform( pntSdvig );
		 plg[1] =plgTEmp ;

		 plg[2].calcBoundBox()  ;
		 pntSdvig = TURPointXY (-plg[2].Box[0]+ 1.,- plg[2].Box[1]+1.);
		 plgTEmp = plg[2].SdvigTransform( pntSdvig );
		 plg[2] =plgTEmp ;
	 for (int i =0; i <  3; i++)
	 {
	 //	 plg[i].WriteSetSHPFiles(L"E:\\�����\\F14\\plgi.shp", &plg[i], 1) ;
		 String wchFileName = wchFoldName ;
		 wchFileName += L"\\Image";
		 wchFileName += i;
		 wchFileName += L".shp";
		 plg[i].WriteSetSHPFiles(wchFileName.w_str(), &plg[i], 1) ;
	 //	 mFight.mTarget.mpArrPlanePolygon[i].mPolygon.
	 }


}
//---------------------------------------------------------------------------

void __fastcall TForm3::Button7Click(TObject *Sender)
{
	mNumShot =StrToInt(Edit3->Text)-1;

  // �������� ����������
  String stringOutDir = mwchOutFold;
  stringOutDir += L"\\AppointPoint_ ";
  stringOutDir +=  mNumShot;
  _wmkdir(stringOutDir.w_str());
  ///

	// 1. ������� ���� � ���
	TURPolygon *pPlgArr = new TURPolygon[mpPlgArrProjection[mNumShot].NumParts];
	for (int i = 0; i < mpPlgArrProjection[mNumShot].NumParts; i++)
	{
	 pPlgArr[i] =  mpPlgArrProjection[mNumShot].extractSimplePolygon(i);
	}
	wchar_t wchFoldName[300] ={0}, wchFileNamePlg[300] ={0}
	, wchFileNamePolygonEll1[300] ={0}, wchFileNamePolygonEll2[300] ={0}
	, wchFileNamePolygonEll3[300] ={0};
	wcscpy(  wchFoldName,  stringOutDir.w_str());
	wcscat(wchFoldName, L"\\");

	wcscpy(  wchFileNamePlg,  wchFoldName);
	wcscat(wchFileNamePlg, L"TargPlg.shp");
	TURPolygon::WriteSetSHPFiles(wchFileNamePlg,pPlgArr, mpPlgArrProjection[mNumShot].NumParts) ;

	delete [] pPlgArr ;
	///

	// 2. ���������� ����������� ���������� � ��������� ���������
	double *parrElK = &mparrCorMtrxCartinSK[ 4 * mNumShot];
	//   arrF - ������� ����������� ��������  ������ �������
	double arrF[4] = {0.} , arrMtrxLamb[4] = {0.};
	CalcProperVectors2(parrElK, arrF , arrMtrxLamb) ;

	 ///

	 // ������� ��������� ��������������, ������������ ������   ��������� ��������� ����� �������
	 // �� ����������� ������� ����������������� ��������� ��������� ���������
	double arrLinTrasf [4] ={0.}, arrLambSq[4] ={0.};
	arrLambSq[0] =  sqrt(arrMtrxLamb[0]);
	arrLambSq[3] =  sqrt(arrMtrxLamb[3]);
	MtrxMultMatrx(arrF,2, 2, arrLambSq,2, arrLinTrasf) ;
	///
	TURPointXY pntSdvig(0.,0.);
	// 2.1 ������������ ������� ��������� �� ������ 1

	TURPolygon plgCircle1 = TURPolygon::fncCreateCircle(pntSdvig,1, 1001) ; // ��� ��������� ����
	TURPolygon plgonEll1 = plgCircle1.fncLinTransform(arrLinTrasf );// ��� ��� �������� ��������������
	wcscpy(  wchFileNamePolygonEll1,  wchFoldName);
	wcscat( wchFileNamePolygonEll1, L"PolygonEll1.shp");

	plgonEll1.WriteSetSHPFiles( wchFileNamePolygonEll1, &plgonEll1,1);

	///

	 // 2.2 ������������ ������� ��������� �� ������ 2
	TURPolygon plgCircle2 = TURPolygon::fncCreateCircle(pntSdvig
	,2, 1001) ;  // ��� ����  ������� 2

	TURPolygon plgonEll2 = plgCircle2.fncLinTransform(arrLinTrasf ); // ��� ��� �������� ��������������

	wcscpy(  wchFileNamePolygonEll2,  wchFoldName);
	wcscat( wchFileNamePolygonEll2, L"PolygonEll2.shp");

	plgonEll2.WriteSetSHPFiles( wchFileNamePolygonEll2, &plgonEll2,1);

	///

	// 2.3 ������������ ������� ��������� �� ������ 3
	TURPolygon plgCircle3 = TURPolygon::fncCreateCircle(pntSdvig
	,3, 1001) ;// ��� ����  ������� 3

	TURPolygon plgonEll3 = plgCircle3.fncLinTransform(arrLinTrasf ); // ��� ��� �������� ��������������
	wcscpy(  wchFileNamePolygonEll3,  wchFoldName);
	wcscat( wchFileNamePolygonEll3, L"PolygonEll3.shp");

	plgonEll3.WriteSetSHPFiles( wchFileNamePolygonEll3, &plgonEll3,1);
	///

	// 3. ���  ���������
	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  wchFoldName);
	wcscat(wchAxesFileName, L"AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-500., 500.
	,-500., 500.,30.) ;
}
//---------------------------------------------------------------------------





void __fastcall TForm3::Button8Click(TObject *Sender)
{
				// �������� ���� � ����� � �������
	if (!mpwchOutFileTraj0)
	{
	 if (I_ENTER_COUNT > 4) Close();


	 ShowMessage(L"������� ���� � ����� � ���������") ;
	 I_ENTER_COUNT++;
	 return;
	}
	wcscpy(mpwchOutFileTraj,mpwchOutFileTraj0);
	wchar_t *pwchr = wcsrchr(mpwchOutFileTraj, L'\\');
	pwchr[0] = 0;
	   //	String wchFoldName = mwchOutFold;

	 fncInputData() ;
	 const double VAlEps0  = StrTo_Dbl_(LabeledEdit6->Text) * M_PI / 180.;

	 // ���������� �������� ���������� ���� � ��������� ���������� �� ����������
	 wchar_t wchFoldName[400] = {0};
	wcscpy(  wchFoldName,  mpwchOutFileTraj);
	wcscat(wchFoldName, L"\\Traject_and_Scatters");
	_wmkdir (wchFoldName);
	double arrVessVelocity_GSK [3] ={0.};
	TMyShellTraj ShellTraj (arrVessVelocity_GSK, mVessel.mShellBody,  VAlEps0,  0. );
	const double VAlStepInt =    0.005;
	double valDHoriz = 1.;
	ShellTraj.fncMoveClass_TO_ZeroAlt_AND_ShowGraphs(mEnvironment, VAlStepInt, wchFoldName,  valDHoriz ) ;

	 LabeledEdit13->Text =   ShellTraj.mTCur;
	 LabeledEdit11->Text =   sqrt(ShellTraj.marrStrSK_VS [0]* ShellTraj.marrStrSK_VS [0] +  ShellTraj.marrStrSK_VS [1] * ShellTraj.marrStrSK_VS [1]);
	 LabeledEdit9->Text =   ShellTraj.marrStrSK_VS [6] ;
	 LabeledEdit14->Text =   ShellTraj.marrStrSK_VS[7] / M_PI * 180. ;

	// ������������ ������� ���������
		  // ���������� �������������� ������� �������� ������� ��������� �������
double arrMtrxShellDisp  [LEN_ARR_SCATTERS * LEN_ARR_SCATTERS];
// 0. marrStrSK_VS [3]-  ���� ���
// 1. marrStrSK_VS [5]-  ������� ��������� ��������� ��
// 2.  marrStrSK_VS [7]-  ���� �����
// 3. marrStrSK_VS [6]-  ������� ��������
// �������������� ��������:
// 4.  �����
// 5.  ����� ����� Cx
// 6. ���� ����� �� ��� Z
// 7. ������ �������������� �������� �����
// 8. ����������� ��������������� �����
// 9. ������ ������������� �����
 memset(arrMtrxShellDisp, 0,  sizeof(double) * LEN_ARR_SCATTERS * LEN_ARR_SCATTERS);
 mFight.fillShellVozmDispMatr_ShootingEarth (arrMtrxShellDisp) ;

 ///
   double step = 0.1;
	const int QUANT_COLS = 10 , QUANT_POINTS  = (ShellTraj.mTCur) / step ;


	wchar_t *pwcharrFileNames ;
	const int lenName =30 ;// ������������ ����� ����� ����������




	double *parrBuff = new double [QUANT_COLS * QUANT_POINTS] ;
	memset (parrBuff, 0, QUANT_COLS * QUANT_POINTS * sizeof(double)) ;
	pwcharrFileNames = new wchar_t [ QUANT_COLS * lenName] ;
	memset (pwcharrFileNames, 0, QUANT_COLS * lenName* sizeof(wchar_t)) ;

	wcscpy( &pwcharrFileNames[ 0 * 30], L"t");
	wcscpy( &pwcharrFileNames[ 1 * 30], L"SigGSK_X");
	wcscpy( &pwcharrFileNames[ 2* 30],  L"SigGSK_Y");
	wcscpy( &pwcharrFileNames[ 3 * 30], L"SigGSK_Z");
	wcscpy( &pwcharrFileNames[ 4 * 30], L"SigGSK_Vx");
	wcscpy( &pwcharrFileNames[ 5 * 30], L"SigGSK_Vy");
	wcscpy( &pwcharrFileNames[ 6 * 30], L"SigGSK_Vz");
	wcscpy( &pwcharrFileNames[ 7 * 30], L"SigSSK_X");
	wcscpy( &pwcharrFileNames[ 8 * 30], L"SigSSK_Y");
	wcscpy( &pwcharrFileNames[ 9 * 30], L"SigSSK_Z");

	double *pscaleY = new double  [QUANT_COLS] ;
	pscaleY[1] = 1;
	pscaleY[2] = 1;
	pscaleY[3] = 1;
	pscaleY[4] = 1.;
	pscaleY[5] = 1.;
	pscaleY[6] = 1.;
	pscaleY[7] = 1.;
	pscaleY[8] = 1.;
	pscaleY[9] = 1.;
   double arrShellScatteringsCorMtarx_GSK [ 6 * 6] = {0.}, arrShellVS_GSK[6] = {0.},  arrShellScatteringsCorMtrxPos_SSK[9] = {0.}
   , arrStrSK_Jacobian [LEN_ARR_SCATTERS * 8] = {0.};

	for (int i = 0; i  < QUANT_POINTS; i++)
	{
	  TMyShellTraj ShellTraj0 (arrVessVelocity_GSK, mVessel.mShellBody,  VAlEps0,  0. );

	 // ���������� �������������� ������� �������� ������� ��������� �������
// � ���
// INPUT:
// Environment -������� �����
// VAlStepInt - ��� ��������������
// VAlTFlight - �������� �����
// arrMtrxShellDisp  - ������������ ������� ���������� ��������� ����������
//OUTPUT:
// arrStrSK_Jacobian [LEN_ARR_SCATTERS * 8] - ������� ������� ����������� �������� �������
// ������� � ������������� �������� ��������� �� ����������
// arrShellScatteringsCorMtarx_GSK [ 6 * 6] - �������������� ������� � ���
// arrShellScatteringsCorMtrxPos_SSK [3*3] - �������������� ������� ������ ��������� ��������� � ���
// arrShellVS_GSK[6] - ������ ��������� ������� � ���

ShellTraj0.calc_VS_GSK_And_ScatteringsCorrMatrx_GSK (mEnvironment, VAlStepInt
	 ,(double(i) )* step, arrMtrxShellDisp,  arrStrSK_Jacobian
	, arrShellScatteringsCorMtarx_GSK, arrShellVS_GSK,  arrShellScatteringsCorMtrxPos_SSK) ;
	 parrBuff[ i *QUANT_COLS] =  ShellTraj0.mTCur;
	 parrBuff[ i *QUANT_COLS + 1] =  sqrt(arrShellScatteringsCorMtarx_GSK[0]);
	 parrBuff[ i *QUANT_COLS + 2] =  sqrt(arrShellScatteringsCorMtarx_GSK[6 + 1]);
	 parrBuff[ i *QUANT_COLS + 3] =  sqrt(arrShellScatteringsCorMtarx_GSK[ 2 *6 + 2]);
	 parrBuff[ i *QUANT_COLS + 4] =  sqrt(arrShellScatteringsCorMtarx_GSK[ 3 *6 + 3]);
	 parrBuff[ i *QUANT_COLS + 5] =  sqrt(arrShellScatteringsCorMtarx_GSK[ 4 *6 + 4]);
	 parrBuff[ i *QUANT_COLS + 6] =  sqrt(arrShellScatteringsCorMtarx_GSK[ 5 *6 + 5]);
	 parrBuff[ i *QUANT_COLS + 7] =  sqrt(arrShellScatteringsCorMtrxPos_SSK[0]);
	 parrBuff[ i *QUANT_COLS + 8] =  sqrt(arrShellScatteringsCorMtrxPos_SSK[4]);
	 parrBuff[ i *QUANT_COLS + 9] =  sqrt(arrShellScatteringsCorMtrxPos_SSK[8]);


	}
   for (int j = 1; j < QUANT_COLS ; j++)
	 {


	TYrWriteShapeFile::WriteOneReport(wchFoldName  // ���� � �����
								  , parrBuff // ������ � ����������� - ������� nBuffRows x nBuffCols
								  ,QUANT_COLS // - �-�� ���������� � ������� ��������� ���������� � ������
								  ,QUANT_POINTS //  - �-�� �����
								  ,pwcharrFileNames //������� � �������� ���������� - ������� nBuffCols x lenName
								  ,lenName // ������������ ����� ����� ����������
								  ,0  // ����� ���������� �� ��� X
								  ,j  // ����� ���������� �� ��� Y
								  ,100 //  ������� �� ��� X
								  ,pscaleY[j]  // ������� �� ��� Y
								   );
	 }


  wchar_t wchFileName4[300] = {0} ;
  wcscpy(wchFileName4, wchFoldName );
  wcscat(wchFileName4, L"\\Axes.shp");
  TYrWriteShapeFile::CreateShpAxes(wchFileName4,-40000.,40000,-40000.,40000) ;



	 LabeledEdit7->Text =   sqrt(arrShellScatteringsCorMtarx_GSK[0]);
	 LabeledEdit8->Text =   sqrt(arrShellScatteringsCorMtarx_GSK[6 + 1]);
	 LabeledEdit12->Text =  sqrt(arrShellScatteringsCorMtarx_GSK[ 2 *6 + 2]);

   delete []parrBuff;
   delete []pwcharrFileNames;
   delete []pscaleY;
	/*
		const int QUANT_COLS = 9 , QUANT_POINTS  = (VAlFixedT - mTCur) / VAlStepInt ;
	double *parrBuff;
	double *pscaleY;
	wchar_t *pwcharrFileNames ;
	const int lenName =30 ;// ������������ ����� ����� ����������
	wchar_t 	wcharrPath[400] = {0};

   if (wcharrPath1)
   {
   wcscpy(wcharrPath, wcharrPath1 );
   if(wcharrPath[wcslen(wcharrPath) -1] != L'\\')
   {
	   wcharrPath[wcslen(wcharrPath) ] = L'\\';
	   wcharrPath[wcslen(wcharrPath) +1] = 0;

   }


	parrBuff = new double [QUANT_COLS * QUANT_POINTS] ;
	memset (parrBuff, 0, QUANT_COLS * QUANT_POINTS * sizeof(double)) ;
	pwcharrFileNames = new wchar_t [ QUANT_COLS * lenName] ;
	memset (pwcharrFileNames, 0, QUANT_COLS * lenName* sizeof(wchar_t)) ;

	wcscpy( &pwcharrFileNames[ 0 * 30], L"t");
	wcscpy( &pwcharrFileNames[ 1 * 30], L"X");
	wcscpy( &pwcharrFileNames[ 2* 30],  L"Y");
	wcscpy( &pwcharrFileNames[ 3 * 30], L"Z");
	wcscpy( &pwcharrFileNames[ 4 * 30], L"Psi");
	wcscpy( &pwcharrFileNames[ 5 * 30], L"Omega");
	wcscpy( &pwcharrFileNames[ 6 * 30], L"PI");
	wcscpy( &pwcharrFileNames[ 7 * 30], L"V");
	wcscpy( &pwcharrFileNames[ 8 * 30], L"Tetta");

	pscaleY = new double  [QUANT_COLS] ;
	pscaleY[1] = 0.1;
	pscaleY[2] = 0.1;
	pscaleY[3] = 0.1;
	pscaleY[4] = 1000.;
	pscaleY[5] = 1.;
	pscaleY[6] = 1.;
	pscaleY[7] = 1.;
	pscaleY[8] = 1000.;


	 }



  int i = 0 ;

  for ( i = 0; i < QUANT_POINTS; i++)
  {
	  fncEilerStep(Environment, VAlStepInt);

	   if (wcharrPath)
	   {
	   double *p = &parrBuff[i *  QUANT_COLS] ;
	   p[0] =  (double)mTCur ;
	   p[1] = (double)marrStrSK_VS [0];
	   p[2] = (double)marrStrSK_VS [1];
	   p[3] = (double)marrStrSK_VS [2];
	   p[4] = (double)marrStrSK_VS [3];
	   p[5] = (double)marrStrSK_VS [4];
	   p[6] = (double)marrStrSK_VS [5];
	   p[7] = (double)marrStrSK_VS [6];
	   p[8] = (double)marrStrSK_VS [7];
		}

  }

	 if (wcharrPath1)
	 {
	 for (int j = 1; j < QUANT_COLS -1; j++)
	 {


	TYrWriteShapeFile::WriteOneReport(wcharrPath  // ���� � �����
								  , parrBuff // ������ � ����������� - ������� nBuffRows x nBuffCols
								  ,QUANT_COLS // - �-�� ���������� � ������� ��������� ���������� � ������
								  ,QUANT_POINTS //  - �-�� �����
								  ,pwcharrFileNames //������� � �������� ���������� - ������� nBuffCols x lenName
								  ,lenName // ������������ ����� ����� ����������
								  ,0  // ����� ���������� �� ��� X
								  ,j  // ����� ���������� �� ��� Y
								  ,100 //  ������� �� ��� X
								  ,pscaleY[j]  // ������� �� ��� Y
								   );
	 }

	 TYrWriteShapeFile::WriteOneReport(wcharrPath  // ���� � �����
								  , parrBuff // ������ � ����������� - ������� nBuffRows x nBuffCols
								  ,QUANT_COLS // - �-�� ���������� � ������� ��������� ���������� � ������
								  ,QUANT_POINTS //  - �-�� �����
								  ,pwcharrFileNames //������� � �������� ���������� - ������� nBuffCols x lenName
								  ,lenName // ������������ ����� ����� ����������
								  ,1  // ����� ���������� �� ��� X
								  ,2  // ����� ���������� �� ��� Y
								  ,1.//  ������� �� ��� X
								  ,1.  // ������� �� ��� Y
								   );
   TYrWriteShapeFile::WriteOneReport(wcharrPath  // ���� � �����
								  , parrBuff // ������ � ����������� - ������� nBuffRows x nBuffCols
								  ,QUANT_COLS // - �-�� ���������� � ������� ��������� ���������� � ������
								  ,QUANT_POINTS //  - �-�� �����
								  ,pwcharrFileNames //������� � �������� ���������� - ������� nBuffCols x lenName
								  ,lenName // ������������ ����� ����� ����������
								  ,1  // ����� ���������� �� ��� X
								  ,3  // ����� ���������� �� ��� Y
								  ,1.//  ������� �� ��� X
								  ,1.  // ������� �� ��� Y
								   );
  wchar_t wchFileName4[300] = {0} ;
  wcscpy(wchFileName4, wcharrPath );
  wcscat(wchFileName4, L"\\Axes.shp");
  TYrWriteShapeFile::CreateShpAxes(wchFileName4,-40000.,40000,-40000.,40000) ;

  delete parrBuff ;
  delete pwcharrFileNames ;
  delete pscaleY ;
  */



}
//---------------------------------------------------------------------------



void __fastcall TForm3::Button9Click(TObject *Sender)
{
   OpenDialog1->Filter = L"����� � ��������� (*.shp)|*.shp" ;

	if (OpenDialog1->Execute())
	{
	mpwchOutFileTraj0 =  (OpenDialog1->FileName).w_str();

	}
	Edit4->Text =  mpwchOutFileTraj0;
}
//---------------------------------------------------------------------------
