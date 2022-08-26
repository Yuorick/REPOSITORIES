// ---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "FormSeaTargs.h"
// ---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"

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
#include "SeaTargPicturte.h"

#include "NeighbourhoodAppPoint.h"
#include "ProbabilityTheory.h"

static int I_ENTER_COUNT1 = 0;
extern const double VAL_C;
TForm4 *Form4;
extern TPlanePolygon PlanePlgDEBUG;
extern double arrDestrUZP [];
extern const int NumColsTab1_DestrUZP;
extern double arrCutterUZP [];
extern const int NumColsTab1_CutterUZP;
// ---------------------------------------------------------------------------
__fastcall TForm4::TForm4(TComponent* Owner) : TForm(Owner)
{
}
// ---------------------------------------------------------------------------

void __fastcall TForm4::fncInputData() {
	String s_22 = LabeledEdit27->Text;
 //	if (wcschr(s_22.w_str(), L',')) {

		// ����
		// ���� �������, ����

		mBearing0 = StrTo_Dbl_(LabeledEdit23->Text) * M_PI / 180.;

		// ���������, �

	 //	mDist0 = StrTo_Dbl_(ComboBox5->Items->Text); // ����
		mDist0 = StrTo_Dbl_(ComboBox5->Text);

		// ���� �����
		mTargCourse0 = StrTo_Dbl_(LabeledEdit29->Text) * M_PI / 180.;

		///

		double valAUDelayT = StrTo_Dbl_(LabeledEdit50->Text)/ 1000.;
		mSigAUDelayT = sqrt(TProbabilityTheory::calcDispRavnomern(valAUDelayT) );

		mSigDrivAY_U = StrTo_Dbl_(LabeledEdit52->Text) / 1000.;
		mSigDrivAY_dU_po_dt = StrTo_Dbl_(LabeledEdit53->Text) / 1000.;

 //	}
 //	else {
		// ���������������
		// mWidthDgr =  StrTo_Dbl_(LabeledEdit3->Text) * M_PI/3000./2. ;
 //	}

	// ������, �
	mElev0 = 0.;
	//
	mTargZenitAng0 = M_PI / 2.;

	switch(ComboBox3->ItemIndex) // ����
	{
	case 0:
		mEnumTargType = DESTROYER;
		mVelocity0 = 10.;
		mTargEPR = 1000.;
		mWSkz = 0.002;
		break;
	case 1:
		mEnumTargType = MOTORBOAT;
		mVelocity0 = 20.;
		mTargEPR = 500.;
		mWSkz = 0.01;
		break;

	default:

		break;

	}
	LabeledEdit25->Text = int(mVelocity0);
	LabeledEdit28->Text = double(int(mTargEPR * 1000.)) / 1000.;
	LabeledEdit27->Text = double(int(mWSkz * 1000.)) / 1000.;
	// LabeledEdit22->Text = mElev0 ;
	// ��. ���������� ShellBody (mEnumShellType). ShellBody ���������� ����� ��������� �����������
	switch(ComboBox1->ItemIndex) {
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
	default:
		mEnumShellType = CALIBRO_UNKNOWN;
		menumCannonType = CANNON_UNKNOWN;
		break;
	}

	///
	// ��� ����������
			// ��� ����������
		switch(ComboBox2->ItemIndex)
		{
			case 0:
			mDetonatorType = D4MRM;
			break;

			case 1:
			mDetonatorType = MFIVU;
			break;

			default:
			mDetonatorType = DETONATOR_UNKNOWN;
			break;
		}


		///


	TDetonator Detonator(mDetonatorType);
	mShellBody = TShellBody(mEnumShellType, Detonator);
	mArtCannon = TArtCannon(menumCannonType, mSigAUDelayT);
	mArtComplex = TArtComplex(mArtCannon, mSigDrivAY_U, mSigDrivAY_dU_po_dt);

	LabeledEdit49->Text = mArtCannon.mRateOfFire;

	mQuantShells = StrToInt(LabeledEdit10->Text);

	mInitTargData = TInitTargData(mBearing0, mTargCourse0, mTargZenitAng0,
		mVelocity0, mDist0, mElev0, 0.);
	const double TCur0 = 0.;
	double arrWSkz[3] = {
		0.
	};
	arrWSkz[0] = mWSkz;
	arrWSkz[1] = mWSkz;
	arrWSkz[2] = 0.;
	TTraject Traj0(TCur0, arrWSkz, mInitTargData);
	TTarget Target(Traj0, mEnumTargType, mTargEPR, NULL);
	mVessel = TVessel(mShellBody, mFar_2D, mTransmitAnt, mDriverSigBet,
		mDriverSigEps, mDriverDynamicSigBet, mDriverDynamicSigEps, mMaxSig_Q,
		mMaxSig_Psi, mMaxSig_Tet, mMaxSig_dQdt, mMaxSig_dPsidt, mMaxSig_dTetdt
		// ����
		, mMaxSig_H, mMaxSig_VH, mK1, mSigV, mEnvironment, mVesselWidth,
		mVesselLength, marrFarParallacs, mMaxQ, mT_Q, mMaxPsi, mT_Psi, mMaxTet,
		mT_Tet, mMaxVert, mQ0, mVVess, mInitTargData, mMaxAmp_AftFlexure,
		mT_AftFlexure, mMaxAmp_BoardFlexure, mT_BoardFlexure, mControlSyst, marrArtParral,
		mArtComplex, NULL);
	mFight = TFight(mVessel, Target, mVessel.mControlSyst.mFiltT, mEtalonSign,
		mEnvironment, NULL);

  //	TFar Far(mVessel.mFar_2D, true);
  //	const double VAlTetta0 = findDiagrWidth
   //			(Far.mFaceta.m_d, Far.m_D, Far.mFaceta.m_n, Far.m_N, mVessel.mFar_2D.mLambda);
   //	int ib = (VAlTetta0 * 100000.);
   //	LabeledEdit34->Text = ((double)ib) / 100.;

	LabeledEdit51->Text = double
			(int(mVessel.mArtComplex.mArtCannon.mAngGroupedFire * 100000.)) / 100.;

}

void __fastcall TForm4::Button4Click(TObject *Sender) {

	OpenDialog1->Filter = L"����� � ��������� (*.shp)|*.shp";

	if (OpenDialog1->Execute()) {
		mpwchOutFile0 = (OpenDialog1->FileName).w_str();

	}
	Edit2->Text = mpwchOutFile0;
}
// ---------------------------------------------------------------------------

void __fastcall TForm4::Button1Click(TObject *Sender) {
	Close();
}
// ---------------------------------------------------------------------------

void __fastcall TForm4::Button3Click(TObject *Sender)
{

//double valOmega = TNeighbourhoodAppPoint::calcOmega(&arrCutterUZP[NumColsTab1_CutterUZP], NumColsTab1_CutterUZP);
/*
 PlanePlgDEBUG;
	 double arrTargV [3] = {0.}, arrVMiss[3] = {0.};
	 double valTargCourse = 225. * M_PI / 180.;
	 arrTargV [0] = 10. * sin(valTargCourse);
	 arrTargV [1] = 10. * cos(valTargCourse);

	 double valTet = -45. * M_PI /180.;
	 double valBet =   45. * M_PI / 180.;
	 double valV0 = 300.;
	 arrVMiss[0] = valV0 * cos (valTet) * sin ( valBet);
	 arrVMiss[1] = valV0 * cos (valTet) * cos ( valBet);
	 arrVMiss[2] = valV0 * sin (valTet) ;



	TURPolygon 	pPllShadow =  PlanePlgDEBUG.createShadowPlg(arrTargV,  arrVMiss);
	pPllShadow.WriteSetSHPFiles(L"E:\\PROJECTS_C++\\TARAN\\New\\tempShadow3.shp", &pPllShadow, 1) ;
	 return;    */
	// �������� ���� � ����� � �������
	if (!mpwchOutFile0) {
		if (I_ENTER_COUNT1 > 4)
			Close();

		ShowMessage(L"������� ���� � ����� � ���������");
		I_ENTER_COUNT1++;
		return;
	}
	wcscpy(mwchOutFold, mpwchOutFile0);
	wchar_t *pwchr = wcsrchr(mwchOutFold, L'\\');
	pwchr[0] = 0;
	String wchFoldName = mwchOutFold;

	fncInputData();

	double  valProb = -1.0, valProb0 = -1.0;
	double valSKZPromach = -1.0, valHAntenna = 20.;
	double arrCoMtrx_1_and_3_Groups [36] = {0.} // ������ ������� ������� ������� ������������� �������� 1 � 3 �����
	  , arrCoMtrx_2_Group [36] = {0.}; // ������ ������� ������� ������� ������������� �������� 1 � 3 �����
	double valKGSKEps = -1.0, valKGSKBet = -1.0, valTFlight = -1.0; // ����, ����, �������� �����
	double arrCorrMatrxVectMiss_GSK [36] = {0.}; // ������� ������� �������� �������
	double arrVectMiss_GSK[6] = {0.}; // ������� ������ ������� � ����� �������

	if (!mFight.calcSuccessProbSeaTarg(mQuantShells, valHAntenna, &valSKZPromach
	,arrCoMtrx_1_and_3_Groups, arrCoMtrx_2_Group, &valKGSKEps , &valKGSKBet, &valTFlight
	, arrVectMiss_GSK, &valProb, &valProb0))
	{

	}

	MtrxSumMatrx(arrCoMtrx_1_and_3_Groups, arrCoMtrx_2_Group, 6, 6, arrCorrMatrxVectMiss_GSK);
	///
	int ia = (valProb * 100);
	LabeledEdit46->Text = ((double)ia) / 100.;

	ia = (valProb0 * 1000);
	LabeledEdit34->Text = ((double)ia) / 1000.;

	createNibourAppointmPointPictureForSeaTarg(wchFoldName.w_str(),arrCoMtrx_1_and_3_Groups, arrCoMtrx_2_Group,arrVectMiss_GSK
	, mFight.mTarget, mQuantShells) ;


}
// ---------------------------------------------------------------------------

void __fastcall TForm4::Button2Click(TObject *Sender) {

	if (!mpwchOutDataFileTaran)
	{
		if (I_ENTER_COUNT1 > 4)
			Close();

		ShowMessage(L"������� ���� � ����� � ���������");
		I_ENTER_COUNT1++;
		return;
	}

	double arrCoMtrx_1_and_3_Groups[36] = {	0.}, arrCoMtrx_2_Group[36] = {0.}
				,arrCorMtrx[36] = {0.}, arrCorMtrxRez[4] = {0.};
	int iNumCols = 3;
	int iNumRows = ComboBox5->GetCount();
	int iLenName = 30;
	double *parrBuff = new double[iNumCols * iNumRows];

	wchar_t *pwcharrColNames = new wchar_t[iNumCols * iLenName];
	memset(pwcharrColNames, 0, iNumCols * iLenName*sizeof(wchar_t));

	wcscpy(pwcharrColNames, L"DIST");
	wcscpy(&pwcharrColNames[iLenName], L"FIRST_SEMIAXIS");
	wcscpy(&pwcharrColNames[2 * iLenName], L"SECOND_SEMIAXIS");

	for (int i = 0; i < iNumRows; i++)
 //	for (int i = 2; i < 3; i++)
	 {
		ComboBox5->ItemIndex = i;
		fncInputData();
		mFight.mVessel.mMaxQ = 0.; // ������� ��������� �� ��������� ������
		memset( mFight.mTarget.mTraject.marrSigW, 0, 3 * sizeof(double));  // ���� ��������� �� ��������� ������
		// ������� ������ � ����� �������

			// 1. ������� ������ � ����� �������
  // ���������� ������� ��������� �� � ����
	double arrPositionAY_KGSK[3] = {0.};
	mFight.mVessel.calcAY_Position(arrPositionAY_KGSK);
	double valKGSKEps = 0., valKGSKBet  = 0.;
	double arrVectAppointmentPointGSK[6] = {0.};
	double valMiss = -1., valTFlight = -1.;

	// ���������� ��������� ������� ��������� ���� � ���� �� �������� ��������
	double arrTargVS_KGSK0[6] ={0.};
	MtrxMinusMatrx(mFight.mTarget.mTraject.marrVectSostGSK_Begin, mFight.mVessel.marrVectSost,1, 6, arrTargVS_KGSK0);
	///
	if (!mFight.calcAppointmentPoint(NULL, NULL, &(mFight.mVessel.marrVectSost[3])
	,arrTargVS_KGSK0, arrPositionAY_KGSK
	, &valKGSKEps, &valKGSKBet, &valTFlight, arrVectAppointmentPointGSK, &valMiss))
	{
		ShowMessage(L"ERROR_ TFight::calcSuccessProbCoast");
		return ;
	}
///
 /////////////////////////////////////////////////////////////////////////
	//////////////  ��� ��� �������!!! ������ ����������������� � ��� ������� !!!!///////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////

//valKGSKEps =	0.35873669492725 ;
//valKGSKBet =	6.27370551180827 ;
//valTFlight =	43.3708000003162 ;
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////


    double arrShellVS_GSK [6] = {0.};
		mFight.calcCorMtrx_First_And_Third_Group(true,valKGSKEps, valKGSKBet, valTFlight,arrCoMtrx_1_and_3_Groups, arrShellVS_GSK);
		const double VAlObservTime = 10.;
		mFight.calcCorMtrx_Second_Group(VAlObservTime, valTFlight, arrCoMtrx_2_Group);
		MtrxSumMatrx(arrCoMtrx_1_and_3_Groups, arrCoMtrx_2_Group, 6, 6, arrCorMtrx);
		arrCorMtrxRez[0] = arrCorMtrx[0];
		arrCorMtrxRez[1] = arrCorMtrx[1];
		arrCorMtrxRez[2] = arrCorMtrx[1];
		arrCorMtrxRez[3] = arrCorMtrx[7];

		double arrV[4] = {0.}, arrLamb[4] = {0.};
		CalcProperVectors2(arrCorMtrxRez, arrV, arrLamb);
		parrBuff[i * iNumCols] = mDist0;
		parrBuff[i * iNumCols + 1] = sqrt(arrLamb[0]);
		parrBuff[i * iNumCols + 2] = sqrt(arrLamb[3]);

	}
	TYrWrite::WriteMassiveInFIleSCV(mpwchOutDataFileTaran, parrBuff, iNumRows, iNumCols,
		NULL, pwcharrColNames, iLenName);
	delete parrBuff;
	delete pwcharrColNames;
}

// ---------------------------------------------------------------------------

void __fastcall TForm4::Button5Click(TObject *Sender)
{
 SaveDialog1->Filter = L"TXT ����� (*.csv)|*.csv" ;


	if (SaveDialog1->Execute())
	{
	 //	ShowMessage( (SaveDialog1->FileName).w_str()) ;
	 mpwchOutDataFileTaran =  (SaveDialog1->FileName).w_str();

	}
	Edit1->Text =mpwchOutDataFileTaran;

}
//---------------------------------------------------------------------------

// ---------------------------------------------------------------------------

void __fastcall TForm4::create5P10()
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
	 //double arrDelt[4] = {0.};


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
void __fastcall TForm4::Button6Click(TObject *Sender)
{

	if (mpwchOutFile0 == NULL) {
		if (I_ENTER_COUNT1 > 4)
			Close();

		ShowMessage(L"������� ���� � ����� � ���������");
		I_ENTER_COUNT1++;
		return;
	}
	wcscpy(mwchOutFold, mpwchOutFile0);
	wchar_t *pwchr = wcsrchr(mwchOutFold, L'\\');
	pwchr[0] = 0;
	String wchFoldName = mwchOutFold;




	int iNumCols = 2;


	int iNumRows = 90./ 5. + 1;
	int iLenName = 30;
	double *parrBuff = new double[iNumCols * iNumRows];

	wchar_t *pwcharrColNames = new wchar_t[iNumCols * iLenName];
	memset(pwcharrColNames, 0, iNumCols * iLenName*sizeof(wchar_t));

	wcscpy(pwcharrColNames, L"TargCourse");
	wcscpy(&pwcharrColNames[iLenName], L"Probability");

	double pValSum  = 0.;
	for (int i = 0; i < iNumRows; i++)
 //	for (int i = 2; i < 3; i++)
	 {
	 LabeledEdit29->Text = 180. +((double)i) * 5.;
	 Button3Click(Sender);
	 double valPCur =  StrTo_Dbl_(LabeledEdit46->Text) ;
	 pValSum   += valPCur;

		parrBuff[i * iNumCols] = mTargCourse0 * 180./ M_PI;
		parrBuff[i * iNumCols + 1] = valPCur * 100.;

	}

	// ���  ���������
	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  wchFoldName.w_str());
	wcscat(wchAxesFileName, L"AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-500., 500.
	,-500., 500.,30.) ;
	// ������

	TYrWriteShapeFile::WriteOneReport(wchFoldName.w_str() // ���� � �����
									,parrBuff // ������ � ����������� - ������� nBuffRows x nBuffCols
									,iNumCols  // - �-�� ���������� � ������� ��������� ���������� � ������
									,iNumRows //  - �-�� �����
									,pwcharrColNames //������� � �������� ���������� - ������� nBuffCols x lenName
									,30 // ������������ ����� ����� ����������
									,0  // ����� ���������� �� ��� Y
									,1  // ����� ���������� �� ��� X
									,1  //  ������� �� ��� Y
								  ,1  // ������� �� ��� X
									 ) ;

	delete parrBuff;
	delete pwcharrColNames;



	int ia = (pValSum / ((double)iNumRows) * 100.);
	LabeledEdit46->Text = ((double)ia) / 100.;
}
//---------------------------------------------------------------------------

