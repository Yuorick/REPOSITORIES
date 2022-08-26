//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "MainForm.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"

#include <math.h>
#include <dir.h>
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
#include "ProbabilityTheory.h"
#include "Circle.h"
#include "NeighbourhoodAppPoint.h"

#include "LinOptimization.h"
#include "GameTheory.h"

#include "Table_1D.h"
#include "Table_2D.h"
#include "Table_3D.h"
#include "DiagrSinX.h"


static int I_ENTER_COUNT = 0;
extern const double DTi;
extern const double tb1;
extern const double VAL_C;

TForm1 *Form1;

_fastcall TForm1::TForm1(TComponent* Owner) : TForm(Owner)
{
	create5P10();


	fillVessel();
	fncInputData();



}

// ---------------------------------------------------------------------------

void __fastcall TForm1::fncInputData()
{


	//if (wcschr(s_22.w_str(), L','))
 //	{

		// ����� �����
		mLambda = VAL_C / StrTo_Dbl_(LabeledEdit15->Text) * 100. / 1000000000.;
		LabeledEdit13->Text = mLambda;

		// ���������� ����� ������������
		mdEmitCol = mLambda * 0.55;
		mdEmitRow = mLambda * 0.55;
		LabeledEdit12->Text = mdEmitCol;

		// ��������� ������ ��� ������� ���������
		// ���������
		mEtalonAmp = StrTo_Dbl_(LabeledEdit4->Text);
		// ���������
		mEtalonDist = StrTo_Dbl_(LabeledEdit10->Text) * 1000.;
		// ���
		mEtalonAPR = StrTo_Dbl_(LabeledEdit16->Text);
		// ��� ����� ���� ��������� ��������� 5�10
		mNoiseSKZ_5P10 = StrTo_Dbl_(LabeledEdit17->Text);
		// ��� �������� ������� �������� ��������� ��������� 5�10
		mEtalonSigAmplFact_5P10 = StrTo_Dbl_(LabeledEdit18->Text);

		// �������� �� ��������
		mEtalonPowerPrd = StrTo_Dbl_(LabeledEdit36->Text);
		// �� �� ��������
		mEtalonKYPrd = StrTo_Dbl_(LabeledEdit38->Text);
		// �������� �� �����
		mEtalonKYPriem = StrTo_Dbl_(LabeledEdit37->Text);


		mEtalonSign = TEtalonSign(mEtalonAmp, mEtalonDist, mEtalonAPR,
			mNoiseSKZ_5P10, mEtalonSigAmplFact_5P10, mEtalonPowerPrd, mEtalonKYPrd
          , mEtalonKYPriem);
		///

		//  ������� �� ��������
        	// �������� �� ��������
		mPowerPrd = StrTo_Dbl_(LabeledEdit40->Text);
		// �� �� ��������
		 mKYPrd = StrTo_Dbl_(LabeledEdit39->Text);

		 mTransmitAnt.mKYPrd  =mKYPrd;
		 mTransmitAnt.mPowerPrd = mPowerPrd;
		///


		///








 //	}
 //	else
 //	{
		// ���������������
		// mWidthDgr =  StrTo_Dbl_(LabeledEdit3->Text) * M_PI/3000./2. ;
 //	}



	// �-�� ��������  �����������

	mNumEmitCols = StrToInt(LabeledEdit9->Text);
	// �-�� ����� �����������
	mNumEmitRows = StrToInt(LabeledEdit11->Text);
	// �-�� �������� ��  (�� ����������� )
	mNumAMCols = StrToInt(LabeledEdit5->Text);
	// �-�� ����� �� (�� ��������� )
	mNumAMRows = StrToInt(LabeledEdit6->Text);
	// ���������� ����� �� �� ���������
	mdAMRow = ((double)mNumEmitRows) * mdEmitRow;
	// ���������� ����� �� �� �����������
	mdAMCol = ((double)mNumEmitCols) * mdEmitCol;
	LabeledEdit8->Text = mdAMRow;
	LabeledEdit7->Text = mdAMCol;
	///

	const double SigEmitNoise = mNoiseSKZ_5P10 / sqrt(((double)(8 * 16 * 28)));
	const double SigEmitAmplFact = mEtalonSigAmplFact_5P10 / sqrt
		(((double)(8 * 16 * 28)));
	mAM_2D = TAM_2D(mNumEmitCols, mNumEmitRows, mdEmitCol, mdEmitRow,
		SigEmitNoise, SigEmitAmplFact);


	mFar_2D = TFar_2D(mNumAMCols, mNumAMRows, mLambda, mdAMCol, mdAMRow,
		mAM_2D, mbarrAM, mSigmaR);

	for (int i = 0; i < mNumAMCols * mNumAMRows; i++) {

		if (!mbarrAM[i]) {
			mFar_2D.mpAm2D[i].mOtklCoefUs = 1.;
			mFar_2D.mpAm2D[i].mSigEmitNoise = 0.;
		}
	}

	mAntpCoef = StrTo_Dbl_(LabeledEdit29->Text);

	// ������ ����������
	 mTrajectH = StrTo_Dbl_(LabeledEdit28->Text);

	// ��� ����
	  mEPR = StrTo_Dbl_(LabeledEdit1->Text);

	  // ������ �������
	  mAntH =  StrTo_Dbl_(LabeledEdit14->Text) ;

	  //
	  LabeledEdit21->Text = mLambda/ mFar_2D.calcAppertVert();
}

//---------------------------------------------------------------------------

void __fastcall TForm1::Button3Click(TObject *Sender)
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


	 fncInputData() ;

	 bool bLat = false;
	 double valSigE = -1.,valSigQ = -1.;
	  double valDistzahv = mFar_2D.calc_TwoTargsZahvatDist(mAntpCoef,mTrajectH, mEPR
	   ,mPowerPrd, mKYPrd, mEtalonSign , mAntH
	   ,&bLat  , &valSigE ,&valSigQ  ) ;

	   int ia = (valDistzahv * 10.) ;
	   LabeledEdit2->Text = ((double)ia)/ 10.;

	   ia = (valSigE * 1000. * 100.) ;
	   LabeledEdit3->Text = ((double)ia)/ 100.;

	   ia = (valSigQ * 1000. * 100.) ;
	   LabeledEdit19->Text = ((double)ia)/ 100.;

	 const int NCols = 21;

	 double  valDMax = 6000.; // ���� ���������
	 double  valDMin = 500.; // ���� ���������
	  double step = 100.;
	 const int NRows = (valDMax - valDMin)/ step;
	 double *arrBuff  = new double [NRows* NCols];

      TFar Far0(mFar_2D, true);
	// TFar Far(Far0, 4) ;
		int quantFalseSign = 1;
	  // �������� �������� ���� ���� � ���������
	  double valBetTargTrue = 0., valBetAntpTrue =0.;

	 TComp *cmparrTrueAmMeasures = new TComp[ mFar_2D.mNumAMCols * mFar_2D.mNumAMRows];
	 TComp *cmparrNoisedAmMeasures = new TComp [ mFar_2D.mNumAMCols * mFar_2D.mNumAMRows];



	 double arrDisp[4]  = {0.};
	 arrDisp[0]=arrDisp[1]=arrDisp[2]=arrDisp[3]=  mNoiseSKZ_5P10* mNoiseSKZ_5P10/4.;
     TFaceta Faceta(mNumEmitRows,mdEmitRow, mLambda) ;
	 TFar Far1( 4, mdAMRow, mLambda  , Faceta, arrDisp) ;
     double valEstEpsTarg = 0.,valEstEpsAntp = 0., arrMtrxCorrEps [4]= {0.};
	 TComp cmpKTarg , cmpKAntp;
	 for (int i =0; i < NRows; i++)
	 {
		double valHorDist =  valDMax - step * ((double)i);
		double valDist = sqrt( valHorDist * valHorDist + mTrajectH * mTrajectH);
		// ���������� ��������� ���� ����� ���� � ��������
		double valEpsTargTrue = atan((mTrajectH - mAntH)/  valHorDist);
		double valEpsAntpTrue = -atan((mTrajectH + mAntH)/  valHorDist);
		///

		// �������� ������� ����
		double valAmpTarg = mFar_2D.calcCurrentSignalAmpl (valDist, mEPR
		,mPowerPrd, mKYPrd, mEtalonSign) ;
		///

		// �������� ������� ��������
		double valAmpAntp =  valAmpTarg  *mAntpCoef;
		///

		// �������� ������� ���� � �������� ��������

		double valTargSignPhase = 2. * M_PI * getRand01( );
	   //	TSingleSign SignTarg(valBetTargTrue, valEpsTargTrue,  valAmpTarg, valTargSignPhase ) ;

		double valAntpSignPhase = 2. * M_PI * getRand01( );
	   //	TSingleSign SignAntp(valBetAntpTrue, valEpsAntpTrue, valAmpAntp, valAntpSignPhase) ;
		///


         	// ���������� ������� ��������� ���������
		//  ���������� �������� ��������  ��������
		////////////////////////////////////////////////////////////
		TComp arrFourRowMeasures[4];
		memset(arrFourRowMeasures, 0, 4 * sizeof(TComp));







	TComp pcmpS[4];


	TComp cmpTrueKTarg (valAmpTarg * cos(valTargSignPhase), valAmpTarg * sin(valTargSignPhase) );
	TComp cmpTrueKAntp (valAmpAntp * cos(valAntpSignPhase), valAmpAntp * sin(valAntpSignPhase) );

	Far1.ImitateMeasureArray(valEpsTargTrue,cmpTrueKTarg, valEpsAntpTrue,cmpTrueKAntp, pcmpS,arrFourRowMeasures);

		// ����� ��������� ��������
		Far1.fncEstimateMsd(arrFourRowMeasures, &valEstEpsTarg, &valEstEpsAntp
		, &cmpKTarg , &cmpKAntp , arrMtrxCorrEps );
		double valEstHMSD = mAntH + sin(valEstEpsTarg) * valDist;
		///


		// ������ ���
		double valRSMDisp = -1.;
		double valEstEpsRSM = Far1.fncEstimateRsmTetta(arrFourRowMeasures,  &valRSMDisp ) ;
		double valEstHRSM = mAntH + sin(valEstEpsRSM) * valDist;

		// ���������� ������
		arrBuff[i *  NCols] =  valHorDist*0.1;
		for (int j = 0; j < 4; j++)
		{
		  arrBuff[i *  NCols + 1 + j *2] = arrFourRowMeasures[j].modul();  // ����
		  arrBuff[i *  NCols + 1 + j *2 +1] = arrFourRowMeasures[j].phase() * 180./ M_PI; // ����
		}
		arrBuff[i *  NCols + 9] = valAmpTarg/ sqrt(Far1.calcNoiseDisp()) ;   // ������/ ���
		arrBuff[i *  NCols + 10] =  valEstEpsRSM * 1000.;// ��� ����������
		arrBuff[i *  NCols + 11] =  sqrt(valRSMDisp) * 1000.;  // ��� ���
		arrBuff[i *  NCols + 12] = valEstHRSM ; // ������ ������ ���

		arrBuff[i *  NCols + 13] =  valEstEpsTarg * 1000.; // ������ �� ���� ���
		arrBuff[i *  NCols + 14] =  valEstEpsAntp * 1000.; // ������ �� ���� ���
		arrBuff[i *  NCols + 15] =  sqrt(arrMtrxCorrEps[0]) * 1000.; // ��� ������ ���� ���� ���
		arrBuff[i *  NCols + 16] =  sqrt(arrMtrxCorrEps[3]) * 1000.; // ��� ������ ���� ���� ���
		arrBuff[i *  NCols + 17] =   valEstHMSD ; // ������ ������ ���
		arrBuff[i *  NCols + 18] =  valEpsTargTrue * 1000.; // �������� �� ����
		arrBuff[i *  NCols + 19] =  valEpsAntpTrue * 1000.; // �������� �� ����
		arrBuff[i *  NCols + 20] =  mTrajectH ; // �������� ������
	 }

	 wchar_t *wcharrFileNames = new wchar_t[NCols * 30];
	 memset(wcharrFileNames, 0,NCols * 30* sizeof (wchar_t));
	 wcscpy(&wcharrFileNames[0],L"D");
	 wcscpy(&wcharrFileNames[30],L"A0");
	 wcscpy(&wcharrFileNames[30 * 2],L"Fi0");
	 wcscpy(&wcharrFileNames[30 * 3],L"A1");
	 wcscpy(&wcharrFileNames[30 * 4],L"Fi1");
	 wcscpy(&wcharrFileNames[30 * 5],L"A2");
	 wcscpy(&wcharrFileNames[30 * 6],L"Fi2");
	 wcscpy(&wcharrFileNames[30 * 7],L"A3");
	 wcscpy(&wcharrFileNames[30 * 8],L"Fi3");
	 wcscpy(&wcharrFileNames[30 * 9],L"SigNoise");
	 wcscpy(&wcharrFileNames[30 * 10],L"EstEps_RSM");
	 wcscpy(&wcharrFileNames[30 * 11],L"SKZ_RSM");
	 wcscpy(&wcharrFileNames[30 * 12],L"EstH_RSM");
	 wcscpy(&wcharrFileNames[30 * 13],L"EstTargEps_MSD");
	 wcscpy(&wcharrFileNames[30 * 14],L"EstAntpEps_MSD");
	 wcscpy(&wcharrFileNames[30 * 15],L"SKZ_TargEps_MSD");
	 wcscpy(&wcharrFileNames[30 * 16],L"SKZ_AntpEps_MSD");
	 wcscpy(&wcharrFileNames[30 * 17],L"EstH_MSD");
	 wcscpy(&wcharrFileNames[30 * 18],L"EpsTarg_True");
	 wcscpy(&wcharrFileNames[30 * 19],L"EpsAntp_True");
	 wcscpy(&wcharrFileNames[30 * 20],L"H_True");


	 double scalex = 1.;
	 double scaley = 1.;
		wchar_t wchFoldName[300] ={0};
	wcscpy(  wchFoldName,  mwchOutFold);
	wcscat(wchFoldName, L"\\");
	 for (int i = 1; i < NCols; i++)
	 {
	  TYrWriteShapeFile::WriteOneReport(wchFoldName  // ���� � �����
								  , arrBuff // ������ � ����������� - ������� nBuffRows x nBuffCols
								  ,NCols // - �-�� ���������� � ������� ��������� ���������� � ������
								  ,NRows //  - �-�� �����
								  ,wcharrFileNames //������� � �������� ���������� - ������� nBuffCols x lenName
								  ,30 // ������������ ����� ����� ����������
								  ,0  // ����� ���������� �� ��� X
								  ,i   // ����� ���������� �� ��� Y
								  , scalex  //  ������� �� ��� X
								  , scaley  // ������� �� ��� Y
								   ) ;
	 }

free(arrBuff);

	///


	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  wchFoldName);
	wcscat(wchAxesFileName, L"AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,0., 10000.
	 ,0.,10000., 2.) ;
     delete [] cmparrTrueAmMeasures;
	 delete [] cmparrNoisedAmMeasures ;


	

}
//-----------------------------------------------------------------------

//---------------------------------------------------------------------------

void __fastcall TForm1::Button4Click(TObject *Sender)
{
OpenDialog1->Filter = L"����� � ��������� (*.shp)|*.shp" ;

	if (OpenDialog1->Execute())
	{
	mpwchOutFile0 =  (OpenDialog1->FileName).w_str();

	}
	Edit2->Text =  mpwchOutFile0;
}
//---------------------------------------------------------------------------


void __fastcall TForm1::Panel7Click(TObject *Sender)
{
fncInputData()  ;



}
//---------------------------------------------------------------------------

void __fastcall TForm1::Panel2Click(TObject *Sender)
{
fncInputData();
}

//---------------------------------------------------------------------------

void __fastcall TForm1::Panel6Click(TObject *Sender)
{
fncInputData();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::Panel3Click(TObject *Sender)
{
fncInputData();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::Panel8Click(TObject *Sender)
{
fncInputData();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::Panel9Click(TObject *Sender)
{
fncInputData();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::Panel10Click(TObject *Sender)
{
fncInputData()  ;

}
//---------------------------------------------------------------------------



//---------------------------------------------------------------------------


void __fastcall TForm1::Button1Click(TObject *Sender)
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


	 fncInputData() ;
	 TFar Far0(mFar_2D, true);
	 TFar Far(Far0, 2) ;

	 double valDiagrWirth = Far.findDiagrWidthApprox();
	 const double VAlNAppert =  mFar_2D.calcAppertVert();

	 	// ���������� ��������� ���� ����� ���� � ��������
	   //	double valEpsTargTrue = 0.;//mLambda/ VAlNAppert-0.005  ;//0.026;//0.026;
		double valEpsTargTrue = StrTo_Dbl_(LabeledEdit20->Text) /1000.;


		double valFDiagr = Far.fncFFarApprox(valEpsTargTrue);

		double valFacetaDiagr = Far.mFaceta.fncFFacetaApprox (valEpsTargTrue);
		double valEpsAntpTrue = 0.;
		///

		// �������� ������� ����
		double valAmpTarg = 500.;//
		 //mFar_2D.calcCurrentSignalAmpl (valDist, mEPR,mPowerPrd, mKYPrd, mEtalonSign) ;
		///

		// �������� ������� ��������
		double valAmpAntp =  0.;
		///

	   double valDiagrNull = mLambda/ mFar_2D.calcAppertVert();
	   TComp cmpTemp = mFar_2D.fncIdealDiagr (0.,valDiagrNull);
	   double valTemp = Far.fncFFarApprox (valDiagrNull);




		 const int NCols = 7;


	 const int NRows = 10000;
	 double *arrBuff  = new double [NRows* NCols];

	 double valSKZ = 0.;



	 double valEstEpsTarg = 0. ,valEstEpsAntp = 0.;
	 TComp cmpKTarg , cmpKAntp;
     double valSumDisp = mFar_2D.calcSumDisp();
	 const double VAlNWaveCur = sqrt(mFar_2D.calcSumDisp())/ valAmpTarg;


	 double dispCalc = -1.;
	 for (int i =0; i < NRows; i++)
	 {
		// �������� ������� ���� � �������� ��������

		double valTargSignPhase = 2. * M_PI * getRand01( );
		double valAntpSignPhase = 2. * M_PI * getRand01( );

        		// ���������� ������� ��������� ���������
		//  ���������� �������� ��������  ��������
		TComp cmpSZv[2], cmpS[2];


	TComp cmpTrueKTarg (valAmpTarg * cos(valTargSignPhase), valAmpTarg * sin(valTargSignPhase) );
	TComp cmpTrueKAntp (valAmpAntp * cos(valAntpSignPhase), valAmpAntp * sin(valAntpSignPhase) );

		 Far.ImitateMeasureArray(valEpsTargTrue,cmpTrueKTarg, valEpsAntpTrue,cmpTrueKAntp
	 , cmpS,cmpSZv);


	 double valDisp = -1.;
	 double valEstEps = Far.fncEstimateRsmTetta(cmpSZv,  &valDisp);

	 valSKZ += 1000. *(valEstEps - valEpsTargTrue) * 1000. *(valEstEps - valEpsTargTrue);
	 /////



	 double valTheorDisp = TFar::calcTheorDisp_RSM(VAlNWaveCur, VAlNAppert
	  , mLambda, valEpsTargTrue);


	 // ����� �� �������������
	 double del = 1.;
	 double deriv[4] = {0.};
	 double valDisp1 = 0.;
	 TComp cmpSZv1[2];
	 double sum = 0.;

	 memcpy(cmpSZv1, cmpS, 2 * sizeof(TComp));
	 cmpSZv1[0].m_Re += del;
	 deriv[0] = (Far.fncEstimateRsmTetta(cmpSZv1,  &valDisp1) - valEpsTargTrue)/ del;
	 sum += deriv[0] * deriv[0];

	 memcpy(cmpSZv1, cmpS, 2 * sizeof(TComp));
	 cmpSZv1[0].m_Im += del;
	 deriv[1] = (Far.fncEstimateRsmTetta(cmpSZv1,  &valDisp1) - valEpsTargTrue)/ del;
	 sum += deriv[1] * deriv[1];

	 memcpy(cmpSZv1, cmpS, 2 * sizeof(TComp));
	 cmpSZv1[1].m_Re += del;
	 deriv[2] = (Far.fncEstimateRsmTetta(cmpSZv1,  &valDisp1) - valEpsTargTrue)/ del;
	 sum += deriv[2] * deriv[2];

	 memcpy(cmpSZv1, cmpS, 2 * sizeof(TComp));
	 cmpSZv1[1].m_Im += del;
	 deriv[3] = (Far.fncEstimateRsmTetta(cmpSZv1,  &valDisp1) - valEpsTargTrue)/ del;
	 sum += deriv[3] * deriv[3];

	 dispCalc =  valSumDisp * sum / 4.;






		// ���������� ������
		arrBuff[i *  NCols] =  ((double)i);
		arrBuff[i *  NCols + 1] =  valEpsTargTrue* 1000.;
		arrBuff[i *  NCols + 2] = valEstEps * 1000.;
		arrBuff[i *  NCols + 3] = (valEpsTargTrue - valEstEps) * 1000.;
		arrBuff[i *  NCols + 4] = 3. *sqrt(valTheorDisp) * 1000.;
		arrBuff[i *  NCols + 5] = -3. *sqrt(valTheorDisp) * 1000.;
		arrBuff[i *  NCols + 6] = sqrt(valDisp) * 1000.;



		TComp pcmpAmMeasures[64], pcmpAmMeasuresNoised[64];
		TSingleSign SingleSign(0., valEpsTargTrue, valAmpTarg, valTargSignPhase);
		mFar_2D.fncImitateMeasuresArray (SingleSign, pcmpAmMeasures, pcmpAmMeasuresNoised);

	 double valEstBet0, valSigBet0, valEstEps0, valSigEps0;
		mFar_2D.calcEstRSM(pcmpAmMeasuresNoised, &valEstBet0, &valSigBet0
	 , &valEstEps0, &valSigEps0) ;
	 int iii=0;

	 }

	 valSKZ = sqrt(valSKZ / ((double)NRows));


	 double valNAppert = VAlNAppert;
	 double valDispNew = calcDispNew(valEpsTargTrue, valSumDisp, valAmpTarg, mLambda,valNAppert/ 2.
	 ,valFacetaDiagr);
	 double temp = sqrt(valDispNew);


	 double valTheorDisp = TFar::calcTheorDisp_RSM(VAlNWaveCur, VAlNAppert
	  , mLambda, valEpsTargTrue);
	  double temp1 = sqrt(valTheorDisp);

	  double temp2 = sqrt(dispCalc);


	 wchar_t *wcharrFileNames = new wchar_t[NCols * 30];
	 memset(wcharrFileNames, 0,NCols * 30* sizeof (wchar_t));
	 wcscpy(&wcharrFileNames[0],L"n");
	 wcscpy(&wcharrFileNames[30],L"eps_True");
	 wcscpy(&wcharrFileNames[30 * 2],L"eps_Est");
	 wcscpy(&wcharrFileNames[30 * 3],L"eps_Del");
	 wcscpy(&wcharrFileNames[30 * 4],L"SKZ_Theor");
	 wcscpy(&wcharrFileNames[30 * 5],L"MinusSKZ_Theor");
	 wcscpy(&wcharrFileNames[30 * 6],L"SKZ_Cur");



	 double scalex = 1.;
	 double scaley = 1.;
		wchar_t wchFoldName[300] ={0};
	wcscpy(  wchFoldName,  mwchOutFold);
	wcscat(wchFoldName, L"\\");
	 for (int i = 1; i < NCols; i++)
	 {
	  TYrWriteShapeFile::WriteOneReport(wchFoldName  // ���� � �����
								  , arrBuff // ������ � ����������� - ������� nBuffRows x nBuffCols
								  ,NCols // - �-�� ���������� � ������� ��������� ���������� � ������
								  ,NRows //  - �-�� �����
								  ,wcharrFileNames //������� � �������� ���������� - ������� nBuffCols x lenName
								  ,30 // ������������ ����� ����� ����������
								  ,0  // ����� ���������� �� ��� X
								  ,i   // ����� ���������� �� ��� Y
								  , scalex  //  ������� �� ��� X
								  , scaley  // ������� �� ��� Y
								   ) ;
	 }

free(arrBuff);

	///


	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  wchFoldName);
	wcscat(wchAxesFileName, L"AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,0., 10000.
	 ,-10000.,10000., 2.) ;












	  int iii=0;
 /*
	 bool bLat = false;
	 double valSigE = -1.,valSigQ = -1.;
	  double valDistzahv = mFar_2D.calc_TwoTargsZahvatDist(mAntpCoef,mTrajectH, mEPR
	   ,mPowerPrd, mKYPrd, mEtalonSign , mAntH
	   ,&bLat  , &valSigE ,&valSigQ  ) ;

	   int ia = (valDistzahv * 10.) ;
	   LabeledEdit2->Text = ((double)ia)/ 10.;

	   ia = (valSigE * 1000. * 100.) ;
	   LabeledEdit3->Text = ((double)ia)/ 100.;

	   ia = (valSigQ * 1000. * 100.) ;
	   LabeledEdit19->Text = ((double)ia)/ 100.;

	 const int NCols = 21;

	 double  valDMax = 6000.; // ���� ���������
	 double  valDMin = 500.; // ���� ���������
	  double step = 100.;
	 const int NRows = (valDMax - valDMin)/ step;
	 double *arrBuff  = new double [NRows* NCols];

	  TFar Far0(mFar_2D, true);
	// TFar Far(Far0, 4) ;
		int quantFalseSign = 1;
	  // �������� �������� ���� ���� � ���������
	  double valBetTargTrue = 0., valBetAntpTrue =0.;

	 TComp *cmparrTrueAmMeasures = new TComp[ mFar_2D.mNumAMCols * mFar_2D.mNumAMRows];
	 TComp *cmparrNoisedAmMeasures = new TComp [ mFar_2D.mNumAMCols * mFar_2D.mNumAMRows];



	 double arrDisp[4]  = {0.};
	 arrDisp[0]=arrDisp[1]=arrDisp[2]=arrDisp[3]=  mNoiseSKZ_5P10* mNoiseSKZ_5P10/4.;
     TFaceta Faceta(mNumEmitRows,mdEmitRow, mLambda) ;
	 TFar Far1( 4, mdAMRow, mLambda  , Faceta, arrDisp) ;
     double valEstEpsTarg = 0.,valEstEpsAntp = 0., arrMtrxCorrEps [4]= {0.};
	 TComp cmpKTarg , cmpKAntp;
	 for (int i =0; i < NRows; i++)
	 {
		double valHorDist =  valDMax - step * ((double)i);
		double valDist = sqrt( valHorDist * valHorDist + mTrajectH * mTrajectH);
		// ���������� ��������� ���� ����� ���� � ��������
		double valEpsTargTrue = atan((mTrajectH - mAntH)/  valHorDist);
		double valEpsAntpTrue = -atan((mTrajectH + mAntH)/  valHorDist);
		///

		// �������� ������� ����
		double valAmpTarg = mFar_2D.calcCurrentSignalAmpl (valDist, mEPR
		,mPowerPrd, mKYPrd, mEtalonSign) ;
		///

		// �������� ������� ��������
		double valAmpAntp =  valAmpTarg  *mAntpCoef;
		///

		// �������� ������� ���� � �������� ��������

		double valTargSignPhase = 2. * M_PI * getRand01( );
	   //	TSingleSign SignTarg(valBetTargTrue, valEpsTargTrue,  valAmpTarg, valTargSignPhase ) ;

		double valAntpSignPhase = 2. * M_PI * getRand01( );
	   //	TSingleSign SignAntp(valBetAntpTrue, valEpsAntpTrue, valAmpAntp, valAntpSignPhase) ;
		///


         	// ���������� ������� ��������� ���������
		//  ���������� �������� ��������  ��������
		////////////////////////////////////////////////////////////
		TComp arrFourRowMeasures[4];
		memset(arrFourRowMeasures, 0, 4 * sizeof(TComp));







	TComp pcmpS[4];


	TComp cmpTrueKTarg (valAmpTarg * cos(valTargSignPhase), valAmpTarg * sin(valTargSignPhase) );
	TComp cmpTrueKAntp (valAmpAntp * cos(valAntpSignPhase), valAmpAntp * sin(valAntpSignPhase) );

	Far1.ImitateMeasureArray(valEpsTargTrue,cmpTrueKTarg, valEpsAntpTrue,cmpTrueKAntp, pcmpS,arrFourRowMeasures);

		// ����� ��������� ��������
		Far1.fncEstimateMsd(arrFourRowMeasures, &valEstEpsTarg, &valEstEpsAntp
		, &cmpKTarg , &cmpKAntp , arrMtrxCorrEps );
		double valEstHMSD = mAntH + sin(valEstEpsTarg) * valDist;
		///


		// ������ ���
		double valRSMDisp = -1.;
		double valEstEpsRSM = Far1.fncEstimateRsmTetta(arrFourRowMeasures,  &valRSMDisp ) ;
		double valEstHRSM = mAntH + sin(valEstEpsRSM) * valDist;

		// ���������� ������
		arrBuff[i *  NCols] =  valHorDist*0.1;
		for (int j = 0; j < 4; j++)
		{
		  arrBuff[i *  NCols + 1 + j *2] = arrFourRowMeasures[j].modul();  // ����
		  arrBuff[i *  NCols + 1 + j *2 +1] = arrFourRowMeasures[j].phase() * 180./ M_PI; // ����
		}
		arrBuff[i *  NCols + 9] = valAmpTarg/ sqrt(Far1.calcNoiseDisp()) ;   // ������/ ���
		arrBuff[i *  NCols + 10] =  valEstEpsRSM * 1000.;// ��� ����������
		arrBuff[i *  NCols + 11] =  sqrt(valRSMDisp) * 1000.;  // ��� ���
		arrBuff[i *  NCols + 12] = valEstHRSM ; // ������ ������ ���

		arrBuff[i *  NCols + 13] =  valEstEpsTarg * 1000.; // ������ �� ���� ���
		arrBuff[i *  NCols + 14] =  valEstEpsAntp * 1000.; // ������ �� ���� ���
		arrBuff[i *  NCols + 15] =  sqrt(arrMtrxCorrEps[0]) * 1000.; // ��� ������ ���� ���� ���
		arrBuff[i *  NCols + 16] =  sqrt(arrMtrxCorrEps[3]) * 1000.; // ��� ������ ���� ���� ���
		arrBuff[i *  NCols + 17] =   valEstHMSD ; // ������ ������ ���
		arrBuff[i *  NCols + 18] =  valEpsTargTrue * 1000.; // �������� �� ����
		arrBuff[i *  NCols + 19] =  valEpsAntpTrue * 1000.; // �������� �� ����
		arrBuff[i *  NCols + 20] =  mTrajectH ; // �������� ������
	 }

	 wchar_t *wcharrFileNames = new wchar_t[NCols * 30];
	 memset(wcharrFileNames, 0,NCols * 30* sizeof (wchar_t));
	 wcscpy(&wcharrFileNames[0],L"D");
	 wcscpy(&wcharrFileNames[30],L"A0");
	 wcscpy(&wcharrFileNames[30 * 2],L"Fi0");
	 wcscpy(&wcharrFileNames[30 * 3],L"A1");
	 wcscpy(&wcharrFileNames[30 * 4],L"Fi1");
	 wcscpy(&wcharrFileNames[30 * 5],L"A2");
	 wcscpy(&wcharrFileNames[30 * 6],L"Fi2");
	 wcscpy(&wcharrFileNames[30 * 7],L"A3");
	 wcscpy(&wcharrFileNames[30 * 8],L"Fi3");
	 wcscpy(&wcharrFileNames[30 * 9],L"SigNoise");
	 wcscpy(&wcharrFileNames[30 * 10],L"EstEps_RSM");
	 wcscpy(&wcharrFileNames[30 * 11],L"SKZ_RSM");
	 wcscpy(&wcharrFileNames[30 * 12],L"EstH_RSM");
	 wcscpy(&wcharrFileNames[30 * 13],L"EstTargEps_MSD");
	 wcscpy(&wcharrFileNames[30 * 14],L"EstAntpEps_MSD");
	 wcscpy(&wcharrFileNames[30 * 15],L"SKZ_TargEps_MSD");
	 wcscpy(&wcharrFileNames[30 * 16],L"SKZ_AntpEps_MSD");
	 wcscpy(&wcharrFileNames[30 * 17],L"EstH_MSD");
	 wcscpy(&wcharrFileNames[30 * 18],L"EpsTarg_True");
	 wcscpy(&wcharrFileNames[30 * 19],L"EpsAntp_True");
	 wcscpy(&wcharrFileNames[30 * 20],L"H_True");


	 double scalex = 1.;
	 double scaley = 1.;
		wchar_t wchFoldName[300] ={0};
	wcscpy(  wchFoldName,  mwchOutFold);
	wcscat(wchFoldName, L"\\");
	 for (int i = 1; i < NCols; i++)
	 {
	  TYrWriteShapeFile::WriteOneReport(wchFoldName  // ���� � �����
								  , arrBuff // ������ � ����������� - ������� nBuffRows x nBuffCols
								  ,NCols // - �-�� ���������� � ������� ��������� ���������� � ������
								  ,NRows //  - �-�� �����
								  ,wcharrFileNames //������� � �������� ���������� - ������� nBuffCols x lenName
								  ,30 // ������������ ����� ����� ����������
								  ,0  // ����� ���������� �� ��� X
								  ,i   // ����� ���������� �� ��� Y
								  , scalex  //  ������� �� ��� X
								  , scaley  // ������� �� ��� Y
								   ) ;
	 }

free(arrBuff);

	///


	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  wchFoldName);
	wcscat(wchAxesFileName, L"AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,0., 10000.
	 ,0.,10000., 2.) ;
     delete [] cmparrTrueAmMeasures;
	 delete [] cmparrNoisedAmMeasures ;


 */
  
}
//---------------------------------------------------------------------------
/*
void __fastcall TForm1::Button3Click(TObject *Sender)
{

}      */
//---------------------------------------------------------------------------
/*
void __fastcall TForm1::Button4Click(TObject *Sender)
{

} */
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------








//---------------------------------------------------------------------------

void __fastcall TForm1::Button2Click(TObject *Sender)
{
 // ������
double parrArgTab1[] = {0.,4.,5.,6.,7.,8.,9.,10.
};
int NumColsTab1 = 8;

// ���������
double parrArgTab2[] = {800., 6000.
};
int NumColsTab2 = 2;

// ������
double parrArgTab3[] = {4.,  20.,50.,100.,300.,1001.
};
int NumColsTab3 =  6;


double parrVal[] = {1.,1.,1.,    0.166, 0.,0.,0.,0. ,    1.,1.,1.,1.,0.997, 0.981,0.929, 0.

									 ,1.,1.,1.,    0.166, 0.,0.,0.,0.,    1.,1.,1.,1.,0.997, 0.981,0.929, 0.
								 ,1.,1.,1.,    0.252, 0.,0.,0.,0.,    1.,1.,1.,1.,0.997, 0.981,0.929, 0.
									 ,1.,1.,1.,    0.317, 0.,0.,0.,0.,    1.,1.,1.,1.,0.997, 0.981,0.929, 0.
									 ,1.,1.,0.804, 0.462, 0.,0.,0.,0.,    1.,1.,1.,1.,0.997, 0.981,0.929, 0.
									 ,1.,1.,0.5, 0.25, 0.,0.,0.,0.,       1.,1.,1.,1.,0.997, 0.981,0.929, 0.
										};

 TTable_3D Table( parrArgTab1,  NumColsTab1, parrArgTab2,  NumColsTab2
 ,parrArgTab3,  NumColsTab3, parrVal ) ;


 double visota = 500;
 double dalnost = 3400.;
 double promah = 6;



double ver = Table.calcValue(visota, dalnost, promah);
String s_22= L"H = ";
 s_22 = s_22 + visota + L"; D = " + dalnost + L" ; Promah = " +  promah + L" ; Ver = " + ver;

ShowMessage(s_22);

return;
/*double p = 50. / 2/ M_PI/ sqrt(9416 * 1757.);
//p = 0.01207;
double x = 1. - exp(1666. * log(1.-p)/4.);
int iii = 0;
*

/*
// ������������ ��� ���
int nvars = 4;
int  nrows =2;
int bvars = 0;

double a[] ={1.,3.,1.,4.
						,2.,1.,4.,0.} ;
for (int i =0; i < nvars * nrows; i++)
{
	a[i] += 1.;
}


double f[] ={ -1.,-1.,-1.,-1.
							};
double b[] =  { 1.,1.,1.,1.
							};

int  nrows_eq = 0;

double a_eq[6] = {1.,1.,1.,0.,0.,0.
								 };
double b_eq[1] = {1.
								 };
double lb[4] ={0.,0.,0.,0.
							 };
double ub[] = {10000.,10000.,100000.,10000.,10000.,100000.
							 };

int ix [1] ={0};
double x [10] ={0.};
double fval = -1.;

	 int irez =  LinNumericalSolver(  nvars, bvars, f,  nrows,
						a, b, nrows_eq,	a_eq, b_eq, lb,
							ub,ix, x, fval) ;

							for (int i = 0; i < nvars; i++)
							{
							 x[i] = -x[i] / fval;
							}
							int iii0 = 0;
const int NUmRows = 3;
const int NUmCols = 3;
double arr_x [100] = {0.}, arr_y [100] = {0.};
double fval1 = -1.;
double arr_a[] = {1.,2.,3.
									,4.,0.,1.
									,2.,3.,0.} ;
TGameTheory::solvMartrxGame(arr_a, NUmRows,  NUmCols
									, arr_x, arr_y, fval1) ;

	*/





TURPolygon *purPlg =  (TURPolygon *)malloc(sizeof( TURPolygon)* 10);
 TURPolygon **ppurPlg = &purPlg;
 int quantPlg =1;
TURPolygon::ReadSHPFile(L"E:\\�����\\F14\\plg2.shp",ppurPlg,  &quantPlg)  ;
TURPolygon plg_temp0 = (*ppurPlg)[0];
for (int i =0; i < plg_temp0.NumPoints; i++)
{
 plg_temp0.Points[i].Y = plg_temp0.Points[i].Y * 0.6961;
}


TURPolygon::WriteSetSHPFiles(L"E:\\�����\\F14\\Image20.shp",&plg_temp0 ,  quantPlg) ;
plg_temp0.WriteToASCII__(L"E:\\�����\\F14\\plg2.txt");

(*ppurPlg)[0].WriteToASCII__(L"E:\\�����\\F14\\AboveWiev.txt");
 TURPolygon plg4 = (*ppurPlg)[0].SimOtragenieTransform();
// TURPolygon plg4 = (*ppurPlg)[0].SimOtragenieY_and_flip();
 TURPolygon::WriteSetSHPFiles(L"E:\\�����\\F14\\Image10.shp",&plg4 ,  quantPlg) ;
TURPointXY pntSdvig2(0.,0.);
double valRastigenie = 1.;
TURPolygon plg5 = (*ppurPlg)[0].LinTransform(-M_PI/2. , pntSdvig2, valRastigenie ) ;
TURPolygon::WriteSetSHPFiles(L"E:\\Ametist\\toVertorGraf\\MyF14\\nos0.shp",&plg5 ,  quantPlg) ;




valRastigenie  = 16.34/ (7.16 -4.6);
TURPolygon plg6 = plg5.LinTransform(0. , pntSdvig2, valRastigenie ) ;
TURPolygon::WriteSetSHPFiles(L"E:\\Ametist\\toVertorGraf\\MyF14\\nos1.shp",&plg6 ,  quantPlg) ;

TURPointXY pntSdvig(-37.5,  16.3 );

TURPolygon plg7 = plg6.LinTransform(0. , pntSdvig, 1 ) ;
TURPolygon::WriteSetSHPFiles(L"E:\\Ametist\\toVertorGraf\\MyF14\\nos2.shp",&plg7 ,  quantPlg) ;
plg7.WriteToASCII__(L"E:\\Ametist\\toVertorGraf\\TEXTs\\F14nos.txt");
TURPointXY pntSdvig1(0.,0.);
/*
TURPolygon plg1 = plg0.LinTransform(-M_PI/2. , pntSdvig1, valRastigenie ) ;

 valRastigenie  = 18.86/ ((*ppurPlg)[0].Points[19].Y - (*ppurPlg)[0].Points[20].Y);
TURPolygon plg2 = plg1.LinTransform(0. , pntSdvig1, valRastigenie ) ;
TURPolygon::WriteSetSHPFiles(L"E:\\Ametist\\toVertorGraf\\MyF14\\MyKrilo1.shp",&plg2,  quantPlg) ;
TURPolygon plg3 = plg2.SimOtragenieTransform();
TURPolygon::WriteSetSHPFiles(L"E:\\Ametist\\toVertorGraf\\MyF14\\MyKrilo2.shp",&plg3 ,  quantPlg) ;
 */
free(purPlg);
}
//---------------------------------------------------------------------------



/*
void __fastcall TForm1::ExchangeData()
{
 //	fncInputData();

	if(RadioGroup1->ItemIndex ==0 )
	{
	// ��� ������� ������ ���� (����� �����)
		Form3->mSigSins = mSigSins;

		// ��� ������� ������ ������ �������� ����� (����� �����)
		Form3->mSig_d_po_dt_Sins = mSig_d_po_dt_Sins;

		Form3->mEtalonSign = mEtalonSign;

		// ������� ���
		Form3->mVesselWidth =mVesselWidth ; // ������(�)
		Form3->mVesselLength = mVesselLength ;
		memcpy(Form3-> marrFarParallacs, marrFarParallacs, 3 * sizeof(double));

		Form3->mMaxQ = mMaxQ ; /// ������������ ���� ���������� �� ������������ �����(��������� ���� ��������)
		Form3->mT_Q = mT_Q; // ������ ��������
		Form3->mMaxPsi = mMaxPsi ;// ������������ ���� ������� �����(���������)
		Form3->mT_Psi = mT_Psi ; // ������ ������� �����
		Form3->mMaxTet = mMaxTet ; //������������ ���� ��������� �����(���������)
		Form3->mT_Tet = mT_Tet; // ������ �������� �����
		Form3->mMaxVert = mMaxVert  ;

		// ���������� ��������  ������� ������
		Form3->mQ0 = mQ0  ; // ����������� ����
		Form3->mVVess = mVVess  ;// �������� ������� ������ 20 �����
		memcpy(Form3-> marrDelt,  marrDelt, 4 * sizeof(double)); //  ��������� ����

		//������������ ��������� ��������� ������ ������� � ��� �� 100 �
		Form3->mMaxAmp_AftFlexure = mMaxAmp_AftFlexure;
		// ������ ��������� ��������� ������
		Form3->mT_AftFlexure = mT_AftFlexure;
		//������������ ��������� ��������� ������ ������� � ��� �� 100 �
		Form3->mMaxAmp_BoardFlexure = mMaxAmp_BoardFlexure;
		// ������ ��������� ��������� ������
		Form3->mT_BoardFlexure = mT_BoardFlexure;

		// 3.1 �������� ����
		TSins mSins ;
		Form3->mMaxSig_Q = mMaxSig_Q ;
		Form3->mMaxSig_Psi = mMaxSig_Psi  ;
		Form3->mMaxSig_Tet = mMaxSig_Tet  ;
		Form3->mMaxSig_dQdt = mMaxSig_dQdt ;
		Form3->mMaxSig_dPsidt = mMaxSig_dPsidt ;
		Form3->mMaxSig_dTetdt = mMaxSig_dTetdt ;
		Form3->mK1 = mK1         ;
		Form3->mSigV = mSigV      ;
		Form3->mSigH = mSigH     ;
		Form3->mMaxSig_H = mMaxSig_H ;
		Form3->mMaxSig_VH = mMaxSig_VH ;

		// ������
		Form3->mDriverSigBet = mDriverSigBet ;// �������� ��������� ���� Bet �������
		Form3->mDriverSigEps = mDriverSigEps ;// �������� ��������� ���� Eps  ������� (���� �����)
		Form3->mDriverDynamicSigBet = mDriverDynamicSigBet ;// �������� ��������� ���� �����  �������
		Form3->mDriverDynamicSigEps = mDriverDynamicSigEps ;// ��������  ������� ��������� ���� �����
		//

		//���� ����������


		Form3->mControlSyst  = mControlSyst;
		memcpy(Form3-> marrArtParral,  marrArtParral, 3 * sizeof(double));  // ������ ���������� ��

		// 	��������� ��������� �������� ��������  �����
		Form3->mSigDrivAY_dU_po_dt = mSigDrivAY_dU_po_dt ;

		Form3->mTransmitAnt = mTransmitAnt;
		Form3->mFar_2D = mFar_2D;
		Form3->mEnvironment = mEnvironment;

	}

	if(RadioGroup1->ItemIndex ==1 )
	{
    	// ��� ������� ������ ���� (����� �����)
		Form4->mSigSins = mSigSins;

		// ��� ������� ������ ������ �������� ����� (����� �����)
		Form4->mSig_d_po_dt_Sins = mSig_d_po_dt_Sins;

	 Form4->mEtalonSign = mEtalonSign;


		// ������� ���
		Form4->mVesselWidth =mVesselWidth ; // ������(�)
		Form4->mVesselLength = mVesselLength ;
		memcpy(Form4-> marrFarParallacs, marrFarParallacs, 3 * sizeof(double));

		Form4->mMaxQ = mMaxQ ; /// ������������ ���� ���������� �� ������������ �����(��������� ���� ��������)
		Form4->mT_Q = mT_Q; // ������ ��������
		Form4->mMaxPsi = mMaxPsi ;// ������������ ���� ������� �����(���������)
		Form4->mT_Psi = mT_Psi ; // ������ ������� �����
		Form4->mMaxTet = mMaxTet ; //������������ ���� ��������� �����(���������)
		Form4->mT_Tet = mT_Tet; // ������ �������� �����
		Form4->mMaxVert = mMaxVert  ;

		// ���������� ��������  ������� ������
		Form4->mQ0 = mQ0  ; // ����������� ����
		Form4->mVVess = mVVess  ;// �������� ������� ������ 20 �����
		memcpy(Form4-> marrDelt,  marrDelt, 4 * sizeof(double)); //  ��������� ����

		//������������ ��������� ��������� ������ ������� � ��� �� 100 �
		Form4->mMaxAmp_AftFlexure = mMaxAmp_AftFlexure;
		// ������ ��������� ��������� ������
		Form4->mT_AftFlexure = mT_AftFlexure;
		//������������ ��������� ��������� ������ ������� � ��� �� 100 �
		Form4->mMaxAmp_BoardFlexure = mMaxAmp_BoardFlexure;
		// ������ ��������� ��������� ������
		Form4->mT_BoardFlexure = mT_BoardFlexure;

		// 3.1 �������� ����
		TSins mSins ;
		Form4->mMaxSig_Q = mMaxSig_Q ;
		Form4->mMaxSig_Psi = mMaxSig_Psi  ;
		Form4->mMaxSig_Tet = mMaxSig_Tet  ;
		Form4->mMaxSig_dQdt = mMaxSig_dQdt ;
		Form4->mMaxSig_dPsidt = mMaxSig_dPsidt ;
		Form4->mMaxSig_dTetdt = mMaxSig_dTetdt ;
		Form4->mK1 = mK1         ;
		Form4->mSigV = mSigV      ;
		Form4->mSigH = mSigH     ;
		Form4->mMaxSig_H = mMaxSig_H ;
		Form4->mMaxSig_VH = mMaxSig_VH ;

		// ������
		Form4->mDriverSigBet = mDriverSigBet ;// �������� ��������� ���� Bet �������
		Form4->mDriverSigEps = mDriverSigEps ;// �������� ��������� ���� Eps  ������� (���� �����)
		Form4->mDriverDynamicSigBet = mDriverDynamicSigBet ;// �������� ��������� ���� �����  �������
		Form4->mDriverDynamicSigEps = mDriverDynamicSigEps ;// ��������  ������� ��������� ���� �����
		//



		Form4->mControlSyst  = mControlSyst;
		memcpy(Form4-> marrArtParral,  marrArtParral, 3 * sizeof(double));  // ������ ���������� ��

		// 	��������� ��������� �������� ��������  �����
		Form4->mSigDrivAY_dU_po_dt = mSigDrivAY_dU_po_dt ;

		Form4->mTransmitAnt = mTransmitAnt;
		Form4->mFar_2D = mFar_2D;
		Form4->mEnvironment = mEnvironment;

	}

	if(RadioGroup1->ItemIndex ==2 )
	{

		Form2->mEtalonSign = mEtalonSign;
		// ��� ������� ������ ���� (����� �����)
		Form2->mSigSins = mSigSins;

		// ��� ������� ������ ������ �������� ����� (����� �����)
		Form2->mSig_d_po_dt_Sins = mSig_d_po_dt_Sins;

		// ������� ���
		Form2->mVesselWidth =mVesselWidth ; // ������(�)
		Form2->mVesselLength = mVesselLength ;
		memcpy(Form2-> marrFarParallacs, marrFarParallacs, 3 * sizeof(double));

		Form2->mMaxQ = mMaxQ ; /// ������������ ���� ���������� �� ������������ �����(��������� ���� ��������)
		Form2->mT_Q = mT_Q; // ������ ��������
		Form2->mMaxPsi = mMaxPsi ;// ������������ ���� ������� �����(���������)
		Form2->mT_Psi = mT_Psi ; // ������ ������� �����
		Form2->mMaxTet = mMaxTet ; //������������ ���� ��������� �����(���������)
		Form2->mT_Tet = mT_Tet; // ������ �������� �����
		Form2->mMaxVert = mMaxVert  ;

		// ���������� ��������  ������� ������
		Form2->mQ0 = mQ0  ; // ����������� ����
		Form2->mVVess = mVVess  ;// �������� ������� ������ 20 �����
		memcpy(Form2-> marrDelt,  marrDelt, 4 * sizeof(double)); //  ��������� ����

		//������������ ��������� ��������� ������ ������� � ��� �� 100 �
		Form2->mMaxAmp_AftFlexure = mMaxAmp_AftFlexure;
		// ������ ��������� ��������� ������
		Form2->mT_AftFlexure = mT_AftFlexure;
		//������������ ��������� ��������� ������ ������� � ��� �� 100 �
		Form2->mMaxAmp_BoardFlexure = mMaxAmp_BoardFlexure;
		// ������ ��������� ��������� ������
		Form2->mT_BoardFlexure = mT_BoardFlexure;

		// 3.1 �������� ����
		TSins mSins ;
		Form2->mMaxSig_Q = mMaxSig_Q ;
		Form2->mMaxSig_Psi = mMaxSig_Psi  ;
		Form2->mMaxSig_Tet = mMaxSig_Tet  ;
		Form2->mMaxSig_dQdt = mMaxSig_dQdt ;
		Form2->mMaxSig_dPsidt = mMaxSig_dPsidt ;
		Form2->mMaxSig_dTetdt = mMaxSig_dTetdt ;
		Form2->mK1 = mK1         ;
		Form2->mSigV = mSigV      ;
		Form2->mSigH = mSigH     ;
		Form2->mMaxSig_H = mMaxSig_H ;
		Form2->mMaxSig_VH = mMaxSig_VH ;

		// ������
		Form2->mDriverSigBet = mDriverSigBet ;// �������� ��������� ���� Bet �������
		Form2->mDriverSigEps = mDriverSigEps ;// �������� ��������� ���� Eps  ������� (���� �����)
		Form2->mDriverDynamicSigBet = mDriverDynamicSigBet ;// �������� ��������� ���� �����  �������
		Form2->mDriverDynamicSigEps = mDriverDynamicSigEps ;// ��������  ������� ��������� ���� �����
		//



		Form2->mControlSyst  = mControlSyst;
		memcpy(Form2-> marrArtParral,  marrArtParral, 3 * sizeof(double));  // ������ ���������� ��

		// 	��������� ��������� �������� ��������  �����
		Form2->mSigDrivAY_dU_po_dt = mSigDrivAY_dU_po_dt ;

		Form2->mTransmitAnt = mTransmitAnt;
		Form2->mFar_2D = mFar_2D;
		Form2->mEnvironment = mEnvironment;


	}
}
*/
//---------------------------------------------------------------------------




// ---------------------------------------------------------------------------

void __fastcall TForm1::create5P10()
{
   	//���������
	 mEtalonAmp = 500.;
	//���������
	 mEtalonDist = 12000.;
	//���
	 mEtalonAPR = 1.;
	//��� ����� ���� ��������� ��������� 5�10
	 mNoiseSKZ_5P10 = 18.735;
	// ��� �������� ������� �������� ��������� ��������� 5�10
	 mEtalonSigAmplFact_5P10 = 0.01;
	//
	// �������� �� ��������
	 mEtalonPowerPrd = 4000.;
	// �� �� ��������
	 mEtalonKYPrd = 840.;
	// �������� �� �����
	 mEtalonKYPriem = 5200.;

	 mEtalonSign = TEtalonSign(mEtalonAmp, mEtalonDist, mEtalonAPR,
			mNoiseSKZ_5P10, mEtalonSigAmplFact_5P10, mEtalonPowerPrd, mEtalonKYPrd
		  , mEtalonKYPriem);

	///

	// �������
	// ������ �� ��
	// �-�� ����������� �� �����������
	 mNumEmitCols = 8;
	// �-�� ����������� �� ���������
	 mNumEmitRows = 8;
		// ����� �����
		mLambda = VAL_C / 8.* 100. / 1000000000.;

			// ���������� ����� ������������
		mdEmitCol = mLambda * 0.55;
		mdEmitRow = mLambda * 0.55;


	// ������ �� ���
	// �-�� ��  �� �����������
	 mNumAMCols = 8;
	// �-�� �� �� ���������
	 mNumAMRows = 8;
		// ���������� ����� �� �� ���������
	mdAMRow = ((double)mNumEmitRows) * mdEmitRow;
	// ���������� ����� �� �� �����������
	mdAMCol = ((double)mNumEmitCols) * mdEmitCol;
	// ��� ���� � ��������� ���������

	 for (int i = 0; i < MAX_QUANT_AM; i++)
	{
		mbarrAM[i] = true;
	}

	 // for (int i = 0; i < 3; i++)  // ��� ������� !!!
	 // {
	 //	marrArtParral[i] = 0.; // ������ ���������� ��
	 //	marrFarParallacs[i]  = 0.; //

	 // }


	// �������� �� ��������
	 mPowerPrd =4000.;
	// �� �� ��������
	 mKYPrd = 840.;
	 mTransmitAnt.mKYPrd  =mKYPrd;
	 mTransmitAnt.mPowerPrd = mPowerPrd;


	///

const double SigEmitNoise = mNoiseSKZ_5P10 / sqrt(((double)(8 * 16 * 28)));
	const double SigEmitAmplFact = mEtalonSigAmplFact_5P10 / sqrt
		(((double)(8 * 16 * 28)));
	mAM_2D = TAM_2D(mNumEmitCols, mNumEmitRows, mdEmitCol, mdEmitRow,
		SigEmitNoise, SigEmitAmplFact);

	mFar_2D = TFar_2D(mNumAMCols, mNumAMRows, mLambda, mdAMCol, mdAMRow,
		mAM_2D, mbarrAM);

	for (int i = 0; i < mNumAMCols * mNumAMRows; i++) {

		if (!mbarrAM[i]) {
			mFar_2D.mpAm2D[i].mOtklCoefUs = 1.;
			mFar_2D.mpAm2D[i].mSigEmitNoise = 0.;
		}
	}





}

//-------------------------------------------------------
void __fastcall TForm1::fillVessel()
{
	   LabeledEdit15->Text = VAL_C /  mLambda / 10000000.;
	   LabeledEdit13->Text = mLambda;


		LabeledEdit4->Text = mEtalonAmp ;
		// ���������
		LabeledEdit10->Text = mEtalonDist / 1000.;
		// ���
		LabeledEdit16->Text = mEtalonAPR ;

		LabeledEdit17->Text = mNoiseSKZ_5P10 ;
		// ��� �������� ������� �������� ��������� ��������� 5�10
		LabeledEdit18->Text = mEtalonSigAmplFact_5P10 ;

		// �������� �� ��������
		LabeledEdit36->Text = mEtalonPowerPrd ;
		// �� �� ��������
		LabeledEdit38->Text = mEtalonKYPrd ;
		// �������� �� �����
		LabeledEdit37->Text = mEtalonKYPriem ;
		///

		LabeledEdit40->Text = mPowerPrd ;
		// �� �� ��������
		 LabeledEdit39->Text = mKYPrd;
		 ///


		///

		// �-�� ��������  �����������

	LabeledEdit9->Text = mNumEmitCols ;
	// �-�� ����� �����������
	LabeledEdit11->Text = mNumEmitRows ;
	// �-�� �������� ��  (�� ����������� )
	LabeledEdit5->Text = mNumAMCols ;
	// �-�� ����� �� (�� ��������� )
	LabeledEdit6->Text = mNumAMRows ;

		// ���������� ����� ������������ �� �����������
	LabeledEdit12->Text = mdEmitCol;


	LabeledEdit8->Text = mdAMRow;
	LabeledEdit7->Text = mdAMCol;
	///



}
void __fastcall TForm1::Button7Click(TObject *Sender)
{
	create5P10();
	fillVessel();
}
//---------------------------------------------------------------------------





void __fastcall TForm1::LabeledEdit15Change(TObject *Sender)
{
fncInputData();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::LabeledEdit15Exit(TObject *Sender)
{
	fncInputData();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::LabeledEdit31Change(TObject *Sender)
{
    fncInputData();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::LabeledEdit21Change(TObject *Sender)
{
  fncInputData();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::LabeledEdit48Change(TObject *Sender)
{
  fncInputData();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::LabeledEdit54Change(TObject *Sender)
{
  fncInputData();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::LabeledEdit35Change(TObject *Sender)
{
  fncInputData();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::LabeledEdit40Change(TObject *Sender)
{
  fncInputData();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::LabeledEdit39Change(TObject *Sender)
{
  fncInputData();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::LabeledEdit29Change(TObject *Sender)
{
fncInputData();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::LabeledEdit28Change(TObject *Sender)
{
fncInputData();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::LabeledEdit1Change(TObject *Sender)
{
fncInputData();
}
//---------------------------------------------------------------------------

double  __fastcall TForm1::calcDispNew(double eps, double dispSum, double A, double lamb,double d
, double valFDiagr)
{

double temp =  lamb/ M_PI/ d;
double fi = 2. * M_PI * d / lamb * eps;
double D = fncDiagrSinx_div_x(fi);
double sigSq = dispSum * cos(fi/2.)* cos(fi/2.)/A/A/valFDiagr/valFDiagr * temp* temp/2.;

return sigSq;

/*
double temp =  lamb/ M_PI/ d;


double sigSq = dispSum /A/A/valFDiagr/valFDiagr * temp* temp/2.;

return sigSq; */
}
