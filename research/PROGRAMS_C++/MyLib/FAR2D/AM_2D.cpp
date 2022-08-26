//---------------------------------------------------------------------------


#pragma hdrstop
#include <vcl.h>
#include "AM_2D.h"
#include "Comp.h"
#include "SingleSign.h"

//---------------------------------------------------------------------------


#pragma hdrstop
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <stdio.h>
#include <stdlib.h>
#include "Gauss.h"
#include "YrWriteShapeFile.h"
#include "URPolyLine.h"
#include "UrPointXY.h"
//#include "DiagrSinX.h"

//extern const double TET0707;

__fastcall  TAM_2D::TAM_2D()
{
 // �-�� ����������� �� �����������
	 mNumEmitCols = 8;
	// �-�� ����������� �� �����������
	 mNumEmitRows =16;
	// ���������� ����� ������������ �� �����������
	 mdCol = 2.;
	// ���������� ����� ������������ �� ���������
	 mdRow = 2.;
	
	// ��� ����  ����������
	 mSigEmitNoise = 10./ 32.;
	 mOtklCoefUs = 1.;
	 mOtklBetta = 0.;
	 mOtklEps =0.;
	 mSigEmitAmplFact = 0.;


}
// ����������� �����������
__fastcall  TAM_2D::TAM_2D (const TAM_2D &R2)
 {
	 mNumEmitCols = R2.mNumEmitCols ;
	 mNumEmitRows = R2.mNumEmitRows;

	 mdCol = R2.mdCol;
	 mdRow = R2.mdRow;
	 mSigEmitNoise = R2.mSigEmitNoise  ;
	 mOtklCoefUs = R2.mOtklCoefUs;
	 mOtklBetta = R2.mOtklBetta;
	 mOtklEps = R2.mOtklEps;
	 mSigEmitAmplFact = R2.mSigEmitAmplFact;
 }
 // ����� ������ 1
 __fastcall TAM_2D::TAM_2D(const int NumEmitCols,const int NumEmitRows
	  ,const double dCol,const double dRow, const double SigEmitNoise)
 {
	 mNumEmitCols = NumEmitCols ;
	 mNumEmitRows = NumEmitRows;

	 mdCol = dCol;
	 mdRow = dRow;
	 mSigEmitNoise = SigEmitNoise;
	 mOtklCoefUs = 1.;
	 mOtklBetta = 0.;
	 mOtklEps = 0.;


 }
 // ����� ������ 2
 __fastcall TAM_2D::TAM_2D(const int NumEmitCols,const int NumEmitRows
	  ,const double dCol,const double dRow, const double SigEmitNoise
	  , const double OtklCoefUs, const double OklBetta
	  , const double OtklEps )
 {
	 mNumEmitCols = NumEmitCols ;
	 mNumEmitRows = NumEmitRows;

	 mdCol = dCol;
	 mdRow = dRow;
	 mSigEmitNoise = SigEmitNoise;
	 mOtklCoefUs = OtklCoefUs;
	 mOtklBetta = mOtklBetta;
	 mOtklEps = OtklEps ;


 }


// �������� ������������
  TAM_2D &TAM_2D::operator=(const TAM_2D  &R2)
 {
	 mNumEmitCols = R2.mNumEmitCols ;
	 mNumEmitRows = R2.mNumEmitRows;

	 mdCol = R2.mdCol;
	 mdRow = R2.mdRow;
	 mSigEmitNoise = R2.mSigEmitNoise ;
	 mOtklCoefUs = R2.mOtklCoefUs;
	 mOtklBetta = R2.mOtklBetta;
	 mOtklEps = R2.mOtklEps;
	 mSigEmitAmplFact = R2.mSigEmitAmplFact;

	 return *this ;
 }

 // ����� ������ 3
 __fastcall TAM_2D::TAM_2D(const int NumEmitCols,const int NumEmitRows
	  ,const double dCol,const double dRow, const double SigEmitNoise
	  ,const double SigEmitAmplFact)
 {
	 mNumEmitCols = NumEmitCols ;
	 mNumEmitRows = NumEmitRows;

	 mdCol = dCol;
	 mdRow = dRow;
	 mSigEmitNoise = SigEmitNoise;
	 mSigEmitAmplFact = SigEmitAmplFact;
	 mOtklCoefUs = 1.;
	 mOtklBetta = 0.;
	 mOtklEps = 0.;


 }
 // ���������� �������
  double TAM_2D::calcSquare()
  {
	  return mNumEmitCols* mNumEmitRows * mdCol * mdRow ;
  }

 TComp TAM_2D::fncIdealDiagr (const double valBetta,const double valEps,const double valLambda)
 {
  TComp cmpRez(0.,0.), cmpRez0(0.,0.);
  for (int m = 0; m < mNumEmitRows; m++)
  {
	for (int l = 0; l < mNumEmitCols; l++)
	{
		double valq = -( mdCol * cos( valEps) * sin(valBetta ) * ( (((double)mNumEmitCols)-1.)/ 2. - (double)l )
			 + mdRow * sin( valEps) * ( (((double)mNumEmitRows)-1.)/ 2. - (double)m ));
		cmpRez += exp_(valq * 2. * M_PI / valLambda);
	}
  }
  TComp cmpTemp ( cos( valEps) * cos(valBetta ) / ((double)mNumEmitCols)/ ((double)mNumEmitRows), 0.);
  cmpRez0 = cmpRez *  cmpTemp  ;
	 return cmpRez0    ;
 }

 // ���������� ��������� ������� ����������� � ����������� valBetta, valEps � �������� ��������   valAmpl
 // � ���� valPhase (������������ �������� ������  ��)
 //OUTPUT:
 // pcmparrEmitMeasures - ������ �������� �������� ������������   ����������� mNumEmitRows*  mNumEmitCols
 // ���������� ��������� ������
 // ���� pcmpEmitMeasures!= NULL, �� � ����� ������� ������������ ������� � �����������
 TComp TAM_2D::fncImitateSingleTargMeasure (const double ValBetta,const double ValEps ,const double valLambda
	 , const double valAmpl, const double valPhase, TComp * pcmpEmitMeasures)
 {
	 const double valBetta = ValBetta - mOtklBetta;
	 const double valEps = ValEps - mOtklEps;
	 TComp cmpNoise(0.,0.), cmpRez0(0.,0.), cmpSdvig;
	 TComp cmpSignal( valAmpl * cos(valPhase), valAmpl * sin(valPhase));// �������� ������ � ����������� �����
	 TComp cmpNorm ( cos( valEps) * cos(valBetta ) / ((double)mNumEmitCols)/ ((double)mNumEmitRows), 0.);
  TComp *pcmparrEmitMeasures = (TComp *)malloc( mNumEmitRows * mNumEmitCols * sizeof(TComp));
  memset(pcmparrEmitMeasures, 0, mNumEmitRows * mNumEmitCols * sizeof(TComp));
  // ���������� ������� �������� � ����������� �����
  TComp cmpOklCoefUs(mOtklCoefUs, 0.);
  for (int m = 0; m < mNumEmitRows; m++)
  {
	for (int l = 0; l < mNumEmitCols; l++)
	{
	 // double Xml = (((double)mNumEmitCols)-1.)/ 2. *mdCol - ((double)l) * mdCol;
	//  double Zml = (((double)mNumEmitRows)-1.)/ 2. *mdRow - ((double)m) * mdRow;
	 // double valq =  -( Xml *cos( valEps)* sin(valBetta ) +  Zml * sin( valEps));


		double valq = -( mdCol * cos( valEps) * sin(valBetta ) * ( (((double)mNumEmitCols)-1.)/ 2. - (double)l )
	 		 + mdRow * sin( valEps) * ( (((double)mNumEmitRows)-1.)/ 2. - (double)m ));

		cmpSdvig = exp_(valq * 2. * M_PI /valLambda);
		cmpNoise.m_Re = getGauss(0., mSigEmitNoise/ sqrt(2.) ) ;
		cmpNoise.m_Im = getGauss(0., mSigEmitNoise/ sqrt(2.) ) ;
		pcmparrEmitMeasures [ m *mNumEmitCols + l] =  (((cmpSignal* cmpSdvig)* cmpNorm)*cmpOklCoefUs) + cmpNoise;
		cmpRez0 +=  pcmparrEmitMeasures [ m *mNumEmitCols + l];
	}
  }
   if  (pcmpEmitMeasures)
   {
	   memcpy(pcmpEmitMeasures, pcmparrEmitMeasures,mNumEmitRows * mNumEmitCols * sizeof(TComp));
   }
   free(pcmparrEmitMeasures);
	 return cmpRez0;
 }


 // ���������� �������� ��������� � ������������ ���������
 void TAM_2D::createVertIdealDiagrGraphs(const double valLambda, wchar_t *wchFoldName1 )
{
	wchar_t wchFoldName[300] ={0};
	wcscpy(  wchFoldName,  wchFoldName1);
	wcscat(wchFoldName, L"\\");
	double step = M_PI / 3000./ 10.;

	const int nBuffRows = 1500 *2;
	const int nBuffCols =2;
	double  *parrBuff = new double [nBuffRows  * nBuffCols] ;
	memset(parrBuff, 0, nBuffRows * nBuffCols * sizeof(double));

	const int lenName = 30;
	wchar_t wcharrFileNames [150] ={0};
	wcscpy( wcharrFileNames, L"Tetta");
	wcscpy( &wcharrFileNames[30], L"F_Vert");



	const double valBetta = 0.;
	for (int i=0 ; i < nBuffRows; i++)
	{
		parrBuff[ i * nBuffCols] = ((double) (-nBuffRows/2 +i)) / 10.;
		double tet = step * ((double) (-nBuffRows/2 +i));
		TComp cmpTemp = fncIdealDiagr (valBetta,tet, valLambda);

		parrBuff[ i * nBuffCols + 1]= cmpTemp.m_Re;


	}

	double scaley = 100.;
	for (int i=1; i < nBuffCols; i++)
	{
		TYrWriteShapeFile::WriteOneReport(                 wchFoldName  // ���� � �����
		, parrBuff // ������ � ����������� - ������� nBuffRows x nBuffCols
		, nBuffCols // - �-�� ���������� � ������� ��������� ���������� � ������
		, nBuffRows //  - �-�� �����
		,wcharrFileNames //������� � �������� ���������� - ������� nBuffCols x lenName
		,lenName // ������������ ����� ����� ����������
		,0  // ����� ���������� �� ��� X
		,i  // ����� ���������� �� ��� Y
		,1. //  ������� �� ��� X
		,scaley // ������� �� ��� Y
		) ;
	}

	delete []parrBuff;
	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  wchFoldName);
	wcscat(wchAxesFileName, L"AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-parrBuff[ (nBuffRows-1) * nBuffCols]-100000.*step, parrBuff[ (nBuffRows-1) * nBuffCols]+100000.*step
	,0,scaley * 1.1, 50000.*step) ;
	TURPointXY pPoints[2];
	pPoints[0] =  TURPointXY (-100., scaley * sqrt(2.)/2.);
	pPoints[1] =  TURPointXY (100., scaley * sqrt(2.)/2.);
	TURPolyLine pln( pPoints,2) ;

	wcscpy(  wchAxesFileName,  wchFoldName);
	wcscat(wchAxesFileName, L"Line0.shp");
	TURPolyLine::WriteSetSHPFiles(wchAxesFileName,&pln, 1) ;
}

 // ���������� ������� ��������� ��������  �������
 // ������������ �������� ������ �������� parrsignFalse ����������� quantFalseSign
 //� �������� ��������, ���������� �� ���� signTarg
 //  INPUT:
 // signTarg -  ������ ����
 // parrsignFalse - ������ �������� ������ �����(���������)
 // quantFalseSign - ����� ����� �������
 // valLambda -  ����� �����
 //OUTPUT:
 // pcmpTrueEmitMeasures - ������ �������� ��������, ����������� �� ���� ����������� ����������� mNumEmitRows*  mNumEmitCols
 // pcmpNoisedEmitMeasures - ������ ����������� ��������, �� ������ �����������   ����������� mNumEmitRows*  mNumEmitCols
 // ���������� ��������� ����������� ������
 // ���� pcmpTrueEmitMeasures == NULL, �� ���� ������ �� �����������
 // ���� pcmpNoisedEmitMeasures == NULL, �� ���� ������ �� �����������
 /*
 TComp TAM_2D::fncImitateMultiTargMeasure ( TSingleSign signTarg, TSingleSign *parrsignFalse, int quantFalseSign
	 ,const double valLambda, TComp * pcmpTrueEmitMeasures, TComp * pcmpNoisedEmitMeasures)
 {
  TComp cmpRez0(0.,0.);
  TComp *pcmparrTrueEmitMeasures = (TComp *)malloc( mNumEmitRows * mNumEmitCols * sizeof(TComp));
  memset(pcmparrTrueEmitMeasures, 0, mNumEmitRows * mNumEmitCols * sizeof(TComp));
  TComp *pcmparrNoisedEmitMeasures = (TComp *)malloc( mNumEmitRows * mNumEmitCols * sizeof(TComp));
  memset(pcmparrNoisedEmitMeasures, 0, mNumEmitRows * mNumEmitCols * sizeof(TComp));

	// ���������� ������� �������� � ����������� �����
   TComp cmpOklCoefUs(mOtklCoefUs, 0.);
   // ������������ ������� �������� ������ �����, ����������� �� ���� �����������
   for (int n = 0; n < quantFalseSign; n++)
   {
	 double valBetta = parrsignFalse[n].mBet - mOtklBetta;  // �������� ���� ����� �������
	 double valEps   = parrsignFalse[n].mEps   - mOtklEps; // �������� ���� ����� �������
	 TComp cmpSignal(parrsignFalse[n].mAmpl * cos(parrsignFalse[n].mPhase)
			  , parrsignFalse[n].mAmpl * sin(parrsignFalse[n].mPhase));// �������� ������ � ����������� �����
	 TComp cmpNorm ( cos( valEps) * cos(valBetta ) / ((double)mNumEmitCols)/ ((double)mNumEmitRows), 0.);


		for (int m = 0; m < mNumEmitRows; m++)
		{
		for (int l = 0; l < mNumEmitCols; l++)
		{
			double valq = -( mdCol * cos( valEps) * sin(valBetta ) * ( (((double)mNumEmitCols)-1.)/ 2. - (double)l )
				 + mdRow * sin( valEps) * ( (((double)mNumEmitRows)-1.)/ 2. - (double)m ));
			TComp cmpSdvig = exp_(valq * 2. * M_PI /valLambda);
			pcmparrTrueEmitMeasures[ m *mNumEmitCols + l] =  (((cmpSignal* cmpSdvig)* cmpNorm)*cmpOklCoefUs) ;//+ cmpNoise;
		}
		}
	}
	///

	// ����������� ��������� �������
	 double valBetta = signTarg.mBet - mOtklBetta;  // �������� ���� ����� �������
	 double valEps   = signTarg.mEps   - mOtklEps; // �������� ���� ����� �������
	 TComp cmpSignal(signTarg.mAmpl * cos(signTarg.mPhase)
			  , signTarg.mAmpl * sin(signTarg.mPhase));// �������� ������ � ����������� �����
	 TComp cmpNorm ( cos( valEps) * cos(valBetta ) / ((double)mNumEmitCols)/ ((double)mNumEmitRows), 0.);


		for (int m = 0; m < mNumEmitRows; m++)
		{
		for (int l = 0; l < mNumEmitCols; l++)
		{
			double valq = -( mdCol * cos( valEps) * sin(valBetta ) * ( (((double)mNumEmitCols)-1.)/ 2. - (double)l )
				 + mdRow * sin( valEps) * ( (((double)mNumEmitRows)-1.)/ 2. - (double)m ));
			TComp cmpSdvig = exp_(valq * 2. * M_PI /valLambda);
			pcmparrTrueEmitMeasures[ m *mNumEmitCols + l] =  (((cmpSignal* cmpSdvig)* cmpNorm)*cmpOklCoefUs) ;//+ cmpNoise;
		}
		}

	///

	// ������������ ������� ����������� ��������
		for (int m = 0; m < mNumEmitRows; m++)
		{
		for (int l = 0; l < mNumEmitCols; l++)
		{
			TComp cmpNoise (0.,0.);
			cmpNoise.m_Re = getGauss(0., mSigEmitNoise/ sqrt(2.) ) ;
			cmpNoise.m_Im = getGauss(0., mSigEmitNoise/ sqrt(2.) ) ;
			pcmpNoisedEmitMeasures[ m *mNumEmitCols + l] =   pcmparrTrueEmitMeasures[ m *mNumEmitCols + l] + cmpNoise;
			cmpRez0 +=  pcmpNoisedEmitMeasures[ m *mNumEmitCols + l];
		}
		}


   if  (pcmpTrueEmitMeasures)
   {
	   memcpy(pcmpTrueEmitMeasures, pcmparrTrueEmitMeasures,mNumEmitRows * mNumEmitCols * sizeof(TComp));
   }
   if  (pcmpNoisedEmitMeasures)
   {
	   memcpy(pcmpNoisedEmitMeasures, pcmparrNoisedEmitMeasures,mNumEmitRows * mNumEmitCols * sizeof(TComp));
   }
   free(pcmparrNoisedEmitMeasures);
   free(pcmparrTrueEmitMeasures);
   return cmpRez0;
 }
  */
 TAM_2D TAM_2D::create_AM_5P10(const double valSigNoise)
 {
	return TAM_2D (8, 16, 2.,  2.,valSigNoise/ sqrt(28.*8.*16.)) ;
 }

 // ���������� ����������� ���� ���������
double TAM_2D::calcSumDisp()
{

  return ((double)mNumEmitCols * mNumEmitRows) * mSigEmitNoise * mSigEmitNoise;
}


 // ���������� ���������  �������� �� ��
double TAM_2D::calcSumAmplFactDisp()
{
  return   mSigEmitAmplFact * mSigEmitAmplFact /((double)mNumEmitCols * mNumEmitRows);
}


//    �������������!!!!!
 // ���������� ������� ���������   �������    ��������� ����

 //  INPUT:
 // signTarg -  ������ ����
 // valLambda -  ����� �����
 //OUTPUT:
 // pcmpTrueEmitMeasures - ������ �������� ��������, ����������� �� ���� ����������� ����������� mNumEmitRows*  mNumEmitCols
 // pcmpNoisedEmitMeasures - ������ ����������� ��������, �� ������ �����������   ����������� mNumEmitRows*  mNumEmitCols
 // ���������� ��������� ����������� ������
 // ���� pcmpTrueEmitMeasures == NULL, �� ���� ������ �� �����������
 // ���� pcmpNoisedEmitMeasures == NULL, �� ���� ������ �� �����������
TComp TAM_2D::fncImitateSingleTargMeasure (const TSingleSign Sign, const double valLambda
	   , TComp * pcmpTrueEmitMeasures, TComp * pcmpNoisedEmitMeasures)
{
  TComp cmpRez0(0.,0.);
  TComp *pcmparrTrueEmitMeasures = (TComp *)malloc( mNumEmitRows * mNumEmitCols * sizeof(TComp));
  memset(pcmparrTrueEmitMeasures, 0, mNumEmitRows * mNumEmitCols * sizeof(TComp));
  TComp *pcmparrNoisedEmitMeasures = (TComp *)malloc( mNumEmitRows * mNumEmitCols * sizeof(TComp));
  memset(pcmparrNoisedEmitMeasures, 0, mNumEmitRows * mNumEmitCols * sizeof(TComp));

	// ���������� ������� �������� � ����������� �����
   TComp cmpOklCoefUs(mOtklCoefUs, 0.);
   // ������������ ������� �������� ������ �����, ����������� �� ���� �����������



	///
	  TSingleSign sign1 =  Sign;
	  sign1.mBet -=  mOtklBetta;  // �������� ���� ����� �������
	  sign1.mEps -=  mOtklEps; // �������� ���� ����� �������
	  TComp cmpSignal(Sign.mAmpl * cos(Sign.mPhase)
			  , Sign.mAmpl * sin(Sign.mPhase));// �������� ������ � ����������� �����
	 TComp cmpNorm ( cos( sign1.mEps) * cos(sign1.mBet ) / ((double)mNumEmitCols)/ ((double)mNumEmitRows), 0.);

	 	for (int m = 0; m < mNumEmitRows; m++)
		{
		for (int l = 0; l < mNumEmitCols; l++)
		{
			double valq = -( mdCol * cos( sign1.mEps ) * sin(sign1.mBet ) * ( (((double)mNumEmitCols)-1.)/ 2. - (double)l )
				 + mdRow * sin( sign1.mEps) * ( (((double)mNumEmitRows)-1.)/ 2. - (double)m ));
		   //	double valq = -( mdCol * cos( sign1.mEps ) * sin(sign1.mBet ) * ( (((double)mNumEmitCols)-1.)/ 2. - (double)l )
			 //	 - mdRow * sin( sign1.mEps) * ( (((double)mNumEmitRows)-1.)/ 2. - (double)m ));

			TComp cmpSdvig = exp_(valq * 2. * M_PI /valLambda);
			pcmparrTrueEmitMeasures[ m *mNumEmitCols + l] =  (((cmpSignal* cmpSdvig)* cmpNorm)*cmpOklCoefUs) ;//+ cmpNoise;
		}
		}

	///

	// ������������ ������� ����������� ��������
		for (int m = 0; m < mNumEmitRows; m++)
		{
		for (int l = 0; l < mNumEmitCols; l++)
		{
			TComp cmpNoise (0.,0.);
			cmpNoise.m_Re = getGauss(0., mSigEmitNoise ) ;
			cmpNoise.m_Im = getGauss(0., mSigEmitNoise ) ;
			pcmparrNoisedEmitMeasures[ m *mNumEmitCols + l] =   pcmparrTrueEmitMeasures[ m *mNumEmitCols + l] + cmpNoise;
			cmpRez0 +=  pcmparrNoisedEmitMeasures[ m *mNumEmitCols + l];
		}
		}







   if  (pcmpTrueEmitMeasures)
   {
	   memcpy(pcmpTrueEmitMeasures, pcmparrTrueEmitMeasures,mNumEmitRows * mNumEmitCols * sizeof(TComp));
   }
   if  (pcmpNoisedEmitMeasures)
   {
	   memcpy(pcmpNoisedEmitMeasures, pcmparrNoisedEmitMeasures,mNumEmitRows * mNumEmitCols * sizeof(TComp));
   }
   free(pcmparrNoisedEmitMeasures);
   free(pcmparrTrueEmitMeasures);
   return cmpRez0;

}

 /*// ������� ��������� ���������
 double  TAM_2D::fncFSource (const double valTetta)
 {
	 return sqrt(fabs(cos(valTetta)));
 }
  //---------------------------------------------------------
 // ������� �����������  ��������� ���������
 double  TAM_2D::fnc_dFSource_po_dTetta (const double valTetta)
 {
	 if ((fabs(valTetta) -M_PI/ 2.) <=0.)
	 {
	  return - sin (valTetta)/fncFSource (valTetta)/2.;
	 }
	 return sin (valTetta)/fncFSource (valTetta)/2.;
 }
  //---------------------------------------------------------

 // �������  ����������� ��������� ������
 double  TAM_2D::fnc_dFAM_2D_po_dTet (const double tet)
 {
  double valk = 2. * M_PI / mLambda;
double alf = valk * m_d * ((double)m_n)/ 2.;
double bet = valk * m_d / 2.;

 if (fabs(tet)< 0.0000001) return 0.;


 double val_F1 = fncFSource (tet);
double val_dF1 = fnc_dFSource_po_dTetta (tet) ;

double f31 = fnc_F3( alf,  tet)  ;
double f32 = fnc_F3( bet,  tet)  ;


double val_F21 = f31/ f32;


double a1 =  fnc_dF3_po_tet(alf,   tet) ;
double a2 =  fnc_dF3_po_tet(bet,   tet) ;


double val_dF21_po_dTet =  fnc_F4(f31 ,f32, a1, a2) ;


return ( val_dF1 *  val_F21  + val_F1 * val_dF21_po_dTet )/ ((double ) m_n);
 }
  //---------------------------------------------------------
double  TAM_2D::fncFAM_2D (const double tet)
{

  if (fabs(tet)< 0.0000001) return 1.;
    double valk = 2. * M_PI / mLambda;
double alf = valk * m_d * ((double)m_n)/ 2.;
double bet = valk * m_d / 2.;
  double temp = sqrt(cos(tet))* sin( alf * sin(tet) )
		   /(((double)m_n)* sin( bet * sin(tet) )) ;

  return temp  ;
}
 //
  //---------------------------------------------------------
 double TAM_2D::fnc_dF3_po_tet(const double val_gam, const double tet)
{
   return cos(val_gam * sin(tet) ) * val_gam * cos(tet);
}
  //---------------------------------------------------------
/// F3 = sin(gam(sin(tet)-mu) -
double TAM_2D::fnc_F3(const double val_gam, const double tet)
{
   return sin(val_gam * sin(tet) );
}

double TAM_2D::fnc_F4(const double v,const  double u, const double v1,const  double u1)
{
   return (v1 * u - u1 * v)/ u / u;
}

//-----------------------------------------------------------------

double TAM_2D::findDiagrWidth()
{
 double tet0 = mLambda/m_d/ m_n * 0.75;

 int i =0;
 double a = sqrt(2.)/2.;
 for ( i =0; i < 10; i++)
 {
  double del = -(fncFAM_2D ( tet0) - a)/ fnc_dFAM_2D_po_dTet ( tet0);
  tet0 += del;
  if (fabs(del) < 0.0000001) break;

 }
 return tet0;
}


void TAM_2D::createDiagrGraphs(wchar_t *wchFoldName1 )
{
	wchar_t wchFoldName[300] ={0};
	wcscpy(  wchFoldName,  wchFoldName1);
	wcscat(wchFoldName, L"\\");
	double step = M_PI / 3000./ 10.;

	const int nBuffRows = 1500 *2;
	const int nBuffCols =3;
	double  *parrBuff = new double [nBuffRows  * nBuffCols] ;
	memset(parrBuff, 0, nBuffRows * nBuffCols * sizeof(double));

	const int lenName = 30;
	wchar_t wcharrFileNames [150] ={0};
	wcscpy( wcharrFileNames, L"Tetta");
	wcscpy( &wcharrFileNames[30], L"F");
	wcscpy( &wcharrFileNames[60], L"SinxF");

	double valWidthDgr = findDiagrWidth() ;

	for (int i=0 ; i < nBuffRows; i++)
	{
		parrBuff[ i * nBuffCols] = ((double) (-nBuffRows/2 +i)) / 10.;
		double tet = step * ((double) (-nBuffRows/2 +i));
		parrBuff[ i * nBuffCols + 1]= fncFAM_2D( tet)   ;
		parrBuff[ i * nBuffCols + 2]= fncDiagrSimple(valWidthDgr, tet * 1000.)   ;

	}

	double scaley = 100.;
	for (int i=1; i < nBuffCols; i++)
	{
		TYrWriteShapeFile::WriteOneReport(                 wchFoldName  // ���� � �����
		, parrBuff // ������ � ����������� - ������� nBuffRows x nBuffCols
		, nBuffCols // - �-�� ���������� � ������� ��������� ���������� � ������
		, nBuffRows //  - �-�� �����
		,wcharrFileNames //������� � �������� ���������� - ������� nBuffCols x lenName
		,lenName // ������������ ����� ����� ����������
		,0  // ����� ���������� �� ��� X
		,i  // ����� ���������� �� ��� Y
		,1. //  ������� �� ��� X
		,scaley // ������� �� ��� Y
		) ;
	}

	delete parrBuff;
	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  wchFoldName);
	wcscat(wchAxesFileName, L"AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-parrBuff[ (nBuffRows-1) * nBuffCols]-100000.*step, parrBuff[ (nBuffRows-1) * nBuffCols]+100000.*step
	,0,scaley * 1.1, 50000.*step) ;
	TURPointXY pPoints[2];
	pPoints[0] =  TURPointXY (-100., scaley * sqrt(2.)/2.);
	pPoints[1] =  TURPointXY (100., scaley * sqrt(2.)/2.);
	TURPolyLine pln( pPoints,2) ;

	wcscpy(  wchAxesFileName,  wchFoldName);
	wcscat(wchAxesFileName, L"Line0.shp");
	TURPolyLine::WriteSetSHPFiles(wchAxesFileName,&pln, 1) ;
}

 // ������������ ������ ������ ��������� ����� sinx/x
 double TAM_2D::findDiagrWidthApprox()
{
 double vala = 2. * M_PI / mLambda * m_d * ((double)m_n)/ 2.;
 double tet1 = TET0707 / vala;
 return tet1;
}

 // ������������ ������  ��������� ����� sinx/x
double  TAM_2D::fncFAM_2DApprox (const double valTetta)
{
	double valWidth = findDiagrWidthApprox();
	return fncDiagrSimple(valWidth, valTetta );

}

 // ������������ ������ �����������  ��������� ����� sinx/x
double  TAM_2D::fnc_dFAM_2DApprox_po_dTet (const double tet)
{
	double valWidth = findDiagrWidthApprox();
	return fncDerivDiagrSimple(valWidth, tet );

}

 // ������������ ������ ������ �����������  ��������� ����� sinx/x
double TAM_2D::fnc_d2FAM_2DApprox_po_dTet2(const double tet)
{
double valWidth = findDiagrWidthApprox();
return fncDeriv2DiagrSimple(valWidth, tet );
} */
#pragma package(smart_init)

