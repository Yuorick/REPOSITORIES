//---------------------------------------------------------------------------


#pragma hdrstop
#include <vcl.h>
#include "Far_2D.h"
#include <dir.h>


#pragma hdrstop
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "Comp.h"
#include <stdio.h>
#include "YrWriteShapeFile.h"
#include "URPolyLine.h"
#include "UrPointXY.h"
#include "DiagrSinX.h"
#include "Diagrams.h"
#include "SingleSign.h"
#include "EtalonSign.h"
#include "Gauss.h"
#include "Far.h"
#include "Equations.h"
#include "Diagrams.h"
#include "MatrixProccess.h"


  class TFar;
//extern const double TET0707;
__fastcall TFar_2D::~TFar_2D()
{
	if(mpAm2D) free(mpAm2D) ;
	mpAm2D = NULL ;

	if(mpbarrAM)
	{
	free(mpbarrAM);
	}
	mpbarrAM = NULL;

}

__fastcall  TFar_2D::TFar_2D()
{
	// �-�� ��  �� �����������
	 mNumAMCols = 0;
	// �-�� �� �� ���������
	 mNumAMRows =0;
	// ���������� ����� �� �� �����������
	 mdAMCol = 0.;
	// ���������� ����� �� �� ���������
	 mdAMRow = 0.;

	mLambda = 0.;
	mpAm2D = NULL;
	mpbarrAM = NULL;
	mDistSKZ = 0.;

}
// ����������� �����������
__fastcall  TFar_2D::TFar_2D (const TFar_2D &R2)
 {
		mNumAMCols = R2.mNumAMCols ;
		mNumAMRows = R2.mNumAMRows;
		mdAMCol = R2.mdAMCol ;
		mdAMRow = R2.mdAMRow;
		mLambda = R2.mLambda;
		mDistSKZ = R2.mDistSKZ;

		//if (mpAm2D != NULL)
	   //	{
	   //	free(mpAm2D);
	   //	mpAm2D = NULL;
	   //	}
		if(R2.mpAm2D != NULL)
		{
		mpAm2D = (TAM_2D*)malloc(mNumAMCols*mNumAMRows* sizeof(TAM_2D)) ;
		if(mpAm2D == NULL)
		{
		ShowMessage(L"Not memory for mpAm2D") ;
		Abort() ;
		}
		memcpy( mpAm2D,R2.mpAm2D, mNumAMCols*mNumAMRows* sizeof(TAM_2D)) ;
		}
		///
	   //	if (mpbarrAM != NULL)
		////{
	 	//free(mpbarrAM);
	   //	mpbarrAM = NULL;
		//}
		if(R2.mpbarrAM != NULL)
		{
		mpbarrAM = (bool*)malloc(mNumAMCols*mNumAMRows* sizeof(bool)) ;
		if(mpbarrAM == NULL)
		{
		ShowMessage(L"Not memory for mpAm2D") ;
		Abort() ;
		}
		memcpy( mpbarrAM,R2.mpbarrAM, mNumAMCols*mNumAMRows* sizeof(bool)) ;
	}

 }



// �������� ������������
  TFar_2D &TFar_2D::operator=(const TFar_2D  &R2)
 {
		mNumAMCols = R2.mNumAMCols ;
		mNumAMRows = R2.mNumAMRows;
		mdAMCol = R2.mdAMCol ;
		mdAMRow = R2.mdAMRow;
		mLambda = R2.mLambda;
		mDistSKZ = R2.mDistSKZ;

	   //	if (mpAm2D != NULL)
	   //	{
	   ///////////	free(mpAm2D);
	   ////	mpAm2D = NULL;
	   //	}
		if(R2.mpAm2D != NULL)
		{
		mpAm2D = (TAM_2D*)malloc(mNumAMCols*mNumAMRows* sizeof(TAM_2D)) ;
		if(mpAm2D == NULL)
		{
		ShowMessage(L"Not memory for mpAm2D") ;
		Abort() ;
		}
		memcpy( mpAm2D,R2.mpAm2D, mNumAMCols*mNumAMRows* sizeof(TAM_2D)) ;
		}
		///
	   //	if (mpbarrAM != NULL)
	   //	{
	   //	free(mpbarrAM);
	   //	mpbarrAM = NULL;
		//}
		if(R2.mpbarrAM != NULL)
		{
		mpbarrAM = (bool*)malloc(mNumAMCols*mNumAMRows* sizeof(bool)) ;
		if(mpbarrAM == NULL)
		{
		ShowMessage(L"Not memory for mpAm2D") ;
		Abort() ;
		}
		memcpy( mpbarrAM,R2.mpbarrAM, mNumAMCols*mNumAMRows* sizeof(bool)) ;
	}

	 return *this ;
}

 // ����� ������
 __fastcall TFar_2D::TFar_2D(int NumAMCols, int NumAMRows, double Lambda ,double dAMCol
	 ,double dAMRow , TAM_2D *pAm2D )
 {
	 mNumAMCols= NumAMCols ;
	 mNumAMRows = NumAMRows;
	 mdAMCol = dAMCol ;
	 mdAMRow = dAMRow;
	 mLambda = Lambda;
	 if (mpAm2D != NULL)
	{
	free(mpAm2D);
	mpAm2D = NULL;
	}
	if(pAm2D != NULL)
	{
		mpAm2D = (TAM_2D*)malloc(mNumAMCols*mNumAMRows* sizeof(TAM_2D)) ;
		if(mpAm2D == NULL) 	   {
		ShowMessage(L"Not memory for mpAm2D") ;
		Abort() ;
	}
	memcpy( mpAm2D,pAm2D, mNumAMCols*mNumAMRows* sizeof(TAM_2D)) ;
	}


 }

 // ����� ������
 __fastcall TFar_2D::TFar_2D(int NumAMCols, int NumAMRows, double Lambda ,double dAMCol
	 ,double dAMRow , TAM_2D am2D )
 {
	 mNumAMCols= NumAMCols ;
	 mNumAMRows = NumAMRows;
	 mdAMCol = dAMCol ;
	 mdAMRow = dAMRow;
	 mLambda = Lambda;
	 if (mpAm2D != NULL)
	{
	free(mpAm2D);
	mpAm2D = NULL;
	}
	mpAm2D = (TAM_2D*)malloc(mNumAMCols*mNumAMRows* sizeof(TAM_2D)) ;
	if(mpAm2D == NULL)
	{
	ShowMessage(L"Not memory for mpAm2D") ;
		Abort() ;
	}
	for (int i =0; i < NumAMCols * NumAMRows; i++)
	{
		mpAm2D[i] =  am2D;
	}

 }

	// ����� ������
 __fastcall TFar_2D::TFar_2D(int NumAMCols, int NumAMRows, double Lambda ,double dAMCol
	 ,double dAMRow , TAM_2D am2D, bool *pbarrAM , const double DistSKZ)
 {
	 mNumAMCols= NumAMCols ;
	 mNumAMRows = NumAMRows;
	 mdAMCol = dAMCol ;
	 mdAMRow = dAMRow;
	 mLambda = Lambda;
	 mDistSKZ = DistSKZ;
	if (mpAm2D != NULL)
	{
	//free(mpAm2D);
	mpAm2D = NULL;
	}
	mpAm2D = (TAM_2D*)malloc(mNumAMCols*mNumAMRows* sizeof(TAM_2D)) ;
	if(mpAm2D == NULL)
	{
	ShowMessage(L"Not memory for mpAm2D") ;
		Abort() ;
	}
	for (int i =0; i < NumAMCols * NumAMRows; i++)
	{
	  mpAm2D[i] =  am2D;
	}
//////
	if (mpbarrAM != NULL)
	{
	//free(mpAm2D);
	mpbarrAM = NULL;
	}
	mpbarrAM = (bool*)malloc(mNumAMCols*mNumAMRows* sizeof(bool)) ;
	if(mpbarrAM == NULL)
	{
	ShowMessage(L"Not memory for mpbarrAM") ;
		Abort() ;
	}
	for (int i =0; i < NumAMCols * NumAMRows; i++)
	{
	  mpbarrAM[i] =  pbarrAM[i];
	}

 }

 	// ����� ������
 __fastcall TFar_2D::TFar_2D(int NumAMCols, int NumAMRows, double Lambda ,double dAMCol
	 ,double dAMRow , TAM_2D am2D, bool *pbarrAM )
 {
	 mNumAMCols= NumAMCols ;
	 mNumAMRows = NumAMRows;
	 mdAMCol = dAMCol ;
	 mdAMRow = dAMRow;
	 mLambda = Lambda;
	if (mpAm2D != NULL)
	{
	//free(mpAm2D);
	mpAm2D = NULL;
	}
	mpAm2D = (TAM_2D*)malloc(mNumAMCols*mNumAMRows* sizeof(TAM_2D)) ;
	if(mpAm2D == NULL)
	{
	ShowMessage(L"Not memory for mpAm2D") ;
		Abort() ;
	}
	for (int i =0; i < NumAMCols * NumAMRows; i++)
	{
	  mpAm2D[i] =  am2D;
	}
//////
	if (mpbarrAM != NULL)
	{
	//free(mpAm2D);
	mpbarrAM = NULL;
	}
	mpbarrAM = (bool*)malloc(mNumAMCols*mNumAMRows* sizeof(bool)) ;
	if(mpbarrAM == NULL)
	{
	ShowMessage(L"Not memory for mpbarrAM") ;
		Abort() ;
	}
	for (int i =0; i < NumAMCols * NumAMRows; i++)
	{
	  mpbarrAM[i] =  pbarrAM[i];
	}

 }




 TComp TFar_2D::fncIdealDiagr (const double valBetta,const double valEps)
 {
  TComp cmpRez(0.,0.), cmpRez0(0.,0.);
  for (int m = 0; m < mNumAMRows; m++)
  {
	for (int l = 0; l < mNumAMCols; l++)
	{
		double valq = -( mdAMCol * cos( valEps) * sin(valBetta ) * ( (((double)mNumAMCols)-1.)/ 2. - (double)l )
			 + mdAMRow * sin( valEps) * ( (((double)mNumAMRows)-1.)/ 2. - (double)m ));
		cmpRez += exp_(valq * 2. * M_PI / mLambda);
	}
  }
  TComp cmpTemp ( 1. / ((double)mNumAMCols)/ ((double)mNumAMRows), 0.);
  cmpRez0 = mpAm2D[0].fncIdealDiagr ( valBetta, valEps, mLambda)* cmpTemp  * cmpRez;

	 return cmpRez0    ;
 }



 // ���������� �������� ��������� � ������������ ���������
 void TFar_2D::createVertIdealDiagrGraphs( wchar_t *wchFoldName1 )
{
	wchar_t wchFoldName[300] ={0};
	wcscpy(  wchFoldName,  wchFoldName1);
	wcscat(wchFoldName, L"\\");


	const int nBuffRows = 1500 *2;
	const int nBuffCols =2;
	double step = M_PI/ 10./ ((double) nBuffRows);
	double  *parrBuff = new double [nBuffRows  * nBuffCols] ;
	memset(parrBuff, 0, nBuffRows * nBuffCols * sizeof(double));

	const int lenName = 30;
	wchar_t wcharrFileNames [150] ={0};
	wcscpy( wcharrFileNames, L"Tetta");
	wcscpy( &wcharrFileNames[30], L"F_Vert");



	const double valBetta = 0.;
	double tet0 = - ((double)nBuffRows )/2. * step;
	for (int i=0 ; i < nBuffRows; i++)
	{
		double tet = tet0 + ((double)i) * step;
		parrBuff[ i * nBuffCols] = tet;
		TComp cmpTemp = fncIdealDiagr (valBetta,tet);

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
		,1000. //  ������� �� ��� X
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

// ������������ ������� ������� �������� �������
// ��������� ����
// INPUT :
// valAmpl - ��������� ������� ����
// valPhase  - ���� ������� ���� � ������ ���
// valBetta, valEps  - ���� ����
// OUTPUT:
// pcmpAmMeasures - ������ ��������� ��

// �� ��������!!!! 17.07.2019
 void TFar_2D::fncImitateMeasuresArray (const double valAmpl,const double valPhase
 ,const double valBetta,const double valEps, TComp *pcmpAmMeasures)
 {
  TComp cmpRez(0.,0.), cmpRez0(0.,0.);
  TComp *pcmpEmitMeasures  = (TComp *)malloc(mpAm2D[0].mNumEmitCols * mpAm2D[0].mNumEmitRows * sizeof(TComp));
  for (int m = 0; m < mNumAMRows; m++)
  {
	for (int l = 0; l < mNumAMCols; l++)
	{
		double valq = -( mdAMCol * cos( valEps) * sin(valBetta ) * ( (((double)mNumAMCols)-1.)/ 2. - (double)l )
			 + mdAMRow * sin( valEps) * ( (((double)mNumAMRows)-1.)/ 2. - (double)m ));


		double valPhSdvig =  valq * 2. * M_PI / mLambda;

	   pcmpAmMeasures[m * mNumAMCols + l]=	mpAm2D[m * mNumAMCols + l].fncImitateSingleTargMeasure (valBetta, valEps , mLambda
	 , valAmpl,  valPhase + valPhSdvig,  pcmpEmitMeasures) ;
	}
  }
 }

 // �������������!!!!
 // ������������ ������� ������� �������� �������
// ��������� ����
// INPUT :
// Sign - ������ ����, ����������� � ����� ���
// OUTPUT:
//pcmpNoisedAmMeasures- ����������� ������ ��������� ��
// pcmpTrueAmMeasures - �������� (��������������) ������ ��������� ��
// ���� pcmpTrueAmMeasures==NULL, �� ���� ������ �� �����������
 void TFar_2D::fncImitateMeasuresArray (const TSingleSign Sign, TComp *pcmpTrueAmMeasures
 , TComp *pcmpNoisedAmMeasures)
 {
  TComp cmpRez(0.,0.), cmpRez0(0.,0.);
  TComp *pcmpTrueEmitMeasures  = (TComp *)malloc(mpAm2D[0].mNumEmitCols * mpAm2D[0].mNumEmitRows * sizeof(TComp));
  TComp *pcmpNoisedEmitMeasures  = (TComp *)malloc(mpAm2D[0].mNumEmitCols * mpAm2D[0].mNumEmitRows * sizeof(TComp));
  TComp *pcmparrE = (TComp *)malloc(mpAm2D[0].mNumEmitCols * mpAm2D[0].mNumEmitRows * sizeof(TComp));
   for (int i=0; i < mpAm2D[0].mNumEmitCols * mpAm2D[0].mNumEmitRows; i++)
   {
	pcmparrE[i] = TComp(1.,0.);
   }
  // ������� � �)����� �����
  // cos(eps)sin(bet) ,  cos(eps)cos(bet)  , sin(eps)
  // ���������� ������ �� � ������� m, l  � ����
  //Xm,l =(Ncols-1)dcol/2 - l* dcol, Ym,l = 0, Zm,l = (Nrows-1)*drow/2 - m* drow
  for (int m = 0; m < mNumAMRows; m++)
  {
	for (int l = 0; l < mNumAMCols; l++)
	{
	 // double Xml = (((double)mNumAMCols)-1.)/ 2. *mdAMCol - ((double)l) * mdAMCol;
	  //double Zml = (((double)mNumAMRows)-1.)/ 2. *mdAMRow - ((double)m) * mdAMRow;
	 // double valq =  -( Xml *cos( Sign.mEps)* sin(Sign.mBet ) +  Zml * sin( Sign.mEps));
		double valq = -( mdAMCol * cos( Sign.mEps) * sin(Sign.mBet ) * ( (((double)mNumAMCols)-1.)/ 2. - (double)l )
		 + mdAMRow * sin( Sign.mEps) * ( (((double)mNumAMRows)-1.)/ 2. - (double)m ));

		double valPhSdvig =  valq * 2. * M_PI / mLambda;

	   //	mpAm2D[m * mNumAMCols + l].fncImitateSingleTargMeasure (valBetta, valEps , mLambda
	// , valAmpl,  valPhase + valPhSdvig,  pcmpEmitMeasures) ;
	TSingleSign sign1 =  Sign;
	sign1.mPhase +=  valPhSdvig;
	sign1.mAmpl = sign1.mAmpl / ((double)mNumAMRows)/ ((double)mNumAMCols);
	pcmpNoisedAmMeasures [m * mNumAMCols + l] =
			mpAm2D[m * mNumAMCols + l].fncImitateSingleTargMeasure (sign1, mLambda, pcmpTrueEmitMeasures, pcmpNoisedEmitMeasures) ;
   // ������������ �������������� ������� ��������� ��
	if (pcmpTrueAmMeasures)
	{
		MtrxMultMatrx(pcmpTrueEmitMeasures,1, mpAm2D[0].mNumEmitCols * mpAm2D[0].mNumEmitRows
		   ,pcmparrE, 1 , &pcmpTrueAmMeasures[m * mNumAMCols + l]);
	}

	}
  }

   free(pcmpTrueEmitMeasures );
   free(pcmpNoisedEmitMeasures) ;
   free(pcmparrE);
 }

// ���������� ���� � ������� �������� �������� ��
// INPUT:
// pcmpTrueAmMeasures -  ������ �������� ��������
// OUTPUT:
// pcmpNoisedAmMeasures -  ������ ����������� ��������
//
 void TFar_2D::fncAddNoise_To_TrueMeasuresArray (TComp *pcmpTrueAmMeasures, TComp *pcmpNoisedAmMeasures)
 {
  TComp cmpRez(0.,0.), cmpRez0(0.,0.);
  TComp *pcmpTrueEmitMeasures  = (TComp *)malloc(mpAm2D[0].mNumEmitCols * mpAm2D[0].mNumEmitRows * sizeof(TComp));
  TComp *pcmpNoisedEmitMeasures  = (TComp *)malloc(mpAm2D[0].mNumEmitCols * mpAm2D[0].mNumEmitRows * sizeof(TComp));
  TComp *pcmparrTemp  = (TComp *)malloc(mNumAMRows * mNumAMCols * sizeof(TComp));
  memset(pcmparrTemp, 0, mNumAMRows * mNumAMCols * sizeof(TComp));
  TComp *pcmparrE = (TComp *)malloc(mpAm2D[0].mNumEmitCols * mpAm2D[0].mNumEmitRows * sizeof(TComp));
   for (int i=0; i < mpAm2D[0].mNumEmitCols * mpAm2D[0].mNumEmitRows; i++)
   {
	pcmparrE[i] = TComp(1.,0.);
   }
	TSingleSign sign1(0.,0.,0.,0.);
  for (int m = 0; m < mNumAMRows; m++)
  {
	for (int l = 0; l < mNumAMCols; l++)
	{
		if (!mpbarrAM)
		{
		continue;
		}
	mpAm2D[m * mNumAMCols + l].fncImitateSingleTargMeasure (sign1, mLambda, pcmpTrueEmitMeasures, pcmpNoisedEmitMeasures) ;
	MtrxMultMatrx(pcmpNoisedEmitMeasures,1, mpAm2D[0].mNumEmitCols * mpAm2D[0].mNumEmitRows
	,pcmparrE, 1 , &pcmparrTemp[m * mNumAMCols + l]);
	}

  }

   MtrxSumMatrx(pcmpTrueAmMeasures, pcmparrTemp,1, mNumAMRows * mNumAMCols, pcmpNoisedAmMeasures) ;
   free(pcmpTrueEmitMeasures );
   free(pcmpNoisedEmitMeasures) ;
   free(pcmparrE);
   free(pcmparrTemp);
 }

 // ���������� �������� �������������� ������������ 1-�� � 4-��
 // ��������� �������� ��� 5�10
 //
 void TFar_2D::createNormCoefGraphs_For_5P10( wchar_t *wchFoldName1 )
{

	double valDagrWidth = findDiagrWidth(mpAm2D[0].mdCol,mdAMCol, mpAm2D[0].mNumEmitCols
  , 6,mLambda);
	int nBuffRows =  valDagrWidth * 10000.-1;
	wchar_t wchFoldName[300] ={0};
	wcscpy(  wchFoldName,  wchFoldName1);
	wcscat(wchFoldName, L"\\");



	const int nBuffCols =3;
	double step = 0.0001  ;
	double  *parrBuff = new double [nBuffRows  * nBuffCols] ;
	memset(parrBuff, 0, nBuffRows * nBuffCols * sizeof(double));

	const int lenName = 30;
	wchar_t wcharrFileNames [150] ={0};
	wcscpy( wcharrFileNames, L"Betta");
	wcscpy( &wcharrFileNames[30], L"Modul");
	wcscpy( &wcharrFileNames[60], L"Phase");
	parrBuff[0] = 0.;
	parrBuff[1] = 8./ 6.;
	parrBuff[2] = 0.;



	 double valBetta = 0.;

	for (int i= 1 ; i < nBuffRows ; i++)
	{
		valBetta = step * ((double)i );
		parrBuff[ i * nBuffCols] = valBetta;
	 //	    TComp cmp1(1.,0.);
	  //	double valMux = 2. * M_PI/ mLambda * mdAMCol * sin(valBetta);
	  //	TComp cmpZx = exp_(-valMux);
	  //	TComp cmpZx6 = exp_(6. *valMux);
	  //	TComp cmpZx8 = exp_(8. *valMux);
	  //	TComp cmpTemp1 =  (cmp1 - cmpZx8) / (cmp1 - cmpZx6);
		TComp cmpKn = calcNormCoeffFor_1_and_4_Rows_5P10(valBetta ) ;
	  //	cmpTemp1  * cmpZx ;
		parrBuff[ i * nBuffCols + 1] = cmpKn.modul();
		parrBuff[ i * nBuffCols + 2] = cmpKn.phase();

	}

	double scaley = 1.;
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
		,100. //  ������� �� ��� X
		,scaley // ������� �� ��� Y
		) ;
	}

	delete []parrBuff;
	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  wchFoldName);
	wcscat(wchAxesFileName, L"AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,0, parrBuff[ (nBuffRows-1) * nBuffCols]* 100.* 1.3
	,0,scaley * 1.4* 1.2, 0.15);

}


 TComp TFar_2D::calcNormCoeffFor_1_and_4_Rows_5P10(const double ValBetta )
{
	TComp cmp1(1.,0.);
	double valMux = 2. * M_PI/ mLambda * mdAMCol * sin(ValBetta);
	TComp cmpZx = exp_(-valMux);
	TComp cmpZx6 = exp_(6. *valMux);
	TComp cmpZx8 = exp_(8. *valMux);
	TComp cmpTemp1 =  (cmp1 - cmpZx8) / (cmp1 - cmpZx6);
	TComp cmpKn = cmpTemp1  * cmpZx ;
	return cmpKn;
}



TFar_2D TFar_2D::create_5P10(const double valSigNoise, const double valLambda )
{
	TAM_2D Am2D (8, 16, 2.,  2.,valSigNoise/ sqrt(28.*8.*16.)) ;
	double dx = 16.3;
	double dz = 32.7;
	 TFar_2D Far5P10(8, 4, valLambda ,dx
	 ,dz, Am2D );
	 Far5P10.mpAm2D[0].mOtklCoefUs = 0.;
	 Far5P10.mpAm2D[7].mOtklCoefUs = 0.;
	 Far5P10.mpAm2D[24].mOtklCoefUs = 0.;
	 Far5P10.mpAm2D[31].mOtklCoefUs = 0.;
	 Far5P10.mpAm2D[0].mSigEmitNoise = 0.;
	 Far5P10.mpAm2D[7].mSigEmitNoise = 0.;
	 Far5P10.mpAm2D[24].mSigEmitNoise = 0.;
	 Far5P10.mpAm2D[31].mSigEmitNoise = 0.;
	 return Far5P10;
}

// ������ ����� ����� ��������� ��������� �������
// INPUT:
// pcmpAmMeasures - ������ ������� ��
// OUTPUT:
// *pEstBet, *pSigBet - ������ �� � ��� ������
// *pEstEps , *pSigEps - ������ �� � ��� ������
//
void TFar_2D::calcEstRSM(TComp *pcmpAmMeasures, double *pEstBet, double *pSigBet
	 , double *pEstEps, double *pSigEps)
{
   calcEstBetRSM(pcmpAmMeasures, pEstBet, pSigBet) ;
   calcEstEpsRSM(pcmpAmMeasures, pEstEps, pSigEps) ;
  // TComp cmpRight = calcRigthDiagr (pcmpAmMeasures);
  // TComp cmpLeft = calcLeftDiagr (pcmpAmMeasures);
 //  TComp cmpUp = calcUpperDiagr (pcmpAmMeasures);
 //  TComp cmpDown = calcDownDiagr (pcmpAmMeasures);

}

void TFar_2D::calcEstBetRSM(TComp *pcmpAmMeasures, double *pEstBet, double *pSigBet)
{
  double valDispRight = 0., valDispLeft = 0.;
  TComp cmpRight = calcRigthDiagr (pcmpAmMeasures, &valDispRight);
  TComp cmpLeft = calcLeftDiagr (pcmpAmMeasures, &valDispLeft);
  TComp cmp_dPel_po_dSUp(0.,0.),cmp_dPel_po_dSDown (0.,0.);
  double valPelFnc = fncPeleng(cmpRight,cmpLeft, &cmp_dPel_po_dSUp, &cmp_dPel_po_dSDown);
  double valDispPelFnc =  cmp_dPel_po_dSUp.modul() * cmp_dPel_po_dSUp.modul() * valDispRight/ 2.
		 + cmp_dPel_po_dSDown.modul() * cmp_dPel_po_dSDown.modul() * valDispLeft/ 2.;


	double valD = mdAMCol * ((double)mNumAMCols) /2.;
	double valCoef =  mLambda/( M_PI * valD );
	*pEstBet = -atan(valPelFnc)* valCoef;
	double valDeriv = 1./ (1. + valPelFnc * valPelFnc );
	*pSigBet = sqrt( valDispPelFnc * valDeriv * valDeriv * valCoef* valCoef);

}
void TFar_2D::calcEstEpsRSM(TComp *pcmpAmMeasures, double *pEstEps, double *pSigEps)
{
  double valDispUp = 0., valDispDown = 0.;
  TComp cmpUp = calcUpperDiagr (pcmpAmMeasures, &valDispUp);
  TComp cmpDown = calcDownDiagr (pcmpAmMeasures, &valDispDown);
  TComp cmp_dPel_po_dSUp(0.,0.),cmp_dPel_po_dSDown (0.,0.);
  double valPelFnc = fncPeleng(cmpUp,cmpDown, &cmp_dPel_po_dSUp, &cmp_dPel_po_dSDown);
  double valDispPelFnc =  cmp_dPel_po_dSUp.modul() * cmp_dPel_po_dSUp.modul() * valDispUp/2.
		 + cmp_dPel_po_dSDown.modul() * cmp_dPel_po_dSDown.modul() * valDispDown/2.;


	double valD = mdAMRow * ((double)mNumAMRows) /2.;
	double valCoef =  mLambda/( M_PI * valD );
	*pEstEps = -atan(valPelFnc)* valCoef;
	double valDeriv = 1./ (1. + valPelFnc * valPelFnc );
	*pSigEps = sqrt( valDispPelFnc * valDeriv * valDeriv * valCoef* valCoef);
}

// ���������� ������ �� ��������� ���� ���������-��������� �������
// � ���������� ������� 1-�� � 4-�� �����  ��� 5�10
void TFar_2D::calcEstEpsRSM_5P10_Correction_1_and_4_Rows(TComp *pcmpAmMeasures
  , const double ValBet, double *pEstEps, double *pSigEps)
{
  TComp cmpCorrectCoeff = calcNormCoeffFor_1_and_4_Rows_5P10( ValBet )  ;

	TFar_2D far2DTemp = *this;
	TComp cmparrAmMeasures [32];
	memcpy( cmparrAmMeasures, pcmpAmMeasures, 32 * sizeof(TComp));
	for (int i =0; i < 8; i++)
	{
	  cmparrAmMeasures[i] *= cmpCorrectCoeff;
	  far2DTemp.mpAm2D[i].mSigEmitNoise *= cmpCorrectCoeff.modul() ;
	  cmparrAmMeasures[24 + i] *= cmpCorrectCoeff;
	  far2DTemp.mpAm2D[24 +  i].mSigEmitNoise *= cmpCorrectCoeff.modul() ;
	}
	far2DTemp.calcEstEpsRSM(cmparrAmMeasures,pEstEps, pSigEps ) ;

}


// ���������� 4-� ���������� ��������
// ��������� ������� ������ ����, ����� �������
// INPUT:
// pcmpAmMeasures -  ������ ��������� ��
// OUTPUT:
// pcmparrQuaterMeasures [4] -   ������ ��������� ���������� ��������
// parrQuaterDisp [4] -  ������ ��������� ������ ��������� ���������� ��������
void TFar_2D::calcQuarterDiagrams(TComp *pcmpAmMeasures,TComp *pcmparrQuaterMeasures, double *parrQuaterDisp)
{
	memset(parrQuaterDisp,0, 4 *sizeof(double));
	memset(pcmparrQuaterMeasures,0, 4 *sizeof(TComp));
	for (int i = 0; i < mNumAMRows /2   ; i++)
	for (int j = 0; j < mNumAMCols / 2; j++)
	{
	  pcmparrQuaterMeasures[0] +=   pcmpAmMeasures[ i * mNumAMCols + j];
	  pcmparrQuaterMeasures[1] +=   pcmpAmMeasures[ i * mNumAMCols + mNumAMCols / 2 + j];
	  pcmparrQuaterMeasures[2] +=   pcmpAmMeasures[ (i + mNumAMRows /2)* mNumAMCols + j];
	  pcmparrQuaterMeasures[3] +=   pcmpAmMeasures[ (i + mNumAMRows /2)* mNumAMCols + mNumAMCols / 2 + j];
	  parrQuaterDisp[0] += mpAm2D[ i * mNumAMCols + j].calcSumDisp();
	  parrQuaterDisp[1] += mpAm2D[ i * mNumAMCols + mNumAMCols / 2 + j].calcSumDisp();
	  parrQuaterDisp[2] += mpAm2D[ (i + mNumAMRows /2)* mNumAMCols + j].calcSumDisp();
	  parrQuaterDisp[3] += mpAm2D[ (i + mNumAMRows /2)* mNumAMCols + mNumAMCols / 2 + j].calcSumDisp();
	}

}


// ������ ��������� �� ����������� �������� ������ ����� ������� � ����� ������ ���������� ��������
TComp TFar_2D::calcRigthDiagr(TComp *pcmpAmMeasures, double *pDisp)
{
  /*	TComp *cmparrTemp = new TComp[mNumAMCols *mNumAMRows];
	memset(cmparrTemp , 0, mNumAMCols *mNumAMRows * sizeof(TComp));

	*pDisp = 0.;
	for (int i = 0; i < mNumAMRows   ; i++)
	for (int j = 0; j < mNumAMCols / 2; j++)
	{
	  cmparrTemp [ i * mNumAMCols + j] = TComp(1.,0.);
	  *pDisp += mpAm2D[ i * mNumAMCols + j].calcSumDisp();
	}

	TComp cmpRez;
	MtrxMultMatrx(pcmpAmMeasures,1, mNumAMCols *mNumAMRows,cmparrTemp,1, &cmpRez);
	delete[]cmparrTemp;  */
	TComp pcmparrQuaterMeasures[4];
	double parrQuaterDisp[4];
	calcQuarterDiagrams(pcmpAmMeasures,pcmparrQuaterMeasures, parrQuaterDisp);
	TComp cmpRez =  pcmparrQuaterMeasures[0] + pcmparrQuaterMeasures[2];
	*pDisp = parrQuaterDisp[0] + parrQuaterDisp[2];
	return cmpRez;
}

TComp TFar_2D::calcLeftDiagr(TComp *pcmpAmMeasures, double *pDisp)
{
  /*	TComp *cmparrTemp = new TComp[mNumAMCols *mNumAMRows];
	memset(cmparrTemp , 0, mNumAMCols *mNumAMRows * sizeof(TComp));
	*pDisp = 0.;
	for (int i = 0; i < mNumAMRows   ; i++)
	for (int j = 0; j < mNumAMCols / 2; j++)
	{
	  cmparrTemp [ i * mNumAMCols + mNumAMCols / 2 + j] = TComp(1.,0.);
	  *pDisp += mpAm2D[i * mNumAMCols + mNumAMCols / 2 + j].calcSumDisp();
	}

	TComp cmpRez;
	MtrxMultMatrx(pcmpAmMeasures,1, mNumAMCols *mNumAMRows,cmparrTemp,1, &cmpRez);
	delete[]cmparrTemp;
	return cmpRez; */
	TComp pcmparrQuaterMeasures[4];
	double parrQuaterDisp[4];
	calcQuarterDiagrams(pcmpAmMeasures,pcmparrQuaterMeasures, parrQuaterDisp);
	TComp cmpRez =  pcmparrQuaterMeasures[1] + pcmparrQuaterMeasures[3];
	*pDisp = parrQuaterDisp[1] + parrQuaterDisp[3];
	return cmpRez;

}


TComp TFar_2D::calcUpperDiagr(TComp *pcmpAmMeasures, double *pDisp)
{
  /*	TComp cmpRez(0.,0.);
	*pDisp = 0.;
	for (int i = 0; i < mNumAMCols *mNumAMRows/ 2; i++)
	{

	  *pDisp += mpAm2D[i ].calcSumDisp();
	  cmpRez += pcmpAmMeasures[i];
	}

	return cmpRez; */
	TComp pcmparrQuaterMeasures[4];
	double parrQuaterDisp[4];
	calcQuarterDiagrams(pcmpAmMeasures,pcmparrQuaterMeasures, parrQuaterDisp);
	TComp cmpRez =  pcmparrQuaterMeasures[0] + pcmparrQuaterMeasures[1];
	*pDisp = parrQuaterDisp[0] + parrQuaterDisp[1];
	return cmpRez;

}


TComp TFar_2D::calcDownDiagr(TComp *pcmpAmMeasures, double *pDisp)
{
  /*	TComp cmpRez(0.,0.);
	*pDisp = 0.;
	for (int i = 0; i < mNumAMCols *mNumAMRows/ 2; i++)
	{

	  *pDisp += mpAm2D[mNumAMCols *mNumAMRows/ 2 + i ].calcSumDisp();
	  cmpRez += pcmpAmMeasures[mNumAMCols *mNumAMRows/ 2 + i];
	}

	return cmpRez; */
	TComp pcmparrQuaterMeasures[4];
	double parrQuaterDisp[4];
	calcQuarterDiagrams(pcmpAmMeasures,pcmparrQuaterMeasures, parrQuaterDisp);
	TComp cmpRez =  pcmparrQuaterMeasures[2] + pcmparrQuaterMeasures[3];
	*pDisp = parrQuaterDisp[2] + parrQuaterDisp[3];
	return cmpRez;

}

// �������������� ������� ���
// INPUT:
// cmpUp, cmpDown -  ������ ������� � ������ ��������
//  OUTPUT
// *pcmp_dPel_po_dSUp , *pcmp_dPel_po_dSDown - ������� ���������� �������� ������� �� ���������� ������� � ������ ���������

 double TFar_2D::fncPeleng(TComp cmpUp,TComp cmpDown, TComp *pcmp_dPel_po_dSUp
   , TComp *pcmp_dPel_po_dSDown)
 {
  TComp cmpRazn = cmpUp -   cmpDown;
  TComp cmpSum  = cmpUp +   cmpDown;
  TComp cmpPel =  cmpRazn/cmpSum;


	TComp cmp2(2.,0.);
	TComp cmpMin2(-2.,0.);

	TComp cmpt0 = cmp2 *cmpDown/ TComp(1000.,0.);
	TComp cmpt2 = cmp2 *cmpUp/ TComp(1000.,0.);
	TComp cmpt1 = (cmpSum*cmpSum)/ TComp(1000.,0.);
	*pcmp_dPel_po_dSUp = cmpt0 / cmpt1 ;
	*pcmp_dPel_po_dSDown =   cmpt2 / cmpt1 ;

  return cmpPel.m_Im;
 }


 void TFar_2D::fncImitateMultiTargMeasure ( TSingleSign SignTarg, TSingleSign *parrsignFalse, int quantFalseSign
	 , TComp * pcmpTrueAmMeasures, TComp * pcmpNoisedAMMeasures)
{
	memset(pcmpTrueAmMeasures, 0, mNumAMCols * mNumAMRows * sizeof(TComp));
	memset(pcmpNoisedAMMeasures, 0, mNumAMCols * mNumAMRows * sizeof(TComp));

	TComp * pcmpAmMeasuresCur = new TComp [mNumAMCols * mNumAMRows];
	memset(pcmpAmMeasuresCur, 0, mNumAMCols * mNumAMRows * sizeof(TComp));

	fncImitateMeasuresArray (SignTarg, pcmpTrueAmMeasures, pcmpAmMeasuresCur);



	TComp * pcmpTrueAmMeasuresCur = new TComp [mNumAMCols * mNumAMRows];
	memset(pcmpTrueAmMeasuresCur, 0, mNumAMCols * mNumAMRows * sizeof(TComp));


	for (int i =0; i < quantFalseSign; i++)
	{
	 fncImitateMeasuresArray (parrsignFalse[i], pcmpTrueAmMeasuresCur , pcmpAmMeasuresCur);
	 MtrxSumMatrx(pcmpTrueAmMeasuresCur, pcmpTrueAmMeasures,1, mNumAMCols * mNumAMRows, pcmpAmMeasuresCur) ;
	 memcpy(pcmpTrueAmMeasures, pcmpAmMeasuresCur, mNumAMCols * mNumAMRows * sizeof(TComp));
	}
	fncAddNoise_To_TrueMeasuresArray (pcmpTrueAmMeasuresCur, pcmpNoisedAMMeasures) ;

}

// ���������� ������ �� ��������� ���� ���������-��������� �������
// � ���������� ������� 1-�� � 4-�� �����  ��� 5�10
bool TFar_2D::calcEstRSM_5P10_Correction_InAccordanceWithModel(TComp *pcmpAmMeasures
  , double *pEstBet, double *pSigBet, double *pEstEps, double *pSigEps)
{
  bool breturn = false;
  double valEstBet, valSigBet, valEstEps,  valSigEps
	 ,valEstBet0, valSigBet0, valEstEps0,  valSigEps0;
  calcEstRSM(pcmpAmMeasures, &valEstBet, &valSigBet
	 ,&valEstEps, &valSigEps) ;
 // ������������ �������
 int i =0;
  for (i =0; i < 100; i++)
  {
  calcEstEpsRSM_5P10_Correction_1_and_4_Rows( pcmpAmMeasures
  , valEstBet, &valEstEps0, &valSigEps0) ;
  valEstBet0 = valEstBet;
  calcEstBetRSM_5P10_StepCorrection_Left_and_Right_Diagr( pcmpAmMeasures
  , valEstEps0,  &valEstBet0, &valSigBet0) ;
  if ( (fabs(valEstBet0- valEstBet) < 0.0001)
	&& (fabs(valEstEps0- valEstEps) < 0.0001))
   {
	 breturn = true;
	 *pEstBet = valEstBet0;
	 *pSigBet = valSigBet0;
	 *pEstEps = valEstEps0;
	 *pSigEps = valSigEps0;
	 break;
   }
   else
   {
	valEstBet = valEstBet0;
	valEstEps = valEstEps0;
   }

  }


 return breturn;
}


// ���� ��� ��������� ������ ��������� ���� ������������ ���  ��� 5�10
//INPUT:
// pcmpAmMeasures - ������ ��������� ��
// ValEstEps -   ������� ������ ���� �����
// OUTPUT:
// pEstBet - ������ ���� ����
// pSigBet -  ��� ���� ����
void TFar_2D::calcEstBetRSM_5P10_StepCorrection_Left_and_Right_Diagr( TComp *pcmpAmMeasures
  , const double ValEstEps,  double *pEstBet, double *pSigBet)
{
  TComp  cmpSumR(0.,0.), cmpSumL(0.,0.);
  double arrc[4] = {3.,4.,4.,3.};
  for (int i = 0; i < 4; i++)
  {
	TComp cmpTemp = exp_(((double)i + 1.) * mdAMRow * sin(ValEstEps) * 2.* M_PI/ mLambda);
	TComp cmpTemp1 = exp_(arrc[i]* mdAMCol * sin(*pEstBet) * 2.* M_PI/ mLambda);
	cmpSumR +=  (TComp(1.,0.) - cmpTemp1 ) * cmpTemp;
	cmpSumL +=  ((TComp(1.,0.) - cmpTemp1 ) * cmpTemp1) * cmpTemp;

  }
  TComp cmpDzittaR =  TComp(1.,0.)/ cmpSumR;

  TComp cmpDzittaL = (exp_(4.* 2.* M_PI * mdAMCol* sin(*pEstBet)/ mLambda) ) / cmpSumL;
   TComp *pcmpAmMeasures0 = new TComp[32];
   TFar_2D Far2DTemp = *this;
   for (int i =0; i < 4; i++)
   for (int j =0; j < 4; j++)
   {
	 pcmpAmMeasures0 [i * 8 + j]     =  pcmpAmMeasures [i * 8 + j] * cmpDzittaR;
	 pcmpAmMeasures0 [i * 8 + 4 + j] =  pcmpAmMeasures [i * 8 + 4 + j] * cmpDzittaL;
	 Far2DTemp.mpAm2D[i * 8 + j].mSigEmitNoise *= cmpDzittaR.modul();
	 Far2DTemp.mpAm2D[i * 8 + 4 + j].mSigEmitNoise *= cmpDzittaL.modul();
   }
   Far2DTemp.calcEstBetRSM(pcmpAmMeasures0 , pEstBet,pSigBet )   ;
   delete []pcmpAmMeasures0 ;
}



// ���������� ������ ����� ��������� ���� ������������  ���������-��������� �������
// ��� ������ ����� ��������� �� ����� �� ��������
// ��������������� ������ ����� �������� ������ ���� �������, ������
// ������ � ����� ���������
// INPUT:
// �mpUp, �mpDown  - ��������� ������� � ������ ��������
// ValDispUp, ValDispDown -  ���������� ��������� ������� � ������ ��������
// �mpRight, �mpLeft  - ��������� ������ � ����� ��������
// ValDispRight, ValDispLeft -  ���������� ��������� ������ � ����� ��������
// OUTPUT:
// *pEstBet, *pSigBet - ������ �� � ��� ������
// *pEstEps , *pSigEps - ������ �� � ��� ������
bool TFar_2D::calcEstIterRSM( const TComp CmpUp	,const TComp CmpDown
	,const double ValDispUp, const double ValDispDown
   ,const TComp CmpRight, const TComp CmpLeft
   ,const double ValDispRight, const double ValDispLeft
  , double *pEstBet, double *pSigBet, double *pXiSq4Bet
  , double *pEstEps, double *pSigEps, double *pXiSq4Eps)
{
 bool breturn = false;
  double valEstBet, valSigBet, valEstEps,  valSigEps
	 ,valEstBet0, valSigBet0, valEstEps0,  valSigEps0;
  //calcEstRSM(pcmpAmMeasures, &valEstBet, &valSigBet
  //	 ,&valEstEps, &valSigEps) ;
	 const double ValVertD  = mdAMRow * ((double)mNumAMRows) /2.;
	 const double ValHorD  = mdAMCol * ((double)mNumAMCols) /2.;
  calcIdealEstRSM( CmpUp, CmpDown, ValDispUp,ValDispDown,  ValVertD, &valEstEps,  &valSigEps, pXiSq4Eps);

  calcIdealEstRSM(CmpRight, CmpLeft, ValDispRight, ValDispLeft,  ValHorD, &valEstBet,  &valSigBet, pXiSq4Bet);
 // ������������ �������
 int i =0;
  for (i =0; i < 100; i++)
  {
  valEstEps0 = valEstEps;
  fncStep_CorrectEps_UpAndDownDiagr_RSM(  CmpUp, CmpDown
  , valEstBet, &valEstEps0, &valSigEps0,pXiSq4Eps) ;
  valEstBet0 = valEstBet;
  fncStep_CorrectBet_RightAndLeftDiagr_RSM(  CmpRight,  CmpLeft
  , valEstEps0,  &valEstBet0, &valSigBet0,pXiSq4Bet) ;

  if ( (fabs(valEstBet0- valEstBet) < 0.0001)
	&& (fabs(valEstEps0- valEstEps) < 0.0001))
   {
	 breturn = true;
	 *pEstBet = valEstBet0;
	 *pSigBet = valSigBet0;
	 *pEstEps = valEstEps0;
	 *pSigEps = valSigEps0;
//  valD = mdAMRow * ((double)mNumAMRows) /2.;   ��� ���� �����
//  valD = mdAMCol * ((double)mNumAMCols) /2.;  ��� ��������� ����
	// ���������� �� �������!!!!
	 break;
   }
   else
   {
	valEstBet = valEstBet0;
	valEstEps = valEstEps0;
   }

  }


 return breturn;
}



// ������ ���� ������������  ��� ����� ��������� ������� � ������ �������� (������ � �����)
// INPUT:
// CmpUp, CmpDown -  ��������� ������� � ������� ��������
// ValDispUp, ValDispDown -  ���������� ��������� ������� � ������ ��������
// ValD -  ���������� ����� �������� �������� ���� � �����. ��������
//  valD = mdAMRow * ((double)mNumAMRows) /2.;   ��� ���� �����
//  valD = mdAMCol * ((double)mNumAMCols) /2.;  ��� ��������� ����
// OUTPUT:
// *pEstEps - ������ ����
// *pSigEps - ��� ������ ����
void TFar_2D::calcIdealEstRSM(const TComp CmpUp,const TComp CmpDown, const double ValDispUp
, const double ValDispDown, const double ValD, double *pEstEps, double *pSigEps, double *pXiSq4Eps)
{

	TComp cmp_dPel_po_dSUp(0.,0.),cmp_dPel_po_dSDown (0.,0.);
	double valPelFnc = fncPeleng(CmpUp,CmpDown, &cmp_dPel_po_dSUp, &cmp_dPel_po_dSDown);
	double valDispPelFnc =  cmp_dPel_po_dSUp.modul() * cmp_dPel_po_dSUp.modul() * ValDispUp
	 + cmp_dPel_po_dSDown.modul() * cmp_dPel_po_dSDown.modul() * ValDispDown;

	double valCoef =  mLambda/( M_PI * ValD );
	*pEstEps = -atan(valPelFnc)* valCoef;
	double valDeriv = 1./ (1. + valPelFnc * valPelFnc );
	*pSigEps = sqrt( valDispPelFnc * valDeriv * valDeriv * valCoef* valCoef);
	*pXiSq4Eps = calcIdealXiSq4_For_RSM( CmpUp, CmpDown,  ValDispUp
	,ValDispDown,  ValD, *pEstEps) ;
}


// ���������� ������� �� ������� 4 ��� ������ ��������� ��������
// �� ����, ����� ������� ������������ � ��� ������ ��������
double TFar_2D::calcIdealXiSq4_For_RSM( TComp CmpUp, TComp CmpDown, const double ValDispUp
, const double ValDispDown, const double ValD, const double ValEstEps)
{
	TComp cmpRotate = exp_(-ValD * sin(ValEstEps)* 2.* M_PI/ mLambda);
	TComp cmpDownZv = cmpRotate * CmpDown;
	TComp cmpRazn = CmpUp -  cmpDownZv ;
	return 2.* cmpRazn.modul() * cmpRazn.modul()/ (ValDispUp + ValDispDown);
}





void TFar_2D::fncStep_CorrectEps_UpAndDownDiagr_RSM (const TComp CmpUp	,const TComp CmpDown
	,const double ValEstBet, double *pvalEstEps0, double *pvalSigEps0, double *pvalXiSq4Eps)
{
  TComp cmp_fUp = fnc_fUp( ValEstBet,*pvalEstEps0);
  TComp cmp_fDown = fnc_fDown( ValEstBet,*pvalEstEps0);;
  TComp  cmpXiUp = TComp(1.,0.)/cmp_fUp;
  TComp  cmpXiDown = exp_( ((double)mNumAMRows)/ 2.*mdAMRow * sin(*pvalEstEps0) * 2. * M_PI/  mLambda )
				/cmp_fDown;
  TComp cmpSZvUp = cmpXiUp * CmpUp;
  TComp cmpSZvDown = cmpXiDown * CmpDown;


  double valDispUp = 0., valDispDown = 0.;
  for (int i =0; i < mNumAMRows/ 2 ; i++)
  for (int j = 0; j < mNumAMCols; j++)
  {
  valDispUp += mpAm2D[ i * mNumAMCols + j].calcSumDisp();
  valDispDown += mpAm2D[ (mNumAMRows/ 2 + i) * mNumAMCols + j].calcSumDisp();

  }
  valDispUp  = valDispUp / cmp_fUp.modul()/ cmp_fUp.modul();
  valDispDown = valDispDown / cmp_fDown.modul()/ cmp_fDown.modul()  ;
  double valD = mdAMRow * ((double)mNumAMRows) /2.;

  calcIdealEstRSM(cmpSZvUp ,cmpSZvDown, valDispUp
	,valDispDown, valD, pvalEstEps0, pvalSigEps0, pvalXiSq4Eps) ;
}

// ���������� ������� �� ������� 4 ��� ������ ���� ����� ��� ����������� ������ ���
// - ����� ��������� ������ ����� �� ��������
// �������� ������ �������� ���������(�������, ������, ������, �����)
// INPUT:
// CmpUp - ����� ������� ���������
// CmpDown - ����� ������ ���������
// ValEstBet - ������ ��
// ValEstEps - ������ ��
// OUTPUT:
// ��4 �������
double TFar_2D::calcRealXi4Sq_Eps_RSM (const TComp CmpUp	,const TComp CmpDown
	,const double ValEstBet,const double ValEstEps)
{
  TComp cmp_fUp = fnc_fUp( ValEstBet,ValEstEps);
  TComp cmp_fDown = fnc_fDown( ValEstBet,ValEstEps);;
  TComp  cmpXiUp = TComp(1.,0.)/cmp_fUp;
 // TComp  cmpXiDown = exp_( ((double)mNumAMRows)/ 2.*mdAMRow * sin(ValEstEps) * 2. * M_PI/  mLambda )
  //				/cmp_fDown;
  TComp  cmpXiDown = TComp(1.,0.)/cmp_fDown;
  TComp cmpSZvUp = cmpXiUp * CmpUp;
  TComp cmpSZvDown = cmpXiDown * CmpDown;


  double valDispUp = 0., valDispDown = 0.;
  for (int i =0; i < mNumAMRows/ 2 ; i++)
  for (int j = 0; j < mNumAMCols; j++)
  {
  valDispUp += mpAm2D[ i * mNumAMCols + j].calcSumDisp();
  valDispDown += mpAm2D[ (mNumAMRows/ 2 + i) * mNumAMCols + j].calcSumDisp();

  }
  valDispUp  = valDispUp / cmp_fUp.modul()/ cmp_fUp.modul();
  valDispDown = valDispDown / cmp_fDown.modul()/ cmp_fDown.modul()  ;

  double valXiSq4Eps = 0.;
 // calcIdealEstRSM(cmpSZvUp ,cmpSZvDown, valDispUp
  //	,valDispDown, valD, pvalEstEps0, pvalSigEps0, &valXiSq4Eps) ;
  TComp cmpDel =  cmpSZvUp - cmpSZvDown;
  valXiSq4Eps = cmpDel.modul() * cmpDel.modul()/ ( valDispUp + valDispDown);
 //valXiSq4Eps =  calcIdealXiSq4_For_RSM( cmpSZvUp ,cmpSZvDown, valDispUp
   //	,valDispDown, valD,   ValEstEps) ;
	return 2. *valXiSq4Eps;
}


void TFar_2D::fncStep_CorrectBet_RightAndLeftDiagr_RSM(  const TComp CmpRight, const TComp CmpLeft
  ,const double ValEstEps, double *pvalEstBet0, double *pvalSigBet0, double *pvalXiSq4Bet)
{
  TComp cmp_fRight = fnc_fRight( *pvalEstBet0,ValEstEps);
  TComp cmp_fLeft = fnc_fLeft (*pvalEstBet0,ValEstEps);
  TComp  cmpXiRight = TComp(1.,0.)/cmp_fRight ;
  TComp  cmpXiLeft = exp_( ((double)mNumAMCols)/ 2.*mdAMCol * sin(*pvalEstBet0) * 2. * M_PI/  mLambda )
				/cmp_fLeft;
  TComp cmpSZvRight = cmpXiRight * CmpRight;
  TComp cmpSZvLeft = cmpXiLeft * CmpLeft;


  double valDispRight = 0., valDispLeft = 0.;
  for (int i =0; i < mNumAMRows ; i++)
  for (int j = 0; j < mNumAMCols/ 2; j++)
  {
   valDispRight += mpAm2D[ i * mNumAMCols + j].calcSumDisp();
   valDispLeft += mpAm2D[ i * mNumAMCols + mNumAMCols/ 2 + j].calcSumDisp();

  }
  valDispRight = valDispRight / cmp_fRight.modul()/ cmp_fRight.modul();
  valDispLeft = valDispLeft / cmp_fLeft.modul()/ cmp_fLeft.modul()  ;
  double valD = mdAMCol * ((double)mNumAMCols) /2.;

  calcIdealEstRSM(cmpSZvRight ,cmpSZvLeft, valDispRight
	,valDispLeft, valD, pvalEstBet0, pvalSigBet0, pvalXiSq4Bet) ;

}



// ���������� ������� �� ������� 4 ��� ������ ��������� ���� ��� ����������� ������ ���
// - ����� ��������� ������ ����� �� ��������
// �������� ������ �������� ���������(�������, ������, ������, �����)
// INPUT:
// CmpUp - ����� ������� ���������
// CmpDown - ����� ������ ���������
// ValEstBet - ������ ��
// ValEstEps - ������ ��
// OUTPUT:
// ��4 �������
double TFar_2D::calcRealXi4Sq_Bet_RSM (const TComp CmpRight	,const TComp CmpLeft
	,const double ValEstBet,const double ValEstEps)
{
  TComp cmp_fRight = fnc_fRight( ValEstBet,ValEstEps);
  TComp cmp_fLeft = fnc_fLeft( ValEstBet,ValEstEps);;
  TComp  cmpXiRight = TComp(1.,0.)/cmp_fRight;
  TComp  cmpXiLeft = TComp(1.,0.)/cmp_fLeft;;
  TComp cmpSZvRight = cmpXiRight * CmpRight;
  TComp cmpSZvLeft = cmpXiLeft * CmpLeft;


  double valDispRight= 0., valDispLeft = 0.;
  for (int i =0; i < mNumAMRows/ 2 ; i++)
  for (int j = 0; j < mNumAMCols; j++)
  {
  valDispRight += mpAm2D[ i * mNumAMCols + j].calcSumDisp();
  valDispLeft += mpAm2D[ (mNumAMRows/ 2 + i) * mNumAMCols + j].calcSumDisp();

  }
  valDispRight  = valDispRight / cmp_fRight.modul()/ cmp_fRight.modul();
  valDispLeft = valDispLeft / cmp_fLeft.modul()/ cmp_fLeft.modul()  ;

  double valXiSq4Bet = 0.;

  TComp cmpDel =  cmpSZvRight - cmpSZvLeft;
  valXiSq4Bet = cmpDel.modul() * cmpDel.modul()/ ( valDispRight +valDispLeft );
  return 2. *valXiSq4Bet;
}

TComp TFar_2D::fnc_fUp( const double ValEstBet, const double ValEstEps)
{
  TComp cmpSum(0.,0.);
  for (int i =0; i < mNumAMRows/ 2 ; i++)
  for (int j = 0; j < mNumAMCols; j++)
  {
  double temp = fnc_qml(ValEstBet,  ValEstEps, i, j);
  cmpSum += TComp(mpAm2D[ i * mNumAMCols + j].mOtklCoefUs , 0.)
	 * exp_(temp * 2. * M_PI / mLambda);

  }
	return cmpSum;
}
TComp TFar_2D::fnc_fDown( const double ValEstBet, const double ValEstEps)
{
  TComp cmpSum(0.,0.);
  for (int i =0; i < mNumAMRows/ 2 ; i++)
  for (int j = 0; j < mNumAMCols; j++)
  {
  cmpSum += TComp(mpAm2D[ (mNumAMRows/ 2 + i) * mNumAMCols + j].mOtklCoefUs , 0.)
	 * exp_(fnc_qml(ValEstBet,  ValEstEps, mNumAMRows/ 2 + i, j) * 2. * M_PI / mLambda);

  }
	return cmpSum;
}


TComp TFar_2D::fnc_fRight( const double ValEstBet, const double ValEstEps)
{
  TComp cmpSum(0.,0.);
  for (int i =0; i < mNumAMRows ; i++)
  for (int j = 0; j < mNumAMCols /2; j++)
  {
  cmpSum += TComp(mpAm2D[ i * mNumAMCols + j].mOtklCoefUs , 0.)
	 * exp_(fnc_qml(ValEstBet,  ValEstEps, i, j) * 2. * M_PI / mLambda);

  }
	return cmpSum;
}
TComp TFar_2D::fnc_fLeft( const double ValEstBet, const double ValEstEps)
{
  TComp cmpSum(0.,0.);
  for (int i =0; i < mNumAMRows ; i++)
  for (int j = 0; j < mNumAMCols /2; j++)
  {
  cmpSum += TComp(mpAm2D[ i * mNumAMCols +  mNumAMCols /2 + j].mOtklCoefUs , 0.)
	 * exp_(fnc_qml(ValEstBet,  ValEstEps, i,  mNumAMCols /2 + j) * 2. * M_PI / mLambda);

  }
	return cmpSum;
}

double  TFar_2D::fnc_qml( const double ValEstBet, const double ValEstEps, const int Num_m, const int Num_l)
{
  double valreturn = -(mdAMCol * ((double)(mNumAMCols + 1.)/ 2. - (double) Num_l) * sin(ValEstBet)
	   + mdAMRow * ((double)(mNumAMRows + 1.)/ 2. - (double) Num_m)* sin(ValEstEps) );
  return valreturn ;
}


 // ���������� �������� ��������� � ������������ ���������
 // ��� ������ �������:
 // 1.���������� �������� ��������� ������ ��������� �� � ��  ��� ������������� �� = 0, 2, 4, 6,8,10, -2, -4, -6,-8,-10 ����
 // � ����������� �� ���� �����  � �������� -10 +10 ����
 // 2.���������� �������� ��������� ������ ��������� �� � ��  ��� ������������� �� = 0, 2, 4, 6,8,10, -2, -4, -6,-8,-10 ����
 // � ����������� �� ��  � �������� -10 +10 ����
 // INPUT:
 // ValBet - ��
 // ValEps - ��
 // wchFoldName1 - ����  � ����� � �������

 void TFar_2D::createSystemErrorGraphs_RSM( wchar_t *wchFoldName1 )
{
	wchar_t wchFoldName[300] ={0}, wchFoldNameBet[300] ={0}, wchFoldNameEps[300] ={0};
	wcscpy(  wchFoldName,  wchFoldName1);
	wcscat(wchFoldName, L"\\");
	wchar_t wstr[10] = {0};
	for (int i = 0; i < 21; i++)
	{

	  int num = -10 +i;
	  wcscpy(  wchFoldNameBet,  wchFoldName1);
	  wcscpy(  wchFoldNameEps,  wchFoldName1);
	  wcscat(wchFoldNameBet, L"\\FixedBet=");
	  wcscat(wchFoldNameBet, wstr);
	  wcscat(wchFoldNameEps, L"\\FixedEps=");
	  wcscat(wchFoldNameEps, wstr);

	  _wmkdir(wchFoldNameBet);
	  _wmkdir(wchFoldNameEps);

	  const double   ValEps = ((double)num)*0.001;
	  const double   ValBet = ((double)num)*0.001;
	  createSystemErrorGraphsFixedBetta_RSM( ValBet,wchFoldNameBet ) ;
	  createSystemErrorGraphsFixedEps_RSM( ValEps,wchFoldNameEps ) ;


	}


	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  wchFoldName);
	wcscat(wchAxesFileName, L"AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-10., 11.
	,-6,6,0.5) ;

}


// ���������� �������� ��������� � ������������ ���������
 // ��� ������ �������:
 // 1.���������� �������� ��������� ������ ��������� �� � ��  ��� ������������� ��  ValBet
 // � ����������� �� ���� �����  � �������� -10 +10 ����
 // INPUT:
 // ValBet - ��
 // ValEps - ��
 // wchFoldName1 - ����  � ����� � �������
  void TFar_2D::createSystemErrorGraphsFixedBetta_RSM(const double ValBet,wchar_t *wchFoldName1 )
{
	wchar_t wchFoldName[300] ={0};
	wcscpy(  wchFoldName,  wchFoldName1);
	wcscat(wchFoldName, L"\\");

	const int nBuffCols = 3;
	double step = 0.1;
	const int nBuffRows = 20./ step + 1;

	double  *parrBuff = new double [nBuffRows  * nBuffCols] ;
	memset(parrBuff, 0, nBuffRows * nBuffCols * sizeof(double));

	const int lenName = 30;
	wchar_t wcharrFileNames [150] ={0};
	wcscpy( wcharrFileNames, L"Eps");
	wcscpy( &wcharrFileNames[30], L"DeltaSystBet");
	wcscpy( &wcharrFileNames[60], L"DeltaSystEps");




	for (int i=0 ; i < nBuffRows; i++)
	{
		const double ValEps = (-10.0 + ((double)i) * step)* 0.001;
		TComp *pcmpTrueAmMeasures  = (TComp *)malloc(mNumAMCols * mNumAMRows * sizeof(TComp));
		TComp *pcmpNoisedAmMeasures  = (TComp *)malloc(mNumAMCols * mNumAMRows * sizeof(TComp));
		TSingleSign targSign(ValBet, ValEps, 800., 0.);

		 fncImitateMeasuresArray (targSign, pcmpTrueAmMeasures, pcmpNoisedAmMeasures) ;
		 double valEstBet = 0.,  valSigBet = 0., valEstEps = 0.,  valSigEps =0.;
		  calcEstRSM(pcmpTrueAmMeasures, &valEstBet, &valSigBet
		 , &valEstEps, &valSigEps) ;

		parrBuff[i * nBuffCols] = (-10. +((double)i) * step);
		parrBuff[i * nBuffCols + 1] = (valEstBet - ValBet) * 1000.;
		parrBuff[i * nBuffCols + 2] = (valEstEps - ValEps) * 1000.;


	}

	double scaley = 5.;
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
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-10., 11.
	,-6,6,0.5) ;


}


 // ��� ������ �������:
 // ���������� �������� ��������� ������ ��������� �� � ��  ��� ������������� ��  ValEps
 // � ����������� �� ��  � �������� -10 +10 ����
 // INPUT:
 // ValEps - ��
 // wchFoldName1 - ����  � ����� � �������


 void TFar_2D::createSystemErrorGraphsFixedEps_RSM( const double ValEps, wchar_t *wchFoldName1 )
{
	wchar_t wchFoldName[300] ={0};
	wcscpy(  wchFoldName,  wchFoldName1);
	wcscat(wchFoldName, L"\\");

	const int nBuffCols = 5;
	double step = 0.1;
	const int nBuffRows = 20./ step + 1;

	double  *parrBuff = new double [nBuffRows  * nBuffCols] ;
	memset(parrBuff, 0, nBuffRows * nBuffCols * sizeof(double));

	const int lenName = 30;
	wchar_t wcharrFileNames [150] ={0};
	wcscpy( wcharrFileNames, L"Bet");
	wcscpy( &wcharrFileNames[30], L"DeltaSystBet");
	wcscpy( &wcharrFileNames[60], L"DeltaSystEps");
	wcscpy( &wcharrFileNames[90], L"IterSigtBet");
	wcscpy( &wcharrFileNames[120], L"IterSigEps");





	for (int i=0 ; i < nBuffRows; i++)
	{
		const double ValBet = (-10.0 + ((double)i) * step)* 0.001;
		TComp *pcmpTrueAmMeasures  = (TComp *)malloc(mNumAMCols * mNumAMRows * sizeof(TComp));
		TComp *pcmpNoisedAmMeasures  = (TComp *)malloc(mNumAMCols * mNumAMRows * sizeof(TComp));
		TSingleSign targSign(ValBet, ValEps, 800., 0.);
	  //TSingleSign targSign(0.0017, 0., 800., 0.);

		 fncImitateMeasuresArray (targSign, pcmpTrueAmMeasures, pcmpNoisedAmMeasures) ;
		 double valEstBet = 0.,  valSigBet = 0., valEstEps = 0.,  valSigEps =0.;

		 double valDispUp = 0., valDispDown = 0.;
		 TComp cmpUp = calcUpperDiagr(pcmpTrueAmMeasures, &valDispUp)  ;
		 TComp cmpDown = calcDownDiagr(pcmpTrueAmMeasures, &valDispDown)  ;

		 double valD =   mdAMRow * ((double)mNumAMRows) /2.;
		 double valXiSq4Eps =0.;
		 calcIdealEstRSM(cmpUp, cmpDown,  valDispUp
		   , valDispDown,  valD, &valEstEps, &valSigEps, &valXiSq4Eps) ;


		  double valDispRight = 0., valDispLeft = 0.;
		 TComp cmpRight = calcRigthDiagr(pcmpTrueAmMeasures, &valDispRight)  ;
		 TComp cmpLeft = calcLeftDiagr(pcmpTrueAmMeasures, &valDispLeft)  ;

		 valD =   mdAMCol * ((double)mNumAMCols) /2.;
		 double valXiSq4Bet = 0.;
		 calcIdealEstRSM(cmpRight, cmpLeft,  valDispRight
		   , valDispLeft,  valD, &valEstBet, &valSigBet, &valXiSq4Bet) ;

	   //	 double valEstBet1 = 0.,  valSigBet1 = 0., valEstEps1 = 0.,  valSigEps1 =0.;
		// calcEstRSM(pcmpTrueAmMeasures, &valEstBet1, &valSigBet1
		// , &valEstEps1, &valSigEps1) ;
		double valIterEstBet, valIterSigBet, valIterEstEps, valIterSigEps ;
		double valIterXiSq4Bet, valIterXiSq4Eps;
		calcEstIterRSM( cmpUp,cmpDown,valDispUp, valDispDown
   ,cmpRight, cmpLeft, valDispRight, valDispLeft
  , &valIterEstBet, &valIterSigBet, &valIterXiSq4Bet
  , &valIterEstEps, &valIterSigEps, &valIterXiSq4Eps)  ;


		parrBuff[i * nBuffCols] = (-10. +((double)i) * step);
		parrBuff[i * nBuffCols + 1] = (valEstBet - ValBet) * 1000.;
		parrBuff[i * nBuffCols + 2] = (valEstEps - ValEps) * 1000.;
		parrBuff[i * nBuffCols + 3] = valIterSigBet * 1000.;
		parrBuff[i * nBuffCols + 4] = valIterSigEps * 1000.;


	}


	double scaley = 5.;
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
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-10., 11.
	,-6.,6.,0.5) ;

}




// ���������� ������� �� ������� 8 ��� ������������� ��� ��� ����������� ������ ���
// - ����� ��������� ������ ����� �� ��������
// �������� 4 ������������ ���������
// INPUT:
// pcmparrQuaterMeasures[4] - ������ ��������� ������������ ��������
// parrQuaterDisp[4] - ������ ��������� ������ ��������� ������������ ��������
// ValEstBet - ������ ��
// ValEstEps - ������ ��
// OUTPUT:
// ��8 �������
double TFar_2D::calcRealXi8Sq_RSM (TComp *pcmparrQuaterMeasures, double *parrQuaterDisp
	,const double ValEstBet,const double ValEstEps)
{
 TComp cmparrRoQuarters[4],cmparrQuaterMeasuresWave[4];
 memset(cmparrRoQuarters, 0, 4 * sizeof(TComp));
 memset(cmparrQuaterMeasuresWave, 0, 4 * sizeof(TComp));
 double  arrQuaterDisp[4] = {0.};
 calcFuncRoQuarters(  ValEstBet,  ValEstEps, cmparrRoQuarters) ;
 for (int i =0; i < 4; i++)
 {
   if(cmparrRoQuarters[i].modul()  < 0.0000001)
   {
	   continue;
   }
   cmparrQuaterMeasuresWave[i] =  pcmparrQuaterMeasures[i]/ cmparrRoQuarters[i];
   arrQuaterDisp[i] = parrQuaterDisp[i]/ cmparrRoQuarters[i].modul()/cmparrRoQuarters[i].modul();

 }

 TComp cmparrLamb[4];
 memset(cmparrLamb, 0, 4 * sizeof(TComp));
 double valSum = 0.;
 for (int i =0; i < 4; i++)
 {
  if(cmparrRoQuarters[i].modul()  < 0.0000001)
   {
	   continue;
   }
  valSum += 1./ arrQuaterDisp[i] ;
 }
 for (int i =0; i < 4; i++)
 {
 if(cmparrRoQuarters[i].modul()  < 0.0000001)
   {
	   continue;
   }
 cmparrLamb[i] = TComp(1./ arrQuaterDisp[i], 0.)/ TComp(valSum , 0.) ;
 }

 TComp cmp_aEst(0.,0.);
 MtrxMultMatrx(cmparrLamb,1, 4, cmparrQuaterMeasuresWave,1, &cmp_aEst);

 double valXi8Sq = 0.;
 for (int i = 0; i < 4; i++)
 {
  if(cmparrRoQuarters[i].modul()  < 0.0000001)
   {
	   continue;
   }
  TComp cmpDel = cmparrQuaterMeasuresWave[i] -  cmp_aEst;
  valXi8Sq += cmpDel.modul() * cmpDel.modul()/arrQuaterDisp[i];
 }

return 2. * valXi8Sq ;
}

//


void TFar_2D::calcFuncRoQuarters( const double ValEstBet, const double ValEstEps, TComp *pcmparrRoQuarters)
{
	memset(pcmparrRoQuarters, 0, 4 * sizeof(TComp));
	for (int i = 0; i < mNumAMRows /2   ; i++)
	for (int j = 0; j < mNumAMCols / 2; j++)
	{
		double temp = fnc_qml(ValEstBet,  ValEstEps, i, j);
	  pcmparrRoQuarters[0] +=   TComp(mpAm2D[ i * mNumAMCols + j].mOtklCoefUs , 0.)
	 * exp_(temp * 2. * M_PI / mLambda);


	  temp = fnc_qml(ValEstBet,  ValEstEps, i, mNumAMCols / 2 + j);
	  pcmparrRoQuarters[1] +=   TComp(mpAm2D[ i * mNumAMCols + mNumAMCols / 2 + j].mOtklCoefUs , 0.)
	 * exp_(temp * 2. * M_PI / mLambda);

	  temp = fnc_qml(ValEstBet,  ValEstEps, i + mNumAMRows /2,  j);
	  pcmparrRoQuarters[2] +=   TComp(mpAm2D[ (i + mNumAMRows /2) * mNumAMCols  + j].mOtklCoefUs , 0.)
	 * exp_(temp * 2. * M_PI / mLambda);


	  temp = fnc_qml(ValEstBet,  ValEstEps, i + mNumAMRows /2,  mNumAMCols / 2 + j);
	  pcmparrRoQuarters[3] +=   TComp(mpAm2D[ (i + mNumAMRows /2) * mNumAMCols  + mNumAMCols / 2 + j].mOtklCoefUs , 0.)
	 * exp_(temp * 2. * M_PI / mLambda);

	}

}



// ���������� �������� ������ ���� ��������� ��� �������������
// ����� ����� � ����� ����
// ��� ������ �������:
 // �� ������ ��������� ����������� ��������� ��
 // ���������� ������, �����, ���� � ������ ���������(��������� ���� ��������)
 // �� ���� ���������� ����������� �� � �� ������� ��� � ������������ ���
 // ����������� ��� � ������� �� �������
 // �������� �������
 // INPUT:
 // SingleSig- ������  �� ���� ��������
 // Nisp - ����� ���������
 // barrRealAM -������ ����������� ������� ����������  ��
 // wchFoldName1 - ����  � ����� � �������
 void TFar_2D::createMonteCarloRSM( const TSingleSign SingleSig,const int Nisp
   , const int QuantMeasXi2, bool *barrRealAM, wchar_t *wchFoldName1 )
 {
	TFar_2D far2DReal = *this;
	for ( int i = 0; i < mNumAMCols * mNumAMRows; i++)
	{
	  if (!barrRealAM[i])
	  {
		far2DReal.mpAm2D[i].mOtklCoefUs =0.;
		far2DReal.mpAm2D[i].mSigEmitNoise = 0.;
      }
	}
	wchar_t wchFoldName[300] ={0};
	wcscpy(  wchFoldName,  wchFoldName1);
	wcscat(wchFoldName, L"\\");

	const int nBuffCols = 21;

	const int nBuffRows = Nisp;

	double  *parrBuff = new double [nBuffRows  * nBuffCols] ;
	memset(parrBuff, 0, nBuffRows * nBuffCols * sizeof(double));

	const int lenName = 30;
	wchar_t *wcharrFileNames  = new wchar_t[lenName * nBuffCols] ;
	memset(wcharrFileNames, 0, sizeof (wchar_t )* lenName * nBuffCols);
	wcscpy( wcharrFileNames, L"n");
	wcscpy( &wcharrFileNames[lenName], L"TrueBet");
	wcscpy( &wcharrFileNames[2 * lenName], L"RsmEstBet");
	wcscpy( &wcharrFileNames[3 * lenName], L"IterRsmEstBet");
	wcscpy( &wcharrFileNames[4 * lenName], L"DelRsmBet");
	wcscpy( &wcharrFileNames[5 * lenName], L"DelIterRsmBet");
	wcscpy( &wcharrFileNames[6 * lenName], L"XiSq4RSMBet");
	wcscpy( &wcharrFileNames[7 * lenName], L"XiSq4IterRSMBet");
	wcscpy( &wcharrFileNames[8 * lenName], L"SigIterRSMBet");

	wcscpy( &wcharrFileNames[9 * lenName], L"TrueEps");
	wcscpy( &wcharrFileNames[10 * lenName], L"RsmEstEps");
	wcscpy( &wcharrFileNames[11 * lenName], L"IterRsmEstEps");
	wcscpy( &wcharrFileNames[12 * lenName], L"DelRsmEps");
	wcscpy( &wcharrFileNames[13 * lenName], L"DelIterRsmEps");
	wcscpy( &wcharrFileNames[14 * lenName], L"XiSq4RSMEps");
	wcscpy( &wcharrFileNames[15 * lenName], L"XiSq4IterRSMEps");
	wcscpy( &wcharrFileNames[16 * lenName], L"SigIterRSMEps");
	wcscpy( &wcharrFileNames[17 * lenName], L"Xi4Sq_0_98");
	wcscpy( &wcharrFileNames[18 * lenName], L"Xi8Sq");
	wcscpy( &wcharrFileNames[19 * lenName], L"Xi8Sq_0_98");
	wcscpy( &wcharrFileNames[20 * lenName], L"Xi8SqIdealRSM");

	for (int i=0 ; i < nBuffRows; i++) 	{

		TComp *pcmpTrueAmMeasures  = (TComp *)malloc(mNumAMCols * mNumAMRows * sizeof(TComp));
		TComp *pcmpNoisedAmMeasures  = (TComp *)malloc(mNumAMCols * mNumAMRows * sizeof(TComp));

		 far2DReal.fncImitateMeasuresArray ( SingleSig, pcmpTrueAmMeasures, pcmpNoisedAmMeasures) ;
		 ///

	   TComp cmparrQuaterMeasures[4];
	   double arrQuaterDisp[4] = {0.};
	   far2DReal.calcQuarterDiagrams(pcmpNoisedAmMeasures, cmparrQuaterMeasures, arrQuaterDisp);
	   ///
		 double valEstBet = 0.,  valSigBet = 0., valEstEps = 0.,  valSigEps =0.;

		 double valDispUp = 0., valDispDown = 0.;
		 TComp cmpUp = far2DReal.calcUpperDiagr(pcmpNoisedAmMeasures, &valDispUp)  ;
		 TComp cmpDown = far2DReal.calcDownDiagr(pcmpNoisedAmMeasures, &valDispDown)  ;

		 double valD =   mdAMRow * ((double)mNumAMRows) /2.;
		 double valXiSq4Eps =0.;
		 calcIdealEstRSM(cmpUp, cmpDown,  valDispUp
		   , valDispDown,  valD, &valEstEps, &valSigEps, &valXiSq4Eps) ;


		  double valDispRight = 0., valDispLeft = 0.;
		 TComp cmpRight = far2DReal.calcRigthDiagr(pcmpNoisedAmMeasures, &valDispRight)  ;
		 TComp cmpLeft = far2DReal.calcLeftDiagr(pcmpNoisedAmMeasures, &valDispLeft)  ;

		 valD =   mdAMCol * ((double)mNumAMCols) /2.;
		 double valXiSq4Bet = 0.;
		 calcIdealEstRSM(cmpRight, cmpLeft,  valDispRight
		   , valDispLeft,  valD, &valEstBet, &valSigBet, &valXiSq4Bet) ;
		// ���������� ��4 ������� ��� ��������� ������ ���
		 valXiSq4Eps = calcRealXi4Sq_Eps_RSM (cmpUp	,cmpDown,valEstBet,valEstEps);
		 valXiSq4Bet = calcRealXi4Sq_Bet_RSM (cmpRight	,cmpLeft,valEstBet,valEstEps);
		 // ���������� ��8 ������� ��� ��������� ������ ���
		double valXi8SqIdeal =  calcRealXi8Sq_RSM (cmparrQuaterMeasures,arrQuaterDisp,valEstBet,valEstEps);
		///

		double valIterEstBet, valIterSigBet, valIterEstEps, valIterSigEps ;
		double valIterXiSq4Bet, valIterXiSq4Eps;
		calcEstIterRSM( cmpUp,cmpDown,valDispUp, valDispDown
   ,cmpRight, cmpLeft, valDispRight, valDispLeft
  , &valIterEstBet, &valIterSigBet, &valIterXiSq4Bet
  , &valIterEstEps, &valIterSigEps, &valIterXiSq4Eps)  ;


	   // ���������� �� ������� 8 ��� ������ ���������� ������������ �������
	   // ����� ����������� ���������
	   double valXi8Sq = calcRealXi8Sq_RSM (cmparrQuaterMeasures, arrQuaterDisp, valIterEstBet,valIterEstEps) ;

		parrBuff[i * nBuffCols] = ((double)i) + 1.;
		parrBuff[i * nBuffCols + 1] = SingleSig.mBet * 1000.;
		parrBuff[i * nBuffCols + 2] = valEstBet * 1000.;
		parrBuff[i * nBuffCols + 3] = valIterEstBet * 1000.;
		parrBuff[i * nBuffCols + 4] = (valEstBet - SingleSig.mBet) * 1000.* 10.;
		parrBuff[i * nBuffCols + 5] = (valIterEstBet - SingleSig.mBet) * 1000.* 10.;
		parrBuff[i * nBuffCols + 6] =  valXiSq4Bet;
		parrBuff[i * nBuffCols + 7] =  valIterXiSq4Bet;
		parrBuff[i * nBuffCols + 8] =  valIterSigBet*1000.* 10.;

		parrBuff[i * nBuffCols + 9] = SingleSig.mEps * 1000.;
		parrBuff[i * nBuffCols + 10] = valEstEps * 1000.;
		parrBuff[i * nBuffCols + 11] = valIterEstEps * 1000.;
		parrBuff[i * nBuffCols + 12] = (valEstEps - SingleSig.mEps) * 1000.* 10.;
		parrBuff[i * nBuffCols + 13] = (valIterEstEps - SingleSig.mEps) * 1000.* 10.;
		parrBuff[i * nBuffCols + 14] =  valXiSq4Eps;
		parrBuff[i * nBuffCols + 15] =  valIterXiSq4Eps;
		parrBuff[i * nBuffCols + 16] =  valIterSigEps * 1000.* 10.;

		parrBuff[i * nBuffCols + 17] =  11.6;
		parrBuff[i * nBuffCols + 18] =  valXi8Sq ;
		parrBuff[i * nBuffCols + 19] = 18.2 ;
		parrBuff[i * nBuffCols + 20] =  valXi8SqIdeal;
	}


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
		,1. // ������� �� ��� Y
		) ;
	}

	// ���������� �������� ������� ������������� ��� ������������ � ������������ ���
	int nBuffCols1 = 4;
	double  *parrBuff1 = new double [(nBuffRows -QuantMeasXi2)  * nBuffCols1] ;
	memset(parrBuff1, 0, (nBuffRows -QuantMeasXi2) * nBuffCols1 * sizeof(double));
	double valTreshold_0_98 = 0.;
	switch(QuantMeasXi2)
	{   case 1:
		valTreshold_0_98= 18.2;
		break;
		case 2:
		valTreshold_0_98 = 29.6;
		break;

		case 3:
		 valTreshold_0_98 = 40.3;
		break;
		default:
		valTreshold_0_98 = 2.1;

    }


	for (int i =0; i < (nBuffRows -QuantMeasXi2); i++)
	{
	   parrBuff1[i *  nBuffCols1 ] = QuantMeasXi2 + i;
	   for (int j = 0; j < QuantMeasXi2; j++)
	   {
		 parrBuff1[i *  nBuffCols1 + 1] += parrBuff[(i + j) * nBuffCols + 18] ;
		 parrBuff1[i *  nBuffCols1 + 2] += parrBuff[(i + j) * nBuffCols + 20] ;
	   }
	   parrBuff1[i *  nBuffCols1 + 3] = valTreshold_0_98;
	   if (QuantMeasXi2 > 3)
	   {
		 parrBuff1[i *  nBuffCols1 + 1] = sqrt(2. * parrBuff1[i *  nBuffCols1 + 1] ) - sqrt(2. * (8. *  ((double)QuantMeasXi2))-1.) ;
		 parrBuff1[i *  nBuffCols1 + 2] = sqrt(2. * parrBuff1[i *  nBuffCols1 + 2] ) - sqrt(2. * (8. *  ((double)QuantMeasXi2)) -1.) ;
	   }
   }

   delete []parrBuff;
	memset(wcharrFileNames, 0, sizeof (wchar_t )* lenName * nBuffCols);
	wcscpy( wcharrFileNames, L"n");
	wcscpy( &wcharrFileNames[lenName], L"XiSqSumIter");
	wcscpy( &wcharrFileNames[2 * lenName], L"XiSqSumIdeal");
	wcscpy( &wcharrFileNames[3 * lenName], L"Treshold_0_98_Sum");

	for (int i=1; i < nBuffCols1; i++)
	{
		TYrWriteShapeFile::WriteOneReport(                 wchFoldName  // ���� � �����
		, parrBuff1 // ������ � ����������� - ������� nBuffRows x nBuffCols
		, nBuffCols1 // - �-�� ���������� � ������� ��������� ���������� � ������
		, nBuffRows -QuantMeasXi2 //  - �-�� �����
		,wcharrFileNames //������� � �������� ���������� - ������� nBuffCols x lenName
		,lenName // ������������ ����� ����� ����������
		,0  // ����� ���������� �� ��� X
		,i  // ����� ���������� �� ��� Y
		,1. //  ������� �� ��� X
		,1. // ������� �� ��� Y
		) ;
	}

	delete []parrBuff1;

	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  wchFoldName);
	wcscat(wchAxesFileName, L"AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,0., ((double)Nisp)*1.1,-20.,20.,1.5) ;


	wcscpy(  wchAxesFileName,  wchFoldName);
	wcscat(wchAxesFileName, L"FishNet.shp");
	TURPolyLine  plnFishNet = TURPolyLine::createFishNet(0., ((double)Nisp)*1.1
   ,-26. ,18., 10., 2. );
   plnFishNet.WriteSetSHPFiles(wchAxesFileName, &plnFishNet, 1 );
	return;
 }

 /////////////////////////////////////////////////////////////////////////////////////////////////

double  TFar_2D::calcSigEtalonRSM_Eps(TEtalonSign EtalonSign, double valEps)
{

	const double VAl_NWaveEtalon = calcNWaveEtalon (EtalonSign);

	const double VAlMuCourse = mLambda / 2./ M_PI /( ((double)(mNumAMRows * mpAm2D[0].mNumEmitRows)) * mpAm2D[0].mdRow / 2.);
	return VAlMuCourse   *VAl_NWaveEtalon ;
}
double  TFar_2D::calcSigEtalonRSM_Bet(TEtalonSign EtalonSign, double valBet)
{

	const double VAl_NWaveEtalon = calcNWaveEtalon (EtalonSign);
	const double VAlMuCourse = mLambda / 2./ M_PI /( ((double)(mNumAMCols * mpAm2D[0].mNumEmitCols)) * mpAm2D[0].mdCol / 2.);
	return VAlMuCourse   *VAl_NWaveEtalon ;
}

int TFar_2D::calcQuantWorkingAM()
{
   int irez = 0;
   for (int i =0; i < (mNumAMCols * mNumAMRows); i++)
   {
	irez += mpbarrAM[i];
   }
  return irez;
}

// ���������� ���������� ��������� ���/������
double  TFar_2D::calcNWaveEtalon (TEtalonSign EtalonSign)
{
		const int ItotalQuantEmit = calcQuantWorkingAM()* mpAm2D[0].mNumEmitCols * mpAm2D[0].mNumEmitRows;
	const double VAl_NWaveEtalon = sqrt( EtalonSign.mNoiseSKZ_5P10 * EtalonSign.mNoiseSKZ_5P10
	/ EtalonSign.mEtalonAmp/ EtalonSign.mEtalonAmp
	* 56.* 64/((double)ItotalQuantEmit));
	 //  + ((double)ItotalQuantEmit) /( 56. * 64)* EtalonSign.mEtalonSigAmplFact_5P10* EtalonSign.mEtalonSigAmplFact_5P10);
   return VAl_NWaveEtalon;
}
/*
// ���� ������ ���������� ����������������� �
 double TFar_2D::calc_TwoTargsResolutionDist(const double VAlTargH, const double VAlTargEPR
  ,const double PowerPrd, const double KYPrd, const TEtalonSign ETalonSign , const double VAlHAntenna, double *pvalSigTarg )
{
   double valHorisDist0 = 100000.;
   double valAntNormalAng = 0.;
   TFar Far(*this, true);

   double valDistCur =  0.;
   double valStep = 50.;
   const int INtNC = valHorisDist0/ valStep -1;

   for (int i = 0; i < INtNC; i++)
   {
		double valHorisDistCur = valHorisDist0 - ((double)i) * valStep;
		valDistCur =  sqrt(VAlTargH * VAlTargH  + valHorisDistCur * valHorisDistCur);

		double alfUMTrg = asin( (VAlTargH- VAlHAntenna)/ valDistCur);
		double alfUMAntp = -atan( (VAlTargH+ VAlHAntenna) / valHorisDistCur);
		alfUMAntp = alfUMAntp - alfUMTrg;
		alfUMTrg = 0.;
		// double temp = calcKYPriem();
		//	double valBCur =  calcCurrentSignalAmpl ( valDistCur ,  VAlTargEPR
		//  , PowerPrd,  KYPrd, ETalonSign) ;

		double valSigTarg1 = 0., valSigQ1 = 0.;
		double valAntpPhaze = getRand01( ) * 2. * M_PI;
		if( calc_SKZ_LAT(valDistCur,VAlTargH ,VAlHAntenna, valAntNormalAng,VAlTargEPR
		,ETalonSign  , PowerPrd, KYPrd, valAntpPhaze, &valSigTarg1 , &valSigQ1  ))
		{
		*pvalSigTarg  =  valSigTarg1;
		return valDistCur ;
		}
   }
   return 0.;

}

 */


// ���� ������ ���������� ����������������� �
 double TFar_2D::calc_TwoTargsResolutionDist(const double VAlCoefAntp,const double VAlTargH, const double VAlTargEPR
  ,const double PowerPrd, const double KYPrd, const TEtalonSign ETalonSign , const double VAlHAntenna, double *pvalSigTarg )
{
   double valHorisDist0 = 100000.;
   double valAntNormalAng = 0.;
   TFar Far(*this, true);

   double valDistCur =  0.;
   double valStep = 50.;
   const int INtNC = valHorisDist0/ valStep -1;

   for (int i = 0; i < INtNC; i++)
   {
		double valHorisDistCur = valHorisDist0 - ((double)i) * valStep;
		valDistCur =  sqrt(VAlTargH * VAlTargH  + valHorisDistCur * valHorisDistCur);

		double alfUMTrg = asin( (VAlTargH- VAlHAntenna)/ valDistCur);
		double alfUMAntp = -atan( (VAlTargH+ VAlHAntenna) / valHorisDistCur);
		alfUMAntp = alfUMAntp - alfUMTrg;
		alfUMTrg = 0.;
		// double temp = calcKYPriem();
		//	double valBCur =  calcCurrentSignalAmpl ( valDistCur ,  VAlTargEPR
		//  , PowerPrd,  KYPrd, ETalonSign) ;

		double valSigTarg1 = 0., valSigQ1 = 0.;
		double valAntpPhaze = getRand01( ) * 2. * M_PI;
		if( calc_SKZ_LAT(VAlCoefAntp,valDistCur,VAlTargH ,VAlHAntenna, valAntNormalAng,VAlTargEPR
		,ETalonSign  , PowerPrd, KYPrd, valAntpPhaze, &valSigTarg1 , &valSigQ1  ))
		{
		*pvalSigTarg  =  valSigTarg1;
		return valDistCur ;
		}
   }
   return 0.;

}



// ���������� ��������� ������� ��� ��� ������� ������ ��������
// � ���������� ������ ��������� ������� ����
// ���� ����� ������� � ������� ����� 0 ���
// OUTPUT:
// *pvalSigE, *pvalSigQ - ��� ������ �� � �� �� ������ �������
// *bLat - ������� ��� �� ������ ������� (= true, ���� �� ���)
//
//
double TFar_2D::calc_TwoTargsZahvatDist(const double VAlCoefAntp,const double VAlTargH, const double VAlTargEPR
  ,const double PowerPrd, const double KYPrd, const TEtalonSign ETalonSign , const double VAlHAntenna
  , bool *bLat  , double *pvalSigE ,double *pvalSigQ  )
{
   double valHorisDist0 = 22150.;//100000.;
   double valAntNormalAng = 0.; // ���� ����� ������� � �������
   TFar Far0(*this, true);
   TFar Far(Far0, 4) ;

	// ������ ��������� �����:
   const double VAlSumDiagrWidthVert = findDiagrWidth(Far.mFaceta.m_d  ,Far.m_D, Far.mFaceta.m_n
  , Far.m_N ,mLambda);

   TFar Far00(*this, false);
   TFar Far01(Far00, 2) ;
   // ������ ��������� ����� � �����. ���������:
   const double VAlSumDiagrWidthHor = findDiagrWidth(Far01.mFaceta.m_d  ,Far01.m_D, Far01.mFaceta.m_n
  , Far01.m_N ,mLambda);

   double valDistCur =  0.;
   double valStep = 50.;
   const int INtNC = valHorisDist0/ valStep -1;

   for (int i = 0; i < INtNC; i++)
   {

	if (i == 401)
	{
      int iii=0;
	}
	 double valHorisDistCur = valHorisDist0 - ((double)i) * valStep;
	 valDistCur =  sqrt(VAlTargH * VAlTargH  + valHorisDistCur * valHorisDistCur);

		// ������ �������
	 const double VAlBCur = calcCurrentSignalAmpl ( valDistCur,  VAlTargEPR
   , PowerPrd,  KYPrd, ETalonSign) ;

	//	const double VAl_NWaveEtalon   =     calcNWaveEtalon ( ETalonSign);

     double temp = sqrt(calcSumDisp());
	 double valNWaveCur = sqrt(calcSumDisp())/ VAlBCur;

	 double valAppertHor = calcAppertHor() ;
	// double valAppertVert = calcAppertVert() ;

	 *pvalSigQ = sqrt(TFar::calcTheorDisp_RSM(valNWaveCur, valAppertHor,mLambda, 0.005));
	 if ((*pvalSigQ ) * 3. >  VAlSumDiagrWidthHor)
	 {
		continue;
	 }


	 double valAntpPhaze = getRand01( ) * 2. * M_PI;  // ����  ������� ��������
	//double valAntpPhaze = M_PI/2.;
	 calc_SKZ_LAT( VAlCoefAntp,valDistCur, VAlTargH, VAlHAntenna, valAntNormalAng, VAlTargEPR
	,ETalonSign , PowerPrd, KYPrd, valAntpPhaze, pvalSigE ,pvalSigQ  );

	if (3. * (*pvalSigE) > VAlSumDiagrWidthVert)
		{
		   continue;
		}
		else
		{
		  return valDistCur ;
		}


   }
   return -1.;

}
 /*
// ���������� ��������� ������� ��� ��� ������� ������ ��������
// � ���������� ������ ��������� ������� ����
// ���� ����� ������� � ������� ����� 0 ���
// OUTPUT:
// *pvalSigE, *pvalSigQ - ��� ������ �� � �� �� ������ �������
// *bLat - ������� ��� �� ������ ������� (= true, ���� �� ���)
//
//
double TFar_2D::calc_TwoTargsZahvatDist(const double VAlTargH, const double VAlTargEPR
  ,const double PowerPrd, const double KYPrd, const TEtalonSign ETalonSign , const double VAlHAntenna
  , bool *bLat  , double *pvalSigE ,double *pvalSigQ  )
{
   double valHorisDist0 = 22150.;//100000.;
   double valAntNormalAng = 0.; // ���� ����� ������� � �������
   TFar Far0(*this, true);
   TFar Far(Far0, 4) ;

	// ������ ��������� �����:
   const double VAlSumDiagrWidthVert = findDiagrWidth(Far.mFaceta.m_d  ,Far.m_D, Far.mFaceta.m_n
  , Far.m_N ,mLambda);

   TFar Far00(*this, false);
   TFar Far01(Far00, 2) ;
   // ������ ��������� ����� � �����. ���������:
   const double VAlSumDiagrWidthHor = findDiagrWidth(Far01.mFaceta.m_d  ,Far01.m_D, Far01.mFaceta.m_n
  , Far01.m_N ,mLambda);

   double valDistCur =  0.;
   double valStep = 50.;
   const int INtNC = valHorisDist0/ valStep -1;

   for (int i = 0; i < INtNC; i++)
   {
	 double valHorisDistCur = valHorisDist0 - ((double)i) * valStep;
	 valDistCur =  sqrt(VAlTargH * VAlTargH  + valHorisDistCur * valHorisDistCur);

		// ������ �������
	 const double VAlBCur = calcCurrentSignalAmpl ( valDistCur,  VAlTargEPR
   , PowerPrd,  KYPrd, ETalonSign) ;

	//	const double VAl_NWaveEtalon   =     calcNWaveEtalon ( ETalonSign);


	 double valNWaveCur = sqrt(calcSumDisp())/ VAlBCur;

	 double valAppertHor = calcAppertHor() ;
	// double valAppertVert = calcAppertVert() ;

	 *pvalSigQ = sqrt(TFar::calcTheorDisp_RSM(valNWaveCur, valAppertHor,mLambda, 0.005));
	 if ((*pvalSigQ ) * 3. >  VAlSumDiagrWidthHor)
	 {
		continue;
	 }


	 double valAntpPhaze = getRand01( ) * 2. * M_PI;  // ����  ������� ��������
	 calc_SKZ_LAT( valDistCur, VAlTargH, VAlHAntenna, valAntNormalAng, VAlTargEPR
	,ETalonSign , PowerPrd, KYPrd, valAntpPhaze, pvalSigE ,pvalSigQ  );

	if (3. * (*pvalSigE) > VAlSumDiagrWidthVert)
		{
		   continue;
		}
		else
		{
		  return valDistCur ;
		}


   }
   return -1.;

}*/

 /*

// ���������� ��� ������ ��������� �� �  �� ����
// ��������������, ��� � ���� ������� ������������ ������� � �������� �� ������ ����������� � �������� ����
// � ��������� �����, ������� ������������� �������� ��������� �����
// INPUT
// VAlTargD
// VAlTargD
// VAlTargD
// VAlAntNormalAng
// VAlTargEPR
// EtalonSign
// PowerPrd
// KYPrd
// valAntpPhaze
// OUTPUT:
// *pvalSigE, *pvalSigQ
// ���������� true, e��� ������ ��������� �� ��������� ������ ��������� ��������
// false, ���� �� ���������-���������� ������
bool TFar_2D::calc_SKZ_LAT(const double VAlTargD,const double VAlTargH
  ,const double VAlHAntenna, const double VAlAntNormalAng, const double VAlTargEPR
  , const TEtalonSign EtalonSign ,const double PowerPrd,const double KYPrd, double valAntpPhaze
	, double *pvalSigE ,double *pvalSigQ  )
{

	bool bRazreshenie = false;

	double alfUMTrg = asin( (VAlTargH - VAlHAntenna)/ VAlTargD);
	double alfUMAntp = -atan( (VAlTargH + VAlHAntenna)
	/ sqrt(VAlTargD * VAlTargD  - VAlTargH * VAlTargH));
	alfUMAntp = alfUMAntp - alfUMTrg + VAlAntNormalAng;
	alfUMTrg = VAlAntNormalAng;
	TFar Far0(*this, true);
	//
	 if (Far0.m_N % 4 != 0)
  {
	return false;
	}
	TFar Far(4) ;
	mLambda =  Far0.mLambda;
 int qRow = Far0.m_N / 4;  // ������� ��������  ����� � ������


 for (int i = 0; i < 4; i++)
 {

   for (int j =0; j < qRow; j++)
   {

	Far.mparrDisp[i] += Far0.mparrDisp[ i * qRow + j];
	Far.mparrAmplFactDisp [i] += Far0.mparrAmplFactDisp[i * qRow + j];
   }
	 Far.mparrAmplFactDisp [i] = Far.mparrAmplFactDisp [i]/((double)qRow);  // ????????????
 }
 ///

 Far.m_D = qRow * Far0.m_D;
 Far.mFaceta = TFaceta (Far0.mFaceta.m_n * qRow  , Far0.mFaceta.m_d, Far0.mFaceta.mLambda, Far.mparrDisp[0], Far.mparrAmplFactDisp [0]);
	///


	// ������ �������
	const double VAlBCur = calcCurrentSignalAmpl ( VAlTargD,  VAlTargEPR
   , PowerPrd,  KYPrd, EtalonSign) ;
	///




	 double valNWaveCur = sqrt(calcSumDisp())/ VAlBCur;

	 double valAppertHor = calcAppertHor() ;

	 *pvalSigQ = sqrt(TFar::calcTheorDisp_RSM(valNWaveCur, valAppertHor,mLambda, 0.005)) ;

	bRazreshenie = true;
	TComp  cmpKTarg(VAlBCur,0.);

	TComp cmpKAntp  = TComp(VAlBCur* cos(valAntpPhaze) ,VAlBCur  * sin(valAntpPhaze));


	*pvalSigE = Far.calc_Guaranted_SKZ_RSM_For_TwoTargs (  alfUMTrg, alfUMAntp, cmpKTarg, cmpKAntp) ;
	 //	*pvalSigE = Far.calc_Mean_SKZ_RSM_For_TwoTargs (  alfUMTrg, alfUMAntp, cmpKTarg, cmpKAntp) ;
		// *pvalSigE = 100000.;/// ������ !"!!!!!!
	TComp cmpTemp(0.25, 0.);
	cmpKTarg = cmpKTarg * cmpTemp;
	cmpKAntp = cmpKAntp * cmpTemp;
	double valSigETemp = 0.;
	Far.calcTheoretical_LATarg_SKZ( &bRazreshenie, &valSigETemp
	,alfUMTrg, alfUMAntp, cmpKTarg, cmpKAntp);
 //  valSigETemp  = 10000000.;// ������ !!!!!!!
	 if (valSigETemp < (*pvalSigE ))
		 {
		*pvalSigE = valSigETemp;
		return true;
	   }
	return false;

}
 */
// �������������!!!
// ���������� ��� ������ ��������� �� �  �� ����
// ��������������, ��� � ���� ������� ������������ ������� � �������� �� ������ �� ����������� � �������� ����
// � ��������� �����, ������� ������������� �������� ��������� �����
// INPUT
// VAlCoefAntp - ����� ��������
// VAlTargD
// VAlAntNormalAng
// VAlTargEPR
// EtalonSign
// PowerPrd
// KYPrd
// valAntpPhaze
// OUTPUT:
// *pvalSigE, *pvalSigQ
// ���������� true, e��� ������ ��������� �� ��������� ������ ��������� ��������
// false, ���� �� ���������-���������� ������
bool TFar_2D::calc_SKZ_LAT(const double VAlCoefAntp,const double VAlTargD,const double VAlTargH
  ,const double VAlHAntenna, const double VAlAntNormalAng, const double VAlTargEPR
  , const TEtalonSign EtalonSign ,const double PowerPrd,const double KYPrd, double valAntpPhaze
	, double *pvalSigE ,double *pvalSigQ  )
{

	bool bRazreshenie = false;

	double alfUMTrg = asin( (VAlTargH - VAlHAntenna)/ VAlTargD);
	double alfUMAntp = -atan( (VAlTargH + VAlHAntenna)
	/ sqrt(VAlTargD * VAlTargD  - VAlTargH * VAlTargH));
	alfUMAntp = alfUMAntp - alfUMTrg + VAlAntNormalAng;
	alfUMTrg = VAlAntNormalAng;
	TFar Far0(*this, true);
	//
	 if (Far0.m_N % 4 != 0)
  {
	return false;
	}
	TFar Far(4) ;
	mLambda =  Far0.mLambda;
 int qRow = Far0.m_N / 4;  // ������� ��������  ����� � ������


 for (int i = 0; i < 4; i++)
 {

   for (int j =0; j < qRow; j++)
   {

	Far.mparrDisp[i] += Far0.mparrDisp[ i * qRow + j];
	Far.mparrAmplFactDisp [i] += Far0.mparrAmplFactDisp[i * qRow + j];
   }
	 Far.mparrAmplFactDisp [i] = Far.mparrAmplFactDisp [i]/((double)qRow);  // ????????????
 }
 ///

 Far.m_D = qRow * Far0.m_D;
 Far.mFaceta = TFaceta (Far0.mFaceta.m_n * qRow  , Far0.mFaceta.m_d, Far0.mFaceta.mLambda, Far.mparrDisp[0], Far.mparrAmplFactDisp [0]);
	///


	// ������ �������
	const double VAlBCur = calcCurrentSignalAmpl ( VAlTargD,  VAlTargEPR
   , PowerPrd,  KYPrd, EtalonSign) ;
	///




	 double valNWaveCur = sqrt(calcSumDisp())/ VAlBCur;

	 double valAppertHor = calcAppertHor() ;
	 double valAppertVert = calcAppertVert() ;

	 *pvalSigQ = sqrt(TFar::calcTheorDisp_RSM(valNWaveCur, valAppertHor,mLambda, 0.004)) ;

	bRazreshenie = true;
	TComp  cmpKTarg(VAlBCur,0.);

	TComp cmpKAntp  = TComp(VAlCoefAntp *VAlBCur* cos(valAntpPhaze) ,VAlCoefAntp * VAlBCur  * sin(valAntpPhaze));

	double valSigETemp = 0.;

	 if (cmpKAntp.modul() < 0.05 * cmpKTarg.modul())
	 {
	   *pvalSigE = sqrt(TFar::calcTheorDisp_RSM(valNWaveCur, valAppertVert,mLambda, 0.004)) ;
	 }
	 else
	 {
		*pvalSigE = Far.calc_Guaranted_SKZ_RSM_For_TwoTargs (  alfUMTrg, alfUMAntp, cmpKTarg, cmpKAntp) ;
		TComp cmpTemp(0.25, 0.);
		cmpKTarg = cmpKTarg * cmpTemp;
		cmpKAntp = cmpKAntp * cmpTemp;

		Far.calcTheoretical_LATarg_SKZ( &bRazreshenie, &valSigETemp
		,alfUMTrg, alfUMAntp, cmpKTarg, cmpKAntp);


		// ����� ������ !!!
		*pvalSigE = valSigETemp;
        return true;

		if (valSigETemp < (*pvalSigE ))
		{
		*pvalSigE = valSigETemp;
		return true;
		}
		return false;
	 }

	 return true;



}



// ���������� ��� ������ ��������� �� �  �� ����
// ��������������, ��� � ���� ������� ������������ ������� � �������� �� ������ ����������� � �������� ����
// � ��������� �����, ������� ������������� �������� ��������� �����
// INPUT
// VAlTargD
// VAlTargD
// VAlTargD
// VAlAntNormalAng
// VAlTargEPR
// EtalonSign
// PowerPrd
// KYPrd
// valAntpPhaze
// OUTPUT:
// *pvalSigE, *pvalSigQ
// ���������� true, e��� ������ ��������� �� ��������� ������ ��������� ��������
// false, ���� �� ���������-���������� ������
bool TFar_2D::calc_SKZ_SingleTarg(const double VAlTargD
 , const double VAlAntNormalAng, const double VAlTargEPR
  , const TEtalonSign EtalonSign ,const double PowerPrd,const double KYPrd
	, double *pvalSigE ,double *pvalSigQ  )
{
	// ������ �������
	const double VAlBCur = calcCurrentSignalAmpl ( VAlTargD,  VAlTargEPR
   , PowerPrd,  KYPrd, EtalonSign) ;
   ///
	 double valNWaveCur = sqrt(calcSumDisp())/ VAlBCur;

	 double valAppertHor = calcAppertHor() ;
	 double valAppertVert = calcAppertVert() ;
	 *pvalSigQ = sqrt(TFar::calcTheorDisp_RSM(valNWaveCur, valAppertHor,mLambda, 0.005)) ;
	 *pvalSigE = sqrt(TFar::calcTheorDisp_RSM(valNWaveCur, valAppertVert ,mLambda, 0.005)) ;
	return true;
 }







double  TFar_2D::calcCurrentSignalAmpl (const double VAlTargD, const double VAlTargEPR
   ,const double PowerPrd, const double KYPrd, const TEtalonSign ETalonSign)
{
return sqrt(VAlTargEPR / ETalonSign.mEtalonEPR * PowerPrd/ ETalonSign.mEtalonPowerPrd *  calcKYPriem()/ ETalonSign.mEtalonKYPriem
		  * KYPrd/ ETalonSign.mEtalonKYPrd)
		* ETalonSign.mEtalonDist *  ETalonSign.mEtalonDist
		/ (VAlTargD * VAlTargD) *  ETalonSign.mEtalonAmp ;
}
 // ���������� ��������� ��������� ����������� ����
double TFar_2D::calcSumDisp()
{
   double valRez = 0.;
   for (int i =0; i < mNumAMCols * mNumAMRows; i++)
   {
	  if(mpbarrAM[i])
	  {
		valRez += mpAm2D[i].calcSumDisp();
      }
   }
  return valRez;
}

double TFar_2D::calcAppertHor()
{
	return ((double)mNumAMCols) * mdAMCol;
}

double TFar_2D::calcAppertVert()
{
	return ((double)mNumAMRows) * mdAMRow;
}



double TFar_2D::calcKYPriem()
{
 double valS = calcQuantWorkingAM()* mpAm2D[0].calcSquare()/ mLambda/mLambda;
 return 4. * M_PI * valS * 0.35;

}

void TFar_2D::fncEstimateCourseAngs_For_TwoTargs (TComp *pcmparrAmMeasures, const double VAlEstEpsTarg, const double VAlEstEpsAntp
	  , TComp CMpKTarg , TComp CMpKAntp , double *arrMtrxCorrEps,
		double *pvalEstBetTarg, double *pvalEstBetAntp , double *arrMtrxCorrBet)
{
	double valMuTarg =  sin(VAlEstEpsTarg ) * calcAppertVert() * 3./ 8. / mLambda * 2. * M_PI;
	TComp cmpZTarg = exp_(valMuTarg);

	double valMuAntp =  sin(VAlEstEpsAntp ) * calcAppertVert() * 3./ 8. / mLambda * 2. * M_PI;
	TComp cmpZAntp = exp_(valMuAntp);


	double valFTarg = fncDiagrPar(mpAm2D[0].mdRow , mdAMRow, mpAm2D[0].mNumEmitRows, mNumAMRows ,mLambda , 0.,VAlEstEpsTarg);
	double valFAntp = fncDiagrPar(mpAm2D[0].mdRow , mdAMRow, mpAm2D[0].mNumEmitRows, mNumAMRows ,mLambda , 0.,VAlEstEpsAntp);
	TComp cmpDiagrTarg (valFTarg,0.);
	TComp cmpDiagrAntp (valFAntp,0.);

	TComp cmp_a =  (CMpKTarg * cmpZTarg) * ( cmpDiagrTarg * TComp(0.,-1));
	TComp cmp_b =  (CMpKAntp * cmpZAntp) * (cmpDiagrAntp * TComp(0.,-1));

	double valLeftDisp = 0., valRightDisp = 0.;
	TComp cmpSRight = calcRigthDiagr (pcmparrAmMeasures, &valLeftDisp ) ;
	TComp cmpSLeft  = calcLeftDiagr(pcmparrAmMeasures, &valRightDisp ) ;
	TComp cmpDif = cmpSRight - cmpSLeft;

	double arr[4] ={0.}, arr_b [2] = {0.}, arrx[2] ={0.};
	arr[0] = cmp_a.m_Re;
	arr[1] = cmp_b.m_Re;
	arr[2] = cmp_a.m_Im;
	arr[3] = cmp_b.m_Im;

	arr_b[0] = cmpDif.m_Re;
	arr_b[1] = cmpDif.m_Im;

	SolvLinEq2(arr, arr_b,arrx);
	///
	double valTargGamma = atan(arrx[0]);
	double valAntpGamma = atan(arrx[1]);

	*pvalEstBetTarg = 2. * valTargGamma *mLambda/ M_PI/ calcAppertHor();
	*pvalEstBetAntp = 2. * valAntpGamma *mLambda/ M_PI/ calcAppertHor();
	return;
}


#pragma package(smart_init)


