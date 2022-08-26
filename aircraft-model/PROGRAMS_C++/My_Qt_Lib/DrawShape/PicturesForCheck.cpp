//---------------------------------------------------------------------------


#pragma hdrstop

#include "PicturesForCheck.h"


// ������� ��� ��������� ��������� ��������� � �� ������������� ��������� 4 �������
void Pict1(wchar_t *Fold, mQuantDiagr , mAppert, mLambda
   ,mAngSdvig, mNoiseSkz*mNoiseSkz, mAmplFactSig  )
{
    // ��������� �������������� �������
	TParAnt ParAnt_5P10_03(mQuantDiagr , mAppert, mLambda
   ,mAngSdvig, mNoiseSkz*mNoiseSkz, mAmplFactSig);
	double valGenAngSdvig = ParAnt_5P10_03.transformAngToGeneralizedAng (mAngSdvig ) ;
	double valStep =0.01;

	double valZero =  ParAnt_5P10_03.findFirstZero_For_SumDiagr();
	const double VAlDiap = valZero + valGenAngSdvig;  // ��������� �������� �� ����� ���������� ���� �����
	const int nBuffRows = 2. *VAlDiap/ valStep ;
	const int nBuffCols = 3;
	double  *parrBuff = new double [nBuffRows  * nBuffCols] ;
	memset(parrBuff, 0, nBuffRows * nBuffCols * sizeof(double));

	const int lenName = 30;
	wchar_t *wcharrFileNames  = new wchar_t [30 * nBuffCols];
	memset(wcharrFileNames, 0, 30 * nBuffCols * sizeof(wchar_t));

	wcscpy( wcharrFileNames, L"GenAng");
	wcscpy( &wcharrFileNames[30], L"FSum");
	wcscpy( &wcharrFileNames[60], L"Polinom4");

	double valF0 = fncSumDiagr( 0.) ;
	double root1 = findFirstZero_For_SumDiagr();
	double root2 = findSecondZero_For_SumDiagr();
   for (int i=0 ; i < nBuffRows; i++)
  {
   double valGenAngCur = VAlDiap - ((double)i ) * valStep;
   parrBuff[ i * nBuffCols] =  valGenAngCur;
   parrBuff[ i * nBuffCols +1] = ParAnt_5P10_03.fncSumDiagr( valGenAngCur)/valF0 ;
   parrBuff[ i * nBuffCols +2] = ( valGenAngCur * valGenAngCur -  root1 * root1)
		* (valGenAngCur * valGenAngCur -  root2 * root2)/ (root1 * root1 * root2 * root2);

  }


  double *pscaley = new double [nBuffCols] ;
  pscaley[0] = 1000.;
  pscaley[1] = 1.;
  pscaley[2] = 1.;

  for (int i=1 ; i < nBuffCols; i++)
  {

  TYrWriteShapeFile::WriteOneReport(Fold  // ���� � �����
								  , parrBuff // ������ � ����������� - ������� nBuffRows x nBuffCols
								  , nBuffCols // - �-�� ���������� � ������� ��������� ���������� � ������
								  , nBuffRows //  - �-�� �����
								  ,wcharrFileNames //������� � �������� ���������� - ������� nBuffCols x lenName
								  ,lenName // ������������ ����� ����� ����������
								  ,0  // ����� ���������� �� ��� X
								  ,i  // ����� ���������� �� ��� Y
								  ,1. //  ������� �� ��� X
								  ,pscaley[i]// ������� �� ��� Y
								   ) ;
  }




	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  Fold);
	wcscat(wchAxesFileName, L"AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-VAlDiap, VAlDiap + 0.5
	 ,-2.,2., 0.25) ;
	 delete pscaley ;
	 delete  parrBuff;
	 delete wcharrFileNames;

}


#pragma package(smart_init)
