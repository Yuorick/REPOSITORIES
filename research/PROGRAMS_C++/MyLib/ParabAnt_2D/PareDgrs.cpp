//---------------------------------------------------------------------------


#pragma hdrstop

#include "PareDgrs.h"
#include <math.h>
#include <string.h>

#include "Comp.h"
#include "Gauss.h"
#include "DiagrSinX.h"
#include "URPolyLine.h"
#include "UrPointXY.h"
#include "YrWriteShapeFile.h"
#include "DiagrSinX.h"
#include "MatrixProccess.h"



extern const double TET0707;
extern const double VAL_C ;



TPareDgrs::TPareDgrs()
{

 // ����� �����
   mLambda = 3. ;
 // ���� ������ ����������� ��������  �������

 mAngSdvig = 0.03;


   mSincDgr = TSincDgr();

}


// ����������� �����������
TPareDgrs::TPareDgrs (const TPareDgrs &R2)
 {
 // ����� �����
   mLambda = R2. mLambda ;
 // ���� ������ ����������� ��������  �������
  mAngSdvig =  R2.mAngSdvig;
 //
  mSincDgr = R2.mSincDgr;

 }

  // �������� ������������
  TPareDgrs TPareDgrs::operator=(TPareDgrs  R2)
{
 // ����� �����
   mLambda = R2. mLambda ;
 // ���� ������ ����������� ��������  �������
  mAngSdvig =  R2.mAngSdvig;
 //
  mSincDgr = R2.mSincDgr;
  return *this ;
}

__fastcall TPareDgrs::TPareDgrs(const double Lambda,double RadAngSdvig, TSincDgr SincDgr)
{
 mLambda = Lambda;
 mAngSdvig = RadAngSdvig;
 mSincDgr = SincDgr;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
///// ��� �������� ���� SIN(X)/X /////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

void TPareDgrs::createGraphs_For_BearingFuncs( wchar_t *wchFoldName1)
{
double valGenAngSdvig = mSincDgr.transformAngToGeneralizedAng (mAngSdvig, mLambda ) ;
 // ������� ��� ����������� ��������
	wchar_t wchFoldName[300] ={0};
	wcscpy(  wchFoldName,  wchFoldName1);
	wcscat(wchFoldName, L"\\");
	double step =0.0001;

	const double VAlDiap = M_PI;  // �������� �������� �� ����� ���������� ���� �����
								   // ��� ������ ����������� ��������
	const int nBuffRows = 2. *VAlDiap/ step ;
	const int nBuffCols = 5;
	double  *parrBuff = new double [nBuffRows  * nBuffCols] ;
	memset(parrBuff, 0, nBuffRows * nBuffCols * sizeof(double));

	const int lenName = 30;
	wchar_t *wcharrFileNames  = new wchar_t [30 * nBuffCols];
	memset(wcharrFileNames, 0, 30 * nBuffCols * sizeof(wchar_t));

	wcscpy( wcharrFileNames, L"GeneralizedTetta");

	wcscpy( &wcharrFileNames[ 30], L"fncBear_Part_ARSM");

	wcscpy( &wcharrFileNames[2 * 30], L"fncBear_Part_Division");
	wcscpy( &wcharrFileNames[3 * 30], L"fncFirstApprox_Bear_Part_ARSM");
	wcscpy( &wcharrFileNames[4* 30], L"fncThirdApprox_Bear_Part_ARSM");

  double valCoefTalor1 =  calcFirstCoefTalorBearFunc_For_PartDiagr();
  double valCoefTalor3 = calcThirdCoefTalorBearFunc_For_PartDiagr() ;
  double tet = -VAlDiap;
  for (int i=0 ; i < nBuffRows; i++)
  {
   tet = -VAlDiap + ((double)i) * step;
   double temp0 = fncDiagrSinx_div_x( tet - valGenAngSdvig/2.);
   double temp1 = fncDiagrSinx_div_x( tet + valGenAngSdvig/2.);

   parrBuff[ i * nBuffCols] = tet ; // ����� ����
   parrBuff[ i * nBuffCols + 1] = (temp0 - temp1)/ (temp0 + temp1);
   parrBuff[ i * nBuffCols + 2] =  temp0/ temp1 -1.;
   parrBuff[ i * nBuffCols + 3] =  -tet * valCoefTalor1;
   parrBuff[ i * nBuffCols + 4] =  -tet * valCoefTalor1 + valCoefTalor3* tet* tet*tet;



  }

 // double scalex = 100.;
  double *pscaley = new double [nBuffCols] ;
  pscaley[0] = 1.;
  pscaley[1] = 1.;
  pscaley[2] = 1.;
  pscaley[3] = 1.;
  pscaley[4] = 1.;
//  pscaley[5] = 10.;
//  pscaley[6] = 10.;
//  pscaley[7] = 10.;
 // pscaley[6] = 1;
 // pscaley[7] = 1;
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
								  ,1.// ������� �� ��� Y
								   ) ;
  }

  TURPointXY  arrPoints[2];
  arrPoints[0] =  TURPointXY ( -VAlDiap, 0.);
  arrPoints[1] =  TURPointXY (  VAlDiap, 0.);
  TURPolyLine plnDiap0(arrPoints,2) ;
  wchar_t wchFileName[300] ={0};
  wcscpy(  wchFileName,  wchFoldName);
  wcscat(wchFileName, L"Diap0.shp");
  TURPolyLine::WriteSetSHPFiles(wchFileName, &plnDiap0, 1) ;
  delete parrBuff;
  delete pscaley;
   ///





  ///

	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  wchFoldName);
	wcscat(wchAxesFileName, L"AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-(VAlDiap + 0.05), VAlDiap + 0.5
	 ,-3. ,3.5, 0.3) ;

  delete wcharrFileNames;
 // delete wcharrFileNames1;
}//



// ���������� �������� ���� ������� ������ ������� ��� ����  ��������
// � ���������� ����������� (������ ����������� � ����)
double  TPareDgrs::calcFirstCoefTalorBearFunc_For_PartDiagr()
{
	double valGenAngSdvig = mSincDgr.transformAngToGeneralizedAng (mAngSdvig, mLambda ) ;
	return  fncDerivDiagrSinx_div_x(valGenAngSdvig/2.)/fncDiagrSinx_div_x(valGenAngSdvig/2.) ;
}

// ���������� ������ �����������  ������ ������� ��� ����  ��������
// � ���������� ����������� � �����  tet
double  TPareDgrs::calcSecondDerivBearFunc_For_PartDiagr(const double tet)
{

	double bet = mSincDgr.transformAngToGeneralizedAng (mAngSdvig, mLambda )/2.;
	double temp0 =  fncDeriv2DiagrSinx_div_x( tet -bet );
	double temp1 =  fncDeriv2DiagrSinx_div_x( tet +bet );

	double temp2 = fncDerivDiagrSinx_div_x( tet -bet );
	double temp3 = fncDerivDiagrSinx_div_x( tet +bet );

	double temp4 = fncDiagrSinx_div_x( tet -bet );
	double temp5 = fncDiagrSinx_div_x( tet +bet );

	double tempSum =  temp4 +  temp5;

	double temp6 =  (temp0* temp5   - temp1 * temp4) * tempSum;

	double temp7 =   -2.*(temp2 * temp5 - temp3 * temp4) * ( temp2 + temp3);

	double temp8 = temp6 + temp7;

	return 2. * temp8 /(tempSum * tempSum * tempSum);
}

double  TPareDgrs::calcThirdCoefTalorBearFunc_For_PartDiagr()
{
	double temp = calcSecondDerivBearFunc_For_PartDiagr(0.01)/0.01/6.;
	return  temp ;
}

 //������ ���� ���� ���� � ���
 //INPUT:
 // cmparrPartSZv[ mQuantDiagr] - ����������� ������ ��������� ����� ��������
 // NumPartRayPare  - ����������1 ����� ���� ������� �������� (��������� ���������� � ����)
 // OUTPUT:
 // *pvalEstPartDiaRadAngTarg - ������ ���� ����, ���
 // *pcmpPartDiaKTarg  - ������ ������� ����
 // *pvalPartDiaDisp - ��������� ������ ������ ���� ����, ���*���
 // *pvalPsi2 - ������������� ������ ����� ��������������� ���������

int TPareDgrs::estimateARSMethPartDiagr( TComp *cmparrPartSZv
		  ,  double *pvalEstPartDiaRadAngTarg,  TComp *pcmpPartDiaKTarg
		   , double *pvalPartDiaDisp, double *pvalPsi2)
{
	double valGenAngSdvig = mSincDgr.transformAngToGeneralizedAng (mAngSdvig, mLambda ) ;
	// ���������� ������� ������� �������������� �������
	double val_a0 =  calcThirdCoefTalorBearFunc_For_PartDiagr(); // ��� ����� � ����
	double val_a2 =  -calcFirstCoefTalorBearFunc_For_PartDiagr(); // ��� �����
	///

	// ���������� ������ ����� ��������������� ���������
	const TComp CmpTemp = (cmparrPartSZv [ 1] - cmparrPartSZv [0] )
					/(cmparrPartSZv [ 1] + cmparrPartSZv [0] );
	double val_a3 = CmpTemp.m_Re;
	///

	double valGenAng = val_a3/val_a2;
	if(abs(valGenAng)>= M_PI)
	{
		return -1;
	}
	// ������ ���� ���� � ����� �����������
	double valF = (val_a0 * valGenAng * valGenAng * valGenAng + val_a2 * valGenAng  -val_a3);
	double valFDeriv = (3. * val_a0 * valGenAng * valGenAng  + val_a2 );
	valGenAng = valGenAng - valF / valFDeriv;
	   ///
		// ������ ���� ���� � ��������
   *pvalEstPartDiaRadAngTarg = mSincDgr.transformGeneralizedAngToAng (valGenAng , mLambda ) ;
	   ///

	   // ������ ������� ����
	   double  temp = fncDiagrSinx_div_x(valGenAng + valGenAngSdvig/2.);
	   if (valGenAng < 0.)
	   {
		 *pcmpPartDiaKTarg =   cmparrPartSZv [0]/ TComp(fncDiagrSinx_div_x(valGenAng + valGenAngSdvig/2.), 0.);
	   }
	   else
	   {
		 *pcmpPartDiaKTarg =   cmparrPartSZv [1]/ TComp(fncDiagrSinx_div_x(valGenAng - valGenAngSdvig/2.), 0.);

	   }
	   ///

	   // ������������ �������� ������� arrMtrxCorr[4] ������ ��������� ������� ��������
	   double arrMtrxCorr[16] = {0.};
	   double valGenAntpEps = 0.;
	   TComp cmpKAntp(0.,0.);
	   int lenAnsamble =  2;
	   createMtrxCorrForMeasures(valGenAng,   valGenAntpEps
	   ,*pcmpPartDiaKTarg ,cmpKAntp, arrMtrxCorr )  ;
		///

		// ���������� ��������� ������   ����������� ����������� ����


		double valS1 = (*pcmpPartDiaKTarg).modul() * fncDiagrSinx_div_x(valGenAng + valGenAngSdvig/2.) ;
		double valS2 = (*pcmpPartDiaKTarg).modul() * fncDiagrSinx_div_x(valGenAng - valGenAngSdvig/2.) ;
		double valSum2= (valS1 + valS2) * (valS1 + valS2);
		double arra[4] = {0.}, arrT[4] ={0.};
		double valDispPsi = 0.;
		arra[0] =  2. *valS2/valSum2;
		arra[2] = -2. * valS1/valSum2;
		MtrxMultMatrx(arra,1, 4, arrMtrxCorr,4, arrT) ;
		MtrxMultMatrx(arrT,1, 4, arra,1, &valDispPsi);

		double valDispGen =  valDispPsi / (valFDeriv*valFDeriv);
		///

		double valSKZ = mSincDgr.transformGeneralizedAngToAng (sqrt(valDispGen), mLambda ) ;
		 *pvalPartDiaDisp =  valSKZ * valSKZ;
		 *pvalPsi2 = fabs(CmpTemp.m_Im)/sqrt(valDispPsi);
	return 0;
}

//������������ �������������� ������� ������ ��������� ���������
// �������� ��������� ��������
// ������ ��������� ������ ����� ��� ...S[i].m_Re S , S[i].m_Im...
// INPUT:
// valGenTargEps  - ����� ���� ����
// valGenAntpEps - ����� ���� ��������
// valGenAngSdvig - ����� ���� ������ ��������
// cmpKTarg, cmpKAntp - �������(����������) ���� � ��������
//  iNumAnsamble   - ����� ������ ��������� � ��������
//          (��� ��������� �������������� � ������� ����������� ���� ������ ������� � �������� ������)
// lenAnsamble  - ���������� �������� � ��������
// OTPUT:
// arrMtrxCorr[2 *lenAnsamble *2 *lenAnsamble] - �������� �������
void TPareDgrs::createMtrxCorrForMeasures(double valGenTargEps,  double valGenAntpEps
  ,TComp cmpKTarg ,TComp cmpKAntp,  double *arrMtrxCorr )
 {
	const double  valGenAngSdvig =  mSincDgr.transformAngToGeneralizedAng (mAngSdvig, mLambda );
	memset(arrMtrxCorr, 0, 16 * sizeof(double));

  // ���������� ������ �������, ������������� ��������� ������� �������� ��������
	for (int i = 0; i < 2; i++)  // ���� �� ���������� � ��������
	{
	   // ���������� ���� ����� ������������ ��� ��������� � ������� i � ���������� �����������
	  double valDiagrGenAng = - valGenAngSdvig/2. + ((double)i) * valGenAngSdvig; // ���� ��� ����� � ������� i
	 // ((double)(-mQuantDiagr/2 + iNumAnsamble + i))*  valGenAngSdvig;
	  double valTargGen =  valGenTargEps - valDiagrGenAng; // ���� ���� � ��������� � ������� i, ���.
	  double valAntpGen =  valGenAntpEps - valDiagrGenAng; // ���� ���� � ��������� � ������� i, ���.
	  ///
	  // ��������� � ����������� ����
	  double valDiagrTarg = fncDiagrSinx_div_x(valTargGen);

	  // ��������� � �������� ��������
	  double valDiagrAntp = fncDiagrSinx_div_x(valAntpGen);

	  // ��������� �������������� ����� �������
	   double valDispRe = ( (valDiagrTarg * cmpKTarg.m_Re * valDiagrTarg * cmpKTarg.m_Re)
						  + (valDiagrAntp * cmpKAntp.m_Re * valDiagrAntp * cmpKAntp.m_Re) )
						   * mSincDgr.mAmplFactSig * mSincDgr.mAmplFactSig;
	   // ��������� ������� ����� �������
	   double valDispIm = ( (valDiagrTarg * cmpKTarg.m_Im * valDiagrTarg * cmpKTarg.m_Im)
						  + (valDiagrAntp * cmpKAntp.m_Im * valDiagrAntp * cmpKAntp.m_Im) )
						   * mSincDgr.mAmplFactSig * mSincDgr.mAmplFactSig;
	 ///
	  // ���������� ������������ ��������� ������� arrMtrxCorr � ��������������� �������
	  // �������������� � �������5���� ����� ������� ��������� � ������� i �������� ������ 2*i � 2*i+1 ��������������
	  arrMtrxCorr[2*2 * (2 * i) + 2 *i] =  valDispRe  + mSincDgr.mNoiseDisp;
	  arrMtrxCorr[2*2 * (2 * i + 1) + 2 *i + 1] =  valDispIm  + mSincDgr.mNoiseDisp;

	}
	return;
 }
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
///////////// ��� ������������ �������� ////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// ��������� �0
 TComp  TPareDgrs::fncDgrF0(const double VAlRadTetta )
 {
	 TComp cmpMod(fncModDgr0( VAlRadTetta), 0.);
	 double  valArg = fncArgDgr0( VAlRadTetta);
	 return cmpMod * exp_(valArg);
 }


// ��������� �1
 TComp  TPareDgrs::fncDgrF1(const double VAlRadTetta )
 {
	 TComp cmpMod(fncModDgr1( VAlRadTetta), 0.);
	 double  valArg = fncArgDgr1( VAlRadTetta);
	 return cmpMod * exp_(valArg);
 }

 double TPareDgrs::fncModDgr0(const double VAlRadTetta )
 {
	 double vala =  M_PI /0.0353;
	 return fncDiagrSinx_div_x(VAlRadTetta * vala );
 }

 double TPareDgrs::fncDerivModDgr0(const double VAlRadTetta )
 {
	 double vala =  M_PI /0.0353;
	 return  vala * fncDerivDiagrSinx_div_x(VAlRadTetta * vala );
 }

 double TPareDgrs::fncArgDgr0(const double VAlRadTetta )
 {
	 return 0.;
 }

 double TPareDgrs::fncModDgr1(const double VAlRadTetta )
 {
	 double vala =  M_PI /0.0353;
	 return fncDiagrSinx_div_x(VAlRadTetta * vala );
 }

 double TPareDgrs::fncArgDgr1(const double VAlRadTetta )
 {
	 return 0.;
 }

 double TPareDgrs::fncDerivModDgr1(const double VAlRadTetta )
 {
	 double vala =  M_PI /0.0353;
	 return  vala * fncDerivDiagrSinx_div_x(VAlRadTetta * vala );
 }

 double TPareDgrs::fncBearingCoeff()
 {
	 return  fncModDgr0(mAngSdvig/2. ) / fncDerivModDgr0(mAngSdvig/2.  ) ;

 }
////
// ��������� �0
 TComp  TPareDgrs::fncDgrF0(const double VAlRadTetta , TURPolyLine &plnModulGraph, TURPolyLine &plnArgGraph)
 {
	 TComp cmpMod(fncModDgr0( VAlRadTetta, plnModulGraph), 0.);
	 double  valArg = fncArgDgr0( VAlRadTetta, plnArgGraph);
	 return cmpMod * exp_(valArg);
 }

  double TPareDgrs::fncModDgr0(const double VAlRadTetta , TURPolyLine &plnGraph)
 {

	 return plnGraph.LinearValueApprox(VAlRadTetta);
 }

 double TPareDgrs::fncArgDgr0(const double VAlRadTetta , TURPolyLine &plnGraph)
 {
	 return plnGraph.LinearValueApprox(VAlRadTetta);
 }

 double TPareDgrs::fncDerivModDgr0(const double VAlRadTetta , TURPolyLine &plnGraph)
 {

	 return  (plnGraph.LinearValueApprox(VAlRadTetta + 0.00001) - plnGraph.LinearValueApprox(VAlRadTetta))/0.00001;
 }

 double TPareDgrs::fncBearingCoeff(TURPolyLine &plnGraph)
 {
	 return  fncModDgr0(mAngSdvig/2., plnGraph ) / fncDerivModDgr0(mAngSdvig/2. , plnGraph ) ;

 }
/*
// ��������� �1
 TComp  TPareDgrs::fncDgrF1(const double VAlRadTetta , TURPolyLine plnGraph)
 {
	 TComp cmpMod(fncModDgr1( VAlRadTetta), 0.);
	 double  valArg = fncArgDgr1( VAlRadTetta);
	 return cmpMod * exp_(valArg);
 }






 double TPareDgrs::fncModDgr1(const double VAlRadTetta , TURPolyLine plnGraph)
 {
	 double vala =  M_PI /0.0353;
	 return fncDiagrSinx_div_x(VAlRadTetta * vala );
 }

 double TPareDgrs::fncArgDgr1(const double VAlRadTetta , TURPolyLine plnGraph)
 {
	 return 0.;
 }

 double TPareDgrs::fncDerivModDgr1(const double VAlRadTetta , TURPolyLine plnGraph)
 {
	 double vala =  M_PI /0.0353;
	 return  vala * fncDerivDiagrSinx_div_x(VAlRadTetta * vala );
 } */
////


 // ���������� ������� ������� ������ ���� ��� ������� 2 �����
 // ������� ��������� ���� 1 �����������, ���  valRadEps
 // ������� ��������� ���� 2 ��������
 // ��� ������� ��������� ���� 2 �������������� ��������������� ������ ���������� ��������� ���� 1���
 // ���������� ����. �������� ������ ������ � ����������� �� ��������� ������ ����
 // INPUT:
 // wchFoldName1 - ���� � ����� � ���������
 // valRadEps - ������� ��������� ����1
 // valRo - ��������� ������� ���� 1 � ���� 2 (������� �������)
 // valRadDiap - �������� ��������� ��������� ���� 2
bool TPareDgrs::createGraphs_GuarantSystError_2Targs( wchar_t *wchFoldName1, double valRadEps, double valRo
  , double valRadDiap)
{
	wchar_t wchFoldName[300] ={0};
	wcscpy(  wchFoldName,  wchFoldName1);
	wcscat(wchFoldName, L"\\");


  double step = 0.001;
  int iC = 2 * valRadDiap / step -1;

	const int nBuffRows = iC ;
	const int nBuffCols = 3;
	double  *parrBuff = new double [nBuffRows  * nBuffCols] ;
	memset(parrBuff, 0, nBuffRows * nBuffCols * sizeof(double));

	const int lenName = 30;
	wchar_t *wcharrFileNames  = new wchar_t [30 * nBuffCols];
	memset(wcharrFileNames, 0, 30 * nBuffCols * sizeof(wchar_t));
	wcscpy( wcharrFileNames, L"x2");
	wcscpy( &wcharrFileNames[ 30], L"SystErr");
	wcscpy( &wcharrFileNames[ 60], L"DeltaFi");

  double alf = mAngSdvig/2.;
  const double VAlBearCoeff = fncBearingCoeff();
  for (int i =0; i < iC; i++)
  {
	double valX2 =  -valRadDiap + step/2. + ((double)i)* step;
	TComp cmpa = fncDgrF0(valX2 + alf ) - fncDgrF1(valX2 - alf );
	TComp cmpb = fncDgrF0(valRadEps + alf ) - fncDgrF1(valRadEps - alf );
	TComp cmpc = fncDgrF0(valX2 + alf ) + fncDgrF1(valX2 - alf );
	TComp cmpd = fncDgrF0(valRadEps + alf ) + fncDgrF1(valRadEps - alf );

	TComp cmpZ0(1.,0.), cmpZ1(0.,1.), cmpZ2(-1.,0.);
	TComp cmp0,cmp1, cmp2;
	bool b0 = fncLinFrac( cmpa,  cmpb,  cmpc,  cmpd,  cmpZ0, &cmp0) ;
	bool b1 = fncLinFrac( cmpa,  cmpb,  cmpc,  cmpd,  cmpZ1, &cmp1) ;
	bool b2 = fncLinFrac( cmpa,  cmpb,  cmpc,  cmpd,  cmpZ2, &cmp2) ;
	if (!(b0&&b1&&b2))
	{
      return false;
	}

	TComp cmpCentre;
	double valRadius;
	bool b3 = findCircleParams (cmp0, cmp1, cmp2, &cmpCentre, &valRadius);
	double valTemp0 = fabs((cmpCentre.m_Re + valRadius) * VAlBearCoeff - valRadEps);
	double valTemp1 = fabs((cmpCentre.m_Re - valRadius) * VAlBearCoeff- valRadEps);

	double valSystErr = valTemp0;
	TComp cmpOmega (cmpCentre.m_Re + valRadius, cmpCentre.m_Im);
	if (valTemp1 > valTemp0)
	{
	  valSystErr = valTemp1;
	  cmpOmega = TComp  (cmpCentre.m_Re - valRadius, cmpCentre.m_Im);
	}
	TComp cmpProImage;
	fncLinFrac( cmpd * TComp(-1.,0.), cmpb,  cmpc,  cmpa * TComp(-1.,0.),  cmpOmega, &cmpProImage) ;
	parrBuff[ i * nBuffCols] = valX2 ; //
	parrBuff[ i * nBuffCols + 1] = valSystErr;
	parrBuff[ i * nBuffCols + 2] = cmpProImage.phase() * 180. / M_PI;

  }

   double *pscaley = new double [nBuffCols] ;
  pscaley[0] = 1000.;
  pscaley[1] = 1000.;
  pscaley[2] = 1.;

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
								  ,pscaley[0] //  ������� �� ��� X
								  ,pscaley[i]// ������� �� ��� Y
								   ) ;
  }

  delete parrBuff;
  delete pscaley;

	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  wchFoldName);
	wcscat(wchAxesFileName, L"AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-(valRadDiap + 0.05) *pscaley[0], (valRadDiap + 0.5) *pscaley[0]
	 ,-10. ,10, 0.3) ;

  delete wcharrFileNames;
}


 // ���������� ������� ������� ������ ���� ��� ������� 2 �����
 // ������� ��������� ���� 1 �����������, ���  valRadEps
 // ������� ��������� ���� 2 ��������
 // ��� ������� ��������� ���� 2 �������������� ��������������� ������ ���������� ��������� ���� 1���
 // ���������� ����. �������� ������ ������ � ����������� �� ��������� ������ ����
 // INPUT:
 // wchFoldName1 - ���� � ����� � ���������
 // valRadEps - ������� ��������� ����1
 // valRo1,valRo2 - ��������� ������jd ���� 1 � ���� 2
 // valRadDiap - �������� ��������� ��������� ���� 2
bool TPareDgrs::createGraphs_GuarantSystError_2Targs( wchar_t *wchFoldName1
  , TURPolyLine &plnModulGraph0, TURPolyLine &plnModulArg0, TURPolyLine &plnModulGraph1, TURPolyLine &plnModulArg1
  ,double valRo1, double valRo2, double valRadEps,  double valRadDiap)
{
	wchar_t wchFoldName[300] ={0};
	wcscpy(  wchFoldName,  wchFoldName1);
	wcscat(wchFoldName, L"\\");


  double step = 0.001;
  int iC = 2 * valRadDiap / step -1;

	const int nBuffRows = iC ;
	const int nBuffCols = 3;
	double  *parrBuff = new double [nBuffRows  * nBuffCols] ;
	memset(parrBuff, 0, nBuffRows * nBuffCols * sizeof(double));

	const int lenName = 30;
	wchar_t *wcharrFileNames  = new wchar_t [30 * nBuffCols];
	memset(wcharrFileNames, 0, 30 * nBuffCols * sizeof(wchar_t));
	wcscpy( wcharrFileNames, L"x2");
	wcscpy( &wcharrFileNames[ 30], L"SystErr");
	wcscpy( &wcharrFileNames[ 60], L"DeltaFi");

  double alf = mAngSdvig/2.;
  const double VAlBearCoeff = fncBearingCoeff();
  int iCur = 0;
  for (int i =0; i < iC; i++)
  {
	double valX2 =  -valRadDiap + step/2. + ((double)i)* step;
	TComp cmpa = TComp(valRo2/ valRo1, 0.)
		 * (fncDgrF0(valX2 + alf,plnModulGraph0,plnModulArg0 ) - fncDgrF0(valX2 - alf ,plnModulGraph1,plnModulArg1 ));
	TComp cmpb = fncDgrF0(valRadEps + alf,plnModulGraph0,plnModulArg0 ) - fncDgrF0(valRadEps - alf,plnModulGraph1,plnModulArg1 );
	TComp cmpc = TComp(valRo2/ valRo1, 0.)
	 * (fncDgrF0(valX2 + alf,plnModulGraph0,plnModulArg0 ) + fncDgrF0(valX2 - alf,plnModulGraph1,plnModulArg1 ));
	TComp cmpd = fncDgrF0(valRadEps + alf,plnModulGraph0,plnModulArg0 ) + fncDgrF0(valRadEps - alf,plnModulGraph1,plnModulArg1 );

	TComp cmpZ0(1.,0.), cmpZ1(0.,1.), cmpZ2(-1.,0.);
	TComp cmp0,cmp1, cmp2;
	bool b0 = fncLinFrac( cmpa,  cmpb,  cmpc,  cmpd,  cmpZ0, &cmp0) ;
	bool b1 = fncLinFrac( cmpa,  cmpb,  cmpc,  cmpd,  cmpZ1, &cmp1) ;
	bool b2 = fncLinFrac( cmpa,  cmpb,  cmpc,  cmpd,  cmpZ2, &cmp2) ;
	if (!(b0&&b1&&b2))
	{
	  continue;
	}

	TComp cmpCentre;
	double valRadius;
	bool b3 = findCircleParams (cmp0, cmp1, cmp2, &cmpCentre, &valRadius);
	double valTemp0 = fabs((cmpCentre.m_Re + valRadius) * VAlBearCoeff - valRadEps);
	double valTemp1 = fabs((cmpCentre.m_Re - valRadius) * VAlBearCoeff- valRadEps);

	double valSystErr = valTemp0;
	TComp cmpOmega (cmpCentre.m_Re + valRadius, cmpCentre.m_Im);
	if (valTemp1 > valTemp0)
	{
	  valSystErr = valTemp1;
	  cmpOmega = TComp  (cmpCentre.m_Re - valRadius, cmpCentre.m_Im);
	}
	TComp cmpProImage;
	if(!fncLinFrac( cmpd * TComp(-1.,0.), cmpb,  cmpc,  cmpa * TComp(-1.,0.),  cmpOmega, &cmpProImage))
	{
	 continue;
	}
	parrBuff[ i * nBuffCols] = valX2 ; //
	parrBuff[ i * nBuffCols + 1] = valSystErr;
	parrBuff[ i * nBuffCols + 2] = cmpProImage.phase() ;
	if (fabs(parrBuff[ i * nBuffCols + 2] + M_PI) < 0.0001 )
	{
	  parrBuff[ i * nBuffCols + 2] = M_PI;
	}
	iCur++;

  }

   double *pscaley = new double [nBuffCols] ;
  pscaley[0] = 1000.;
  pscaley[1] = 1000.;
  pscaley[2] = 1.;

  for (int i=1; i < nBuffCols; i++)
  {

  TYrWriteShapeFile::WriteOneReport(                 wchFoldName  // ���� � �����
								  , parrBuff // ������ � ����������� - ������� nBuffRows x nBuffCols
								  , nBuffCols // - �-�� ���������� � ������� ��������� ���������� � ������
								  , iCur //  - �-�� �����
								  ,wcharrFileNames //������� � �������� ���������� - ������� nBuffCols x lenName
								  ,lenName // ������������ ����� ����� ����������
								  ,0  // ����� ���������� �� ��� X
								  ,i  // ����� ���������� �� ��� Y
								  ,pscaley[0] //  ������� �� ��� X
								  ,pscaley[i]// ������� �� ��� Y
								   ) ;
  }

  delete parrBuff;
  delete pscaley;

	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  wchFoldName);
	wcscat(wchAxesFileName, L"AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-(valRadDiap + 0.05) *pscaley[0], (valRadDiap + 0.5) *pscaley[0]
	 ,-500. ,500., 0.3) ;

  delete wcharrFileNames;
}

/*
// ���������� ������� ������� ������ ���� ������ ��� ������� 2 �����
 // ������� ��������� ���� 1 �����������, ���  valRadEps
 // ������� ��������� ���� 2 ��������
 // ��� ������� ��������� ���� 2 �������������� ��������������� ������ ���������� ��������� ���� 1���
 // ���������� ����. �������� ������ ������ � ����������� �� ��������� ������ ����
 // INPUT:
 // wchFoldName1 - ���� � ����� � ���������
 // valRadEps - ������� ��������� ����1
 // valRo1,valRo2 - ��������� ������jd ���� 1 � ���� 2
 // valRadDiap - �������� ��������� ��������� ���� 2
bool TPareDgrs::createGraphs_GuarantSystErrARSM_Mod_2Targs( wchar_t *wchFoldName1
  , TURPolyLine &plnModulGraph0, TURPolyLine &plnModulArg0, TURPolyLine &plnModulGraph1, TURPolyLine &plnModulArg1
  ,double valRo1, double valRo2, double valRadEps,  double valRadDiap)
{
	wchar_t wchFoldName[300] ={0};
	wcscpy(  wchFoldName,  wchFoldName1);
	wcscat(wchFoldName, L"\\");


  double step = 0.001;
  int iC = 2 * valRadDiap / step -1;

	const int nBuffRows = iC ;
	const int nBuffCols = 3;
	double  *parrBuff = new double [nBuffRows  * nBuffCols] ;
	memset(parrBuff, 0, nBuffRows * nBuffCols * sizeof(double));

	const int lenName = 30;
	wchar_t *wcharrFileNames  = new wchar_t [30 * nBuffCols];
	memset(wcharrFileNames, 0, 30 * nBuffCols * sizeof(wchar_t));
	wcscpy( wcharrFileNames, L"x2");
	wcscpy( &wcharrFileNames[ 30], L"SystErr");
	wcscpy( &wcharrFileNames[ 60], L"DeltaFi");

  double alf = mAngSdvig/2.;
  const double VAlBearCoeff = fncBearingCoeff();
  int iCur = 0;
  for (int i =0; i < iC; i++)
  {
	double valX2 =  -valRadDiap + step/2. + ((double)i)* step;
	TComp cmpa = TComp(valRo1/ valRo1, 0.) * fncDgrF0(valRadEps - alf ,plnModulGraph1,plnModulArg1 );
		// * (fncDgrF0(valX2 + alf,plnModulGraph0,plnModulArg0 ) - fncDgrF0(valX2 - alf ,plnModulGraph1,plnModulArg1 ));


	TComp cmpb = TComp(valRo2/ valRo1, 0.) * fncDgrF0(valX2 - alf ,plnModulGraph1,plnModulArg1 );
	//fncDgrF0(valRadEps + alf,plnModulGraph0,plnModulArg0 ) - fncDgrF0(valRadEps - alf,plnModulGraph1,plnModulArg1 );
	TComp cmpc =  TComp(valRo1/ valRo1, 0.) * fncDgrF0(valRadEps + alf ,plnModulGraph0,plnModulArg0 );
	//TComp(valRo2/ valRo1, 0.)
	 //* (fncDgrF0(valX2 + alf,plnModulGraph0,plnModulArg0 ) + fncDgrF0(valX2 - alf,plnModulGraph1,plnModulArg1 ));
	TComp cmpd =  TComp(valRo2/ valRo1, 0.) * fncDgrF0(valX2 + alf ,plnModulGraph0,plnModulArg0 );
	//fncDgrF0(valRadEps + alf,plnModulGraph0,plnModulArg0 ) + fncDgrF0(valRadEps - alf,plnModulGraph1,plnModulArg1 );

	TComp cmpZ0(1.,0.), cmpZ1(0.,1.), cmpZ2(-1.,0.);
	TComp cmp0,cmp1, cmp2;
	bool b0 = fncLinFrac( cmpa,  cmpb,  cmpc,  cmpd,  cmpZ0, &cmp0) ;
	bool b1 = fncLinFrac( cmpa,  cmpb,  cmpc,  cmpd,  cmpZ1, &cmp1) ;
	bool b2 = fncLinFrac( cmpa,  cmpb,  cmpc,  cmpd,  cmpZ2, &cmp2) ;
	if (!(b0&&b1&&b2))
	{
	  continue;
	}

	TComp cmpCentre;
	double valRadius;
	bool b3 = findCircleParams (cmp0, cmp1, cmp2, &cmpCentre, &valRadius);
	double valTemp0 = fabs((cmpCentre.m_Re + valRadius) * VAlBearCoeff - valRadEps);
	double valTemp1 = fabs((cmpCentre.m_Re - valRadius) * VAlBearCoeff- valRadEps);

	double valSystErr = valTemp0;
	TComp cmpOmega (cmpCentre.m_Re + valRadius, cmpCentre.m_Im);
	if (valTemp1 > valTemp0)
	{
	  valSystErr = valTemp1;
	  cmpOmega = TComp  (cmpCentre.m_Re - valRadius, cmpCentre.m_Im);
	}
	TComp cmpProImage;
	if(!fncLinFrac( cmpd * TComp(-1.,0.), cmpb,  cmpc,  cmpa * TComp(-1.,0.),  cmpOmega, &cmpProImage))
	{
	 continue;
	}
	parrBuff[ i * nBuffCols] = valX2 ; //
	parrBuff[ i * nBuffCols + 1] = valSystErr;
	parrBuff[ i * nBuffCols + 2] = cmpProImage.phase() ;
	if (fabs(parrBuff[ i * nBuffCols + 2] + M_PI) < 0.0001 )
	{
	  parrBuff[ i * nBuffCols + 2] = M_PI;
	}
	iCur++;

  }

   double *pscaley = new double [nBuffCols] ;
  pscaley[0] = 1000.;
  pscaley[1] = 1000.;
  pscaley[2] = 1.;

  for (int i=1; i < nBuffCols; i++)
  {

  TYrWriteShapeFile::WriteOneReport(                 wchFoldName  // ���� � �����
								  , parrBuff // ������ � ����������� - ������� nBuffRows x nBuffCols
								  , nBuffCols // - �-�� ���������� � ������� ��������� ���������� � ������
								  , iCur //  - �-�� �����
								  ,wcharrFileNames //������� � �������� ���������� - ������� nBuffCols x lenName
								  ,lenName // ������������ ����� ����� ����������
								  ,0  // ����� ���������� �� ��� X
								  ,i  // ����� ���������� �� ��� Y
								  ,pscaley[0] //  ������� �� ��� X
								  ,pscaley[i]// ������� �� ��� Y
								   ) ;
  }

  delete parrBuff;
  delete pscaley;

	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  wchFoldName);
	wcscat(wchAxesFileName, L"AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-(valRadDiap + 0.05) *pscaley[0], (valRadDiap + 0.5) *pscaley[0]
	 ,-500. ,500., 0.3) ;

  delete wcharrFileNames;
}
 */
double TPareDgrs::max_(double a, double b)
{
	return (a>b)?a:b;
}

void TPareDgrs::imitateMeasureArr(TComp cmpK, double valEps, TComp *cmparrMeas)
{
  cmparrMeas[0] =  cmpK * fncDgrF0(valEps + mAngSdvig/2. );
  cmparrMeas[1] =  cmpK * fncDgrF1(valEps - mAngSdvig/2. );
}

double TPareDgrs::estimateARSM(TComp *cmparrMeas)
{
   double valBearCoeff =  fncBearingCoeff() ;
   TComp cmpQZv = (cmparrMeas[0] - cmparrMeas[1])/ (cmparrMeas[0] + cmparrMeas[1]) ;
   double valReal = cmpQZv.m_Re;
   return  valReal * valBearCoeff;
}

#pragma package(smart_init)
