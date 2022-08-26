//---------------------------------------------------------------------------


#pragma hdrstop
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "YrWriteShapeFile.h"
#include "URPolyLine.h"
#include "UrPointXY.h"
#include "Diagrams.h"
#include "Equations.h"
#include "MatrixProccess.h"
#include "Comp.h"
#include "DiagrSinX.h"
extern const double TET0707 ;


// ���������� �� ���


// ���������� �� ���   �������������
// val_d -���������� ����� ������������ � ������,
// val_D - ���������� ����� ��������
// ival_n - ���������� ����������� � ������
// ival_N  � ���������� �������
// valTetScan  � ���� ������������.
// valLamb  - ����� �����
// tet - ���� ��������
double fncDiagrPar(const double val_d,const  double val_D, const int ival_n
  , const  int ival_N,const  double valLamb, const double valTetScan,const double tet)
{
double valk = 2. * M_PI / valLamb;
double alf = valk * val_d * ((double)ival_n)/ 2.;
double bet = valk * val_d / 2.;
double gam = valk * val_D * ((double)ival_N)/ 2.;
double del = valk * val_D / 2.;
//double mu = sin(valTetScan);

  if (valTetScan > 0.000001)
  {
  if (fabs(tet-valTetScan)< 0.000001)
  {
   double bb = sqrt(cos(valTetScan)) * sin(alf * sin(valTetScan))/sin(bet * sin(valTetScan)) /((double)ival_n) ;
   return bb;
  }
  if (fabs(tet)< 0.0000001)
  {
   double aa = sin( gam * ( sin( valTetScan)))/ sin( del * ( sin( valTetScan)))/((double)ival_N);
   return aa ;
  }
  }
  else
  {
     if (fabs(tet)< 0.0000001) return 1.;
  }

  double temp = sqrt(cos(tet))* sin( alf * sin(tet) )
		   /(((double)ival_n)* sin( bet * sin(tet) )) ;
  double temp1 =  sin( gam * (sin(tet) - sin( valTetScan)) )
		   / sin( del* (sin(tet) - sin(valTetScan)) )/((double)ival_N) ;
  return temp * temp1 ;
}

// ����������� ������� ���������� �� ���   �� �����
// val_d -���������� ����� ������������ � ������,
// val_D - ���������� ����� ��������
// ival_n - ���������� ����������� � ������
// ival_N  � ���������� �������
// valLamb � ���� ������������.
// tet - ���� ��������
double fnc_dDiagrPar_po_dTet(const double val_d,const  double val_D, const int ival_n
  , const  int ival_N,const  double valLamb, const double valTetScan,const double tet)
{
  if (fabs(tet)< 0.0000001) return 0.;
  double valk = 2. * M_PI / valLamb;
double alf = valk * val_d * ((double)ival_n)/ 2.;
double bet = valk * val_d / 2.;
double gam = valk * val_D * ((double)ival_N)/ 2.;
double del = valk * val_D / 2.;
double mu = sin(valTetScan);

double val_F1 = sqrt(cos(tet));
double val_dF1 = - sin (tet)/val_F1/2.;
double f31 = fnc_F3( alf,0.,  tet)  ;
double f32 = fnc_F3( bet,0.,  tet)  ;
double f33 = fnc_F3( gam, mu,  tet)  ;
double f34 = fnc_F3( del, mu,  tet)  ;

double val_F21 = f31/ f32;
double val_F22 = f33/ f34;

double a1 =  fnc_dF3_po_tet(alf, 0. ,  tet) ;
double a2 =  fnc_dF3_po_tet(bet, 0. ,  tet) ;
double b1 =  fnc_dF3_po_tet(gam, mu ,  tet) ;
double b2 =  fnc_dF3_po_tet(del, mu ,  tet) ;

double val_dF21_po_dTet =  fnc_F4(f31 ,f32, a1, a2) ;
double val_dF22_po_dTet =  fnc_F4(f33 ,f34, b1, b2) ;

return ( val_dF1 *  val_F21 * val_F22 + val_F1 * val_dF21_po_dTet * val_F22
	 + val_F1 * val_F21 * val_dF22_po_dTet)/ ((double ) ival_n)/ ((double ) ival_N);

}
/// dF3 = dF3/dTet
double fnc_dF3_po_tet(const double val_gam,const  double val_mu, const double tet)
{
   return cos(val_gam * (sin(tet)-val_mu) ) * val_gam * cos(tet);
}

/// F3 = sin(gam(sin(tet)-mu) -
double fnc_F3(const double val_gam,const  double val_mu, const double tet)
{
   return sin(val_gam * (sin(tet)-val_mu) );
}
//-----------------------------------------------------------------
/// F3 = sin(gam(sin(tet)-mu) -
double fnc_F4(const double v,const  double u, const double v1,const  double u1)
{
   return (v1 * u - u1 * v)/ u / u;
}
//-----------------------------------------------------------------

double findDiagrWidth(const double val_d,const  double val_D, const int ival_n
  , const  int ival_N,const  double valLamb)
{
 //double tet0 = 0.02;
 double tet0 = valLamb/val_D/ ival_N * 0.7;

 int i =0;
 double a = sqrt(2.)/2.;
 for ( i =0; i < 10; i++)
 {
  double del = -(fncDiagrPar( val_d, val_D,  ival_n, ival_N, valLamb, 0., tet0) - a)/ fnc_dDiagrPar_po_dTet( val_d, val_D,  ival_n, ival_N, valLamb, 0., tet0);
  tet0 += del;
  if (fabs(del) < 0.0000001) break;

 }
 return tet0;
}
//-----------------------------------------------------------------
/*
void createDiagrGraphs(double valLamb,double valWidthDgr,wchar_t *wchFoldName1 )
{
  wchar_t wchFoldName[300] ={0};
  wcscpy(  wchFoldName,  wchFoldName1);
  wcscat(wchFoldName, L"\\");
  double step = M_PI / 3000./ 10.;
  const int nBuffRows = 470 *2;
  double  *parrBuff = new double [nBuffRows  * 3] ;
  memset(parrBuff, 0, nBuffRows * 3 * sizeof(double));
  const int nBuffCols =3;
  const int lenName = 30;
  wchar_t wcharrFileNames [90] ={0};
  wcscpy( wcharrFileNames, L"Tetta");
  wcscpy( &wcharrFileNames[30], L"SimpleDiagr");
  wcscpy( &wcharrFileNames[60], L"ParDiagr");

  for (int i=0 ; i < nBuffRows; i++)
  {
   parrBuff[ i * nBuffCols] = ((double) (-nBuffRows/2 +i)) / 10.;
   double tet = step * ((double) (-nBuffRows/2 +i));
   parrBuff[ i * nBuffCols + 1]=fncDiagrSimple(valWidthDgr, tet * 1000.) ;
   parrBuff[ i * nBuffCols + 2]=fncDiagrPar( valLamb,0., tet);

  }

  double scaley = 100.;
  TYrWriteShapeFile::WriteOneReport(                 wchFoldName  // ���� � �����
								  , parrBuff // ������ � ����������� - ������� nBuffRows x nBuffCols
								  , nBuffCols // - �-�� ���������� � ������� ��������� ���������� � ������
								  , nBuffRows //  - �-�� �����
								  ,wcharrFileNames //������� � �������� ���������� - ������� nBuffCols x lenName
								  ,lenName // ������������ ����� ����� ����������
								  ,0  // ����� ���������� �� ��� X
								  ,1  // ����� ���������� �� ��� Y
								  ,1. //  ������� �� ��� X
								  ,scaley // ������� �� ��� Y
								   ) ;
TYrWriteShapeFile::WriteOneReport(                 wchFoldName  // ���� � �����
								  , parrBuff // ������ � ����������� - ������� nBuffRows x nBuffCols
								  , nBuffCols // - �-�� ���������� � ������� ��������� ���������� � ������
								  , nBuffRows //  - �-�� �����
								  ,wcharrFileNames //������� � �������� ���������� - ������� nBuffCols x lenName
								  ,lenName // ������������ ����� ����� ����������
								  ,0  // ����� ���������� �� ��� X
								  ,2  // ����� ���������� �� ��� Y
								  ,1. //  ������� �� ��� X
								  ,scaley  // ������� �� ��� Y
								   ) ;
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

} */
//-----------------------------------------------------------------
void createDiagrGraphs_v1(const double val_d,const  double val_D, const int ival_n
  , const  int ival_N,const  double valLamb,const double valTetScan,wchar_t *wchFoldName1 )
{
  wchar_t wchFoldName[300] ={0};
  wcscpy(  wchFoldName,  wchFoldName1);
  wcscat(wchFoldName, L"\\");
  double step = M_PI / 3000./ 10.;
  const int nBuffRows = 470 *2;
  const int nBuffCols =4;
  double  *parrBuff = new double [nBuffRows  * nBuffCols] ;
  memset(parrBuff, 0, nBuffRows * nBuffCols * sizeof(double));

  const int lenName = 30;
  wchar_t wcharrFileNames [120] ={0};
  wcscpy( wcharrFileNames, L"Tetta");
  wcscpy( &wcharrFileNames[30], L"SimpleDiagr");
  wcscpy( &wcharrFileNames[60], L"ParDiagr");
  wcscpy( &wcharrFileNames[90], L"ParDiagr_with_TetScan");

  double valWidthDgr = findDiagrWidth( val_d, val_D,  ival_n
  ,  ival_N, valLamb);
  for (int i=0 ; i < nBuffRows; i++)
  {
   double tet = step * ((double) (-(double)nBuffRows/2. + (double)i));
   parrBuff[ i * nBuffCols] = tet ;
   parrBuff[ i * nBuffCols + 1]=fncDiagrSimple(valWidthDgr, tet ) ;
   parrBuff[ i * nBuffCols + 2]=fncDiagrPar( val_d, val_D,  ival_n
  ,  ival_N,valLamb,0., tet);
   parrBuff[ i * nBuffCols + 3]=fncDiagrPar( val_d, val_D,  ival_n
  ,  ival_N,valLamb,valTetScan, tet);


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
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////   ����� ������� ��� ������������ ������ � 6 ������������   ////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
// ������� �������� ��������� ������������ ����� ���� � �������� ������� �������
// ��� ������������� �����
// INPUT:
// quantDiagr  - �-�� ��������
// valWidthDgr  - ������ ��������� �� ������ 0707 � ���
// pcmpSZv   - ������ ��������� �� ���������� (����������� �����)
// pAlf     - ���� ���������� ��������� �� ������������
// fncF    - ������� ���������
// fncdF_po_dTet - ������� ������ ����������� ���������
// fncd2F_po_dTet2  - ������� ������  ����������� ���������
// OUTPUT:
// pcmpZ1  - ����������� ����������� ��������� ����
// pcmpZ2  - ����������� ����������� ��������� ��������
// palfTrg - ���� ����� ����
// palfAnp - ���� ����� ��������
// ����������:
//  -1 - ���� �-�� �������� quantDiagr< 3
// -2 - ����� ������� ����� ���������
//
int find_Targ_and_Antp_Newton6(const int quantDiagr,double valWidthDgr,TComp *pcmpSZv, double *pAlf,
 double(*pfncF)(double valWidthDgr, double tet)
 , double(*pfncdF_po_dTet)(double valWidthDgr, double tet), double(*pfncd2F_po_dTet2)(double valWidthDgr, double tet),
 TComp *pcmpZ1, TComp *pcmpZ2, double *palfTrg, double *palfAnp )
{
if (quantDiagr < 3)
{
return -1;
}
// ��������� �����������:
double arrX[6] ={0.};

arrX[0] =  *palfTrg ;
arrX[1] = *palfAnp ;
arrX[2] = (*pcmpZ1).m_Re  ;   // ����� ��������� ����
arrX[3] = (*pcmpZ2).m_Re  ;   // ����� ��������� ����
arrX[4] = (*pcmpZ1).m_Im  ; // ����� ��������� ��������
arrX[5] = (*pcmpZ2).m_Im  ; // ����� ��������� ��������

//
// ��������� ������� ���������
TComp pcmpZ[2];
if(!findOtimal_Z1_and_Z2(quantDiagr,valWidthDgr,pcmpSZv, pAlf,
pfncF, pcmpZ, arrX )) return -3 ;
arrX[2] = pcmpZ[0].m_Re  ;   // ����� ��������� ����
arrX[3] = pcmpZ[1].m_Re  ;   // ����� ��������� ����
arrX[4] = pcmpZ[0].m_Im  ; // ����� ��������� ��������
arrX[5] = pcmpZ[1].m_Im  ; // ����� ��������� ��������


double arrDelX[6] ={0.};
int i = 0;
double arr_FGreek[6] = {0.}, arr_dFGreek[36] ={0.}, arr_dFGreekInv[36] ={0.};
int ireturn = -2;
double del = -2.;
for (i = 0; i < 100; i++)
{
	if(!findOtimal_Z1_and_Z2(quantDiagr,valWidthDgr,pcmpSZv, pAlf,
	pfncF, pcmpZ, arrX )) return -3 ;
	arrX[2] = pcmpZ[0].m_Re  ;   // ����� ��������� ����
	arrX[3] = pcmpZ[1].m_Re  ;   // ����� ��������� ����
	arrX[4] = pcmpZ[0].m_Im  ; // ����� ��������� ��������
	arrX[5] = pcmpZ[1].m_Im  ; // ����� ��������� ��������

	///
	calcF_and_dF_Newton6( quantDiagr, valWidthDgr,pcmpSZv, pAlf,
	pfncF,pfncdF_po_dTet, pfncd2F_po_dTet2, arrX, arr_FGreek, arr_dFGreek ) ;

	// ��� ��������
	double arrX1[6] ={0.}, arr_dFGreek1[36] ={0.}, arr_FGreek1[6] ={0.};
	memcpy(arrX1,arrX, 6 * sizeof(double));
	arrX1 [5] += 0.0001;
	calcF_and_dF_Newton6( quantDiagr, valWidthDgr,pcmpSZv, pAlf,
	pfncF,pfncdF_po_dTet, pfncd2F_po_dTet2, arrX1, arr_FGreek1, arr_dFGreek1 ) ;
	double arrDelX1 [6] ={0.}, arrDelX2 [6] ={0.};
	MtrxMinusMatrx(arr_FGreek1, arr_FGreek,6, 1, arrDelX1);
	MatrxMultScalar(arrDelX1, 6, 1, 1. /0.0001,arrDelX2);

	//
	if(!InverseMtrx6(arr_dFGreek, arr_dFGreekInv)) return -2;
	MtrxMultMatrx(arr_dFGreekInv ,6, 6, arr_FGreek,1, arrDelX1) ;
	del = NormVect(arrDelX1, 2);
	MatrxMultScalar(arrDelX1, 1, 6, 1. ,arrDelX);
	double arrT[6] ={0.};
	MtrxMinusMatrx(arrX, arrDelX,1, 6, arrT) ;
	memcpy( arrX, arrT, 6 * sizeof(double));
	if (fabs(arrX[0])> 12.) arrX[0] = 12. * SIGN__(arrX[0]);
	if (fabs(arrX[1])> 12.) arrX[1] = 12. * SIGN__(arrX[1]);


	if (del< 0.001)
	{
		ireturn = 0;
		if (arrX[0] > arrX[1])
		{
		(*pcmpZ1).m_Re = arrX[2];
		(*pcmpZ1).m_Im = arrX[4];
		(*pcmpZ2).m_Re = arrX[3];
		(*pcmpZ2).m_Im = arrX[5];
		*palfTrg =  arrX[0];
		*palfAnp =  arrX[1];
		}
		else
		{
		(*pcmpZ2).m_Re = arrX[2];
		(*pcmpZ2).m_Im = arrX[4];
		(*pcmpZ1).m_Re = arrX[3];
		(*pcmpZ1).m_Im = arrX[5];
		*palfTrg =  arrX[1];
		*palfAnp =  arrX[0];
		}

		break ;
	}
}
return ireturn;
}

void calcF_and_dF_Newton6(const int quantDiagr,double valWidthDgr,TComp *pcmpSZv, double *pAlf,
 double(*pfncF)(double valWidthDgr, double tet)
 , double(*pfncdF_po_dTet)(double valWidthDgr, double tet), double(*pfncd2F_po_dTet2)(double valWidthDgr, double tet),
 double *arrX, double * arr_F, double * arr_dF )
 {
  double arrFT[6] ={0.}, arrdFT[36] ={0.}, arrFT0[6] ={0.}, arrdFT0[36] ={0.};
  memset( arr_F, 0, 6 * sizeof(double));
  memset( arr_dF, 0, 36 * sizeof(double));
  for (int i = 0; i < quantDiagr; i++)
  {
	calcPartial_F_and_dF_Newton6( valWidthDgr,pcmpSZv[i], pAlf[i],
	 pfncF, pfncdF_po_dTet,pfncd2F_po_dTet2,
	   arrX,arrFT, arrdFT ) ;
	 MtrxSumMatrx(arr_F, arrFT,1, 6, arrFT0) ;
	 memcpy(arr_F, arrFT0, 6 * sizeof(double));
	 MtrxSumMatrx(arr_dF, arrdFT,6, 6, arrdFT0) ;
	 memcpy(arr_dF, arrdFT0, 36 * sizeof(double));
  }
 }

 // ���������� ������� ������� ������� �������� ���������  ��� ���������� �����
 // ���� � �������� ��� ������������� ������� ���������
 // INPUT:
// quantDiagr  - �-�� ��������
// valWidthDgr  - ������ ��������� �� ������ 0707 � ���
// pcmpSZv   - ������ ��������� �� ���������� (����������� �����)
// pAlf     - ���� ���������� ��������� �� ������������
// fncF    - ������� ���������
// fncdF_po_dTet - ������� ������ ����������� ���������
// fncd2F_po_dTet2  - ������� ������  ����������� ���������


 // arrX[0] - ���� ����� ����
 // arrX[1] -  ���� ����� ��������
 // arrX[2] -k11 �������������� �����  ������������ ��������� ����
 // arrX[3] -k21 �������������� �����  ������������ ��������� ���������
 // arrX[4] - k12 ������ �����  ������������ ��������� ����
 // arrX[5] - k22 ������ �����  ������������ ��������� ���������
 // OUTPUT:
 // arr_F[2]  - F(X)
 // arr_dF[4] -dF(X)/dX

 void calcF_and_dF_AnglesRelax_Newton6(const int quantDiagr,double valWidthDgr,TComp *pcmpSZv, double *pAlf,
 double(*pfncF)(double valWidthDgr, double tet)
 , double(*pfncdF_po_dTet)(double valWidthDgr, double tet), double(*pfncd2F_po_dTet2)(double valWidthDgr, double tet),
 double *arrX, double * arr_F, double * arr_dF )
 {
  double arrFT[2] ={0.}, arrdFT[4] ={0.}, arrFT0[2] ={0.}, arrdFT0[4] ={0.};
  memset( arr_F, 0, 2 * sizeof(double));
  memset( arr_dF, 0, 4 * sizeof(double));
  for (int i = 0; i < quantDiagr; i++)
  {
	calcPartial_F_and_dF_AnglesRelax_Newton6( valWidthDgr,pcmpSZv[i], pAlf[i],
	 pfncF, pfncdF_po_dTet,pfncd2F_po_dTet2,
	   arrX,arrFT, arrdFT ) ;
	 MtrxSumMatrx(arr_F, arrFT,1, 2, arrFT0) ;
	 memcpy(arr_F, arrFT0, 2 * sizeof(double));
	 MtrxSumMatrx(arr_dF, arrdFT,2, 2, arrdFT0) ;
	 memcpy(arr_dF, arrdFT0, 4 * sizeof(double));
  }
 }
 // ������ ������ ������� F(X)=dFGreek(X)/dX � ������� ����� dF(X)/dX ��� ���� ���������
 // ��� ���������� ����������� ����� ��� ������������� ����� ���������
 // INPUT:
 // valWidthDgr - ������ ��������� �� ������ 0707 � ���
 // cmpSZv  - ���������, ����������� �����
 // valAlf  - ���� ���������� ��������� �� �����������(���� �����)
 // fncF    - ������� ���������
// fncdF_po_dTet - ������� ������ ����������� ���������
// fncd2F_po_dTet2  - ������� ������  ����������� ���������
 // arrX[6] -  ������ ����������, � ������� ������������ ����������
 // arrX[0] - ���� ����� ����
 // arrX[1] -  ���� ����� ��������
 // arrX[2] -k11 �������������� �����  ������������ ��������� ����
 // arrX[3] -k21 �������������� �����  ������������ ��������� ���������
 // arrX[4] - k12 ������ �����  ������������ ��������� ����
 // arrX[5] - k22 ������ �����  ������������ ��������� ���������
 // OUTPUT:
 // arr_F[2]  - F(X)
 // arr_dF[4] -dF(X)/dX
 void calcPartial_F_and_dF_AnglesRelax_Newton6(double valWidthDgr,TComp cmpSZv, double valAlf0,
 double(*pfncF)(double valWidthDgr, double tet)
 , double(*pfncdF_po_dTet)(double valWidthDgr, double tet), double(*pfncd2F_po_dTet2)(double valWidthDgr, double tet),
 double *arrX, double * arr_F, double * arr_dF )
 {
  double valAlf = valAlf0 ;

  memset( arr_F, 0, 2 * sizeof(double));
  memset( arr_dF, 0, 4 * sizeof(double));
  double val_F1 = (*pfncF)(valWidthDgr, valAlf  - arrX[0]);
  double val_dF1_po_dTetta = (*pfncdF_po_dTet)(valWidthDgr, valAlf - arrX[0]);
  double val_d2F1_po_dTetta2 = (*pfncd2F_po_dTet2)(valWidthDgr, valAlf - arrX[0]);


  double val_F2 = (*pfncF)(valWidthDgr, valAlf - arrX[1]);
  double val_dF2_po_dTetta = (*pfncdF_po_dTet)(valWidthDgr, valAlf - arrX[1]);
  double val_d2F2_po_dTetta2 = (*pfncd2F_po_dTet2)(valWidthDgr, valAlf - arrX[1]);

  double modSqK1 = arrX[2]* arrX[2] + arrX[4]* arrX[4];
  double modSqK2 = arrX[3]* arrX[3] + arrX[5]* arrX[5];

	// ����������  arr_F
  double temp0 =  modSqK1 * val_F1 + arrX[2] *(arrX[3] *val_F2 - cmpSZv.m_Re)+ arrX[4] *(arrX[5] *val_F2 - cmpSZv.m_Im) ;
  arr_F[0] =- temp0  *  val_dF1_po_dTetta;
  double temp1 =  modSqK2 * val_F2 + arrX[3] *(arrX[2] *val_F1 - cmpSZv.m_Re)+ arrX[5] *(arrX[4] *val_F1 - cmpSZv.m_Im) ;
  arr_F[1] = -temp1 *  val_dF2_po_dTetta;
  //

  // ���������� dF
 arr_dF[0] = modSqK1 * val_dF1_po_dTetta * val_dF1_po_dTetta + temp0 * val_d2F1_po_dTetta2;
  // + ((arrX[2] * arrX[3] + arrX[4] * arrX[5] )*  val_F2 -arrX[2] * cmpSZv.m_Re -arrX[4]* cmpSZv.m_Im )* val_d2F1_po_dTetta2;

arr_dF[1] = (arrX[2] * arrX[3] + arrX[4] * arrX[5] ) *  val_dF1_po_dTetta *  val_dF2_po_dTetta;

arr_dF[2] =arr_dF[1] ;
arr_dF[3] = modSqK2 * val_dF2_po_dTetta * val_dF2_po_dTetta  + temp0 * val_d2F2_po_dTetta2;
//  + ((arrX[2] * arrX[3] + arrX[4] * arrX[5] )*  val_F1 -arrX[3] * cmpSZv.m_Re -arrX[5]* cmpSZv.m_Im )* val_d2F2_po_dTetta2;






  ///
 /*
 // 1. ���������� ������ ������ arr_dF
 arr_dF[0] = ( arrX[2]*arrX[2] +  arrX[4]* arrX[4]) * val_dF1_po_dTetta * val_dF1_po_dTetta
   +(arrX[2] * val_F1 + arrX[3]* val_F2 - cmpSZv.m_Re + arrX[4]* val_F1 + arrX[5]* val_F2 - cmpSZv.m_Im) *  val_d2F1_po_dTetta2 ;

 arr_dF[1] = (arrX[2] * arrX[3] + arrX[4]*arrX[5]) *val_dF1_po_dTetta *val_dF2_po_dTetta ;

 arr_dF[2] = -2.* arrX[2] *val_F1 * val_dF1_po_dTetta - ( arrX[3] * val_F2 - cmpSZv.m_Re  )*val_dF1_po_dTetta ;

 arr_dF[3] = - arrX[2] * val_F2 * val_dF1_po_dTetta ;

 arr_dF[4] = -2. * arrX[4] *val_F1 * val_dF1_po_dTetta  -( arrX[5]*val_F2 - cmpSZv.m_Im ) *  val_dF1_po_dTetta ;

 arr_dF[5] = -arrX[4]* val_F2 * val_dF1_po_dTetta ;
 ///
 // 2. ���������� ������  ������ arr_dF
 arr_dF[6] = arr_dF[1];

 arr_dF[7] = ( arrX[3]*arrX[3] +  arrX[5]* arrX[5]) * val_dF2_po_dTetta * val_dF2_po_dTetta
  + (valG1 *arrX[3] + valG2 *arrX[5]) *val_d2F2_po_dTetta2 ;

 arr_dF[8] =  -val_F1 * val_dF2_po_dTetta *arrX[3] ;


 arr_dF[9] = -2.  *val_F2 * val_dF2_po_dTetta * arrX[3]  -( arrX[2]*val_F1 - cmpSZv.m_Re ) *  val_dF2_po_dTetta ;



 arr_dF[10] =  -val_F1* val_dF2_po_dTetta * arrX[5] ;


 arr_dF[11] =  -2. * arrX[5] *val_F2 * val_dF2_po_dTetta  -( arrX[4] *val_F1 - cmpSZv.m_Im ) *  val_dF2_po_dTetta ;
 ///
 // 3. ���������� �������  ������ arr_dF
 arr_dF[12] =  -arrX[2]*val_dF1_po_dTetta *val_F1 -  valG1 * val_dF1_po_dTetta ;

 arr_dF[13] = -arrX[3] *val_F1 *val_dF2_po_dTetta ;

 arr_dF[14] =  val_F1* val_F1;

 arr_dF[15] = val_F1* val_F2;

 arr_dF[16] =   0.;

 arr_dF[17] = 0;
 ///
 // 4. ���������� 4-��   ������ arr_dF
 arr_dF[18] =  -arrX[2]* val_F2 * val_dF1_po_dTetta ;

 arr_dF[19] = -arrX[3]*val_F2 *val_dF2_po_dTetta - valG1 * val_dF2_po_dTetta;

 arr_dF[20] =  val_F1* val_F2;

 arr_dF[21] = val_F1* val_F1;

 arr_dF[22] =   0.;

 arr_dF[23] = 0;
 ///

 // 5. ���������� 5-��   ������ arr_dF
 arr_dF[24] =  -arrX[4] * val_F1 * val_dF1_po_dTetta - valG2 * val_dF1_po_dTetta;


 arr_dF[25] = -arrX[5] * val_F1* val_dF2_po_dTetta ;

 arr_dF[26] = 0.;


 arr_dF[27] =  0.;


 arr_dF[28] =    val_F1* val_F1;;

 arr_dF[29] = val_F1* val_F2;
 ///
// 6. ���������� 6-��   ������ arr_dF
 arr_dF[30] =  -arrX[4] * val_F2 * val_dF1_po_dTetta ;


 arr_dF[31] = -arrX[5] * val_dF2_po_dTetta* val_F2 - valG2 * val_dF2_po_dTetta;

 arr_dF[32] = 0.;


 arr_dF[33] =  0.;


 arr_dF[34] =   val_F1* val_F2;

 arr_dF[35] = val_F2* val_F2; //
   */
 }


 // ������ ������ ������� F(X)=dFGreek(X)/dX � ������� ����� dF(X)/dX ��� ���� ���������
 // INPUT:
 // valWidthDgr - ������ ��������� �� ������ 0707 � ���
 // cmpSZv  - ���������, ����������� �����
 // valAlf  - ���� ���������� ��������� �� �����������(���� �����)
 // fncF    - ������� ���������
// fncdF_po_dTet - ������� ������ ����������� ���������
// fncd2F_po_dTet2  - ������� ������  ����������� ���������
 // arrX[6] -  ������ ����������, � ������� ������������ ����������
 // arrX[0] - ���� ����� ����
 // arrX[1] -  ���� ����� ��������
 // arrX[2] -k11 �������������� �����  ������������ ��������� ����
 // arrX[3] -k21 �������������� �����  ������������ ��������� ���������
 // arrX[4] - k12 ������ �����  ������������ ��������� ����
 // arrX[5] - k22 ������ �����  ������������ ��������� ���������
 // OUTPUT:
 // arr_F[6]  - F(X)
 // arr_dF[36] -dF(X)/dX
 void calcPartial_F_and_dF_Newton6 (double valWidthDgr,TComp cmpSZv, double valAlf0,
 double(*pfncF)(double valWidthDgr, double tet)
 , double(*pfncdF_po_dTet)(double valWidthDgr, double tet), double(*pfncd2F_po_dTet2)(double valWidthDgr, double tet),
 double *arrX, double * arr_F, double * arr_dF )
 {
  double valAlf = valAlf0 ;

  memset( arr_F, 0, 6 * sizeof(double));
  memset( arr_dF, 0, 36 * sizeof(double));
  double val_F1 = (*pfncF)(valWidthDgr, valAlf  - arrX[0]);
  double val_dF1_po_dTetta = (*pfncdF_po_dTet)(valWidthDgr, valAlf - arrX[0]);
  double val_d2F1_po_dTetta2 = (*pfncd2F_po_dTet2)(valWidthDgr, valAlf - arrX[0]);

  double val_F2 = (*pfncF)(valWidthDgr, valAlf - arrX[1]);
  double val_dF2_po_dTetta = (*pfncdF_po_dTet)(valWidthDgr, valAlf - arrX[1]);
  double val_d2F2_po_dTetta2 = (*pfncd2F_po_dTet2)(valWidthDgr, valAlf - arrX[1]);

  double valG1 = arrX[2] * val_F1 +  arrX[3] * val_F2 - cmpSZv.m_Re ;
  double valG2 = arrX[4] * val_F1 +  arrX[5] * val_F2 - cmpSZv.m_Im ;

  // ����������  arr_F
  arr_F[0] =  - valG1 * arrX[2] *val_dF1_po_dTetta -  valG2 * arrX[4] * val_dF1_po_dTetta;
  arr_F[1] =  - valG1 * arrX[3] *val_dF2_po_dTetta -  valG2 * arrX[5] *val_dF2_po_dTetta ;
  arr_F[2] =   valG1 *val_F1;
  arr_F[3] =   valG1 *val_F2;
  arr_F[4] =   valG2 *val_F1;
  arr_F[5] =   valG2 *val_F2 ;
  ///

 // 1. ���������� ������ ������ arr_dF
 arr_dF[0] = ( arrX[2]*arrX[2] +  arrX[4]* arrX[4]) * val_dF1_po_dTetta * val_dF1_po_dTetta
   +(arrX[2] * val_F1 + arrX[3]* val_F2 - cmpSZv.m_Re + arrX[4]* val_F1 + arrX[5]* val_F2 - cmpSZv.m_Im) *  val_d2F1_po_dTetta2 ;

 arr_dF[1] = (arrX[2] * arrX[3] + arrX[4]*arrX[5]) *val_dF1_po_dTetta *val_dF2_po_dTetta ;

 arr_dF[2] = -2.* arrX[2] *val_F1 * val_dF1_po_dTetta - ( arrX[3] * val_F2 - cmpSZv.m_Re  )*val_dF1_po_dTetta ;

 arr_dF[3] = - arrX[2] * val_F2 * val_dF1_po_dTetta ;

 arr_dF[4] = -2. * arrX[4] *val_F1 * val_dF1_po_dTetta  -( arrX[5]*val_F2 - cmpSZv.m_Im ) *  val_dF1_po_dTetta ;

 arr_dF[5] = -arrX[4]* val_F2 * val_dF1_po_dTetta ;
 ///
 // 2. ���������� ������  ������ arr_dF
 arr_dF[6] = arr_dF[1];

 arr_dF[7] = ( arrX[3]*arrX[3] +  arrX[5]* arrX[5]) * val_dF2_po_dTetta * val_dF2_po_dTetta
  + (valG1 *arrX[3] + valG2 *arrX[5]) *val_d2F2_po_dTetta2 ;

 arr_dF[8] =  -val_F1 * val_dF2_po_dTetta *arrX[3] ;


 arr_dF[9] = -2.  *val_F2 * val_dF2_po_dTetta * arrX[3]  -( arrX[2]*val_F1 - cmpSZv.m_Re ) *  val_dF2_po_dTetta ;



 arr_dF[10] =  -val_F1* val_dF2_po_dTetta * arrX[5] ;


 arr_dF[11] =  -2. * arrX[5] *val_F2 * val_dF2_po_dTetta  -( arrX[4] *val_F1 - cmpSZv.m_Im ) *  val_dF2_po_dTetta ;
 ///
 // 3. ���������� �������  ������ arr_dF
 arr_dF[12] =  -arrX[2]*val_dF1_po_dTetta *val_F1 -  valG1 * val_dF1_po_dTetta ;

 arr_dF[13] = -arrX[3] *val_F1 *val_dF2_po_dTetta ;

 arr_dF[14] =  val_F1* val_F1;

 arr_dF[15] = val_F1* val_F2;

 arr_dF[16] =   0.;

 arr_dF[17] = 0;
 ///
 // 4. ���������� 4-��   ������ arr_dF
 arr_dF[18] =  -arrX[2]* val_F2 * val_dF1_po_dTetta ;

 arr_dF[19] = -arrX[3]*val_F2 *val_dF2_po_dTetta - valG1 * val_dF2_po_dTetta;

 arr_dF[20] =  val_F1* val_F2;

 arr_dF[21] = val_F1* val_F1;

 arr_dF[22] =   0.;

 arr_dF[23] = 0;
 ///

 // 5. ���������� 5-��   ������ arr_dF
 arr_dF[24] =  -arrX[4] * val_F1 * val_dF1_po_dTetta - valG2 * val_dF1_po_dTetta;


 arr_dF[25] = -arrX[5] * val_F1* val_dF2_po_dTetta ;

 arr_dF[26] = 0.;


 arr_dF[27] =  0.;


 arr_dF[28] =    val_F1* val_F1;;

 arr_dF[29] = val_F1* val_F2;
 ///
// 6. ���������� 6-��   ������ arr_dF
 arr_dF[30] =  -arrX[4] * val_F2 * val_dF1_po_dTetta ;


 arr_dF[31] = -arrX[5] * val_dF2_po_dTetta* val_F2 - valG2 * val_dF2_po_dTetta;

 arr_dF[32] = 0.;


 arr_dF[33] =  0.;


 arr_dF[34] =   val_F1* val_F2;

 arr_dF[35] = val_F2* val_F2; //

 }

 // ������� �������� ��������� ������������ ����� ���� � �������� ������� �������
// ��� ������������� �����
// INPUT:
// quantDiagr  - �-�� ��������
// valWidthDgr  - ������ ��������� �� ������ 0707 � ���
// pcmpSZv   - ������ ��������� �� ���������� (����������� �����)
// pAlf     - ������ ����� ���������� ��������� �� ������������
// fncF    - ������� ���������
// cmpZ1  - ����������� ����������� ��������� ����
// cmpZ2  - ����������� ����������� ��������� ��������
// alfTrg - ���� ����� ����
// alfAnp - ���� ����� ��������
double calcFncObj_Newton6(const int quantDiagr,double valWidthDgr,TComp *pcmpSZv, double *pAlf,
 double(*pfncF)(double valWidthDgr, double tet)
, TComp cmpZ1, TComp cmpZ2, double alfTrg, double alfAnp )
{

 double arrX[6] ={0.};

 arrX[0] =  alfTrg ;
 arrX[1] = alfAnp ;
 arrX[2] = cmpZ1.m_Re  ;   // ����� ��������� ����
 arrX[3] = cmpZ2.m_Re  ;   // ����� ��������� ����
 arrX[4] = cmpZ1.m_Im  ; // ����� ��������� ��������
 arrX[5] = cmpZ2.m_Im  ; // ����� ��������� ��������
 ///

 double valRez = 0.;
 for (int i =0; i < quantDiagr; i++)
 {
   double valPartTemp =  calcFncPartialObj_Newton6( valWidthDgr, pcmpSZv[i], pAlf[i],
 pfncF, arrX ) ;
 valRez +=  valPartTemp ;
 }
return valRez;
}


double calcFncPartialObj_Newton6(double valWidthDgr,TComp cmpSZv, double valAlf0,
 double(*pfncF)(double valWidthDgr, double tet), double *arrX )
{
  double valAlf = valAlf0 ;
  double val_F1 = (*pfncF)(valWidthDgr, valAlf  - arrX[0]);
  double val_F2 = (*pfncF)(valWidthDgr, valAlf - arrX[1]);
  double valG1 = arrX[2] * val_F1 +  arrX[3] * val_F2 - cmpSZv.m_Re ;
  double valG2 = arrX[4] * val_F1 +  arrX[5] * val_F2 - cmpSZv.m_Im ;

  return valG1* valG1 + valG2*valG2;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////  ����� ������ ������� ��� ������������ ������ � 6 ������������////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////  ����� ���������� ��������������� ������ ��� ������������ ������ � 6 ������������///////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

// ������� �������� ��������� ������������ ����� ���� � �������� ������� �������
// ��� ������������� �����
// INPUT:
// quantDiagr  - �-�� ��������
// valWidthDgr  - ������ ��������� �� ������ 0707 � ���
// pcmpSZv   - ������ ��������� �� ���������� (����������� �����)
// pAlf     - ���� ���������� ��������� �� ������������
// fncF    - ������� ���������
// fncdF_po_dTet - ������� ������ ����������� ���������
// fncd2F_po_dTet2  - ������� ������  ����������� ���������
// OUTPUT:
// pcmpZ1  - ����������� ����������� ��������� ����
// pcmpZ2  - ����������� ����������� ��������� ��������
// palfTrg - ���� ����� ����
// palfAnp - ���� ����� ��������
// ����������:
//  -1 - ���� �-�� �������� quantDiagr< 3
// -2 - ����� ������� ����� ���������
// -3 -���� ������������� ������� ������ ������� ���������������� =0
int find_Targ_and_Antp_GroupRelax6(const int quantDiagr,double valWidthDgr,TComp *pcmpSZv, double *pAlf,
 double(*pfncF)(double valWidthDgr, double tet)
 , double(*pfncdF_po_dTet)(double valWidthDgr, double tet), double(*pfncd2F_po_dTet2)(double valWidthDgr, double tet),
 TComp *pcmpZ1, TComp *pcmpZ2, double *palfTrg, double *palfAnp )
{
	if (quantDiagr < 3)
	{
	return -1;
	}
	// ��������� �����������:
	double arrX0[6] ={0.}, arrX[6] ={0.};

	arrX[0] =  *palfTrg ;
	arrX[1] = *palfAnp ;
	arrX[2] = (*pcmpZ1).m_Re  ;   // ����� ��������� ���� Re
	arrX[3] = (*pcmpZ2).m_Re  ;   // ����� ��������� ��������  Im
	arrX[4] = (*pcmpZ1).m_Im  ; // ����� ��������� ���� Re
	arrX[5] = (*pcmpZ2).m_Im  ; // ����� ��������� �������� Im
	memcpy(arrX0, arrX, 6 * sizeof(double));
	///
	double arrDelX[6] ={0.};
	int i = 0;

	int ireturn = -2;
	double del = -2.;
	TComp pcmpZ[2];
	pcmpZ[0] = *pcmpZ1;
	pcmpZ[1] = *pcmpZ2;
	double  valObj0,valObj1,valObj00;
	for (i = 0; i < 1000; i++)
	{
	// ���������� ������������� ��������� ��� ������������� ����� ����� �������
	// 2 ����� ������� ���������������� ��� ����������� 2-�� �������

		calcFncObj_Newton6(quantDiagr,valWidthDgr,pcmpSZv, pAlf,
		pfncF, pcmpZ[0], pcmpZ[1], arrX[0] , arrX[1] ) ;

		if(! findOtimal_Z1_and_Z2( quantDiagr,valWidthDgr,pcmpSZv, pAlf,
		pfncF, pcmpZ,  arrX ) ) return -3;

		calcFncObj_Newton6(quantDiagr,valWidthDgr,pcmpSZv, pAlf,
		pfncF, pcmpZ[0], pcmpZ[1], arrX[0] , arrX[1] ) ;
		arrX[2] = pcmpZ[0].m_Re  ;   // ����� ��������� ���� Re
		arrX[3] = pcmpZ[1].m_Re  ;   // ����� ��������� ��������  Im
		arrX[4] = pcmpZ[0].m_Im  ; // ����� ��������� ���� Re
		arrX[5] = pcmpZ[1].m_Im  ; // ����� ��������� �������� Im

		///
		int irez= Relax_Angles_GroupRelax6( quantDiagr,valWidthDgr,pcmpSZv, pAlf,
		pfncF, pfncdF_po_dTet, pfncd2F_po_dTet2, &arrX[2], arrX );

		calcFncObj_Newton6(quantDiagr,valWidthDgr,pcmpSZv, pAlf,
		pfncF, pcmpZ[0], pcmpZ[1], arrX[0] , arrX[1] ) ;

		if (irez !=0)
		{
		return irez ;
		}

		MtrxMinusMatrx(arrX, arrX0,6, 1, arrDelX);
		del = NormVect(arrDelX, 2);
		memcpy(arrX0 ,arrX, 6 * sizeof(double));
		if (del< 0.01)
		{
		ireturn = 0;


		break ;
		}

  }
	if (arrX[0] > arrX[1])
	{
	*pcmpZ1 = pcmpZ[0];

	*pcmpZ2 = pcmpZ[1];

	*palfTrg =  arrX[0];
	*palfAnp =  arrX[1];
	}
	else
	{
	*pcmpZ1 = pcmpZ[1];

	*pcmpZ2 = pcmpZ[0];

	*palfTrg =  arrX[1];
	*palfAnp =  arrX[0];

	}
	return ireturn;
}

// ���������� ����������� ������������ ���������  ���� � �������� ��� ������������� �����
// ��������� ���� � ��������
// ��� ������������� �����
// INPUT:
// quantDiagr  - �-�� ��������
// valWidthDgr  - ������ ��������� �� ������ 0707 � ���
// pcmpSZv   - ������ ��������� �� ���������� (����������� �����)
// pAlf     - ���� ���������� ��������� �� ������������
// fncF    - ������� ���������
// fncdF_po_dTet - ������� ������ ����������� ���������
// fncd2F_po_dTet2  - ������� ������  ����������� ���������
// arrMilradAngles[0] - ���� ����� ����  � ������������
// arrMilradAngles[1]  - ���� ����� �������� � ������������

// OUTPUT:
// pcmpSZv [2] - ������ ����������� ����� ������� ���������
// ����������:
//  -1 - ���� �-�� �������� quantDiagr< 3
// -2 - ����� ������� ����� ���������
// -3 - ���� ����� �� �������
bool findOtimal_Z1_and_Z2(const int quantDiagr,double valWidthDgr,TComp *pcmpSZv, double *pAlf,
 double(*pfncF)(double valWidthDgr, double tet), TComp *pcmpZ, double *arrMilradAngles )
{
// ���������� ������������� ��������� ��� ������������� ����� ����� �������
	// 2 ����� ������� ���������������� ��� ����������� 2-�� �������
	double arrZ[4] ={0.};
	double arrFSum[4] ={0.}, arrB1[2]={0.}, arrB2[2]={0.}, arrF[2]={0.}, arrT0[2] ={0.}, arrT1[4]={0.}, arrT2[4]={0.};
	for (int i =0; i < quantDiagr; i++)
	{
	 arrF[0] =  pfncF( valWidthDgr, pAlf[i] - arrMilradAngles[0]) ;
	 arrF[1] =  pfncF( valWidthDgr, pAlf[i]  - arrMilradAngles[1]) ;
	 arrB1[0] +=  pcmpSZv[i].m_Re * arrF[0];
	 arrB1[1] +=  pcmpSZv[i].m_Re * arrF[1];
	 arrB2[0] +=  pcmpSZv[i].m_Im * arrF[0];
	 arrB2[1] +=  pcmpSZv[i].m_Im * arrF[1];
	 memcpy(arrT0,arrF,  2 * sizeof(double));
	 MtrxMultMatrxTransp(arrF, 2, 1, arrT0,2, arrT1) ;
	 MtrxSumMatrx(arrFSum, arrT1,2, 2, arrT2) ;
	 memcpy(arrFSum,arrT2, 4 * sizeof(double));
	}
   if( !(SolvLinEq2(arrFSum, arrB1,arrZ)&&SolvLinEq2(arrFSum, arrB2, &arrZ[2]))) return false;

	pcmpZ[0].m_Re = arrZ[0] ;
	pcmpZ[0].m_Im = arrZ[2] ;
	pcmpZ[1].m_Re = arrZ[1] ;
	pcmpZ[1].m_Im = arrZ[3] ;
	return true;
}
// ���������� ����������� ����� ���� � �������� ��� ������������� ������� ���������
// ��� ������������� �����
// INPUT:
// quantDiagr  - �-�� ��������
// valWidthDgr  - ������ ��������� �� ������ 0707 � ���
// pcmpSZv   - ������ ��������� �� ���������� (����������� �����)
// pAlf     - ���� ���������� ��������� �� ������������
// fncF    - ������� ���������
// fncdF_po_dTet - ������� ������ ����������� ���������
// fncd2F_po_dTet2  - ������� ������  ����������� ���������
// arrZ[4] - ������ � �������������� ��������� ���� � ��������  ZTarg = k11 + i *k12, ZAntp = k21 + i *k22,
//    arrZ[0] = k11
//    arrZ[2] = k21
//    arrZ[0] = k12
//    arrZ[0] = k22
// arrMilradAngles[0] - ���� ����� ����  � ������������
// arrMilradAngles[1]  - ���� ����� �������� � ������������

// OUTPUT:
// arrMilradAngles[0] - ���� ����� ����   � ������������
// arrMilradAngles[1]  - ���� ����� �������� � ������������
// ����������:
//  -1 - ���� �-�� �������� quantDiagr< 3
// -2 - ����� ������� ����� ���������
// -3 - ���� ����� �� �������
int Relax_Angles_GroupRelax6(const int quantDiagr,double valWidthDgr,TComp *pcmpSZv, double *pAlf,
 double(*pfncF)(double valWidthDgr, double tet)
 , double(*pfncdF_po_dTet)(double valWidthDgr, double tet), double(*pfncd2F_po_dTet2)(double valWidthDgr, double tet),
 double *arrZ, double *arrAngles )
{
  if (quantDiagr < 3)
  {
   return -1;
  }
  // ��������� �����������:
  double arrX[6] ={0.};

 arrX[0] =  arrAngles[0];
 arrX[1] = arrAngles[1];
  arrX[2] = arrZ[0]  ;   // ����� ��������� ����
  arrX[3] = arrZ[1]  ;   // ����� ��������� ��������
  arrX[4] = arrZ[2]  ; // ����� ���������  ����
  arrX[5] = arrZ[3]  ; // ����� ��������� ��������
  ///
  double arrDelX[6] ={0.};
  int i = 0;
  double arr_FGreek[6] = {0.}, arr_dFGreek[36] ={0.};
  double arr_FGreekDim2[2] ={0.},arr_dFGreekDim2[4] ={0.},arr_dFGreekDim2Inv[4] ={0.} ;

  double del = -2.;
  for (i = 0; i < 100; i++)
  {
	calcF_and_dF_Newton6( quantDiagr, valWidthDgr,pcmpSZv, pAlf,
	   pfncF,pfncdF_po_dTet, pfncd2F_po_dTet2, arrX, arr_FGreek, arr_dFGreek ) ;
	arr_FGreekDim2[0] = arr_FGreek[0];
	arr_FGreekDim2[1] = arr_FGreek[1];
	arr_dFGreekDim2[0] = arr_dFGreek[0];
	arr_dFGreekDim2[1] = arr_dFGreek[1];
	arr_dFGreekDim2[2] = arr_dFGreekDim2[1];
	arr_dFGreekDim2[3] = arr_dFGreek[6 + 1];

	// ���������!!!
	double arr_F00[2] ={0.},arr_dF00[4] ={0.} ;
	calcF_and_dF_AnglesRelax_Newton6 ( quantDiagr, valWidthDgr,pcmpSZv, pAlf,
	   pfncF,pfncdF_po_dTet, pfncd2F_po_dTet2, arrX, arr_F00, arr_dF00 ) ;
	arr_dFGreekDim2[0] = arr_dF00 [0];
	arr_dFGreekDim2[1] = arr_dF00 [1];
	arr_dFGreekDim2[2] = arr_dFGreekDim2[1];
	arr_dFGreekDim2[3] = arr_dF00 [3];

	///

	if(!InverseMtrx2(arr_dFGreekDim2, arr_dFGreekDim2Inv)) return -2;
	MtrxMultMatrx(arr_dFGreekDim2Inv ,2, 2, arr_FGreekDim2,1, arrDelX) ;
	del = NormVect(arrDelX, 2);
	double arrT[2] ={0.};
	MtrxMinusMatrx(arrX, arrDelX,1, 2, arrT) ;
	memcpy( arrX, arrT, 2 * sizeof(double));
	if (del< 0.00001)
	{

	  arrAngles[0] = arrX[0]   ;
	  arrAngles[1] = arrX[1]  ;
	  return 0;
	}
 }
  return -3;
}

double SIGN__(double a)
{
	return (a > 0.)?1:-1;
}
#pragma package(smart_init)
