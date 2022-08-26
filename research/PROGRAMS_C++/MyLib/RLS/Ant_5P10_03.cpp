//---------------------------------------------------------------------------


#pragma hdrstop

#include "Ant_5P10_03.h"
#include <string.h>
#include <math.h>
#include <float.h>
#include "Comp.h"
#include "DiagrSinX.h"

extern const double TET0707;
extern const double VAL_C ;


TAnt_5P10_03::TAnt_5P10_03()
{
 // �-�� ����������� ��������� � �����
  mQuantDiagr = 7;
 // �������� �������
  mAppert = 120.;
 // ����� �����
   mLambda = 3. ;
 // ���� ������ ����������� ��������  �������
  mAngSdvig = 1.5026 *M_PI/ 180.;
 // ��������z ���� � ����������� ����������
  mNoiseDisp = 400;

 // ������� ������� �������� k � ����������� ���������
  mAmplFactSig = 0.01;
}
// ����������� �����������
TAnt_5P10_03::TAnt_5P10_03 (const TAnt_5P10_03 &R2)
 {
// �-�� ����������� ��������� � �����
  mQuantDiagr = R2.mQuantDiagr;
 // �������� �������
  mAppert = R2.mAppert;
 // ����� �����
   mLambda = R2. mLambda ;
 // ���� ������ ����������� ��������  �������
  mAngSdvig = R2.mAngSdvig;
 // ��������z ���� � ����������� ����������
  mNoiseDisp = R2.mNoiseDisp;

 // ������� ������� �������� k � ����������� ���������
  mAmplFactSig = R2. mAmplFactSig;
 }

  // �������� ������������
  TAnt_5P10_03 TAnt_5P10_03::operator=(TAnt_5P10_03  R2)
{
// �-�� ����������� ��������� � �����
  mQuantDiagr = R2.mQuantDiagr;
 // �������� �������
  mAppert = R2.mAppert;
 // ����� �����
   mLambda = R2. mLambda ;
 // ���� ������ ����������� ��������  �������
  mAngSdvig = R2.mAngSdvig;
 // ��������z ���� � ����������� ����������
  mNoiseDisp = R2.mNoiseDisp;

 // ������� ������� �������� k � ����������� ���������
  mAmplFactSig = R2.mAmplFactSig;

  return *this ;
}



// ����� ������  1
 TAnt_5P10_03::TAnt_5P10_03(const int quantDiagr,const double Appert,const double Lambda
   ,const double AngSdvig, const double NoiseDisp, const double AmplFactSig)
 {
	 mAppert = Appert ;
	 mLambda = Lambda ;
	 mNoiseDisp = NoiseDisp;
	 mAngSdvig = AngSdvig;
	 mAmplFactSig = AmplFactSig;
	 mQuantDiagr = quantDiagr;
 }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// ������� �����/////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// �������� ���� � ��� � ���������� ����
double TAnt_5P10_03::transformAngToGeneralizedAng__ (const double  valAng  )
{
	return  mAppert * M_PI/ mLambda  * sin (valAng);
}

// �������� ����������� ���� � ������� � ���
double TAnt_5P10_03::transformGeneralizedAngToAng__ (const double  GeneralizedAng  )
{
	return GeneralizedAng  * mLambda/ mAppert/ M_PI;//asin ( GeneralizedAng  * mLambda/ mAppert/ M_PI);
}




// ���������� ������� b
// INPUT:
// cmparrS[3] - ������ ��������� ������ ��������
// OUTPUT:
// arr_b[2] -  ������������ ���������
bool TAnt_5P10_03::calcVect_b__(TComp *cmparrS, double *arr_b)
{
 // ����������� ������� b
   double arrS[4] = {0.};
   arrS[0] = cmparrS[0].m_Re;
   arrS[1] = cmparrS[2].m_Re;
   arrS[2] = cmparrS[0].m_Im;
   arrS[3] = cmparrS[2].m_Im;
   double arrs[2] = {0.};
   arrs[0] = cmparrS[1].m_Re;
   arrs[1] = cmparrS[1].m_Im;

   double  arrSInv[4] = {0.};
   bool brez =   InverseMtrx2__(arrS, arrSInv);
   MtrxMultMatrx__(arrSInv,2, 2, arrs,1, arr_b) ;
   return brez;
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
void TAnt_5P10_03::createMtrxCorrForMeasuresForAnsambleOfPartDiagrs__(double valGenTargEps,  double valGenAntpEps
  ,TComp cmpKTarg ,TComp cmpKAntp, int iNumAnsamble, int lenAnsamble
  , double *arrMtrxCorr )
 {
	const double  valGenAngSdvig =  transformAngToGeneralizedAng__ (mAngSdvig );
	memset(arrMtrxCorr, 0, 4 * lenAnsamble * lenAnsamble * sizeof(double));

  // ���������� ������ �������, ������������� ��������� ������� �������� ��������
	for (int i = 0; i < lenAnsamble; i++)  // ���� �� ���������� � ��������
	{
	   // ���������� ���� ����� ������������ ��� ��������� � ������� i � ���������� �����������
	  double valDiagrGenAng = calGenAngSdvigaPartDiagr__(iNumAnsamble + i); // ���� ��� ����� � ������� i
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
						   * mAmplFactSig * mAmplFactSig;
	   // ��������� ������� ����� �������
	   double valDispIm = ( (valDiagrTarg * cmpKTarg.m_Im * valDiagrTarg * cmpKTarg.m_Im)
						  + (valDiagrAntp * cmpKAntp.m_Im * valDiagrAntp * cmpKAntp.m_Im) )
						   * mAmplFactSig * mAmplFactSig;
	 ///
	  // ���������� ������������ ��������� ������� arrMtrxCorr � ��������������� �������
	  // �������������� � �������5���� ����� ������� ��������� � ������� i �������� ������ 2*i � 2*i+1 ��������������
	  arrMtrxCorr[2*lenAnsamble * (2 * i) + 2 *i] =  valDispRe  + mNoiseDisp;
	  arrMtrxCorr[2*lenAnsamble * (2 * i + 1) + 2 *i + 1] =  valDispIm  + mNoiseDisp;

	}
	return;
 }



// ���������� ������ ������� ������ ���������� ���������� ����� ���� � ��������
// �������������� �������� ����� ��������� ��������
// INPUT:
// valGenTargEps, valGenAntpEps - ������ ���������� ����� ���� � �������� ������������ ����������� ��������� ������
// OUTPUT:
// arrMtrxCorrGenAngs[4] - ������ ������� ������ ���������� ����� ���� � ���������
bool TAnt_5P10_03::calcMtrxCorrGenAngs_AlternativeMeth3PartDiagr__(TComp *cmparrS,  double valGenTargEps,  double valGenAntpEps
  ,TComp cmpKTarg ,TComp cmpKAntp ,const int NumAnsamble, double *arrMtrxCorrGenAngs )
{
	 double valGenAngSdvig = transformAngToGeneralizedAng__ (mAngSdvig ) ;

	 // ��������
	 double arr_b[2] ={0.};
	 calcVect_b__(&cmparrS[NumAnsamble], arr_b);
	 double val_b0 = arr_b[0];
	 double val_b1 = arr_b[1];

  // ������������ �������� ������� ����������� dB/dS
  // ������ ������ ��� �������� ������� b0 �� S
  // ������ ������ ��� �������� ������� b1 �� S
  double arrdB_po_dS[12] = {0.};
  calc_dB0_po_dS__(&cmparrS[NumAnsamble], valGenAngSdvig,  val_b0,  val_b1, arrdB_po_dS);
  calc_dB1_po_dS__(&cmparrS[NumAnsamble], valGenAngSdvig,  val_b0,  val_b1, &arrdB_po_dS[6]);
  ///

  double arrGenAngs[2] = {0.};
  arrGenAngs[0] =  valGenTargEps;
  arrGenAngs[1] =  valGenAntpEps;
  ///

  // ������������ ������� ������� ����������� ����� �� ����������
  double arr_dMu_po_dS[12] = {0.};
  for (int i =0; i < 2; i++)
  {

	 double coeff = fncDerivDiagrSinx_div_x(arrGenAngs[i])
				   -  val_b0 *  fncDerivDiagrSinx_div_x(arrGenAngs[i] + valGenAngSdvig)
				   -  val_b1 *  fncDerivDiagrSinx_div_x(arrGenAngs[i] - valGenAngSdvig) ;
	 double arrF[2] = {0.};
	 arrF[0] = -fncDiagrSinx_div_x(arrGenAngs[i] + valGenAngSdvig) / coeff;
	 arrF[1] = -fncDiagrSinx_div_x(arrGenAngs[i] - valGenAngSdvig) / coeff;
	// MtrxTranspMultMatrx(arrdB_po_dS,2, 6, arrF, 1, &arr_dMu_po_dS[ i * 6]) ;
	 MtrxMultMatrx__(arrF,1, 2, arrdB_po_dS, 6, &arr_dMu_po_dS[ i * 6]) ;
  }

  ///

  // ������������ �������� ������� ���������
  double arrCorrPartDiagr[36] = {0.};

  createMtrxCorrForMeasuresForAnsambleOfPartDiagrs__( valGenTargEps,  valGenAntpEps
  , cmpKTarg , cmpKAntp ,NumAnsamble, 3
  , arrCorrPartDiagr );
  ///

  // ���������� ������ ������� ������� �������
  double  arrtemp0[12] = {0.};
  MtrxMultMatrx__(arr_dMu_po_dS,2, 6, arrCorrPartDiagr, 6, arrtemp0) ;
  MtrxMultMatrxTransp__(arrtemp0,2, 6, arr_dMu_po_dS,2, arrMtrxCorrGenAngs) ;
  ///
  return true;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//////  �������������� ������� ��� ������ 3 ����������� ��������  /////////////////////////////////////////////////////////////////////////
////////  ������� ��� ������������� ���������� ///////////////////////////////////////////////////////////////////////
////////// ������� ���� /////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// ���������� ���������� ����� ���� � �������� �� ���������� 3 ����������� ��������
//INPUT:
// cmparrS[7] -  ������ �������
// NumRayTriple-  ����� ������ ������� ����� (��������� ����� ���� ����� ������� � ����)
//                ����� ������ ��� ���������� ����� ������� ����
//  OUTPUT:
//  *valEstAngTarg -   ���� ����, ���
// *valEstAngAntp -    ���� ��������, ���
// cmpKTarg - ���� ��������� ����
// cmpKAntp - ���� ��������� ��������
// arrMtrxCorr  - ������ ������� ������ ��������� ���� ���� � ��������
int TAnt_5P10_03::estimateMethThreePartDiagr__( TComp *cmparrS
  , int iNumRayTriple , double *valEstAngTarg, double *valEstAngAntp, TComp *cmpKTarg , TComp *cmpKAntp , double *arrMtrxCorr
  , double *pval_b0, double *pval_b1 )
{
	double valGenAngSdvig = transformAngToGeneralizedAng__ (mAngSdvig ) ;
	// ���������� ���������� ����� ���� � ��������
	double  valGenTargEps = -100.,  valGenAntpEps = -100.0;  // ����� ���� ���� � ��������
	if (EstGenAngsThreePartDiagr__(  &cmparrS[iNumRayTriple]
	,  &valGenTargEps,  &valGenAntpEps, pval_b0, pval_b1) !=2)
	{
	   return -1;
	}



	//
	// ���������� ������� ���������
	TComp cmparrA[4],cmparrAInv[4], cmparrK[2], cmparrSTemp[2];

	cmparrA[0] = TComp(fncDiagrSinx_div_x(valGenTargEps + valGenAngSdvig), 0.);
	cmparrA[1] = TComp(fncDiagrSinx_div_x (valGenAntpEps + valGenAngSdvig ), 0.);
	cmparrA[2] = TComp(fncDiagrSinx_div_x (valGenTargEps - valGenAngSdvig ), 0.);
	cmparrA[3] = TComp(fncDiagrSinx_div_x (valGenAntpEps - valGenAngSdvig ), 0.);

	bool brez =   InverseMtrx2( cmparrA, cmparrAInv);
	if (!brez)
	{
	 return -2;
	}

	cmparrSTemp[0] = cmparrS [iNumRayTriple ];
	cmparrSTemp[1] = cmparrS [iNumRayTriple + 2 ];
	MtrxMultMatrx(cmparrAInv,2, 2, cmparrSTemp,1, cmparrK);
	*cmpKTarg = cmparrK[0] ;
	*cmpKAntp = cmparrK[1] ;

	///

	// ���������� ����� ���� � �������� � ���
		  // �������������� ���� ����������� ���� ������ � ��������
			double valRadCentreRayPos = calcRadAngSdvigaPartDiagr__(iNumRayTriple + 1);
			// ( - ((double)mQuantDiagr -1.) /2. + 1. + ((double)iNumRayTriple))				*mAngSdvig ;
		  // ���������� ���� ���� � �������� � ���
			const double ValRadTargEps =  transformGeneralizedAngToAng__ ( valGenTargEps ) ;// � �������
			*valEstAngTarg = ValRadTargEps + valRadCentreRayPos;
		  // ���������� ���� �������� � �������� � ����
			const double ValRadAntp�Eps =  transformGeneralizedAngToAng__ (valGenAntpEps ) ;// � �������
			*valEstAngAntp = ValRadAntp�Eps + valRadCentreRayPos;
	///

	// ���������� ������ ������� ������ ����������� ���������� ����� ���� � ��������
	double arrMtrxCorrMu[4] = {0.};
   //	calcMtrxCorrGenAngs_Meth3PartDiagr(iNumRayTriple ,valGenTargEps,  valGenAntpEps
	 // ,*cmpKTarg,*cmpKAntp, arrMtrxCorrMu) ;

   calcMtrxCorrGenAngs_AlternativeMeth3PartDiagr__(cmparrS,   valGenTargEps,  valGenAntpEps
  , *cmpKTarg , *cmpKAntp ,iNumRayTriple, arrMtrxCorrMu );

	///

	// ���������� �������� ������� ������ ����������� ����� ���� � �������� � ���������
	double arrQ[4] = {0.}, arrTemp[4] = {0.};
	double coef = mLambda / mAppert / M_PI;
	arrQ[0] =  coef / sqrt( 1.- coef * valGenTargEps * coef * valGenTargEps);
	arrQ[3] =  coef / sqrt( 1.- coef * valGenAntpEps * coef * valGenAntpEps);
	MtrxMultMatrx__(arrQ,2, 2, arrMtrxCorrMu,2, arrTemp)  ;
	MtrxMultMatrxTransp__(arrTemp,2, 2, arrQ,2, arrMtrxCorr) ;
	///
	return 0;
}


// ���������� ���������� ����� ���� � �������� �� ���������� 3 ����������� ��������
//INPUT:
// wchFoldName1 - ���� � ������ � ���������
// cmparrS[3] -  ������ �������
//  OUTPUT:
//  *pvalGenTargEps -   ���� ����, ����� �����
// *pvalGenAntpEps -    ���� ��������, ����� �����
int TAnt_5P10_03::EstGenAngsThreePartDiagr__(TComp *cmparrS ,  double  *pvalGenTargEps
   ,  double  *pvalGenAntpEps, double *pval_b0, double *pval_b1)
{

  // ����������� ������� b
   double arr_b[2] = {0.}  ;
   bool brez = calcVect_b__(cmparrS,arr_b) ;
   ///
   *pval_b0 = arr_b[0];
   *pval_b1 = arr_b[1];


   double arrRoots[2] ={0.};
   int iNUmRoots = findRootsFgr_For_3PartDiagr__(arr_b,  arrRoots);
   *pvalGenTargEps =  arrRoots[1];
   *pvalGenAntpEps =  arrRoots[0];
   return iNUmRoots ;
}


// ���������� ������ ������� FGr ���� ������ 3 ��������� ��������
// ���������� �-�� ������
 int TAnt_5P10_03::findRootsFgr_For_3PartDiagr__(double *arr_b,  double *arrRoots)
 {
	double valGenAngSdvig = transformAngToGeneralizedAng__ (mAngSdvig ) ;
   int iNUmRoots = 0;

   double valMu0 =  M_PI + valGenAngSdvig - 0.3;
   double valStep = 0.25;
	const int lenBuff = 2. *valMu0/ valStep ;
   double *parrBuff  = new double [lenBuff];

  double valmu = -valMu0;

  for (int i=0 ; i < lenBuff; i++)
  {
   valmu = -valMu0 + ((double)i) * valStep;
   parrBuff[ i ] = fncFGr_PartDiagr__(arr_b, valmu);

  }

  for (int i = 0; i < (lenBuff-1); i++)
  {
   if (parrBuff[ i ] * parrBuff[ i + 1] <= 0.)
   {
	 if (iNUmRoots >1)
	 {
	   break;
	 }
	 double valX0 = -valMu0 + ((double)i) * valStep;
	 double valX1 = -valMu0 + ((double) (i +1)) * valStep;
	 arrRoots [iNUmRoots ] = findRootMethChord_For_Fgr_For_3PartDiagr__(  arr_b,  valX0, valX1) ;
	 iNUmRoots++;


   }
  }


  delete parrBuff;
  return iNUmRoots;

 }


double TAnt_5P10_03::fncFGr_PartDiagr__(double *arr_b, double valmu)
{
  double valGenAngSdvig = transformAngToGeneralizedAng__ (mAngSdvig ) ;
  return (fncDiagrSinx_div_x(valmu) - arr_b[0] * fncDiagrSinx_div_x(valmu +valGenAngSdvig)
		 - arr_b[1] * fncDiagrSinx_div_x( valmu - valGenAngSdvig));
}


double TAnt_5P10_03::findRootMethChord_For_Fgr_For_3PartDiagr__(double *arr_b, double  valX0, double valX1)
{
  double valMu0 = valX0 - 0.00001;
  double valMu1 = valX1 + 0.00001;
  double valMuRez =  valMu0 ;
  double eps = 0.001;
  for (int i = 0; i < 100; i++)
  {
	double f0 = fncFGr_PartDiagr__(arr_b, valMu0 ) ;
	double f1 = fncFGr_PartDiagr__(arr_b, valMu1 ) ;
	double valMuRez0 = valMu0 - f0 * ( valMu1 - valMu0)/ (f1 - f0);
	if (fabs(valMuRez0 - valMuRez) < eps)
	{
	 return valMuRez0;
	}
	double f2 = fncFGr_PartDiagr__(arr_b, valMuRez0 ) ;
	if( f2* f0 < 0.)
	{
	   valMu1 = valMuRez0;

	}
	else
	{
      valMu0 = valMuRez0;
    }
	 valMuRez = valMuRez0;
  }
	return -100000.;
}

// ���������� ���� ������������ ��� ����������� ��������� � �������� � ������� numPartDiagr
// ����� ��������� ����������� �� 0 �� (mQuantDiagr -1)
double  TAnt_5P10_03::calcRadAngSdvigaPartDiagr__(const int numPartDiagr)
{
	return (- (((double) mQuantDiagr)-1.)/2. + (double)numPartDiagr) *   mAngSdvig;
}

// ���������� ���� ������������ ��� ����������� ��������� � ���������� ���� � ������� numPartDiagr
// ����� ��������� ����������� �� 0 �� (mQuantDiagr -1)
double  TAnt_5P10_03::calGenAngSdvigaPartDiagr__(const int numPartDiagr)
{
	double valRad = calcRadAngSdvigaPartDiagr__( numPartDiagr);
	return transformAngToGeneralizedAng__ (valRad  ) ;;
}

// ���������� ���� ������������ ��� ���������� ��������� � �������� � ������� numSumDiagr
// ����� ��������� ����������� �� 0 �� (mQuantDiagr -2)
double  TAnt_5P10_03::calcRadAngSdvigaSumDiagr__(const int numSumDiagr)
{
    // ���� ����������� ��� ������� ����������� ���������
	double valPart0_Rad = calcRadAngSdvigaPartDiagr__( 0);
	///
	// ���� ����������� ��� ������� ���������  ���������
	double valSum0_Rad =valPart0_Rad + mAngSdvig/ 2.;
	///
	return valSum0_Rad + ((double)numSumDiagr) *   mAngSdvig;
}

// ���������� ���� ������������ ��� ��������� ��������� � ���������� ���� � ������� numSumDiagr
// ����� ��������� ����������� �� 0 �� (mQuantDiagr -2)
double  TAnt_5P10_03::calGenAngSdvigaSumDiagr__(const int numSumDiagr)
{
	double valRad = calcRadAngSdvigaSumDiagr__( numSumDiagr);
	return transformAngToGeneralizedAng__ (valRad  ) ;
}




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////  ����� ���������  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////



// ���������� ���������� ����� ���� � �������� �� ���������� 3 ����������� ��������
// ������� ���������
//!!!!!!!!!!!!!   �������� �������  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//INPUT:
// cmparrS[6] -  ������ �������
// NumRayTriple-  ����� ������ ������� ����� (��������� ����� ���� ����� ������� � ����)
//                ����� ������ ��� ���������� ����� ������� ����
//  OUTPUT:
//  *valEstAngTarg -   ���� ����, ���
// *valEstAngAntp -    ���� ��������, ���
// cmpKTarg - ���� ��������� ����
// cmpKAntp - ���� ��������� ��������
// arrMtrxCorr  - ������ ������� ������ ��������� ���� ���� � ��������
int TAnt_5P10_03::tabulatedSolution_5P10_03__( TComp *cmparrS
  , int iNumRayTriple , double *valEstAngTarg, double *valEstAngAntp, TComp *cmpKTarg , TComp *cmpKAntp , double *arrMtrxCorr)
{
	double valGenAngSdvig = transformAngToGeneralizedAng__ (mAngSdvig ) ;
	// ���������� ���������� ����� ���� � ��������
	double  valGenTargEps = -100.,  valGenAntpEps = -100.0;  // ����� ���� ���� � ��������
	double pval_b0 =0., pval_b1 =0.;
	int irez = tabulatedEstimationGenAngs_5P10_03__(  &cmparrS[iNumRayTriple]
	,  &valGenTargEps,  &valGenAntpEps);
	if (irez < 1)
	{
	   return -1;
	}
	//
	// ���������� ������� ���������
	TComp cmparrA[4],cmparrAInv[4], cmparrK[2], cmparrSTemp[2];

	cmparrA[0] = TComp(fncDiagrSinx_div_x(valGenTargEps + valGenAngSdvig), 0.);
	cmparrA[1] = TComp(fncDiagrSinx_div_x (valGenAntpEps + valGenAngSdvig ), 0.);
	cmparrA[2] = TComp(fncDiagrSinx_div_x (valGenTargEps - valGenAngSdvig ), 0.);
	cmparrA[3] = TComp(fncDiagrSinx_div_x (valGenAntpEps - valGenAngSdvig ), 0.);

	bool brez =   InverseMtrx2( cmparrA, cmparrAInv);
	if (!brez)
	{
	 return -2;
	}

	cmparrSTemp[0] = cmparrS [iNumRayTriple ];
	cmparrSTemp[1] = cmparrS [iNumRayTriple + 2 ];
	MtrxMultMatrx(cmparrAInv,2, 2, cmparrSTemp,1, cmparrK);
	*cmpKTarg = cmparrK[0] ;
	*cmpKAntp = cmparrK[1] ;

	///

	// ���������� ����� ���� � �������� � ���
		  // �������������� ���� ����������� ���� ������ � ��������
			double valRadCentreRayPos = calcRadAngSdvigaPartDiagr__(iNumRayTriple + 1);
			// ( - ((double)mQuantDiagr -1.) /2. + 1. + ((double)iNumRayTriple))				*mAngSdvig ;
		  // ���������� ���� ���� � �������� � ���
			const double ValRadTargEps =  transformGeneralizedAngToAng__( valGenTargEps ) ;// � �������
			*valEstAngTarg = ValRadTargEps + valRadCentreRayPos;
		  // ���������� ���� �������� � �������� � ����
			const double ValRadAntp�Eps =  transformGeneralizedAngToAng__ (valGenAntpEps ) ;// � �������
			*valEstAngAntp = ValRadAntp�Eps + valRadCentreRayPos;
	///

	// ���������� ������ ������� ������ ����������� ���������� ����� ���� � ��������
	double arrMtrxCorrMu[4] = {0.};



   calcMtrxCorrGenAngs_AlternativeMeth3PartDiagr__(cmparrS,   valGenTargEps,  valGenAntpEps
  , *cmpKTarg , *cmpKAntp ,iNumRayTriple, arrMtrxCorrMu );

	///

	// ���������� �������� ������� ������ ����������� ����� ���� � �������� � ���������
	double arrQ[4] = {0.}, arrTemp[4] = {0.};
	double coef = mLambda / mAppert / M_PI;
	arrQ[0] =  coef / sqrt( 1.- coef * valGenTargEps * coef * valGenTargEps);
	arrQ[3] =  coef / sqrt( 1.- coef * valGenAntpEps * coef * valGenAntpEps);
	MtrxMultMatrx__(arrQ,2, 2, arrMtrxCorrMu,2, arrTemp)  ;
	MtrxMultMatrxTransp__(arrTemp,2, 2, arrQ,2, arrMtrxCorr) ;
	///
	return irez;
}
//
// ���������� ���������� ����� ���� � �������� �� ���������� 3 ����������� ��������
//INPUT:
// cmparrS[3] -  ������ ������� 3 ���������������� ��������
// ��������� ���� ����� �����

//  OUTPUT:
//  *pvalGenTargEps -   ���� ����, ����� �����
// *pvalGenAntpEps -    ���� ��������, ����� �����

//��� ��������
//-1 -  ������� �� ������������� ��������
int TAnt_5P10_03::tabulatedEstimationGenAngs_5P10_03__(TComp *cmparrS ,  double  *pvalGenTargEps
   ,double  *pvalGenAntpEps)
{

  // ����������� ������� b
   double arr_b[2] = {0.}  ;
   bool brez = calcVect_b__(cmparrS,arr_b) ;
   ///

   // ����������� ����������� ������� �������, ��������� � 0
   double valFreq= VAL_C * 100./ mLambda/ 1000000000.;
	int  numFreq = -1;
   for (int i = 0; i < (QUantCasesHyperbolicType__ + QUantCasesEllipticType__); i++)
   {
	 if (fabs(valFreq - constArrMtrxTransf__[i * 5]) < 0.01)
	 {
	   numFreq = i;
	   break;
	 }
   }
   if (numFreq == -1)
   {
	return -1;
   }
   ///
  // �������� �������������� ������� b
	// �����
  double arrbT[2] = {0.}, arrTransf_b[2] = {0.};
  arrbT[0] = arr_b[0] + constArrPointsSdvig__[ 3 * numFreq +1];
  arrbT[1] = arr_b[1] + constArrPointsSdvig__[ 3 * numFreq +2];
	//��������
	double arrMtrxTransf[4] ;
	memcpy(arrMtrxTransf, &constArrMtrxTransf__[5 * numFreq +1], 4 * sizeof(double));
  MtrxMultMatrx__(arrMtrxTransf, 2, 2, arrbT, 1, arrTransf_b) ;
 ///

 //
	double arrGenAngs[2] ={0.};
	int ireturn = 0;
 if (mLambda < 3.14)
 { // ������������� ������
	// ���������� ����� ����� ������� ������, ���������� ����� ����� arrTransf_b
	//   ���������  ����������
	double valFi0 = atan2(arrTransf_b[1], arrTransf_b[0]);
	double valD = sqrt(arrTransf_b[0] * arrTransf_b[0] + arrTransf_b[1] * arrTransf_b[1]);
	if (valD < 1.)
	{
	  return 0;
	}
	double gam  =  acos(1. / valD);
	double arrt[2]= {0.};
	arrt[0] = valFi0  + gam;
	arrt[1] = valFi0  - gam;
	fncCorrectArg__(arrt);
	fncCorrectArg__(&arrt[1]);

	///

	//


	ireturn = 0;
	// ���������� ������ ������ � ��������    ������������� ������//  ������� ������������� ��������� ������� 6
	int iNumRow1 = findNumRow__(constArrCoeffPolinEllipticCase__,  QUantCasesEllipticType__
		, 4 + LEnArrCoefPolinEllipticCase__,  valFreq) ;
	for (int i =0; i < 2; i++)
	{
	  if (fabs(arrt[ireturn])> constArrCoeffPolinEllipticCase__[iNumRow1 * (4 + LEnArrCoefPolinEllipticCase__) +1])
	  {
		 continue;
	  }
		 double temp = fabs( arrt[i]);
		 double valPolinom = fncPolinom__(&constArrCoeffPolinEllipticCase__[iNumRow1 * (4 + LEnArrCoefPolinEllipticCase__) + 4]
		  , LEnArrCoefPolinEllipticCase__, temp);
		 if (fabs(valPolinom) > M_PI/2.)
		 {
		   int iii =0;
		 }
		 arrGenAngs [ireturn] = -tan(valPolinom)*constArrCoeffPolinEllipticCase__[iNumRow1 * 11 + 3];

	 if (arrt[i] >0.)
	 {
	   arrGenAngs[ireturn] = -arrGenAngs [ireturn];
	 }
	 ireturn++;

	}



 }
 else
 {
   // ��������������� ������
	double pntOutX[2]  ={0.}, pntOutY [2] = {0.};

   int irez = findTangencyHyperbolaPoints(arrTransf_b[0], arrTransf_b[1]
 ,pntOutX[0], pntOutY[0], pntOutX[1], pntOutY[1]);
   if (irez == 0)
   {
	 return 0;
   }
	ireturn = 0;
	for (int i = 0; i < 2; i++)
	{
	  bool brez = findGenAng_UsingApproximation_ForHyperbolicCase__(pntOutX[i],pntOutY[i], &arrGenAngs[ireturn]);

	  if(brez)
	  {
		ireturn++ ;
	  }
	}
 }

   if (ireturn ==2)
   {
	 *pvalGenTargEps =  arrGenAngs[0];
	*pvalGenAntpEps =  arrGenAngs[1];
	if ((*pvalGenTargEps) < (*pvalGenAntpEps) )
	 {
	 double temp = *pvalGenAntpEps;
	 (*pvalGenAntpEps) = (*pvalGenTargEps);
	 (*pvalGenTargEps) = temp;
	 }
   }

   if (ireturn == 1)
   {
        *pvalGenTargEps =  arrGenAngs[0];
   }

   return ireturn ;
}


///


bool TAnt_5P10_03::findGenAng_UsingApproximation_ForHyperbolicCase__(double pntInpX, double pntInpY,  double *pGenAng)
{
	double valFreq= VAL_C * 100./ mLambda/ 1000000000.;
	double tang = pntInpY/pntInpX;
	double valt = log((1.+ tang)/(1. - tang))/ 2.;
	if (pntInpX <=0.)
	{// ����� ����� ���������

	 int iNumRow0 = findNumRow__(constArrCoeffPolinHypCaseLeft__, QUantCasesHyperbolicType__, LEnArrCoefPolinHypCaseLeft__ + 1,  valFreq) ;
	 double *arrCoeff = new double  [LEnArrCoefPolinHypCaseLeft__];
	 memcpy(arrCoeff, &constArrCoeffPolinHypCaseLeft__[(LEnArrCoefPolinHypCaseLeft__ + 1)*iNumRow0 +1]
		, LEnArrCoefPolinHypCaseLeft__ * sizeof(double));
	 *pGenAng = fncOddPolinom__(arrCoeff, LEnArrCoefPolinHypCaseLeft__, valt);
	 delete arrCoeff ;
	 return true;
	}
	else
	{ // ������ ����� ���������
	  int iNumRow0 = findNumRow__(constArrLoranCoeffPolinHypCaseRight__, QUantCasesHyperbolicType__, LEnArrLoranCoeff__ + 1,  valFreq) ;
	   double *arrCoeff = new double  [LEnArrLoranCoeff__];
	   memcpy(arrCoeff, &constArrLoranCoeffPolinHypCaseRight__[(LEnArrLoranCoeff__ + 1)*iNumRow0 +1], LEnArrLoranCoeff__ * sizeof(double));

	   int iNumRow1= findNumRow__(constArrAreaCoeffHypCaseRight__, QUantCasesHyperbolicType__,5,  valFreq) ;
	   double bound = constArrAreaCoeffHypCaseRight__[iNumRow1 * 5 + 2];
	   if (valt < bound)
	   {
		*pGenAng = fncPolinom__(arrCoeff, LEnArrLoranCoeff__, valt) / valt/ valt;
		 delete arrCoeff ;
		 return true;
	   }
	   else
	   {
		   if (valt > -bound)
		   {
			 *pGenAng = -fncPolinom__(arrCoeff, LEnArrLoranCoeff__, -valt) / valt/ valt;
			 delete arrCoeff ;
			 return true;
		   }
		   else
		   {
			 delete arrCoeff ;
			 return false;
		   }
	   }
	  delete arrCoeff ;
	}
}

///
int  findNumRow__(const double *arr, int  numRows, int  numCols, double valFreq)
{
for (int i = 0; i < numRows; i++)
   {
	 if (fabs(valFreq - arr[i * numCols]) < 0.01)
	 {
	   return i;
	   break;
	 }
   }
 return -1;
}
///


void fncCorrectArg__(double *t)
{
	 if ((*t) > M_PI)
	{
	  (*t)  -= 2. *M_PI;
	}

	if ((*t)  < -M_PI)
	{
	  (*t)  += 2. *M_PI;
	}

}
///
double fncPoinomApproximatedEstimationEllipticCase(double valArgTemp ,double valSlop
	 ,double *arrCoeff, int lenArrCoeff)
{
  if (valArgTemp >0)
	  {
	   return tan(fncPolinom__(arrCoeff, lenArrCoeff, valArgTemp))*valSlop;
	  }
	  else
	  {
	   return -tan(fncPolinom__(arrCoeff, lenArrCoeff, valArgTemp))*valSlop;
	  }
}



// ���������� ����� ������� ������ ����������� �� ��������� ������  pntInpX, pntInpY
// � ��������� x*x -y*y = 1
// pntInp0X, pntInp0Y -  �������� ����� �� ���������
// OUTPUT:
// pntOut0X, pntOut0Y,  pntOut1X, pntOut1Y -   ����� �������
// ��������� �-�� ����� �������
//
int findTangencyHyperbolaPoints(double  pntInp0X,double  pntInp0Y
 ,double &pntOut0X,double &pntOut0Y, double &pntOut1X,double &pntOut1Y)
{


  ///
  if ((pntInp0X * pntInp0X - pntInp0Y * pntInp0Y -1.) > 0.)
  {
	return 0;
  }
  // ���������� ����� ������� ��������� ��������� x*x - y*y =1
  double vala = pntInp0Y *  pntInp0Y  -   pntInp0X * pntInp0X;
  double valb = 2. * pntInp0X;
  double valc = -(1. +  pntInp0Y *  pntInp0Y);
  TComp x1,x2;
  // ������� ����������� ��������� a*x*x + b*x + c =0
// ����������
// 0 - 2 �������������� ��������� �����
// 1 - 2 �������������� ������� �����
// 2- ������� �� ������� ���� ���� ������� ������ a!= 0, c=0
// 3 - 2 ���������� ����������� �����
// 4 -  1 �������������� ������ (a =0)
// 5  - �������������� a=b=0, c!=0
// 6 - ��������� )a=b=c=0
int irez =  SolvEq2__(vala, valb, valc,x1,x2) ;
switch (irez)
{
	case 0:
	case 2:
	 pntOut0X =  x1.m_Re;
	 pntOut0Y = sqrt( pntOut0X * pntOut0X  -1);

	 pntOut1X =  x2.m_Re;
	 pntOut1Y = sqrt( pntOut1X * pntOut1X  -1);

	 return 2;
	break;
	case 1:
	 pntOut0X =  x1.m_Re;
	 pntOut0Y = sqrt( pntOut0X * pntOut0X  -1);

	 return 1;
	default:
	return 0;
	break;
}
}

  // ���������� �������� �������������� ������� � ��������� ���������
// ������� � 1
 // arrCoeff - ������ �������������
 // arrCoeff [0] - ������� ��� x, arrCoeff [1]  - ������� ��� x*x*x, � �.�.
 double fncOddPolinom__(double *arrCoeff, int lenarr, double valArg)
 {

   double valTemp =  valArg;
   double valSum = 0.;
   for (int i = 0; i < lenarr; i++)
   {
	 valSum +=  valTemp *  arrCoeff[i];
	 valTemp = valTemp * valArg * valArg;
   }
   return valSum;
 }

 ///


 // ���������� �������� �������������� �������
 // arrCoeff - ������ �������������
 // arrCoeff [0] - ������� ��� 1, arrCoeff [1]  - ������� ��� x,arrCoeff [2]  ��� x*x � �.�.
 double fncPolinom__(double *arrCoeff, int lenarrCoeff, double valArg)
 {
   double valTemp =  1.;
   double valSum = 0.;
   for (int i = 0; i < lenarrCoeff; i++)
   {
	 valSum +=  valTemp *  arrCoeff[i];
	 valTemp = valTemp * valArg ;
   }
   return valSum;
 }


// ��������� ������� 2-�� �������
// Returns:
// true , ���� det(A) != 0
// false - � ��������� ������
bool   InverseMtrx2__(double *arrA, double *arrOut)
{  double EPS = 0.00000000001;
   double det = arrA [ 0] * arrA [ 3] - arrA [ 1] * arrA [ 2] ;
   if (fabs(det) < EPS) return false ;
   arrOut [ 0] =  arrA [ 3] / det ;
   arrOut [ 1] = -arrA [ 1] / det ;
   arrOut [ 2] = -arrA [ 2] / det ;
   arrOut [ 3] =  arrA [ 0] / det ;
  return true;
}

// ��������� ������� �� �������
void MtrxMultMatrx__(double *parrA,int nRowsA, int nColsA, double * parrB,int nColsB, double *parrRez)
{
  memset (parrRez,0,sizeof(double) * nRowsA *  nColsB) ;
 double *pRez = parrRez ;
  double *pA = parrA ;
  double *pB = parrB ;

  for (int i = 0; i < nRowsA; i++)
  {
	double *pAt = pA;
   double *pBt =pB ;
  for (int j = 0; j < nColsB; j++)
  {
	pAt = pA ;
	pBt = parrB + j ;

	for (int k = 0; k < nColsA ; k++)
	{
	  *pRez += (*pAt) * (*pBt);
	  pAt ++;
	  pBt += nColsB ;
	}
	pRez++;
  }
  pA += nColsA ;
  }
}

// ��������� ������� �� ����������������� �������
void MtrxMultMatrxTransp__(double *parrA,int nRowsA, int nColsA, double * parrB,int nRowsB, double *parrRez)
{
   memset (parrRez,0,sizeof(double) * nRowsA *  nRowsB) ;
  double *pRez = parrRez ;
  double *pA = parrA ;
  double *pB = parrB ;

  for (int i = 0; i < nRowsA; i++)
  {
	double *pAt = pA;
   double *pBt =pB ;
  for (int j = 0; j < nRowsB; j++)
  {
	pAt = pA ;


	for (int k = 0; k < nColsA ; k++)
	{
	  *pRez += (*pAt) * (*pBt);
	  pAt ++;
	  pBt ++;
	}
	pRez++;
  }
  pA += nColsA ;
  }
}


void calc_dB0_po_dS__(TComp *cmparrS, double valGenAngSdvig, double val_b0, double val_b1, double *arrdB_po_dS)
{


 double temp = cmparrS[0].m_Re * cmparrS[2].m_Im - cmparrS[2].m_Re * cmparrS[0].m_Im;
 arrdB_po_dS[0] = -val_b0 * cmparrS[2].m_Im/ temp;
 arrdB_po_dS[1] =  val_b0 * cmparrS[2].m_Re/ temp;
 arrdB_po_dS[2] =  cmparrS[2].m_Im/ temp;
 arrdB_po_dS[3] =  -cmparrS[2].m_Re/ temp;
 arrdB_po_dS[4] =  -cmparrS[2].m_Im * val_b1/ temp;
 arrdB_po_dS[5] =  cmparrS[2].m_Re * val_b1/ temp;

}


void calc_dB1_po_dS__(TComp *cmparrS, double valGenAngSdvig, double val_b0, double val_b1, double *arrdB_po_dS)
{
 double temp = cmparrS[0].m_Re * cmparrS[2].m_Im - cmparrS[2].m_Re * cmparrS[0].m_Im;
 arrdB_po_dS[0] = val_b0 * cmparrS[0].m_Im/ temp;
 arrdB_po_dS[1] =  -val_b0 * cmparrS[0].m_Re/ temp;
 arrdB_po_dS[2] = - cmparrS[0].m_Im/ temp;
 arrdB_po_dS[3] =  cmparrS[0].m_Re/ temp;
 arrdB_po_dS[4] =  cmparrS[0].m_Im * val_b1/ temp;
 arrdB_po_dS[5] =  -cmparrS[0].m_Re * val_b1/ temp;
}


// ������� ����������� ��������� a*x*x + b*x + c =0
// ����������
// 0 - 2 �������������� ��������� �����
// 1 - 2 �������������� ������� �����
// 2- ������� �� ������� ���� ���� ������� ������ a!= 0, c=0
// 3 - 2 ���������� ����������� �����
// 4 -  1 �������������� ������ (a =0)
// 5  - �������������� a=b=0, c!=0
// 6 - ��������� )a=b=c=0
int SolvEq2__(const double a,const double b,const double c,TComp &x1,TComp &x2)
{
	x1.m_Re = 0;
	x1.m_Im = 0;
	x2.m_Re = 0;
	x2.m_Im = 0;
	if (fabs(a) < 0.000000000001)
	{
	  if (fabs(b) < 0.000000000001)
	  {
		if (fabs(c) < 0.000000000001) return 6; //��������� a=b=c=0

		return 5 ;  // �������������� a=b=0, c!=0
	  }
	  x1.m_Re = -c/b ;  // 1 �������������� ������ (a =0)
	  return 4 ;
	}
	if (fabs(c) < 0.000000000001)
	{
	  x1.m_Re = - b/a;
	 return 2;       // ������� �� ������� ���� ���� ������� ������ a!= 0, c=0
	}

	if(fabs(c) > DBL_MAX/10)
	{
	   //	ShowMessage(L"fabs(c) > DBL_MAX/10");
		return 10;
    }
	double d2 = b * b - 4 * c * a ;
	if (d2 > DBL_MIN)
	{
	 double d = sqrt(d2);
	 x1.m_Re = (-b + d)/ 2/a ;
	 x2.m_Re = (-b - d)/ 2/a ;
	 return 0;

	}
	if (d2 < - DBL_MIN)
	{
	  x1.m_Re = x2.m_Re = -b/2/a ;
	  x1.m_Im = sqrt(-d2);
	  x2.m_Im = -x1.m_Im;
	  return 3;
	}
	x1.m_Re = x2.m_Re =  -b/2/a ;
	return 1;

}

#pragma package(smart_init)
