//---------------------------------------------------------------------------


#pragma hdrstop

#include "BlastSituation.h"

#include <math.h>
#include <string.h>
#include "UrPointXY.h"
#include "URPolygon.h"
#include "URPolyLine.h"
#include "MatrixProccess.h"
#include <stdlib.h>
#include <string.h>
#include <vcl.h>
#include "YrWriteShapeFile.h"



TBlastSituation::TBlastSituation()
{
 	//���� ����������� �� ������� ������� � ����
	 mAlfVis = M_PI / 4.;
	// ������ ��������� �
	 mRDefeat = 467.;
	// ���� ��� ������� ������ ����������
	 mFiDefeat = 4.7 / 180. * M_PI;
	// ��� ������ �� ���������
	mSigR = 33.;
	// ��� ������ �� ����  �����������
	mSigAlf = 0.03;
	// ����� ������� ����� ����
	mLDanger = 20.;



}

//---------------------------------------------------------------------------


// ����������� �����������
 TBlastSituation ::TBlastSituation (const TBlastSituation &R)
 {
	 mAlfVis = R.mAlfVis ;
	 mFiDefeat = R.mFiDefeat;
	 mRDefeat = R.mRDefeat;
	 mLDanger = R.mLDanger ;
	 mSigR = R. mSigR ;
	 mSigAlf = R. mSigAlf ;
 }

 // �������� ������������
 TBlastSituation TBlastSituation::operator=(TBlastSituation  R)
 {
	 mAlfVis = R.mAlfVis ;
	 mFiDefeat = R.mFiDefeat;
	 mRDefeat = R.mRDefeat;
	 mLDanger = R.mLDanger ;
	 mSigR = R. mSigR ;
	 mSigAlf = R. mSigAlf ;
	return *this ;
 }
 // ����� �����������
 TBlastSituation::TBlastSituation (const double AlfVis, const double RDefeat,const double FiDefeat
	 ,const double LDanger ,const double SigR ,const double SigAlf)
 {
	mAlfVis = AlfVis;
	mRDefeat = RDefeat ;
	mFiDefeat = FiDefeat ;
	mLDanger = LDanger ;
	mSigR = SigR;
	mSigAlf = SigAlf ;
}

 // ���������� ������� ������������ ����������� ��������� �� ���� ���������
// ��������� ��������� * ���� ��������� = const
// �������� ����� �� ���� ��������� � �����  valAngStep ������� � valAngMin  �� 0.5 *M_PI
// ��������� ����������� � ����� mpwcharrFoldReport
// � ��������� Probab_FROM_FiDefeat.shp
bool  TBlastSituation::fncCreateGraph_Probab_FROM_FiDiagr(wchar_t  *pwchFileName, const double valAngStep,const double valAngMin,const double valAngMax)
{

  double valAngCur = valAngMin;
  int numP = ((valAngMax  - valAngCur)/ valAngStep)  ;
  //numP = 40;
  TURPolyLine  plnGraph  = TURPolyLine(1, numP) ;
  double valDTresh = 0., valProbTresh =0.;
  for (int i = 0; i < numP; i++)

  {
	 valAngCur =  valAngMin + ((double)i) * valAngStep ;

	 double valRDeafCur = mRDefeat * sin(mFiDefeat  / 4.) /sin( valAngCur / 4.) ;


TBlastSituation BlSitCur   ( mAlfVis ,valRDeafCur,valAngCur,mLDanger , mSigR, mSigAlf );

 BlSitCur.findOptBlastTreshold_and_createProbabGraph(valDTresh , valProbTresh, NULL) ;
   ///
	 plnGraph.Points[i].X = valAngCur * 180./ M_PI;
	 plnGraph.Points[i].Y = valProbTresh * 100.;

  }

	  plnGraph.WriteSetSHPFiles(pwchFileName, &plnGraph,1 ) ;
 return true ;
}

// ���������� ����������� ��������� �������
// OUTPUT:
// valDTresh - �����������������������
// valProbTresh - ����������� ����������� ���������
// pwchFileName -����� �������� ������
// ���� pwchFileName == NULL, �� ������ �� ��������

void TBlastSituation::findOptBlastTreshold_and_createProbabGraph(double &valDTresh
	,double &valProbTresh, wchar_t *pwchFileName)
{
   int ii = (mRDefeat * 10./ mSigR);
   double valDMin = mRDefeat - ((double)(ii - 1)) *mSigR /10. ;
   double valDMax =  mRDefeat + 3. * mSigR;
   int iCircle = (valDMax - valDMin)/mSigR * 10.;
   double valDCur = valDMin ;

   TURPointXY *pPoints = new TURPointXY[iCircle];

   valProbTresh = -1.;
   for (int i = 0; i < iCircle; i++)
   {
	  valDCur =  valDMin + ((double)i) *mSigR /10. ;
	  double valProbCur = calcProb( valDCur);
	  if (valProbCur > valProbTresh)
	  {
		valProbTresh  =valProbCur ;
		valDTresh =  valDCur ;
	  }
	  if (pwchFileName)
	  {
		pPoints[i].X = valDCur;
		pPoints[i].Y = valProbCur * 100.;
	  }
   }

   if (pwchFileName)
   {
	 TURPolyLine polyLine( pPoints,iCircle)  ;
	 polyLine.WriteSetSHPFiles(pwchFileName, &polyLine,1) ;

   }
	 delete []pPoints  ;

}

double TBlastSituation::calcProb(const double valDpod)
{

	const int iInt = 100;
	double discr = sqrt(mRDefeat * mRDefeat -mLDanger * mLDanger * sin(mAlfVis)* sin(mAlfVis)/ 4.);
	double val_d1 = -mLDanger * cos(mAlfVis) / 2. - discr;

	  val_d1 = 0.;

	double val_d2 = -mLDanger * cos(mAlfVis) / 2. + discr;
	if (val_d2 < 0.)
	{
	  return 0. ;
	}
	double step = val_d2  / ((double)iInt);
	double sum = 0.;
	for (int i =0; i < iInt; i++)
	{
	 // double tau = val_d1 - valDpod + ((double)i) * step ;
	 double tau =  ((double)i) * step ;
	 // double vald = valDpod + tau ;
	  sum += exp(- (tau-valDpod) * (tau-valDpod) /mSigR / mSigR / 2.)* fncUslovnP(tau) ;
	}
	sum = sum * step /sqrt( M_PI * 2.) / mSigR ;
	return sum;
}

double TBlastSituation::fncUslovnP(const double vald)
{

	double valr1 = sqrt( vald  * vald   + mLDanger *mLDanger/ 4. +  mLDanger * vald * cos ( mAlfVis) );
	double valr2 = sqrt( vald  * vald   + mLDanger *mLDanger/ 4. - mLDanger * vald * cos ( mAlfVis ));
	double valPsi1 = asin( mLDanger * sin(fabs( mAlfVis ))/ valr1/2.) ;
	double valPsi2 = asin( mLDanger * sin(fabs( mAlfVis ))/ valr2/ 2.) ;
	double valLim1 = -mFiDefeat/ 2. +  valPsi1;
	double valLim2 = mFiDefeat/ 2. -  valPsi2;
	double vala = valLim1 / mSigAlf;
	double valb = valLim2 / mSigAlf;
	double p = fncNormRaspr( vala, valb) ;
	return p ;
}

void TBlastSituation::createGrapUslovnP_from_D( wchar_t *pwchFileName)
{
	const int numP = 1000;
	 TURPolyLine  plnGraph  = TURPolyLine(1, numP) ;
	 double step = 1.;
	for (int i =0; i < numP; i++)
	{

	 plnGraph.Points[i].X =  ((double)i) * step ;
	 plnGraph.Points[i].Y =  fncUslovnP(plnGraph.Points[i].X) * 100. ;
	 }
	 plnGraph.WriteSetSHPFiles(pwchFileName, &plnGraph,1) ;
}



double TBlastSituation::fncNormRaspr( const double a, const double b)
{
	if (a > b) return 0;
	if (a > 5. ) return 0.;
	if (b < -5. ) return 0.;
	double vala = ( a < -5.)? -5.:a;
	double valb= ( b > 5.)? 5.:b;
	double valSum = 0.;
	double valSt = 0.001;
	int Num = (valb - vala)/ valSt ;
	double temp = vala;
	for (int i = 0; i < Num; i++)
	{
	  valSum += exp(- temp * temp / 2.);
	  temp +=  valSt;
	}
	return valSum * valSt / sqrt( 2. * M_PI) ;
}

#pragma package(smart_init)



