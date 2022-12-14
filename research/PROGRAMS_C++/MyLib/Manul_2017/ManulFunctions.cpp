//---------------------------------------------------------------------------


#pragma hdrstop

#include "ManulFunctions.h"
#include <float.h>
#include <stdlib.h>
#include <math.h>
 #include <string.h>
#include "YrWriteShapeFile.h"
#include "URPolyLine.h"
#include "UrPointXY.h"
#include "URPolygon.h"
#include "Comp.h"
#include "Far.h"
#include "StabSyst2.h"
#include "Gauss.h"

#include "Far_2D.h"

void __fastcall TManulFunctions::createSKZ_Graphs(wchar_t* wchOutFold, const double VAlTargH, const double VAlTargEPR
  , const TEtalonSign ETalonSign ,const double PowerPrd,const double KYPrd, const double VAlHAntenna, const double VAlZAhvatDist
 , TFar_2D Far_2D, const double VAlSigW, const double VAlSigEpsMMO, const double VAlMeasT
  , const double VAlVelocity0)
{
  double arrK[4]  = {1000000,-150000.,-150000.,22500},  arrKExtr[4] = {0.};
  double valAntNormalAng = 0.;
  TFar Far0(Far_2D, true);
  TFar Far(Far0, 4);
  const int INc = ((VAlZAhvatDist -1000. )/ VAlVelocity0 /VAlMeasT);
   TURPolyLine plnMeasSKZ( 1, INc) ;
   TURPolyLine plnExtrErrSKZ( 1, INc) ;
   TURPolyLine plnBoolMSD( 1, INc) ;
  for (int i =0; i < INc; i++)
  {
		 double valHorisDistCur = VAlZAhvatDist - ((double)i) * VAlVelocity0 * VAlMeasT ;
		 double valDistCur =  sqrt(VAlTargH * VAlTargH  + valHorisDistCur * valHorisDistCur);
		double valSigTarg = 0.;

  double valSigQ  =0.;
  double valAntpPhaze = getRand01( ) * 2. * M_PI;



	  bool brez =  	Far_2D.calc_SKZ_LAT(1., valDistCur,VAlTargH,VAlHAntenna, valAntNormalAng , VAlTargEPR
		,ETalonSign  , PowerPrd, KYPrd,  valAntpPhaze, &valSigTarg , &valSigQ  ) * 100. ;
	   plnBoolMSD.Points[i].Y = (brez)?200.:0.;
		TStabSyst2 StabSyst2(VAlMeasT,VAlSigW, VAlSigEpsMMO * valDistCur
	,valSigTarg * valDistCur ) ;
	double arrP [2] ={0.};
		  StabSyst2.OneStepGolubev(arrK, arrP) ;
		  StabSyst2.OneStepGolubevExtrapolation(arrK, VAlMeasT, arrKExtr) ;
		 plnMeasSKZ.Points[i].X = valHorisDistCur;
		 plnExtrErrSKZ.Points[i].X = valHorisDistCur;
		 plnBoolMSD.Points[i].X = valHorisDistCur;
		 plnMeasSKZ.Points[i].Y = valSigTarg * 1000. * 100.;
		 plnExtrErrSKZ.Points[i].Y = sqrt(arrKExtr[0])/ valDistCur * 1000.* 100.;
  }

	wchar_t wchMeasSKZ[300] ={0};
	wcscpy(  wchMeasSKZ,  wchOutFold);
	wcscat(wchMeasSKZ, L"\\plnMeasSKZ.shp");
	plnMeasSKZ.WriteSetSHPFiles(wchMeasSKZ, &plnMeasSKZ, 1) ;

	wchar_t wchExtrErrSKZ[300] ={0};
	wcscpy(  wchExtrErrSKZ,  wchOutFold);
	wcscat(wchExtrErrSKZ, L"\\plnExtrErrSKZ.shp");
	plnExtrErrSKZ.WriteSetSHPFiles(wchExtrErrSKZ, &plnExtrErrSKZ, 1) ;

	wchar_t wchBoolMSD[300] ={0};
	wcscpy( wchBoolMSD,  wchOutFold);
	wcscat(wchBoolMSD, L"\\plnBoolMSD.shp");
	plnBoolMSD.WriteSetSHPFiles(wchBoolMSD, &plnBoolMSD, 1) ;

	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName, wchOutFold);
	wcscat(wchAxesFileName, L"\\AxesArr.shp");
	TYrWriteShapeFile::CreateShpAxes(wchAxesFileName,-100000., 100000., -1000., 10000.) ;

	// ????? 1,5
  TURPolyLine pln1P5( 1, 2) ;
  pln1P5.Points[0].X = 0.;
  pln1P5.Points[0].Y = 1.5 * 100;;
  pln1P5.Points[1].X = 100000.;
  pln1P5.Points[1].Y = 1.5 * 100;;
	wchar_t wchLine1P5[300] ={0};
	wcscpy(  wchLine1P5,  wchOutFold);
	wcscat(wchLine1P5, L"\\plnLine1P5.shp");
	pln1P5.WriteSetSHPFiles(wchLine1P5, &pln1P5, 1) ;

}
#pragma package(smart_init)
