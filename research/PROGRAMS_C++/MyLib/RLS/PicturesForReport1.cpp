//---------------------------------------------------------------------------


#pragma hdrstop
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dir.h>
#include "SeaPlg.h"
#include "PicturesForReport1.h"
#include "URPolygon.h"
#include "UrPointXY.h"
#include "UABPlg.h"
#include "CapMathewVess.h"
#include "URPolyLine.h"
#include "YrWriteShapeFile.h"
// ���������� ���������������  ��������  ��� ������ ���� 1
// �������� 4 ��������� ���������
// ������� �����������
// ������ ���, ���������� ���, �����������
// ����, �������
void  _fastcall createPicture1(wchar_t *wchFoldName )
{
   // �������� ���� ����� ����������� ����
	wchar_t wchFileSeaSinName[300] ={0};
	wcscpy(  wchFileSeaSinName,  wchFoldName);
	wcscat(wchFileSeaSinName, L"\\SeaSinName.shp");
	const double valAmpl = 5./*20.*/, valOmega =0.05, valPh0 =0.,  valXMin = -1000.
	 ,valXMax =10000., valDeep =-10000.;
	 const int numPoints = 100000;
	drawSeaSin( wchFileSeaSinName, valAmpl , valOmega,  valPh0,  valXMin
   ,  valXMax,  valDeep,  numPoints );
   ///

   // �������� ���� ����� ����
   TURPointXY pntTagr(3000., 300.);
	 wchar_t wchFileTargName[300] ={0};
	wcscpy( wchFileTargName,  wchFoldName);
	wcscat(wchFileTargName, L"\\Targ.shp");

	 const double  valAngRotation = M_PI;
	 const double valRastigenie = 14.;
	 drawUAB(wchFileTargName,  valAngRotation, pntTagr,valRastigenie) ;
   ///

	 // �������� ���� ����� �������

	wchar_t wchFileVessPath[300] ={0};
	wcscpy( wchFileVessPath,  wchFoldName);
	wcscat(wchFileVessPath, L"\\Vessel");
	_wmkdir(wchFileVessPath);
	const TURPointXY pntSdvig(0.,0.);
	const double  valAngPsi = -0.02;
	const double  valAngEps = 0.03;
	const double  valCentreVessX = 0.;
	double valRast = 1.;//3.;
	TURPointXY pntRadar(0.,0.), pntEndCentralAxe(0.,0.);
	drawVesselMathewPicture( wchFileVessPath,  valAngPsi,  valAngEps
	 ,   valCentreVessX, valRast, pntRadar, pntEndCentralAxe);
	 ///

	// ���������� ����� ������� ����������� � ������� ���������� ���������  pntMirrow
	TURPointXY pntMirrow(0.,0.), pntAntp;
	bool brez =   findMirrowPnt(valAmpl , valOmega,  valPh0,pntRadar,pntTagr,pntTagr.X -200., &pntMirrow);
	///


// �������� ������  ���� ������ ��� ����� �������������� ���������� �����
		// ������ ������� �����������
		//������ ����������� �����������
		// ����� �������
		// ����� ���� - �����
		// ����� ���� - �������
		// �����  �����  - �������
		// ����� ��������
	wchar_t wchFoldName1[300] ={0};

	wcscpy( wchFoldName1,  wchFoldName);
	wcscat(wchFoldName1, L"\\Set1");
	_wmkdir(wchFoldName1);
	createSetGraphs(wchFoldName1, valAmpl ,  valOmega,   valPh0, pntRadar, pntTagr, pntMirrow, pntAntp) ;
 ///

 // �������� ���� ����� ��������
	wchar_t wchFileAntpName[300] ={0};
	wcscpy( wchFileAntpName,  wchFoldName);
	wcscat(wchFileAntpName, L"\\Antp.shp");

	 drawUAB(wchFileAntpName,  valAngRotation, pntAntp,valRastigenie) ;
   ///

   // �������� ���� ����� ����������� ��� ���������
	TURPolyLine plnCentral( pntRadar, pntEndCentralAxe) ;
	wcscpy(wchFileAntpName, wchFoldName);
	wcscat(wchFileAntpName, L"\\CentralLine.shp");
	 plnCentral.WriteSetSHPFiles(wchFileAntpName, & plnCentral ,1);
	 ///

   //�������� ���� ����� ��� ������ ����
  wcscpy(wchFileAntpName, wchFoldName);
	wcscat(wchFileAntpName, L"\\Sect1.shp");
 TYrWriteShapeFile::CreateAngleMarks(wchFileAntpName, pntRadar
	 , pntEndCentralAxe, pntTagr,500.,520. ) ;
 ///

 //�������� ���� ����� ��� ������ ��������
  wcscpy(wchFileAntpName, wchFoldName);
	wcscat(wchFileAntpName, L"\\Sect2.shp");
 TYrWriteShapeFile::CreateAngleMarks(wchFileAntpName, pntRadar
	 , pntEndCentralAxe, pntAntp,600.,620. ) ;
 ///


// ���������� ������ ����� ������� ����������� � ������� ���������� ���������  pntMirrow

  //	brez =   findMirrowPnt(valAmpl , valOmega,  valPh0,pntRadar,pntTagr,pntTagr.X -5000., &pntMirrow);
  brez =   findMirrowPnt(valAmpl , valOmega,  valPh0,pntRadar,pntTagr,pntTagr.X -1200., &pntMirrow);
	///

	// �������� ������  ���� ������ ��� ����� �������������� ���������� �����
		// ������ ������� �����������
		//������ ����������� �����������
		// ����� �������
		// ����� ���� - �����
		// ����� ���� - �������
		// �����  �����  - �������
		// ����� ��������
	wchar_t wchFoldName2[300] ={0};

	wcscpy( wchFoldName2,  wchFoldName);
	wcscat(wchFoldName2, L"\\Set2");
	_wmkdir(wchFoldName2);
	createSetGraphs(wchFoldName2, valAmpl ,  valOmega
	,   valPh0, pntRadar, pntTagr, pntMirrow, pntAntp) ;

	///

	// ���������� ���� ���������� �����
	 TURPointXY arrPntAntipods[10000];
	 int quantAntipodPoints = 0;
	 findAllMirrPnts(valAmpl , valOmega,  valPh0,pntRadar,pntTagr,arrPntAntipods, quantAntipodPoints );
	 wchar_t wchFilePoints[300] ={0};

	wcscpy( wchFilePoints,  wchFoldName);
	wcscat(wchFilePoints, L"\\AntipodPointsArr.shp");
	TURPointXY::WriteSetSHPFiles(wchFilePoints,arrPntAntipods, quantAntipodPoints) ;
}

 // �������������
// ���������� ���������������  ��������  ��� ������ ���� 1
// �������� 4 ��������� ���������
// ������� �����������
// ������ ���, ���������� ���, �����������
// ����, �������
// valVessRast  - ������� ���������� �������
void  _fastcall createPicture1(wchar_t *wchFoldName, double valVessRast )
{
   // �������� ���� ����� ����������� ����
	wchar_t wchFileSeaSinName[300] ={0};
	wcscpy(  wchFileSeaSinName,  wchFoldName);
	wcscat(wchFileSeaSinName, L"\\SeaSinName.shp");
	const double valAmpl = 5./*20.*/, valOmega =0.05, valPh0 =0.,  valXMin = -1000.
	 ,valXMax =10000., valDeep =-10000.;
	 const int numPoints = 100000;
	drawSeaSin( wchFileSeaSinName, valAmpl , valOmega,  valPh0,  valXMin
   ,  valXMax,  valDeep,  numPoints );
   ///

   // �������� ���� ����� ����
   TURPointXY pntTagr(3000., 300.);
	 wchar_t wchFileTargName[300] ={0};
	wcscpy( wchFileTargName,  wchFoldName);
	wcscat(wchFileTargName, L"\\Targ.shp");

	 const double  valAngRotation = M_PI;
	 const double valRastigenie = 14.;
	 drawUAB(wchFileTargName,  valAngRotation, pntTagr,valRastigenie) ;
   ///

	 // �������� ���� ����� �������

	wchar_t wchFileVessPath[300] ={0};
	wcscpy( wchFileVessPath,  wchFoldName);
	wcscat(wchFileVessPath, L"\\Vessel");
	_wmkdir(wchFileVessPath);
	const TURPointXY pntSdvig(0.,0.);
	const double  valAngPsi = -0.02;
	const double  valAngEps = 0.03;
	const double  valCentreVessX = 0.;

	TURPointXY pntRadar(0.,0.), pntEndCentralAxe(0.,0.);
	drawVesselMathewPicture( wchFileVessPath,  valAngPsi,  valAngEps
	 ,   valCentreVessX, valVessRast, pntRadar, pntEndCentralAxe);
	 ///

	// ���������� ����� ������� ����������� � ������� ���������� ���������  pntMirrow
	TURPointXY pntMirrow(0.,0.), pntAntp;
	bool brez =   findMirrowPnt(valAmpl , valOmega,  valPh0,pntRadar,pntTagr,pntTagr.X -200., &pntMirrow);
	///


// �������� ������  ���� ������ ��� ����� �������������� ���������� �����
		// ������ ������� �����������
		//������ ����������� �����������
		// ����� �������
		// ����� ���� - �����
		// ����� ���� - �������
		// �����  �����  - �������
		// ����� ��������
	wchar_t wchFoldName1[300] ={0};

	wcscpy( wchFoldName1,  wchFoldName);
	wcscat(wchFoldName1, L"\\Set1");
	_wmkdir(wchFoldName1);
	createSetGraphs(wchFoldName1, valAmpl ,  valOmega,   valPh0, pntRadar, pntTagr, pntMirrow, pntAntp) ;
 ///

 // �������� ���� ����� ��������
	wchar_t wchFileAntpName[300] ={0};
	wcscpy( wchFileAntpName,  wchFoldName);
	wcscat(wchFileAntpName, L"\\Antp.shp");

	 drawUAB(wchFileAntpName,  valAngRotation, pntAntp,valRastigenie) ;
   ///

   // �������� ���� ����� ����������� ��� ���������
	TURPolyLine plnCentral( pntRadar, pntEndCentralAxe) ;
	wcscpy(wchFileAntpName, wchFoldName);
	wcscat(wchFileAntpName, L"\\CentralLine.shp");
	 plnCentral.WriteSetSHPFiles(wchFileAntpName, & plnCentral ,1);
	 ///

   //�������� ���� ����� ��� ������ ����
  wcscpy(wchFileAntpName, wchFoldName);
	wcscat(wchFileAntpName, L"\\Sect1.shp");
 TYrWriteShapeFile::CreateAngleMarks(wchFileAntpName, pntRadar
	 , pntEndCentralAxe, pntTagr,500.,520. ) ;
 ///

 //�������� ���� ����� ��� ������ ��������
  wcscpy(wchFileAntpName, wchFoldName);
	wcscat(wchFileAntpName, L"\\Sect2.shp");
 TYrWriteShapeFile::CreateAngleMarks(wchFileAntpName, pntRadar
	 , pntEndCentralAxe, pntAntp,600.,620. ) ;
 ///


// ���������� ������ ����� ������� ����������� � ������� ���������� ���������  pntMirrow

  //	brez =   findMirrowPnt(valAmpl , valOmega,  valPh0,pntRadar,pntTagr,pntTagr.X -5000., &pntMirrow);
  brez =   findMirrowPnt(valAmpl , valOmega,  valPh0,pntRadar,pntTagr,pntTagr.X -1200., &pntMirrow);
	///

	// �������� ������  ���� ������ ��� ����� �������������� ���������� �����
		// ������ ������� �����������
		//������ ����������� �����������
		// ����� �������
		// ����� ���� - �����
		// ����� ���� - �������
		// �����  �����  - �������
		// ����� ��������
	wchar_t wchFoldName2[300] ={0};

	wcscpy( wchFoldName2,  wchFoldName);
	wcscat(wchFoldName2, L"\\Set2");
	_wmkdir(wchFoldName2);
	createSetGraphs(wchFoldName2, valAmpl ,  valOmega
	,   valPh0, pntRadar, pntTagr, pntMirrow, pntAntp) ;

	///

	// ���������� ���� ���������� �����
	 TURPointXY arrPntAntipods[10000];
	 int quantAntipodPoints = 0;
	 findAllMirrPnts(valAmpl , valOmega,  valPh0,pntRadar,pntTagr,arrPntAntipods, quantAntipodPoints );
	 wchar_t wchFilePoints[300] ={0};

	wcscpy( wchFilePoints,  wchFoldName);
	wcscat(wchFilePoints, L"\\AntipodPointsArr.shp");
	TURPointXY::WriteSetSHPFiles(wchFilePoints,arrPntAntipods, quantAntipodPoints) ;
}

// ���������� ������ �������� ��� ����� ���������� �����
// ������ �������
//������ �����������
// ����� �������
// ����� ���� - �����
// ����� ���� - �������
 // �����  �����  - �������
  // ����� ��������
void   createSetGraphs(wchar_t *wchFoldName,const double valAmpl ,const double  valOmega
	, const double  valPh0,const TURPointXY pntRadar, const TURPointXY pntTagr
	, TURPointXY pntMirrow, TURPointXY &pntAntipod)
{

	// ������ �������
	TURPointXY pntN(- valAmpl*valOmega* cos(valOmega * pntMirrow.X + valPh0), 1.);
	TURPointXY::DoNorm(pntN) ;
	double valL = 500.;
	TURPointXY pointEnd(pntMirrow.X + valL* pntN.X,pntMirrow.Y + valL* pntN.Y);
	TURPolyLine arrowN = TURPolyLine::fncCreateArrow(pntMirrow, pointEnd
					,70.,10. * M_PI/ 180.);
	 wchar_t wchFileNormArrName[300] ={0};
	wcscpy( wchFileNormArrName,  wchFoldName);
	wcscat(wchFileNormArrName, L"\\NormArrow.shp");
	arrowN.WriteSetSHPFiles(wchFileNormArrName, &arrowN, 1);
	//������ �����������
	TURPointXY vectCasat(pntN.Y, -pntN.X);
	TURPointXY pointBegin(pntMirrow.X + valL* vectCasat.X,pntMirrow.Y + valL* vectCasat.Y);
	pointEnd = TURPointXY(pntMirrow.X - valL* vectCasat.X,pntMirrow.Y - valL* vectCasat.Y);

	wchar_t wchFileCasatName[300] ={0};
	wcscpy( wchFileCasatName,  wchFoldName);
	wcscat(wchFileCasatName, L"\\Casat.shp");
	TURPolyLine plnCasat = TURPolyLine::fncCreateArrow(pointBegin, pointEnd
					,0.01,10. * M_PI/ 180.);
	plnCasat.WriteSetSHPFiles(wchFileCasatName, &plnCasat, 1);

	// ����� �������
	wchar_t wchFileMirrowName[300] ={0};
	wcscpy( wchFileMirrowName,  wchFoldName);
	wcscat(wchFileMirrowName, L"\\MirrowPnt.shp");
	pntMirrow.WriteSetSHPFiles(wchFileMirrowName, &pntMirrow, 1);
	// ����� ���� - �����
	wchar_t wchFileplnTarg_Radar[300] ={0};
	wcscpy( wchFileplnTarg_Radar,  wchFoldName);
	wcscat(wchFileplnTarg_Radar, L"\\Targ_radar.shp");
	TURPolyLine plnTarg_Radar = TURPolyLine::fncCreateArrow(pntTagr, pntRadar
					,0.01,0.01 * M_PI/ 180.);
	plnTarg_Radar.WriteSetSHPFiles(wchFileplnTarg_Radar, &plnTarg_Radar, 1);

	//
	// ����� ���� - �������
	wchar_t wchFileplnTarg_Mirrow[300] ={0};
	wcscpy( wchFileplnTarg_Mirrow,  wchFoldName);
	wcscat(wchFileplnTarg_Mirrow, L"\\Targ_Mirrow.shp");
	TURPolyLine plnTarg_Mirrow = TURPolyLine::fncCreateArrow(pntTagr, pntMirrow
					,0.01,0.01 * M_PI/ 180.);
	plnTarg_Mirrow.WriteSetSHPFiles(wchFileplnTarg_Mirrow, &plnTarg_Mirrow, 1);
	//
	// ���������� ����� ��������
	  // ���������� �� ������ �� �������
	double val_l0 = TURPointXY::dist( pntMirrow, pntRadar) ;
	  // ���������� �� ���� �� �������
	double val_l1 = TURPointXY::dist(pntMirrow, pntTagr) ;
	  // ���������� �� ������ �� ��������
	double val_l =  val_l0 +  val_l1;
	  // ����� ��������
	 pntAntipod =  TURPointXY::ParamPoint(pntRadar,pntMirrow,val_l /val_l0 ); // p1 + alf*(p2-p1)
	  // �����  �����  - �������
	wchar_t wchFileplnAntip_Radar[300] ={0};
	wcscpy(  wchFileplnAntip_Radar,  wchFoldName);
	wcscat( wchFileplnAntip_Radar, L"\\Antip_radar.shp");
	TURPolyLine plnAtip_Radar = TURPolyLine::fncCreateArrow(pntAntipod, pntRadar
					,0.01,0.01 * M_PI/ 180.);
	plnAtip_Radar.WriteSetSHPFiles( wchFileplnAntip_Radar, &plnAtip_Radar, 1);

	// ����� ��������
	wchar_t wchFileAntipodName[300] ={0};
	wcscpy( wchFileAntipodName,  wchFoldName);
	wcscat(wchFileAntipodName, L"\\AntipodPnt.shp");
	pntAntipod.WriteSetSHPFiles(wchFileMirrowName, &pntAntipod, 1);
}


// ���������� �������� ��� ����� ���������� �����
void   calcAntipod(const double valAmpl ,const double  valOmega
	, const double  valPh0,const TURPointXY pntRadar, const TURPointXY pntTagr
	, TURPointXY pntMirrow, TURPointXY *pntAntipod)
{

	// ������ �������
	TURPointXY pntN(- valAmpl*valOmega* cos(valOmega * pntMirrow.X + valPh0), 1.);
	TURPointXY::DoNorm(pntN) ;
	double valL = 500.;
	TURPointXY pointEnd(pntMirrow.X + valL* pntN.X,pntMirrow.Y + valL* pntN.Y);
	TURPolyLine arrowN = TURPolyLine::fncCreateArrow(pntMirrow, pointEnd
					,70.,10. * M_PI/ 180.);

	//������ �����������
	TURPointXY vectCasat(pntN.Y, -pntN.X);
	TURPointXY pointBegin(pntMirrow.X + valL* vectCasat.X,pntMirrow.Y + valL* vectCasat.Y);
	pointEnd = TURPointXY(pntMirrow.X - valL* vectCasat.X,pntMirrow.Y - valL* vectCasat.Y);

	//

	// ���������� ����� ��������
	  // ���������� �� ������ �� �������
	double val_l0 = TURPointXY::dist( pntMirrow, pntRadar) ;
	  // ���������� �� ���� �� �������
	double val_l1 = TURPointXY::dist(pntMirrow, pntTagr) ;
	  // ���������� �� ������ �� ��������
	double val_l =  val_l0 +  val_l1;
	  // ����� ��������
	*pntAntipod =  TURPointXY::ParamPoint(pntRadar,pntMirrow,val_l /val_l0 ); // p1 + alf*(p2-p1)


}

// ���������� ����� ������� ��������� ��������
// ������� ����������� - ���������
// INPUT:
// valAmpl, valOmega , valPh0 - ��������, �������, ��� ���� ������� �����
// pntRadar - ����� ������������� ������
// pntTagr - ����� ����
// OUTPUT:
// pntTagr - ����� �������
bool findMirrowPnt(const double valAmpl ,const double  valOmega
	, const double  valPh0,const TURPointXY pntRadar, const TURPointXY pntTagr
	, const double valx0, TURPointXY *pntMirrow)
{

  for (int i =0; i < 10000000; i++)
  {
	double valx = valx0 - ((double)i) * 0.1;
	if (valx < (pntRadar.X + 400.))
	{
	   break;
	}
	TURPointXY pntCur(valx,valAmpl * sin(valOmega * valx + valPh0));
	TURPointXY pntTarg0(pntTagr.X -  pntCur.X,pntTagr.Y -  pntCur.Y);
	TURPointXY pntRadar0(pntRadar.X -  pntCur.X,pntRadar.Y -  pntCur.Y);
	TURPointXY pntN(- valAmpl * valOmega * cos(valOmega * valx + valPh0), 1.);
	double angTarg  = fabs(TURPointXY::calcAng( pntTarg0,  pntN)) ;
	double angRadar = fabs(TURPointXY::calcAng( pntRadar0, pntN)) ;
	if (fabs(angTarg - angRadar)< 2. * M_PI/ 180.)
	{
	  *pntMirrow = pntCur;
	  return true;
	}

  }
  return  false;
}
 void	 findAllMirrPnts(const double valAmpl ,const double  valOmega
	, const double  valPh0,const TURPointXY pntRadar, const TURPointXY pntTagr
	, TURPointXY *arrPntAntipod, int &quantMirrPoints)
{
  double valx0 =   pntTagr.X - 100.;
  quantMirrPoints = 0;
  TURPointXY pntMirrow(0.,0.);
   for (int i = 0; i < 10000000; i++)
   {
	 if (quantMirrPoints >100000) break;
	 if(!findMirrowPnt(valAmpl , valOmega, valPh0, pntRadar,  pntTagr
	,  valx0, &pntMirrow)  ) break  ;
	valx0 =  pntMirrow.X - M_PI/valOmega/100. ;

	 calcAntipod( valAmpl ,  valOmega
	,   valPh0, pntRadar,  pntTagr
	,  pntMirrow, &arrPntAntipod[quantMirrPoints]) ;
	 quantMirrPoints++;


   }
}

#pragma package(smart_init)
