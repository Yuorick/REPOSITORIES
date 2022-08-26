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
#include "JetPlg.h"

#include "PicturesForReport2.h"

// ���������� ���������������  ��������  ��� ������ ���� 2
// ��������:
// ������� �����������
// �������
// �����
// ���������� �����  (��������)
// ��������� �������� (������ ��������������)
// �������
// ���������
// ��� ����� - �������
// ���� ����� ����� ������
// INPUT:
// valJetAlt0 - ��� ������ ��������
// valJetX0 - ��� ���������� �������� �� ��� X
// valJetVelo - �������� ��������
// valJetVeloAng - ���� ������� �������� �������� � ���������
// valT - �����. ���������� �������� �� ��������� ������� [0;  valT]
void  _fastcall createPicture2(wchar_t *wchFoldName, double valJetAlt0, double valJetX0,  double valJetVelo,  double valJetVeloAng
	, double valT )
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

   // �������� ���� ������ ���������� �������� � �����
   int quantPointsTraj = 1000;
   TURPointXY *parrPntTrajJet = new TURPointXY[ quantPointsTraj] ;
   TURPointXY *parrPntTrajBomb = new TURPointXY[ quantPointsTraj] ;
   double valDelT = valT / ((double) (quantPointsTraj -1));
   for (int i = 0; i < quantPointsTraj; i++)
   {
	 double valTCur = ((int)i) * valDelT;
	 parrPntTrajJet[i] =  TURPointXY(valJetX0 +  valTCur * valJetVelo*cos(valJetVeloAng)
		,valJetAlt0 +  valTCur * valJetVelo*sin(valJetVeloAng));
	 parrPntTrajBomb[i] = TURPointXY(valJetX0 +  valTCur * valJetVelo*cos(valJetVeloAng)
		,valJetAlt0 +  valTCur * valJetVelo*sin(valJetVeloAng) - 9.8 *valTCur * valTCur /2. );
   }
	TURPolyLine plnJetTraj(parrPntTrajJet,quantPointsTraj) ;
	TURPolyLine plnBombTraj(parrPntTrajBomb,quantPointsTraj) ;
	 // ����������  ��������
	wchar_t wchFileName[300] ={0};
	wcscpy(  wchFileName,  wchFoldName);
	wcscat(wchFileName, L"\\TrajJet.shp");
	plnJetTraj.WriteSetSHPFiles(wchFileName, &plnJetTraj, 1);
		// ����������  �����
	wcscpy(  wchFileName,  wchFoldName);
	wcscat(wchFileName, L"\\TrajBomb.shp");
	plnBombTraj.WriteSetSHPFiles(wchFileName, &plnBombTraj, 1);
	///

	// �������� ���� ����� �����
	const double  valAngRotBomb = atan2( valJetVelo*sin(valJetVeloAng) -9.8 *valT , valJetVelo * cos(valJetVeloAng));
	wcscpy(  wchFileName,  wchFoldName);
	wcscat(wchFileName, L"\\Bomb.shp");
	double valRastigenie = 30.;
	drawUAB(wchFileName,  valAngRotBomb, parrPntTrajBomb[quantPointsTraj - 1] ,valRastigenie) ;
   ///

   // �������� ���� ����� ��������
	wcscpy(  wchFileName,  wchFoldName);
	wcscat(wchFileName, L"\\Jet.shp");
	valRastigenie = 5.;
	drawJet(wchFileName ,valJetVeloAng, parrPntTrajJet[quantPointsTraj - 1] , valRastigenie, true);
   ///

	 // �������� ���� ����� �������

	wchar_t wchFileVessPath[300] ={0};
	wcscpy( wchFileVessPath,  wchFoldName);
	wcscat(wchFileVessPath, L"\\Vessel");
	_wmkdir(wchFileVessPath);
	const TURPointXY pntSdvig(0.,0.);
	const double  valAngPsi = 0.;
	const double  valAngEps = atan2(parrPntTrajJet[quantPointsTraj - 1].Y, parrPntTrajJet[quantPointsTraj - 1].X);
	const double  valCentreVessX = 0.;
	double valRast = 3.;
	TURPointXY pntRadar(0.,0.), pntEndCentralAxe(0.,0.);
   double valGiagrR = 6000.;
   double valDiagrWidth = 4./ 180* M_PI;
	drawVesselMathewPicture( wchFileVessPath,  valAngPsi,  valAngEps
	 ,   valCentreVessX, valRast, valGiagrR, valDiagrWidth, pntRadar, pntEndCentralAxe);
	 ///

   // �������� ���� ����� ����������� ��� ���������
	TURPolyLine plnCentral( pntRadar, pntEndCentralAxe) ;
	wcscpy(wchFileName, wchFoldName);
	wcscat(wchFileName, L"\\CentralLine.shp");
	 plnCentral.WriteSetSHPFiles(wchFileName, & plnCentral ,1);
	 ///

	 // �������� ���� ����� ����� ����� - �������
	TURPolyLine plnRadar_Jet( pntRadar, parrPntTrajJet[quantPointsTraj - 1]) ;
	wcscpy(wchFileName, wchFoldName);
	wcscat(wchFileName, L"\\Radar_Jet.shp");
	plnRadar_Jet.WriteSetSHPFiles(wchFileName, &plnRadar_Jet ,1);
	 ///

	  // �������� ���� ����� ����� ����� - �����
	TURPolyLine plnRadar_Bomb( pntRadar, parrPntTrajBomb[quantPointsTraj - 1]) ;
	wcscpy(wchFileName, wchFoldName);
	wcscat(wchFileName, L"\\Radar_Bomb.shp");
	plnRadar_Bomb.WriteSetSHPFiles(wchFileName, &plnRadar_Bomb ,1);
	 ///

   //�������� ���� ����� ��� �������
  wcscpy(wchFileName, wchFoldName);
	wcscat(wchFileName, L"\\Sect1.shp");
	double valDistJet = sqrt( parrPntTrajJet[quantPointsTraj - 1].X *  parrPntTrajJet[quantPointsTraj - 1].X
		 +  parrPntTrajJet[quantPointsTraj - 1].Y *  parrPntTrajJet[quantPointsTraj - 1].Y);
 TYrWriteShapeFile::CreateAngleMarks(wchFileName, pntRadar
	 , parrPntTrajBomb[quantPointsTraj - 1], parrPntTrajJet[quantPointsTraj - 1],2.*valDistJet / 5.,2.*valDistJet / 5.* 1.025 ) ;
 ///

	delete []parrPntTrajJet;
	delete []parrPntTrajBomb;
}

// �������������, � ��� ������ ������ ���� ���������� �������
// ���������� ���������������  ��������  ��� ������ ���� 2
// ��������:
// ������� �����������
// �������
// �����
// ���������� �����  (��������)
// ��������� �������� (������ ��������������)
// �������
// ���������
// ��� ����� - �������
// ���� ����� ����� ������
// INPUT:
// valJetAlt0 - ��� ������ ��������
// valJetX0 - ��� ���������� �������� �� ��� X
// valJetVelo - �������� ��������
// valJetVeloAng - ���� ������� �������� �������� � ���������
// valT - �����. ���������� �������� �� ��������� ������� [0;  valT]
void  _fastcall createPicture2(wchar_t *wchFoldName, double valJetAlt0, double valJetX0,  double valJetVelo,  double valJetVeloAng
	, double valT , double valVessRast )
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

   // �������� ���� ������ ���������� �������� � �����
   int quantPointsTraj = 1000;
   TURPointXY *parrPntTrajJet = new TURPointXY[ quantPointsTraj] ;
   TURPointXY *parrPntTrajBomb = new TURPointXY[ quantPointsTraj] ;
   double valDelT = valT / ((double) (quantPointsTraj -1));
   for (int i = 0; i < quantPointsTraj; i++)
   {
	 double valTCur = ((int)i) * valDelT;
	 parrPntTrajJet[i] =  TURPointXY(valJetX0 +  valTCur * valJetVelo*cos(valJetVeloAng)
		,valJetAlt0 +  valTCur * valJetVelo*sin(valJetVeloAng));
	 parrPntTrajBomb[i] = TURPointXY(valJetX0 +  valTCur * valJetVelo*cos(valJetVeloAng)
		,valJetAlt0 +  valTCur * valJetVelo*sin(valJetVeloAng) - 9.8 *valTCur * valTCur /2. );
   }
	TURPolyLine plnJetTraj(parrPntTrajJet,quantPointsTraj) ;
	TURPolyLine plnBombTraj(parrPntTrajBomb,quantPointsTraj) ;
	 // ����������  ��������
	wchar_t wchFileName[300] ={0};
	wcscpy(  wchFileName,  wchFoldName);
	wcscat(wchFileName, L"\\TrajJet.shp");
	plnJetTraj.WriteSetSHPFiles(wchFileName, &plnJetTraj, 1);
		// ����������  �����
	wcscpy(  wchFileName,  wchFoldName);
	wcscat(wchFileName, L"\\TrajBomb.shp");
	plnBombTraj.WriteSetSHPFiles(wchFileName, &plnBombTraj, 1);
	///

	// �������� ���� ����� �����
	const double  valAngRotBomb = atan2( valJetVelo*sin(valJetVeloAng) -9.8 *valT , valJetVelo * cos(valJetVeloAng));
	wcscpy(  wchFileName,  wchFoldName);
	wcscat(wchFileName, L"\\Bomb.shp");
	double valRastigenie = 30.;
	drawUAB(wchFileName,  valAngRotBomb, parrPntTrajBomb[quantPointsTraj - 1] ,valRastigenie) ;
   ///

   // �������� ���� ����� ��������
	wcscpy(  wchFileName,  wchFoldName);
	wcscat(wchFileName, L"\\Jet.shp");
	valRastigenie = 5.;
	drawJet(wchFileName ,valJetVeloAng, parrPntTrajJet[quantPointsTraj - 1] , valRastigenie, true);
   ///

	 // �������� ���� ����� �������

	wchar_t wchFileVessPath[300] ={0};
	wcscpy( wchFileVessPath,  wchFoldName);
	wcscat(wchFileVessPath, L"\\Vessel");
	_wmkdir(wchFileVessPath);
	const TURPointXY pntSdvig(0.,0.);
	const double  valAngPsi = 0.;
	const double  valAngEps = atan2(parrPntTrajJet[quantPointsTraj - 1].Y, parrPntTrajJet[quantPointsTraj - 1].X);
	const double  valCentreVessX = 0.;
	double valRast = 3.;
	TURPointXY pntRadar(0.,0.), pntEndCentralAxe(0.,0.);
   double valGiagrR = 6000.;
   double valDiagrWidth = 4./ 180* M_PI;
	drawVesselMathewPicture( wchFileVessPath,  valAngPsi,  valAngEps
	 ,   valCentreVessX, valVessRast, valGiagrR, valDiagrWidth, pntRadar, pntEndCentralAxe );
	 ///

   // �������� ���� ����� ����������� ��� ���������
	TURPolyLine plnCentral( pntRadar, pntEndCentralAxe) ;
	wcscpy(wchFileName, wchFoldName);
	wcscat(wchFileName, L"\\CentralLine.shp");
	 plnCentral.WriteSetSHPFiles(wchFileName, & plnCentral ,1);
	 ///

	 // �������� ���� ����� ����� ����� - �������
	TURPolyLine plnRadar_Jet( pntRadar, parrPntTrajJet[quantPointsTraj - 1]) ;
	wcscpy(wchFileName, wchFoldName);
	wcscat(wchFileName, L"\\Radar_Jet.shp");
	plnRadar_Jet.WriteSetSHPFiles(wchFileName, &plnRadar_Jet ,1);
	 ///

	  // �������� ���� ����� ����� ����� - �����
	TURPolyLine plnRadar_Bomb( pntRadar, parrPntTrajBomb[quantPointsTraj - 1]) ;
	wcscpy(wchFileName, wchFoldName);
	wcscat(wchFileName, L"\\Radar_Bomb.shp");
	plnRadar_Bomb.WriteSetSHPFiles(wchFileName, &plnRadar_Bomb ,1);
	 ///

   //�������� ���� ����� ��� �������
  wcscpy(wchFileName, wchFoldName);
	wcscat(wchFileName, L"\\Sect1.shp");
	double valDistJet = sqrt( parrPntTrajJet[quantPointsTraj - 1].X *  parrPntTrajJet[quantPointsTraj - 1].X
		 +  parrPntTrajJet[quantPointsTraj - 1].Y *  parrPntTrajJet[quantPointsTraj - 1].Y);
 TYrWriteShapeFile::CreateAngleMarks(wchFileName, pntRadar
	 , parrPntTrajBomb[quantPointsTraj - 1], parrPntTrajJet[quantPointsTraj - 1],2.*valDistJet / 5.,2.*valDistJet / 5.* 1.025 ) ;
 ///

	delete []parrPntTrajJet;
	delete []parrPntTrajBomb;
}


#pragma package(smart_init)
